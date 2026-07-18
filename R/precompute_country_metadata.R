#!/usr/bin/env Rscript

if (sys.nframe() == 0L) {

script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(script_argument)) stop("Could not determine the script path.", call. = FALSE)
script_path <- normalizePath(sub("^--file=", "", script_argument[[1]]), mustWork = TRUE)
project_root <- dirname(dirname(script_path))
setwd(project_root)

source(file.path(project_root, "R", "config.R"), local = FALSE)
source(file.path(project_root, "R", "helpers.R"), local = FALSE)
assert_packages(c("terra", "sf", "rnaturalearth", "rnaturalearthdata"))

output_path <- file.path(project_root, "outputs", "country_prediction_metadata.csv")

message("Loading and validating model rasters...")
stacks <- load_model_stacks(check_values = TRUE)

reference <- stacks[[1]]
for (model_id in names(stacks)[-1]) {
  same_geometry <- terra::compareGeom(
    reference, stacks[[model_id]],
    lyrs = FALSE, crs = TRUE, ext = TRUE, rowcol = TRUE, res = TRUE,
    stopOnError = FALSE
  )
  if (!isTRUE(same_geometry)) {
    stop(model_id, ": raster geometry differs from the other models.", call. = FALSE)
  }
}

message("Loading African country boundaries...")
countries_wgs84 <- load_african_countries()
countries_raster_crs <- sf::st_transform(countries_wgs84, terra::crs(reference))

message("Calculating cell areas for area-aware medians...")
cell_area_raster <- terra::cellSize(reference[[1]], unit = "km")
cell_areas <- terra::values(cell_area_raster, mat = FALSE)

summarise_country <- function(raster, country, years) {
  names(raster) <- paste0("year_", years)
  layer_columns <- names(raster)
  polygon <- terra::vect(country)

  extracted <- tryCatch(
    terra::extract(raster, polygon, cells = TRUE, exact = TRUE),
    error = function(error) NULL
  )
  coverage_column <- if (!is.null(extracted)) {
    intersect(c("fraction", "weight"), names(extracted))
  } else {
    character()
  }

  if (is.null(extracted) || !length(coverage_column)) {
    extracted <- terra::extract(raster, polygon, cells = TRUE, weights = TRUE)
    coverage_column <- intersect(c("fraction", "weight"), names(extracted))
  }
  if (!length(coverage_column)) {
    stop("terra extraction did not return boundary-cell coverage weights.", call. = FALSE)
  }
  if (!"cell" %in% names(extracted)) {
    stop("terra extraction did not return cell identifiers.", call. = FALSE)
  }

  if (!nrow(extracted)) {
    return(data.frame(
      year = years,
      median_prediction = NA_real_,
      n_valid_cells = 0L
    ))
  }

  coverage <- extracted[[coverage_column[[1]]]]
  area_weights <- coverage * cell_areas[extracted$cell]
  medians <- numeric(length(years))
  counts <- integer(length(years))

  for (i in seq_along(years)) {
    values <- extracted[[layer_columns[[i]]]]
    valid <- is.finite(values)
    counts[[i]] <- sum(valid)
    medians[[i]] <- weighted_median(values, area_weights)
  }

  data.frame(
    year = years,
    median_prediction = medians,
    n_valid_cells = counts
  )
}

results <- vector("list", nrow(MODEL_CATALOG) * nrow(countries_raster_crs))
result_index <- 0L

for (model_index in seq_len(nrow(MODEL_CATALOG))) {
  model <- MODEL_CATALOG[model_index, , drop = FALSE]
  raster <- stacks[[model$model_id]]
  years <- extract_layer_years(names(raster))
  message(sprintf("[%d/%d] %s", model_index, nrow(MODEL_CATALOG), model$model_label))

  for (country_index in seq_len(nrow(countries_raster_crs))) {
    country <- countries_raster_crs[country_index, ]
    summary <- summarise_country(raster, country, years)
    result_index <- result_index + 1L
    results[[result_index]] <- data.frame(
      model_id = model$model_id,
      model_label = model$model_label,
      year = summary$year,
      country_name = country$country_name,
      country_iso3 = country$country_iso3,
      median_prediction = summary$median_prediction,
      n_valid_cells = summary$n_valid_cells,
      stringsAsFactors = FALSE
    )
  }
}

metadata <- do.call(rbind, results)
metadata <- metadata[order(metadata$model_id, metadata$year, metadata$country_iso3), ]
row.names(metadata) <- NULL

key <- paste(metadata$model_id, metadata$year, metadata$country_iso3, sep = "|")
expected_rows <- nrow(MODEL_CATALOG) * length(EXPECTED_YEARS) * nrow(countries_raster_crs)
if (nrow(metadata) != expected_rows) {
  stop("Incomplete output: expected ", expected_rows, " rows but created ", nrow(metadata), ".", call. = FALSE)
}
if (anyDuplicated(key)) stop("Duplicate model/year/country keys were generated.", call. = FALSE)
if (!setequal(metadata$model_id, MODEL_CATALOG$model_id)) stop("Incomplete model coverage.", call. = FALSE)
if (!setequal(metadata$year, EXPECTED_YEARS)) stop("Incomplete year coverage.", call. = FALSE)

dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
temporary_path <- tempfile("country_prediction_metadata_", tmpdir = dirname(output_path), fileext = ".csv")
on.exit(if (file.exists(temporary_path)) unlink(temporary_path), add = TRUE)
utils::write.csv(metadata, temporary_path, row.names = FALSE, na = "NA")

backup_path <- paste0(output_path, ".bak")
if (file.exists(backup_path)) unlink(backup_path)
if (file.exists(output_path) && !file.rename(output_path, backup_path)) {
  stop("Could not prepare the existing metadata CSV for replacement.", call. = FALSE)
}
if (!file.rename(temporary_path, output_path)) {
  if (file.exists(backup_path)) file.rename(backup_path, output_path)
  stop("Could not move the completed metadata CSV into outputs/.", call. = FALSE)
}
if (file.exists(backup_path)) unlink(backup_path)

message(
  "Done: wrote ", nrow(metadata), " rows for ", length(unique(metadata$model_id)),
  " models, ", length(unique(metadata$year)), " years, and ",
  length(unique(metadata$country_iso3)), " countries to ", output_path, "."
)

}
