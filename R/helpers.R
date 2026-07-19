# Fail early with one actionable message if any runtime packages are unavailable.
assert_packages <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
  if (length(missing)) {
    stop(
      "Missing required R packages: ", paste(missing, collapse = ", "),
      ". Install them before continuing.",
      call. = FALSE
    )
  }
}

# Resolve one catalog record from the stable machine-readable model ID.
model_row <- function(model_id) {
  row <- MODEL_CATALOG[MODEL_CATALOG$model_id == model_id, , drop = FALSE]
  if (nrow(row) != 1L) stop("Unknown model ID: ", model_id, call. = FALSE)
  row
}

# Extract the four-digit year embedded in each raster layer name.
extract_layer_years <- function(layer_names) {
  hits <- regexpr("(?:19|20)[0-9]{2}", layer_names, perl = TRUE)
  if (any(hits < 0L)) {
    stop("Every raster layer name must contain one four-digit year.", call. = FALSE)
  }

  years <- as.integer(regmatches(layer_names, hits))
  remaining <- sub("(?:19|20)[0-9]{2}", "", layer_names, perl = TRUE)
  if (any(grepl("(?:19|20)[0-9]{2}", remaining, perl = TRUE))) {
    stop("A raster layer name contains more than one four-digit year.", call. = FALSE)
  }
  years
}

# Check layer count, year order, CRS, and value availability for one stack.
validate_raster_stack <- function(raster, model_id, check_values = TRUE) {
  if (!inherits(raster, "SpatRaster")) stop(model_id, ": not a SpatRaster.", call. = FALSE)
  if (terra::nlyr(raster) != length(EXPECTED_YEARS)) {
    stop(model_id, ": expected 29 layers but found ", terra::nlyr(raster), ".", call. = FALSE)
  }
  years <- extract_layer_years(names(raster))
  if (anyDuplicated(years) || !identical(years, EXPECTED_YEARS)) {
    stop(
      model_id, ": layer years must be unique and ordered from 2000 through 2028; found ",
      paste(years, collapse = ", "), ".", call. = FALSE
    )
  }
  raster_crs <- terra::crs(raster, proj = TRUE)
  if (!nzchar(raster_crs)) stop(model_id, ": raster CRS is missing.", call. = FALSE)
  if (!terra::hasValues(raster)) stop(model_id, ": raster contains no values.", call. = FALSE)

  if (check_values) {
    ranges <- terra::minmax(raster)
    invalid <- !is.finite(ranges[1, ]) | !is.finite(ranges[2, ])
    if (any(invalid)) {
      stop(model_id, ": all-NA or invalid layer(s): ", paste(years[invalid], collapse = ", "), ".", call. = FALSE)
    }
  }
  invisible(years)
}

# Open all raster stacks lazily after validating their metadata.
load_model_stacks <- function(check_values = TRUE) {
  stacks <- setNames(vector("list", nrow(MODEL_CATALOG)), MODEL_CATALOG$model_id)
  for (i in seq_len(nrow(MODEL_CATALOG))) {
    path <- MODEL_CATALOG$relative_path[[i]]
    if (!file.exists(path)) stop("Missing model file: ", path, call. = FALSE)
    raster <- terra::rast(path)
    validate_raster_stack(raster, MODEL_CATALOG$model_id[[i]], check_values = check_values)
    stacks[[MODEL_CATALOG$model_id[[i]]]] <- raster
  }
  stacks
}

# Load the supported country geometries and establish stable ISO3 join keys.
load_african_countries <- function() {
  countries <- rnaturalearth::ne_countries(
    scale = "medium", continent = "Africa", returnclass = "sf"
  )
  countries <- sf::st_make_valid(countries)

  choose_iso <- function(primary, fallback) {
    primary <- as.character(primary)
    fallback <- as.character(fallback)
    invalid <- is.na(primary) | primary == "" | primary == "-99"
    primary[invalid] <- fallback[invalid]
    primary
  }

  countries$country_iso3 <- choose_iso(countries$iso_a3, countries$adm0_a3)
  countries$country_name <- as.character(countries$name_long)
  countries <- countries[
    !is.na(countries$country_iso3) & grepl("^[A-Z]{3}$", countries$country_iso3),
    c("country_iso3", "country_name", "geometry")
  ]
  countries <- countries[!countries$country_iso3 %in% NORTHERN_AFR_ISO3, ]
  countries <- countries[order(countries$country_name), ]

  if (anyDuplicated(countries$country_iso3)) {
    duplicates <- unique(countries$country_iso3[duplicated(countries$country_iso3)])
    stop("Duplicate African country ISO3 codes: ", paste(duplicates, collapse = ", "), call. = FALSE)
  }
  if (!nrow(countries)) stop("No African country boundaries were returned.", call. = FALSE)
  sf::st_transform(countries, 4326)
}

# Return the fixed colour domain for the selected model family.
model_domain <- function(raster, is_k13) {
  if (is_k13) {
    return(list(min = K13_COLOUR_DOMAIN[["min"]], max = K13_COLOUR_DOMAIN[["max"]]))
  }
  list(min = PARTNER_DRUG_COLOUR_DOMAIN[["min"]], max = PARTNER_DRUG_COLOUR_DOMAIN[["max"]])
}

# Build a clipping palette function while leaving source prediction values intact.
make_palette <- function(domain, is_k13) {
  if (is_k13) {
    sqrt_palette <- leaflet::colorNumeric(
      palette = viridisLite::viridis(256),
      domain = sqrt(c(domain$min, domain$max)),
      na.color = "transparent"
    )
    function(values) {
      missing <- is.na(values)
      adjusted <- pmin(pmax(values, domain$min), domain$max)
      colours <- sqrt_palette(sqrt(adjusted))
      colours[missing] <- "transparent"
      colours
    }
  } else {
    partner_palette <- leaflet::colorNumeric(
      palette = as.vector(iddoPal::iddo_palettes_sequential$BlGyRd),
      domain = c(domain$min, domain$max),
      na.color = "transparent"
    )
    function(values) {
      missing <- is.na(values)
      adjusted <- pmin(pmax(values, domain$min), domain$max)
      colours <- partner_palette(adjusted)
      colours[missing] <- "transparent"
      colours
    }
  }
}

# Format raw prediction values for country tooltips.
format_prediction <- function(value) {
  if (!length(value) || is.na(value) || !is.finite(value)) return("No prediction available")
  absolute <- abs(value)
  if (absolute > 0 && absolute < 0.001) {
    formatC(value, format = "e", digits = 2)
  } else {
    formatC(value, format = "f", digits = 4)
  }
}

# Construct a continuous HTML colour bar with family-specific tick positions.
legend_html <- function(model_label, legend_note, domain, palette, is_k13) {
  gradient_positions <- seq(0, 1, length.out = 64)
  if (is_k13) {
    breaks <- K13_LEGEND_BREAKS
    gradient_values <- (gradient_positions^2) * domain$max
    tick_positions <- sqrt(breaks / domain$max) * 100
  } else {
    breaks <- seq(domain$min, domain$max, length.out = 3)
    gradient_values <- domain$min + gradient_positions * (domain$max - domain$min)
    tick_positions <- seq(0, 100, length.out = length(breaks))
  }

  gradient_stops <- paste0(
    palette(gradient_values), " ", formatC(gradient_positions * 100, format = "f", digits = 1), "%",
    collapse = ""
  )
  gradient_stops <- gsub("%(?=#[0-9A-Fa-f])", "%, ", gradient_stops, perl = TRUE)
  tick_labels <- paste0(
    '<span class="prediction-legend__tick',
    ifelse(tick_positions == 0, " is-start", ifelse(tick_positions == 100, " is-end", "")),
    '" style="left:',
    formatC(tick_positions, format = "f", digits = 2), '%">',
    if (is_k13) K13_LEGEND_LABELS else PARTNER_DRUG_LEGEND_LABELS,
    "</span>",
    collapse = ""
  )

  paste0(
    '<div class="prediction-legend"><div class="prediction-legend__title">',
    htmltools::htmlEscape(model_label),
    '</div><div class="prediction-legend__bar" style="background:linear-gradient(to right,',
    gradient_stops, ')"></div><div class="prediction-legend__ticks">',
    tick_labels, '</div><div class="prediction-legend__note">',
    htmltools::htmlEscape(legend_note),
    "</div></div>"
  )
}

# Read and validate the complete model/year/country metadata grid.
load_country_metadata <- function(path, countries) {
  if (!file.exists(path)) {
    stop(
      "Country metadata is missing. Run: Rscript R/precompute_country_metadata.R",
      call. = FALSE
    )
  }
  metadata <- utils::read.csv(path, stringsAsFactors = FALSE, na.strings = c("NA", ""))
  required <- c(
    "model_id", "model_label", "year", "country_name", "country_iso3",
    "median_prediction", "n_valid_cells"
  )
  absent <- setdiff(required, names(metadata))
  if (length(absent)) stop("Country metadata is missing columns: ", paste(absent, collapse = ", "), call. = FALSE)
  metadata$year <- as.integer(metadata$year)
  metadata$n_valid_cells <- as.integer(metadata$n_valid_cells)

  key <- paste(metadata$model_id, metadata$year, metadata$country_iso3, sep = "|")
  if (anyDuplicated(key)) stop("Country metadata contains duplicate model/year/country keys.", call. = FALSE)

  expected <- expand.grid(
    model_id = MODEL_CATALOG$model_id,
    year = EXPECTED_YEARS,
    country_iso3 = countries$country_iso3,
    stringsAsFactors = FALSE
  )
  expected_key <- paste(expected$model_id, expected$year, expected$country_iso3, sep = "|")
  missing_keys <- setdiff(expected_key, key)
  extra_keys <- setdiff(key, expected_key)
  if (length(missing_keys) || length(extra_keys)) {
    stop(
      "Country metadata coverage is invalid (", length(missing_keys), " missing and ",
      length(extra_keys), " unexpected keys). Re-run preprocessing.", call. = FALSE
    )
  }
  metadata
}

# Load and validate the preprocessed observational records used by the point layer.
load_observation_data <- function(path) {
  if (!file.exists(path)) {
    stop(
      "Observation data are missing. Run: Rscript R/link_analysis_data_to_metadata.R",
      call. = FALSE
    )
  }

  observations <- utils::read.csv(
    path, stringsAsFactors = FALSE, na.strings = c("NA", "")
  )
  required <- c(
    "Longitude", "Latitude", "year", "model", "Marker", "Present", "Tested",
    "Prevalence", "Authors", "PubMedID", "Marker_details"
  )
  absent <- setdiff(required, names(observations))
  if (length(absent)) {
    stop("Observation data are missing columns: ", paste(absent, collapse = ", "), call. = FALSE)
  }

  numeric_columns <- c("Longitude", "Latitude", "year", "Present", "Tested", "Prevalence")
  for (column in numeric_columns) observations[[column]] <- as.numeric(observations[[column]])

  expected_models <- MODEL_CATALOG$model_id
  if (!setequal(unique(observations$model), expected_models)) {
    stop("Observation data do not contain exactly the expected model IDs.", call. = FALSE)
  }
  if (anyNA(observations[c("year", "model")])) {
    stop("Observation data contain missing years or model IDs.", call. = FALSE)
  }
  if (max(observations$year, na.rm = TRUE) != OBSERVATION_MAX_YEAR) {
    stop("The most recent observation year must be ", OBSERVATION_MAX_YEAR, ".", call. = FALSE)
  }

  valid_tested <- is.finite(observations$Tested) & observations$Tested > 0
  expected_prevalence <- observations$Present[valid_tested] / observations$Tested[valid_tested]
  invalid_prevalence <- abs(observations$Prevalence[valid_tested] - expected_prevalence) > 1e-12
  if (any(invalid_prevalence, na.rm = TRUE)) {
    stop("Observation prevalence must equal Present / Tested.", call. = FALSE)
  }

  k13_all <- observations$model == "k13_all"
  k13_key <- with(
    observations[k13_all, ],
    paste(Longitude, Latitude, year, PubMedID, sep = "|")
  )
  if (anyDuplicated(k13_key)) {
    stop("Aggregate K13 observations contain duplicate site records.", call. = FALSE)
  }
  if (any(is.na(observations$Marker_details[k13_all]) | observations$Marker_details[k13_all] == "")) {
    stop("Aggregate K13 observations are missing marker-level details.", call. = FALSE)
  }

  missing_coordinates <- !is.finite(observations$Longitude) | !is.finite(observations$Latitude)
  attr(observations, "missing_coordinate_count") <- sum(missing_coordinates)
  observations
}

# Format integer-like sample counts without unnecessary decimal places.
format_count <- function(value) {
  if (!length(value) || is.na(value) || !is.finite(value)) return("Not available")
  if (abs(value - round(value)) < 1e-9) {
    format(round(value), big.mark = ",", scientific = FALSE, trim = TRUE)
  } else {
    formatC(value, format = "f", digits = 1, big.mark = ",")
  }
}

# Build an area-based radius function and stable representative counts for one model.
observation_size_scale <- function(tested) {
  tested <- tested[is.finite(tested) & tested > 0]
  if (!length(tested)) stop("No positive sample sizes are available for this model.", call. = FALSE)

  maximum <- max(tested)
  radius <- function(values) {
    scaled <- sqrt(pmax(values, 0) / maximum) * OBSERVATION_RADIUS_RANGE[["max"]]
    pmin(
      OBSERVATION_RADIUS_RANGE[["max"]],
      pmax(OBSERVATION_RADIUS_RANGE[["min"]], scaled)
    )
  }
  breaks <- unique(as.numeric(stats::quantile(
    tested, probs = c(0.25, 0.5, 1), names = FALSE, type = 1
  )))
  list(radius = radius, breaks = breaks, maximum = maximum)
}

# Render the custom sample-size legend using the same radii as the Leaflet circles.
observation_size_legend_html <- function(size_scale) {
  rows <- paste0(
    '<div class="size-legend__row"><span class="size-legend__symbol" style="width:',
    formatC(2 * size_scale$radius(size_scale$breaks), format = "f", digits = 1),
    'px;height:',
    formatC(2 * size_scale$radius(size_scale$breaks), format = "f", digits = 1),
    'px"></span><span>',
    vapply(size_scale$breaks, format_count, character(1)),
    "</span></div>",
    collapse = ""
  )
  paste0(
    '<div class="size-legend"><div class="size-legend__title">Number tested</div>',
    rows, "</div>"
  )
}

# Construct a safe point popup, adding marker-level rows for aggregate K13 sites.
observation_popup_html <- function(observation) {
  value <- function(column, fallback = "Not available") {
    candidate <- observation[[column]]
    if (!length(candidate) || is.na(candidate) || !nzchar(trimws(as.character(candidate)))) {
      fallback
    } else {
      as.character(candidate)
    }
  }

  rows <- c(
    paste0("<dt>Marker</dt><dd>", htmltools::htmlEscape(value("Marker")), "</dd>"),
    paste0("<dt>Present</dt><dd>", format_count(observation$Present), "</dd>"),
    paste0("<dt>Tested</dt><dd>", format_count(observation$Tested), "</dd>"),
    paste0("<dt>Year</dt><dd>", as.integer(observation$year), "</dd>"),
    paste0("<dt>Author</dt><dd>", htmltools::htmlEscape(value("Authors")), "</dd>")
  )
  if (!is.na(observation$PubMedID) && nzchar(trimws(as.character(observation$PubMedID)))) {
    rows <- c(
      rows,
      paste0("<dt>PubMedID</dt><dd>", htmltools::htmlEscape(observation$PubMedID), "</dd>")
    )
  }

  marker_details <- ""
  if (identical(observation$model, "k13_all")) {
    details <- strsplit(value("Marker_details"), " | ", fixed = TRUE)[[1]]
    marker_details <- paste0(
      '<div class="observation-popup__details"><strong>Marker-level counts</strong>',
      paste0("<span>", htmltools::htmlEscape(details), "</span>", collapse = ""),
      "</div>"
    )
  }

  htmltools::HTML(paste0(
    '<div class="observation-popup"><dl>', paste0(rows, collapse = ""), "</dl>",
    marker_details, "</div>"
  ))
}

# Calculate the first value at which cumulative area weight reaches 50%.
weighted_median <- function(values, weights) {
  valid <- is.finite(values) & is.finite(weights) & weights > 0
  if (!any(valid)) return(NA_real_)
  values <- values[valid]
  weights <- weights[valid]
  ordering <- order(values)
  values <- values[ordering]
  weights <- weights[ordering]
  values[which(cumsum(weights) >= sum(weights) / 2)[1]]
}
