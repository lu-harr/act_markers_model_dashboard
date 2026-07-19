# Shiny auto-sources files under R/, so run this standalone workflow only via Rscript.
if (sys.nframe() == 0L) {

library(dplyr)

# for some reason Chatty G loves to mess around with paths and roots..
# it's ok sir! I'll behave and run this script from the right place!
data_path <- function(filename) file.path("data", filename)

# Read a CSV without altering source column names, then verify its required schema.
read_analysis_set <- function(filename, required_columns) {
  path <- data_path(filename)
  if (!file.exists(path)) stop("Missing input file: ", path, call. = FALSE)

  data <- utils::read.csv(
    path,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    na.strings = c("", "NA")
  )
  blank_names <- which(is.na(names(data)) | names(data) == "")
  names(data)[blank_names] <- paste0("unnamed_", seq_along(blank_names))
  names(data) <- make.unique(names(data))
  missing_columns <- setdiff(required_columns, names(data))
  if (length(missing_columns)) {
    stop(
      filename, " is missing required columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
  data
}

# Treat blank and literal "null" publication identifiers as missing join keys.
clean_pubmed_id <- function(values) {
  values <- trimws(as.character(values))
  values[is.na(values) | values == "" | tolower(values) == "null"] <- NA_character_
  values
}

# Convert a full comma-separated author list to a compact first-author citation.
format_first_author <- function(authors) {
  authors <- trimws(as.character(authors))
  missing <- is.na(authors) | authors == "" | tolower(authors) == "null"
  first_author <- trimws(sub(",.*$", "", authors))
  citation <- paste0(first_author, " et al.")
  citation[missing] <- NA_character_
  citation
}

# Return the first usable value in a group while preserving a missing value if none exist.
first_non_missing <- function(values) {
  usable <- values[!is.na(values) & trimws(as.character(values)) != ""]
  if (length(usable)) usable[[1]] else NA_character_
}

# Sum or maximise count fields without turning an entirely missing group into zero.
sum_if_present <- function(values) {
  if (all(is.na(values))) NA_real_ else sum(values, na.rm = TRUE)
}

max_if_present <- function(values) {
  if (all(is.na(values))) NA_real_ else max(values, na.rm = TRUE)
}

# Read the complete MOLDM export and retain one row for each unique publication.
message("Reading unique MOLDM articles...")
moldm_articles <- read_analysis_set(
  "moldm_data.csv",
  c("Title", "Journal", "Authors", "Year Published", "PubMedID", "uniq_id_publication")
) |>
  dplyr::select(
    Title, Journal, Authors, `Year Published`, PubMedID, uniq_id_publication
  ) |>
  dplyr::mutate(
    uniq_id_publication = as.character(uniq_id_publication),
    PubMedID = clean_pubmed_id(PubMedID)
  ) |>
  dplyr::distinct()

if (anyDuplicated(moldm_articles$uniq_id_publication)) {
  stop("MOLDM article metadata contains duplicate uniq_id_publication values.", call. = FALSE)
}

# Build one deterministic metadata row per PubMed ID so pfmdr joins cannot multiply rows.
pubmed_article_lookup <- moldm_articles |>
  dplyr::filter(!is.na(PubMedID)) |>
  dplyr::arrange(uniq_id_publication) |>
  dplyr::distinct(PubMedID, .keep_all = TRUE) |>
  dplyr::select(PubMedID, Authors)

# Link crt76 analysis rows to publication metadata using the public publication ID.
message("Linking crt76 analysis data...")
crt76 <- read_analysis_set(
  "moldm_crt76.csv",
  c("Longitude", "Latitude", "year", "Present", "Tested", "uniq_id_publication")
) |>
  dplyr::mutate(
    uniq_id_publication = as.character(uniq_id_publication),
    model = "crt76",
    Marker = "K76T",
    Marker_details = NA_character_
  ) |>
  dplyr::left_join(
    moldm_articles |>
      dplyr::select(uniq_id_publication, Authors, PubMedID),
    by = "uniq_id_publication"
  )

# Link all three pfmdr loci to publication metadata using PubMed ID.
message("Linking pfmdr analysis data...")
pfmdr <- read_analysis_set(
  "pfmdr_single_locus.csv",
  c("Longitude", "Latitude", "year", "Present", "Tested", "PubMedID", "loc")
) |>
  dplyr::mutate(
    PubMedID = clean_pubmed_id(PubMedID),
    Marker = as.character(loc),
    model = dplyr::recode(
      Marker,
      N86Y = "mdr1_86",
      Y184F = "mdr1_184",
      D1246Y = "mdr1_1246"
    ),
    Marker_details = NA_character_
  ) |>
  dplyr::left_join(pubmed_article_lookup, by = "PubMedID")

# Recreate the aggregate K13 analysis set while keeping one record per sampled site.
message("Preparing aggregate K13 analysis data...")
k13_source <- read_analysis_set(
  "moldm_marcse_with_markers.csv",
  c(
    "Longitude", "Latitude", "year", "Marker", "Authors", "PubMedID",
    "mutant", "Present", "Tested"
  )
) |>
  dplyr::mutate(PubMedID = clean_pubmed_id(PubMedID)) |>
  dplyr::filter(mutant | Marker == "wildtype")

k13_mutants <- k13_source |>
  dplyr::filter(mutant, Present > 0)

k13_wildtypes <- k13_source |>
  dplyr::filter(Marker == "wildtype") |>
  dplyr::anti_join(
    k13_mutants |>
      dplyr::distinct(Longitude, Latitude, year, PubMedID),
    by = c("Longitude", "Latitude", "year", "PubMedID")
  ) |>
  dplyr::mutate(
    Marker = "Any WHO-listed",
    Present = 0
  )

# First consolidate repeated rows for the same marker, then calculate site-level totals.
k13_marker_level <- dplyr::bind_rows(k13_mutants, k13_wildtypes) |>
  dplyr::group_by(Longitude, Latitude, year, PubMedID, Marker) |>
  dplyr::summarise(
    Present = sum_if_present(Present),
    Tested = max_if_present(Tested),
    Authors = first_non_missing(Authors),
    .groups = "drop"
  ) |>
  dplyr::arrange(Longitude, Latitude, year, PubMedID, Marker) |>
  dplyr::mutate(
    Marker_detail = paste0(Marker, " - Present: ", Present, "; Tested: ", Tested)
  )

k13_all <- k13_marker_level |>
  dplyr::group_by(Longitude, Latitude, year, PubMedID) |>
  dplyr::summarise(
    model = "k13_all",
    Marker = paste(Marker, collapse = "; "),
    Present = sum_if_present(Present),
    Tested = max_if_present(Tested),
    Authors = first_non_missing(Authors),
    Marker_details = paste(Marker_detail, collapse = " | "),
    .groups = "drop"
  )

# Read each K13 SNP analysis set directly because it already contains article metadata.
k13_snp_files <- c(
  A675V = "moldm_marcse_k13snp_A675V.csv",
  C469Y = "moldm_marcse_k13snp_C469Y.csv",
  P441L = "moldm_marcse_k13snp_P441L.csv",
  R561H = "moldm_marcse_k13snp_R561H.csv",
  R622I = "moldm_marcse_k13snp_R622I.csv"
)

k13_snps <- lapply(names(k13_snp_files), function(marker) {
  read_analysis_set(
    k13_snp_files[[marker]],
    c("Longitude", "Latitude", "year", "Present", "Tested", "Authors", "PubMedID")
  ) |>
    dplyr::mutate(
      PubMedID = clean_pubmed_id(PubMedID),
      model = paste0("k13_", marker),
      Marker = marker,
      Marker_details = NA_character_
    )
})

# Combine every analysis set, calculate prevalence, and shorten author lists.
message("Combining and formatting analysis metadata...")
analysis_sets <- c(list(crt76, pfmdr, k13_all), k13_snps)
analysis_sets <- lapply(analysis_sets, function(data) {
  dplyr::select(
    data, Longitude, Latitude, year, model, Marker, Present, Tested,
    Authors, PubMedID, Marker_details
  )
})

moldm_reorged <- dplyr::bind_rows(analysis_sets) |>
  dplyr::mutate(
    Prevalence = dplyr::if_else(
      !is.na(Tested) & Tested > 0,
      as.numeric(Present) / as.numeric(Tested),
      NA_real_
    ),
    Authors = format_first_author(Authors),
    PubMedID = clean_pubmed_id(PubMedID)
  ) |>
  dplyr::select(
    Longitude, Latitude, year, model, Marker, Present, Tested, Prevalence,
    Authors, PubMedID, Marker_details
  )

# Validate model coverage, point coordinates, ratios, and aggregate-site uniqueness.
expected_models <- c(
  "k13_all", "k13_A675V", "k13_C469Y", "k13_P441L", "k13_R561H",
  "k13_R622I", "crt76", "mdr1_86", "mdr1_184", "mdr1_1246"
)
if (!setequal(unique(moldm_reorged$model), expected_models)) {
  stop("The reorganised data do not contain exactly the ten expected model IDs.", call. = FALSE)
}
if (anyNA(moldm_reorged[c("year", "model")])) {
  stop("The reorganised data contain missing years or model IDs.", call. = FALSE)
}
missing_coordinates <- is.na(moldm_reorged$Longitude) | is.na(moldm_reorged$Latitude)
if (any(missing_coordinates)) {
  message(sprintf(
    "Retaining %d source records with missing coordinates; the dashboard must not map them.",
    sum(missing_coordinates)
  ))
}
valid_tested <- !is.na(moldm_reorged$Tested) & moldm_reorged$Tested > 0
expected_prevalence <- moldm_reorged$Present[valid_tested] / moldm_reorged$Tested[valid_tested]
if (any(abs(moldm_reorged$Prevalence[valid_tested] - expected_prevalence) > 1e-12, na.rm = TRUE)) {
  stop("At least one prevalence does not equal Present / Tested.", call. = FALSE)
}
k13_site_key <- with(
  moldm_reorged[moldm_reorged$model == "k13_all", ],
  paste(Longitude, Latitude, year, PubMedID, sep = "|")
)
if (anyDuplicated(k13_site_key)) {
  stop("Aggregate K13 data contain duplicate site/year/publication records.", call. = FALSE)
}

# Write atomically so a failed run cannot leave a partial output CSV.
output_path <- data_path("moldm_reorged.csv")
temporary_path <- tempfile("moldm_reorged_", tmpdir = dirname(output_path), fileext = ".csv")
on.exit(if (file.exists(temporary_path)) unlink(temporary_path), add = TRUE)
utils::write.csv(moldm_reorged, temporary_path, row.names = FALSE, na = "")
if (!file.rename(temporary_path, output_path)) {
  stop("Could not replace output file: ", output_path, call. = FALSE)
}

message(sprintf("Wrote %s rows to %s", nrow(moldm_reorged), output_path))

}
