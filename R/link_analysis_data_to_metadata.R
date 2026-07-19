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
  c("Longitude", "Latitude", "year", "uniq_id_publication")
) |>
  dplyr::mutate(
    uniq_id_publication = as.character(uniq_id_publication),
    Marker = "K76T"
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
  c("Longitude", "Latitude", "year", "PubMedID", "loc")
) |>
  dplyr::mutate(
    PubMedID = clean_pubmed_id(PubMedID),
    Marker = as.character(loc)
  ) |>
  dplyr::left_join(pubmed_article_lookup, by = "PubMedID")

# Recreate the aggregate K13 analysis set while keeping one record per sampled site.
message("Preparing aggregate K13 analysis data...")
k13_source <- read_analysis_set(
  "moldm_marcse_with_markers.csv",
  c("Longitude", "Latitude", "year", "Marker", "Authors", "PubMedID", "mutant", "Present")
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
  )

k13_all <- dplyr::bind_rows(k13_mutants, k13_wildtypes) |>
  dplyr::mutate(Marker = "Any WHO-listed")

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
    c("Longitude", "Latitude", "year", "Authors", "PubMedID")
  ) |>
    dplyr::mutate(
      PubMedID = clean_pubmed_id(PubMedID),
      Marker = marker
    )
})

# Combine every analysis set, retain the requested fields, and shorten author lists.
message("Combining and formatting analysis metadata...")
analysis_sets <- c(list(crt76, pfmdr, k13_all), k13_snps)
analysis_sets <- lapply(analysis_sets, function(data) {
  dplyr::select(data, Longitude, Latitude, year, Marker, Authors, PubMedID)
})

moldm_reorged <- dplyr::bind_rows(analysis_sets) |>
  dplyr::mutate(
    Authors = format_first_author(Authors),
    PubMedID = clean_pubmed_id(PubMedID)
  )

# Write atomically so a failed run cannot leave a partial output CSV.
output_path <- data_path("moldm_reorged.csv")
temporary_path <- tempfile("moldm_reorged_", tmpdir = dirname(output_path), fileext = ".csv")
on.exit(if (file.exists(temporary_path)) unlink(temporary_path), add = TRUE)
utils::write.csv(moldm_reorged, temporary_path, row.names = FALSE, na = "")
if (!file.rename(temporary_path, output_path)) {
  stop("Could not replace output file: ", output_path, call. = FALSE)
}

message(sprintf("Wrote %s rows to %s", nrow(moldm_reorged), output_path))

