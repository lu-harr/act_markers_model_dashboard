# Years represented by the 29 layers in every model raster.
EXPECTED_YEARS <- 2000:2028

# Single source of truth for model labels, file locations, and palette family.
MODEL_CATALOG <- data.frame(
  model_id = c(
    "k13_all", "k13_A675V", "k13_C469Y", "k13_P441L", "k13_R561H",
    "k13_R622I", "crt76", "mdr1_86", "mdr1_184", "mdr1_1246"
  ),
  model_label = c(
    "Kelch 13 (all ART-associated mutations)",
    "Kelch 13 A675V",
    "Kelch 13 C469Y",
    "Kelch 13 P441L",
    "Kelch 13 R561H",
    "Kelch 13 R622I",
    "Pfcrt K76T",
    "Pfmdr1 N86Y",
    "Pfmdr1 Y184F",
    "Pfmdr1 D1246Y"
  ),
  legend_note = c(
    "0 = 0% ART-associated mutations; 0.6 = 60% ART-associated mutations",
    "0 = 0% 675V; 0.6 = 60% 675V",
    "0 = 0% 469Y; 0.6 = 60% 469Y",
    "0 = 0% 441L; 0.6 = 60% 441L",
    "0 = 0% 561H; 0.6 = 60% 561H",
    "0 = 0% 622I; 0.6 = 60% 622I",
    "0 = 100% K76; 1 = 100% 76T",
    "0 = 100% N86; 1 = 100% 86Y",
    "0 = 100% Y184; 1 = 100% 184F",
    "0 = 100% D1246; 1 = 100% 1246Y"
  ),
  relative_path = file.path(
    "outputs",
    c(
      "k13_all", "k13_A675V", "k13_C469Y", "k13_P441L", "k13_R561H",
      "k13_R622I", "crt76", "mdr1_86", "mdr1_184", "mdr1_1246"
    ),
    "preds_medians.tif"
  ),
  is_k13 = c(rep(TRUE, 6), rep(FALSE, 4)),
  stringsAsFactors = FALSE
)

# Packages required when the interactive dashboard starts.
REQUIRED_PACKAGES <- c(
  "shiny", "bslib", "leaflet", "terra", "sf", "rnaturalearth",
  "rnaturalearthdata", "viridisLite", "htmltools", "iddoPal", "dplyr"
)

# Initial camera position; users can zoom out to the unrestricted global map.
INITIAL_MAP_VIEW <- c(lng = 18, lat = 1, zoom = 3)

# Countries intentionally omitted from shapes, selectors, and summaries.
NORTHERN_AFR_ISO3 <- c("DZA", "EGY", "ESH", "LBY", "MAR", "TUN")

# Fixed domains ensure model maps remain comparable across years and SNPs.
K13_COLOUR_DOMAIN <- c(min = 0, max = 0.6)
K13_LEGEND_BREAKS <- c(0.02, 0.1, 0.2, 0.4, 0.6)
K13_LEGEND_LABELS <- c("0.02", "0.1", "0.2", "0.4", "0.6")
PARTNER_DRUG_COLOUR_DOMAIN <- c(min = 0, max = 1)
PARTNER_DRUG_LEGEND_LABELS <- c("0", "0.5", "1")
