EXPECTED_YEARS <- 2000:2028

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

REQUIRED_PACKAGES <- c(
  "shiny", "bslib", "leaflet", "terra", "sf", "rnaturalearth",
  "rnaturalearthdata", "viridisLite", "htmltools"
)

AFRICA_BOUNDS <- c(west = -20, south = -36, east = 55, north = 39)

NORTHERN_AFR_ISO3 <- c("DZA", "EGY", "ESH", "LBY", "MAR", "TUN")
K13_COLOUR_DOMAIN <- c(min = 0, max = 0.6)
K13_LEGEND_BREAKS <- c(0.02, 0.1, 0.2, 0.4, 0.6)
K13_LEGEND_LABELS <- c("0.02", "0.1", "0.2", "0.4", "0.6")

MAJOR_CITIES <- data.frame(
  name = c(
    "Cairo", "Lagos", "Kinshasa", "Johannesburg", "Nairobi", "Addis Ababa",
    "Dar es Salaam", "Luanda", "Khartoum", "Abidjan", "Accra", "Dakar",
    "Kampala", "Maputo", "Antananarivo", "Casablanca"
  ),
  lon = c(31.24, 3.38, 15.31, 28.05, 36.82, 38.76, 39.21, 13.23, 32.56, -4.01, -0.19, -17.45, 32.58, 32.59, 47.51, -7.59),
  lat = c(30.04, 6.52, -4.32, -26.20, -1.29, 9.01, -6.79, -8.84, 15.50, 5.36, 5.56, 14.69, 0.35, -25.97, -18.88, 33.57),
  stringsAsFactors = FALSE
)

LANDMARKS <- data.frame(
  name = c(
    "Pyramids of Giza", "Victoria Falls", "Mount Kilimanjaro", "Serengeti",
    "Okavango Delta", "Table Mountain", "Timbuktu", "Virunga National Park"
  ),
  lon = c(31.13, 25.86, 37.36, 34.83, 22.91, 18.41, -3.00, 29.20),
  lat = c(29.98, -17.92, -3.07, -2.33, -19.28, -33.96, 16.77, -0.70),
  stringsAsFactors = FALSE
)
