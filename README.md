# ACT molecular markers model dashboard

~Courtesy of ChatGPT 5.6~

An R Shiny and Leaflet dashboard for exploring predictions described in the ten spatiotemporal models of molecular markers of artemisinin and ACT partner drug resistance/reduced susceptibility/drug pressure described in Harrison et al. (<https://doi.org/10.64898/2026.03.03.26347488>).

## Requirements

Use R 4.1 or newer. Install the required packages once:

```r
install.packages(c(
  "shiny", "bslib", "leaflet", "terra", "sf", "rnaturalearth",
  "rnaturalearthdata", "viridisLite", "htmltools", "iddoPal"
))
```

The app uses public basemap tiles and Google Fonts, so those visual features require an internet connection. The pinned `html2canvas` 1.4.1 dependency is bundled under its MIT licence in `www/`, allowing the PNG-export logic itself to load without a CDN.

## Precompute country summaries

From the project root, run:

```sh
Rscript R/precompute_country_metadata.R
```

This validates every model stack and writes `outputs/country_prediction_metadata.csv`. It calculates an area-weighted median using the fraction of each boundary cell covered by the country. Algeria, Egypt, Libya, Morocco, Tunisia, and Western Sahara are deliberately excluded from both the output and the dashboard boundary layer. The calculation reads large raster stacks and can take several minutes.

The existing CSV is replaced only after the new output has been fully calculated and validated.

## Run the dashboard

From the project root:

```r
shiny::runApp()
```

The dashboard validates all raster stacks and the country-metadata CSV at startup. If the CSV is missing or invalid, predictions still render, but country tooltips report that no prediction is available and the sidebar shows the preprocessing command.
