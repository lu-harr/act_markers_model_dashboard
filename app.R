source(file.path("R", "config.R"), local = FALSE)
source(file.path("R", "helpers.R"), local = FALSE)

countries <- load_african_countries()
model_stacks <- load_model_stacks(check_values = TRUE)

model_domains <- setNames(vector("list", nrow(MODEL_CATALOG)), MODEL_CATALOG$model_id)
for (i in seq_len(nrow(MODEL_CATALOG))) {
  model_domains[[MODEL_CATALOG$model_id[[i]]]] <- model_domain(
    model_stacks[[MODEL_CATALOG$model_id[[i]]]], MODEL_CATALOG$is_k13[[i]]
  )
}

metadata_state <- tryCatch(
  list(
    data = load_country_metadata(
      file.path("outputs", "country_prediction_metadata.csv"), countries
    ),
    error = NULL
  ),
  error = function(error) list(data = NULL, error = conditionMessage(error))
)

model_choices <- stats::setNames(MODEL_CATALOG$model_id, MODEL_CATALOG$model_label)
country_choices <- c(
  "All Africa" = "ALL",
  stats::setNames(countries$country_iso3, countries$country_name)
)

theme <- bslib::bs_theme(
  version = 5,
  bg = "#f5f4ef",
  fg = "#17231f",
  primary = "#176b57",
  secondary = "#d8e2dc",
  base_font = bslib::font_google("DM Sans"),
  heading_font = bslib::font_google("DM Sans")
)

ui <- shiny::fluidPage(
  theme = theme,
  shiny::tags$head(
    shiny::tags$meta(name = "viewport", content = "width=device-width, initial-scale=1"),
    shiny::includeCSS(file.path("www", "styles.css")),
    shiny::includeScript(file.path("www", "html2canvas.min.js")),
    shiny::includeScript(file.path("www", "html2canvas-browser-bridge.js")),
    shiny::includeScript(file.path("www", "map-screenshot.js"))
  ),
  shiny::div(
    class = "app-shell",
    shiny::tags$aside(
      class = "control-panel",
      shiny::div(
        class = "brand-block",
        shiny::div(class = "eyebrow", "MARCSE · MODEL EXPLORER"),
        shiny::h1("Resistance predictions"),
        shiny::p("Explore modelled mutation prevalence across Africa, 2000–2028.")
      ),
      shiny::div(
        class = "controls",
        shiny::selectInput(
          "model", "Model", choices = model_choices, selected = "k13_all",
          width = "100%"
        ),
        shiny::sliderInput(
          "year", "Prediction year", min = 2000, max = 2028,
          value = 2026, step = 1, sep = "", width = "100%"
        ),
        shiny::selectizeInput(
          "country", "Zoom to country", choices = country_choices,
          selected = "ALL", options = list(placeholder = "Choose a country"),
          width = "100%"
        ),
        shiny::selectInput(
          "basemap", "Basemap",
          choices = c("Light map" = "CartoDB.PositronNoLabels", "Street map" = "OpenStreetMap.Mapnik"),
          selected = "CartoDB.PositronNoLabels", width = "100%"
        ),
        shiny::div(
          class = "overlay-controls",
          shiny::checkboxInput("show_place_labels", "Place labels", value = TRUE),
          shiny::checkboxInput("show_cities", "Major cities", value = FALSE),
          shiny::checkboxInput("show_landmarks", "Landmarks", value = FALSE)
        ),
        shiny::actionButton(
          "download_map", "Save map as PNG",
          icon = shiny::icon("camera"), class = "btn-screenshot", width = "100%"
        )
      ),
      shiny::div(
        class = "panel-footer",
        shiny::uiOutput("metadata_notice"),
        shiny::p("Country values are area-aware medians from the precomputed metadata file.")
      )
    ),
    shiny::div(
      class = "map-panel",
      shiny::div(
        id = "map-shell",
        class = "map-shell",
        leaflet::leafletOutput("prediction_map", width = "100%", height = "100%"),
        shiny::uiOutput("map_status")
      )
    )
  )
)

server <- function(input, output, session) {
  raster_cache <- new.env(parent = emptyenv())
  cache_order <- character()
  cache_limit <- 4L

  get_raster_layer <- function(model_id, year) {
    key <- paste(model_id, year, sep = "|")
    if (exists(key, envir = raster_cache, inherits = FALSE)) {
      cache_order <<- c(setdiff(cache_order, key), key)
      return(get(key, envir = raster_cache, inherits = FALSE))
    }

    raster <- model_stacks[[model_id]]
    years <- extract_layer_years(names(raster))
    layer_index <- which(years == year)
    if (length(layer_index) != 1L) {
      stop("Year ", year, " does not match exactly one layer for ", model_id, ".", call. = FALSE)
    }
    layer <- raster[[layer_index]]
    ranges <- terra::minmax(layer)
    if (!all(is.finite(ranges))) stop("The selected raster layer contains no valid values.", call. = FALSE)
    if (!terra::is.lonlat(layer)) {
      layer <- terra::project(layer, "EPSG:4326", method = "bilinear")
    }

    assign(key, layer, envir = raster_cache)
    cache_order <<- c(cache_order, key)
    while (length(cache_order) > cache_limit) {
      remove(list = cache_order[[1]], envir = raster_cache)
      cache_order <<- cache_order[-1]
    }
    layer
  }

  active_model <- shiny::reactive({
    shiny::req(input$model)
    model_row(input$model)
  })

  tooltip_labels <- function(model_id, year, model_label) {
    values <- rep(NA_real_, nrow(countries))
    if (!is.null(metadata_state$data)) {
      subset <- metadata_state$data[
        metadata_state$data$model_id == model_id & metadata_state$data$year == year,
        , drop = FALSE
      ]
      values <- subset$median_prediction[match(countries$country_iso3, subset$country_iso3)]
    }
    lapply(seq_len(nrow(countries)), function(i) {
      htmltools::HTML(paste0(
        '<div class="country-tooltip"><strong>',
        htmltools::htmlEscape(countries$country_name[[i]]),
        "</strong><span>", htmltools::htmlEscape(model_label),
        " · ", year, "</span><span>Median prediction: <b>",
        format_prediction(values[[i]]), "</b></span></div>"
      ))
    })
  }

  add_country_borders <- function(proxy, model_id, year, model_label) {
    proxy |>
      leaflet::clearGroup("National borders") |>
      leaflet::addPolygons(
        data = countries,
        layerId = ~country_iso3,
        group = "National borders",
        color = "#263d35",
        weight = 1,
        opacity = 0.9,
        fill = TRUE,
        fillColor = "transparent",
        fillOpacity = 0.01,
        label = tooltip_labels(model_id, year, model_label),
        labelOptions = leaflet::labelOptions(
          direction = "auto", opacity = 1, sticky = TRUE,
          className = "country-label"
        ),
        highlightOptions = leaflet::highlightOptions(
          color = "#fff6d5", weight = 3, opacity = 1,
          fillOpacity = 0.08, bringToFront = TRUE
        ),
        options = leaflet::pathOptions(pane = "border-pane")
      )
  }

  output$prediction_map <- leaflet::renderLeaflet({
    leaflet::leaflet(options = leaflet::leafletOptions(
      minZoom = 2, worldCopyJump = TRUE, preferCanvas = TRUE
    )) |>
      leaflet::addMapPane("prediction-pane", zIndex = 250) |>
      leaflet::addMapPane("border-pane", zIndex = 410) |>
      leaflet::addMapPane("label-pane", zIndex = 425) |>
      leaflet::addMapPane("marker-pane", zIndex = 450) |>
      leaflet::addProviderTiles(
        leaflet::providers$CartoDB.PositronNoLabels,
        group = "Basemap",
        options = leaflet::providerTileOptions(noWrap = TRUE)
      ) |>
      leaflet::fitBounds(
        AFRICA_BOUNDS[["west"]], AFRICA_BOUNDS[["south"]],
        AFRICA_BOUNDS[["east"]], AFRICA_BOUNDS[["north"]]
      ) |>
      leaflet::addScaleBar(position = "bottomleft", options = leaflet::scaleBarOptions(imperial = FALSE))
  })

  shiny::observeEvent(input$basemap, {
    provider <- leaflet::providers[[input$basemap]]
    leaflet::leafletProxy("prediction_map") |>
      leaflet::clearGroup("Basemap") |>
      leaflet::addProviderTiles(
        provider, group = "Basemap",
        options = leaflet::providerTileOptions(noWrap = TRUE)
      )
  }, ignoreInit = TRUE)

  shiny::observeEvent(input$show_place_labels, {
    proxy <- leaflet::leafletProxy("prediction_map") |> leaflet::clearGroup("Place labels")
    if (isTRUE(input$show_place_labels)) {
      proxy |>
        leaflet::addProviderTiles(
          leaflet::providers$CartoDB.PositronOnlyLabels,
          group = "Place labels",
          options = leaflet::providerTileOptions(noWrap = TRUE, pane = "label-pane")
        )
    }
  }, ignoreInit = FALSE)

  shiny::observeEvent(list(input$model, input$year), {
    model <- active_model()
    year <- as.integer(input$year)

    tryCatch(
      shiny::withProgress(message = "Loading prediction surface", value = 0.2, {
        layer <- get_raster_layer(model$model_id, year)
        shiny::incProgress(0.35, detail = "Preparing colours")
        domain <- model_domains[[model$model_id]]
        palette <- make_palette(domain, model$is_k13)
        legend <- legend_html(model$model_label, domain, palette, model$is_k13)

        proxy <- leaflet::leafletProxy("prediction_map") |>
          leaflet::clearGroup("Prediction") |>
          leaflet::removeControl("prediction-legend") |>
          leaflet::addRasterImage(
            layer,
            colors = palette,
            opacity = 0.82,
            project = TRUE,
            method = "bilinear",
            group = "Prediction",
            layerId = "active-prediction",
            maxBytes = 24 * 1024^2,
            options = leaflet::gridOptions(pane = "prediction-pane")
          ) |>
          leaflet::addControl(
            html = htmltools::HTML(legend), position = "bottomright",
            layerId = "prediction-legend"
          )

        shiny::incProgress(0.35, detail = "Updating country summaries")
        add_country_borders(proxy, model$model_id, year, model$model_label)
      }),
      error = function(error) {
        shiny::showNotification(conditionMessage(error), type = "error", duration = NULL)
      }
    )
  }, ignoreInit = FALSE)

  shiny::observeEvent(input$country, {
    proxy <- leaflet::leafletProxy("prediction_map")
    if (identical(input$country, "ALL")) {
      proxy |>
        leaflet::fitBounds(
          AFRICA_BOUNDS[["west"]], AFRICA_BOUNDS[["south"]],
          AFRICA_BOUNDS[["east"]], AFRICA_BOUNDS[["north"]]
        )
    } else {
      selected <- countries[countries$country_iso3 == input$country, ]
      shiny::req(nrow(selected) == 1L)
      bounds <- sf::st_bbox(selected)
      proxy |>
        leaflet::fitBounds(bounds[["xmin"]], bounds[["ymin"]], bounds[["xmax"]], bounds[["ymax"]])
    }
  }, ignoreInit = TRUE)

  shiny::observeEvent(input$show_cities, {
    proxy <- leaflet::leafletProxy("prediction_map") |> leaflet::clearGroup("Major cities")
    if (isTRUE(input$show_cities)) {
      proxy |>
        leaflet::addCircleMarkers(
          data = MAJOR_CITIES, lng = ~lon, lat = ~lat, label = ~name,
          group = "Major cities", radius = 3.5, stroke = TRUE, weight = 1,
          color = "#ffffff", fillColor = "#17231f", fillOpacity = 0.95,
          options = leaflet::pathOptions(pane = "marker-pane")
        )
    }
  }, ignoreInit = FALSE)

  shiny::observeEvent(input$show_landmarks, {
    proxy <- leaflet::leafletProxy("prediction_map") |> leaflet::clearGroup("Landmarks")
    if (isTRUE(input$show_landmarks)) {
      proxy |>
        leaflet::addCircleMarkers(
          data = LANDMARKS, lng = ~lon, lat = ~lat, label = ~name,
          group = "Landmarks", radius = 5, stroke = TRUE, weight = 1.5,
          color = "#17231f", fillColor = "#f4c95d", fillOpacity = 1,
          options = leaflet::pathOptions(pane = "marker-pane")
        )
    }
  }, ignoreInit = FALSE)

  output$map_status <- shiny::renderUI({
    model <- active_model()
    shiny::div(
      class = "map-status",
      shiny::span(class = "map-status__dot"),
      shiny::strong(model$model_label),
      shiny::span(" · ", as.integer(input$year))
    )
  })

  output$metadata_notice <- shiny::renderUI({
    if (is.null(metadata_state$error)) return(NULL)
    shiny::div(
      class = "metadata-warning",
      shiny::icon("triangle-exclamation"),
      shiny::span(metadata_state$error)
    )
  })

  shiny::observeEvent(input$screenshot_error, {
    shiny::showNotification(input$screenshot_error, type = "error", duration = 8)
  })
}

shiny::shinyApp(ui, server)
