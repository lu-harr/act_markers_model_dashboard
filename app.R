# Load shared constants and helper functions before constructing the app.
source(file.path("R", "config.R"), local = FALSE)
source(file.path("R", "helpers.R"), local = FALSE)

# Validate static spatial inputs once at process startup.
countries <- load_african_countries()
model_stacks <- load_model_stacks(check_values = TRUE)

# Cache one fixed display domain per model for consistent temporal comparison.
model_domains <- setNames(vector("list", nrow(MODEL_CATALOG)), MODEL_CATALOG$model_id)
for (i in seq_len(nrow(MODEL_CATALOG))) {
  model_domains[[MODEL_CATALOG$model_id[[i]]]] <- model_domain(
    model_stacks[[MODEL_CATALOG$model_id[[i]]]], MODEL_CATALOG$is_k13[[i]]
  )
}

# Country metadata is optional at startup so raster exploration can still load.
metadata_state <- tryCatch(
  list(
    data = load_country_metadata(
      file.path("outputs", "country_prediction_metadata.csv"), countries
    ),
    error = NULL
  ),
  error = function(error) list(data = NULL, error = conditionMessage(error))
)

# Load the prepared observation table once and retain a readable startup error if invalid.
observation_state <- tryCatch(
  {
    data <- load_observation_data(file.path("data", "moldm_reorged.csv"))
    missing_coordinates <- attr(data, "missing_coordinate_count")
    if (missing_coordinates > 0) {
      message(
        "Observation data: ", missing_coordinates,
        " record(s) with missing coordinates will not be mapped."
      )
    }
    list(data = data, error = NULL)
  },
  error = function(error) list(data = NULL, error = conditionMessage(error))
)

# Named vectors give Shiny readable labels while retaining stable IDs as values.
model_choices <- stats::setNames(MODEL_CATALOG$model_id, MODEL_CATALOG$model_label)
country_choices <- c(
  "All Africa" = "ALL",
  stats::setNames(countries$country_iso3, countries$country_name)
)

# CARTO provides a transparent label-only layer that can sit above either basemap.
label_providers <- list(
  "CartoDB.PositronNoLabels" = leaflet::providers$CartoDB.PositronOnlyLabels,
  # bug here ... todo:
  "OpenStreetMap.Mapnik" = leaflet::providers$CartoDB.PositronOnlyLabels
)

# Apply the visual theme consistently to inputs, typography, and notifications.
theme <- bslib::bs_theme(
  version = 5,
  bg = "#F3EEED",
  fg = "#675957",
  primary = "#9D2123",
  secondary = "#D36F5E",
  base_font = bslib::font_google("DM Sans"),
  heading_font = bslib::font_google("DM Sans")
)

# The UI uses a compact control sidebar and a map that fills the remaining area.
ui <- shiny::fluidPage(
  theme = theme,
  shiny::tags$head(
    shiny::tags$meta(name = "viewport", content = "width=device-width, initial-scale=1"),
    shiny::includeCSS(file.path("www", "styles.css")),
    shiny::includeScript(file.path("www", "html2canvas.min.js")),
    shiny::includeScript(file.path("www", "html2canvas-browser-bridge.js")),
    shiny::includeScript(file.path("www", "map-screenshot.js")),
    shiny::includeScript(file.path("www", "observation-slider.js"))
  ),
  shiny::div(
    class = "app-shell",
    shiny::tags$aside(
      class = "control-panel",
      shiny::div(
        class = "brand-block",
        shiny::div(class = "eyebrow", "MARC SE-AFRICA · MODEL EXPLORER"),
        shiny::h1("Molecular surveillance data & model predictions"),
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
        shiny::actionButton(
          "toggle_data", "Show data",
          icon = shiny::icon("circle-dot"), class = "btn-data-toggle", width = "100%"
        ),
        shiny::conditionalPanel(
          condition = "input.toggle_data % 2 === 1",
          class = "observation-year-controls",
          shiny::sliderInput(
            "observation_years", "Observation years",
            min = OBSERVATION_YEAR_SENTINEL, max = OBSERVATION_MAX_YEAR,
            value = c(OBSERVATION_YEAR_SENTINEL, OBSERVATION_MAX_YEAR),
            step = 1, sep = "", width = "100%"
          )
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
          class = "map-label-control",
          shiny::checkboxInput(
            "labels_on_top", "Map labels on top", value = FALSE, width = "100%"
          )
        ),
        shiny::actionButton(
          "download_map", "Save map as PNG",
          icon = shiny::icon("camera"), class = "btn-screenshot", width = "100%"
        )
      ) #,
      # shiny::div(
      #   class = "panel-footer",
      #   shiny::uiOutput("metadata_notice"),
      #   shiny::p("Country values are area-aware medians from the precomputed metadata file.")
      # )
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
  # Keep a small per-session LRU cache to avoid repeatedly reopening recent layers.
  raster_cache <- new.env(parent = emptyenv())
  cache_order <- character()
  cache_limit <- 4L

  # Resolve and, if necessary, reproject one model/year raster layer.
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

  # Expose the selected catalog row to all dependent reactives and observers.
  active_model <- shiny::reactive({
    shiny::req(input$model)
    model_row(input$model)
  })

  # Interpret the action-button count as a persistent Show data / Hide data toggle.
  data_visible <- shiny::reactive({
    toggle_count <- if (is.null(input$toggle_data)) 0L else input$toggle_data
    toggle_count %% 2L == 1L
  })

  shiny::observeEvent(input$toggle_data, {
    visible <- data_visible()
    shiny::updateActionButton(
      session, "toggle_data",
      label = if (visible) "Hide data" else "Show data",
      icon = shiny::icon(if (visible) "circle-xmark" else "circle-dot")
    )
    if (visible && !is.null(observation_state$error)) {
      shiny::showNotification(observation_state$error, type = "error", duration = NULL)
    }
  }, ignoreInit = TRUE)

  # Filter prepared records independently of the selected prediction year.
  active_observations <- shiny::reactive({
    shiny::req(data_visible(), is.null(observation_state$error), input$observation_years)
    years <- as.integer(input$observation_years)
    observations <- observation_state$data
    keep <- observations$model == input$model &
      is.finite(observations$Longitude) & is.finite(observations$Latitude) &
      observations$year <= years[[2]]
    if (years[[1]] != OBSERVATION_YEAR_SENTINEL) {
      keep <- keep & observations$year >= years[[1]]
    }
    observations[keep, , drop = FALSE]
  })

  # Join precomputed medians to map polygons and build safe hover labels.
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

  # Redraw supported national borders above the active prediction raster.
  add_country_borders <- function(proxy, model_id, year, model_label) {
    proxy |>
      leaflet::clearGroup("National borders") |>
      leaflet::addPolygons(
        data = countries,
        layerId = ~country_iso3,
        group = "National borders",
        color = "#675957",
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
          color = "#FDEBE8", weight = 3, opacity = 1,
          fillOpacity = 0.08, bringToFront = TRUE
        ),
        options = leaflet::pathOptions(pane = "border-pane")
      )
  }

  # Build an unrestricted world map, initially centred on Africa.
  output$prediction_map <- leaflet::renderLeaflet({
    leaflet::leaflet(options = leaflet::leafletOptions(
      worldCopyJump = TRUE, preferCanvas = TRUE
    )) |>
      leaflet::addMapPane("prediction-pane", zIndex = 250) |>
      leaflet::addMapPane("border-pane", zIndex = 410) |>
      leaflet::addMapPane("labels", zIndex = 416) |>
      leaflet::addMapPane("observation-pane", zIndex = 440) |>
      leaflet::addProviderTiles(
        leaflet::providers$CartoDB.PositronNoLabels,
        group = "Basemap"
      ) |>
      leaflet::setView(
        lng = INITIAL_MAP_VIEW[["lng"]],
        lat = INITIAL_MAP_VIEW[["lat"]],
        zoom = INITIAL_MAP_VIEW[["zoom"]]
      ) |>
      leaflet::addScaleBar(position = "bottomleft", options = leaflet::scaleBarOptions(imperial = FALSE))
  })

  # Replace only the basemap; the street option explicitly uses OSM Mapnik.
  shiny::observeEvent(input$basemap, {
    provider <- switch(
      input$basemap,
      "OpenStreetMap.Mapnik" = leaflet::providers$OpenStreetMap.Mapnik,
      leaflet::providers$CartoDB.PositronNoLabels
    )
    leaflet::leafletProxy("prediction_map") |>
      leaflet::clearGroup("Basemap") |>
      leaflet::addProviderTiles(provider, group = "Basemap")
  }, ignoreInit = TRUE)

  # Optionally redraw a transparent place-label layer above the prediction raster.
  shiny::observeEvent(list(input$labels_on_top, input$basemap), {
    proxy <- leaflet::leafletProxy("prediction_map") |>
      leaflet::clearGroup("Map labels")

    if (isTRUE(input$labels_on_top)) {
      provider <- label_providers[[input$basemap]]
      proxy |>
        leaflet::addProviderTiles(
          provider,
          group = "Map labels",
          options = leaflet::providerTileOptions(pane = "labels", zIndex = 416)
        )
    }
  }, ignoreInit = TRUE)

  # Replace the raster, legend, and country summaries when model or year changes.
  shiny::observeEvent(list(input$model, input$year), {
    model <- active_model()
    year <- as.integer(input$year)

    tryCatch(
      shiny::withProgress(message = "Loading prediction surface", value = 0.2, {
        layer <- get_raster_layer(model$model_id, year)
        shiny::incProgress(0.35, detail = "Preparing colours")
        domain <- model_domains[[model$model_id]]
        palette <- make_palette(domain, model$is_k13)
        legend <- legend_html(
          model$model_label, model$legend_note, domain, palette, model$is_k13
        )

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

  # Replace observational circles and their size legend when data controls change.
  shiny::observe({
    visible <- data_visible()
    model <- active_model()
    proxy <- leaflet::leafletProxy("prediction_map") |>
      leaflet::clearGroup("Observations") |>
      leaflet::removeControl("observation-size-legend")

    if (!visible || !is.null(observation_state$error)) return(invisible(NULL))

    points <- active_observations()
    all_model_points <- observation_state$data[
      observation_state$data$model == model$model_id &
        is.finite(observation_state$data$Longitude) &
        is.finite(observation_state$data$Latitude),
      , drop = FALSE
    ]
    if (!nrow(points) || !nrow(all_model_points)) return(invisible(NULL))

    domain <- model_domains[[model$model_id]]
    palette <- make_palette(domain, model$is_k13)
    size_scale <- observation_size_scale(all_model_points$Tested)
    popups <- lapply(seq_len(nrow(points)), function(index) {
      observation_popup_html(points[index, , drop = FALSE])
    })

    proxy |>
      leaflet::addCircleMarkers(
        data = points,
        lng = ~Longitude,
        lat = ~Latitude,
        radius = size_scale$radius(points$Tested),
        layerId = paste0("observation-", seq_len(nrow(points))),
        group = "Observations",
        stroke = TRUE,
        color = "#F3EEED",
        weight = 1.25,
        opacity = 1,
        fillColor = palette(points$Prevalence),
        fillOpacity = 0.9,
        popup = popups,
        popupOptions = leaflet::popupOptions(maxWidth = 340, maxHeight = 360),
        options = leaflet::pathOptions(pane = "observation-pane")
      ) |>
      leaflet::addControl(
        html = htmltools::HTML(observation_size_legend_html(size_scale)),
        position = "bottomright",
        layerId = "observation-size-legend"
      )
  })

  # Zoom to a country or return to the standard Africa-centred camera.
  shiny::observeEvent(input$country, {
    proxy <- leaflet::leafletProxy("prediction_map")
    if (identical(input$country, "ALL")) {
      proxy |>
        leaflet::setView(
          lng = INITIAL_MAP_VIEW[["lng"]],
          lat = INITIAL_MAP_VIEW[["lat"]],
          zoom = INITIAL_MAP_VIEW[["zoom"]]
        )
    } else {
      selected <- countries[countries$country_iso3 == input$country, ]
      shiny::req(nrow(selected) == 1L)
      bounds <- sf::st_bbox(selected)
      proxy |>
        leaflet::fitBounds(bounds[["xmin"]], bounds[["ymin"]], bounds[["xmax"]], bounds[["ymax"]])
    }
  }, ignoreInit = TRUE)

  # Keep the floating model/year status synchronized with the active controls.
  output$map_status <- shiny::renderUI({
    model <- active_model()
    shiny::div(
      class = "map-status",
      shiny::span(class = "map-status__dot"),
      shiny::strong(model$model_label),
      shiny::span(" · ", as.integer(input$year))
    )
  })

  # Surface preprocessing problems in the sidebar without crashing the map.
  output$metadata_notice <- shiny::renderUI({
    if (is.null(metadata_state$error)) return(NULL)
    shiny::div(
      class = "metadata-warning",
      shiny::icon("triangle-exclamation"),
      shiny::span(metadata_state$error)
    )
  })

  # Forward client-side screenshot errors through Shiny notifications.
  shiny::observeEvent(input$screenshot_error, {
    shiny::showNotification(input$screenshot_error, type = "error", duration = 8)
  })
}

shiny::shinyApp(ui, server)
