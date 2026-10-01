# Interactive 3D waterfall plot for the shrimp SD-DRM

library(dplyr)
library(purrr)
library(plotly)

source("model_functions.R")

# Model inputs

set.seed(123)
exposure <- read.csv("raw_shrimp_mc2d_simulation_data.csv") |> arrange(inner_id)
n_inner <- nrow(exposure)
n_outer <- 500
dose <- rep(exposure$D_raw, times = n_outer)

alpha <- 0.60
beta <- 1.31e6
tfs <- 1
MIC <- 0.13

pd_mean <- c(Emax_per_day = 81.587790, EC50_ug_mL = 0.166877)
pd_cov <- matrix(c(32.631766, 0.124149033, 0.124149033, 0.002101326),
                 nrow = 2, byrow = TRUE)

draw_pd <- function(n) {
  draws <- MASS::mvrnorm(n, pd_mean, pd_cov)
  invalid <- draws[, 1] <= 0 | draws[, 2] <= 0
  
  while (any(invalid)) {
    draws[invalid, ] <- MASS::mvrnorm(sum(invalid), pd_mean, pd_cov)
    invalid <- draws[, 1] <= 0 | draws[, 2] <= 0
  }
  
  data.frame(Emax = draws[, 1], EC50 = draws[, 2])
}

set.seed(123)
pd <- draw_pd(n_outer)

ab_lookup <- make_ab_lookup(
  alpha = alpha, beta = beta, Emax = pd$Emax, EC50 = pd$EC50,
  C_max = 0.10 * MIC, tfs = tfs, n_grid = 250, seed = 0, nsim = 1000
)

# Risk summaries

run_scenario <- function(fr, pMIC) {
  result <- calc_sd_drm(
    dose = dose, fr = fr, C = pMIC / 100 * MIC,
    alpha = alpha, beta = beta,
    Emax = rep(pd$Emax, each = n_inner),
    EC50 = rep(pd$EC50, each = n_inner),
    tfs = tfs, ab_lookup = ab_lookup
  )
  
  list(
    risk = matrix(result$risk_total, nrow = n_inner),
    delta = matrix(result$delta_p_treat, nrow = n_inner)
  )
}

summarise_risk <- function(vary_by, values, fixed_value) {
  map_dfr(values, function(value) {
    fr <- if (vary_by == "fr") value else fixed_value
    pMIC <- if (vary_by == "pMIC") value else fixed_value
    scenario <- run_scenario(fr, pMIC)
    expected_risk <- colMeans(scenario$risk)
    limits <- quantile(expected_risk, c(0.025, 0.50, 0.975), names = FALSE)
    prop_more <- mean(colMeans(scenario$delta) > 0)
    
    tibble(
      value = value,
      lower95 = limits[1], median = limits[2], upper95 = limits[3],
      prop_more = prop_more,
      outcome = ifelse(prop_more > 0.5,
                       "More likely treatable", "Less likely treatable")
    )
  })
}

find_cutoff <- function(data, vary_by, fixed_value, tolerance) {
  data <- arrange(data, value)
  transition <- which(head(data$prop_more > 0.5, -1) &
                        !tail(data$prop_more > 0.5, -1))[1]
  if (is.na(transition)) return(NA_real_)
  
  uniroot(
    function(value) summarise_risk(vary_by, value, fixed_value)$prop_more - 0.5,
    interval = data$value[c(transition, transition + 1)],
    tol = tolerance
  )$root
}

make_view_data <- function(vary_by, slices, grid) {
  map_dfr(slices, function(slice) {
    summarise_risk(vary_by, grid, slice) |>
      mutate(
        fr = if (vary_by == "fr") value else slice,
        pMIC = if (vary_by == "pMIC") value else slice,
        log10_lower95 = log10(lower95),
        log10_median = log10(median),
        log10_upper95 = log10(upper95),
        hover = paste0(
          "Resistant fraction: ", sprintf("%.3f", fr),
          "<br>pMIC: ", sprintf("%.2f", pMIC), "% MIC",
          "<br>Median log10 risk: ", sprintf("%.3f", log10_median),
          "<br>95% interval: ", sprintf("%.3f", log10_lower95),
          " to ", sprintf("%.3f", log10_upper95),
          "<br>Treatment outcome: ", outcome
        )
      )
  })
}

make_cutoffs <- function(data, vary_by, slices, tolerance) {
  x_name <- if (vary_by == "pMIC") "pMIC" else "fr"
  slice_name <- if (vary_by == "pMIC") "fr" else "pMIC"
  
  map_dfr(slices, function(slice) {
    scenario <- data[abs(data[[slice_name]] - slice) < 1e-10, ]
    cutoff <- find_cutoff(scenario, vary_by, slice, tolerance)
    if (!is.finite(cutoff)) return(tibble())
    
    cutoff_risk <- summarise_risk(vary_by, cutoff, slice)$median
    values <- setNames(c(cutoff, slice), c(x_name, slice_name))
    
    tibble(
      fr = values["fr"], pMIC = values["pMIC"],
      log10_median = log10(cutoff_risk),
      hover = paste0(
        "Treatment cut-point",
        "<br>Resistant fraction: ", sprintf("%.3f", values["fr"]),
        "<br>pMIC: ", sprintf("%.3f", values["pMIC"]), "% MIC"
      )
    )
  })
}

fr_slices <- seq(0, 0.5, by = 0.05)
pMIC_slices <- 0:10

view_A <- make_view_data("pMIC", fr_slices, seq(0, 10, by = 0.1))
view_B <- make_view_data("fr", pMIC_slices, seq(0, 0.5, by = 0.025))
cutoffs_A <- make_cutoffs(view_A, "pMIC", fr_slices, 0.001)
cutoffs_B <- make_cutoffs(view_B, "fr", pMIC_slices, 0.0001)

# Interactive plot

make_mesh <- function(data, x_name, slice_name) {
  data <- data[order(data[[x_name]]), ]
  n <- nrow(data)
  interval <- 0:(n - 2)
  
  list(
    x = rep(data[[x_name]], 2),
    y = rep(data[[slice_name]], 2),
    z = c(data$log10_lower95, data$log10_upper95),
    i = c(interval, interval + 1),
    j = c(interval + 1, n + interval + 1),
    k = c(n + interval, n + interval)
  )
}

split_lines <- function(data, cutoff, x_name) {
  if (nrow(cutoff) > 0) {
    cutoff <- select(cutoff, fr, pMIC, log10_median, hover)
    more <- bind_rows(filter(data, outcome == "More likely treatable"), cutoff)
    less <- bind_rows(cutoff, filter(data, outcome == "Less likely treatable"))
  } else {
    more <- filter(data, outcome == "More likely treatable")
    less <- filter(data, outcome == "Less likely treatable")
  }
  
  list(
    more = more |> distinct(fr, pMIC, .keep_all = TRUE) |>
      arrange(.data[[x_name]]),
    less = less |> distinct(fr, pMIC, .keep_all = TRUE) |>
      arrange(.data[[x_name]])
  )
}

add_view <- function(plot, data, cutoffs, x_name, slice_name,
                     slices, view_id, visible) {
  trace_ids <- character()
  more_legend <- less_legend <- FALSE
  
  for (slice in slices) {
    scenario <- data[abs(data[[slice_name]] - slice) < 1e-10, ]
    scenario <- scenario[order(scenario[[x_name]]), ]
    cutoff <- cutoffs[abs(cutoffs[[slice_name]] - slice) < 1e-10, ]
    lines <- split_lines(scenario, cutoff, x_name)
    mesh <- make_mesh(scenario, x_name, slice_name)
    
    plot <- add_trace(
      plot, x = mesh$x, y = mesh$y, z = mesh$z,
      i = mesh$i, j = mesh$j, k = mesh$k, type = "mesh3d",
      color = "#6BAED6", facecolor = rep("#6BAED6", length(mesh$i)),
      opacity = 0.40, name = "95% uncertainty interval",
      legendrank = 4, showlegend = slice == slices[1],
      visible = visible, hoverinfo = "skip"
    )
    trace_ids <- c(trace_ids, view_id)
    
    if (nrow(lines$more) > 1) {
      plot <- add_trace(
        plot, x = lines$more[[x_name]], y = lines$more[[slice_name]],
        z = lines$more$log10_median, text = lines$more$hover,
        type = "scatter3d", mode = "lines", hoverinfo = "text",
        line = list(color = "black", width = 5),
        name = "More likely treatable", legendrank = 2,
        showlegend = !more_legend, visible = visible
      )
      trace_ids <- c(trace_ids, view_id)
      more_legend <- TRUE
    }
    
    if (nrow(lines$less) > 1) {
      plot <- add_trace(
        plot, x = lines$less[[x_name]], y = lines$less[[slice_name]],
        z = lines$less$log10_median, text = lines$less$hover,
        type = "scatter3d", mode = "lines", hoverinfo = "text",
        line = list(color = "red", width = 5),
        name = "Less likely treatable", legendrank = 1,
        showlegend = !less_legend, visible = visible
      )
      trace_ids <- c(trace_ids, view_id)
      less_legend <- TRUE
    }
  }
  
  if (nrow(cutoffs) > 0) {
    plot <- add_trace(
      plot, x = cutoffs[[x_name]], y = cutoffs[[slice_name]],
      z = cutoffs$log10_median, text = cutoffs$hover,
      type = "scatter3d", mode = "markers", hoverinfo = "text",
      marker = list(
        symbol = "circle", color = "white", size = 6,
        line = list(color = "black", width = 2)
      ),
      name = "Treatment cut-point", legendrank = 3,
      showlegend = TRUE, visible = visible
    )
    trace_ids <- c(trace_ids, view_id)
  }
  
  list(plot = plot, trace_ids = trace_ids)
}

plot_3d <- plot_ly()
result_A <- add_view(plot_3d, view_A, cutoffs_A, "pMIC", "fr",
                     fr_slices, "A", TRUE)
result_B <- add_view(result_A$plot, view_B, cutoffs_B, "fr", "pMIC",
                     pMIC_slices, "B", FALSE)
plot_3d <- result_B$plot
trace_ids <- c(result_A$trace_ids, result_B$trace_ids)

z_limits <- range(c(view_A$log10_lower95, view_A$log10_upper95,
                    view_B$log10_lower95, view_B$log10_upper95), finite = TRUE)
z_ticks <- seq(floor(z_limits[1] * 5) / 5,
               ceiling(z_limits[2] * 5) / 5, by = 0.2)
z_limits <- range(z_ticks) + c(-0.05, 0.05)

make_scene <- function(x_title, x_ticks, y_title, y_ticks, aspect, eye) {
  list(
    xaxis = list(title = x_title, tickmode = "array", tickvals = x_ticks[x_ticks > 0],
                 range = range(x_ticks)),
    yaxis = list(title = y_title, tickmode = "array", tickvals = y_ticks,
                 range = range(y_ticks)),
    zaxis = list(title = "log<sub>10</sub>(P<sub>illness</sub>)",
                 tickmode = "array", tickvals = z_ticks,
                 ticktext = sprintf("%.1f", z_ticks), range = z_limits),
    aspectmode = "manual", aspectratio = aspect,
    camera = list(
      eye = eye,
      center = list(x = -0.18, y = 0, z = 0)
    )
  )
}

scene_A <- make_scene(
  "Residual doxycycline concentration (pMIC)", pMIC_slices,
  "Resistant fraction (fᵣ)", fr_slices,
  list(x = 1.5, y = 0.9, z = 1),
  list(x = 1.5, y = -1.8, z = 1.1)
)

scene_B <- make_scene(
  "Resistant fraction (fᵣ)", fr_slices,
  "Residual doxycycline concentration (pMIC)", pMIC_slices,
  list(x = 1.5, y = 1.1, z = 1),
  list(x = 1.9, y = -1.4, z = 1.1)
)

view_button <- function(label, view_id, scene) {
  list(
    label = label, method = "update",
    args = list(
      list(visible = trace_ids == view_id),
      list(title = list(text = label, x = 0.5), scene = scene)
    )
  )
}

plot_3d <- plot_3d |>
  layout(
    title = list(text = "A. Fixed fᵣ and varying pMIC", x = 0.5),
    scene = scene_A,
    updatemenus = list(list(
      type = "buttons", direction = "right",
      x = 0.5, xanchor = "center", y = 1.16,
      buttons = list(
        view_button("A. Fixed fᵣ and varying pMIC", "A", scene_A),
        view_button("B. Fixed pMIC and varying fᵣ", "B", scene_B)
      )
    )),
    legend = list(
      orientation = "h", x = 0.5, xanchor = "center", y = 1.06,
      itemsizing = "constant", itemclick = FALSE, itemdoubleclick = FALSE
    ),
    margin = list(l = 0, r = 0, b = 0, t = 100)
  ) |>
  config(displayModeBar = FALSE, responsive = TRUE) |>
  htmlwidgets::onRender(
    "function(el) { document.title = 'SD-DRM shrimp model'; }"
  )


# Save the interactive plot

htmlwidgets::saveWidget(
  plot_3d, "images/Interactive 3D plot.html", selfcontained = TRUE
)

unlink(
  "images/Interactive 3D plot_files",
  recursive = TRUE
)

