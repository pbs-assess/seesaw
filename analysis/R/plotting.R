library(dplyr, warn.conflicts = FALSE)
library(tidyr)
library(purrr)
library(ggplot2)
theme_set(ggsidekick::theme_sleek())
dir.create(here::here("figs"), showWarnings = FALSE)

# Okabe-Ito (colour-blind safe), dropping black and yellow
highlight_colours <- function(n) {
  unname(palette.colors(palette = "Okabe-Ito"))[c(2, 3, 4, 6, 7, 8)][seq_len(n)]
}

# Case studies -------------------------------------------------------------------

# max_A: one row per species x model with the maximum A across windows
plot_A_moving_window <- function(max_A, window, top_species) {
  dat <- mutate(
    max_A,
    model = reorder(model, A, FUN = mean),
    species_hl = ifelse(species %in% top_species, as.character(species), "other")
  )
  n_hl <- length(top_species)
  pal <- setNames(
    c(highlight_colours(n_hl), "grey50", "black"),
    c(top_species, "other", "Mean")
  )
  ggplot(dat, aes(model, A)) +
    coord_flip() +
    geom_violin(scale = "width", colour = "grey85", fill = "grey85") +
    geom_line(aes(group = species, colour = species_hl), alpha = 0.6) +
    geom_point(aes(colour = species_hl), alpha = 0.8) +
    geom_point(stat = "summary", fun = mean, aes(colour = "Mean"), fill = "black", shape = 23, size = 2.5) +
    scale_colour_manual(values = pal, breaks = c(top_species, "Mean"), labels = tools::toTitleCase, name = NULL) +
    guides(colour = guide_legend(override.aes = list(
      shape = c(rep(16, n_hl), 23), linetype = c(rep(1, n_hl), 0), alpha = 1, fill = "black"
    ))) +
    scale_y_sqrt(
      limits = c(0, NA), breaks = c(0, 2, 10, 25, 50, 100, 200, 400, 600, 800),
      expand = expansion(mult = c(0, 0.05))
    ) +
    ylab(paste0("Maximum biennial amplitude (%)\nacross ", window, "-year windows")) +
    theme(
      axis.title.y = element_blank(),
      panel.grid.major = element_line(colour = "grey90", linewidth = 0.3),
      panel.grid.minor = element_line(colour = "grey90", linewidth = 0.3),
      legend.position = "inside", legend.position.inside = c(0.98, 0.02), legend.justification = c(1, 0),
      legend.background = element_rect(fill = "white", colour = NA)
    )
}

# Self-tests -----------------------------------------------------------------------

summary_metric_labels <- c(
  rmse_ratio = "RMSE / full-data RMSE",
  coverage_50 = "50% CI coverage",
  seesaw_abs_error = "Seesaw error (%)",
  variability = "Year-to-year SD\n(est. / true)"
)
summary_nominal <- c(rmse_ratio = 1, coverage_50 = 0.5, seesaw_abs_error = 0, variability = 1)
# Summarised across replicates by the median (the rest by the mean)
summary_median_metrics <- c("seesaw_abs_error", "variability")

# Violins of all replicates, one point per species, crossbar at the median
# across species. Each panel zooms to the widest 10-90% range of the
# replicates, widened to keep the species points and reference lines in view;
# oob_keep so the violins still use all replicates.
#
# y: factor in top-to-bottom order. seesaw_max: species beyond it are pinned to
# it and drawn as triangles. rmse_pad: extra room in the RMSE panel (fraction
# of its range), e.g. for an inset legend.
plot_metric_summary <- function(rep_metrics, metrics, y, floor, row = NULL, connect = FALSE,
                                seesaw_max = NULL, rmse_pad = 0) {
  labels <- summary_metric_labels[metrics]
  groups <- c(y, row)
  rep_long <- rep_metrics |>
    select(species, rep, all_of(c(groups, metrics))) |>
    pivot_longer(all_of(metrics), names_to = "metric") |>
    drop_na(value)
  species_long <- rep_long |>
    summarise(
      value = if (first(metric) %in% summary_median_metrics) median(value) else mean(value),
      .by = all_of(c("species", groups, "metric"))
    ) |>
    mutate(beyond = metric == "seesaw_abs_error" & value > (seesaw_max %||% Inf)) |>
    mutate(value = if_else(beyond, seesaw_max %||% Inf, value)) |>
    # geom_path joins points in row order
    arrange(species, metric, .data[[y]])
  nominal <- tibble(metric = metrics, value = summary_nominal[metrics])
  floors <- tibble(metric = intersect("seesaw_abs_error", metrics), value = floor)

  limits <- bind_rows(
    summarise(rep_long, lo = quantile(value, 0.1), hi = quantile(value, 0.9), .by = all_of(c("metric", groups))),
    transmute(species_long, metric, lo = value, hi = value),
    transmute(nominal, metric, lo = value, hi = value),
    transmute(floors, metric, lo = value, hi = value)
  ) |>
    summarise(lo = min(lo), hi = max(hi), .by = metric)
  limits <- limits[match(metrics, limits$metric), ]
  if (!is.null(seesaw_max)) limits$hi[limits$metric == "seesaw_abs_error"] <- seesaw_max
  i_rmse <- limits$metric == "rmse_ratio"
  limits$hi[i_rmse] <- limits$hi[i_rmse] + rmse_pad * (limits$hi[i_rmse] - limits$lo[i_rmse])
  # The SD ratio is on a log scale (too smooth and too wiggly by the same
  # factor are equally far from 1)
  x_scales <- pmap(limits, \(metric, lo, hi) {
    if (metric == "variability") {
      scale_x_log10(limits = c(lo, hi), oob = scales::oob_keep)
    } else {
      scale_x_continuous(limits = c(lo, hi), oob = scales::oob_keep)
    }
  })

  as_metric <- \(d) mutate(d, metric = factor(metric, levels = metrics, labels = labels))
  ggplot(as_metric(species_long), aes(value, .data[[y]])) +
    geom_violin(data = as_metric(rep_long), fill = "grey92", colour = NA, scale = "width") +
    geom_vline(aes(xintercept = value), data = as_metric(nominal), linetype = 2, colour = "grey50") +
    geom_vline(aes(xintercept = value), data = as_metric(floors), linetype = 3, colour = "grey50") +
    (if (connect) geom_path(aes(colour = species, group = species))) +
    geom_point(
      aes(colour = species, shape = beyond),
      position = if (connect) "identity" else position_jitter(height = 0.15, seed = 1)
    ) +
    stat_summary(fun = median, geom = "crossbar", fun.min = median, fun.max = median, width = 0.6) +
    scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17), guide = "none") +
    (if (is.null(row)) {
      facet_wrap(~metric, nrow = 1, scales = "free_x")
    } else {
      facet_grid(reformulate("metric", row), scales = "free_x")
    }) +
    ggh4x::facetted_pos_scales(x = x_scales) +
    # Paired has 12 colours (Dark2 has only 8)
    scale_colour_brewer(palette = "Paired", labels = tools::toTitleCase) +
    scale_y_discrete(limits = rev) +
    labs(x = NULL, y = NULL, colour = NULL) +
    theme(
      legend.position = "bottom",
      panel.grid.major = element_line(colour = "grey90", linewidth = rel(0.5)),
      panel.spacing.x = unit(1, "lines")
    ) +
    guides(colour = guide_legend(nrow = 2))
}
