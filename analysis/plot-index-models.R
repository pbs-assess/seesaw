# Shared plots for the index model comparisons in the survey-specific scripts.

#' Standardized index time series, one panel per model x species, with model
#' rows ordered by mean moving-window amplitude.
#' @param out Output of `fit_index_models()` with a `species` column.
#' @param lu Year lookup with `year` and the `colour` column.
#' @param seesaw_mw Moving-window amplitudes with `species`, `model`, `A`.
#' @param colour Name of the column in `lu` to colour points by.
plot_indexes <- function(out, lu, seesaw_mw, colour, colour_lab = "Survey", .ylab = "Index") {
  seesaw_summary <- seesaw_mw |>
    summarise(mean_A = mean(A), .by = c(species, model))

  out |>
    left_join(lu, by = "year") |>
    left_join(seesaw_summary, by = c("species", "model")) |>
    group_by(species, model) |>
    mutate(geomean = exp(mean(log(est))), est = est / geomean, lwr = lwr / geomean, upr = upr / geomean) |>
    ggplot(aes(year, log(est), ymin = log(lwr), ymax = log(upr))) +
    geom_ribbon(fill = "grey90") +
    geom_linerange(aes(colour = .data[[colour]])) +
    geom_point(aes(colour = .data[[colour]])) +
    scale_colour_brewer(palette = "Dark2") +
    facet_grid(forcats::fct_reorder(model, mean_A) ~ species) +
    ylab(.ylab) +
    xlab("Year") +
    labs(colour = colour_lab) +
    ggsidekick::theme_sleek()
}

#' Distribution of moving-window biennial amplitudes by model.
#' @param include_all_data Show the "all data" models in a second colour;
#'   if FALSE they are dropped.
plot_A_moving_window <- function(seesaw_mw, window, include_all_data = FALSE) {
  if (!include_all_data) {
    dat_mw <- seesaw_mw |>
      filter(!grepl("all data", model)) |>
      mutate(all_data = FALSE)
  } else {
    dat_mw <- seesaw_mw |> mutate(all_data = grepl("all data", model))
  }

  dat_mw <- dat_mw |>
    mutate(model = reorder(model, A, FUN = mean))

  blue <- RColorBrewer::brewer.pal(8, "Blues")[3]
  orange <- RColorBrewer::brewer.pal(8, "Oranges")[3]

  g <- dat_mw |>
    ggplot(aes(model, A)) +
    coord_flip(ylim = c(0, 200)) +
    ylab(paste0("Estimated biennial amplitude (%)\nacross ", window, "-survey windows")) +
    ggsidekick::theme_sleek() +
    theme(axis.title.y = element_blank(), panel.grid.major = element_line(colour = "grey90", linewidth = 0.3), panel.grid.minor = element_line(colour = "grey90", linewidth = 0.3))

  if (include_all_data) {
    g <- g + geom_violin(scale = "width", mapping = aes(colour = all_data, fill = all_data)) +
      scale_colour_manual(values = c(blue, orange)) +
      scale_fill_manual(values = c(blue, orange)) +
      guides(colour = "none", fill = "none")
  } else {
    g <- g + geom_violin(scale = "width", colour = blue, fill = blue)
  }

  g +
    geom_point(position = position_jitter(width = 0.1), colour = "grey25", alpha = 0.3) +
    geom_point(stat = "summary", fun = mean, colour = "black") +
    scale_y_sqrt(limits = c(0, NA), expand = expansion(mult = c(0, 0.05)))
}
