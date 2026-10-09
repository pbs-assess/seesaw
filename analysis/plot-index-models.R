# Shared plots for the index model comparisons in the survey-specific scripts.

# Okabe-Ito (colour-blind safe), dropping black and yellow
highlight_colours <- function(n) {
  unname(grDevices::palette.colors(palette = "Okabe-Ito"))[c(2, 3, 4, 6, 7, 8)][seq_len(n)]
}

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
#' @param connect_stocks Summarise each stock by its maximum amplitude across
#'   windows (for the violins, points, and mean) and connect stocks by a line
#'   across models.
#' @param y_max Upper limit of the amplitude axis; NULL shows the full range.
plot_A_moving_window <- function(seesaw_mw, window, include_all_data = FALSE, connect_stocks = FALSE,
                                 y_max = if (connect_stocks) NULL else 200, n_highlight = 6) {
  if (!include_all_data) {
    dat_mw <- seesaw_mw |>
      filter(!grepl("all data", model)) |>
      mutate(all_data = FALSE)
  } else {
    dat_mw <- seesaw_mw |> mutate(all_data = grepl("all data", model))
  }

  if (connect_stocks) {
    dat_mw <- dat_mw |>
      summarise(A = max(A), .by = c(species, model, all_data))
  }

  dat_mw <- dat_mw |>
    mutate(model = reorder(model, A, FUN = mean))

  light_grey <- "grey85"
  orange <- RColorBrewer::brewer.pal(8, "Oranges")[3]

  if (connect_stocks) {
    y_breaks <- c(0, 2, 10, 25, 50, 100, 200, 400, 600, 800)
    y_lab <- paste0("Maximum biennial amplitude (%)\nacross ", window, "-survey windows")
  } else {
    y_breaks <- c(0, 2, 10, seq(25, 200, 25))
    y_lab <- paste0("Estimated biennial amplitude (%)\nacross ", window, "-survey windows")
  }

  g <- dat_mw |>
    ggplot(aes(model, A)) +
    coord_flip(ylim = if (!is.null(y_max)) c(0, y_max)) +
    ylab(y_lab) +
    ggsidekick::theme_sleek() +
    theme(axis.title.y = element_blank(), panel.grid.major = element_line(colour = "grey90", linewidth = 0.3), panel.grid.minor = element_line(colour = "grey90", linewidth = 0.3))

  if (include_all_data) {
    g <- g + geom_violin(scale = "width", mapping = aes(colour = all_data, fill = all_data)) +
      scale_colour_manual(values = c(light_grey, orange)) +
      scale_fill_manual(values = c(light_grey, orange)) +
      guides(colour = "none", fill = "none")
  } else {
    g <- g + geom_violin(scale = "width", colour = light_grey, fill = light_grey)
  }

  if (connect_stocks) {
    # highlight the top n stocks (by maximum amplitude) in the reference model
    top_species <- dat_mw |>
      filter(model == "IID RF, factor(year)") |>
      slice_max(A, n = n_highlight) |>
      pull(species)
    dat_mw <- dat_mw |>
      mutate(species_hl = ifelse(species %in% top_species, as.character(species), "other"))
    # Okabe-Ito (colour-blind safe), dropping black and yellow
    pal <- setNames(
      c(highlight_colours(n_highlight), "grey50", "black"),
      c(top_species, "other", "Mean")
    )
    n_hl <- length(top_species)
    g <- g +
      geom_line(data = dat_mw, aes(group = species, colour = species_hl), alpha = 0.6) +
      geom_point(data = dat_mw, aes(colour = species_hl), alpha = 0.8) +
      geom_point(stat = "summary", fun = mean, aes(colour = "Mean"), fill = "black", shape = 23, size = 2.5) +
      scale_colour_manual(values = pal, breaks = c(top_species, "Mean"), labels = tools::toTitleCase, name = NULL) +
      guides(colour = guide_legend(override.aes = list(
        shape = c(rep(16, n_hl), 23), linetype = c(rep(1, n_hl), 0), alpha = 1, fill = "black"
      ))) +
      theme(legend.position = "inside", legend.position.inside = c(0.98, 0.02), legend.justification = c(1, 0), legend.background = element_rect(fill = "white", colour = NA))
  } else {
    g <- g +
      geom_point(position = position_jitter(width = 0.1), colour = "grey25", alpha = 0.3) +
      geom_point(stat = "summary", fun = mean, colour = "black", fill = "black", shape = 23, size = 2.5)
  }

  g +
    scale_y_sqrt(limits = c(0, NA), breaks = y_breaks, expand = expansion(mult = c(0, 0.05)))
}

#' Index time series for the stocks with the largest amplitude in a reference
#' model: rows are stocks, columns are models.
#' @param n_top Number of stocks to show (ranked by maximum amplitude in the
#'   reference model).
#' @param models Optional subset of models to show as columns (in order); the
#'   reference model is always shown first.
#' @param ref_model Reference model used to choose the stocks.
plot_top_stock_indexes <- function(out, seesaw_mw, n_top = 6, ref_model = "IID RF, factor(year)",
                                   include_all_data = FALSE, models = NULL, lu = NULL,
                                   .ylab = "Centered biomass index") {
  if (!include_all_data) {
    out <- filter(out, !grepl("all data", model))
    seesaw_mw <- filter(seesaw_mw, !grepl("all data", model))
  }

  top_species <- seesaw_mw |>
    filter(model == ref_model) |>
    summarise(A = max(A), .by = species) |>
    slice_max(A, n = n_top) |>
    pull(species)

  # columns ordered by mean amplitude across all stocks, reference first
  model_order <- seesaw_mw |>
    summarise(mean_A = mean(A), .by = model) |>
    arrange(mean_A) |>
    pull(model)
  model_order <- c(ref_model, setdiff(model_order, ref_model))
  if (!is.null(models)) model_order <- c(ref_model, setdiff(models, ref_model))

  dat <- out |>
    filter(species %in% top_species, model %in% model_order) |>
    group_by(species, model) |>
    mutate(geomean = exp(mean(log(est))), est = est / geomean, lwr = lwr / geomean, upr = upr / geomean) |>
    ungroup() |>
    mutate(
      species = factor(species, levels = top_species),
      model = factor(model, levels = model_order)
    )

  # per stock, the year parity (even/odd) with the higher mean index in the
  # reference model; if `lu` (year, phase) is supplied, use its phase instead
  if (is.null(lu)) {
    dat <- mutate(dat, even = year %% 2 == 0)
  } else {
    dat <- left_join(dat, distinct(lu, year, phase), by = "year") |>
      mutate(even = phase > 0)
  }
  upper_parity <- dat |>
    filter(model == ref_model) |>
    summarise(m = mean(log(est)), .by = c(species, even)) |>
    slice_max(m, n = 1, by = species) |>
    select(species, upper_even = even)
  dat <- dat |>
    left_join(upper_parity, by = "species") |>
    mutate(upper = even == upper_even)

  ggplot(dat, aes(year, est, ymin = lwr, ymax = upr)) +
    geom_ribbon(fill = "grey90") +
    geom_linerange(aes(colour = species)) +
    geom_point(aes(colour = species, shape = upper)) +
    scale_colour_manual(values = setNames(highlight_colours(n_top), top_species)) +
    # open = the parity with the higher mean index in the reference model
    scale_shape_manual(values = c("TRUE" = 1, "FALSE" = 16)) +
    guides(colour = "none", shape = "none") +
    facet_grid(species ~ model, scales = "free_y", labeller = labeller(species = \(x) label_wrap_gen(14)(tools::toTitleCase(x)))) +
    scale_y_log10() +
    scale_x_continuous(breaks = seq(2005, 2025, 5)) +
    ylab(.ylab) +
    ggsidekick::theme_sleek() +
    theme(axis.title.x = element_blank(), panel.grid.major = element_line(colour = "grey90", linewidth = 0.3), panel.grid.minor = element_blank())
}
