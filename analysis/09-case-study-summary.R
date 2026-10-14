# Combined case-study figure: moving-window seesaw amplitude by model for the
# BC synoptic trawl, HBLL OUT, and transboundary AFSC GOA + PBS trawl surveys
# (one row each), with sparklines of the most seesawing stocks beside each.
#
# Uses the cached index fits from 03-fit-syn-models.R,
# 04-fit-afsc-pbs-models.R, and 06-fit-hbll-models.R; run those first.

library(dplyr)
library(ggplot2)
library(patchwork)

source(here::here("analysis/metric-functions.R"))
source(here::here("analysis/plot-index-models.R"))

seesaw_window <- 10L
n_top <- 5L
ref_model <- "IID RF, factor(year)"
spark_models <- c(ref_model, "IID RF, RW year")

read_indexes <- function(file) {
  readRDS(here::here("data-generated", file)) |>
    filter(!grepl("depth", model), !grepl("all data", model), model != "Spatial only, RW year", model != "AR1 RF") |>
    mutate(species = replace(species, species == "north pacific spiny dogfish", "pacific spiny dogfish"))
}

# HBLL N and S alternate, but the parity flips after the 2013 gap, so phase
# is coded by region (as in 06-fit-hbll-models.R); the trawl case studies use
# calendar-year parity
hbll_phase <- readRDS(here::here("data-raw/hbll-out-sets.rds")) |>
  distinct(year, survey_abbrev) |>
  mutate(phase = ifelse(survey_abbrev == "HBLL OUT N", 0.5, -0.5)) |>
  select(year, phase)
stopifnot(!any(duplicated(hbll_phase$year)))

case_studies <- list(
  "BC synoptic trawl" = list(out = read_indexes("bc-indexes4.rds"), phase = NULL),
  "HBLL outside" = list(out = read_indexes("hbll-indexes.rds"), phase = hbll_phase),
  "AFSC GOA + PBS synoptic trawl" = list(out = read_indexes("transboundary-afsc-pbs-indexes.rds"), phase = NULL)
)

calc_mw <- function(out, phase) {
  if (!is.null(phase)) out <- left_join(out, phase, by = "year")
  out |>
    arrange(species, model, year) |>
    group_by(species, model) |>
    group_modify(\(.x, .y) {
      moving_window(.x$est, year = .x$year, window = seesaw_window, phase = if ("phase" %in% names(.x)) .x$phase)
    }) |>
    ungroup()
}
seesaw_mw <- purrr::map(case_studies, \(cs) calc_mw(cs$out, cs$phase))

# Common model order (mean of the per-case-study mean max A) so the rows are
# comparable
max_A <- purrr::map_dfr(seesaw_mw, \(mw) summarise(mw, A = max(A), .by = c(species, model)), .id = "case_study")
model_order <- max_A |>
  summarise(A = mean(A), .by = c(case_study, model)) |>
  summarise(A = mean(A), .by = model) |>
  arrange(A) |>
  pull(model)

# Same ranking as plot_A_moving_window() so the colours match
top_stocks <- function(mw) {
  mw |>
    filter(model == ref_model) |>
    summarise(A = max(A), .by = species) |>
    slice_max(A, n = n_top) |>
    pull(species)
}

plot_A_panel <- function(mw, title, show_x_title) {
  plot_A_moving_window(mw, seesaw_window, connect_stocks = TRUE, n_highlight = n_top) +
    scale_x_discrete(limits = model_order) +
    ggtitle(title) +
    theme(
      legend.position = "none",
      plot.title = element_text(size = rel(1)),
      axis.title.x = if (show_x_title) element_text() else element_blank()
    )
}

# Tufte-style sparklines: centred log index for the top stocks under the
# offending and fixing models, same y scale within each stock, with the
# stock's maximum moving-window A at the end of each line. Model headers are
# drawn above the top panels (in the space of the A panel's title) rather than
# as strips, which patchwork would align by pushing down the A panel.
plot_sparklines <- function(out, mw, show_header) {
  top <- top_stocks(mw)
  dat <- out |>
    filter(species %in% top, model %in% spark_models) |>
    mutate(log_est = log(est) - mean(log(est)), .by = c(species, model)) |>
    mutate(
      species = factor(species, levels = top),
      model = factor(model, levels = spark_models)
    )
  labs_A <- mw |>
    filter(species %in% top, model %in% spark_models) |>
    summarise(A = max(A), .by = c(species, model)) |>
    left_join(
      dat |> slice_max(year, n = 1, by = c(species, model)) |> select(species, model, year, log_est),
      by = c("species", "model")
    ) |>
    mutate(
      species = factor(species, levels = top),
      model = factor(model, levels = spark_models),
      label = paste0(round(A), "%")
    )
  headers <- tibble::tibble(
    species = factor(top[1], levels = top),
    model = factor(spark_models, levels = spark_models),
    year = mean(range(dat$year)),
    label = spark_models
  )

  ggplot(dat, aes(year, log_est, colour = species)) +
    geom_line() +
    geom_text(aes(label = label), data = labs_A, hjust = -0.15, size = 2.5, colour = "grey40") +
    geom_text(
      aes(label = label), data = if (show_header) headers else headers[0, ],
      y = Inf, vjust = -0.6, size = 2.8, colour = "grey20"
    ) +
    facet_grid(
      species ~ model,
      scales = "free_y", switch = "y",
      labeller = labeller(species = \(x) tools::toTitleCase(x))
    ) +
    scale_colour_manual(values = setNames(highlight_colours(n_top), top), guide = "none") +
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.22))) +
    coord_cartesian(clip = "off") +
    theme_void() +
    theme(
      strip.text.y.left = element_text(angle = 0, hjust = 1, size = rel(0.85)),
      strip.text.x = element_blank(),
      panel.spacing.y = unit(2, "pt"),
      panel.spacing.x = unit(6, "pt"),
      plot.margin = margin(5.5, 5.5, 5.5, 5.5)
    )
}

# One 3 x 2 grid (not nested rows) so the columns align across case studies
panels <- purrr::imap(names(case_studies), \(cs, i) {
  list(
    plot_A_panel(seesaw_mw[[cs]], cs, show_x_title = i == length(case_studies)),
    plot_sparklines(case_studies[[cs]]$out, seesaw_mw[[cs]], show_header = i == 1L)
  )
}) |>
  purrr::list_flatten()

wrap_plots(panels, ncol = 2, widths = c(1.5, 1))
ggsave(here::here("figs/case-study-A-moving-window.pdf"), width = 7.5, height = 8)
