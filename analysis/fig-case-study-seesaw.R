# Norway uses the fits to all the data (phase as if the alternation had
# continued), not the alternating subset.

library(dplyr)
library(patchwork)
source(here::here("analysis/R/seesaw-metrics.R"))
source(here::here("analysis/R/plotting.R"))

seesaw_window <- 10L # calendar years
n_top <- 5L
ref_model <- "IID RF, factor(year)"
spark_models <- c(ref_model, "IID RF, RW year")

read_case_study <- \(name) readRDS(here::here("data-generated", paste0("case-study-", name, ".rds")))$indexes

# Order, titles and panel tags match the maps in fig-case-study-maps.R
case_studies <- list(
  "Gulf of Alaska (GOA) + BC trawl" = read_case_study("goa-bc-trawl"),
  "Norwegian gillnet–fyke net" = read_case_study("norway") |> filter(design == "all"),
  "BC trawl" = read_case_study("bc-trawl"),
  "BC longline" = read_case_study("bc-longline")
)

max_A <- map(case_studies, \(out) {
  out |>
    arrange(species, model, year) |>
    summarise(A = max(moving_window(est, year, phase, window = seesaw_window)$A), .by = c(species, model))
})

top_stocks <- map(max_A, \(d) {
  d |>
    filter(model == ref_model) |>
    slice_max(A, n = n_top) |>
    pull(species)
})

# One model order so the rows are comparable
model_order <- bind_rows(max_A, .id = "case_study") |>
  summarise(A = mean(A), .by = c(case_study, model)) |>
  summarise(A = mean(A), .by = model) |>
  arrange(A) |>
  pull(model)

plot_A_panel <- function(d, top, title, tag, show_x_title) {
  g <- plot_A_moving_window(d, seesaw_window, top)
  # The default breaks crowd together on the sqrt scale when A runs into the
  # thousands (Norway)
  if (max(d$A) > 1000) {
    g <- suppressMessages(g + scale_y_sqrt(
      limits = c(0, NA), breaks = c(0, 10, 50, 200, 500, 1000, 2000),
      expand = expansion(mult = c(0, 0.05))
    ))
  }
  g +
    scale_x_discrete(limits = model_order) +
    labs(title = title, tag = tag) +
    theme(
      legend.position = "none",
      plot.title = element_text(size = rel(1)),
      axis.title.x = if (show_x_title) element_text() else element_blank()
    )
}

# Model headers are drawn above the top panels rather than as strips, which
# patchwork would align by pushing down the A panel
plot_sparklines <- function(out, d_A, top, show_header) {
  dat <- out |>
    filter(species %in% top, model %in% spark_models) |>
    mutate(log_est = log(est) - mean(log(est)), .by = c(species, model)) |>
    mutate(
      species = factor(species, levels = top),
      model = factor(model, levels = spark_models)
    )
  labs_A <- d_A |>
    filter(species %in% top, model %in% spark_models) |>
    left_join(
      dat |> slice_max(year, n = 1, by = c(species, model)) |> select(species, model, year, log_est) |>
        mutate(across(c(species, model), as.character)),
      by = c("species", "model")
    ) |>
    mutate(
      species = factor(species, levels = top),
      model = factor(model, levels = spark_models),
      label = paste0(round(A), "%")
    )
  headers <- tibble(
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
    scale_colour_manual(values = set_names(highlight_colours(n_top), top), guide = "none") +
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

# One 4 x 2 grid (not nested rows) so the columns align across case studies
panels <- imap(names(case_studies), \(cs, i) {
  list(
    plot_A_panel(max_A[[cs]], top_stocks[[cs]], cs, letters[i], show_x_title = i == length(case_studies)),
    plot_sparklines(case_studies[[cs]], max_A[[cs]], top_stocks[[cs]], show_header = i == 1L)
  )
}) |>
  list_flatten()

g <- wrap_plots(panels, ncol = 2, widths = c(1.5, 1))
ggsave(here::here("figs/case-study-seesaw.pdf"), g, width = 7.5, height = 10.5, device = cairo_pdf)
