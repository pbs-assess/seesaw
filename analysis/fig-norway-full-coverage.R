library(dplyr)
library(patchwork)
source(here::here("analysis/R/seesaw-metrics.R"))
source(here::here("analysis/R/plotting.R"))

norway_models <- c("IID RF, factor(year)", "IID RF, RW year")
norway_window <- 10L # calendar years
norway_data <- c(alternating = "Alternating only", all = "With full 2024-2026")
norway_colours <- set_names(c("#D55E00", "#0072B2"), norway_data)
inside_legend_x <- 0.71

norway <- readRDS(here::here("data-generated/case-study-norway.rds"))$indexes |>
  filter(model %in% norway_models) |>
  mutate(data = norway_data[design]) |>
  # only species with every model fit under both designs, so both panels show
  # the same species
  filter(n_distinct(data, model) == length(norway_data) * length(norway_models), .by = species) |>
  mutate(
    model = factor(model, levels = norway_models),
    data = factor(data, levels = norway_data)
  )
norway_A <- norway |>
  arrange(year) |>
  summarise(
    A_mean = mean(moving_window(est, year, phase, window = norway_window)$A),
    .by = c(species, model, data)
  )
norway_order <- norway_A |>
  filter(model == norway_models[1], data == norway_data[1]) |>
  arrange(A_mean) |>
  pull(species)
norway_A <- mutate(norway_A, species = factor(species, levels = norway_order))
norway_arrows <- pivot_wider(norway_A, names_from = data, values_from = A_mean)

g_norway_A <- ggplot(norway_A, aes(A_mean, species)) +
  geom_segment(
    aes(x = .data[[norway_data[1]]], xend = .data[[norway_data[2]]], yend = species),
    data = norway_arrows, colour = "grey60",
    arrow = arrow(length = unit(4, "pt"), type = "closed")
  ) +
  geom_point(aes(colour = data)) +
  facet_wrap(~model) +
  scale_x_sqrt(breaks = c(0, 10, 50, 200, 500, 1000, 2000)) +
  scale_colour_manual(values = norway_colours) +
  # in the empty right side of the RW year panel
  theme(
    legend.position = "inside",
    legend.position.inside = c(inside_legend_x, 0.5),
    legend.justification = c(0, 0.5),
    legend.background = element_blank(),
    legend.key = element_blank()
  ) +
  labs(
    x = paste0("Mean A (%), ", norway_window, "-year windows"), y = NULL, colour = NULL,
    title = "Norwegian gillnet–fyke net: adding full coverage in 2024–2026"
  )

# The three most abundant fitted stocks: a typical seesaw that the full years
# shrink (pollack), one they overcorrect (cod; the South got relatively poorer
# by 2024-2026), and one with no seesaw to fix (saithe)
norway_ts_species <- c("Pollack", "Atlantic cod", "Saithe")
g_norway_ts <- norway |>
  filter(model == norway_models[1], species %in% norway_ts_species) |>
  mutate(
    species = factor(species, levels = norway_ts_species),
    log_est = log(est) - mean(log(est)),
    .by = c(species, data)
  ) |>
  ggplot(aes(year, log_est, colour = data)) +
  annotate("rect", xmin = 2023.5, xmax = 2026.5, ymin = -Inf, ymax = Inf, fill = "grey92") +
  geom_line() +
  geom_point() +
  facet_wrap(~species, scales = "free_y") +
  scale_colour_manual(values = norway_colours) +
  scale_x_continuous(breaks = c(2014, 2018, 2022, 2026)) +
  labs(
    x = NULL, y = "Centred log index", colour = NULL,
    title = "IID RF, factor(year); shaded: both regions sampled"
  ) +
  theme(legend.position = "none")

# free() so the time-series y-axis title isn't pushed left to line up with
# the species labels above
g <- (g_norway_A / free(g_norway_ts, type = "label", side = "l")) +
  plot_layout(heights = c(0.95, 0.85)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.title = element_text(size = rel(0.9)))
ggsave(here::here("figs/norway-full-coverage.pdf"), g, width = 7.5, height = 5, device = cairo_pdf)
