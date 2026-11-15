library(dplyr)
library(patchwork)
source(here::here("analysis/R/seesaw-metrics.R"))
source(here::here("analysis/R/plotting.R"))
source(here::here("analysis/R/gradient-sim.R"))

example <- readRDS(here::here("data-generated/gradient-example.rds"))
s <- example$s
fits <- example$fits
# as in fit-gradient-example.R
delta <- 2
sigma_O <- 1
seed <- 1

field_levels <- c("Truth", "IID RF, RW year", "IID RF, factor(year)")

# Spatial effects --------------------------------------------------------------

# The true time-invariant spatial component (omega + gradient) is eta minus
# epsilon in any year; the year-1 effect is a constant removed by centring
g1 <- gradient_grid[gradient_grid$year == 1L, ]
spatial <- bind_rows(
  tibble(field = "Truth", s |> filter(year == 1L) |> transmute(X, Y, omega = eta - epsilon_st)),
  imap_dfr(fits, \(f, m) {
    tibble(field = m, transmute(predict(f, newdata = g1), X, Y, omega = omega_s))
  })
) |>
  mutate(
    field = factor(field, levels = field_levels),
    region = if_else(Y >= region_cutoff, "North", "South")
  ) |>
  group_by(field) |>
  mutate(omega = omega - mean(omega)) |> # only the pattern is identified
  ungroup()

reg_means <- spatial |>
  group_by(field, region) |>
  summarise(omega = mean(omega), .groups = "drop")

# for the caption
step <- reg_means |>
  group_by(field) |>
  summarise(step = omega[region == "North"] - omega[region == "South"])
print(step)

field_colours <- setNames(c("black", rev(highlight_colours(2))), field_levels)

profile <- spatial |>
  mutate(band = cut(Y, breaks = seq(0, 1, by = 0.05), include.lowest = TRUE)) |>
  group_by(field, band) |>
  summarise(Y = mean(Y), omega = mean(omega), .groups = "drop")

# Year effects -----------------------------------------------------------------

# True: grid mean of eta - epsilon (year effect plus a constant spatial mean).
# Estimated: the non-random-field part of the linear predictor (factor(year)
# coefficients or the random walk). All centred across years.
sampled <- tibble(
  year = seq_len(gradient_n_year),
  region = if_else(year %in% north_years, "North", "South")
)
year_effects <- bind_rows(
  s |>
    group_by(year) |>
    summarise(effect = mean(eta - epsilon_st)) |>
    mutate(field = "Truth"),
  imap_dfr(fits, \(f, m) {
    predict(f, newdata = g1[rep(1L, gradient_n_year), ] |> mutate(year = seq_len(gradient_n_year))) |>
      transmute(year, effect = est_non_rf, field = m)
  })
) |>
  group_by(field) |>
  mutate(effect = effect - mean(effect)) |>
  ungroup() |>
  left_join(sampled, by = "year") |>
  mutate(field = factor(field, levels = field_levels))

# Gradient grid ------------------------------------------------------------------

m <- gradient_grid_metrics(readRDS(here::here("data-generated/gradient-grid.rds")), window = 8L)
grid_cells <- filter(m$by_cell, sigma_O == !!sigma_O)
this_rep <- filter(m$by_seed, delta == !!delta, sigma_O == !!sigma_O, seed == !!seed)

# Plot -------------------------------------------------------------------------

# Squishes the few most extreme fine-scale cells of the truth (~2-3%) so the
# regional pattern dominates
fill_lim <- log(16)

# Effects are centred log values; label them as fold differences from the mean
# (geometric mean across space or years); below 1 as fractions (1/4x, not 0.25x)
fold_labels <- \(x) {
  fold <- exp(x)
  paste0(if_else(fold < 1, paste0("1/", round(1 / fold)), as.character(round(fold))), "x")
}

g_maps <- ggplot(spatial, aes(X, Y, fill = omega)) +
  geom_raster() +
  geom_hline(yintercept = region_cutoff, linetype = 2) +
  scale_fill_distiller(
    palette = "RdBu", direction = -1,
    limits = c(-fill_lim, fill_lim), oob = scales::squish,
    breaks = log(c(1 / 16, 1 / 4, 1, 4, 16)), labels = fold_labels
  ) +
  coord_equal(expand = FALSE) +
  ggh4x::facet_wrap2(
    ~field,
    strip = ggh4x::strip_themed(
      text_x = map(field_colours, \(x) element_text(colour = x))
    )
  ) +
  labs(x = "X", y = "Y", fill = "Spatial\neffect") +
  theme(axis.text = element_blank(), axis.ticks = element_blank())

# The facet strip lines this panel up with the maps
g_profile <- ggplot(profile, aes(omega, Y, colour = field)) +
  geom_vline(xintercept = 0, colour = "grey80") +
  geom_hline(yintercept = region_cutoff, linetype = 2) +
  geom_path() +
  scale_colour_manual(values = field_colours) +
  scale_x_continuous(breaks = log(c(0.25, 1, 4)), labels = fold_labels) +
  scale_y_continuous(expand = c(0, 0)) +
  facet_wrap(~"South-north profile") +
  labs(x = "Mean spatial effect", y = NULL, colour = NULL) +
  theme(aspect.ratio = 1, axis.text.y = element_blank(), axis.ticks.y = element_blank())

g_years <- ggplot(year_effects, aes(year, effect, colour = field)) +
  geom_hline(yintercept = 0, colour = "grey80") +
  geom_line() +
  geom_point(aes(shape = region), show.legend = c(colour = FALSE)) + # one colour legend
  scale_colour_manual(values = field_colours) +
  scale_shape_manual(values = c(North = 16, South = 1)) +
  scale_x_continuous(breaks = seq(2, gradient_n_year, 2)) +
  scale_y_continuous(breaks = log(c(0.5, 1, 2)), labels = fold_labels) +
  labs(x = "Year", y = "Year effect", colour = NULL, shape = "Region\nsampled")

dodge <- position_dodge(width = 0.12)
g_grid <- ggplot(grid_cells, aes(delta, seesaw_med, colour = model)) +
  geom_hline(yintercept = 0, colour = "grey80") +
  geom_linerange(aes(ymin = seesaw_lo, ymax = seesaw_hi), position = dodge, show.legend = FALSE) +
  geom_line(position = dodge, show.legend = FALSE) +
  geom_point(position = dodge, show.legend = FALSE) +
  geom_point(aes(y = seesaw_error), data = this_rep, shape = 4, position = dodge, show.legend = FALSE) +
  scale_colour_manual(values = field_colours) +
  scale_x_continuous(
    breaks = log(c(1, 2, 4, 8)),
    labels = fold_labels,
    expand = expansion(mult = 0.08)
  ) +
  labs(x = "North:south density ratio", y = "Signed seesaw\nerror (%)")

fig <- (g_maps + g_profile + plot_layout(widths = c(3, 1))) /
  (g_years + g_grid + plot_layout(widths = c(1.5, 1))) +
  plot_layout(heights = c(1, 0.85), guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(
    legend.justification = "left",
    strip.text = element_text(size = rel(0.9)),
    legend.text = element_text(size = rel(0.85))
  )

ggsave(here::here("figs/gradient-mechanism.pdf"), plot = fig, width = 8, height = 4.2)
