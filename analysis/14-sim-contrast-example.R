# Example of what causes the seesaw: one simulated population with a
# north-south contrast in its spatial field, sampled one region per year
# (alternating north / south) and fit with two estimation models. With
# factor(year), the year effects absorb the regional contrast (they alternate
# with the region sampled) and the estimated spatial field loses most of it.
# With a random walk on the year effects, the contrast stays in the spatial
# field. The last panel places this example in the contrast grid from
# 10-sim-contrast-grid.R, whose setup, functions, and metrics are reused here
# (run that script first for its cached results).

lines <- readLines(here::here("analysis/10-sim-contrast-grid.R"))
eval(parse(text = lines[1:(grep("^# Run [(]one cache", lines) - 1)]))
library(patchwork)
source(here::here("analysis/plot-index-models.R")) # highlight_colours()

# One replicate of one cell of the contrast grid
delta <- 2
sigma_O <- 1
seed <- 1

# Simulate and fit -------------------------------------------------------------

# Cached; delete the file to resimulate and refit
cache_file <- here::here(sprintf(
  "data-generated/sim-contrast-example__delta-%s__sigma-O-%s__seed-%03d.rds",
  delta, sigma_O, seed
))
if (!file.exists(cache_file)) {
  s <- simulate_cell(delta, sigma_O, seed, make_square_mesh(grid, max_edge = 0.02))
  d <- observe(
    s,
    sample_n = obs_n, seed = seed, gap = gap, region_cutoff = region_cutoff,
    north_yrs = north_yrs, south_yrs = south_yrs
  )
  mesh <- make_square_mesh(d, max_edge = 0.05)
  fits <- purrr::map(setNames(models, models), \(m) fit_model(m, d, mesh))
  stopifnot(!purrr::some(fits, is.null))
  saveRDS(list(s = s, fits = fits), cache_file)
}
cached <- readRDS(cache_file)
s <- cached$s
fits <- cached$fits

field_levels <- c("Truth", "IID RF, RW year", "IID RF, factor(year)")

# Spatial effects --------------------------------------------------------------

# The true time-invariant spatial component (omega + gradient) is eta minus
# epsilon in any year; the year-1 effect is a constant removed by centring
g1 <- grid[grid$year == 1L, ]
spatial <- bind_rows(
  tibble::tibble(field = "Truth", s |> filter(year == 1L) |> transmute(X, Y, omega = eta - epsilon_st)),
  purrr::imap_dfr(fits, \(f, m) {
    tibble::tibble(field = m, transmute(predict(f, newdata = g1), X, Y, omega = omega_s))
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

# North - south difference in each field (for the caption)
step <- reg_means |>
  group_by(field) |>
  summarise(step = omega[region == "North"] - omega[region == "South"])
print(step)

# Colour-blind safe (Okabe-Ito), as in the case-study figures
field_colours <- setNames(c("black", rev(highlight_colours(2))), field_levels)

# Profile from south to north: mean across X within bands of Y
profile <- spatial |>
  mutate(band = cut(Y, breaks = seq(0, 1, by = 0.05), include.lowest = TRUE)) |>
  group_by(field, band) |>
  summarise(Y = mean(Y), omega = mean(omega), .groups = "drop")

# Year effects -----------------------------------------------------------------

# True: grid mean of eta - epsilon (year effect plus a constant spatial mean).
# Estimated: the non-random-field part of the linear predictor (factor(year)
# coefficients or the random walk). All centred across years.
sampled <- tibble::tibble(
  year = seq_len(n_year),
  region = if_else(year %in% north_yrs, "North", "South")
)
year_effects <- bind_rows(
  s |>
    group_by(year) |>
    summarise(effect = mean(eta - epsilon_st)) |>
    mutate(field = "Truth"),
  purrr::imap_dfr(fits, \(f, m) {
    predict(f, newdata = g1[rep(1L, n_year), ] |> mutate(year = seq_len(n_year))) |>
      transmute(year, effect = est_non_rf, field = m)
  })
) |>
  group_by(field) |>
  mutate(effect = effect - mean(effect)) |>
  ungroup() |>
  left_join(sampled, by = "year") |>
  mutate(field = factor(field, levels = field_levels))

# Contrast grid -----------------------------------------------------------------

# Metrics section of 10-sim-contrast-grid.R (excess seesaw A by cell and seed)
out <- readRDS(here::here("data-generated/contrast-grid.rds"))
eval(parse(text = lines[grep("^# Metrics", lines):(grep("^# Plots", lines) - 1)]))
grid_cells <- filter(cell_summary, sigma_O == !!sigma_O)
this_rep <- filter(metrics, delta == !!delta, sigma_O == !!sigma_O, seed == !!seed)

# Plot -------------------------------------------------------------------------

# Diverging scale centred on zero, symmetric limits at 16-fold; squishes the
# few most extreme fine-scale cells of the truth (~2-3%) so the regional
# pattern dominates
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
  facet_wrap(~field) +
  labs(x = "X", y = "Y", fill = "Spatial\neffect") +
  theme(axis.text = element_blank(), axis.ticks = element_blank())

# Beside the maps with Y vertical, as a marginal profile; the facet strip
# lines its panel up with the maps
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
  scale_x_continuous(breaks = seq(2, n_year, 2)) +
  scale_y_continuous(breaks = log(c(0.5, 1, 2)), labels = fold_labels) +
  labs(x = "Year", y = "Year effect", colour = NULL, shape = "Region\nsampled")

# Excess seesaw across the grid at this sigma_O (median and IQR across seeds);
# crosses: this example. Colours match the legend, so no legend of its own.
dodge <- position_dodge(width = 0.12)
g_grid <- ggplot(grid_cells, aes(delta, excess_med, colour = model)) +
  geom_hline(yintercept = 0, colour = "grey80") +
  geom_linerange(aes(ymin = excess_lo, ymax = excess_hi), position = dodge, show.legend = FALSE) +
  geom_line(position = dodge, show.legend = FALSE) +
  geom_point(position = dodge, show.legend = FALSE) +
  geom_point(aes(y = excess_A), data = this_rep, shape = 4, position = dodge, show.legend = FALSE) +
  scale_colour_manual(values = field_colours) +
  # Delta is a log difference: label as a fold difference at round ratios
  scale_x_continuous(
    breaks = log(c(1, 2, 4, 8)),
    labels = fold_labels,
    expand = expansion(mult = 0.08)
  ) +
  labs(x = "North:south density ratio", y = "Excess\nseesaw A (%)")

# Top: maps and the profile at the same size; bottom: the consequences. All
# legends in one left-justified column.
fig <- (g_maps + g_profile + plot_layout(widths = c(3, 1))) /
  (g_years + g_grid + plot_layout(widths = c(1.5, 1))) +
  plot_layout(heights = c(1, 0.85), guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(
    legend.justification = "left",
    strip.text = element_text(size = rel(0.9)),
    legend.text = element_text(size = rel(0.85))
  )

ggsave(here::here("figs/sim-contrast-example.pdf"), plot = fig, width = 8, height = 4.2)
ggsave(here::here("figs/sim-contrast-example.png"), plot = fig, width = 8, height = 4.2, dpi = 200)
