library(dplyr)
library(patchwork)
library(sf)
source(here::here("analysis/R/seesaw-metrics.R"))
source(here::here("analysis/R/plotting.R"))

norway_models <- c("IID RF, factor(year)", "IID RF, RW year")
norway_window <- 10L # calendar years
norway_data <- c(alternating = "Alternating only", all = "With full 2024-2026")
norway_colours <- set_names(c("#CC79A7", "#009E73"), norway_data)

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
    A_max = max(moving_window(est, year, phase, window = norway_window)$A),
    .by = c(species, model, data)
  )
norway_order <- norway_A |>
  filter(model == norway_models[1], data == norway_data[1]) |>
  arrange(A_max) |>
  pull(species)
norway_A <- mutate(norway_A, species = factor(species, levels = norway_order))
# Arrows stop short of the end point, leaving a gap before the arrowhead. The
# gap is in sqrt units so it is the same on screen for every arrow (the
# facets share the x scale); arrows too short for it are dropped.
arrow_gap <- 1.35
norway_arrows <- pivot_wider(norway_A, names_from = data, values_from = A_max) |>
  mutate(
    from = .data[[norway_data[1]]],
    to = .data[[norway_data[2]]],
    sqrt_shift = arrow_gap * sign(sqrt(from) - sqrt(to)),
    to_short = (sqrt(to) + sqrt_shift)^2
  ) |>
  filter(abs(sqrt(from) - sqrt(to)) > 2 * arrow_gap)

g_norway_A <- ggplot(norway_A, aes(A_max, species)) +
  geom_segment(
    aes(x = from, xend = to_short, yend = species),
    data = norway_arrows, colour = "grey60",
    arrow = arrow(length = unit(4, "pt"), type = "closed")
  ) +
  geom_point(aes(colour = data)) +
  facet_wrap(~model) +
  scale_x_sqrt(breaks = c(0, 50, 200, 500, 1000, 2000)) +
  # all species in factor order, though some have no arrow
  scale_y_discrete(drop = FALSE) +
  scale_colour_manual(values = norway_colours) +
  labs(
    x = paste0("Maximum biennial amplitude (%) across ", norway_window, "-year windows"), y = NULL, colour = NULL
  )

# The three most abundant fitted stocks: a typical seesaw that the full years
# shrink (pollack), one they overcorrect (cod; the South got relatively poorer
# by 2024-2026), and one with no seesaw to fix (saithe)
norway_ts_species <- c("Pollack", "Atlantic cod", "Saithe")

# Rows b and c plot log10 by hand rather than use scale_y_log10(), which would
# turn the shading's -Inf into NaN; the axis is labelled on the original scale
# Round values within each panel's range, thinned when the range is wide
relative_breaks <- function(lims) {
  b <- c(0.3, 0.5, 0.7, 1, 1.5, 2, 3)
  b <- b[b >= 10^lims[1] & b <= 10^lims[2]]
  if (length(b) > 5) b <- intersect(b, c(0.3, 0.5, 1, 2, 3))
  log10(b)
}
log10_axis <- scale_y_continuous(
  breaks = relative_breaks,
  labels = \(b) as.character(signif(10^b, 2))
)
both_sampled <- list(
  geom_rect(
    aes(xmin = xmin, xmax = xmax, fill = fill),
    data = tibble(xmin = 2023.5, xmax = 2026.5, fill = "Both regions sampled"),
    ymin = -Inf, ymax = Inf, inherit.aes = FALSE
  ),
  scale_fill_manual(values = c("Both regions sampled" = "grey92")),
  guides(colour = guide_legend(order = 1), fill = guide_legend(order = 2)),
  theme(legend.spacing.y = unit(0, "pt"), legend.margin = margin(2, 0, 2, 0))
)

# Each series relative to its geometric mean
g_norway_ts <- norway |>
  filter(model == norway_models[1], species %in% norway_ts_species) |>
  mutate(
    species = factor(species, levels = norway_ts_species),
    rel_est = est / exp(mean(log(est))),
    .by = c(species, data)
  ) |>
  ggplot(aes(year, log10(rel_est), colour = data)) +
  both_sampled +
  geom_line() +
  geom_point() +
  facet_wrap(~species, scales = "free_y") +
  log10_axis +
  scale_colour_manual(values = norway_colours) +
  scale_x_continuous(breaks = c(2015, 2020, 2025)) +
  labs(x = NULL, y = "Biomass index", colour = NULL, fill = NULL)

# Only the years each subarea was sampled. Relative to each species' geometric
# mean across both subareas, so the North-South gap is kept.
# as in fig-case-study-maps.R
subarea_colours <- c(North = "#0072B2", South = "#D55E00")
g_subarea <- readRDS(here::here("data-generated/norway-subarea.rds")) |>
  mutate(species = factor(species, levels = norway_ts_species)) |>
  mutate(across(c(est, lwr, upr), \(x) x / exp(mean(log(est)))), .by = species) |>
  ggplot(aes(year, log10(est), colour = subarea)) +
  both_sampled +
  geom_linerange(aes(ymin = log10(lwr), ymax = log10(upr)), position = position_dodge(width = 0.4), alpha = 0.6) +
  geom_line(position = position_dodge(width = 0.4)) +
  geom_point(position = position_dodge(width = 0.4)) +
  facet_wrap(~species, scales = "free_y") +
  log10_axis +
  scale_colour_manual(values = subarea_colours) +
  scale_x_continuous(breaks = c(2015, 2020, 2025)) +
  labs(x = NULL, y = "Biomass index", colour = NULL, fill = NULL)

# Mini map of the sites, coloured by subarea as in row c, in the same
# projection as the Norway map in fig-case-study-maps.R
crs_no <- "+proj=aea +lat_0=60 +lat_1=62.5 +lat_2=67 +lon_0=10 +datum=WGS84 +units=m +no_defs"
sf_use_s2(FALSE)
coast_no <- rnaturalearth::ne_countries(scale = 10, returnclass = "sf") |>
  st_geometry() |>
  st_make_valid() |>
  st_crop(st_bbox(c(xmin = -5, xmax = 30, ymin = 55, ymax = 75), crs = 4326)) |>
  suppressWarnings() |>
  st_transform(crs_no)
norway_sites <- readRDS(here::here("data-raw/garn_ruse_predgrid.rds")) |>
  ungroup() |>
  distinct(site, subarea, lon, lat) |>
  st_as_sf(coords = c("lon", "lat"), crs = 4326) |>
  st_transform(crs_no)
map_lims <- st_bbox(c(xmin = 3.4, xmax = 16.3, ymin = 61.6, ymax = 68.7), crs = 4326) |>
  st_as_sfc() |>
  st_segmentize(0.1) |>
  st_transform(crs_no) |>
  st_bbox()
# Labels as in fig-case-study-maps.R
map_labels <- tribble(
  ~label, ~lon, ~lat, ~kind,
  "North", 10.9, 65.9, "survey",
  "South", 5.7, 63.6, "survey",
  "Norway", 13.8, 64.4, "place"
)
map_xy <- st_coordinates(st_transform(st_as_sf(map_labels, coords = c("lon", "lat"), crs = 4326), crs_no))
map_labels <- mutate(map_labels, X = map_xy[, 1], Y = map_xy[, 2])
g_map <- ggplot() +
  geom_sf(data = coast_no, fill = "grey85", colour = NA) +
  geom_sf(aes(colour = subarea), data = norway_sites, size = 1) +
  scale_colour_manual(values = subarea_colours, guide = "none") +
  ggrepel::geom_text_repel(
    aes(X, Y, label = label),
    data = filter(map_labels, kind == "survey"), size = 2.4, fontface = "bold", colour = "grey15",
    bg.colour = "white", bg.r = 0.15, force = 0, force_pull = 0, min.segment.length = Inf
  ) +
  geom_text(
    aes(X, Y, label = label),
    data = filter(map_labels, kind == "place"), size = 2.2, fontface = "italic", colour = "grey35"
  ) +
  coord_sf(
    crs = crs_no, xlim = map_lims[c("xmin", "xmax")], ylim = map_lims[c("ymin", "ymax")],
    expand = FALSE
  ) +
  scale_x_continuous(breaks = seq(4, 16, 4)) +
  scale_y_continuous(breaks = seq(62, 68, 2)) +
  labs(x = NULL, y = NULL) +
  theme(
    axis.text = element_text(size = rel(0.7)),
    panel.grid = element_line(colour = "grey92")
  )

# Nesting the map and panel b keeps row a's long species names from setting the
# left edge of rows b and c
row_a <- (g_map | g_norway_A) +
  plot_layout(widths = c(1.33, 4.5))
# Space after the y title moves the panels below right, to line up with the map
y_title_gap <- theme(axis.title.y = element_text(margin = margin(r = 3)))
g <- (row_a / (g_norway_ts + y_title_gap) / (g_subarea + y_title_gap)) +
  plot_layout(heights = c(0.75, 0.85, 0.85)) +
  plot_annotation(tag_levels = "a") &
  # tags inside each plot's area rather than in an extra
  # row above and column to the left
  theme(
    axis.title = element_text(size = rel(0.9)),
    plot.tag.location = "plot",
    plot.tag.position = c(0, 1),
    plot.tag = element_text(hjust = 0, vjust = 1),
    plot.margin = margin(5.5, 5.5, 5.5, 2)
  )
ggsave(here::here("figs/norway-full-coverage.pdf"), g, width = 7.83, height = 5.99, device = cairo_pdf)
