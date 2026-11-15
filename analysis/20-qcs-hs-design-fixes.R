# Can survey design fix the seesaw? Observed QCS + HS data, several species.
#
# As in analysis/07-qcs-hs-biennial-model-comparison.R, downgrade the full
# QCS + HS data to a pretend alternating design (one region per occasion) and
# compare against the index from the full data, but now with design tweaks:
#
# - full coverage in 1 or 2 occasions in the middle or at the end of an
#   otherwise alternating series
# - full coverage in the first 1 or 2 occasions, then alternating
# - overlap: each occasion, also sample the 10% or 25% of the off-region's tows
#   closest to the other region (i.e., the surveys overlap in the middle)
#
# No known truth is needed: the question is whether seesaw A and max A of each
# design's index move back to those of the full-data index.

source(here::here("analysis/self-test-funcs.R"))
theme_set(ggsidekick::theme_sleek())

surveyjoin::cache_data()
surveyjoin::load_sql_data()

# As in analysis/08-qcs-hs-biennial-self-test.R
species_to_fit <- c(
  "pacific cod",
  "spotted ratfish",
  "pacific ocean perch",
  "pacific halibut",
  "pacific spiny dogfish",
  "lingcod",
  "sablefish",
  "dover sole",
  "arrowtooth flounder",
  "shortspine thornyhead"
)
# Moving window (in survey occasions) for max A; as in analysis/07
seesaw_window <- 8L
models <- c(
  "IID RF, factor(year)",
  "IID RF, RW year"
)
n_workers <- max(1L, min(20L, future::availableCores()/2))

# Data and designs ---------------------------------------------------------------

qcs_hs_base_grid <- qcs_hs_grid()

design_group <- c(
  "Full data" = "Reference", "Alternating" = "Reference",
  "Full first 1" = "Full first", "Full first 2" = "Full first",
  "Full middle 1" = "Full added", "Full middle 2" = "Full added",
  "Full last 1" = "Full added", "Full last 2" = "Full added",
  "Overlap boundary 10%" = "Overlap", "Overlap boundary 25%" = "Overlap"
)
design_names <- names(design_group)

# Data, grid, mesh, and the tows kept by each design for one species
setup_designs <- function(species) {
  dat <- qcs_hs_data(species)
  grid <- replicate_df(qcs_hs_base_grid, "year", sort(unique(dat$year)))
  n_occ <- max(dat$year)
  mid <- ceiling(n_occ / 2)

  # Distance (km) from each tow to the nearest grid cell of the other region
  grid_xy <- distinct(qcs_hs_base_grid, X, Y, survey)
  dist_to_other <- numeric(nrow(dat))
  for (s in qcs_hs_surveys) {
    i <- dat$survey_name == s
    other <- as.matrix(grid_xy[grid_xy$survey != s, c("X", "Y")])
    dist_to_other[i] <- FNN::get.knnx(other, as.matrix(dat[i, c("X", "Y")]), k = 1L)$nn.dist[, 1L]
  }
  # Rank of each tow within its occasion and region: 0 = closest to the other
  # region, 1 = farthest
  dat$boundary_rank <- ave(dist_to_other, dat$year, dat$survey_name, FUN = \(x) (rank(x) - 1) / length(x))

  primary <- dat$survey_name == dat$primary_survey
  full_in <- \(occ) primary | dat$year %in% occ
  overlap <- \(rank, p) primary | rank < p

  designs <- list(
    "Full data" = rep(TRUE, nrow(dat)),
    "Alternating" = primary,
    "Full first 1" = full_in(1),
    "Full first 2" = full_in(1:2),
    "Full middle 1" = full_in(mid),
    "Full middle 2" = full_in(mid + 0:1),
    "Full last 1" = full_in(n_occ),
    "Full last 2" = full_in(n_occ - 1:0),
    "Overlap boundary 10%" = overlap(dat$boundary_rank, 0.1),
    "Overlap boundary 25%" = overlap(dat$boundary_rank, 0.25)
  )
  stopifnot(identical(names(designs), design_names))
  list(
    dat = dat, grid = grid, mid = mid, designs = designs,
    mesh = make_mesh(dat, c("X", "Y"), cutoff = 10)
  )
}

setups <- purrr::map(rlang::set_names(species_to_fit), setup_designs)

# Proportion of all tows kept by each design
purrr::map_dfr(setups, \(x) as.list(sapply(x$designs, mean)), .id = "species") |>
  mutate(across(where(is.numeric), \(x) round(x, 2))) |>
  print(width = Inf)

# Sampled tows by occasion, region, and design (one page per species)
pdf(here::here("figs/qcs-hs-design-fixes-tows.pdf"), width = 10, height = 7)
for (sp in species_to_fit) {
  x <- setups[[sp]]
  g <- purrr::imap_dfr(x$designs, \(keep, design) {
    x$dat[keep, ] |>
      count(year, survey_name) |>
      mutate(design = design)
  }) |>
    mutate(design = factor(design, levels = design_names)) |>
    ggplot(aes(year, n, fill = survey_name)) +
    geom_col() +
    facet_wrap(~design) +
    labs(x = "Survey occasion", y = "Tows", fill = "Region", title = stringr::str_to_sentence(sp))
  print(g)
}
dev.off()

# Map of one occasion's sampled tows for the overlap designs (tow locations
# barely differ among species, so one species is enough)
x <- setups[[1]]
x$dat |>
  filter(year == x$mid) |>
  mutate(
    `Overlap boundary 10%` = x$designs[["Overlap boundary 10%"]][x$dat$year == x$mid],
    `Overlap boundary 25%` = x$designs[["Overlap boundary 25%"]][x$dat$year == x$mid]
  ) |>
  tidyr::pivot_longer(starts_with("Overlap"), names_to = "design", values_to = "sampled") |>
  ggplot(aes(X, Y, colour = survey_name, alpha = sampled)) +
  geom_point() +
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.1)) +
  coord_fixed() +
  facet_wrap(~design) +
  labs(colour = "Region", alpha = "Sampled", title = paste("Occasion", x$mid))
ggsave(here::here("figs/qcs-hs-design-fixes-overlap-map.pdf"), width = 9, height = 6)

# Fit ------------------------------------------------------------------------------

fit_design <- function(species, design, model, keep, dat, grid, mesh) {
  single_thread()
  d <- dat[keep, ]
  # same knots for every design
  m <- make_mesh(d, c("X", "Y"), mesh = mesh$mesh)
  fit <- fit_ok(fit_model(model, d, m))
  index <- if (!is.null(fit)) get_index_ok(fit, grid)
  if (is.null(index)) {
    return(tibble::tibble(species = species, design = design, model = model, converged = FALSE))
  }
  mutate(index, species = species, design = design, model = model, converged = TRUE, .before = 1L)
}

# One cache file per design and model; delete a file to refit it
fits_dir <- here::here("data-generated/qcs-hs-design-fixes")
dir.create(fits_dir, showWarnings = FALSE)
tasks <- tidyr::expand_grid(species = species_to_fit, design = design_names, model = models) |>
  filter(design != "Full data" | model == "IID RF, factor(year)") |>
  mutate(file = file.path(fits_dir, paste0(slug(species), "__", slug(design), "__", slug(model), ".rds")))

to_fit <- filter(tasks, !file.exists(file))
if (nrow(to_fit) > 0L) {
  future::plan(future::multisession, workers = min(n_workers, nrow(to_fit)))
  # Pass each task's own data so workers don't get every species' setup
  furrr::future_pwalk(
    list(
      to_fit$species, to_fit$design, to_fit$model, to_fit$file,
      purrr::map2(to_fit$species, to_fit$design, \(sp, d) setups[[sp]]$designs[[d]]),
      purrr::map(to_fit$species, \(sp) setups[[sp]]$dat),
      purrr::map(to_fit$species, \(sp) setups[[sp]]$grid),
      purrr::map(to_fit$species, \(sp) setups[[sp]]$mesh)
    ),
    \(species, design, model, file, keep, dat, grid, mesh) {
      saveRDS(fit_design(species, design, model, keep, dat, grid, mesh), file)
    },
    .options = furrr::furrr_options(seed = TRUE, scheduling = Inf, packages = c("dplyr", "sdmTMB"))
  )
  future::plan(future::sequential)
}

occasions <- purrr::map_dfr(setups, \(x) distinct(x$dat, year, cal_year, primary_survey, phase), .id = "species")
# Fits cached before species was a column don't have it
indexes <- purrr::map2_dfr(tasks$file, tasks$species, \(f, sp) mutate(readRDS(f), species = sp, .before = 1L)) |>
  left_join(occasions, by = c("species", "year"))

filter(indexes, !converged) |> distinct(species, design, model)
indexes <- filter(indexes, converged)

# Seesaw metrics -------------------------------------------------------------------

# Phase follows which region the alternating design samples (also in occasions
# where a design samples both)
seesaw_A <- function(d) {
  d <- arrange(d, year)
  mw <- moving_window(d$est, year = d$year, window = seesaw_window, phase = d$phase)
  as_tibble(t(period2_metric(d$est, year = d$year, phase = d$phase))) |>
    mutate(A_max = max(mw$A), A_mean = mean(mw$A))
}

reference <- filter(indexes, design == "Full data")

metrics <- indexes |>
  group_by(species, design, model) |>
  group_modify(\(.x, .y) seesaw_A(.x)) |>
  ungroup()
metrics <- metrics |>
  left_join(
    indexes |>
      left_join(select(reference, species, year, ref_log_est = log_est), by = c("species", "year")) |>
      group_by(species, design, model) |>
      summarise(rmse_vs_full = sqrt(mean((log_est - ref_log_est)^2)), .groups = "drop"),
    by = c("species", "design", "model")
  ) |>
  left_join(
    filter(metrics, design == "Full data") |> select(species, ref_A = A, ref_A_max = A_max, ref_A_mean = A_mean),
    by = "species"
  ) |>
  mutate(
    group = design_group[design],
    species = factor(species, levels = species_to_fit),
    design = factor(design, levels = design_names),
    model = factor(model, levels = models)
  ) |>
  arrange(species, model, design)

saveRDS(list(indexes = indexes, metrics = metrics), here::here("data-generated/qcs-hs-design-fixes.rds"))

metrics |>
  select(species, model, design, beta, A, A_lower, A_upper, A_max, A_mean, rmse_vs_full) |>
  mutate(across(where(is.numeric), \(x) round(x, 2))) |>
  print(n = Inf, width = Inf)

# Plots ------------------------------------------------------------------------------

metric_labels <- c(A = "Seesaw A (%)", A_max = "Max A (%)", A_mean = "Mean A (%)")
metrics_long <- metrics |>
  filter(design != "Full data") |>
  tidyr::pivot_longer(c(A, A_max, A_mean), names_to = "metric") |>
  mutate(
    ref = case_when(metric == "A" ~ ref_A, metric == "A_max" ~ ref_A_max, metric == "A_mean" ~ ref_A_mean),
    metric = factor(metric, levels = names(metric_labels), labels = metric_labels),
    design = factor(design, levels = rev(design_names))
  )

# A of each design's index by species; dashed: the full-data index
for (m in metric_labels["A"]) {
  d <- filter(metrics_long, metric == m)
  ggplot(d, aes(value, design, colour = model)) +
    geom_vline(aes(xintercept = ref), data = distinct(d, species, ref), linetype = 2) +
    geom_point() +
    facet_wrap(~species, ncol = 5, scales = "free_x") +
    labs(
      x = m, y = NULL, colour = "Model",
      caption = paste0(
        "Observed QCS + HS data. Dashed: IID RF, factor(year) fit to the full data. ",
        "Max A over ", seesaw_window, "-occasion windows."
      )
    ) +
    theme(legend.position = "bottom")
  ggsave(here::here("figs/qcs-hs-design-fixes-A.pdf"), width = 14, height = 7)
}

# Max (or mean) A over seesaw_window-occasion windows of each design's index
# minus that of the full-data index: one point per species, violin across
# species, crossbar at the median (as in the summary figure of
# analysis/08-qcs-hs-biennial-self-test.R). Pseudo-log x so the large
# alternating seesaws don't squash the rest; also a linear x with values above
# x_max pinned just inside it (triangles). Mean A is less noisy than max A
# (the max of only a few windows picks the worst by chance), and each window
# is detrended locally, unlike whole-series A.
plot_A_max <- function(x_max = Inf, metric_name = "A_max", connect = FALSE, n_highlight = 3L) {
  d <- metrics_long |>
    filter(metric == metric_labels[[metric_name]]) |>
    mutate(species = factor(species, levels = species_to_fit), diff = value - ref)
  ref_model <- d |> filter(model == "IID RF, factor(year)")
  # With connect, highlight the species furthest from the full data under
  # factor(year), plus arrowtooth flounder as a control with almost no seesaw;
  # the rest are light grey
  highlight <- ref_model |>
    summarise(dev = mean(abs(diff)), .by = species) |>
    slice_max(dev, n = n_highlight) |>
    pull(species) |>
    as.character() |>
    c("arrowtooth flounder") |>
    unique()
  d <- d |>
    mutate(
      pinned = diff > x_max,
      # just inside x_max so the triangles aren't hidden by the panel border
      diff = if_else(pinned, x_max * 0.985, diff),
      design = droplevels(design),
      # with connect, a fixed vertical offset per species (instead of random
      # jitter) so each species' line stays in its own lane
      y = as.integer(design) + 0.4 * ((as.integer(species) - 1) / (nlevels(species) - 1) - 0.5),
      species_hl = factor(if_else(species %in% highlight, as.character(species), "other"), levels = c(highlight, "other"))
    ) |>
    arrange(species, design)
  colour_scale <- if (connect) {
    # Okabe-Ito
    scale_colour_manual(
      values = c(setNames(c("#D55E00", "#0072B2", "#009E73", "#E69F00", "#CC79A7")[seq_along(highlight)], highlight), other = "grey80"),
      breaks = highlight, labels = tools::toTitleCase
    )
  } else {
    # Paired has 12 colours (Dark2 has only 8)
    scale_colour_brewer(palette = "Paired", labels = tools::toTitleCase)
  }
  points <- if (connect) {
    list(
      geom_path(aes(y = y, colour = species_hl, group = species), data = filter(d, species_hl == "other")),
      geom_point(aes(y = y, colour = species_hl, shape = pinned), data = filter(d, species_hl == "other")),
      geom_path(aes(y = y, colour = species_hl, group = species), data = filter(d, species_hl != "other")),
      geom_point(aes(y = y, colour = species_hl, shape = pinned), data = filter(d, species_hl != "other"))
    )
  } else {
    geom_point(aes(colour = species, shape = pinned), position = position_jitter(height = 0.2, width = 0, seed = 1))
  }
  ggplot(d, aes(diff, design)) +
    geom_violin(fill = "grey90", colour = NA, scale = "width") +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
    points +
    stat_summary(fun = median, geom = "crossbar", fun.min = median, fun.max = median, width = 0.6) +
    facet_wrap(~model) +
    colour_scale +
    scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17), guide = "none") +
    labs(
      x = paste(
        if (metric_name == "A_max") "Max A minus full-data max A" else "Mean A minus full-data mean A",
        "(percentage points)"
      ),
      y = NULL, colour = NULL
    ) +
    theme(legend.position = "bottom") +
    guides(colour = guide_legend(nrow = if (connect) 1 else 2))
}
plot_A_max() +
  scale_x_continuous(
    transform = scales::pseudo_log_trans(sigma = 5),
    breaks = c(-50, -20, -5, 0, 5, 20, 50, 200, 800)
  )
ggsave(here::here("figs/qcs-hs-design-fixes-A-max.pdf"), width = 9, height = 5.5)
# No padding above x_max; the pinned values sit just inside the edge
plot_A_max(x_max = 200) +
  scale_x_continuous(limits = c(NA, 200), expand = expansion(mult = c(0.05, 0)))
ggsave(here::here("figs/qcs-hs-design-fixes-A-max-linear.pdf"), width = 9, height = 5.5)
# Left edge of the legends inside the RW year panels of a and b (same x so the
# keys line up)
inside_legend_x <- 0.71
g_qcs_hs <- plot_A_max(x_max = 200, metric_name = "A_mean") +
  scale_x_continuous(limits = c(NA, 200), expand = expansion(mult = c(0.05, 0))) +
  ggtitle("QCS + HS synoptic trawl: partial-coverage designs vs. full coverage") +
  # species legend in the empty right side of the RW year panel
  guides(colour = guide_legend(ncol = 1)) +
  theme(
    legend.position = "inside",
    legend.position.inside = c(inside_legend_x, 0.5),
    legend.justification = c(0, 0.5),
    legend.background = element_blank(),
    legend.key = element_blank()
  )

# What actually happened in the Norwegian gillnet-fyke net survey when both
# regions were sampled in 2024-2026 (analysis/19-fit-norway-models.R; run that
# first). Mean A over 10-year windows, as in that script. Phase is the region
# the alternating design samples, as in analysis/09-case-study-summary.R.
norway_models <- c("IID RF, factor(year)", "IID RF, RW year")
norway_window <- 10L
norway_data <- c("Alternating only", "With full 2024-2026")
norway_colours <- setNames(c("#D55E00", "#0072B2"), norway_data)
norway_sites <- readRDS(here::here("data-raw/garn_ruse_predgrid.rds")) |>
  ungroup() |>
  distinct(site, subarea)
norway_phase <- readRDS(here::here("data-raw/garn_ruse_catch_21sp-2026.rds")) |>
  inner_join(norway_sites, by = "site") |>
  filter(year < 2024) |>
  distinct(year, subarea) |>
  mutate(subarea = as.character(subarea)) |>
  bind_rows(tibble::tibble(year = 2024:2026, subarea = c("North", "South", "North"))) |>
  mutate(phase = ifelse(subarea == "North", 0.5, -0.5)) |>
  select(year, phase)
norway <- bind_rows(
  readRDS(here::here("data-generated/norway-fit-indexes.rds")) |> mutate(data = norway_data[1]),
  readRDS(here::here("data-generated/norway-fit-indexes-all-data.rds")) |> mutate(data = norway_data[2])
) |>
  filter(model %in% norway_models) |>
  left_join(norway_phase, by = "year") |>
  # only species fit under both designs (lemon sole fails on the alternating data)
  filter(n_distinct(data) == 2L, .by = c(species, model)) |>
  mutate(
    model = factor(model, levels = norway_models),
    data = factor(data, levels = norway_data),
    species = replace(species, species == "Pollack", "Pollock")
  )
norway_A <- norway |>
  arrange(year) |>
  summarise(
    A_mean = mean(moving_window(est, year = year, window = norway_window, phase = phase)$A),
    .by = c(species, model, data)
  )
norway_order <- norway_A |>
  filter(model == norway_models[1], data == norway_data[1]) |>
  arrange(A_mean) |>
  pull(species)
norway_A <- mutate(norway_A, species = factor(species, levels = norway_order))
norway_arrows <- tidyr::pivot_wider(norway_A, names_from = data, values_from = A_mean)

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
  # legend (shared with c) in the empty right side of the RW year panel
  theme(
    legend.position = "inside",
    legend.position.inside = c(inside_legend_x, 0.5),
    legend.justification = c(0, 0.5),
    legend.background = element_blank(),
    legend.key = element_blank()
  ) +
  labs(
    x = paste0("Mean A (%), ", norway_window, "-year windows"), y = NULL, colour = NULL,
    title = "Norwegian gillnet\u2013fyke net: adding full coverage in 2024\u20132026"
  )

# Time series for the three most abundant fitted stocks (by total catch),
# ordered to tell the story: a typical seesaw that the full years shrink
# (pollock), one they overcorrect (cod; the South got relatively poorer by
# 2024-2026), and one with no seesaw to fix (saithe)
norway_ts_species <- c("Pollock", "Atlantic cod", "Saithe")
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

library(patchwork)
# free() the y-axis title in (c) so it isn't pushed left to line up with the
# long design labels in (a)
(g_qcs_hs / g_norway_A / free(g_norway_ts, type = "label", side = "l")) +
  plot_layout(heights = c(1.7, 0.95, 0.85)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.title = element_text(size = rel(0.9)))
ggsave(here::here("figs/qcs-hs-design-fixes-A-mean-linear.pdf"), width = 7.5, height = 9, device = cairo_pdf)

# Norwegian panels only (b and c above, retagged a and b)
(g_norway_A / free(g_norway_ts, type = "label", side = "l")) +
  plot_layout(heights = c(0.95, 0.85)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.title = element_text(size = rel(0.9)))
ggsave(here::here("figs/qcs-hs-design-fixes-A-mean-linear-norway.pdf"), width = 7.5, height = 5, device = cairo_pdf)

plot_A_max(x_max = 200, metric_name = "A_mean", connect = TRUE) +
  scale_x_continuous(limits = c(NA, 200), expand = expansion(mult = c(0.05, 0)))
ggsave(here::here("figs/qcs-hs-design-fixes-A-mean-linear-connected.pdf"), width = 7.5, height = 4.8)

# Across species: each design's A minus that of the full-data index
metrics_long |>
  ggplot(aes(value - ref, design, colour = model)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.7), width = 0.6) +
  geom_point(position = position_jitterdodge(jitter.height = 0.1, jitter.width = 0, dodge.width = 0.7), alpha = 0.6) +
  facet_wrap(~metric, scales = "free_x") +
  labs(
    x = "Design minus full data (percentage points)", y = NULL, colour = "Model",
    caption = paste0("One point per species (", length(species_to_fit), " species), observed QCS + HS data.")
  )
ggsave(here::here("figs/qcs-hs-design-fixes-A-diff.pdf"), width = 10, height = 5)

# Indexes (one page per species)
pdf(here::here("figs/qcs-hs-design-fixes-indexes.pdf"), width = 8, height = 16)
for (sp in species_to_fit) {
  reference_line <- filter(reference, species == sp) |>
    select(cal_year, ref_est = est, ref_lwr = lwr, ref_upr = upr)
  g <- indexes |>
    filter(species == sp, design != "Full data") |>
    left_join(reference_line, by = "cal_year") |>
    mutate(design = factor(design, levels = design_names), model = factor(model, levels = models)) |>
    ggplot(aes(cal_year, est)) +
    geom_ribbon(aes(ymin = ref_lwr, ymax = ref_upr), fill = "grey85") +
    geom_line(aes(y = ref_est), colour = "grey40") +
    geom_line(colour = "grey50", linetype = 2) +
    geom_point(aes(colour = primary_survey)) +
    facet_grid(design ~ model) +
    labs(
      x = "Year", y = "Biomass index", colour = "Alternating region",
      title = stringr::str_to_sentence(sp),
      caption = "Grey: IID RF, factor(year) fit to the full data (95% CI)."
    )
  print(g)
}
dev.off()
