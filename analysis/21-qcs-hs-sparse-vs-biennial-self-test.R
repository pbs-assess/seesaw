# Simulation self-test: is it better to sample the full QCS + HS domain
# sparsely every year, or to alternate between HS and QCS every other occasion?
#
# As in analysis/08-qcs-hs-biennial-self-test.R, but with two designs of equal
# effort compared against the known truth:
#
# 1. Fit the operating model (IID RF, factor(year)) to the full QCS + HS data.
# 2. For each replicate, draw the spatial fields from their approximate
#    posterior and entirely new spatiotemporal fields and observation error at
#    the observed tows; the true index is the area-weighted sum of the expected
#    density on the grid under the same field draws.
# 3. Two designs, with the same number of tows in each occasion (so the same
#    average effort across the series):
#    - alternating: all tows in one region per occasion (HS / QCS alternating),
#      fit with IID RF + RW year (with the sigma_V prior) and with IID RF +
#      factor(year)
#    - overlap 10% / 25%: the alternating design, but with the 10% or 25% of
#      the other region's tows closest to the sampled region added and the
#      same number of tows dropped at random from the sampled region, so the
#      surveys overlap at the boundary; each fit with IID RF + factor(year)
#      and with IID RF + RW year
#    - sparse annual: a random subsample of the tows of the whole domain
#      (both regions) in every occasion, a different draw for each replicate,
#      fit with IID RF + factor(year) and with IID RF + RW year
#    Also the operating model refit to all tows, as the baseline.
# 4. Score RMSE, 50% (and 95%) CI coverage of the bias-corrected index, seesaw
#    A (and max A), and the signed seesaw error against the truth of each
#    replicate (self_test_metrics(), as in 08).
#
# The pipeline is in analysis/self-test-funcs.R and everything is cached per
# replicate, so n_reps can be raised later and only the new replicates are run.

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
n_reps <- 20L
# Moving window (in survey occasions) for max A; as in 08
seesaw_window <- 8L
n_workers <- max(1L, min(60L, future::availableCores() - 2L))

# Operating model; IID fields for the reasons given in analysis/08
om_model <- "IID RF, factor(year)"
sparse_models <- c("IID RF, factor(year)", "IID RF, RW year")
alternating_models <- c("IID RF, factor(year)", "IID RF, RW year")
# Baseline for the within-replicate RMSE ratios
alternating_model <- "IID RF, RW year"
# Overlap designs: percentage of the off region's tows sampled at the boundary
overlap_pct <- c(10, 25)
overlap_models <- c("IID RF, factor(year)", "IID RF, RW year")
# All tows, the baseline (same as the operating model)
full_model <- "IID RF, factor(year)"

# Labels, in the order plotted; the alternating label has no suffix from
# model_label()
sparse_levels <- paste0(sparse_models, " (sparse annual)")
alternating_levels <- paste0(alternating_models, " (alternating)")
alternating_level <- paste0(alternating_model, " (alternating)")
full_level <- paste0(full_model, " (full data)")
overlap_levels <- paste0(rep(overlap_models, times = length(overlap_pct)), " (overlap ", rep(overlap_pct, each = length(overlap_models)), "%)")
design_levels <- c(full_level, sparse_levels, overlap_levels, alternating_levels)

# Run ------------------------------------------------------------------------

qcs_hs_base_grid <- qcs_hs_grid()

setup_species <- function(species, dat) {
  setup <- qcs_hs_setup(species, dat, qcs_hs_base_grid)
  om <- fit_operating_model(om_model, dat, setup$mesh)
  if (isTRUE(om$om_failed)) return(c(list(species = species), om))
  c(setup, list(fit_om = om$fit, ran_pars = om$ran_pars))
}

dir <- here::here("data-generated/qcs-hs-sparse-vs-biennial-self-test")
results <- run_self_test(
  dir = dir,
  species = species_to_fit,
  n_reps = n_reps,
  prepare_data = qcs_hs_data,
  setup_species = setup_species,
  simulate_rep = simulate_one_arm("spatiotemporal"),
  models = tibble::tibble(
    design = c(
      "full", rep("sparse", length(sparse_models)),
      rep(paste0("overlap", overlap_pct), each = length(overlap_models)),
      rep("alternating", length(alternating_models))
    ),
    model = c(full_model, sparse_models, rep(overlap_models, times = length(overlap_pct)), alternating_models)
  ),
  n_workers = n_workers
)

# Metrics --------------------------------------------------------------------

# self_test_metrics() groups by the model label; the alternating fits are the
# ones with no design suffix, so add one
label_design <- \(model) {
  factor(if_else(grepl(" \\(.*\\)$", model), model, paste0(model, " (alternating)")), levels = design_levels)
}
m <- self_test_metrics(compare_to_truth(results), results, seesaw_window)
compared <- mutate(m$compared, design = label_design(model))
metrics_rep <- mutate(m$metrics_rep, design = label_design(model))

metrics <- metrics_rep |>
  group_by(species, design) |>
  summarise(
    n_reps = n(),
    across(c(rmse, bias, coverage_50, coverage_95), mean),
    across(c(A, A_max, seesaw_error), median),
    .groups = "drop"
  )

saveRDS(
  list(
    compared = compared, metrics_rep = metrics_rep, metrics = metrics, A_truth = m$A_truth,
    failure_rates = results$failure_rates
  ),
  here::here("data-generated/qcs-hs-sparse-vs-biennial-self-test-metrics.rds")
)

metrics |>
  mutate(across(where(is.numeric), \(x) round(x, 2))) |>
  print(n = Inf, width = Inf)

# Within-replicate (same simulated population) RMSE of each sparse annual and
# overlap design relative to the alternating design: geometric mean ratio
# across replicates, and the proportion of replicates in which it has the
# lower RMSE. Ratio < 1: better than alternating.
metrics_rep |>
  select(species, rep, design, rmse) |>
  tidyr::pivot_wider(names_from = design, values_from = rmse) |>
  tidyr::pivot_longer(all_of(c(sparse_levels, overlap_levels)), names_to = "design", values_to = "rmse_design") |>
  summarise(
    rmse_ratio = exp(mean(log(rmse_design / .data[[alternating_level]]))),
    prop_reps_better = mean(rmse_design < .data[[alternating_level]]),
    .by = c(species, design)
  ) |>
  mutate(across(where(is.numeric), \(x) round(x, 2))) |>
  print(n = Inf, width = Inf)

# Plots ----------------------------------------------------------------------

metric_labels <- c(
  rmse = "RMSE (log)", coverage_50 = "50% CI coverage", A = "Seesaw A (%)",
  A_max = paste0("Max A, ", seesaw_window, "-survey windows (%)"),
  seesaw_error = "Seesaw error (%)"
)
# Nominal values; for A and max A, the median of the true index by species
nominal <- bind_rows(
  tidyr::expand_grid(
    species = unique(metrics_rep$species),
    tibble::tibble(metric = c("coverage_50", "seesaw_error"), value = c(0.5, 0))
  ),
  m$A_truth |>
    group_by(species) |>
    summarise(A = median(A), A_max = median(A_max)) |>
    tidyr::pivot_longer(c(A, A_max), names_to = "metric")
) |>
  mutate(metric = factor(metric, levels = names(metric_labels), labels = metric_labels))

rep_long <- metrics_rep |>
  tidyr::pivot_longer(all_of(names(metric_labels)), names_to = "metric") |>
  mutate(metric = factor(metric, levels = names(metric_labels), labels = metric_labels))

# Summary across species (the signed seesaw error is estimated minus true
# seesaw, positive when the index reads high in occasions the richer region is
# sampled; see self_test_metrics()). One point per species (mean across
# replicates; median for the seesaw error), connected between designs, over the
# violins of all replicates pooled across species.
#
# With split = TRUE, rows are the year model (RW year, factor(year)); the full
# data fit is shown in both rows
split_year_model <- function(d) {
  full <- filter(d, design == full_level)
  design_short <- c("full data", "sparse annual", "overlap 10%", "overlap 25%", "alternating")
  bind_rows(
    filter(d, design != full_level) |>
      mutate(year_model = if_else(grepl("RW year", design), "RW year", "factor(year)")),
    mutate(full, year_model = "RW year"),
    mutate(full, year_model = "factor(year)")
  ) |>
    mutate(
      design = factor(
        if_else(design == full_level, "full data", sub(".*\\((.*)\\)$", "\\1", design)),
        levels = design_short
      ),
      year_model = factor(year_model, levels = c("factor(year)", "RW year"))
    )
}

plot_summary <- function(metric_names, file, width, height, split = FALSE) {
  summary_metrics <- metric_labels[metric_names]
  summary_nominal <- filter(nominal, metric %in% summary_metrics, species == first(species)) |>
    select(-species)
  rep_summary <- filter(rep_long, metric %in% summary_metrics)
  if (split) rep_summary <- split_year_model(rep_summary)
  species_long <- rep_summary |>
    summarise(
      value = if_else(first(metric) == metric_labels[["seesaw_error"]], median(value), mean(value)),
      .by = all_of(c("species", "design", "metric", if (split) "year_model"))
    ) |>
    # geom_path joins points in row order
    arrange(species, metric, design)

  ggplot(species_long, aes(value, design)) +
    geom_violin(data = rep_summary, fill = "grey90", colour = NA, scale = "width") +
    geom_vline(aes(xintercept = value), data = summary_nominal, linetype = 2, colour = "grey50") +
    stat_summary(fun = median, geom = "crossbar", fun.min = median, fun.max = median, width = 0.6) +
    geom_path(aes(colour = species, group = species)) +
    geom_point(aes(colour = species)) +
    (if (split) facet_grid(year_model ~ metric, scales = "free_x") else facet_wrap(~metric, nrow = 1, scales = "free_x")) +
    scale_colour_brewer(palette = "Paired", labels = tools::toTitleCase) +
    # first design at the top
    scale_y_discrete(limits = rev) +
    labs(x = NULL, y = NULL, colour = NULL) +
    theme(legend.position = "bottom") +
    guides(colour = guide_legend(nrow = 2))
  ggsave(here::here(file), width = width, height = height)
}
plot_summary(
  c("rmse", "coverage_50", "seesaw_error"),
  "figs/qcs-hs-sparse-vs-biennial-self-test-summary.pdf", width = 12, height = 5
)
# Without coverage, for the many designs
plot_summary(
  c("rmse", "seesaw_error"),
  "figs/qcs-hs-sparse-vs-biennial-self-test-summary-no-coverage.pdf", width = 9, height = 7,
  split = TRUE
)

# By species, all replicates
ggplot(rep_long, aes(value, design)) +
  geom_vline(aes(xintercept = value), data = nominal, linetype = 2, colour = "grey50") +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter(height = 0.15), alpha = 0.4) +
  facet_grid(species ~ metric, scales = "free_x") +
  scale_y_discrete(limits = rev) +
  labs(x = NULL, y = NULL)
ggsave(
  here::here("figs/qcs-hs-sparse-vs-biennial-self-test.pdf"),
  width = 14, height = 1 + 1.6 * length(unique(rep_long$species))
)

# Indexes of the first three replicates
compared |>
  filter(rep <= 3) |>
  ggplot(aes(cal_year, est, colour = design)) +
  geom_line(aes(y = true_est), colour = "black") +
  geom_line() +
  facet_grid(species ~ rep, scales = "free_y", labeller = labeller(rep = label_both)) +
  labs(x = "Year", y = "Biomass index", colour = NULL, caption = "Black: true index") +
  theme(legend.position = "bottom")
ggsave(
  here::here("figs/qcs-hs-sparse-vs-biennial-self-test-indexes.pdf"),
  width = 11, height = 1 + 2 * length(unique(compared$species))
)
