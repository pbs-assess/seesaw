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
#      fit with IID RF + RW year (with the sigma_V prior)
#    - sparse annual: a random subsample of the tows of the whole domain
#      (both regions) in every occasion, a different draw for each replicate,
#      fit with IID RF + factor(year)
# 4. Score RMSE and 50% (and 95%) CI coverage of the bias-corrected index
#    against the truth of each replicate.
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
n_workers <- max(1L, min(60L, future::availableCores() - 2L))

# Operating model; IID fields for the reasons given in analysis/08
om_model <- "IID RF, factor(year)"
sparse_model <- "IID RF, factor(year)"
alternating_model <- "IID RF, RW year"

designs <- c(
  sparse = paste0(sparse_model, " (sparse annual)"),
  alternating = paste0(alternating_model, " (alternating)")
)

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
    design = c("sparse", "alternating"),
    model = c(sparse_model, alternating_model)
  ),
  n_workers = n_workers
)

# Metrics --------------------------------------------------------------------

compared <- compare_to_truth(results) |>
  mutate(design = factor(design, levels = names(designs), labels = designs))

metrics_rep <- compared |>
  group_by(species, rep, design) |>
  summarise(
    rmse = sqrt(mean(log_error^2)),
    bias = mean(log_error),
    coverage_50 = mean(covered_50),
    coverage_95 = mean(covered_95),
    .groups = "drop"
  )

metrics <- metrics_rep |>
  group_by(species, design) |>
  summarise(
    n_reps = n(),
    across(c(rmse, bias, coverage_50, coverage_95), mean),
    .groups = "drop"
  )

saveRDS(
  list(compared = compared, metrics_rep = metrics_rep, metrics = metrics, failure_rates = results$failure_rates),
  here::here("data-generated/qcs-hs-sparse-vs-biennial-self-test-metrics.rds")
)

metrics |>
  mutate(across(where(is.numeric), \(x) round(x, 2))) |>
  print(n = Inf, width = Inf)

# Across species and replicates; the paired RMSE ratio is within a replicate
# (same simulated population), < 1 when the sparse annual design is better
metrics_rep |>
  select(species, rep, design, rmse) |>
  tidyr::pivot_wider(names_from = design, values_from = rmse) |>
  summarise(
    rmse_ratio_sparse_over_alternating = exp(mean(log(.data[[designs[["sparse"]]]] / .data[[designs[["alternating"]]]]))),
    prop_reps_sparse_better = mean(.data[[designs[["sparse"]]]] < .data[[designs[["alternating"]]]]),
    .by = species
  ) |>
  mutate(across(where(is.numeric), \(x) round(x, 2))) |>
  print(n = Inf, width = Inf)

# Plots ----------------------------------------------------------------------

metric_labels <- c(rmse = "RMSE (log)", coverage_50 = "50% CI coverage")
nominal <- tibble::tibble(
  metric = factor(metric_labels[["coverage_50"]], levels = metric_labels),
  value = 0.5
)

rep_long <- metrics_rep |>
  tidyr::pivot_longer(all_of(names(metric_labels)), names_to = "metric") |>
  mutate(metric = factor(metric, levels = names(metric_labels), labels = metric_labels))

# One point per species (mean across replicates) over the violins of all
# replicates pooled across species
species_long <- rep_long |>
  summarise(value = mean(value), .by = c(species, design, metric))

ggplot(species_long, aes(value, design)) +
  geom_violin(data = rep_long, fill = "grey90", colour = NA, scale = "width") +
  geom_vline(aes(xintercept = value), data = nominal, linetype = 2, colour = "grey50") +
  stat_summary(fun = median, geom = "crossbar", fun.min = median, fun.max = median, width = 0.6) +
  geom_point(aes(colour = species), position = position_jitter(height = 0.15, seed = 1)) +
  facet_wrap(~metric, nrow = 1, scales = "free_x") +
  scale_colour_brewer(palette = "Paired", labels = tools::toTitleCase) +
  labs(
    x = NULL, y = NULL, colour = NULL,
    caption = paste0(
      "QCS + HS self-test, ", n_reps, " replicates per species, OM: ", om_model, ". ",
      "Same tows per occasion in both designs. Dashed: nominal coverage."
    )
  ) +
  theme(legend.position = "bottom") +
  guides(colour = guide_legend(nrow = 2))
ggsave(here::here("figs/qcs-hs-sparse-vs-biennial-self-test-summary.pdf"), width = 8, height = 3.8)

# By species, all replicates
ggplot(rep_long, aes(value, design)) +
  geom_vline(aes(xintercept = value), data = nominal, linetype = 2, colour = "grey50") +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter(height = 0.15), alpha = 0.4) +
  facet_grid(species ~ metric, scales = "free_x") +
  labs(x = NULL, y = NULL)
ggsave(
  here::here("figs/qcs-hs-sparse-vs-biennial-self-test.pdf"),
  width = 9, height = 1 + 1.6 * length(unique(rep_long$species))
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
