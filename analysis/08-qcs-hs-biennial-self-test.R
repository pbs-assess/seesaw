# Simulation self-test of QCS + HS biennial sampling with a known truth.
#
# analysis/07-qcs-hs-biennial-model-comparison.R treats the full-data index as
# the truth, but that index is itself an estimate that shares data with the
# biennial fits. Here the truth is known exactly:
#
# 1. Fit the operating model (IID RF, factor(year)) to the full QCS + HS data.
# 2. For each replicate, draw the spatial fields from their approximate
#    posterior and entirely new spatiotemporal fields and observation error at
#    the observed tows and, on the same field draws, the expected density on
#    the prediction grid. The true index is the area-weighted sum of that
#    expected density.
# 3. Keep the biennial subset of the simulated tows (one region per occasion,
#    alternating HS / QCS) and fit the biennial models. Also refit the
#    operating model to the full simulated data as a calibration check.
# 4. Score RMSE, bias, CI coverage, seesaw A, and the signed seesaw error
#    against the known truth.
#
# The pipeline and metrics are in analysis/self-test-funcs.R. Everything is
# cached per replicate, so n_reps can be raised later and only the new
# replicates are run.

source(here::here("analysis/self-test-funcs.R"))
theme_set(ggsidekick::theme_sleek())

surveyjoin::cache_data()
surveyjoin::load_sql_data()

# Strong to moderate seesaws under biennial QCS / HS sampling across a range of
# taxa, plus arrowtooth flounder as a control with almost none; see
# analysis/qcs-hs-biennial-screen.R
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
# Moving window (in survey occasions) for max A; as in
# analysis/05-qcs-hs-experimental-biennial.R. With 11 occasions, 8 gives 4
# windows (10, as in the annual-survey scripts, would give only 2).
seesaw_window <- 8L
n_workers <- max(1L, min(60L, future::availableCores() - 2L))

# Operating model, also refit to the full simulated data for calibration. IID
# spatiotemporal fields: with AR1, rho can go to ~1, so each replicate's new
# epsilon_st is in effect a second, persistent spatial field drawn from the
# prior, swamping the posterior omega_s. With IID, persistent structure stays
# in omega_s and the new epsilon_st are year-to-year noise.
om_model <- "IID RF, factor(year)"
# FALSE: spatial fields are a new draw from their approximate posterior each
# replicate (keeps the fitted gradient across QCS / HS, with uncertainty);
# TRUE: spatial fields are entirely new draws from the GMRF prior
resample_spatial <- FALSE

biennial_models <- c(
  "IID RF, factor(year)",
  "IID RF, RW year",
  "RW RF",
  "RW RF, RW year",
  "AR1 RF, RW year"
)

# Run ------------------------------------------------------------------------

qcs_hs_base_grid <- qcs_hs_grid()

setup_species <- function(species, dat) {
  setup <- qcs_hs_setup(species, dat, qcs_hs_base_grid)
  om <- fit_operating_model(om_model, dat, setup$mesh)
  if (isTRUE(om$om_failed)) return(c(list(species = species), om))
  c(setup, list(fit_om = om$fit, ran_pars = om$ran_pars))
}

results <- run_self_test(
  dir = here::here("data-generated/qcs-hs-biennial-self-test"),
  species = species_to_fit,
  n_reps = n_reps,
  prepare_data = qcs_hs_data,
  setup_species = setup_species,
  simulate_rep = simulate_one_arm(
    if (resample_spatial) c("spatial", "spatiotemporal") else "spatiotemporal"
  ),
  models = tibble::tibble(
    design = c("full", rep("alternating", length(biennial_models))),
    model = c(om_model, biennial_models)
  ),
  n_workers = n_workers
)

# Metrics and plots ------------------------------------------------------------

m <- self_test_metrics(compare_to_truth(results), results, seesaw_window)
saveRDS(m, here::here("data-generated/qcs-hs-biennial-self-test-metrics.rds"))

m$metrics |>
  mutate(across(where(is.numeric), \(x) round(x, 2))) |>
  print(n = Inf, width = Inf)

m$A_truth |>
  group_by(species) |>
  summarise(A_median = median(A), A_max_median = median(A_max))

plot_self_test(
  m,
  full_model = paste0(om_model, " (full data)"),
  window = seesaw_window,
  fig_prefix = here::here("figs/qcs-hs-biennial-self-test"),
  caption = paste0(
    "Biennial QCS/HS self-test, ", n_reps, " replicates per species, OM: ", om_model,
    if (resample_spatial) " (spatial fields new from prior). " else " (spatial fields posterior draws). ",
    "Dashed lines: nominal values (A, max A: median of the true index). ",
    "Seesaw error: estimated minus true phase effect, positive when the index reads high in occasions the richer region is sampled."
  )
)
