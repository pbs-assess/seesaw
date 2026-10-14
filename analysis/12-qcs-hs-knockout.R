# Contrast knockout for the QCS + HS biennial self-test: is the HS / QCS
# density contrast what causes the spurious seesaw?
#
# As analysis/08-qcs-hs-biennial-self-test.R, except the operating model
# includes an HS fixed effect so the mean contrast between the regions is one
# coefficient (per delta component) rather than sitting in omega. Each
# replicate is simulated twice from the same parameter draws and seed:
#
#   k = 1: the fitted contrast (as observed)
#   k = 0: the HS coefficients set to 0 (contrast knocked out)
#
# Everything else (spatial field posterior draws, new spatiotemporal fields,
# and the observation draws in QCS) is identical between the pair, so any
# difference in the excess seesaw is caused by the contrast. Only the models
# needed for the comparison are fit.
#
# The pipeline is in analysis/self-test-funcs.R. Everything is cached per
# replicate, so n_reps can be raised later and only the new replicates are
# run.

source(here::here("analysis/self-test-funcs.R"))
theme_set(ggsidekick::theme_sleek())

surveyjoin::cache_data()
surveyjoin::load_sql_data()

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
knockout_k <- c(0, 1)
# As in 08: with 11 occasions, 8 gives 4 windows
seesaw_window <- 8L
n_workers <- max(1L, min(60L, future::availableCores() - 2L))

# Operating model (plus the HS term), also refit without it to the full
# simulated data for calibration. IID spatiotemporal fields, as in 08: with
# AR1, rho can go to ~1, so each replicate's new epsilon_st is a near-static
# field that sets a random HS / QCS contrast, confounding the knockout.
om_model <- "IID RF, factor(year)"
biennial_models <- c("IID RF, factor(year)", "IID RF, RW year")

# Run ------------------------------------------------------------------------

qcs_hs_base_grid <- qcs_hs_grid()

setup_species <- function(species, dat) {
  setup <- qcs_hs_setup(species, dat, qcs_hs_base_grid)
  om <- fit_operating_model(om_model, dat, setup$mesh, formula = catch_weight ~ 0 + factor(year) + hs)
  if (isTRUE(om$om_failed)) return(c(list(species = species), om))

  # Positions of the HS coefficient for each delta component in the joint
  # parameter vector (b_j: component 1, b_j2: component 2)
  pars <- om$fit$tmb_obj$env$last.par.best
  i_hs <- which(colnames(om$fit$tmb_data$X_ij[[1]]) == "hs")
  idx_hs <- c(which(names(pars) == "b_j")[i_hs], which(names(pars) == "b_j2")[i_hs])

  c(setup, list(fit_om = om$fit, ran_pars = om$ran_pars, idx_hs = idx_hs, hs_coef = pars[idx_hs]))
}

# One parameter draw per replicate, shared by both k; the HS coefficients are
# scaled by k. The same seed gives both the same new spatiotemporal fields.
simulate_rep <- function(setup, sampler, rep) {
  pars <- draw_pars(sampler, rep)
  reports <- lapply(setNames(knockout_k, knockout_k), \(k) {
    pars_k <- pars
    pars_k[setup$idx_hs, ] <- k * pars_k[setup$idx_hs, ]
    simulate_pars(setup, pars_k, rep)
  })
  list(
    catches = lapply(reports, \(r) rep_catch(setup, r)),
    truth = purrr::map_dfr(reports, \(r) rep_truth(setup, r), .id = "arm")
  )
}

results <- run_self_test(
  dir = here::here("data-generated/qcs-hs-knockout"),
  species = species_to_fit,
  n_reps = n_reps,
  prepare_data = qcs_hs_data,
  setup_species = setup_species,
  simulate_rep = simulate_rep,
  models = tibble::tibble(
    design = c("full", rep("alternating", length(biennial_models))),
    model = c(om_model, biennial_models)
  ),
  n_workers = n_workers
)

if (FALSE) {
  s <- setup_species("spotted ratfish", qcs_hs_data("spotted ratfish"))
  s$hs_coef
  sim <- simulate_rep(s, sdmTMB:::.posterior_re_sampler(s$fit_om), 1L)
  sim$truth |> group_by(arm) |> summarise(contrast = mean(contrast))
}

# Metrics ----------------------------------------------------------------------

# Arms are the k values
compared <- mutate(compare_to_truth(results), k = as.numeric(arm))

metrics_rep <- compared |>
  arrange(species, k, rep, model, year) |>
  group_by(species, k, rep, model) |>
  summarise(
    rmse = sqrt(mean(log_error^2)),
    coverage_50 = mean(covered_50),
    coverage_95 = mean(covered_95),
    A_max = seesaw_A_max(est, phase, year, seesaw_window),
    .groups = "drop"
  )

A_truth <- results$true_indexes |>
  mutate(k = as.numeric(arm)) |>
  left_join(results$occasions, by = c("species", "year")) |>
  arrange(species, k, rep, year) |>
  group_by(species, k, rep) |>
  summarise(
    A_max_true = seesaw_A_max(true_est, phase, year, seesaw_window),
    contrast = mean(contrast),
    .groups = "drop"
  )

metrics_rep <- metrics_rep |>
  left_join(A_truth, by = c("species", "k", "rep")) |>
  mutate(A_max_excess = A_max - A_max_true)

metrics <- metrics_rep |>
  group_by(species, model, k) |>
  summarise(
    n_reps = n(),
    contrast = mean(contrast),
    rmse = mean(rmse),
    coverage_50 = mean(coverage_50),
    A_max_excess = median(A_max_excess),
    .groups = "drop"
  )

# Paired effect of the contrast on the excess seesaw (k = 1 minus k = 0, same
# replicate)
paired <- metrics_rep |>
  select(species, model, rep, k, A_max_excess) |>
  tidyr::pivot_wider(names_from = k, values_from = A_max_excess, names_prefix = "k") |>
  mutate(effect = k1 - k0) |>
  group_by(species, model) |>
  summarise(
    n_pairs = sum(!is.na(effect)),
    effect_median = median(effect, na.rm = TRUE),
    prop_positive = mean(effect > 0, na.rm = TRUE),
    .groups = "drop"
  )

saveRDS(
  list(
    compared = compared, metrics_rep = metrics_rep, metrics = metrics, paired = paired,
    failure_rates = results$failure_rates
  ),
  here::here("data-generated/qcs-hs-knockout-metrics.rds")
)

metrics |>
  mutate(across(where(is.numeric), \(x) round(x, 2))) |>
  arrange(model, species, k) |>
  print(n = Inf, width = Inf)
print(paired, n = Inf)

# Plots ----------------------------------------------------------------------

k_labels <- c("0" = "Contrast removed", "1" = "Observed contrast")

metrics_rep |>
  mutate(k = factor(k, levels = c(1, 0), labels = k_labels[c("1", "0")])) |>
  ggplot(aes(k, A_max_excess)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey60") +
  geom_line(aes(group = rep), colour = "grey80") +
  geom_point(colour = "grey60") +
  stat_summary(aes(group = 1), fun = median, geom = "line", colour = "black") +
  stat_summary(fun = median, geom = "point", colour = "black") +
  facet_grid(model ~ species, labeller = labeller(species = \(x) tools::toTitleCase(x))) +
  labs(
    x = NULL,
    y = paste0("Max A - true max A (%), ", seesaw_window, "-survey windows"),
    caption = paste0(
      "QCS/HS biennial self-test with the HS / QCS contrast knocked out. Grey: paired replicates ",
      "(same fields and seed); black: median."
    )
  ) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(here::here("figs/qcs-hs-knockout.pdf"), width = 13, height = 7)
