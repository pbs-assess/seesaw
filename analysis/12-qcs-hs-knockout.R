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

library(dplyr)
library(ggplot2)
library(sdmTMB)

theme_set(ggsidekick::theme_sleek())

source(here::here("analysis/fit-index-models.R"))
source(here::here("analysis/metric-functions.R"))

surveyjoin::cache_data()
surveyjoin::load_sql_data()

survey_names <- c("SYN QCS", "SYN HS")
species_to_fit <- c(
  "pacific cod",
  "spotted ratfish",
  "pacific ocean perch",
  "pacific halibut",
  "pacific spiny dogfish",
  "lingcod",
  "shortspine thornyhead"
)
n_reps <- 20L
knockout_k <- c(0, 1)
ci_levels <- c(0.5, 0.95)
# As in 08: with 11 occasions, 8 gives 4 windows
seesaw_window <- 8L

# Operating model (plus the HS term), also refit without it to the full
# simulated data for calibration
om_model <- "AR1 RF, factor(year)"
biennial_models <- c("IID RF, factor(year)", "IID RF, RW year")

family <- delta_gamma(type = "poisson-link")
control <- sdmTMBcontrol(collapse_spatial_variance = TRUE)
rw_prior <- sdmTMBpriors(sigma_V = lognormal_prior(log(0.3), 0.5))

# Data -----------------------------------------------------------------------

# Model time as survey occasion (see 08)
prepare_species_data <- function(species) {
  dat <- surveyjoin::get_data(species, regions = "pbs") |>
    mutate(year = lubridate::year(lubridate::ymd(date))) |>
    filter(survey_name %in% survey_names, effort > 0) |>
    select(survey_name, year, lon_start, lat_start, depth_m, effort, catch_weight) |>
    tidyr::drop_na()

  complete_years <- dat |>
    distinct(year, survey_name) |>
    count(year) |>
    filter(n == length(survey_names)) |>
    pull(year)

  dat |>
    filter(year %in% complete_years) |>
    mutate(
      hs = as.integer(survey_name == "SYN HS"),
      occasion = match(year, sort(unique(year))),
      primary_survey = if_else(occasion %% 2L == 1L, "SYN HS", "SYN QCS"),
      cal_year = year,
      year = occasion
    ) |>
    sdmTMB::add_utm_columns(ll_names = c("lon_start", "lat_start"), utm_crs = 3156)
}

base_grid <- surveyjoin::dfo_synoptic_grid |>
  filter(survey %in% survey_names) |>
  mutate(hs = as.integer(survey == "SYN HS")) |>
  sdmTMB::add_utm_columns(c("lon", "lat"), utm_crs = 3156)

# Models ---------------------------------------------------------------------

# As in 08
fit_model <- function(model, d, mesh, formula = NULL) {
  args <- list(
    formula = catch_weight ~ 1,
    data = d,
    mesh = mesh,
    offset = log(d$effort),
    family = family,
    time = "year",
    spatial = "on",
    share_range = TRUE,
    anisotropy = TRUE,
    extra_time = seq(min(d$year), max(d$year)),
    control = control,
    silent = TRUE
  )
  factor_year <- list(formula = catch_weight ~ 0 + factor(year), extra_time = NULL)
  rw_year <- list(time_varying = ~1, time_varying_type = "rw0", priors = rw_prior)
  extra <- switch(model,
    "IID RF, factor(year)" = c(spatiotemporal = "iid", factor_year),
    "AR1 RF, factor(year)" = c(spatiotemporal = "ar1", factor_year),
    "IID RF, RW year" = c(spatiotemporal = "iid", rw_year)
  )
  args[names(extra)] <- extra
  if (!is.null(formula)) args$formula <- formula
  do.call(sdmTMB, args)
}

# Operating model and paired simulation ----------------------------------------

setup_species <- function(species) {
  RhpcBLASctl::blas_set_num_threads(1L)
  RhpcBLASctl::omp_set_num_threads(1L)

  dat <- prepare_species_data(species)
  is_biennial <- dat$survey_name == dat$primary_survey
  grid <- sdmTMB::replicate_df(base_grid, "year", sort(unique(dat$year)))
  mesh <- make_mesh(dat, c("X", "Y"), cutoff = 10)
  mesh_biennial <- make_mesh(dat[is_biennial, ], c("X", "Y"), mesh = mesh$mesh)

  fit_om <- fit_ok(fit_model(om_model, dat, mesh, formula = catch_weight ~ 0 + factor(year) + hs))
  if (is.null(fit_om)) {
    cli::cli_warn("Operating model did not pass sanity checks for {species}")
    return(NULL)
  }

  # Positions of the HS coefficient for each delta component in the joint
  # parameter vector (b_j: component 1, b_j2: component 2)
  pars <- fit_om$tmb_obj$env$last.par.best
  i_hs <- which(colnames(fit_om$tmb_data$X_ij[[1]]) == "hs")
  idx_hs <- c(which(names(pars) == "b_j")[i_hs], which(names(pars) == "b_j2")[i_hs])

  # One set of parameter draws shared by both k: fixed effects at their MLEs,
  # random effects from their approximate posterior (as mle_mvn_samples =
  # "multiple"); the spatiotemporal fields are then replaced by new draws
  set.seed(42)
  par_draws <- sdmTMB:::.posterior_re_sampler(fit_om)(n_reps)

  sim_newdata <- bind_rows(
    select(dat, X, Y, year, hs, effort),
    mutate(select(grid, X, Y, year, hs), effort = 1)
  )
  is_tow <- seq_len(nrow(sim_newdata)) <= nrow(dat)
  grid_hs <- grid$hs == 1L

  sims <- purrr::map(setNames(knockout_k, knockout_k), \(k) {
    p <- par_draws
    p[idx_hs, ] <- k * p[idx_hs, ]
    reports <- simulate(
      fit_om,
      nsim = n_reps,
      mcmc_samples = p,
      simulate_re = "spatiotemporal",
      newdata = sim_newdata,
      offset = log(sim_newdata$effort),
      return_tmb_report = TRUE,
      seed = 42,
      silent = TRUE
    )
    truth <- purrr::imap_dfr(reports, \(r, rep) {
      eta <- r$eta_i[!is_tow, 1L] + r$eta_i[!is_tow, 2L]
      grid |>
        mutate(density = exp(eta) * area, eta = eta) |>
        group_by(year) |>
        summarise(
          true_est = sum(density),
          # realized log density contrast (HS minus QCS, area-weighted)
          contrast = weighted.mean(eta[hs == 1L], area[hs == 1L]) -
            weighted.mean(eta[hs == 0L], area[hs == 0L])
        ) |>
        mutate(rep = rep, .before = 1L)
    })
    list(
      catches = lapply(reports, \(r) r$y_i[is_tow, 1L] * r$y_i[is_tow, 2L]),
      truth = truth
    )
  })

  list(
    species = species,
    dat = dat,
    is_biennial = is_biennial,
    grid = grid,
    mesh = mesh,
    mesh_biennial = mesh_biennial,
    hs_coef = pars[idx_hs],
    occasions = distinct(dat, year, cal_year, primary_survey),
    catches = purrr::map(sims, "catches"),
    true_indexes = purrr::map_dfr(sims, "truth", .id = "k") |> mutate(k = as.numeric(k))
  )
}

# Fit replicates --------------------------------------------------------------

fit_rep <- function(setup, k, rep) {
  RhpcBLASctl::blas_set_num_threads(1L)
  RhpcBLASctl::omp_set_num_threads(1L)

  d <- setup$dat
  d$catch_weight <- setup$catches[[as.character(k)]][[rep]]

  index_ok <- function(model, d, mesh) {
    fit <- fit_ok(fit_model(model, d, mesh))
    if (is.null(fit)) tibble::tibble() else get_index_ok(fit, setup$grid)
  }

  bind_rows(
    full = index_ok(om_model, d, setup$mesh) |> mutate(model = om_model),
    biennial = purrr::map(
      biennial_models,
      \(m) index_ok(m, d[setup$is_biennial, ], setup$mesh_biennial)
    ) |>
      setNames(biennial_models) |>
      bind_rows(.id = "model"),
    .id = "design"
  ) |>
    mutate(
      model = if_else(design == "full", paste0(model, " (full data)"), model),
      species = setup$species,
      k = k,
      rep = rep,
      .before = 1L
    )
}

if (FALSE) {
  tictoc::tic()
  s <- setup_species("spotted ratfish")
  tictoc::toc()
  s$hs_coef
  s$true_indexes |> group_by(k) |> summarise(contrast = mean(contrast))
  x <- fit_rep(s, 0, 1L)
}

# One cache file per species; delete a file to refit that species
fits_dir <- here::here("data-generated/qcs-hs-knockout")
dir.create(fits_dir, showWarnings = FALSE)
fits_file <- \(species) file.path(fits_dir, paste0(gsub(" ", "-", species), ".rds"))

to_fit <- species_to_fit[!file.exists(fits_file(species_to_fit))]
if (length(to_fit) > 0L) {
  future::plan(future::multicore, workers = max(1L, min(25L, future::availableCores() - 2L)))
  setups <- furrr::future_map(
    to_fit, setup_species,
    .options = furrr::furrr_options(seed = TRUE, scheduling = Inf)
  ) |>
    purrr::compact()

  # One species x k x replicate per future
  tasks <- tidyr::expand_grid(i = seq_along(setups), k = knockout_k, rep = seq_len(n_reps))
  rep_indexes <- furrr::future_pmap(
    tasks, \(i, k, rep) fit_rep(setups[[i]], k, rep),
    .options = furrr::furrr_options(seed = TRUE, scheduling = Inf)
  )
  future::plan(future::sequential)

  for (i in seq_along(setups)) {
    saveRDS(
      list(
        indexes = bind_rows(rep_indexes[tasks$i == i]),
        true_indexes = mutate(setups[[i]]$true_indexes, species = setups[[i]]$species, .before = 1L),
        occasions = mutate(setups[[i]]$occasions, species = setups[[i]]$species, .before = 1L),
        hs_coef = setups[[i]]$hs_coef
      ),
      fits_file(setups[[i]]$species)
    )
  }
}

knockout <- lapply(fits_file(species_to_fit), \(f) if (file.exists(f)) readRDS(f))
indexes <- purrr::map_dfr(knockout, "indexes")
true_indexes <- purrr::map_dfr(knockout, "true_indexes")
occasions <- purrr::map_dfr(knockout, "occasions")

# Metrics ----------------------------------------------------------------------

compared <- indexes |>
  left_join(true_indexes, by = c("species", "k", "rep", "year")) |>
  left_join(occasions, by = c("species", "year")) |>
  mutate(log_error = log_est - log(true_est))

for (level in ci_levels) {
  q <- qnorm(1 - (1 - level) / 2)
  compared[[paste0("covered_", level * 100)]] <- abs(compared$log_error) <= q * compared$se
}

phase_of <- \(primary_survey) if_else(primary_survey == "SYN HS", 0.5, -0.5)
seesaw_A_max <- function(est, primary_survey, year) {
  max(moving_window(est, year = year, window = seesaw_window, phase = phase_of(primary_survey))$A)
}

metrics_rep <- compared |>
  arrange(species, k, rep, model, year) |>
  group_by(species, k, rep, model) |>
  summarise(
    rmse = sqrt(mean(log_error^2)),
    coverage_50 = mean(covered_50),
    coverage_95 = mean(covered_95),
    A_max = seesaw_A_max(est, primary_survey, year),
    .groups = "drop"
  )

A_truth <- true_indexes |>
  left_join(occasions, by = c("species", "year")) |>
  arrange(species, k, rep, year) |>
  group_by(species, k, rep) |>
  summarise(
    A_max_true = seesaw_A_max(true_est, primary_survey, year),
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
  list(compared = compared, metrics_rep = metrics_rep, metrics = metrics, paired = paired),
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
