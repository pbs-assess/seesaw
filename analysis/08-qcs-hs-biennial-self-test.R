# Simulation self-test of QCS + HS biennial sampling with a known truth.
#
# analysis/07-qcs-hs-biennial-model-comparison.R treats the full-data index as
# the truth, but that index is itself an estimate that shares data with the
# biennial fits. Here the truth is known exactly:
#
# 1. Fit the operating model (AR1 RF, factor(year)) to the full QCS + HS data.
# 2. For each replicate, draw the spatial fields from their approximate
#    posterior and entirely new spatiotemporal fields and observation error at
#    the observed tows and, on the same field draws, the expected density on
#    the prediction grid. The true index is the area-weighted sum of that
#    expected density.
# 3. Keep the biennial subset of the simulated tows (one region per occasion,
#    alternating HS / QCS) and fit the biennial models. Also refit the
#    operating model to the full simulated data as a calibration check.
# 4. Score RMSE, bias, CI coverage, and seesaw A against the known truth.

library(dplyr)
library(ggplot2)
library(sdmTMB)

theme_set(ggsidekick::theme_sleek())

source(here::here("analysis/fit-index-models.R"))
source(here::here("analysis/metric-functions.R"))

surveyjoin::cache_data()
surveyjoin::load_sql_data()

survey_names <- c("SYN QCS", "SYN HS")
species <- "pacific cod"
n_reps <- 10L
ci_levels <- c(0.5, 0.95)

# Operating model, also refit to the full simulated data for calibration
om_model <- "AR1 RF, factor(year)"
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

family <- delta_gamma(type = "poisson-link")
# Drop random fields with SDs collapsing to zero (e.g., the encounter
# spatiotemporal field with biennial data)
control <- sdmTMBcontrol(collapse_spatial_variance = TRUE)
rw_prior <- sdmTMBpriors(sigma_V = gamma_cv(0.3, 0.5))

# Data -----------------------------------------------------------------------

# Model time as survey occasion: QCS + HS years are all odd, so on the
# calendar-year scale extra_time puts an empty year between every survey and
# AR1 rho gets stuck at its starting value of 0 (only rho^2 is identified).
prepare_species_data <- function(species) {
  dat <- surveyjoin::get_data(species, regions = "pbs") |>
    mutate(year = lubridate::year(lubridate::ymd(date))) |>
    filter(survey_name %in% survey_names, effort > 0) |>
    select(survey_name, year, lon_start, lat_start, depth_m, effort, catch_weight) |>
    tidyr::drop_na()

  # QCS was sampled alone in 2003 and 2004; keep years with both surveys
  complete_years <- dat |>
    distinct(year, survey_name) |>
    count(year) |>
    filter(n == length(survey_names)) |>
    pull(year)

  dat |>
    filter(year %in% complete_years) |>
    mutate(
      occasion = match(year, sort(unique(year))),
      primary_survey = if_else(occasion %% 2L == 1L, "SYN HS", "SYN QCS"),
      cal_year = year,
      year = occasion
    ) |>
    sdmTMB::add_utm_columns(ll_names = c("lon_start", "lat_start"), utm_crs = 3156)
}

dat <- prepare_species_data(species)
occasions <- distinct(dat, year, cal_year, primary_survey)
is_biennial <- dat$survey_name == dat$primary_survey

grid <- surveyjoin::dfo_synoptic_grid |>
  filter(survey %in% survey_names) |>
  sdmTMB::add_utm_columns(c("lon", "lat"), utm_crs = 3156) |>
  clamp_depth(dat) |>
  sdmTMB::replicate_df("year", sort(unique(dat$year)))

mesh <- make_mesh(dat, c("X", "Y"), cutoff = 10)
mesh_biennial <- make_mesh(dat[is_biennial, ], c("X", "Y"), mesh = mesh$mesh)

# Operating model ------------------------------------------------------------

# Same specifications as the corresponding models in fit_index_models(), fit
# directly to skip the factor(year) template fit. Every occasion has data, so
# extra_time only matters for documenting the time steps.
fit_model <- function(model, d, mesh) {
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
    "IID RF, RW year" = c(spatiotemporal = "iid", rw_year),
    "RW RF" = list(spatiotemporal = "rw"),
    "RW RF, RW year" = c(spatiotemporal = "rw", rw_year),
    "AR1 RF, RW year" = c(spatiotemporal = "ar1", rw_year)
  )
  args[names(extra)] <- extra
  do.call(sdmTMB, args)
}

fit_om <- fit_model(om_model, dat, mesh)
sanity(fit_om)
print(fit_om)

# Simulate -------------------------------------------------------------------

# Tows and grid cells in one newdata so both share the same field draw.
# eta_i excludes the offset, so grid rows give expected density per unit area;
# the simulated tow catches include the effort offset.
sim_newdata <- bind_rows(
  select(dat, X, Y, year, effort),
  mutate(select(grid, X, Y, year), effort = 1)
)
is_tow <- seq_len(nrow(sim_newdata)) <= nrow(dat)

# Fixed effects at their MLEs. Each replicate takes a new draw of the random
# effects from their approximate posterior (mle_mvn_samples = "multiple");
# those in simulate_re are then replaced by entirely new draws from the GMRF
# prior. So the spatial fields are posterior draws (unless resample_spatial)
# and the spatiotemporal fields are new.
sim_reports <- simulate(
  fit_om,
  nsim = n_reps,
  type = "mle-mvn",
  mle_mvn_samples = "multiple",
  simulate_re = if (resample_spatial) c("spatial", "spatiotemporal") else "spatiotemporal",
  newdata = sim_newdata,
  offset = log(sim_newdata$effort),
  return_tmb_report = TRUE,
  seed = 42,
  silent = TRUE
)

sims <- purrr::imap(sim_reports, \(r, rep) {
  true_index <- grid |>
    mutate(density = exp(r$eta_i[!is_tow, 1L] + r$eta_i[!is_tow, 2L]) * area) |>
    group_by(year) |>
    summarise(true_est = sum(density)) |>
    mutate(rep = rep, .before = 1L)

  list(
    catch = r$y_i[is_tow, 1L] * r$y_i[is_tow, 2L],
    true_index = true_index
  )
})
true_indexes <- purrr::map_dfr(sims, "true_index")

# Fit biennial models ---------------------------------------------------------

fit_rep <- function(rep, catch) {
  RhpcBLASctl::blas_set_num_threads(1L)
  RhpcBLASctl::omp_set_num_threads(1L)

  d <- dat
  d$catch_weight <- catch

  index_ok <- function(model, d, mesh) {
    fit <- fit_ok(fit_model(model, d, mesh))
    if (is.null(fit)) tibble::tibble() else get_index_ok(fit, grid)
  }

  bind_rows(
    full = index_ok(om_model, d, mesh) |> mutate(model = om_model),
    biennial = purrr::map(biennial_models, \(m) index_ok(m, d[is_biennial, ], mesh_biennial)) |>
      setNames(biennial_models) |>
      bind_rows(.id = "model"),
    .id = "design"
  ) |>
    mutate(
      model = if_else(design == "full", paste0(model, " (full data)"), model),
      rep = rep,
      .before = 1L
    )
}

fits_file <- here::here("data-generated/qcs-hs-biennial-self-test.rds")
if (!file.exists(fits_file)) {
  future::plan(future::multisession, workers = min(n_reps, future::availableCores() / 2))
  indexes <- furrr::future_map2_dfr(
    seq_len(n_reps), purrr::map(sims, "catch"), fit_rep,
    # One replicate per future so a slow fit doesn't hold up a whole chunk
    .options = furrr::furrr_options(seed = TRUE, scheduling = Inf)
  )
  future::plan(future::sequential)
  saveRDS(list(indexes = indexes, true_indexes = true_indexes), fits_file)
}
self_test <- readRDS(fits_file)
indexes <- self_test$indexes
true_indexes <- self_test$true_indexes

# Performance metrics ----------------------------------------------------------

compared <- indexes |>
  left_join(true_indexes, by = c("rep", "year")) |>
  left_join(occasions, by = "year") |>
  mutate(log_error = log_est - log(true_est))

for (level in ci_levels) {
  q <- qnorm(1 - (1 - level) / 2)
  compared[[paste0("covered_", level * 100)]] <- abs(compared$log_error) <= q * compared$se
}

# Phase follows which region was sampled, not calendar-year parity
seesaw_A <- function(est, primary_survey, year) {
  as_tibble(t(period2_metric(
    est,
    year = year,
    phase = if_else(primary_survey == "SYN HS", 0.5, -0.5)
  )))
}

metrics_rep <- compared |>
  arrange(rep, model, year) |>
  group_by(rep, model) |>
  summarise(
    rmse = sqrt(mean(log_error^2)),
    bias = mean(log_error),
    coverage_50 = mean(covered_50),
    coverage_95 = mean(covered_95),
    A = seesaw_A(est, primary_survey, year)$A,
    .groups = "drop"
  )

A_truth <- true_indexes |>
  left_join(occasions, by = "year") |>
  arrange(rep, year) |>
  group_by(rep) |>
  summarise(A = seesaw_A(true_est, primary_survey, year)$A)

metrics <- metrics_rep |>
  group_by(model) |>
  summarise(
    n_reps = n(),
    rmse = mean(rmse),
    bias = mean(bias),
    coverage_50 = mean(coverage_50),
    coverage_95 = mean(coverage_95),
    A_median = median(A),
    .groups = "drop"
  ) |>
  arrange(rmse)

saveRDS(
  list(compared = compared, metrics_rep = metrics_rep, metrics = metrics, A_truth = A_truth),
  here::here("data-generated/qcs-hs-biennial-self-test-metrics.rds")
)

metrics |>
  mutate(across(where(is.numeric), \(x) round(x, 2))) |>
  print(width = Inf)
cat("Median truth A:", round(median(A_truth$A), 1), "\n")

# Plots ----------------------------------------------------------------------

metric_levels <- c("rmse", "bias", "coverage_50", "coverage_95", "A")
metric_labels <- c("RMSE (log)", "Bias (log)", "50% CI coverage", "95% CI coverage", "Seesaw A (%)")
nominal <- tibble::tibble(
  metric = factor(c("bias", "coverage_50", "coverage_95", "A"), levels = metric_levels, labels = metric_labels),
  value = c(0, 0.5, 0.95, median(A_truth$A))
)

metrics_rep |>
  tidyr::pivot_longer(all_of(metric_levels), names_to = "metric") |>
  mutate(
    model = factor(model, levels = rev(metrics$model)),
    metric = factor(metric, levels = metric_levels, labels = metric_labels)
  ) |>
  ggplot(aes(value, model)) +
  geom_vline(aes(xintercept = value), data = nominal, linetype = 2, colour = "grey50") +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter(height = 0.15), alpha = 0.4) +
  facet_wrap(~metric, scales = "free_x", nrow = 1) +
  labs(
    x = NULL, y = NULL,
    caption = paste0(
      "Biennial QCS/HS self-test, ", species, ", ", n_reps, " replicates, OM: ", om_model,
      if (resample_spatial) " (spatial fields new from prior). " else " (spatial fields posterior draws). ",
      "Dashed lines: nominal values (A: median of the true index)."
    )
  )
ggsave(here::here("figs/qcs-hs-biennial-self-test.pdf"), width = 13, height = 3.5)

compared |>
  filter(rep <= 6) |>
  ggplot(aes(cal_year, est, colour = model)) +
  geom_line(aes(y = true_est), colour = "black") +
  geom_line() +
  facet_wrap(~rep, scales = "free_y", labeller = label_both) +
  labs(x = "Year", y = "Biomass index", colour = "Biennial model", caption = "Black: true index")
ggsave(here::here("figs/qcs-hs-biennial-self-test-indexes.pdf"), width = 10, height = 6)
