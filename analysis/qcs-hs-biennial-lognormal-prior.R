# Lognormal vs. gamma priors on sigma_V for the IID RF, RW year model, Pacific
# cod under pretend biennial QCS / HS sampling.
#
# The lognormal prior is a custom prior (normal on log(sigma_V)); see
# https://sdmtmb.github.io/sdmTMB/articles/custom-priors.html. No Jacobian is
# needed because the fits are not Bayesian. The gamma fits and the full-data
# reference come from analysis/qcs-hs-biennial-prior-sensitivity.R (run that
# first). Lognormal medians match the gamma means and sdlog values match the
# gamma CVs (a lognormal CV is ~sdlog for small sdlog). Time is survey
# occasion, so sigma_V is per occasion (2 calendar years).

library(dplyr)
library(ggplot2)
library(sdmTMB)

theme_set(ggsidekick::theme_sleek())

source(here::here("analysis/fit-index-models.R"))
source(here::here("analysis/metric-functions.R"))

gamma_file <- here::here("data-generated/qcs-hs-biennial-prior-sensitivity.rds")
stopifnot(file.exists(gamma_file))

surveyjoin::cache_data()
surveyjoin::load_sql_data()

species <- "pacific cod"
survey_names <- c("SYN QCS", "SYN HS")
ci_level <- 0.5
seesaw_window <- 8L

prior_grid <- tidyr::crossing(mean = c(0.1, 0.3, 0.5), spread = c(0.1, 0.3, 0.5, 0.8))

# Data (as in analysis/qcs-hs-biennial-prior-sensitivity.R) ------------------

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

dat <- dat |>
  filter(year %in% complete_years) |>
  mutate(
    occasion = match(year, sort(unique(year))),
    primary_survey = if_else(occasion %% 2L == 1L, "SYN HS", "SYN QCS"),
    cal_year = year,
    year = occasion
  ) |>
  sdmTMB::add_utm_columns(ll_names = c("lon_start", "lat_start"), utm_crs = 3156)
dat_biennial <- filter(dat, survey_name == primary_survey)
occasions <- distinct(dat, year, cal_year, primary_survey)

grid <- surveyjoin::dfo_synoptic_grid |>
  filter(survey %in% survey_names) |>
  sdmTMB::add_utm_columns(c("lon", "lat"), utm_crs = 3156) |>
  clamp_depth(dat) |>
  sdmTMB::replicate_df("year", sort(unique(dat$year)))

mesh <- make_mesh(dat, c("X", "Y"), cutoff = 10)
mesh_biennial <- make_mesh(dat_biennial, c("X", "Y"), mesh = mesh$mesh)

# Fits -----------------------------------------------------------------------

# Same prior on both delta components
lognormal_prior <- function(median, sdlog) {
  sdmTMBpriors(custom = function(par, theta) {
    RTMB::dnorm(log(theta$sigma_V[1, ]), log(median), sdlog, log = TRUE)
  })
}

fit_rw_year <- function(priors) {
  sdmTMB(
    catch_weight ~ 1,
    data = dat_biennial,
    mesh = mesh_biennial,
    offset = log(dat_biennial$effort),
    family = delta_gamma(type = "poisson-link"),
    time = "year",
    spatial = "on",
    spatiotemporal = "iid",
    share_range = TRUE,
    anisotropy = TRUE,
    time_varying = ~1,
    time_varying_type = "rw0",
    priors = priors,
    control = sdmTMBcontrol(collapse_spatial_variance = TRUE),
    silent = TRUE
  )
}

summarise_fit <- function(fit) {
  if (inherits(fit, "error")) {
    return(list(sanity_ok = FALSE, index = NULL, sigma_V = NULL))
  }
  # tidy() omits sigma_V; one column per delta component
  list(
    sanity_ok = isTRUE(all(unlist(sanity(fit, gradient_thresh = 0.01, silent = TRUE)))),
    index = get_index_ok(fit, grid),
    sigma_V = tibble::tibble(
      component = c("encounter", "positive"),
      estimate = as.numeric(as.list(fit$sd_report, "Estimate", report = TRUE)$sigma_V),
      std.error = as.numeric(as.list(fit$sd_report, "Std. Error", report = TRUE)$sigma_V)
    )
  )
}

prior_name <- \(family, mean, spread) paste0(family, ": ", mean, ", ", spread)

fits_file <- here::here("data-generated/qcs-hs-biennial-lognormal-prior.rds")
if (!file.exists(fits_file)) {
  future::plan(future::multisession, workers = min(nrow(prior_grid), future::availableCores() / 2))
  res <- furrr::future_map2(prior_grid$mean, prior_grid$spread, \(m, s) {
    RhpcBLASctl::blas_set_num_threads(1L)
    RhpcBLASctl::omp_set_num_threads(1L)
    summarise_fit(tryCatch(fit_rw_year(lognormal_prior(m, s)), error = \(e) e))
  }, .options = furrr::furrr_options(seed = TRUE, scheduling = Inf)) |>
    setNames(prior_name("lognormal", prior_grid$mean, prior_grid$spread))
  future::plan(future::sequential)
  saveRDS(res, fits_file)
}
res_lognormal <- readRDS(fits_file)

# Gamma fits named to match
res_gamma <- readRDS(gamma_file)
truth <- res_gamma$truth$index |>
  select(year, true_est = est, true_lwr = lwr, true_upr = upr) |>
  left_join(occasions, by = "year")
res_gamma <- res_gamma[paste0("mean = ", prior_grid$mean, ", CV = ", prior_grid$spread)] |>
  setNames(prior_name("gamma", prior_grid$mean, prior_grid$spread))

res <- c(res_gamma, res_lognormal)
fit_info <- bind_rows(
  mutate(prior_grid, family = "gamma"),
  mutate(prior_grid, family = "lognormal")
) |>
  mutate(prior = prior_name(family, mean, spread)) |>
  mutate(sanity_ok = purrr::map_lgl(res[prior], "sanity_ok"))

indexes <- purrr::imap_dfr(res, \(r, name) {
  if (!is.null(r$index)) mutate(r$index, prior = name)
}) |>
  left_join(truth, by = "year") |>
  left_join(fit_info, by = "prior") |>
  mutate(log_error = log_est - log(true_est))

sigma_V <- purrr::imap_dfr(res, \(r, name) {
  if (!is.null(r$sigma_V)) mutate(r$sigma_V, prior = name)
}) |>
  left_join(fit_info, by = "prior")

# Metrics ----------------------------------------------------------------------

phase_of <- \(p) if_else(p == "SYN HS", 0.5, -0.5)
seesaw <- function(est, primary_survey, year) {
  mw <- moving_window(est, year = year, window = seesaw_window, phase = phase_of(primary_survey))
  tibble::tibble(
    A = period2_metric(est, year = year, phase = phase_of(primary_survey))[["A"]],
    A_max = max(mw$A)
  )
}

q <- qnorm(1 - (1 - ci_level) / 2)
metrics <- indexes |>
  arrange(prior, year) |>
  group_by(family, mean, spread, sanity_ok) |>
  summarise(
    rmse = sqrt(mean(log_error^2)),
    bias = mean(log_error),
    coverage = mean(abs(log_error) <= q * se),
    seesaw(est, primary_survey, year),
    .groups = "drop"
  ) |>
  left_join(
    sigma_V |>
      select(family, mean, spread, component, estimate) |>
      tidyr::pivot_wider(names_from = component, values_from = estimate, names_prefix = "sigma_V_"),
    by = c("family", "mean", "spread")
  ) |>
  arrange(mean, spread, family)

cat("\nTruth A:", round(seesaw(truth$true_est, truth$primary_survey, truth$year)$A, 1),
  " A_max:", round(seesaw(truth$true_est, truth$primary_survey, truth$year)$A_max, 1), "\n\n")
metrics |>
  mutate(across(where(is.numeric), \(x) round(x, 3))) |>
  print(n = Inf, width = Inf)

# Plots ----------------------------------------------------------------------

spread_lab <- \(family, spread) paste0(if_else(family == "gamma", "CV = ", "sdlog = "), spread)

indexes |>
  mutate(spread_lab = paste("Spread =", spread)) |>
  ggplot(aes(cal_year, est, colour = factor(mean))) +
  geom_ribbon(aes(cal_year, ymin = true_lwr, ymax = true_upr), data = truth, fill = "grey85", inherit.aes = FALSE) +
  geom_line(aes(y = true_est), data = truth, colour = "black") +
  geom_line() +
  geom_point(aes(shape = primary_survey)) +
  facet_grid(family ~ spread_lab) +
  scale_colour_viridis_d(end = 0.85) +
  scale_shape_manual(values = c(`SYN HS` = 19, `SYN QCS` = 21)) +
  labs(x = "Year", y = "Biomass index", colour = "Prior mean\n(gamma) or\nmedian\n(lognormal)", shape = "Region sampled")
ggsave(here::here("figs/qcs-hs-biennial-lognormal-prior-indexes.pdf"), width = 11, height = 5.5)

metric_levels <- c("rmse", "bias", "coverage", "A", "A_max", "sigma_V_encounter", "sigma_V_positive")
metrics |>
  tidyr::pivot_longer(all_of(metric_levels), names_to = "metric") |>
  mutate(metric = factor(metric, levels = metric_levels)) |>
  ggplot(aes(factor(spread), value, colour = factor(mean), linetype = family, group = interaction(mean, family))) +
  geom_line() +
  geom_point() +
  facet_wrap(~metric, scales = "free_y", nrow = 1) +
  scale_colour_viridis_d(end = 0.85) +
  labs(x = "Prior CV (gamma) or sdlog (lognormal)", y = NULL, colour = "Prior mean\nor median", linetype = "Prior")
ggsave(here::here("figs/qcs-hs-biennial-lognormal-prior-metrics.pdf"), width = 14, height = 3)

# Prior densities with the estimated sigma_V for each component
densities <- bind_rows(
  mutate(prior_grid, family = "gamma"),
  mutate(prior_grid, family = "lognormal")
) |>
  tidyr::crossing(x = seq(0.001, 1.5, length.out = 300)) |>
  mutate(density = if_else(
    family == "gamma",
    dgamma(x, shape = 1 / spread^2, rate = 1 / (spread^2 * mean)),
    dlnorm(x, meanlog = log(mean), sdlog = spread)
  ))

ggplot(densities, aes(x, density, colour = family)) +
  geom_line() +
  geom_vline(aes(xintercept = estimate, colour = family, linetype = component), data = sigma_V) +
  facet_grid(paste("Mean/median =", mean) ~ paste("Spread =", spread), scales = "free_y") +
  coord_cartesian(xlim = c(0, 1.2)) +
  scale_colour_brewer(palette = "Dark2") +
  labs(x = "sigma_V (per occasion)", y = "Prior density", colour = "Prior", linetype = "Estimate")
ggsave(here::here("figs/qcs-hs-biennial-lognormal-prior-sigma.pdf"), width = 9, height = 5.5)
