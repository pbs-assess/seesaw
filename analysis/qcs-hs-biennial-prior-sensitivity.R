# Sensitivity of the IID RF, RW year model to the gamma prior on sigma_V
# (the SD of the random walk intercept) for Pacific cod under pretend biennial
# QCS / HS sampling.
#
# As in analysis/07-qcs-hs-biennial-model-comparison.R, the 'truth' is IID RF,
# factor(year) fit to the full data (both regions every year). Each prior in a
# grid of means x CVs is fit to the biennial subset, plus a fit with no prior.
# Time is survey occasion, so sigma_V is the SD of the change in the intercept
# per occasion (2 calendar years).

library(dplyr)
library(ggplot2)
library(sdmTMB)

theme_set(ggsidekick::theme_sleek())

source(here::here("analysis/fit-index-models.R"))
source(here::here("analysis/metric-functions.R"))

surveyjoin::cache_data()
surveyjoin::load_sql_data()

species <- "pacific cod"
survey_names <- c("SYN QCS", "SYN HS")
ci_level <- 0.5
seesaw_window <- 8L

prior_grid <- tidyr::crossing(mean = c(0.1, 0.3, 0.5), cv = c(0.1, 0.3, 0.5, 0.8)) |>
  mutate(prior = paste0("mean = ", mean, ", CV = ", cv))

# Data -----------------------------------------------------------------------

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

# Model time as survey occasion (see analysis/07-qcs-hs-biennial-model-comparison.R)
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

family <- delta_gamma(type = "poisson-link")
control <- sdmTMBcontrol(collapse_spatial_variance = TRUE)

# Fits -----------------------------------------------------------------------

fit_rw_year <- function(priors) {
  sdmTMB(
    catch_weight ~ 1,
    data = dat_biennial,
    mesh = mesh_biennial,
    offset = log(dat_biennial$effort),
    family = family,
    time = "year",
    spatial = "on",
    spatiotemporal = "iid",
    share_range = TRUE,
    anisotropy = TRUE,
    time_varying = ~1,
    time_varying_type = "rw0",
    priors = priors,
    control = control,
    silent = TRUE
  )
}

# Keep fits that fail sanity checks so they can be flagged rather than dropped
summarise_fit <- function(fit) {
  if (inherits(fit, "error")) {
    return(list(sanity_ok = FALSE, index = NULL, sigma_V = NULL))
  }
  # tidy() omits sigma_V; one column per delta component
  est <- as.list(fit$sd_report, "Estimate", report = TRUE)$sigma_V
  sigma_V <- if (length(est) > 0L) {
    tibble::tibble(
      component = c("encounter", "positive"),
      estimate = as.numeric(est),
      std.error = as.numeric(as.list(fit$sd_report, "Std. Error", report = TRUE)$sigma_V)
    )
  }
  list(
    sanity_ok = isTRUE(all(unlist(sanity(fit, gradient_thresh = 0.01, silent = TRUE)))),
    index = get_index_ok(fit, grid),
    sigma_V = sigma_V
  )
}

fits_file <- here::here("data-generated/qcs-hs-biennial-prior-sensitivity.rds")
if (!file.exists(fits_file)) {
  jobs <- c(
    list(truth = NULL, "no prior" = sdmTMBpriors()),
    setNames(
      purrr::map2(prior_grid$mean, prior_grid$cv, \(m, cv) sdmTMBpriors(sigma_V = gamma_cv(m, cv))),
      prior_grid$prior
    )
  )

  future::plan(future::multisession, workers = min(length(jobs), future::availableCores() / 2))
  res <- furrr::future_imap(jobs, \(priors, name) {
    RhpcBLASctl::blas_set_num_threads(1L)
    RhpcBLASctl::omp_set_num_threads(1L)
    fit <- tryCatch(
      if (name == "truth") {
        sdmTMB(
          catch_weight ~ 0 + factor(year),
          data = dat,
          mesh = mesh,
          offset = log(dat$effort),
          family = family,
          time = "year",
          spatial = "on",
          spatiotemporal = "iid",
          share_range = TRUE,
          anisotropy = TRUE,
          control = control,
          silent = TRUE
        )
      } else {
        fit_rw_year(priors)
      },
      error = \(e) e
    )
    summarise_fit(fit)
  }, .options = furrr::furrr_options(seed = TRUE, scheduling = Inf))
  future::plan(future::sequential)
  saveRDS(res, fits_file)
}
res <- readRDS(fits_file)

stopifnot(res$truth$sanity_ok)
truth <- res$truth$index |>
  select(year, true_est = est, true_lwr = lwr, true_upr = upr) |>
  left_join(occasions, by = "year")

fit_info <- tibble::tibble(prior = names(res)[-1], sanity_ok = purrr::map_lgl(res[-1], "sanity_ok")) |>
  left_join(prior_grid, by = "prior")

indexes <- purrr::imap_dfr(res[-1], \(r, name) {
  if (!is.null(r$index)) mutate(r$index, prior = name)
}) |>
  left_join(truth, by = "year") |>
  left_join(fit_info, by = "prior") |>
  mutate(log_error = log_est - log(true_est))

sigma_V <- purrr::imap_dfr(res[-1], \(r, name) {
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
  group_by(prior, mean, cv, sanity_ok) |>
  summarise(
    rmse = sqrt(mean(log_error^2)),
    bias = mean(log_error),
    coverage = mean(abs(log_error) <= q * se),
    seesaw(est, primary_survey, year),
    .groups = "drop"
  ) |>
  left_join(
    sigma_V |>
      select(prior, component, estimate) |>
      tidyr::pivot_wider(names_from = component, values_from = estimate, names_prefix = "sigma_V_"),
    by = "prior"
  ) |>
  arrange(rmse)

# Empirical anchor: SD of the change in the log truth between occasions
truth_sd_change <- sd(diff(log(truth$true_est)))

cat("\nTruth A:", round(seesaw(truth$true_est, truth$primary_survey, truth$year)$A, 1),
  " A_max:", round(seesaw(truth$true_est, truth$primary_survey, truth$year)$A_max, 1),
  " SD of log change per occasion:", round(truth_sd_change, 3), "\n\n")
metrics |>
  mutate(across(where(is.numeric), \(x) round(x, 3))) |>
  print(n = Inf, width = Inf)

# Plots ----------------------------------------------------------------------

prior_labels <- \(d) mutate(d, mean_lab = paste("Prior mean =", mean), cv_lab = paste("CV =", cv))

indexes |>
  filter(prior != "no prior") |>
  prior_labels() |>
  ggplot(aes(cal_year, colour = factor(mean))) +
  geom_ribbon(aes(cal_year, ymin = true_lwr, ymax = true_upr), data = truth, fill = "grey85", inherit.aes = FALSE) +
  geom_line(aes(y = true_est), data = truth, colour = "black") +
  geom_line(aes(y = est)) +
  geom_point(aes(y = est, shape = primary_survey)) +
  facet_wrap(~cv_lab, nrow = 1) +
  scale_colour_viridis_d(end = 0.85) +
  scale_shape_manual(values = c(`SYN HS` = 19, `SYN QCS` = 21)) +
  labs(
    x = "Year", y = "Biomass index", colour = "Prior mean", shape = "Region sampled"
  )
ggsave(here::here("figs/qcs-hs-biennial-prior-sensitivity-indexes.pdf"), width = 11, height = 3.5)

# Prior densities with the estimated sigma_V for each component
x <- seq(0.001, 1.5, length.out = 300)
densities <- prior_grid |>
  tidyr::crossing(x = x) |>
  mutate(
    shape = 1 / cv^2,
    density = dgamma(x, shape = shape, rate = shape / mean)
  ) |>
  prior_labels()

ggplot(densities, aes(x, density)) +
  geom_line(colour = "grey40") +
  geom_vline(
    aes(xintercept = estimate, colour = component),
    data = prior_labels(filter(sigma_V, prior != "no prior"))
  ) +
  facet_grid(mean_lab ~ cv_lab, scales = "free_y") +
  coord_cartesian(xlim = c(0, 1.5)) +
  labs(x = "sigma_V (per occasion)", y = "Prior density", colour = "Estimate")
ggsave(here::here("figs/qcs-hs-biennial-prior-sensitivity-sigma.pdf"), width = 7.5, height = 5)

metrics |>
  filter(prior != "no prior") |>
  tidyr::pivot_longer(c(rmse, bias, coverage, A, A_max), names_to = "metric") |>
  mutate(metric = factor(metric, levels = c("rmse", "bias", "coverage", "A", "A_max"))) |>
  ggplot(aes(factor(cv), value, colour = factor(mean), group = factor(mean))) +
  geom_hline(
    aes(yintercept = value),
    data = metrics |>
      filter(prior == "no prior") |>
      tidyr::pivot_longer(c(rmse, bias, coverage, A, A_max), names_to = "metric") |>
      mutate(metric = factor(metric, levels = c("rmse", "bias", "coverage", "A", "A_max"))),
    linetype = 2, colour = "grey50"
  ) +
  geom_line() +
  geom_point() +
  facet_wrap(~metric, scales = "free_y", nrow = 1) +
  labs(x = "Prior CV", y = NULL, colour = "Prior mean", caption = "Dashed: no prior.")
ggsave(here::here("figs/qcs-hs-biennial-prior-sensitivity-metrics.pdf"), width = 12, height = 3)
