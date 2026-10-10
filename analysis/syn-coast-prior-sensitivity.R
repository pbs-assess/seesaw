# Sensitivity of the IID RF, RW year model to the gamma prior on sigma_V
# (the SD of the random walk intercept) for a coastwide synoptic index of
# walleye pollock. The synoptic surveys alternate QCS + HS (odd years) and
# WCHG + WCVI (even years), so there is no full-data reference here; IID RF,
# factor(year) is shown for comparison only.
#
# Data as in analysis/03-fit-syn-models.R. Time is calendar year, so sigma_V is
# the SD of the change in the intercept per year (2020 is missing and filled
# with extra_time).

library(dplyr)
library(ggplot2)
library(sdmTMB)

theme_set(ggsidekick::theme_sleek())

source(here::here("analysis/fit-index-models.R"))
source(here::here("analysis/metric-functions.R"))

surveyjoin::cache_data()
surveyjoin::load_sql_data()

species <- "walleye pollock"
seesaw_window <- 10L

prior_grid <- tidyr::crossing(mean = c(0.1, 0.3, 0.5), cv = c(0.1, 0.3, 0.5, 0.8)) |>
  mutate(prior = paste0("mean = ", mean, ", CV = ", cv))

# Data -----------------------------------------------------------------------

dat <- surveyjoin::get_data(species, regions = "pbs") |>
  mutate(year = lubridate::year(lubridate::ymd(date))) |>
  select(survey_name, year, lon_start, lat_start, depth_m, effort, catch_weight) |>
  sdmTMB::add_utm_columns(ll_names = c("lon_start", "lat_start"), utm_crs = 3156) |>
  # Use only complete N/S sampling years
  filter(!(year %in% c(2003, 2004, 2020))) |>
  tidyr::drop_na(effort, catch_weight, depth_m) |>
  # Drop these surveys to be perfectly biennial
  filter(!(year == 2007 & survey_name == "SYN WCHG")) |>
  filter(!(year == 2021 & survey_name == "SYN WCVI"))

survey_group <- \(year) if_else(year %% 2 == 0, "WCHG + WCVI", "QCS + HS")

grid <- surveyjoin::dfo_synoptic_grid |>
  sdmTMB::add_utm_columns(c("lon", "lat"), utm_crs = 3156) |>
  clamp_depth(dat) |>
  sdmTMB::replicate_df("year", sort(unique(dat$year)))

mesh <- make_mesh(dat, c("X", "Y"), cutoff = 10)

family <- delta_gamma(type = "poisson-link")
control <- sdmTMBcontrol(collapse_spatial_variance = TRUE)

# Fits -----------------------------------------------------------------------

fit_model <- function(name, priors) {
  args <- list(
    formula = catch_weight ~ 1,
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
  extra <- if (name == "factor(year)") {
    list(formula = catch_weight ~ 0 + factor(year))
  } else {
    list(
      time_varying = ~1,
      time_varying_type = "rw0",
      priors = priors,
      extra_time = seq(min(dat$year), max(dat$year))
    )
  }
  args[names(extra)] <- extra
  do.call(sdmTMB, args)
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

fits_file <- here::here("data-generated/syn-coast-prior-sensitivity.rds")
if (!file.exists(fits_file)) {
  jobs <- c(
    list("factor(year)" = NULL, "no prior" = sdmTMBpriors()),
    setNames(
      purrr::map2(prior_grid$mean, prior_grid$cv, \(m, cv) sdmTMBpriors(sigma_V = gamma_cv(m, cv))),
      prior_grid$prior
    )
  )

  future::plan(future::multisession, workers = min(length(jobs), future::availableCores() / 2))
  res <- furrr::future_imap(jobs, \(priors, name) {
    RhpcBLASctl::blas_set_num_threads(1L)
    RhpcBLASctl::omp_set_num_threads(1L)
    summarise_fit(tryCatch(fit_model(name, priors), error = \(e) e))
  }, .options = furrr::furrr_options(seed = TRUE, scheduling = Inf))
  future::plan(future::sequential)
  saveRDS(res, fits_file)
}
res <- readRDS(fits_file)

fit_info <- tibble::tibble(prior = names(res), sanity_ok = purrr::map_lgl(res, "sanity_ok")) |>
  left_join(prior_grid, by = "prior")

indexes <- purrr::imap_dfr(res, \(r, name) {
  if (!is.null(r$index)) mutate(r$index, prior = name)
}) |>
  left_join(fit_info, by = "prior") |>
  mutate(survey_group = survey_group(year))

sigma_V <- purrr::imap_dfr(res, \(r, name) {
  if (!is.null(r$sigma_V)) mutate(r$sigma_V, prior = name)
}) |>
  left_join(fit_info, by = "prior")

# Metrics ----------------------------------------------------------------------

seesaw <- function(est, year) {
  mw <- moving_window(est, year = year, window = seesaw_window)
  tibble::tibble(A = period2_metric(est, year = year)[["A"]], A_max = max(mw$A))
}

metrics <- indexes |>
  arrange(prior, year) |>
  group_by(prior, mean, cv, sanity_ok) |>
  summarise(mean_cv = mean(sqrt(exp(se^2) - 1)), seesaw(est, year), .groups = "drop") |>
  left_join(
    sigma_V |>
      select(prior, component, estimate) |>
      tidyr::pivot_wider(names_from = component, values_from = estimate, names_prefix = "sigma_V_"),
    by = "prior"
  ) |>
  arrange(mean, cv)

metrics |>
  mutate(across(where(is.numeric), \(x) round(x, 3))) |>
  print(n = Inf, width = Inf)

# Plots ----------------------------------------------------------------------

prior_labels <- \(d) mutate(d, mean_lab = paste("Prior mean =", mean), cv_lab = paste("CV =", cv))

factor_year <- filter(indexes, prior == "factor(year)") |>
  select(year, est, lwr, upr, survey_group)

rw_indexes <- filter(indexes, !prior %in% c("factor(year)", "no prior"))

rw_indexes |>
  prior_labels() |>
  ggplot(aes(year, est, colour = factor(mean))) +
  geom_ribbon(aes(year, ymin = lwr, ymax = upr), data = factor_year, fill = "grey90", inherit.aes = FALSE) +
  geom_line(data = factor_year, colour = "grey50") +
  geom_line() +
  geom_point(aes(shape = survey_group)) +
  facet_wrap(~cv_lab, nrow = 1) +
  # Clip the factor(year) ribbon and line to the range of the RW year indexes
  coord_cartesian(ylim = c(0, max(rw_indexes$est))) +
  scale_colour_viridis_d(end = 0.85) +
  scale_shape_manual(values = c(`QCS + HS` = 19, `WCHG + WCVI` = 21)) +
  labs(x = "Year", y = "Biomass index", colour = "Prior mean", shape = "Region sampled")
ggsave(here::here("figs/syn-coast-prior-sensitivity-indexes.pdf"), width = 11, height = 3.5)

# Prior densities with the estimated sigma_V for each component
densities <- prior_grid |>
  tidyr::crossing(x = seq(0.001, 1.5, length.out = 300)) |>
  mutate(
    shape = 1 / cv^2,
    density = dgamma(x, shape = shape, rate = shape / mean)
  ) |>
  prior_labels()

ggplot(densities, aes(x, density)) +
  geom_line(colour = "grey40") +
  geom_vline(
    aes(xintercept = estimate, colour = component),
    data = prior_labels(filter(sigma_V, !prior %in% c("factor(year)", "no prior")))
  ) +
  facet_grid(mean_lab ~ cv_lab, scales = "free_y") +
  coord_cartesian(xlim = c(0, 1.5)) +
  labs(x = "sigma_V (per year)", y = "Prior density", colour = "Estimate")
ggsave(here::here("figs/syn-coast-prior-sensitivity-sigma.pdf"), width = 7.5, height = 5)
