# Downgrade QCS + HS to a pretend biennial design and compare index models
# against the full-domain index.
#
# 1. Take years where both SYN QCS and SYN HS were sampled.
# 2. 'Truth': IID RF, factor(year) fit to the full data (both regions every
#    year).
# 3. Biennial: keep only one region per occasion, alternating HS / QCS.
# 4. Fit the shared set of models (analysis/fit-index-models.R) to the
#    biennial data.
# 5. For each model on the biennial data calculate:
#    - seesaw amplitude A (analysis/metric-functions.R)
#    - RMSE and bias of log(index) vs. truth
#    - 50% CI coverage of the truth

library(dplyr)
library(ggplot2)
library(sdmTMB)

theme_set(ggsidekick::theme_sleek())

source(here::here("analysis/fit-index-models.R"))
source(here::here("analysis/metric-functions.R"))

surveyjoin::cache_data()
surveyjoin::load_sql_data()

survey_names <- c("SYN QCS", "SYN HS")
ci_level <- 0.5
# Moving window (in survey occasions) for max A; as in
# analysis/05-qcs-hs-experimental-biennial.R. With 11 occasions, 8 gives 4
# windows (10, as in the annual-survey scripts, would give only 2).
seesaw_window <- 8L

# Strong seesaws under biennial QCS / HS sampling across a range of taxa; see
# analysis/qcs-hs-biennial-screen.R
species_to_fit <- c(
  "pacific cod",
  "spotted ratfish",
  "pacific ocean perch",
  "pacific halibut",
  "pacific spiny dogfish",
  "lingcod",
  "shortspine thornyhead"
)

# Data -----------------------------------------------------------------------

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
      primary_survey = if_else(occasion %% 2L == 1L, "SYN HS", "SYN QCS")
    ) |>
    sdmTMB::add_utm_columns(ll_names = c("lon_start", "lat_start"), utm_crs = 3156)
}

# Fitting --------------------------------------------------------------------

fit_species <- function(species) {
  RhpcBLASctl::blas_set_num_threads(1L)
  RhpcBLASctl::omp_set_num_threads(1L)

  # Model time as survey occasion: QCS + HS years are all odd, so on the
  # calendar-year scale extra_time puts an empty year between every survey and
  # AR1 rho gets stuck at its starting value of 0 (only rho^2 is identified).
  dat <- prepare_species_data(species) |>
    mutate(cal_year = year, year = occasion)
  dat_biennial <- filter(dat, survey_name == primary_survey)

  grid <- surveyjoin::dfo_synoptic_grid |>
    filter(survey %in% survey_names) |>
    sdmTMB::add_utm_columns(c("lon", "lat"), utm_crs = 3156) |>
    clamp_depth(dat) |>
    sdmTMB::replicate_df("year", sort(unique(dat$year)))

  # Same knots for both designs
  mesh <- sdmTMB::make_mesh(dat, c("X", "Y"), cutoff = 10)
  mesh_biennial <- sdmTMB::make_mesh(dat_biennial, c("X", "Y"), mesh = mesh$mesh)

  family <- sdmTMB::delta_gamma(type = "poisson-link")
  # Drop random fields with SDs collapsing to zero (e.g., the encounter
  # spatiotemporal field with biennial data)
  control <- sdmTMB::sdmTMBcontrol(collapse_spatial_variance = TRUE)

  # Truth: same specification as the base model in fit_index_models()
  fit_truth <- fit_ok(sdmTMB::sdmTMB(
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
    silent = FALSE
  ))
  if (is.null(fit_truth)) {
    cli::cli_warn("Full-data reference model did not converge for {species}")
    return(tibble::tibble())
  }

  occasions <- distinct(dat, occasion, cal_year, primary_survey)

  bind_rows(
    full = get_index_ok(fit_truth, grid) |>
      mutate(model = "IID RF, factor(year)", .before = 1L),
    biennial = fit_index_models(
      dat = dat_biennial,
      grid = grid,
      mesh = mesh_biennial,
      response = "catch_weight",
      family = family,
      offset = log(dat_biennial$effort),
      control = control
    ),
    .id = "design"
  ) |>
    rename(occasion = year) |>
    left_join(occasions, by = "occasion") |>
    rename(year = cal_year) |>
    mutate(species = species, .before = 1L)
}

# One cache file per species; delete a file to refit that species
fits_dir <- here::here("data-generated/qcs-hs-biennial-models")
dir.create(fits_dir, showWarnings = FALSE)
fits_file <- \(species) file.path(fits_dir, paste0(gsub(" ", "-", species), ".rds"))

to_fit <- species_to_fit[!file.exists(fits_file(species_to_fit))]
if (length(to_fit) > 0L) {
  fit_and_save <- function(species) {
    out <- fit_species(species)
    saveRDS(out, fits_file(species))
    out
  }
  # Run sequentially for one species so the fitting progress prints
  if (length(to_fit) == 1L) {
    fit_and_save(to_fit)
  } else {
    future::plan(future::multisession, workers = min(length(to_fit), future::availableCores() / 2))
    furrr::future_walk(to_fit, fit_and_save, .options = furrr::furrr_options(seed = TRUE, scheduling = Inf))
    future::plan(future::sequential)
  }
}
indexes <- purrr::map_dfr(fits_file(species_to_fit), readRDS)

# Performance metrics ----------------------------------------------------------

truth <- indexes |>
  filter(design == "full") |>
  select(species, year, true_est = est, true_lwr = lwr, true_upr = upr)

biennial <- filter(indexes, design == "biennial")

compared <- biennial |>
  left_join(truth, by = c("species", "year")) |>
  mutate(
    q = qnorm(1 - (1 - ci_level) / 2),
    covered = log(true_est) >= log_est - q * se & log(true_est) <= log_est + q * se,
    log_error = log_est - log(true_est)
  )

# Phase follows which region was sampled, not calendar-year parity. A_max is
# the maximum A across moving windows.
seesaw_A <- function(d) {
  d <- arrange(d, year)
  phase <- if_else(d$primary_survey == "SYN HS", 0.5, -0.5)
  mw <- moving_window(d$est, year = d$year, window = seesaw_window, phase = phase)
  as_tibble(t(period2_metric(d$est, year = d$year, phase = phase))) |>
    mutate(A_max = max(mw$A))
}

metrics <- compared |>
  group_by(species, model) |>
  summarise(
    rmse = sqrt(mean(log_error^2)),
    bias = mean(log_error),
    coverage = mean(covered),
    n_years = n(),
    .groups = "drop"
  )

A_biennial <- biennial |>
  group_by(species, model) |>
  group_modify(\(.x, .y) seesaw_A(.x)) |>
  ungroup()

A_truth <- indexes |>
  filter(design == "full") |>
  group_by(species) |>
  group_modify(\(.x, .y) seesaw_A(.x)) |>
  ungroup()

metrics <- metrics |>
  left_join(select(A_biennial, species, model, A, A_lower, A_upper, A_max), by = c("species", "model"))

saveRDS(
  list(indexes = indexes, metrics = metrics, A_biennial = A_biennial, A_truth = A_truth),
  here::here("data-generated/qcs-hs-biennial-metrics.rds")
)

metrics |>
  arrange(species, rmse) |>
  print(n = Inf)

A_truth

# Plots ----------------------------------------------------------------------

model_order <- metrics |>
  group_by(model) |>
  summarise(rmse = mean(rmse)) |>
  arrange(rmse) |>
  pull(model)

biennial |>
  mutate(model = factor(model, levels = model_order)) |>
  ggplot(aes(year, est)) +
  geom_ribbon(aes(y = true_est, ymin = true_lwr, ymax = true_upr), data = truth, fill = "grey80") +
  geom_line(aes(y = true_est), data = truth, colour = "grey30") +
  geom_pointrange(aes(ymin = exp(log_est - qnorm(0.75) * se), ymax = exp(log_est + qnorm(0.75) * se), colour = primary_survey)) +
  geom_line(colour = "grey50", linetype = 2) +
  facet_grid(species ~ model, scales = "free_y") +
  labs(
    x = "Year", y = "Biomass index",
    colour = "Sampled region",
    caption = "Grey: full-data reference index (95% CI). Points: biennial index (50% CI)."
  )
ggsave(here::here("figs/qcs-hs-biennial-indexes.pdf"), width = 16, height = 3 + 2.5 * length(unique(biennial$species)))

metrics |>
  tidyr::pivot_longer(c(rmse, bias, coverage, A, A_max), names_to = "metric") |>
  mutate(
    model = factor(model, levels = rev(model_order)),
    metric = factor(metric, levels = c("rmse", "bias", "coverage", "A", "A_max"),
      labels = c(
        "RMSE (log)", "Bias (log)", "50% CI coverage", "Seesaw A (%)",
        paste0("Max A, ", seesaw_window, "-survey windows (%)")
      ))
  ) |>
  ggplot(aes(value, model, colour = species)) +
  geom_point() +
  facet_wrap(~metric, scales = "free_x", nrow = 1) +
  labs(x = NULL, y = NULL, colour = "Species")
ggsave(here::here("figs/qcs-hs-biennial-metrics.pdf"), width = 14, height = 4)
