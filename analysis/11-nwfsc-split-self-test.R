# Simulation self-test of alternating north / south sampling of the NWFSC
# West Coast groundfish bottom trawl survey (WCGBTS) with a known truth.
#
# As analysis/08-qcs-hs-biennial-self-test.R, but the split is artificial: the
# WCGBTS samples the whole coast every year, so here the coast is cut at the
# area-weighted median latitude of the survey grid (fixed by geometry before
# looking at any results) and only one half is kept per survey occasion,
# alternating north / south. This gives a second, independent geometry (a
# long, narrow, along-coast domain like the BC surveys) with the same
# species as far as possible.
#
# 1. Fit the operating model (AR1 RF, factor(year)) to the full WCGBTS data.
# 2. For each replicate, draw the spatial fields from their approximate
#    posterior and entirely new spatiotemporal fields and observation error at
#    the observed tows and, on the same field draws, the expected density on
#    the prediction grid. The true index is the area-weighted sum of that
#    expected density.
# 3. Keep one half of the simulated tows per occasion and fit the alternating
#    models. Also refit the operating model to the full simulated data as a
#    calibration check.
# 4. Score RMSE, bias, CI coverage, and seesaw A against the known truth.

library(dplyr)
library(ggplot2)
library(sdmTMB)

theme_set(ggsidekick::theme_sleek())

source(here::here("analysis/fit-index-models.R"))
source(here::here("analysis/metric-functions.R"))

surveyjoin::cache_data()
surveyjoin::load_sql_data()

survey_name <- "NWFSC.Combo"
utm_crs <- 32610 # UTM 10N
# Present in both halves (pacific cod, pacific ocean perch, and pacific
# halibut are nearly absent south of the split); spans a range of N-S
# contrasts and overlaps the QCS / HS self-test species where possible
species_to_fit <- c(
  "spotted ratfish",
  "pacific spiny dogfish",
  "lingcod",
  "shortspine thornyhead",
  "arrowtooth flounder",
  "sablefish",
  "dover sole",
  "petrale sole"
)
n_reps <- 20L
ci_levels <- c(0.5, 0.95)
# Moving window (in survey occasions) for max A, as in the annual-survey
# scripts; 21 occasions gives 12 windows
seesaw_window <- 10L
# Coarsen the 2 km prediction grid to cut prediction cost
grid_res_km <- 6
# Predict the index in chunks of years to limit peak memory per worker
index_nsplit <- 4L

# Operating model, also refit to the full simulated data for calibration
om_model <- "AR1 RF, factor(year)"
# FALSE: spatial fields are a new draw from their approximate posterior each
# replicate (keeps the fitted gradient across the split, with uncertainty);
# TRUE: spatial fields are entirely new draws from the GMRF prior
resample_spatial <- FALSE

alternating_models <- c(
  "IID RF, factor(year)",
  "IID RF, RW year",
  "RW RF",
  "RW RF, RW year",
  "AR1 RF, RW year"
)

family <- delta_gamma(type = "poisson-link")
# Drop random fields with SDs collapsing to zero (e.g., the encounter
# spatiotemporal field with alternating data)
control <- sdmTMBcontrol(collapse_spatial_variance = TRUE)
rw_prior <- sdmTMBpriors(sigma_V = lognormal_prior(log(0.3), 0.5))

# Grid and split ------------------------------------------------------------

grid_fine <- surveyjoin::nwfsc_grid |>
  filter(survey == survey_name, area > 0)

# Area-weighted median latitude of the survey domain
split_lat <- with(grid_fine, {
  o <- order(lat)
  lat[o][which(cumsum(area[o]) >= sum(area) / 2)[1L]]
})
cli::cli_inform("Splitting north / south at {round(split_lat, 2)} degrees N")

# Assign the half at 2 km, then aggregate within each half so no coarse cell
# straddles the split
base_grid <- grid_fine |>
  mutate(region = if_else(lat >= split_lat, "North", "South")) |>
  sdmTMB::add_utm_columns(c("lon", "lat"), utm_crs = utm_crs) |>
  mutate(X = round(X / grid_res_km) * grid_res_km, Y = round(Y / grid_res_km) * grid_res_km) |>
  group_by(region, X, Y) |>
  summarise(area = sum(area), .groups = "drop")

# Data -----------------------------------------------------------------------

# Model time as survey occasion (no survey in 2020) and alternate the sampled
# half by occasion
prepare_species_data <- function(species) {
  surveyjoin::get_data(species, regions = "nwfsc") |>
    filter(survey_name == .env$survey_name, effort > 0) |>
    mutate(year = lubridate::year(lubridate::ymd(date))) |>
    select(year, lon_start, lat_start, depth_m, effort, catch_weight) |>
    tidyr::drop_na() |>
    mutate(
      region = if_else(lat_start >= split_lat, "North", "South"),
      occasion = match(year, sort(unique(year))),
      primary_region = if_else(occasion %% 2L == 1L, "North", "South"),
      cal_year = year,
      year = occasion
    ) |>
    sdmTMB::add_utm_columns(ll_names = c("lon_start", "lat_start"), utm_crs = utm_crs)
}

# Models ---------------------------------------------------------------------

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

# Operating model and simulation ---------------------------------------------

# Fit the operating model to one species and simulate its replicates. Returns
# only what the replicate fits and metrics need (not the fitted model).
setup_species <- function(species) {
  RhpcBLASctl::blas_set_num_threads(1L)
  RhpcBLASctl::omp_set_num_threads(1L)

  dat <- prepare_species_data(species)
  is_sampled <- dat$region == dat$primary_region
  grid <- sdmTMB::replicate_df(base_grid, "year", sort(unique(dat$year)))
  mesh <- make_mesh(dat, c("X", "Y"), cutoff = 15)
  mesh_alt <- make_mesh(dat[is_sampled, ], c("X", "Y"), mesh = mesh$mesh)

  fit_om <- fit_ok(fit_model(om_model, dat, mesh))
  if (is.null(fit_om)) {
    cli::cli_warn("Operating model did not pass sanity checks for {species}")
    return(NULL)
  }

  # Tows and grid cells in one newdata so both share the same field draws.
  # eta_i excludes the offset, so grid rows give expected density per unit
  # area; the simulated tow catches include the effort offset.
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

  true_indexes <- purrr::imap_dfr(sim_reports, \(r, rep) {
    grid |>
      mutate(density = exp(r$eta_i[!is_tow, 1L] + r$eta_i[!is_tow, 2L]) * area) |>
      group_by(year) |>
      summarise(true_est = sum(density)) |>
      mutate(rep = rep, .before = 1L)
  })

  list(
    species = species,
    dat = dat,
    is_sampled = is_sampled,
    grid = grid,
    mesh = mesh,
    mesh_alt = mesh_alt,
    occasions = distinct(dat, year, cal_year, primary_region),
    catches = lapply(sim_reports, \(r) r$y_i[is_tow, 1L] * r$y_i[is_tow, 2L]),
    true_indexes = true_indexes
  )
}

# Fit replicates --------------------------------------------------------------

fit_rep <- function(setup, rep) {
  RhpcBLASctl::blas_set_num_threads(1L)
  RhpcBLASctl::omp_set_num_threads(1L)

  d <- setup$dat
  d$catch_weight <- setup$catches[[rep]]

  index_ok <- function(model, d, mesh) {
    fit <- fit_ok(fit_model(model, d, mesh))
    if (is.null(fit)) tibble::tibble() else get_index_ok(fit, setup$grid, nsplit = index_nsplit)
  }

  bind_rows(
    full = index_ok(om_model, d, setup$mesh) |> mutate(model = om_model),
    alternating = purrr::map(
      alternating_models,
      \(m) index_ok(m, d[setup$is_sampled, ], setup$mesh_alt)
    ) |>
      setNames(alternating_models) |>
      bind_rows(.id = "model"),
    .id = "design"
  ) |>
    mutate(
      model = if_else(design == "full", paste0(model, " (full data)"), model),
      species = setup$species,
      rep = rep,
      .before = 1L
    )
}

if (FALSE) {
  tictoc::tic()
  s <- setup_species("lingcod")
  tictoc::toc()
  tictoc::tic()
  x <- fit_rep(s, 1L)
  tictoc::toc()
}

# One cache file per species; delete a file to refit that species
fits_dir <- here::here("data-generated/nwfsc-split-self-test")
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

  # One species x replicate per future so a slow fit doesn't hold up a chunk
  tasks <- tidyr::expand_grid(i = seq_along(setups), rep = seq_len(n_reps))
  rep_indexes <- furrr::future_map2(
    tasks$i, tasks$rep, \(i, rep) fit_rep(setups[[i]], rep),
    .options = furrr::furrr_options(seed = TRUE, scheduling = Inf)
  )
  future::plan(future::sequential)

  for (i in seq_along(setups)) {
    saveRDS(
      list(
        indexes = bind_rows(rep_indexes[tasks$i == i]),
        true_indexes = mutate(setups[[i]]$true_indexes, species = setups[[i]]$species, .before = 1L),
        occasions = mutate(setups[[i]]$occasions, species = setups[[i]]$species, .before = 1L)
      ),
      fits_file(setups[[i]]$species)
    )
  }
}

self_test <- lapply(fits_file(species_to_fit), \(f) if (file.exists(f)) readRDS(f))
indexes <- purrr::map_dfr(self_test, "indexes")
true_indexes <- purrr::map_dfr(self_test, "true_indexes")
occasions <- purrr::map_dfr(self_test, "occasions")

# Performance metrics ----------------------------------------------------------

compared <- indexes |>
  left_join(true_indexes, by = c("species", "rep", "year")) |>
  left_join(occasions, by = c("species", "year")) |>
  mutate(log_error = log_est - log(true_est))

for (level in ci_levels) {
  q <- qnorm(1 - (1 - level) / 2)
  compared[[paste0("covered_", level * 100)]] <- abs(compared$log_error) <= q * compared$se
}

# Phase follows which half was sampled
seesaw_A <- function(est, primary_region, year) {
  as_tibble(t(period2_metric(
    est,
    year = year,
    phase = if_else(primary_region == "North", 0.5, -0.5)
  )))
}

# Maximum A across moving windows
seesaw_A_max <- function(est, primary_region, year) {
  mw <- moving_window(
    est,
    year = year,
    window = seesaw_window,
    phase = if_else(primary_region == "North", 0.5, -0.5)
  )
  max(mw$A)
}

metrics_rep <- compared |>
  arrange(species, rep, model, year) |>
  group_by(species, rep, model) |>
  summarise(
    rmse = sqrt(mean(log_error^2)),
    bias = mean(log_error),
    coverage_50 = mean(covered_50),
    coverage_95 = mean(covered_95),
    A = seesaw_A(est, primary_region, year)$A,
    A_max = seesaw_A_max(est, primary_region, year),
    .groups = "drop"
  )

A_truth <- true_indexes |>
  left_join(occasions, by = c("species", "year")) |>
  arrange(species, rep, year) |>
  group_by(species, rep) |>
  summarise(
    A = seesaw_A(true_est, primary_region, year)$A,
    A_max = seesaw_A_max(true_est, primary_region, year),
    .groups = "drop"
  )

metrics <- metrics_rep |>
  group_by(species, model) |>
  summarise(
    n_reps = n(),
    rmse = mean(rmse),
    bias = mean(bias),
    coverage_50 = mean(coverage_50),
    coverage_95 = mean(coverage_95),
    A_median = median(A),
    A_max_median = median(A_max),
    .groups = "drop"
  ) |>
  arrange(species, rmse)

saveRDS(
  list(compared = compared, metrics_rep = metrics_rep, metrics = metrics, A_truth = A_truth),
  here::here("data-generated/nwfsc-split-self-test-metrics.rds")
)

metrics |>
  mutate(across(where(is.numeric), \(x) round(x, 2))) |>
  print(n = Inf, width = Inf)

# Plots ----------------------------------------------------------------------

metric_levels <- c("rmse", "bias", "coverage_50", "coverage_95", "A", "A_max")
metric_labels <- c(
  "RMSE (log)", "Bias (log)", "50% CI coverage", "95% CI coverage",
  "Seesaw A (%)", paste0("Max A, ", seesaw_window, "-survey windows (%)")
)
# Nominal values; for A and max A, the median of the true index by species
nominal <- bind_rows(
  tidyr::expand_grid(
    species = unique(metrics_rep$species),
    tibble::tibble(metric = c("bias", "coverage_50", "coverage_95"), value = c(0, 0.5, 0.95))
  ),
  A_truth |>
    group_by(species) |>
    summarise(A = median(A), A_max = median(A_max)) |>
    tidyr::pivot_longer(c(A, A_max), names_to = "metric")
) |>
  mutate(metric = factor(metric, levels = metric_levels, labels = metric_labels))

# Models ordered by mean RMSE across species
model_order <- metrics |>
  group_by(model) |>
  summarise(rmse = mean(rmse)) |>
  arrange(rmse) |>
  pull(model)

metrics_rep |>
  tidyr::pivot_longer(all_of(metric_levels), names_to = "metric") |>
  mutate(
    model = factor(model, levels = rev(model_order)),
    metric = factor(metric, levels = metric_levels, labels = metric_labels)
  ) |>
  ggplot(aes(value, model)) +
  geom_vline(aes(xintercept = value), data = nominal, linetype = 2, colour = "grey50") +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter(height = 0.15), alpha = 0.4) +
  facet_grid(species ~ metric, scales = "free_x") +
  labs(
    x = NULL, y = NULL,
    caption = paste0(
      "WCGBTS alternating N / S self-test (split at ", round(split_lat, 2), " N), ",
      n_reps, " replicates per species, OM: ", om_model,
      if (resample_spatial) " (spatial fields new from prior). " else " (spatial fields posterior draws). ",
      "Dashed lines: nominal values (A, max A: median of the true index)."
    )
  )
ggsave(
  here::here("figs/nwfsc-split-self-test.pdf"),
  width = 16, height = 1 + 2.2 * length(unique(metrics_rep$species))
)

# Summary across species: one point per species, each a summary across
# replicates, scaled so species are comparable; violins show all replicates
# pooled across species. RMSE is relative to the species' mean full-data
# operating model RMSE (how much is lost by sampling one half per occasion).
# Max A is the excess over the true index in the same replicate (spurious
# seesaw from the model, not from the simulated population).
summary_levels <- c("rmse_ratio", "coverage_50", "A_max_excess")
summary_labels <- c("RMSE / full-data RMSE", "50% CI coverage", "Max A - true max A (%)")
full_model <- paste0(om_model, " (full data)")

rep_summary <- metrics_rep |>
  left_join(select(A_truth, species, rep, A_max_true = A_max), by = c("species", "rep")) |>
  group_by(species) |>
  mutate(
    rmse_ratio = rmse / mean(rmse[model == full_model]),
    A_max_excess = A_max - A_max_true
  ) |>
  ungroup() |>
  tidyr::pivot_longer(all_of(summary_levels), names_to = "metric") |>
  mutate(
    model = factor(model, levels = rev(model_order)),
    metric = factor(metric, levels = summary_levels, labels = summary_labels)
  )

species_summary <- rep_summary |>
  group_by(species, model, metric) |>
  summarise(
    value = if (grepl("A", first(metric))) median(value) else mean(value),
    .groups = "drop"
  )

summary_nominal <- tibble::tibble(
  metric = factor(summary_labels, levels = summary_labels),
  value = c(1, 0.5, 0)
)

# Zoom each panel to the widest 10-90% quantile range of the replicates across
# models, extended to keep the species points and nominal values in view.
# oob_keep so the violins are still computed from all replicates (as with
# coord_cartesian(), which can't be set per panel).
summary_limits <- bind_rows(
  rep_summary |>
    group_by(metric, model) |>
    summarise(lo = quantile(value, 0.1), hi = quantile(value, 0.9), .groups = "drop"),
  rename(species_summary, lo = value) |> mutate(hi = lo),
  rename(summary_nominal, lo = value) |> mutate(hi = lo)
) |>
  group_by(metric) |>
  summarise(lo = min(lo), hi = max(hi)) |>
  arrange(metric)
summary_x_scales <- purrr::map2(
  summary_limits$lo, summary_limits$hi,
  \(lo, hi) scale_x_continuous(limits = c(lo, hi), oob = scales::oob_keep)
)

ggplot(species_summary, aes(value, model)) +
  geom_violin(data = rep_summary, fill = "grey90", colour = NA, scale = "width") +
  geom_vline(aes(xintercept = value), data = summary_nominal, linetype = 2, colour = "grey50") +
  stat_summary(fun = median, geom = "crossbar", fun.min = median, fun.max = median, width = 0.6) +
  geom_point(aes(colour = species), position = position_jitter(height = 0.15, seed = 1)) +
  facet_wrap(~metric, nrow = 1, scales = "free_x") +
  ggh4x::facetted_pos_scales(x = summary_x_scales) +
  scale_colour_brewer(palette = "Dark2", labels = tools::toTitleCase) +
  labs(x = NULL, y = NULL, colour = NULL) +
  theme(legend.position = "bottom", plot.margin = margin(5.5, 12, 5.5, 5.5)) +
  guides(colour = guide_legend(nrow = 2))
ggsave(here::here("figs/nwfsc-split-self-test-summary.pdf"), width = 9, height = 4)

compared |>
  filter(rep <= 3) |>
  ggplot(aes(cal_year, est, colour = model)) +
  geom_line(aes(y = true_est), colour = "black") +
  geom_line() +
  facet_grid(species ~ rep, scales = "free_y", labeller = labeller(rep = label_both)) +
  labs(x = "Year", y = "Biomass index", colour = "Model", caption = "Black: true index")
ggsave(
  here::here("figs/nwfsc-split-self-test-indexes.pdf"),
  width = 11, height = 1 + 2 * length(unique(compared$species))
)
