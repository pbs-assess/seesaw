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
# 1. Fit the operating model (IID RF, factor(year)) to the full WCGBTS data.
# 2. For each replicate, draw the spatial fields from their approximate
#    posterior and entirely new spatiotemporal fields and observation error at
#    the observed tows and, on the same field draws, the expected density on
#    the prediction grid. The true index is the area-weighted sum of that
#    expected density.
# 3. Keep one half of the simulated tows per occasion and fit the alternating
#    models. Also refit the operating model to the full simulated data as a
#    calibration check.
# 4. Score RMSE, bias, CI coverage, seesaw A, and the signed seesaw error
#    against the known truth.
#
# The pipeline and metrics are in analysis/self-test-funcs.R. Everything is
# cached per replicate, so n_reps can be raised later and only the new
# replicates are run. Monitor with:
#   tail -f data-generated/nwfsc-split-self-test/progress.log

source(here::here("analysis/self-test-funcs.R"))
source(here::here("analysis/fit-covariates.R"))
theme_set(ggsidekick::theme_sleek())

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
# Moving window (in survey occasions) for max A, as in the annual-survey
# scripts; 21 occasions gives 12 windows
seesaw_window <- 10L
# Coarsen the 2 km prediction grid to cut prediction cost
grid_res_km <- 6
# Predict the index in chunks of years to limit peak memory per worker
index_nsplit <- 4L
# Mesh cutoff (km); 20 gives ~400 nodes for the WCGBTS (15 gave ~600)
mesh_cutoff <- 20
# Peak memory was ~3.5 GB per worker with mesh_cutoff = 15, so lower this if
# the machine has less than ~3.5 GB x workers of RAM
n_workers <- max(1L, min(60L, future::availableCores() - 2L))

# Operating model, also refit to the full simulated data for calibration. IID
# spatiotemporal fields: with AR1, rho went to ~1 for several species, so each
# replicate's new epsilon_st was in effect a second, persistent spatial field
# drawn from the prior, swamping the posterior omega_s. With IID, persistent
# structure stays in omega_s and the new epsilon_st are year-to-year noise.
om_model <- "IID RF, factor(year)"
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
  add_utm_columns(c("lon", "lat"), utm_crs = utm_crs) |>
  mutate(X = round(X / grid_res_km) * grid_res_km, Y = round(Y / grid_res_km) * grid_res_km) |>
  group_by(region, X, Y) |>
  summarise(area = sum(area), .groups = "drop") |>
  # the region sampled in phase +0.5 occasions, for the true N - S contrast
  mutate(phase_region = as.integer(region == "North"))

# Data -----------------------------------------------------------------------

# Model time as survey occasion (no survey in 2020) and alternate the sampled
# half by occasion
prepare_species_data <- function(species) {
  surveyjoin::get_data(species, regions = "nwfsc") |>
    filter(.data$survey_name == .env$survey_name, effort > 0) |>
    mutate(year = lubridate::year(lubridate::ymd(date))) |>
    select(year, lon_start, lat_start, depth_m, effort, catch_weight) |>
    tidyr::drop_na() |>
    mutate(
      region = if_else(lat_start >= split_lat, "North", "South"),
      occasion = match(year, sort(unique(year))),
      primary_region = if_else(occasion %% 2L == 1L, "North", "South"),
      # seesaw phase follows which half was sampled
      phase = if_else(primary_region == "North", 0.5, -0.5),
      cal_year = year,
      year = occasion
    ) |>
    add_utm_columns(ll_names = c("lon_start", "lat_start"), utm_crs = utm_crs)
}

# Run ------------------------------------------------------------------------

setup_species <- function(species, dat) {
  is_sampled <- dat$region == dat$primary_region
  mesh <- make_mesh(dat, c("X", "Y"), cutoff = mesh_cutoff)
  setup <- sim_setup(
    species = species,
    dat = dat,
    is_sampled = is_sampled,
    grid = replicate_df(base_grid, "year", sort(unique(dat$year))),
    mesh = mesh,
    mesh_alt = make_mesh(dat[is_sampled, ], c("X", "Y"), mesh = mesh$mesh),
    occasions = distinct(dat, year, cal_year, primary_region, phase)
  )
  om <- fit_operating_model(om_model, dat, mesh)
  if (isTRUE(om$om_failed)) return(c(list(species = species), om))
  c(setup, list(fit_om = om$fit, ran_pars = om$ran_pars))
}

results <- run_self_test(
  dir = here::here("data-generated/nwfsc-split-self-test"),
  species = species_to_fit,
  n_reps = n_reps,
  prepare_data = prepare_species_data,
  setup_species = setup_species,
  simulate_rep = simulate_one_arm(
    if (resample_spatial) c("spatial", "spatiotemporal") else "spatiotemporal"
  ),
  models = tibble::tibble(
    design = c("full", rep("alternating", length(alternating_models))),
    model = c(om_model, alternating_models)
  ),
  n_workers = n_workers,
  index_nsplit = index_nsplit
)

# Metrics and plots ------------------------------------------------------------

m <- self_test_metrics(compare_to_truth(results), results, seesaw_window)
# Per-species fit summaries (random-effect parameters, regional omega / density,
# gradient) from the cached operating models, for relating the seesaw error to
# fit properties (analysis/18-seesaw-correlates.R)
m$fit_covariates <- cached_fit_covariates(
  here::here("data-generated/nwfsc-split-self-test"), species_to_fit,
  grid_group = "region", dat_group = "region"
)
saveRDS(m, here::here("data-generated/nwfsc-split-self-test-metrics.rds"))

m$metrics |>
  mutate(across(where(is.numeric), \(x) round(x, 2))) |>
  print(n = Inf, width = Inf)

plot_self_test(
  m,
  full_model = paste0(om_model, " (full data)"),
  window = seesaw_window,
  fig_prefix = here::here("figs/nwfsc-split-self-test"),
  caption = paste0(
    "WCGBTS alternating N / S self-test (split at ", round(split_lat, 2), " N), ",
    n_reps, " replicates per species, OM: ", om_model,
    if (resample_spatial) " (spatial fields new from prior). " else " (spatial fields posterior draws). ",
    "Dashed lines: nominal values (A, max A: median of the true index). ",
    "Seesaw error: estimated minus true phase effect, positive when the index reads high in occasions the richer region is sampled."
  )
)
