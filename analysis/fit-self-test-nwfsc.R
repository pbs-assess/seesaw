# NWFSC West Coast groundfish bottom trawl survey (WCGBTS) self-test. The
# split is artificial: the survey samples the whole coast every year, so the
# coast is cut at the median latitude of the tows (fixed by the design before
# looking at any results) and one half is kept per occasion.
# This gives a second, independent geometry (a long, narrow, along-coast
# domain like the BC surveys) with the same species as far as possible.
#
# Monitor with: tail -f data-generated/self-test-nwfsc/progress.log

library(dplyr)
source(here::here("analysis/R/self-test.R"))
source(here::here("analysis/R/survey-data.R"))
load_surveyjoin()

survey_name <- "NWFSC.Combo"
utm_crs <- 32610 # UTM 10N
# Present in both halves (pacific cod, pacific ocean perch, and pacific
# halibut are nearly absent south of the split); spans a range of N-S
# contrasts and overlaps the QCS + HS self-test species where possible
nwfsc_species <- c(
  "spotted ratfish",
  "pacific spiny dogfish",
  "lingcod",
  "shortspine thornyhead",
  "arrowtooth flounder",
  "sablefish",
  "dover sole",
  "petrale sole"
)
n_reps <- 30L
grid_res_km <- 6 # coarsened from 2 km to cut prediction cost
index_nsplit <- 6L # to limit peak memory per worker
mesh_cutoff <- 15 # ~600 nodes
# Peak memory was ~3.5 GB per worker with mesh_cutoff = 15
n_workers <- max(1L, min(40L, parallelly::availableCores() - 2L))

# The operating model has IID spatiotemporal fields; with AR1, rho went to ~1
# for several species here (see fit-self-test-qcs-hs.R)
fits <- bind_rows(
  tibble(design = "full", model = om_model),
  tibble(design = "alternating", model = names(paper_models))
)

# Grid and split ------------------------------------------------------------

get_nwfsc_data <- function(species) {
  surveyjoin::get_data(species, regions = "nwfsc") |>
    filter(.data$survey_name == .env$survey_name, effort > 0) |>
    mutate(year = lubridate::year(lubridate::ymd(date))) |>
    select(date, year, lon_start, lat_start, effort, catch_weight) |>
    drop_na()
}

# The median tow latitude gives each half about the same number of tows per
# occasion (~310 North, ~335 South), so the alternation is in region only.
# Splitting at the area-weighted median of the grid (37.7 N) instead gave ~400
# North and ~240 South, so sample size alternated with region. The halves
# then differ in area (North ~38%), as the BC survey regions do.
split_lat <- map(nwfsc_species, get_nwfsc_data) |>
  bind_rows() |>
  distinct(date, lon_start, lat_start) |>
  pull(lat_start) |>
  median()
cli::cli_inform("Splitting north / south at {round(split_lat, 2)} degrees N")

grid_fine <- surveyjoin::nwfsc_grid |>
  filter(survey == survey_name, area > 0)

# Assign the half at 2 km, then aggregate within each half so no coarse cell
# straddles the split
base_grid <- grid_fine |>
  mutate(region = if_else(lat >= split_lat, "North", "South")) |>
  add_utm_columns(c("lon", "lat"), utm_crs = utm_crs) |>
  mutate(X = round(X / grid_res_km) * grid_res_km, Y = round(Y / grid_res_km) * grid_res_km) |>
  summarise(area = sum(area), .by = c(region, X, Y))

# Data -----------------------------------------------------------------------

# Time is survey occasion because there was no survey in 2020
prepare_species_data <- function(species) {
  get_nwfsc_data(species) |>
    select(-date) |>
    mutate(
      region = if_else(lat_start >= split_lat, "North", "South"),
      occasion = match(year, sort(unique(year))),
      primary_region = if_else(occasion %% 2L == 1L, "North", "South"),
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
    occasions = distinct(dat, year, cal_year, phase)
  )
  om <- fit_operating_model(dat, mesh)
  if (isTRUE(om$om_failed)) return(c(list(species = species), om))
  c(setup, list(fit_om = om$fit))
}

results <- run_self_test(
  dir = here::here("data-generated/self-test-nwfsc"),
  species = nwfsc_species,
  n_reps = n_reps,
  prepare_data = prepare_species_data,
  setup_species = setup_species,
  fits = fits,
  n_workers = n_workers,
  index_nsplit = index_nsplit
)
saveRDS(results, here::here("data-generated/self-test-nwfsc.rds"))
