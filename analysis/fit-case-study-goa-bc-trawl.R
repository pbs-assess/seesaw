# AFSC Gulf of Alaska and DFO synoptic trawl surveys fit together. GOA and
# QCS + HS are surveyed in odd years and WCHG + WCVI in even years. Unlike the
# BC trawl case study, the irregular years are kept.

library(dplyr)
source(here::here("analysis/R/models.R"))
source(here::here("analysis/R/survey-data.R"))
load_surveyjoin()

# +0.5 = WCHG + WCVI years
goa_bc_trawl_phase <- \(year) ifelse(year %% 2 == 0, 0.5, -0.5)

fit_goa_bc_trawl <- function(species) {
  single_thread()
  dat <- surveyjoin::get_data(species, regions = c("pbs", "afsc")) |>
    mutate(year = lubridate::year(lubridate::ymd(date))) |>
    select(survey_name, year, lon_start, lat_start, effort, catch_weight) |>
    drop_na(lon_start, lat_start, effort, catch_weight) |>
    filter(
      survey_name %in% c("Gulf of Alaska", "SYN HS", "SYN QCS", "SYN WCHG", "SYN WCVI"),
      year >= 2003 # PBS surveys start
    ) |>
    mutate(survey_group = ifelse(grepl("syn", survey_name, ignore.case = TRUE), "pbs", "goa")) |>
    add_utm_columns(ll_names = c("lon_start", "lat_start"), utm_crs = 3156)

  prop_positive <- summarise(dat, prop_positive = mean(catch_weight > 0), .by = survey_group)
  if (!all(c("goa", "pbs") %in% prop_positive$survey_group) || any(prop_positive$prop_positive < 0.25)) {
    return(NULL)
  }

  grid <- surveyjoin::dfo_synoptic_grid |>
    bind_rows(filter(surveyjoin::afsc_grid, survey == "Gulf of Alaska Bottom Trawl Survey")) |>
    add_utm_columns(c("lon", "lat"), utm_crs = 3156) |>
    mutate(survey_group = "pbs") |> # index at the PBS catchability
    replicate_df("year", sort(unique(dat$year)))
  # coarser than the BC-only fits (10 km) for the much larger domain; 20 km
  # keeps ~200 nodes over BC (40 km gave ~80)
  mesh <- make_mesh(dat, c("X", "Y"), cutoff = 20)

  fit_index_models(
    data = dat, grid = grid, mesh = mesh,
    response = "catch_weight",
    family = delta_gamma(type = "poisson-link"),
    offset = log(dat$effort),
    covariates = "factor(survey_group)",
    nsplit = 6L # large grid
  ) |>
    mutate(phase = goa_bc_trawl_phase(year))
}

run_case_study(
  trawl_species, fit_goa_bc_trawl,
  file = here::here("data-generated/case-study-goa-bc-trawl.rds"),
  n_workers = min(length(trawl_species), parallelly::availableCores() / 2)
)
