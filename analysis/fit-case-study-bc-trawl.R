# BC synoptic trawl. QCS + HS are surveyed in odd years and WCHG + WCVI in even
# years.

library(dplyr)
source(here::here("analysis/R/models.R"))
source(here::here("analysis/R/survey-data.R"))
load_surveyjoin()

# +0.5 = WCHG + WCVI years
bc_trawl_phase <- \(year) ifelse(year %% 2 == 0, 0.5, -0.5)

fit_bc_trawl <- function(species) {
  single_thread()
  dat <- surveyjoin::get_data(species, regions = "pbs") |>
    mutate(year = lubridate::year(lubridate::ymd(date))) |>
    select(survey_name, year, lon_start, lat_start, effort, catch_weight) |>
    add_utm_columns(ll_names = c("lon_start", "lat_start"), utm_crs = 3156) |>
    drop_na(effort, catch_weight) |>
    # Keep a cleanly biennial design: drop the years without both halves
    # (2003-2004, QCS only; 2020, no survey) and the off-pattern surveys
    filter(
      !(year %in% c(2003, 2004, 2020)),
      !(year == 2007 & survey_name == "SYN WCHG"),
      !(year == 2021 & survey_name == "SYN WCVI")
    )

  grid <- surveyjoin::dfo_synoptic_grid |>
    add_utm_columns(c("lon", "lat"), utm_crs = 3156) |>
    replicate_df("year", sort(unique(dat$year)))
  mesh <- make_mesh(dat, c("X", "Y"), cutoff = 10)

  fit_index_models(
    data = dat, grid = grid, mesh = mesh,
    response = "catch_weight",
    family = delta_gamma(type = "poisson-link"),
    offset = log(dat$effort)
  ) |>
    mutate(phase = bc_trawl_phase(year))
}

run_case_study(
  trawl_species, fit_bc_trawl,
  file = here::here("data-generated/case-study-bc-trawl.rds"),
  n_workers = min(length(trawl_species), 60L, parallelly::availableCores())
)
