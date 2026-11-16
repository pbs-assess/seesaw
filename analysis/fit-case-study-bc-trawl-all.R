# Alternative to fit-case-study-bc-trawl.R: fit the BC synoptic trawl data as
# collected, without filtering to a cleanly biennial schedule. This keeps the
# QCS-only 2003-2004, the off-pattern WCHG 2007 and WCVI 2021, and the
# WCHG-only 2020.

library(dplyr)
source(here::here("analysis/R/models.R"))
source(here::here("analysis/R/survey-data.R"))
load_surveyjoin()

# Phase follows year parity (+0.5 = WCHG + WCVI years) as if the alternation
# had continued, so it is the same as in the biennial fits
bc_trawl_phase <- \(year) ifelse(year %% 2 == 0, 0.5, -0.5)

fit_bc_trawl_all <- function(species, model) {
  single_thread()
  dat_all <- surveyjoin::get_data(species, regions = "pbs") |>
    mutate(year = lubridate::year(lubridate::ymd(date))) |>
    select(survey_name, year, lon_start, lat_start, effort, catch_weight) |>
    add_utm_columns(ll_names = c("lon_start", "lat_start"), utm_crs = 3156) |>
    drop_na(effort, catch_weight)
  dat_biennial <- filter(
    dat_all,
    !(year %in% c(2003, 2004, 2020)),
    !(year == 2007 & survey_name == "SYN WCHG"),
    !(year == 2021 & survey_name == "SYN WCVI")
  )

  # Same mesh as the biennial fit so the two are comparable
  mesh_biennial <- make_mesh(dat_biennial, c("X", "Y"), cutoff = 10)
  mesh_all <- make_mesh(dat_all, c("X", "Y"), mesh = mesh_biennial$mesh)

  grid <- surveyjoin::dfo_synoptic_grid |>
    add_utm_columns(c("lon", "lat"), utm_crs = 3156) |>
    replicate_df("year", sort(unique(dat_all$year)))

  fit_index_models(
    models = model,
    data = dat_all, grid = grid, mesh = mesh_all,
    response = "catch_weight",
    family = delta_gamma(type = "poisson-link"),
    offset = log(dat_all$effort)
  ) |>
    mutate(phase = bc_trawl_phase(year))
}

run_case_study(
  trawl_species, fit_bc_trawl_all,
  file = here::here("data-generated/case-study-bc-trawl-all.rds"),
  n_workers = min(70L, parallelly::availableCores())
)
