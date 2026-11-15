library(dplyr, warn.conflicts = FALSE)

# Data prep used by a single case study lives in its fit script

trawl_species <- c(
  "arrowtooth flounder",
  "dover sole",
  "english sole",
  "flathead sole",
  "lingcod",
  "longnose skate",
  "pacific cod",
  "pacific halibut",
  "pacific ocean perch",
  "pacific spiny dogfish",
  "petrale sole",
  "redbanded rockfish",
  "rex sole",
  "sablefish",
  "shortspine thornyhead",
  "slender sole",
  "spotted ratfish",
  "walleye pollock"
)

# The QCS + HS self-test species: strong to moderate seesaws under biennial
# QCS / HS sampling across a range of taxa, plus arrowtooth flounder as a
# control with almost none
qcs_hs_species <- c(
  "pacific cod",
  "spotted ratfish",
  "pacific ocean perch",
  "pacific halibut",
  "pacific spiny dogfish",
  "lingcod",
  "sablefish",
  "dover sole",
  "arrowtooth flounder",
  "shortspine thornyhead"
)

# QCS + HS ---------------------------------------------------------------------

qcs_hs_surveys <- c("SYN QCS", "SYN HS")

# Time is survey occasion: QCS + HS years are all odd, so on the calendar-year
# scale extra_time would put an empty year between every survey and AR1 rho
# gets stuck at its starting value of 0 (only rho^2 is identified). QCS was
# sampled alone in 2003 and 2004, so only years with both surveys are kept.
qcs_hs_data <- function(species) {
  dat <- surveyjoin::get_data(species, regions = "pbs") |>
    mutate(year = lubridate::year(lubridate::ymd(date))) |>
    filter(survey_name %in% qcs_hs_surveys, effort > 0) |>
    select(survey_name, year, lon_start, lat_start, effort, catch_weight) |>
    drop_na()

  complete_years <- dat |>
    distinct(year, survey_name) |>
    count(year) |>
    filter(n == length(qcs_hs_surveys)) |>
    pull(year)

  dat |>
    filter(year %in% complete_years) |>
    mutate(
      occasion = match(year, sort(unique(year))),
      primary_survey = if_else(occasion %% 2L == 1L, "SYN HS", "SYN QCS"),
      phase = if_else(primary_survey == "SYN HS", 0.5, -0.5),
      cal_year = year,
      year = occasion
    ) |>
    add_utm_columns(ll_names = c("lon_start", "lat_start"), utm_crs = 3156)
}

# One year; replicate it over the occasions
qcs_hs_grid <- function() {
  surveyjoin::dfo_synoptic_grid |>
    filter(survey %in% qcs_hs_surveys) |>
    add_utm_columns(c("lon", "lat"), utm_crs = 3156)
}

load_surveyjoin <- function() {
  surveyjoin::cache_data()
  surveyjoin::load_sql_data()
}
