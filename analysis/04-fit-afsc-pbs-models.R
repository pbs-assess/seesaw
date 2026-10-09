# Transboundary index models: DFO synoptic trawl (PBS) + AFSC Gulf of Alaska
# Same setup as 03-fit-syn-models.R; assess the seesawness of the various
# index models.

library(sdmTMB)
library(ggplot2)
library(dplyr)

surveyjoin::cache_data()
surveyjoin::load_sql_data()

source(here::here("analysis/fit-index-models.R"))
source(here::here("analysis/metric-functions.R"))
source(here::here("analysis/plot-index-models.R"))

do_fit_afsc_pbs <- function(.sp) {
  RhpcBLASctl::blas_set_num_threads(1L)
  RhpcBLASctl::omp_set_num_threads(1L)

  dat <- surveyjoin::get_data(.sp, regions = c("pbs", "afsc")) |>
    dplyr::mutate(year = lubridate::year(lubridate::ymd(date))) |>
    dplyr::select(survey_name, year, lon_start, lat_start, depth_m, effort, catch_weight, common_name) |>
    tidyr::drop_na(lon_start, lat_start) |>
    dplyr::filter(survey_name %in% c("Gulf of Alaska", "SYN HS", "SYN QCS", "SYN WCHG", "SYN WCVI")) |>
    dplyr::filter(year >= 2003) |> # PBS surveys start
    tidyr::drop_na(effort, catch_weight, depth_m) |>
    dplyr::mutate(survey_group = ifelse(grepl("syn", survey_name, ignore.case = TRUE), "pbs", "goa"))

  dat <- sdmTMB::add_utm_columns(dat, ll_names = c("lon_start", "lat_start"), utm_crs = 3156)

  # Require enough positive sets in both the GOA and PBS surveys
  prop_positive <- dat |>
    dplyr::summarise(prop_positive = mean(catch_weight > 0), .by = survey_group)
  if (!all(c("goa", "pbs") %in% prop_positive$survey_group) || any(prop_positive$prop_positive < 0.25)) {
    return(dplyr::tibble())
  }

  grid <- surveyjoin::dfo_synoptic_grid |>
    dplyr::bind_rows(dplyr::filter(surveyjoin::afsc_grid, survey == "Gulf of Alaska Bottom Trawl Survey")) |>
    sdmTMB::add_utm_columns(c("lon", "lat"), utm_crs = 3156) |>
    clamp_depth(dat) |>
    dplyr::mutate(survey_group = "pbs") |> # index is for the PBS catchability level
    sdmTMB::replicate_df("year", sort(unique(dat$year)))

  mesh <- sdmTMB::make_mesh(dat, c("X", "Y"), cutoff = 40)

  fit_index_models(
    dat = dat,
    grid = grid,
    mesh = mesh,
    response = "catch_weight",
    family = sdmTMB::delta_gamma(type = "poisson-link"),
    offset = log(dat$effort),
    covariates = "factor(survey_group)"
  ) |>
    dplyr::mutate(species = .sp)
}

spp_to_fit_afsc_pbs <- c(
  "arrowtooth flounder",
  "dover sole",
  "english sole",
  "flathead sole",
  # "greenstriped rockfish",
  "lingcod",
  "longnose skate",
  # "north pacific hake",
  "pacific cod",
  "pacific halibut",
  "pacific ocean perch",
  "pacific spiny dogfish",
  "petrale sole",
  "redbanded rockfish",
  "rex sole",
  "sablefish",
  # "sharpchin rockfish",
  "shortspine thornyhead",
  "slender sole",
  "spotted ratfish",
  "walleye pollock"
)

# Delete the .rds to force a refit
fits_file <- here::here("data-generated/transboundary-afsc-pbs-indexes.rds")
if (!file.exists(fits_file)) {
  RhpcBLASctl::blas_set_num_threads(1L)
  RhpcBLASctl::omp_set_num_threads(1L)

  future::plan(future::multisession, workers = min(c(length(spp_to_fit_afsc_pbs), future::availableCores() / 2)))
  out <- furrr::future_map_dfr(spp_to_fit_afsc_pbs, do_fit_afsc_pbs, .options = furrr::furrr_options(seed = TRUE))
  future::plan(future::sequential)
  saveRDS(out, file = fits_file)
}

out <- readRDS(fits_file) |>
  filter(!grepl("depth", model), model != "Spatial only, RW year")

# Phase is coded by calendar year: even vs odd
lu <- distinct(out, year) |>
  mutate(survey_group = ifelse(year %% 2 == 0, "Even years", "Odd years"))

seesaw_window <- 10L
seesaw_mw <- out |>
  arrange(species, model, year) |>
  group_by(species, model) |>
  group_modify(\(.x, .y) moving_window(.x$est, year = .x$year, window = seesaw_window)) |>
  ungroup()

plot_indexes(out, lu, seesaw_mw, colour = "survey_group", colour_lab = "Survey\ngrouping", .ylab = "Biomass index")
ggsave(here::here("figs/transboundary-testing-afsc-pbs.pdf"), width = 30, height = 15)

out |>
  group_by(model) |>
  summarise(n = n())

plot_A_moving_window(seesaw_mw, seesaw_window, connect_stocks = TRUE, n_highlight = 5)
ggsave(here::here("figs/transboundary-trawl-A-moving-window-afsc-pbs.pdf"), width = 5, height = 3.5)

plot_top_stock_indexes(out, seesaw_mw, n_top = 5, models = c("IID RF, factor(year)", "IID RF, RW year", "RW RF"))
ggsave(here::here("figs/transboundary-trawl-top-stock-indexes-afsc-pbs.pdf"), width = 6.2, height = 5)
