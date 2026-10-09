# 2026-06-26
# Running across 21 synoptic groundfish
# Could do for HBLL OUT too
# And maybe HBLL inside??
# Assess the seesawness for the various models
# Could probably also look at correlates in this case (e.g., estimate things like north-south gradient or difference in mean density by subregion)
#
# Then do this for the Norwegian survey?
# Can summarize the overall findings and show some case study panels in paper
#
# This could be more informative than the simulation study on which models work... although only that can get at the coverage and RMSE

library(sdmTMB)
library(ggplot2)
library(dplyr)

surveyjoin::cache_data()
surveyjoin::load_sql_data()

source(here::here("analysis/fit-index-models.R"))
source(here::here("analysis/metric-functions.R"))
source(here::here("analysis/plot-index-models.R"))

do_fit_syn <- function(.sp) {
  RhpcBLASctl::blas_set_num_threads(1L)
  RhpcBLASctl::omp_set_num_threads(1L)

  dat0 <- surveyjoin::get_data(.sp, regions = "pbs") |>
    dplyr::mutate(year = lubridate::year(lubridate::ymd(date))) |>
    dplyr::select(survey_name, year, lon_start, lat_start, depth_m, effort, catch_weight, common_name)

  dat0 <- sdmTMB::add_utm_columns(dat0, ll_names = c("lon_start", "lat_start"), utm_crs = 3156)
  dat <- dat0 |>
    # Use only complete N/S sampling years
    dplyr::filter(!(year %in% c(2003, 2004, 2020))) |>
    tidyr::drop_na(effort, catch_weight, depth_m) |>
    # Drop these surveys to be perfectly bienniel
    dplyr::filter(!(year == 2007 & survey_name == "SYN WCHG")) |>
    dplyr::filter(!(year == 2021 & survey_name == "SYN WCVI"))

  dat_all <- dat0 |>
    tidyr::drop_na(effort, catch_weight, depth_m)

  grid <- surveyjoin::dfo_synoptic_grid |>
    sdmTMB::add_utm_columns(c("lon", "lat"), utm_crs = 3156) |>
    dplyr::mutate(survey_name = "SYN WCVI") |>
    clamp_depth(dat) |>
    sdmTMB::replicate_df("year", sort(unique(dat$year)))

  mesh <- sdmTMB::make_mesh(dat, c("X", "Y"), cutoff = 10)
  mesh_all <- sdmTMB::make_mesh(dat_all, c("X", "Y"), mesh = mesh$mesh)

  fit_index_models(
    dat = dat,
    grid = grid,
    mesh = mesh,
    response = "catch_weight",
    family = sdmTMB::delta_gamma(type = "poisson-link"),
    offset = log(dat$effort),
    all_data = list(data = dat_all, mesh = mesh_all, offset = log(dat_all$effort))
  ) |>
    dplyr::mutate(species = .sp)
}

spp_to_fit_syn <- c(
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

RhpcBLASctl::blas_set_num_threads(1L)
RhpcBLASctl::omp_set_num_threads(1L)

future::plan(future::multisession, workers = min(c(length(spp_to_fit_syn), future::availableCores() / 2)))
out <- furrr::future_map_dfr(spp_to_fit_syn, do_fit_syn, .options = furrr::furrr_options(seed = TRUE))
future::plan(future::sequential)
saveRDS(out, file = here::here("data-generated/bc-indexes4.rds"))
#
out <- readRDS(here::here("data-generated/bc-indexes4.rds")) |>
  filter(!species %in% c("redbanded rockfish", "shortspine thornyhead"))

# Phase is coded by calendar year: even = WCHG + WCVI, odd = QCS + HS
lu <- distinct(out, year) |>
  mutate(survey_group = ifelse(year %% 2 == 0, "WCHG + WCVI", "QCS + HS"))

seesaw_window <- 10L
seesaw_mw <- out |>
  arrange(species, model, year) |>
  group_by(species, model) |>
  group_modify(\(.x, .y) moving_window(.x$est, year = .x$year, window = seesaw_window)) |>
  ungroup()

plot_indexes(out, lu, seesaw_mw, colour = "survey_group", colour_lab = "Survey\ngrouping", .ylab = "Biomass index")
ggsave(here::here("figs/bc-testing2.pdf"), width = 30, height = 15)

out |>
  group_by(model) |>
  summarise(n = n())

plot_A_moving_window(seesaw_mw, seesaw_window)
ggsave(here::here("figs/bc-trawl-A-moving-window.pdf"), width = 5, height = 3.5)

plot_A_moving_window(seesaw_mw, seesaw_window, include_all_data = TRUE)
ggsave(here::here("figs/bc-trawl-A-moving-window2.pdf"), width = 5.5, height = 3.5)
