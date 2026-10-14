# HBLL OUT (N and S) version of 03-fit-bc-models.R
# N and S alternate years, but the alternation flips calendar-year parity
# after the 2013 gap (N: 2006-2012 even, 2015-2025 odd), so the seesaw phase
# is coded by region rather than by year parity.

library(sdmTMB)
library(ggplot2)
library(dplyr)

source(here::here("analysis/fit-index-models.R"))
source(here::here("analysis/metric-functions.R"))
source(here::here("analysis/plot-index-models.R"))

hbll_surveys <- c("HBLL OUT N", "HBLL OUT S")

# Cache the HBLL OUT subset of the gfsynopsis-2025 data pull
sets_file <- here::here("data-raw/hbll-out-sets.rds")
if (!file.exists(sets_file)) {
  readRDS("~/src/gfsynopsis-2025/report/data-cache-2026-07/survey-sets.rds") |>
    filter(survey_abbrev %in% hbll_surveys) |>
    select(
      species_common_name, survey_abbrev, year, fishing_event_id,
      longitude, latitude, depth_m, lglsp_hook_count, count_bait_only, catch_count
    ) |>
    saveRDS(sets_file)
}
# Hook competition adjustment as in gfsynopsis-2025 (prep_stitch_dat()):
# offset = log(hooks / adjust), adjust = -log(p) / (1 - p), p = prop. baited hooks
sets <- readRDS(sets_file) |>
  tidyr::drop_na(lglsp_hook_count, count_bait_only, depth_m, catch_count) |>
  mutate(
    count_bait_only = replace(count_bait_only, count_bait_only == 0, 1),
    prop_bait_hooks = count_bait_only / lglsp_hook_count,
    hook_adjust_factor = -log(prop_bait_hooks) / (1 - prop_bait_hooks),
    offset = log(lglsp_hook_count / hook_adjust_factor)
  ) |>
  filter(is.finite(offset)) |>
  mutate(species_common_name = replace(
    species_common_name, species_common_name == "north pacific spiny dogfish", "pacific spiny dogfish"
  ))

grid0 <- bind_rows(
  gfplot::hbll_n_grid$grid |> mutate(survey_abbrev = "HBLL OUT N"),
  gfplot::hbll_s_grid$grid |> mutate(survey_abbrev = "HBLL OUT S")
) |>
  rename(lon = X, lat = Y, depth_m = depth) |>
  mutate(area = gfplot::hbll_n_grid$cell_area) |>
  sdmTMB::add_utm_columns(c("lon", "lat"), utm_crs = 3156)

do_fit_hbll <- function(.sp) {
  RhpcBLASctl::blas_set_num_threads(1L)
  RhpcBLASctl::omp_set_num_threads(1L)

  dat <- sets |>
    filter(species_common_name == .sp) |>
    sdmTMB::add_utm_columns(c("longitude", "latitude"), utm_crs = 3156)

  grid <- grid0 |>
    clamp_depth(dat) |>
    sdmTMB::replicate_df("year", sort(unique(dat$year)))

  mesh <- make_mesh(dat, c("X", "Y"), cutoff = 10)

  fit_index_models(
    dat = dat,
    grid = grid,
    mesh = mesh,
    response = "catch_count",
    family = list(nbinom1(), nbinom2()),
    offset = dat$offset
  ) |>
    mutate(species = .sp)
}

spp_to_fit_hbll <- c(
  "rougheye/blackspotted rockfish complex",
  "china rockfish",
  "copper rockfish",
  "pacific spiny dogfish",
  "tiger rockfish",
  "lingcod",
  "canary rockfish",
  "quillback rockfish",
  "yelloweye rockfish",
  "silvergray rockfish",
  "spotted ratfish",
  "big skate",
  "rosethorn rockfish",
  "southern rock sole",
  "longnose skate",
  "pacific cod",
  "arrowtooth flounder",
  "pacific halibut"
)
stopifnot(all(spp_to_fit_hbll %in% sets$species_common_name))

# Delete the .rds to force a refit
fits_file <- here::here("data-generated/hbll-indexes.rds")
if (!file.exists(fits_file)) {
  RhpcBLASctl::blas_set_num_threads(1L)
  RhpcBLASctl::omp_set_num_threads(1L)

  future::plan(future::multisession, workers = min(c(length(spp_to_fit_hbll), future::availableCores())))
  out <- furrr::future_map_dfr(spp_to_fit_hbll, do_fit_hbll, .options = furrr::furrr_options(seed = TRUE))
  future::plan(future::sequential)
  saveRDS(out, file = fits_file)
}

out <- readRDS(fits_file) |>
  filter(!grepl("depth", model), model != "Spatial only, RW year") |>
  mutate(species = replace(species, species == "north pacific spiny dogfish", "pacific spiny dogfish"))

# Which region was sampled in each year; +0.5 = N, -0.5 = S
lu <- sets |>
  distinct(year, survey_abbrev) |>
  mutate(phase = ifelse(survey_abbrev == "HBLL OUT N", 0.5, -0.5))
stopifnot(!any(duplicated(lu$year)))

seesaw_window <- 10L
seesaw_mw <- out |>
  left_join(lu, by = "year") |>
  arrange(species, model, year) |>
  group_by(species, model) |>
  group_modify(\(.x, .y) moving_window(.x$est, year = .x$year, window = seesaw_window, phase = .x$phase)) |>
  ungroup()

plot_indexes(out, lu, seesaw_mw, colour = "survey_abbrev", .ylab = "Abundance index")
ggsave(here::here("figs/hbll-testing.pdf"), width = 30, height = 15)

plot_A_moving_window(seesaw_mw, seesaw_window)
ggsave(here::here("figs/hbll-A-moving-window.pdf"), width = 5, height = 3.5)

plot_A_moving_window(seesaw_mw, seesaw_window, connect_stocks = TRUE, n_highlight = 5)
ggsave(here::here("figs/hbll-A-moving-window-connected.pdf"), width = 5, height = 3.5)

plot_top_stock_indexes(out, seesaw_mw, n_top = 5, models = c("IID RF, factor(year)", "IID RF, RW year", "RW RF"),
  lu = lu, .ylab = "Centered abundance index")
ggsave(here::here("figs/hbll-top-stock-indexes.pdf"), width = 6.2, height = 5)
