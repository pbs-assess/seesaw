# Norwegian coastal net survey (garn-ruse tokt) version of 03-fit-syn-models.R
# Based on Brian Stock's analysis/03-norway-testing.R.
#
# North and south alternate years 2013-2023 (no 2014 survey), then both are
# sampled 2024-2026. The main models use an alternating subset that continues
# the pattern (2024 N, 2025 S, 2026 N); the "all data" models add the
# overlapping coverage. As with HBLL, the seesaw phase is coded by region.
#
# Biomass by gear set; effort is one set, and gear is a covariate.
# Crabs are dropped (counts only).

library(sdmTMB)
library(ggplot2)
library(dplyr)

source(here::here("analysis/fit-index-models.R"))
source(here::here("analysis/metric-functions.R"))
source(here::here("analysis/plot-index-models.R"))

spp_lu <- readxl::read_excel(here::here("data-raw/norway-catchsummary.xlsx")) |>
  select(commonname, englishname)

grid0 <- readRDS(here::here("data-raw/garn_ruse_predgrid.rds")) |>
  ungroup() |>
  mutate(gear = as.character(gear))
site_lu <- distinct(grid0, site, subarea)
stopifnot(!any(duplicated(site_lu$site)))

# Keep the sites in the prediction grid; this drops sites north of ~67 N that
# were only sampled in earlier years
sets <- readRDS(here::here("data-raw/garn_ruse_catch_21sp-2026.rds")) |>
  filter(!commonname %in% c("taskekrabbe", "strandkrabbe")) |>
  inner_join(site_lu, by = "site") |>
  tidyr::drop_na(lon, lat, depth, biomass_kg) |>
  rename(depth_m = depth, catch_weight = biomass_kg) |>
  left_join(spp_lu, by = "commonname") |>
  sdmTMB::add_utm_columns(c("lon", "lat"), utm_crs = 32632, units = "km")
stopifnot(!any(is.na(sets$englishname)))

# Alternating design: the region sampled in each year 2024+ continues the
# 2013-2023 alternation
alt_region <- tibble(year = 2024:2026, subarea = c("North", "South", "North"))
is_alternating <- function(d) {
  d$year < 2024 | paste(d$year, d$subarea) %in% paste(alt_region$year, alt_region$subarea)
}

# Some grid cells are missing depth (only used by the depth model); fill with
# the mean sampled depth for that site and gear
grid0 <- grid0 |>
  left_join(
    summarise(distinct(sets, haulid, site, gear, depth_m), mean_depth = mean(depth_m), .by = c(site, gear)),
    by = c("site", "gear")
  ) |>
  mutate(depth_m = ifelse(is.na(depth), mean_depth, depth)) |>
  select(-mean_depth)

do_fit_norway <- function(.sp) {
  RhpcBLASctl::blas_set_num_threads(1L)
  RhpcBLASctl::omp_set_num_threads(1L)

  dat_all <- sets |>
    filter(englishname == .sp) |>
    mutate(survey_name = subarea)
  dat <- filter(dat_all, is_alternating(dat_all))

  grid <- grid0 |>
    clamp_depth(dat) |>
    sdmTMB::replicate_df("year", sort(unique(dat$year)))

  mesh <- sdmTMB::make_mesh(dat, c("X", "Y"), cutoff = 3)
  mesh_all <- sdmTMB::make_mesh(dat_all, c("X", "Y"), mesh = mesh$mesh)

  fit_index_models(
    dat = dat,
    grid = grid,
    mesh = mesh,
    response = "catch_weight",
    family = sdmTMB::delta_gamma(type = "poisson-link"),
    offset = rep(0, nrow(dat)),
    covariates = "factor(gear)",
    all_data = list(data = dat_all, mesh = mesh_all, offset = rep(0, nrow(dat_all))),
    covariate_groups = c(grid = "subarea", dat = "subarea"),
    # e.g., the cod binomial spatiotemporal SD goes to ~0; drop such fields
    control = sdmTMB::sdmTMBcontrol(collapse_spatial_variance = TRUE),
    # Sites lie along a narrow coastal strip, so anisotropy is poorly identified
    anisotropy = FALSE
  ) |>
    with_species(.sp)
}

# Species with >= 10% positive sets and at least one positive set every year
# in the alternating design
spp_to_fit_norway <- sets |>
  filter(is_alternating(sets)) |>
  summarise(n_pos = sum(catch_weight > 0), n = n(), .by = c(englishname, year)) |>
  summarise(prop_pos = sum(n_pos) / sum(n), min_pos = min(n_pos), .by = englishname) |>
  filter(prop_pos >= 0.1, min_pos >= 1) |>
  pull(englishname) |>
  sort()

# Delete the .rds files to force a refit (the covariates file has per-species
# fit summaries; see 03-fit-syn-models.R)
fits_file <- here::here("data-generated/norway-fit-indexes.rds")
covariates_file <- here::here("data-generated/norway-fit-covariates.rds")
if (!file.exists(fits_file) || !file.exists(covariates_file)) {
  RhpcBLASctl::blas_set_num_threads(1L)
  RhpcBLASctl::omp_set_num_threads(1L)

  future::plan(future::multisession, workers = min(c(length(spp_to_fit_norway), future::availableCores() / 2)))
  fitted <- furrr::future_map(setNames(spp_to_fit_norway, spp_to_fit_norway), do_fit_norway, .options = furrr::furrr_options(seed = TRUE))
  future::plan(future::sequential)
  saveRDS(dplyr::bind_rows(lapply(fitted, `[[`, "index")), file = fits_file)
  saveRDS(purrr::compact(lapply(fitted, `[[`, "covariates")), file = covariates_file)
}

out <- readRDS(fits_file) |>
  filter(!grepl("depth", model), model != "Spatial only, RW year")

# Which region was sampled in each year of the alternating design; +0.5 = N, -0.5 = S
lu <- sets |>
  filter(is_alternating(sets)) |>
  distinct(year, subarea) |>
  mutate(phase = ifelse(subarea == "North", 0.5, -0.5))
stopifnot(!any(duplicated(lu$year)))

seesaw_window <- 10L
seesaw_mw <- out |>
  left_join(lu, by = "year") |>
  arrange(species, model, year) |>
  group_by(species, model) |>
  group_modify(\(.x, .y) moving_window(.x$est, year = .x$year, window = seesaw_window, phase = .x$phase)) |>
  ungroup()

plot_indexes(out, lu, seesaw_mw, colour = "subarea", colour_lab = "Region\n(alternating)", .ylab = "Biomass index")
ggsave(here::here("figs/norway-index-models.pdf"), width = 30, height = 15)

plot_A_moving_window(seesaw_mw, seesaw_window, include_all_data = TRUE)
ggsave(here::here("figs/norway-A-moving-window.pdf"), width = 5.5, height = 3.5)

plot_A_moving_window(seesaw_mw, seesaw_window, connect_stocks = TRUE, n_highlight = 5)
ggsave(here::here("figs/norway-A-moving-window-connected.pdf"), width = 5, height = 3.5)

plot_top_stock_indexes(out, seesaw_mw, n_top = 5, models = c("IID RF, factor(year)", "IID RF, RW year", "RW RF"),
  lu = lu, .ylab = "Centered biomass index")
ggsave(here::here("figs/norway-top-stock-indexes.pdf"), width = 6.2, height = 5)
