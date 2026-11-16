# Norwegian coastal gillnet-fyke net survey. North and south alternate years
# 2013-2023 (no 2014 survey), then both are sampled 2024-2026. The
# "alternating" design continues the alternation (2024 N, 2025 S, 2026 N).
#
# Crabs are dropped (counts only).

library(dplyr)
library(sdmTMB)

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
  tidyr::drop_na(lon, lat, biomass_kg) |>
  rename(catch_weight = biomass_kg) |>
  left_join(spp_lu, by = "commonname") |>
  add_utm_columns(c("lon", "lat"), utm_crs = 32632, units = "km")
stopifnot(!any(is.na(sets$englishname)))

alt_region <- tibble(year = 2024:2026, subarea = c("North", "South", "North"))
is_alternating <- function(d) {
  d$year < 2024 | paste(d$year, d$subarea) %in% paste(alt_region$year, alt_region$subarea)
}
norway_phase <- sets |>
  filter(is_alternating(sets)) |>
  distinct(year, subarea) |>
  mutate(phase = ifelse(subarea == "North", 0.5, -0.5)) |>
  select(year, phase)
stopifnot(!any(duplicated(norway_phase$year)))

# Anisotropy is off: the sites lie along a narrow coastal strip, so it is
# poorly identified
fit_norway_model <- function(model, dat, mesh, ...) {
  fit_model(
    model,
    data = dat, mesh = mesh,
    response = "catch_weight",
    family = delta_gamma(type = "poisson-link"),
    covariates = "factor(gear)",
    anisotropy = FALSE,
    ...
  )
}
