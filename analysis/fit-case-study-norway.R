# Norwegian coastal gillnet-fyke net survey. North and south alternate years
# 2013-2023 (no 2014 survey), then both are sampled 2024-2026. The
# "alternating" design continues the alternation (2024 N, 2025 S, 2026 N);
# "all" uses all the data. Phase follows the alternating design in both.
#
# Crabs are dropped (counts only). Anisotropy is off: the sites lie along a
# narrow coastal strip, so it is poorly identified.

library(dplyr)
source(here::here("analysis/R/models.R"))

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
  drop_na(lon, lat, biomass_kg) |>
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

fit_norway <- function(species) {
  single_thread()
  dat_all <- filter(sets, englishname == species)
  dat_alt <- filter(dat_all, is_alternating(dat_all))
  mesh_alt <- make_mesh(dat_alt, c("X", "Y"), cutoff = 3)
  mesh_all <- make_mesh(dat_all, c("X", "Y"), mesh = mesh_alt$mesh)
  designs <- list(
    alternating = list(dat = dat_alt, mesh = mesh_alt),
    all = list(dat = dat_all, mesh = mesh_all)
  )
  map(designs, \(x) {
    fit_index_models(
      data = x$dat, mesh = x$mesh,
      grid = replicate_df(grid0, "year", sort(unique(x$dat$year))),
      response = "catch_weight",
      family = delta_gamma(type = "poisson-link"),
      covariates = "factor(gear)",
      anisotropy = FALSE
    )
  }) |>
    bind_rows(.id = "design") |>
    left_join(norway_phase, by = "year")
}

# Species with >= 10% positive sets and at least one positive set every year
# in the alternating design
norway_species <- sets |>
  filter(is_alternating(sets)) |>
  summarise(n_pos = sum(catch_weight > 0), n = n(), .by = c(englishname, year)) |>
  summarise(prop_pos = sum(n_pos) / sum(n), min_pos = min(n_pos), .by = englishname) |>
  filter(prop_pos >= 0.1, min_pos >= 1) |>
  pull(englishname) |>
  sort()

run_case_study(
  norway_species, fit_norway,
  file = here::here("data-generated/case-study-norway.rds"),
  n_workers = min(length(norway_species), parallelly::availableCores() / 2),
  designs = c("alternating", "all")
)
