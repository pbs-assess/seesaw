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

# design = "alternating": the main models on the alternating subset, plus the
# "all data" factor(year) models; design = "all": the main models on all the
# data (including both regions in 2024-2026), with the same mesh knots
do_fit_norway <- function(.sp, design = c("alternating", "all")) {
  design <- match.arg(design)
  RhpcBLASctl::blas_set_num_threads(1L)
  RhpcBLASctl::omp_set_num_threads(1L)

  dat_all <- sets |>
    filter(englishname == .sp) |>
    mutate(survey_name = subarea)
  dat_alt <- filter(dat_all, is_alternating(dat_all))
  mesh_alt <- sdmTMB::make_mesh(dat_alt, c("X", "Y"), cutoff = 3)
  mesh_all <- sdmTMB::make_mesh(dat_all, c("X", "Y"), mesh = mesh_alt$mesh)

  dat <- if (design == "alternating") dat_alt else dat_all
  mesh <- if (design == "alternating") mesh_alt else mesh_all
  grid <- grid0 |>
    clamp_depth(dat) |>
    sdmTMB::replicate_df("year", sort(unique(dat$year)))

  fit_index_models(
    dat = dat,
    grid = grid,
    mesh = mesh,
    response = "catch_weight",
    family = sdmTMB::delta_gamma(type = "poisson-link"),
    offset = rep(0, nrow(dat)),
    covariates = "factor(gear)",
    all_data = if (design == "alternating") list(data = dat_all, mesh = mesh_all, offset = rep(0, nrow(dat_all))),
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

# Which fits the main figures use: "alternating" (the alternating subset) or
# "all" (all the data, including both regions in 2024-2026)
index_design <- "alternating"
fig_suffix <- if (index_design == "all") "-all-data" else ""

# Delete the .rds files to force a refit (the covariates file has per-species
# fit summaries; see 03-fit-syn-models.R)
fits_files <- c(
  alternating = here::here("data-generated/norway-fit-indexes.rds"),
  all = here::here("data-generated/norway-fit-indexes-all-data.rds")
)
covariates_files <- c(
  alternating = here::here("data-generated/norway-fit-covariates.rds"),
  all = here::here("data-generated/norway-fit-covariates-all-data.rds")
)
for (design in names(fits_files)) {
  if (file.exists(fits_files[[design]]) && file.exists(covariates_files[[design]])) next
  RhpcBLASctl::blas_set_num_threads(1L)
  RhpcBLASctl::omp_set_num_threads(1L)

  future::plan(future::multisession, workers = min(c(length(spp_to_fit_norway), future::availableCores() / 2)))
  fitted <- furrr::future_map(
    setNames(spp_to_fit_norway, spp_to_fit_norway), do_fit_norway, design = design,
    .options = furrr::furrr_options(seed = TRUE)
  )
  future::plan(future::sequential)
  saveRDS(dplyr::bind_rows(lapply(fitted, `[[`, "index")), file = fits_files[[design]])
  saveRDS(purrr::compact(lapply(fitted, `[[`, "covariates")), file = covariates_files[[design]])
}

read_norway <- function(design) {
  readRDS(fits_files[[design]]) |>
    filter(!grepl("depth", model), model != "Spatial only, RW year")
}
out <- read_norway(index_design)

# Which region was sampled in each year of the alternating design; +0.5 = N, -0.5 = S
lu <- sets |>
  filter(is_alternating(sets)) |>
  distinct(year, subarea) |>
  mutate(phase = ifelse(subarea == "North", 0.5, -0.5))
stopifnot(!any(duplicated(lu$year)))

seesaw_window <- 10L
calc_mw <- function(out) {
  out |>
    left_join(lu, by = "year") |>
    arrange(species, model, year) |>
    group_by(species, model) |>
    group_modify(\(.x, .y) moving_window(.x$est, year = .x$year, window = seesaw_window, phase = .x$phase)) |>
    ungroup()
}
seesaw_mw <- calc_mw(out)

plot_indexes(out, lu, seesaw_mw, colour = "subarea", colour_lab = "Region\n(alternating)", .ylab = "Biomass index")
ggsave(here::here(paste0("figs/norway-index-models", fig_suffix, ".pdf")), width = 30, height = 15)

plot_A_moving_window(seesaw_mw, seesaw_window, include_all_data = TRUE)
ggsave(here::here(paste0("figs/norway-A-moving-window", fig_suffix, ".pdf")), width = 5.5, height = 3.5)

plot_A_moving_window(seesaw_mw, seesaw_window, connect_stocks = TRUE, n_highlight = 5)
ggsave(here::here(paste0("figs/norway-A-moving-window-connected", fig_suffix, ".pdf")), width = 5, height = 3.5)

plot_top_stock_indexes(out, seesaw_mw, n_top = 5, models = c("IID RF, factor(year)", "IID RF, RW year", "RW RF"),
  lu = lu, .ylab = "Centered biomass index")
ggsave(here::here(paste0("figs/norway-top-stock-indexes", fig_suffix, ".pdf")), width = 6.2, height = 5)

# How much does adding full coverage in 2024-2026 shrink the seesaw? A over
# the whole series and max A over the moving windows (only windows starting in
# 2017 contain all three full years). Phase is still the region the
# alternating design samples. The "with full" fits are from the "all" design
# (its IID RF, factor(year) is the same model as the alternating design's
# "IID RF, factor(year), all data").
full_years_data <- c("Alternating", "With full 2024-2026")
full_years_out <- bind_rows(
  read_norway("alternating") |>
    filter(!grepl("all data", model)) |>
    mutate(data = full_years_data[1]),
  read_norway("all") |>
    mutate(data = full_years_data[2])
)
A_full_years <- full_years_out |>
  left_join(lu, by = "year") |>
  arrange(species, data, model, year) |>
  group_by(species, data, model) |>
  group_modify(\(.x, .y) {
    tibble(
      A = period2_metric(.x$est, year = .x$year, phase = .x$phase)[["A"]],
      A_max = max(moving_window(.x$est, year = .x$year, window = seesaw_window, phase = .x$phase)$A)
    )
  }) |>
  ungroup() |>
  tidyr::pivot_longer(c(A, A_max), names_to = "metric") |>
  mutate(
    data = factor(data, levels = full_years_data),
    metric = factor(metric, levels = c("A", "A_max"), labels = c("A, 2013-2026", paste0("Max A, ", seesaw_window, "-year windows")))
  )
species_order <- A_full_years |>
  filter(data == "Alternating", model == "IID RF, factor(year)", metric == levels(metric)[1]) |>
  arrange(value) |>
  pull(species)
# Species with no alternating factor(year) fit (lemon sole, which only
# converges with the full years) go at the bottom
species_order <- c(setdiff(sort(unique(A_full_years$species)), species_order), species_order)

plot_full_years <- function(d) {
  d |>
    mutate(species = factor(species, levels = species_order)) |>
    ggplot(aes(value, species)) +
    geom_line(aes(group = species), colour = "grey70") +
    geom_point(aes(colour = data)) +
    scale_x_sqrt(breaks = c(0, 50, 200, 500, 1000, 2000)) +
    labs(x = "Biennial amplitude (%)", y = NULL, colour = NULL) +
    ggsidekick::theme_sleek() +
    theme(legend.position = "bottom")
}

A_full_years |>
  filter(model == "IID RF, factor(year)") |>
  plot_full_years() +
  facet_wrap(~metric, scales = "free_x") +
  ggtitle("IID RF, factor(year)")
ggsave(here::here("figs/norway-A-full-years.pdf"), width = 6.5, height = 4)

# Same for every model
A_full_years |>
  mutate(model = factor(model, levels = A_full_years |>
    filter(data == "Alternating", metric == levels(metric)[1]) |>
    summarise(A = mean(value), .by = model) |>
    arrange(-A) |>
    pull(model))) |>
  plot_full_years() +
  facet_grid(model ~ metric, scales = "free_x")
ggsave(here::here("figs/norway-A-full-years-models.pdf"), width = 6.5, height = 11)
