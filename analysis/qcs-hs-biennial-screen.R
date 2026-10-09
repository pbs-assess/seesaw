# Screen species for strong seesaws under pretend biennial QCS / HS sampling.
#
# For each species, fit only IID RF, factor(year) to:
# - the full data (both regions every occasion; the 'truth'), and
# - the biennial subset (one region per occasion, alternating HS / QCS),
# and compare seesaw A and max A across moving windows. Used to pick species
# for analysis/07-qcs-hs-biennial-model-comparison.R and
# analysis/08-qcs-hs-biennial-self-test.R.

library(dplyr)
library(ggplot2)
library(sdmTMB)

theme_set(ggsidekick::theme_sleek())

source(here::here("analysis/fit-index-models.R"))
source(here::here("analysis/metric-functions.R"))

surveyjoin::cache_data()
surveyjoin::load_sql_data()

survey_names <- c("SYN QCS", "SYN HS")
seesaw_window <- 8L

# Candidate species from analysis/03-fit-syn-models.R
species_to_screen <- c(
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

# Model time as survey occasion (see analysis/07-qcs-hs-biennial-model-comparison.R)
prepare_species_data <- function(species) {
  dat <- surveyjoin::get_data(species, regions = "pbs") |>
    mutate(year = lubridate::year(lubridate::ymd(date))) |>
    filter(survey_name %in% survey_names, effort > 0) |>
    select(survey_name, year, lon_start, lat_start, depth_m, effort, catch_weight) |>
    tidyr::drop_na()

  # QCS was sampled alone in 2003 and 2004; keep years with both surveys
  complete_years <- dat |>
    distinct(year, survey_name) |>
    count(year) |>
    filter(n == length(survey_names)) |>
    pull(year)

  dat |>
    filter(year %in% complete_years) |>
    mutate(
      occasion = match(year, sort(unique(year))),
      primary_survey = if_else(occasion %% 2L == 1L, "SYN HS", "SYN QCS"),
      cal_year = year,
      year = occasion
    ) |>
    sdmTMB::add_utm_columns(ll_names = c("lon_start", "lat_start"), utm_crs = 3156)
}

screen_species <- function(species) {
  RhpcBLASctl::blas_set_num_threads(1L)
  RhpcBLASctl::omp_set_num_threads(1L)

  dat <- prepare_species_data(species)
  dat_biennial <- filter(dat, survey_name == primary_survey)

  grid <- surveyjoin::dfo_synoptic_grid |>
    filter(survey %in% survey_names) |>
    sdmTMB::add_utm_columns(c("lon", "lat"), utm_crs = 3156) |>
    clamp_depth(dat) |>
    sdmTMB::replicate_df("year", sort(unique(dat$year)))

  mesh <- make_mesh(dat, c("X", "Y"), cutoff = 10)

  fit_design <- function(d) {
    fit <- tryCatch(sdmTMB(
      catch_weight ~ 0 + factor(year),
      data = d,
      mesh = make_mesh(d, c("X", "Y"), mesh = mesh$mesh),
      offset = log(d$effort),
      family = delta_gamma(type = "poisson-link"),
      time = "year",
      spatial = "on",
      spatiotemporal = "iid",
      share_range = TRUE,
      anisotropy = TRUE,
      control = sdmTMBcontrol(collapse_spatial_variance = TRUE),
      silent = TRUE
    ), error = \(e) NULL)
    if (is.null(fit)) {
      return(tibble::tibble())
    }
    sanity_ok <- isTRUE(tryCatch(
      all(unlist(sanity(fit, gradient_thresh = 0.01, silent = TRUE))),
      error = \(e) FALSE
    ))
    # Bias correction skipped for speed; it barely affects the year-to-year
    # pattern that A measures
    index <- tryCatch(
      get_index(fit, newdata = grid, area = grid$area, offset = rep(0, nrow(grid)), bias_correct = FALSE),
      error = \(e) NULL
    )
    if (is.null(index)) {
      return(tibble::tibble())
    }
    mutate(index, sanity_ok = sanity_ok)
  }

  occasions <- distinct(dat, year, cal_year, primary_survey)

  bind_rows(
    full = fit_design(dat),
    biennial = fit_design(dat_biennial),
    .id = "design"
  ) |>
    left_join(occasions, by = "year") |>
    mutate(
      species = species,
      n_tows = nrow(dat),
      prop_pos = mean(dat$catch_weight > 0),
      .before = 1L
    )
}

screen_file <- here::here("data-generated/qcs-hs-biennial-screen.rds")
if (!file.exists(screen_file)) {
  future::plan(future::multisession, workers = min(length(species_to_screen), future::availableCores() / 2))
  indexes <- furrr::future_map_dfr(
    species_to_screen, screen_species,
    .options = furrr::furrr_options(seed = TRUE, scheduling = Inf)
  )
  future::plan(future::sequential)
  saveRDS(indexes, screen_file)
}
indexes <- readRDS(screen_file)

# Seesaw metrics ---------------------------------------------------------------

# Phase follows which region was sampled, not calendar-year parity
seesaw <- function(d) {
  d <- arrange(d, year)
  phase <- if_else(d$primary_survey == "SYN HS", 0.5, -0.5)
  mw <- moving_window(d$est, year = d$year, window = seesaw_window, phase = phase)
  as_tibble(t(period2_metric(d$est, year = d$year, phase = phase))) |>
    mutate(A_max = max(mw$A))
}

screen <- indexes |>
  group_by(species, design, n_tows, prop_pos) |>
  summarise(sanity_ok = all(sanity_ok), seesaw(pick(everything())), .groups = "drop") |>
  select(species, design, n_tows, prop_pos, sanity_ok, beta, A, A_lower, A_upper, A_max) |>
  tidyr::pivot_wider(
    names_from = design,
    values_from = c(sanity_ok, beta, A, A_lower, A_upper, A_max)
  ) |>
  mutate(
    A_added = A_biennial - A_full,
    A_max_added = A_max_biennial - A_max_full,
    # Positive beta: HS-sampled occasions higher than QCS-sampled ones
    high_region = if_else(beta_biennial > 0, "HS", "QCS")
  ) |>
  arrange(desc(A_max_biennial))

saveRDS(screen, here::here("data-generated/qcs-hs-biennial-screen-metrics.rds"))

screen |>
  select(
    species, prop_pos, sanity_ok_full, sanity_ok_biennial,
    A_full, A_biennial, A_lower_biennial, A_max_full, A_max_biennial, A_max_added, high_region
  ) |>
  mutate(across(where(is.numeric), \(x) round(x, 2))) |>
  print(n = Inf, width = Inf)

# Plots ----------------------------------------------------------------------

screen |>
  mutate(species = factor(species, levels = rev(screen$species))) |>
  ggplot(aes(y = species)) +
  geom_segment(aes(x = A_max_full, xend = A_max_biennial), colour = "grey60") +
  geom_point(aes(x = A_max_full, colour = "Full data")) +
  geom_point(aes(x = A_max_biennial, colour = "Biennial")) +
  labs(
    x = paste0("Max A across ", seesaw_window, "-survey windows (%)"),
    y = NULL, colour = "QCS + HS data",
    caption = "IID RF, factor(year) fit to each design."
  )
ggsave(here::here("figs/qcs-hs-biennial-screen.pdf"), width = 6, height = 5)

indexes |>
  group_by(species, design) |>
  mutate(est_scaled = est / exp(mean(log(est)))) |>
  ungroup() |>
  mutate(species = factor(species, levels = screen$species)) |>
  ggplot(aes(cal_year, est_scaled, colour = design)) +
  geom_line() +
  geom_point(aes(shape = primary_survey)) +
  facet_wrap(~species, scales = "free_y") +
  labs(x = "Year", y = "Index / geometric mean", colour = "Design", shape = "Sampled region (biennial)")
ggsave(here::here("figs/qcs-hs-biennial-screen-indexes.pdf"), width = 13, height = 9)
