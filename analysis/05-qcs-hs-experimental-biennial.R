# Fit QCS-HS under various scenarios
# 1. all sampled - null case
# 2. all sampled at half effort
# 3. pure biennial
# 4. 25% partial overlap, same, biennial effort
# 5. year of overlap with full effort
# 6. two years of overlap with full effort

# Do this across a largish set of species
# Characterize the see-saw effect magnitude with IID fields and factor(years)

# Main question answered:
# With real data, can we induce it? Yes
# With real data, can simple sampling strategies ameliorate it? Not really? without going to a full coverage

library(dplyr)
library(furrr)
library(ggplot2)
library(sdmTMB)

theme_set(ggsidekick::theme_sleek())

Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1"
)

# Experiment controls --------------------------------------------------------

set.seed(42)

survey_names <- c("SYN QCS", "SYN HS")
seesaw_window <- 8L

# Same synoptic species as analysis/03-fit-bc-models.R.
species_to_fit <- c(
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

scenario_specs <- tibble::tribble(
  ~scenario, ~design, ~n_overlap_years,
  "All sampled", "all", 0L,
  "All sampled at half effort", "half", 0L,
  "Pure biennial", "biennial", 0L,
  "One year 25% partial overlap", "partial", 1L,
  "One year full-effort overlap", "full_overlap", 1L,
  "Two years full-effort overlap", "full_overlap", 2L
) |>
  mutate(scenario = factor(scenario, levels = scenario))

# Overlap is added at the end of the time series (2025, then 2023).

# Data and scenario construction --------------------------------------------

surveyjoin::cache_data()
surveyjoin::load_sql_data()

prepare_species_data <- function(species) {
  dat <- surveyjoin::get_data(species, regions = "pbs") |>
    mutate(year = lubridate::year(lubridate::ymd(date))) |>
    filter(survey_name %in% survey_names) |>
    select(
      survey_name, event_id, year, lon_start, lat_start, depth_m,
      effort, catch_weight, common_name
    ) |>
    filter(effort > 0) |>
    tidyr::drop_na(
      survey_name, event_id, year, lon_start, lat_start, depth_m,
      effort, catch_weight
    )

  # QCS was sampled alone in 2003 and 2004. Keep only occasions for which both
  # surveys are available so all experimental designs share the same years.
  complete_years <- dat |>
    distinct(year, survey_name) |>
    count(year, name = "n_surveys") |>
    filter(n_surveys == length(survey_names)) |>
    pull(year)

  dat |>
    filter(year %in% complete_years) |>
    mutate(
      occasion = match(year, sort(unique(year))),
      primary_survey = if_else(occasion %% 2L == 1L, "SYN HS", "SYN QCS")
    ) |>
    arrange(year, survey_name, event_id) |>
    sdmTMB::add_utm_columns(
      ll_names = c("lon_start", "lat_start"),
      utm_crs = 3156
    )
}

sample_rows <- function(dat, n, seed) {
  if (n >= nrow(dat)) {
    return(dat)
  }
  set.seed(seed)
  dat[sample.int(nrow(dat), size = n, replace = FALSE), , drop = FALSE]
}

sample_half_effort <- function(dat, seed) {
  dat |>
    group_by(year, survey_name) |>
    group_modify(\(.x, .y) {
      group_seed <- seed + .y$year[[1L]] + match(.y$survey_name[[1L]], survey_names)
      sample_rows(.x, n = floor(nrow(.x) / 2), seed = group_seed)
    }) |>
    ungroup()
}

build_scenario_data <- function(dat, design, n_overlap_years = 0L, seed = 42L) {
  stopifnot(
    design %in% c("all", "half", "biennial", "partial", "full_overlap"),
    n_overlap_years >= 0L
  )

  if (design == "all") {
    return(dat)
  }

  if (design == "half") {
    return(sample_half_effort(dat, seed = seed))
  }

  overlap_years <- if (n_overlap_years == 0L) {
    integer(0)
  } else {
    tail(sort(unique(dat$year)), n_overlap_years)
  }

  if (design == "full_overlap") {
    return(dat |>
      filter(survey_name == primary_survey | year %in% overlap_years))
  }

  biennial_dat <- dat |>
    filter(survey_name == primary_survey)

  if (design == "biennial") {
    return(biennial_dat)
  }

  # In the partial-overlap year, retain the same total number of tows as the
  # primary survey would have contributed: 75% primary and 25% other survey.
  stopifnot(design == "partial", length(overlap_years) == 1L)
  overlap_year <- overlap_years[[1L]]
  overlap_dat <- filter(dat, year == overlap_year)
  primary_name <- unique(overlap_dat$primary_survey)
  stopifnot(length(primary_name) == 1L)

  primary_dat <- filter(overlap_dat, survey_name == primary_name)
  secondary_dat <- filter(overlap_dat, survey_name != primary_name)
  n_total <- nrow(primary_dat)
  n_secondary <- min(round(0.25 * n_total), nrow(secondary_dat))
  n_primary <- n_total - n_secondary

  partial_dat <- bind_rows(
    sample_rows(primary_dat, n_primary, seed = seed + overlap_year),
    sample_rows(secondary_dat, n_secondary, seed = seed + overlap_year + 1L)
  )

  bind_rows(
    filter(biennial_dat, year != overlap_year),
    partial_dat
  ) |>
    arrange(year, survey_name, event_id)
}

prediction_grid <- surveyjoin::dfo_synoptic_grid |>
  filter(survey %in% survey_names) |>
  sdmTMB::add_utm_columns(c("lon", "lat"), utm_crs = 3156)

# Model fitting and index calculation ---------------------------------------

fit_one_scenario <- function(dat, grid, mesh, species, scenario, design,
    n_overlap_years, seed = 42L) {
  scenario_dat <- build_scenario_data(
    dat,
    design = design,
    n_overlap_years = n_overlap_years,
    seed = seed
  )

  sampling <- scenario_dat |>
    count(year, occasion, survey_name, primary_survey, name = "n_tows") |>
    mutate(species = species, scenario = as.character(scenario), .before = 1L)

  status <- tibble::tibble(
    species = species,
    scenario = as.character(scenario),
    converged = FALSE,
    stage = "fit",
    message = NA_character_
  )

  mesh0 <- make_mesh(scenario_dat, c("X", "Y"), mesh = mesh$mesh)
  fit <- tryCatch(
    sdmTMB(
      catch_weight ~ 0 + factor(year),
      data = scenario_dat,
      mesh = mesh0,
      offset = log(scenario_dat$effort),
      family = delta_gamma(type = "poisson-link"),
      time = "year",
      spatial = "on",
      spatiotemporal = "iid",
      share_range = TRUE,
      anisotropy = TRUE,
      silent = TRUE
    ),
    error = function(e) e
  )

  if (inherits(fit, "error")) {
    status$message <- conditionMessage(fit)
    return(list(index = tibble::tibble(), sampling = sampling, status = status))
  }

  sanity_ok <- tryCatch(
    isTRUE(all(unlist(sanity(fit, gradient_thresh = 0.01)))),
    error = function(e) FALSE
  )
  if (!sanity_ok) {
    status$stage <- "sanity"
    status$message <- "Model failed one or more sanity checks."
    return(list(index = tibble::tibble(), sampling = sampling, status = status))
  }

  years <- sort(unique(dat$year))
  newdata <- sdmTMB::replicate_df(grid, "year", years)
  index <- tryCatch(
    get_index(
      fit,
      newdata = newdata,
      offset = rep(0, nrow(newdata)),
      bias_correct = TRUE
    ),
    error = function(e) e
  )
  if (inherits(index, "error")) {
    status$stage <- "index"
    status$message <- conditionMessage(index)
    return(list(index = tibble::tibble(), sampling = sampling, status = status))
  }

  index <- index |>
    mutate(
      species = species,
      scenario = as.character(scenario),
      occasion = match(year, years),
      .before = 1L
    )
  status$converged <- TRUE
  status$stage <- "complete"

  list(index = index, sampling = sampling, status = status)
}

fit_species <- function(species, scenario_specs, grid, seed = 42L) {
  RhpcBLASctl::blas_set_num_threads(1L)
  RhpcBLASctl::omp_set_num_threads(1L)

  cli::cli_alert_info("Preparing {species}")
  dat <- prepare_species_data(species)
  # keep knots consistent across scenarios:
  mesh <- make_mesh(dat, c("X", "Y"), cutoff = 10)

  purrr::pmap(
    scenario_specs,
    \(scenario, design, n_overlap_years) {
      cli::cli_alert_info("Fitting {species}: {scenario}")
      fit_one_scenario(
        dat = dat,
        grid = grid,
        mesh = mesh,
        species = species,
        scenario = scenario,
        design = design,
        n_overlap_years = n_overlap_years,
        seed = seed
      )
    }
  )
}

# Seesaw metrics from analysis/09-seesaw-vignette.Rmd -----------------------

period2_metric <- function(x, year = seq_along(x), conf = 0.95) {
  stopifnot(
    length(x) == length(year),
    all(x > 0)
  )

  y <- log(x)

  # Alternating sampling phases, coded so beta is the phase difference.
  phase <- ifelse(year %% 2 == 0, 0.5, -0.5)
  fit <- lm(y ~ year + phase)

  beta <- coef(fit)[["phase"]]
  beta_se <- coef(summary(fit))["phase", "Std. Error"]
  ci <- confint(fit, "phase", level = conf)
  A_fun <- \(b) 100 * (exp(abs(b)) - 1)
  A_ci <- A_fun(ci)
  A_bounds <- if (ci[1L] <= 0 && ci[2L] >= 0) {
    c(0, max(A_ci))
  } else {
    range(A_ci)
  }

  c(
    beta = beta,
    beta_se = beta_se,
    beta_lower = ci[1L],
    beta_upper = ci[2L],
    A = A_fun(beta),
    A_lower = A_bounds[1L],
    A_upper = A_bounds[2L]
  )
}

moving_window <- function(x, year, window = 10L, ...) {
  n <- length(x)
  if (n < window) {
    return(data.frame())
  }

  out <- lapply(seq_len(n - window + 1L), function(i) {
    ind <- i:(i + window - 1L)
    res <- period2_metric(x = x[ind], year = year[ind], ...)

    data.frame(
      start = year[ind[1L]],
      end = year[ind[length(ind)]],
      t(res)
    )
  })

  do.call(rbind, out)
}

calculate_seesaw_metrics <- function(indexes, window = 8L) {
  full_series <- indexes |>
    arrange(species, scenario, occasion) |>
    group_by(species, scenario) |>
    group_modify(\(.x, .y) {
      as.data.frame(t(period2_metric(.x$est, year = .x$occasion)))
    }) |>
    ungroup()

  windows <- indexes |>
    arrange(species, scenario, occasion) |>
    group_by(species, scenario) |>
    group_modify(\(.x, .y) {
      out <- moving_window(.x$est, year = .x$occasion, window = window)
      if (nrow(out) > 0L) {
        year_lookup <- setNames(.x$year, .x$occasion)
        out$start_year <- unname(year_lookup[as.character(out$start)])
        out$end_year <- unname(year_lookup[as.character(out$end)])
      }
      out
    }) |>
    ungroup()

  list(full_series = full_series, moving_window = windows)
}

# Run and save ---------------------------------------------------------------

workers <- min(length(species_to_fit), future::availableCores()/2)
future::plan(future::multisession, workers = workers)

model_results <- furrr::future_map(
# model_results <- purrr::map(
  species_to_fit,
  fit_species,
  scenario_specs = scenario_specs,
  grid = prediction_grid,
  seed = 42L,
  .options = furrr::furrr_options(seed = TRUE)
)

future::plan(future::sequential)

flat_results <- unlist(model_results, recursive = FALSE)
indexes <- purrr::map_dfr(flat_results, "index") |>
  mutate(scenario = factor(scenario, levels = levels(scenario_specs$scenario)))
sampling_counts <- purrr::map_dfr(flat_results, "sampling") |>
  mutate(scenario = factor(scenario, levels = levels(scenario_specs$scenario)))
fit_status <- purrr::map_dfr(flat_results, "status") |>
  mutate(scenario = factor(scenario, levels = levels(scenario_specs$scenario)))

seesaw_metrics <- calculate_seesaw_metrics(indexes, window = seesaw_window)

dir.create(here::here("data-generated"), showWarnings = FALSE, recursive = TRUE)
saveRDS(
  list(
    indexes = indexes,
    full_series_seesaw = seesaw_metrics$full_series,
    moving_window_seesaw = seesaw_metrics$moving_window,
    sampling_counts = sampling_counts,
    fit_status = fit_status,
    scenarios = scenario_specs
  ),
  here::here("data-generated/qcs-hs-experimental-biennial.rds")
)

# Quick checks ---------------------------------------------------------------

out <- readRDS(here::here("data-generated/qcs-hs-experimental-biennial.rds"))
fit_status <- out$fit_status
seesaw_metrics <- list()
seesaw_metrics$full_series <- out$full_series_seesaw
seesaw_metrics$moving_window_seesaw <- out$moving_window_seesaw
indexes <- out$indexes

fit_status |>
  count(scenario, converged)

col <- RColorBrewer::brewer.pal(3, "Blues")[3]

col2 <- RColorBrewer::brewer.pal(3, "Blues")[1]

moving_window_seesaw <- seesaw_metrics$moving_window_seesaw |>
    mutate(scenario = gsub("Pure biennial", "Biennial", scenario)) |>
  mutate(scenario = reorder(scenario, A, FUN = mean))

moving_window_seesaw |>
  filter(scenario != "All sampled at half effort") |>
  ggplot(aes(scenario, A)) +
  geom_violin(quantile.colour = "black", quantile.linetype = 1, quantiles = c(0.5), colour = col, fill = col2) +
  # geom_boxplot() +
  geom_point(position = position_jitter(width = 0.1), size = 0.7, alpha = 0.3, colour = col) +
  geom_point(stat = "summary", fun = mean, size = 3, colour = "black") +
  scale_y_continuous(limits = c(0, NA), expand = c(0, 11)) +
  coord_flip(ylim = c(0, 300), expand = FALSE) +
  labs(
    x = NULL,
    y = "Estimated biennial amplitude (%)"
  )

indexes |>
  group_by(species, scenario) |>
  mutate(est_scaled = est / exp(mean(log(est)))) |>
  ungroup() |>
  ggplot(aes(year, est_scaled, group = scenario, colour = scenario)) +
  geom_line() +
  geom_point() +
  facet_wrap(~species, scales = "free_y") +
  labs(
    x = "Year",
    y = "Index / geometric mean",
    colour = "Sampling scenario"
  )
