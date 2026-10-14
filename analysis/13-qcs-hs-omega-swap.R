# Omega swap for the QCS + HS biennial self-test: does the real spatial
# structure (the fitted omega, including the HS / QCS contrast) cause the
# spurious seesaw?
#
# As analysis/08-qcs-hs-biennial-self-test.R (operating model without an HS
# term), but each replicate is simulated four times from the same parameter
# draws and seed, changing only omega:
#
#   conditional: omega from its approximate posterior (the real geography)
#   conditional_flat: conditional with each region's (HS, QCS) mean shifted to
#     the overall mean: the real within-region pattern without the HS / QCS
#     contrast
#   new_contrast: new with each region's mean shifted to that of the paired
#     conditional draw: the real HS / QCS contrast with a random within-region
#     pattern
#   new: omega replaced by an unconditional draw with the same range, sigma_O,
#     and anisotropy (no real geography)
#
# The new spatiotemporal fields and observation draws are identical across
# the four, so any difference in the excess seesaw is caused by omega. Each
# replicate's realized HS / QCS contrast is kept so the excess seesaw can be
# compared with contrasts that arise by chance in the new fields.

library(dplyr)
library(ggplot2)
library(sdmTMB)

theme_set(ggsidekick::theme_sleek())

source(here::here("analysis/fit-index-models.R"))
source(here::here("analysis/metric-functions.R"))

surveyjoin::cache_data()
surveyjoin::load_sql_data()

survey_names <- c("SYN QCS", "SYN HS")
species_to_fit <- c(
  "pacific cod",
  "spotted ratfish",
  "pacific ocean perch",
  "pacific halibut",
  "pacific spiny dogfish",
  "lingcod",
  "sablefish",
  "dover sole",
  # control: almost no seesaw in analysis/qcs-hs-biennial-screen.R
  "arrowtooth flounder"
  # shortspine thornyhead dropped: its OM encounter AR1 field sits at rho = 1,
  # so each replicate's new epsilon is a near-static field
)
n_reps <- 50L
# Quick check: fewer species and replicates, cached and plotted separately
test_run <- FALSE
if (test_run) {
  species_to_fit <- c("pacific ocean perch", "pacific halibut", "pacific spiny dogfish", "dover sole")
  n_reps <- 10L
}
out_suffix <- if (test_run) "-test" else ""
# Species for the example omega maps
map_species <- "pacific halibut"
n_workers <- max(1L, min(20L, future::availableCores() - 2L))
# Any subset of "conditional", "conditional_flat", "new_contrast", "new" (the
# random-pattern pair is the 2 x 2 check; see the header)
omega_types <- c("conditional", "conditional_flat")
ci_levels <- c(0.5, 0.95)
# As in 08: with 11 occasions, 8 gives 4 windows
seesaw_window <- 8L

om_model <- "AR1 RF, factor(year)"
# The only biennial model with a substantial excess seesaw in 08 / 12
est_model <- "IID RF, factor(year)"

family <- delta_gamma(type = "poisson-link")
control <- sdmTMBcontrol(collapse_spatial_variance = TRUE)

# Data -----------------------------------------------------------------------

# Model time as survey occasion (see 08)
prepare_species_data <- function(species) {
  dat <- surveyjoin::get_data(species, regions = "pbs") |>
    mutate(year = lubridate::year(lubridate::ymd(date))) |>
    filter(survey_name %in% survey_names, effort > 0) |>
    select(survey_name, year, lon_start, lat_start, depth_m, effort, catch_weight) |>
    tidyr::drop_na()

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

base_grid <- surveyjoin::dfo_synoptic_grid |>
  filter(survey %in% survey_names) |>
  mutate(hs = as.integer(survey == "SYN HS")) |>
  sdmTMB::add_utm_columns(c("lon", "lat"), utm_crs = 3156)

# Models ---------------------------------------------------------------------

# As in 08
fit_model <- function(model, d, mesh) {
  args <- list(
    formula = catch_weight ~ 0 + factor(year),
    data = d,
    mesh = mesh,
    offset = log(d$effort),
    family = family,
    time = "year",
    spatial = "on",
    spatiotemporal = switch(model,
      "IID RF, factor(year)" = "iid",
      "AR1 RF, factor(year)" = "ar1"
    ),
    share_range = TRUE,
    anisotropy = TRUE,
    control = control,
    silent = TRUE
  )
  do.call(sdmTMB, args)
}

# Operating model and paired simulation ----------------------------------------

# Shift omega (one delta component, at the mesh nodes) so each region's
# area-weighted mean on the grid moves to a target. Nodes take the region of
# their nearest grid cell; the projection smooths the shift across the HS /
# QCS boundary, so the realized means are close to, not exactly, the targets
# (the realized contrast is kept in the true indexes)
make_region_shifter <- function(mesh, grid) {
  xy <- as.matrix(grid[, c("X", "Y")])
  A <- fmesher::fm_basis(mesh$mesh, loc = xy)
  node_hs <- grid$hs[RANN::nn2(xy, mesh$mesh$loc[, 1:2], k = 1)$nn.idx[, 1]]
  region_area <- c(sum(grid$area[grid$hs == 0L]), sum(grid$area[grid$hs == 1L]))
  # rows: QCS, HS
  W <- rbind(
    (grid$hs == 0L) * grid$area / region_area[1],
    (grid$hs == 1L) * grid$area / region_area[2]
  )
  WA <- as.matrix(W %*% A)
  means <- \(omega) as.vector(WA %*% omega)
  list(
    project = \(omega) as.vector(A %*% omega),
    means = means,
    overall = \(omega) sum(means(omega) * region_area) / sum(region_area),
    shift = \(omega, target) omega + (target - means(omega))[node_hs + 1L]
  )
}

# dat is prepared in the main session (surveyjoin's database is not safe to
# load from parallel workers)
setup_species <- function(species, dat) {
  RhpcBLASctl::blas_set_num_threads(1L)
  RhpcBLASctl::omp_set_num_threads(1L)

  is_biennial <- dat$survey_name == dat$primary_survey
  grid <- sdmTMB::replicate_df(base_grid, "year", sort(unique(dat$year)))
  mesh <- make_mesh(dat, c("X", "Y"), cutoff = 10)
  mesh_biennial <- make_mesh(dat[is_biennial, ], c("X", "Y"), mesh = mesh$mesh)

  fit_om <- fit_ok(fit_model(om_model, dat, mesh))
  if (is.null(fit_om)) {
    cli::cli_warn("Operating model did not pass sanity checks for {species}")
    return(NULL)
  }

  # One set of parameter draws shared by all omega types: fixed effects at
  # their MLEs, random effects from their approximate posterior (as
  # mle_mvn_samples = "multiple")
  set.seed(42)
  par_draws <- sdmTMB:::.posterior_re_sampler(fit_om)(n_reps)
  # (the draws are unnamed; rows follow the joint parameter vector)
  idx_omega <- which(names(fit_om$tmb_obj$env$last.par.best) == "omega_s")
  stopifnot(length(idx_omega) > 0L)

  # Unconditional omega for each draw (its own range, sigma_O, anisotropy).
  # omega_s is nodes x delta components, column-major as in the par vector;
  # a component collapsed by collapse_spatial_variance is absent from the par
  # vector, so keep only the components that are present
  new_omega <- simulate(
    fit_om,
    nsim = n_reps,
    mcmc_samples = par_draws,
    simulate_re = "spatial",
    return_tmb_report = TRUE,
    seed = 99,
    silent = TRUE
  ) |>
    lapply(\(r) r$omega_s)
  n_nodes <- nrow(new_omega[[1]])
  keep <- seq_len(length(idx_omega) / n_nodes)
  stopifnot(length(idx_omega) %% n_nodes == 0L)
  par_new <- par_draws
  for (i in seq_len(n_reps)) par_new[idx_omega, i] <- as.vector(new_omega[[i]][, keep])
  stopifnot(!isTRUE(all.equal(par_new[idx_omega, 1L], par_draws[idx_omega, 1L])))

  # Regional shifts, per delta component
  shifter <- make_region_shifter(mesh, base_grid)
  stopifnot(nrow(mesh$mesh$loc) == n_nodes)
  par_flat <- par_draws
  par_new_contrast <- par_new
  for (i in seq_len(n_reps)) {
    for (j in split(idx_omega, rep(keep, each = n_nodes))) {
      omega_i <- par_draws[j, i]
      par_flat[j, i] <- shifter$shift(omega_i, rep(shifter$overall(omega_i), 2L))
      par_new_contrast[j, i] <- shifter$shift(par_new[j, i], shifter$means(omega_i))
    }
  }

  # Example maps: the first draw of each omega type on the grid, summed over
  # the delta components (the log density contribution)
  project_draw <- function(p) {
    Reduce(`+`, lapply(split(idx_omega, rep(keep, each = n_nodes)), \(j) shifter$project(p[j, 1L])))
  }
  omega_maps <- purrr::map_dfr(
    list(conditional = par_draws, conditional_flat = par_flat, new_contrast = par_new_contrast, new = par_new),
    \(p) mutate(select(base_grid, X, Y, hs, area), value = project_draw(p)),
    .id = "omega"
  )

  sim_newdata <- bind_rows(
    select(dat, X, Y, year, effort),
    mutate(select(grid, X, Y, year), effort = 1)
  )

  list(
    species = species,
    dat = dat,
    is_biennial = is_biennial,
    grid = grid,
    mesh_biennial = mesh_biennial,
    fit_om = fit_om,
    draws = list(
      conditional = par_draws,
      conditional_flat = par_flat,
      new_contrast = par_new_contrast,
      new = par_new
    )[omega_types],
    sim_newdata = sim_newdata,
    is_tow = seq_len(nrow(sim_newdata)) <= nrow(dat),
    om_pars = bind_rows(tidy(fit_om, "ran_pars", model = 1), tidy(fit_om, "ran_pars", model = 2)),
    occasions = distinct(dat, year, cal_year, primary_survey),
    omega_maps = mutate(omega_maps, species = species, .before = 1L)
  )
}

# All replicates of one omega type for one species
simulate_omega_type <- function(fit_om, draws, sim_newdata, is_tow, grid) {
  RhpcBLASctl::blas_set_num_threads(1L)
  RhpcBLASctl::omp_set_num_threads(1L)

  # One call (and seed) per replicate: the number of gamma draws depends on
  # the encounters, so a single stream would desynchronize the omega types
  # after the first replicate
  reports <- lapply(seq_len(ncol(draws)), \(i) {
    simulate(
      fit_om,
      nsim = 1L,
      mcmc_samples = draws[, i, drop = FALSE],
      simulate_re = "spatiotemporal",
      newdata = sim_newdata,
      offset = log(sim_newdata$effort),
      return_tmb_report = TRUE,
      seed = 42 + i,
      silent = TRUE
    )[[1L]]
  })
  truth <- purrr::imap_dfr(reports, \(r, rep) {
    eta <- r$eta_i[!is_tow, 1L] + r$eta_i[!is_tow, 2L]
    omega <- r$omega_s_A[!is_tow, 1L] + r$omega_s_A[!is_tow, 2L]
    grid |>
      mutate(density = exp(eta) * area, eta = eta, omega = omega) |>
      group_by(year) |>
      summarise(
        true_est = sum(density),
        # realized log density contrast (HS minus QCS, area-weighted)
        contrast = weighted.mean(eta[hs == 1L], area[hs == 1L]) -
          weighted.mean(eta[hs == 0L], area[hs == 0L]),
        # the part of it from omega
        contrast_omega = weighted.mean(omega[hs == 1L], area[hs == 1L]) -
          weighted.mean(omega[hs == 0L], area[hs == 0L])
      ) |>
      mutate(rep = rep, .before = 1L)
  })
  list(
    catches = lapply(reports, \(r) r$y_i[is_tow, 1L] * r$y_i[is_tow, 2L]),
    epsilon = lapply(reports, \(r) r$epsilon_st),
    truth = truth
  )
}

# Fit replicates --------------------------------------------------------------

fit_rep <- function(setup, omega, rep) {
  RhpcBLASctl::blas_set_num_threads(1L)
  RhpcBLASctl::omp_set_num_threads(1L)

  d <- setup$dat[setup$is_biennial, ]
  d$catch_weight <- setup$catches[[omega]][[rep]][setup$is_biennial]

  fit <- fit_ok(fit_model(est_model, d, setup$mesh_biennial))
  if (is.null(fit)) return(tibble::tibble())
  get_index_ok(fit, setup$grid) |>
    mutate(species = setup$species, omega = omega, rep = rep, .before = 1L)
}

if (FALSE) {
  tictoc::tic()
  s <- setup_species("spotted ratfish", prepare_species_data("spotted ratfish"))
  tictoc::toc()
  sims <- purrr::map(s$draws, simulate_omega_type,
    fit_om = s$fit_om, sim_newdata = s$sim_newdata, is_tow = s$is_tow, grid = s$grid
  )
  purrr::map_dfr(sims, "truth", .id = "omega") |>
    group_by(omega) |>
    summarise(across(c(contrast, contrast_omega), mean))
}

# One cache file per species; delete a file to refit that species (files
# missing an omega type or replicate are refit)
fits_dir <- here::here(paste0("data-generated/qcs-hs-omega-swap", out_suffix))
dir.create(fits_dir, showWarnings = FALSE)
fits_file <- \(species) file.path(fits_dir, paste0(gsub(" ", "-", species), ".rds"))
maps_file <- here::here(paste0("data-generated/qcs-hs-omega-swap", out_suffix, "-maps.rds"))
cache_ok <- function(f) {
  if (!file.exists(f)) return(FALSE)
  x <- readRDS(f)$true_indexes
  all(omega_types %in% x$omega) && max(x$rep) == n_reps
}

to_fit <- species_to_fit[!vapply(fits_file(species_to_fit), cache_ok, logical(1))]
if (length(to_fit) > 0L) {
  dats <- lapply(to_fit, prepare_species_data)
  future::plan(future::multisession, workers = n_workers)
  opts <- furrr::furrr_options(seed = TRUE, scheduling = Inf)

  # 1. Operating model and omega draws: one future per species
  setups <- furrr::future_map2(to_fit, dats, setup_species, .options = opts, .progress = TRUE) |>
    purrr::compact()

  # 2. Simulation: one future per species x omega type, each sent only what
  # it needs
  sim_tasks <- tidyr::expand_grid(i = seq_along(setups), omega = omega_types)
  sim_out <- furrr::future_pmap(
    list(
      fit_om = purrr::map(setups[sim_tasks$i], "fit_om"),
      draws = purrr::map2(setups[sim_tasks$i], sim_tasks$omega, \(s, o) s$draws[[o]]),
      sim_newdata = purrr::map(setups[sim_tasks$i], "sim_newdata"),
      is_tow = purrr::map(setups[sim_tasks$i], "is_tow"),
      grid = purrr::map(setups[sim_tasks$i], "grid")
    ),
    simulate_omega_type,
    .options = opts,
    .progress = TRUE
  )
  for (i in seq_along(setups)) {
    sims <- rlang::set_names(sim_out[sim_tasks$i == i], omega_types)
    # Pairing check: only omega differs between the omega types
    for (o in omega_types[-1]) stopifnot(identical(sims[[1]]$epsilon, sims[[o]]$epsilon))
    setups[[i]]$fit_om <- NULL
    setups[[i]]$draws <- NULL
    setups[[i]]$catches <- purrr::map(sims, "catches")
    setups[[i]]$true_indexes <- purrr::map_dfr(sims, "truth", .id = "omega")
  }

  # 3. Fits: one future per species x omega type x replicate
  tasks <- tidyr::expand_grid(i = seq_along(setups), omega = omega_types, rep = seq_len(n_reps))
  rep_indexes <- furrr::future_pmap(
    list(setup = setups[tasks$i], omega = tasks$omega, rep = tasks$rep), fit_rep,
    .options = opts,
    .progress = TRUE
  )
  future::plan(future::sequential)

  for (i in seq_along(setups)) {
    saveRDS(
      list(
        indexes = bind_rows(rep_indexes[tasks$i == i]),
        true_indexes = mutate(setups[[i]]$true_indexes, species = setups[[i]]$species, .before = 1L),
        occasions = mutate(setups[[i]]$occasions, species = setups[[i]]$species, .before = 1L),
        om_pars = mutate(setups[[i]]$om_pars, species = setups[[i]]$species, .before = 1L)
      ),
      fits_file(setups[[i]]$species)
    )
    if (setups[[i]]$species == map_species) saveRDS(setups[[i]]$omega_maps, maps_file)
  }
}

# Example maps (only the operating model fit is needed if they are missing)
if (!file.exists(maps_file)) {
  saveRDS(setup_species(map_species, prepare_species_data(map_species))$omega_maps, maps_file)
}
omega_maps <- readRDS(maps_file)

swap <- lapply(fits_file(species_to_fit), \(f) if (file.exists(f)) readRDS(f))
# (a cache file may hold more omega types than requested)
indexes <- filter(purrr::map_dfr(swap, "indexes"), omega %in% omega_types)
true_indexes <- filter(purrr::map_dfr(swap, "true_indexes"), omega %in% omega_types)
occasions <- purrr::map_dfr(swap, "occasions")

# Metrics ----------------------------------------------------------------------

compared <- indexes |>
  left_join(true_indexes, by = c("species", "omega", "rep", "year")) |>
  left_join(occasions, by = c("species", "year")) |>
  mutate(log_error = log_est - log(true_est))

for (level in ci_levels) {
  q <- qnorm(1 - (1 - level) / 2)
  compared[[paste0("covered_", level * 100)]] <- abs(compared$log_error) <= q * compared$se
}

phase_of <- \(primary_survey) if_else(primary_survey == "SYN HS", 0.5, -0.5)
seesaw_A_max <- function(est, primary_survey, year) {
  max(moving_window(est, year = year, window = seesaw_window, phase = phase_of(primary_survey))$A)
}
# Signed seesaw error: beta (HS-year minus QCS-year log difference, detrended)
# of the estimated index minus that of the true index, averaged over windows.
# A uses |beta| and a max over windows, so noise alone inflates it; noise
# averages to zero here
seesaw_beta <- function(est, primary_survey, year) {
  moving_window(est, year = year, window = seesaw_window, phase = phase_of(primary_survey))$beta
}

metrics_rep <- compared |>
  arrange(species, omega, rep, year) |>
  group_by(species, omega, rep) |>
  summarise(
    rmse = sqrt(mean(log_error^2)),
    coverage_50 = mean(covered_50),
    coverage_95 = mean(covered_95),
    A_max = seesaw_A_max(est, primary_survey, year),
    beta_error = mean(
      seesaw_beta(est, primary_survey, year) - seesaw_beta(true_est, primary_survey, year)
    ),
    .groups = "drop"
  )

A_truth <- true_indexes |>
  left_join(occasions, by = c("species", "year")) |>
  arrange(species, omega, rep, year) |>
  group_by(species, omega, rep) |>
  summarise(
    A_max_true = seesaw_A_max(true_est, primary_survey, year),
    contrast = mean(contrast),
    contrast_omega = mean(contrast_omega),
    .groups = "drop"
  )

metrics_rep <- metrics_rep |>
  left_join(A_truth, by = c("species", "omega", "rep")) |>
  mutate(
    A_max_excess = A_max - A_max_true,
    # as a percentage: how much higher HS years are than QCS years in the
    # estimated index, relative to the truth
    seesaw_error = 100 * (exp(beta_error) - 1)
  )

metrics <- metrics_rep |>
  group_by(species, omega) |>
  summarise(
    n_reps = n(),
    abs_contrast = mean(abs(contrast)),
    abs_contrast_omega = mean(abs(contrast_omega)),
    rmse = mean(rmse),
    coverage_50 = mean(coverage_50),
    A_max_excess = median(A_max_excess),
    seesaw_error = median(seesaw_error),
    .groups = "drop"
  )

# Paired effects on the excess seesaw (same replicate): of the real contrast
# given the real within-region pattern (fitted minus flat), of the real omega
# (fitted minus new), and of the real contrast given a random within-region
# pattern (new + contrast minus new); only those whose omega types were run
# (on the log scale for the signed error, so effects are symmetric)
comparisons <- list(
  fitted_minus_flat = c("conditional", "conditional_flat"),
  fitted_minus_new = c("conditional", "new"),
  new_contrast_minus_new = c("new_contrast", "new")
) |>
  purrr::keep(\(x) all(x %in% omega_types))
paired_wide <- metrics_rep |>
  select(species, rep, omega, A_max_excess, beta_error) |>
  tidyr::pivot_longer(c(A_max_excess, beta_error), names_to = "metric") |>
  tidyr::pivot_wider(names_from = omega, values_from = value)
paired <- purrr::imap_dfr(comparisons, \(x, name) {
  mutate(paired_wide, comparison = name, effect = .data[[x[1]]] - .data[[x[2]]])
}) |>
  mutate(comparison = factor(comparison, levels = names(comparisons))) |>
  group_by(metric, species, comparison) |>
  summarise(
    n_pairs = sum(!is.na(effect)),
    effect_median = median(effect, na.rm = TRUE),
    prop_positive = mean(effect > 0, na.rm = TRUE),
    .groups = "drop"
  )

saveRDS(
  list(compared = compared, metrics_rep = metrics_rep, metrics = metrics, paired = paired),
  here::here(paste0("data-generated/qcs-hs-omega-swap", out_suffix, "-metrics.rds"))
)

purrr::map_dfr(swap, "om_pars") |>
  select(species, model, term, estimate, std.error) |>
  print(n = Inf)
metrics |>
  mutate(across(where(is.numeric), \(x) round(x, 2))) |>
  print(n = Inf, width = Inf)
print(paired, n = Inf)

# Plots ----------------------------------------------------------------------

omega_labels <- c(
  conditional = "Fitted",
  conditional_flat = "Fitted, contrast removed",
  new_contrast = "New, fitted contrast added",
  new = "New"
)
species_labeller <- labeller(species = \(x) stringr::str_wrap(tools::toTitleCase(x), 12))
seesaw_label <- "Seesaw error (ratio of index error, HS vs. QCS survey years)"
# beta_error is a log ratio; label it as a ratio
ratio_breaks <- c(0.125, 0.25, 0.5, 1, 2, 4, 8)
ratio_labels <- paste0("×", ratio_breaks)

# Orient each species so its fitted-omega seesaw error is positive (the same
# flip applies to every omega type)
orientation <- metrics_rep |>
  filter(omega == "conditional") |>
  group_by(species) |>
  summarise(
    fitted_contrast = mean(contrast_omega),
    orientation = sign(median(beta_error, na.rm = TRUE))
  )
metrics_oriented <- metrics_rep |>
  left_join(orientation, by = "species") |>
  mutate(beta_error_oriented = orientation * beta_error)

# All four omega types (supplement)
metrics_oriented |>
  mutate(omega = factor(omega, levels = names(omega_labels), labels = omega_labels)) |>
  ggplot(aes(omega, beta_error_oriented)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey60") +
  geom_line(aes(group = rep), colour = "grey80") +
  geom_point(colour = "grey60", alpha = 0.5) +
  stat_summary(aes(group = 1), fun = median, geom = "line", colour = "black") +
  stat_summary(fun = median, geom = "point", colour = "black") +
  scale_y_continuous(breaks = log(ratio_breaks), labels = ratio_labels) +
  facet_wrap(~species, nrow = 1, labeller = species_labeller, scales = "free_y") +
  labs(
    x = NULL,
    y = sub(" (", "\n(", seesaw_label, fixed = TRUE),
    caption = paste0(
      "QCS/HS biennial self-test (", est_model, ") with omega modified or replaced. ",
      "Grey: paired replicates (same epsilon and seed); black: median. ",
      "Oriented so each species' fitted-omega median is above 1."
    )
  ) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(here::here(paste0("figs/qcs-hs-omega-swap", out_suffix, ".pdf")), width = 13, height = 4.5, device = cairo_pdf)

# Seesaw error against the realized contrast (not oriented): its sign follows
# the contrast's
metrics_rep |>
  mutate(omega = factor(omega, levels = names(omega_labels), labels = omega_labels)) |>
  ggplot(aes(contrast, beta_error, colour = omega)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey60") +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey60") +
  geom_point() +
  scale_colour_manual(values = setNames(c("black", "#E69F00", "#56B4E9", "grey60"), omega_labels)) +
  scale_y_continuous(breaks = log(ratio_breaks), labels = ratio_labels) +
  facet_wrap(~species, nrow = 1, labeller = species_labeller, scales = "free_y") +
  labs(
    x = "Realized HS - QCS log-density contrast (replicate mean)",
    y = "Seesaw error (ratio of index error, HS / QCS years)",
    colour = NULL
  ) +
  theme(legend.position = "bottom")
ggsave(here::here(paste0("figs/qcs-hs-omega-swap", out_suffix, "-contrast.pdf")), width = 13, height = 4.5, device = cairo_pdf)

# Main figure: (a) example omega with and without the HS / QCS contrast, (b)
# seesaw error for the same two

main_arms <- c(conditional = "Fitted omega", conditional_flat = "HS-QCS contrast removed")
main_colours <- c("#1F6FB4", "grey55")

maps <- omega_maps |>
  filter(omega %in% names(main_arms)) |>
  select(omega, X, Y, hs, area, value) |>
  # centred on the coast-wide mean (the same in both by construction)
  mutate(value = value - weighted.mean(value, area), .by = omega)
# Fitted = contrast removed + difference: the difference is the removed
# regional contrast (two near-flat blocks, smoothed across the boundary)
map_labels <- c(conditional = "Fitted omega", conditional_flat = "Contrast removed", difference = "Difference")
maps <- maps |>
  tidyr::pivot_wider(names_from = omega, values_from = value) |>
  mutate(difference = conditional - conditional_flat) |>
  tidyr::pivot_longer(names(map_labels), names_to = "omega") |>
  mutate(omega = factor(omega, levels = names(map_labels), labels = map_labels))

# Each region's mean, labelled at the region's centre
region_means <- maps |>
  group_by(omega, hs) |>
  summarise(X = mean(X), Y = mean(Y), value = weighted.mean(value, area), .groups = "drop") |>
  mutate(
    value = round(value, 1) + 0, # no "-0.0"
    label = paste0(if_else(hs == 1L, "HS", "QCS"), ": ", if_else(value == 0, "0.0", sprintf("%+.1f", value)))
  )

fill_limit <- quantile(abs(maps$value), 0.99)
p_maps <- ggplot(maps, aes(X, Y)) +
  geom_tile(aes(fill = value), width = 2, height = 2) +
  geom_label(aes(label = label), data = region_means, size = 3) +
  scale_fill_gradient2(
    low = "#2166AC", mid = "grey95", high = "#B2182B",
    limits = c(-1, 1) * fill_limit, oob = scales::squish
  ) +
  facet_wrap(~omega) +
  coord_fixed() +
  labs(
    title = paste0(tools::toTitleCase(map_species), ", one posterior draw"),
    fill = "Spatial\neffect\n(log)", x = NULL, y = NULL
  ) +
  theme(axis.text = element_blank(), axis.ticks = element_blank())

dots <- metrics_oriented |>
  filter(omega %in% names(main_arms), !is.na(beta_error)) |>
  group_by(species, omega, fitted_contrast) |>
  summarise(
    lwr = quantile(beta_error_oriented, 0.25),
    med = median(beta_error_oriented),
    upr = quantile(beta_error_oriented, 0.75),
    .groups = "drop"
  ) |>
  mutate(
    species_label = paste0(tools::toTitleCase(species), " (", sprintf("%+.1f", fitted_contrast), ")"),
    species_label = forcats::fct_reorder(species_label, abs(fitted_contrast)),
    y = as.integer(species_label) + if_else(omega == "conditional", 0.13, -0.13),
    omega = factor(omega, levels = names(main_arms), labels = main_arms)
  )

p_dots <- ggplot(dots, aes(med, y, colour = omega)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey60") +
  geom_linerange(aes(xmin = lwr, xmax = upr)) +
  geom_point() +
  expand_limits(x = log(0.5)) +
  scale_colour_manual(values = main_colours) +
  scale_x_continuous(breaks = log(ratio_breaks), labels = ratio_labels) +
  scale_y_continuous(
    breaks = seq_along(levels(dots$species_label)),
    labels = levels(dots$species_label),
    expand = expansion(add = 0.5)
  ) +
  labs(x = seesaw_label, y = NULL, colour = NULL) +
  theme(legend.position = "bottom")

patchwork::wrap_plots(p_maps, p_dots, ncol = 1, heights = c(2.4, 0.6 + 0.4 * n_distinct(dots$species))) +
  patchwork::plot_annotation(tag_levels = "a")
ggsave(
  here::here(paste0("figs/qcs-hs-omega-swap", out_suffix, "-main.pdf")),
  width = 7, height = 4.2 + 0.4 * n_distinct(dots$species), device = cairo_pdf
)
