# Shared code for the simulation self-tests:
#   analysis/08-qcs-hs-biennial-self-test.R
#   analysis/11-nwfsc-split-self-test.R
#   analysis/12-qcs-hs-knockout.R
#   analysis/13-qcs-hs-omega-swap.R
#
# Each script supplies its data, operating model setup, and how one replicate
# is simulated; run_self_test() does the rest. Every random draw for a
# replicate is seeded by its replicate number, and every stage is cached per
# species or per replicate, so raising n_reps later only simulates and fits
# the new replicates (earlier ones are reused unchanged).

library(dplyr)
library(ggplot2)
library(sdmTMB)

source(here::here("analysis/fit-index-models.R"))
source(here::here("analysis/metric-functions.R"))

# Models ---------------------------------------------------------------------

st_family <- delta_gamma(type = "poisson-link")
# Drop random fields with SDs collapsing to zero (e.g., the encounter
# spatiotemporal field with alternating data)
st_control <- sdmTMBcontrol(collapse_spatial_variance = TRUE)
st_rw_prior <- sdmTMBpriors(sigma_V = sigma_V_prior(0.3, 0.5))

# Same specifications as the corresponding models in fit_index_models(), fit
# directly to skip the factor(year) template fit. Every occasion has data, so
# extra_time only matters for documenting the time steps. `formula` replaces
# the model's formula (e.g., to add a covariate to the operating model).
fit_model <- function(model, d, mesh, formula = NULL) {
  args <- list(
    formula = catch_weight ~ 1,
    data = d,
    mesh = mesh,
    offset = log(d$effort),
    family = st_family,
    time = "year",
    spatial = "on",
    share_range = TRUE,
    anisotropy = TRUE,
    extra_time = seq(min(d$year), max(d$year)),
    control = st_control,
    silent = TRUE
  )
  factor_year <- list(formula = catch_weight ~ 0 + factor(year), extra_time = NULL)
  rw_year <- list(time_varying = ~1, time_varying_type = "rw0", priors = st_rw_prior)
  extra <- switch(model,
    "IID RF, factor(year)" = c(spatiotemporal = "iid", factor_year),
    "AR1 RF, factor(year)" = c(spatiotemporal = "ar1", factor_year),
    "IID RF, RW year" = c(spatiotemporal = "iid", rw_year),
    "RW RF" = list(spatiotemporal = "rw"),
    "RW RF, RW year" = c(spatiotemporal = "rw", rw_year),
    "AR1 RF, RW year" = c(spatiotemporal = "ar1", rw_year),
    cli::cli_abort("Unknown model: {model}")
  )
  args[names(extra)] <- extra
  if (!is.null(formula)) args$formula <- formula
  do.call(sdmTMB, args)
}

# Fit the operating model. Returns list(fit, ran_pars) or, if it errors or
# fails sanity checks, list(om_failed = TRUE, ...) with diagnostics so the
# cause can be seen without refitting.
fit_operating_model <- function(model, d, mesh, formula = NULL) {
  fit <- tryCatch(fit_model(model, d, mesh, formula), error = \(e) e)
  if (!inherits(fit, "sdmTMB")) {
    return(list(om_failed = TRUE, error = conditionMessage(fit)))
  }
  checks <- tryCatch(
    unlist(sanity(fit, gradient_thresh = 0.01, silent = TRUE)),
    error = \(e) c(sanity_error = FALSE)
  )
  tidy_ok <- \(m) tryCatch(mutate(tidy(fit, "ran_pars", model = m), model = m), error = \(e) NULL)
  ran_pars <- bind_rows(tidy_ok(1), tidy_ok(2))
  if (all(checks)) {
    return(list(fit = fit, ran_pars = ran_pars))
  }
  list(
    om_failed = TRUE,
    sanity = checks,
    pd_hessian = fit$pos_def_hessian,
    convergence = fit$model$message,
    gradients = tibble::tibble(
      par = names(fit$model$par),
      est = unname(fit$model$par),
      gradient = as.numeric(fit$gradients)
    ),
    ran_pars = ran_pars
  )
}

# Fit one model to one replicate's simulated catches: to all tows ("full") or
# to the sampled half of each occasion ("alternating"). Returns the index, or
# an empty tibble if the fit fails sanity checks.
fit_index <- function(setup, catch, design, model, index_nsplit = 2L) {
  d <- setup$dat
  d$catch_weight <- catch
  mesh <- setup$mesh
  if (design == "alternating") {
    d <- d[setup$is_sampled, ]
    mesh <- setup$mesh_alt
  }
  fit <- fit_ok(fit_model(model, d, mesh))
  index <- if (!is.null(fit)) get_index_ok(fit, setup$grid, nsplit = index_nsplit)
  if (is.null(index)) tibble::tibble() else index
}

# QCS + HS data (08, 12, 13) -------------------------------------------------

qcs_hs_surveys <- c("SYN QCS", "SYN HS")

# Model time as survey occasion: QCS + HS years are all odd, so on the
# calendar-year scale extra_time puts an empty year between every survey and
# AR1 rho gets stuck at its starting value of 0 (only rho^2 is identified).
# Occasions alternate between HS (odd) and QCS (even) as the sampled region.
qcs_hs_data <- function(species) {
  dat <- surveyjoin::get_data(species, regions = "pbs") |>
    mutate(year = lubridate::year(lubridate::ymd(date))) |>
    filter(survey_name %in% qcs_hs_surveys, effort > 0) |>
    select(survey_name, year, lon_start, lat_start, depth_m, effort, catch_weight) |>
    tidyr::drop_na()

  # QCS was sampled alone in 2003 and 2004; keep years with both surveys
  complete_years <- dat |>
    distinct(year, survey_name) |>
    count(year) |>
    filter(n == length(qcs_hs_surveys)) |>
    pull(year)

  dat |>
    filter(year %in% complete_years) |>
    mutate(
      hs = as.integer(survey_name == "SYN HS"),
      occasion = match(year, sort(unique(year))),
      primary_survey = if_else(occasion %% 2L == 1L, "SYN HS", "SYN QCS"),
      # seesaw phase follows which region was sampled, not calendar-year parity
      phase = if_else(primary_survey == "SYN HS", 0.5, -0.5),
      cal_year = year,
      year = occasion
    ) |>
    add_utm_columns(ll_names = c("lon_start", "lat_start"), utm_crs = 3156)
}

qcs_hs_grid <- function() {
  surveyjoin::dfo_synoptic_grid |>
    filter(survey %in% qcs_hs_surveys) |>
    mutate(hs = as.integer(survey == "SYN HS"), phase_region = hs) |>
    add_utm_columns(c("lon", "lat"), utm_crs = 3156)
}

# The parts of a setup common to every QCS + HS script: the sampled half, the
# prediction grid, meshes, and the newdata for simulating tows and grid cells
qcs_hs_setup <- function(species, dat, grid) {
  is_sampled <- dat$survey_name == dat$primary_survey
  mesh <- make_mesh(dat, c("X", "Y"), cutoff = 10)
  sim_setup(
    species = species,
    dat = dat,
    is_sampled = is_sampled,
    grid = replicate_df(grid, "year", sort(unique(dat$year))),
    mesh = mesh,
    mesh_alt = make_mesh(dat[is_sampled, ], c("X", "Y"), mesh = mesh$mesh),
    occasions = distinct(dat, year, cal_year, primary_survey, phase),
    # for the HS term in the operating model of 12
    covariates = "hs"
  )
}

# Simulation -----------------------------------------------------------------

# Bundle what the simulation and fits need. Tows and grid cells go in one
# newdata so both share the same field draws; eta_i excludes the offset, so
# grid rows give expected density per unit area, and the simulated tow catches
# include the effort offset.
sim_setup <- function(species, dat, is_sampled, grid, mesh, mesh_alt, occasions, covariates = NULL) {
  cols <- c("X", "Y", "year", covariates)
  sim_newdata <- bind_rows(
    dat[, c(cols, "effort")],
    mutate(grid[, cols], effort = 1)
  )
  list(
    species = species,
    dat = dat,
    is_sampled = is_sampled,
    grid = grid,
    mesh = mesh,
    mesh_alt = mesh_alt,
    occasions = occasions,
    sim_newdata = sim_newdata,
    is_tow = seq_len(nrow(sim_newdata)) <= nrow(dat)
  )
}

# Seeds: the parameter draws, the simulation, and any extra draws (e.g., the
# new omega in 13) use separate seed ranges so their random streams don't
# overlap; each is offset by the replicate number
seed_pars <- 1e6L
seed_sim <- 2e6L
seed_extra <- 3e6L

# One replicate's parameters: fixed effects at their MLEs, random effects drawn
# from their approximate posterior (as mle_mvn_samples = "multiple"). `sampler`
# is sdmTMB:::.posterior_re_sampler(fit), built once per fit. A one-column
# matrix for simulate(mcmc_samples = ).
draw_pars <- function(sampler, rep) {
  set.seed(seed_pars + rep)
  sampler(1L)
}

# Simulate one replicate at the tows and grid cells of setup$sim_newdata from
# parameter draws `pars`; the random effects in simulate_re are replaced by
# entirely new draws from the GMRF prior. The same rep and simulate_re give
# the same new fields, so arms that differ only in `pars` are paired.
simulate_pars <- function(setup, pars, rep, simulate_re = "spatiotemporal") {
  simulate(
    setup$fit_om,
    nsim = 1L,
    mcmc_samples = pars,
    simulate_re = simulate_re,
    newdata = setup$sim_newdata,
    offset = log(setup$sim_newdata$effort),
    return_tmb_report = TRUE,
    seed = seed_sim + rep,
    silent = TRUE
  )[[1L]]
}

# Simulated catches at the tows
rep_catch <- function(setup, report) {
  report$y_i[setup$is_tow, 1L] * report$y_i[setup$is_tow, 2L]
}

# True index by year: the area-weighted sum of the expected density on the
# grid. With a `phase_region` grid column (1 for the region sampled in phase
# +0.5 occasions: HS, North), also the realized contrast of that region minus
# the other (area-weighted) in mean log density, the part of it from omega,
# and in log mean density (which also reflects each region's within-region
# variance).
rep_truth <- function(setup, report) {
  grid_rows <- !setup$is_tow
  out <- setup$grid |>
    mutate(
      eta = report$eta_i[grid_rows, 1L] + report$eta_i[grid_rows, 2L],
      omega = report$omega_s_A[grid_rows, 1L] + report$omega_s_A[grid_rows, 2L],
      density = exp(eta) * area
    ) |>
    group_by(year)
  if (!"phase_region" %in% names(setup$grid)) {
    return(summarise(out, true_est = sum(density)))
  }
  a_minus_b <- \(x, r, area) weighted.mean(x[r == 1L], area[r == 1L]) - weighted.mean(x[r == 0L], area[r == 0L])
  summarise(
    out,
    true_est = sum(density),
    contrast = a_minus_b(eta, phase_region, area),
    contrast_omega = a_minus_b(omega, phase_region, area),
    contrast_density = log(sum(density[phase_region == 1L]) / sum(area[phase_region == 1L])) -
      log(sum(density[phase_region == 0L]) / sum(area[phase_region == 0L]))
  )
}

# simulate_rep for a self-test with a single arm: one simulation per replicate
simulate_one_arm <- function(simulate_re = "spatiotemporal") {
  # evaluate now: the returned function runs in workers
  force(simulate_re)
  function(setup, sampler, rep) {
    report <- simulate_pars(setup, draw_pars(sampler, rep), rep, simulate_re)
    list(
      catches = list(base = rep_catch(setup, report)),
      truth = mutate(rep_truth(setup, report), arm = "base", .before = 1L)
    )
  }
}

# Pipeline -------------------------------------------------------------------

slug <- \(x) gsub("^-|-$", "", gsub("[^a-z0-9]+", "-", tolower(x)))
model_label <- \(model, design) if_else(design == "full", paste0(model, " (full data)"), model)
minutes_since <-\(t0) round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)
single_thread <- function() {
  RhpcBLASctl::blas_set_num_threads(1L)
  RhpcBLASctl::omp_set_num_threads(1L)
}

log_progress <- function(log_file, ...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S "), ..., "\n", sep = "", file = log_file, append = TRUE)
}

# Evaluate expr (which returns a note for the log, or NULL), logging when it's
# done or, on error, the error instead of stopping so other tasks continue
run_logged <- function(log_file, label, expr) {
  t0 <- Sys.time()
  tryCatch(
    log_progress(log_file, "done ", label, expr, " [", minutes_since(t0), " min]"),
    error = \(e) log_progress(log_file, "ERROR ", label, ": ", conditionMessage(e))
  )
  invisible(NULL)
}

# Workers for the three stages of run_self_test(). They're defined here, not
# inside run_self_test(), and take everything as arguments: a closure sent to
# a worker carries its enclosing environment, which would ship every species'
# setup to every task.
setup_worker <- function(species, dat, file, setup_species, log_file) {
  single_thread()
  run_logged(log_file, paste0("setup: ", species), {
    setup <- setup_species(species, dat)
    saveRDS(setup, file)
    if (isTRUE(setup$om_failed)) {
      why <- if (!is.null(setup$error)) setup$error else paste(names(setup$sanity)[!setup$sanity], collapse = ", ")
      paste0(" (OM FAILED: ", why, ")")
    }
  })
}

sim_worker <- function(setup, rep, file, simulate_rep, log_file) {
  single_thread()
  run_logged(log_file, paste0("sim: ", setup$species, " rep ", rep), {
    # A fit read from disk needs its TMB object rebuilt
    setup$fit_om <- sdmTMB:::reload_model(setup$fit_om)
    sampler <- sdmTMB:::.posterior_re_sampler(setup$fit_om)
    saveRDS(simulate_rep(setup, sampler, rep), file)
    NULL
  })
}

fit_worker <- function(setup, catch, rep, arm, design, model, file, index_nsplit, log_file) {
  single_thread()
  run_logged(log_file, paste0("fit: ", setup$species, " rep ", rep, " ", arm, " ", model, " (", design, ")"), {
    index <- fit_index(setup, catch, design, model, index_nsplit)
    saveRDS(index, file)
    if (nrow(index) == 0L) " FAILED SANITY"
  })
}

# Run a self-test in three cached, parallel stages; only what's missing is
# run:
#
#   1. Setup, one per species: setup_species(species, dat) fits the operating
#      model and returns a sim_setup() list plus fit_om (or om_failed = TRUE
#      with diagnostics). Cached in <dir>/setups/.
#   2. Simulation, one per species x replicate: simulate_rep(setup, sampler,
#      rep) returns list(catches = named list of tow catches, one per arm,
#      truth = true index with an `arm` column). Arms are versions of the same
#      replicate (e.g., the k of 12 or omega types of 13). Cached in
#      <dir>/sims/.
#   3. Fits, one per species x replicate x arm x model: each row of `models`
#      (design = "full" or "alternating", model) is fit to each arm in `arms`
#      (default all). Cached in <dir>/fits/.
#
# Errors are logged and not cached, so a rerun retries them. A failed
# operating model or a fit that fails sanity checks is cached so it isn't
# rerun; delete a file to redo it. dat is prepared in the main session
# (surveyjoin's database is not safe to load from parallel workers).
# Monitor with: tail -f <dir>/progress.log
#
# Returns indexes (successful fits only) and true indexes (species, rep, arm,
# year) for replicates 1:n_reps; fit_status, the outcome of every requested
# fit, and failure_rates by species x arm x model; occasions; the setups
# without their fits and data (for anything else a script stored there); and
# the species whose operating model failed.
run_self_test <- function(dir, species, n_reps, prepare_data, setup_species, simulate_rep,
                          models, arms = NULL, n_workers = 1L, index_nsplit = 2L) {
  for (d in file.path(dir, c("setups", "sims", "fits"))) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  setup_file <- \(s) file.path(dir, "setups", paste0(slug(s), ".rds"))
  sim_file <- \(s, rep) file.path(dir, "sims", sprintf("%s-rep%03d.rds", slug(s), rep))
  fit_file <- \(s, rep, arm, design, model) {
    file.path(dir, "fits", sprintf("%s-rep%03d-%s-%s-%s.rds", slug(s), rep, slug(arm), design, slug(model)))
  }
  log_file <- file.path(dir, "progress.log")

  future::plan(future::multisession, workers = n_workers)
  on.exit(future::plan(future::sequential), add = TRUE)
  opts <- furrr::furrr_options(seed = TRUE, scheduling = Inf)
  log_progress(log_file, "=== run started: ", length(species), " species, ", n_reps, " reps, ", n_workers, " workers ===")

  # 1. Setup
  todo <- species[!file.exists(setup_file(species))]
  if (length(todo) > 0L) {
    log_progress(log_file, length(todo), " setups to run")
    furrr::future_pwalk(
      list(species = todo, dat = lapply(todo, prepare_data), file = setup_file(todo)),
      setup_worker,
      setup_species = setup_species, log_file = log_file,
      .options = opts
    )
  }
  setups <- lapply(setNames(species, species), \(s) if (file.exists(setup_file(s))) readRDS(setup_file(s)))
  failed <- names(purrr::keep(setups, \(x) is.null(x) || isTRUE(x$om_failed)))
  if (length(failed) > 0L) cli::cli_warn("No operating model for: {failed}; see {log_file}")
  setups <- setups[setdiff(species, failed)]

  # 2. Simulation
  sims <- tidyr::expand_grid(species = names(setups), rep = seq_len(n_reps)) |>
    mutate(file = sim_file(species, rep))
  todo <- filter(sims, !file.exists(file))
  if (nrow(todo) > 0L) {
    log_progress(log_file, nrow(todo), " simulations to run")
    furrr::future_pwalk(
      list(setup = setups[todo$species], rep = todo$rep, file = todo$file),
      sim_worker,
      simulate_rep = simulate_rep, log_file = log_file,
      .options = opts
    )
  }
  sims <- filter(sims, file.exists(file))
  sim_out <- lapply(sims$file, readRDS)

  # 3. Fits; each task gets only its species' data, grid, and meshes plus its
  # replicate's catches for one arm
  fit_tasks <- purrr::imap(sim_out, \(sim, i) {
    sim_arms <- names(sim$catches)
    tidyr::expand_grid(i_sim = i, arm = if (is.null(arms)) sim_arms else intersect(arms, sim_arms), models)
  }) |>
    bind_rows(tibble::tibble(i_sim = integer(), arm = character(), design = character(), model = character())) |>
    mutate(
      species = sims$species[i_sim],
      rep = sims$rep[i_sim],
      file = fit_file(species, rep, arm, design, model)
    )
  todo <- filter(fit_tasks, !file.exists(file))
  if (nrow(todo) > 0L) {
    log_progress(log_file, nrow(todo), " of ", nrow(fit_tasks), " fits to run")
    setups_lite <- lapply(setups, \(x) x[c("species", "dat", "is_sampled", "grid", "mesh", "mesh_alt")])
    furrr::future_pwalk(
      list(
        setup = setups_lite[todo$species],
        catch = purrr::map2(todo$i_sim, todo$arm, \(i, arm) sim_out[[i]]$catches[[arm]]),
        rep = todo$rep, arm = todo$arm, design = todo$design, model = todo$model, file = todo$file
      ),
      fit_worker,
      index_nsplit = index_nsplit, log_file = log_file,
      .options = opts
    )
  }
  log_progress(log_file, "=== run finished ===")

  # The outcome of every requested fit, so failures are reported rather than
  # silently dropped: "ok", "failed sanity" (cached empty index), "fit error"
  # (not cached), "sim error" (no simulation for the replicate or arm), or "om
  # failed"
  done <- filter(fit_tasks, file.exists(file))
  done_indexes <- lapply(done$file, readRDS)
  done$status <- if_else(vapply(done_indexes, nrow, integer(1L)) > 0L, "ok", "failed sanity")
  requested_arms <- if (is.null(arms)) unique(unlist(lapply(sim_out, \(sim) names(sim$catches)))) else arms
  if (length(requested_arms) == 0L) requested_arms <- NA_character_
  fit_status <- tidyr::expand_grid(species = species, rep = seq_len(n_reps), arm = requested_arms, models) |>
    left_join(select(fit_tasks, species, rep, arm, design, model, file), by = c("species", "rep", "arm", "design", "model")) |>
    left_join(select(done, file, done_status = status), by = "file") |>
    mutate(
      status = case_when(
        species %in% failed ~ "om failed",
        is.na(file) ~ "sim error",
        !is.na(done_status) ~ done_status,
        .default = "fit error"
      ),
      model = model_label(model, design)
    ) |>
    select(species, rep, arm, design, model, status)
  failure_rates <- fit_status |>
    group_by(species, arm, model) |>
    summarise(
      n_requested = n(),
      n_ok = sum(status == "ok"),
      failure_rate = 1 - n_ok / n_requested,
      n_failed_sanity = sum(status == "failed sanity"),
      n_fit_error = sum(status == "fit error"),
      n_sim_error = sum(status == "sim error"),
      n_om_failed = sum(status == "om failed"),
      .groups = "drop"
    )
  if (any(failure_rates$failure_rate > 0)) {
    cli::cli_warn("{sum(fit_status$status != 'ok')} of {nrow(fit_status)} requested fit(s) failed; see failure_rates and {log_file}")
    print(filter(failure_rates, failure_rate > 0), n = Inf, width = Inf)
  }

  ok <- done$status == "ok"
  if (!any(ok)) cli::cli_abort("No fits to summarise; see {log_file}")
  list(
    indexes = purrr::pmap_dfr(
      list(done_indexes[ok], done$species[ok], done$rep[ok], done$arm[ok], done$design[ok], done$model[ok]),
      \(index, s, rep, arm, design, model) {
        mutate(index, species = s, rep = rep, arm = arm, design = design, model = model_label(model, design), .before = 1L)
      }
    ),
    fit_status = fit_status,
    failure_rates = failure_rates,
    true_indexes = purrr::imap_dfr(sim_out, \(sim, i) {
      mutate(sim$truth, species = sims$species[i], rep = sims$rep[i], .before = 1L)
    }),
    occasions = purrr::map_dfr(setups, "occasions", .id = "species"),
    setups = lapply(setups, \(x) x[setdiff(names(x), c("fit_om", "dat", "grid", "mesh", "mesh_alt", "sim_newdata", "is_tow"))]),
    failed = failed
  )
}

# Metrics --------------------------------------------------------------------

# Estimated and true indexes joined, with the log error and CI coverage. With
# common_reps, only replicates in which every requested fit (all models and
# arms) succeeded, so models and arms are compared on the same replicates and
# performance is conditional on success (see results$failure_rates).
compare_to_truth <- function(results, ci_levels = c(0.5, 0.95), common_reps = TRUE) {
  indexes <- results$indexes
  if (common_reps) {
    keep <- results$fit_status |>
      group_by(species, rep) |>
      filter(all(status == "ok")) |>
      ungroup() |>
      distinct(species, rep)
    if (nrow(keep) == 0L) {
      cli::cli_abort("No replicate has every requested fit succeeding; see results$failure_rates")
    }
    dropped <- setdiff(unique(indexes$species), keep$species)
    if (length(dropped) > 0L) cli::cli_warn("No replicate with every fit succeeding for: {dropped}")
    n_dropped <- nrow(distinct(anti_join(indexes, keep, by = c("species", "rep")), species, rep))
    if (n_dropped > 0L) cli::cli_inform("Dropped {n_dropped} species x replicate(s) with a failed fit")
    indexes <- semi_join(indexes, keep, by = c("species", "rep"))
  }
  out <- indexes |>
    left_join(results$true_indexes, by = c("species", "rep", "arm", "year")) |>
    left_join(results$occasions, by = c("species", "year")) |>
    mutate(log_error = log_est - log(true_est))
  for (level in ci_levels) {
    q <- qnorm(1 - (1 - level) / 2)
    out[[paste0("covered_", level * 100)]] <- abs(out$log_error) <= q * out$se
  }
  out
}

# Seesaw A (%) over the whole series; phase (+/- 0.5) marks which half was
# sampled
seesaw_A <- function(est, phase, year) {
  period2_metric(est, year = year, phase = phase)[["A"]]
}

# Maximum A across moving windows of `window` survey occasions
seesaw_A_max <- function(est, phase, year, window) {
  max(moving_window(est, year = year, window = window, phase = phase)$A)
}

# Seesaw beta (sampled-half minus other-half log difference, detrended) in
# each moving window
seesaw_beta <- function(est, phase, year, window) {
  moving_window(est, year = year, window = window, phase = phase)$beta
}

# Self-test summary (08, 11) -------------------------------------------------

# RMSE, bias, coverage, seesaw A, and the seesaw error by species x replicate
# x model, the true index's A, and summaries across replicates.
#
# A uses |beta| (and max A a max over windows), so estimation noise alone
# inflates it: an estimated A above the true A is not evidence of a spurious
# seesaw. The seesaw error is signed instead: the estimated minus true beta
# (phase +0.5 minus phase -0.5 occasions, detrended), averaged over windows,
# oriented by the sign of the replicate's true omega contrast between the
# regions (from the simulation, not the fits), and as a percentage. Positive:
# the index reads too high, relative to the truth, in occasions when the
# richer region is sampled. Noise averages to zero.
self_test_metrics <- function(compared, results, window) {
  metrics_rep <- compared |>
    arrange(species, rep, model, year) |>
    group_by(species, rep, model) |>
    summarise(
      rmse = sqrt(mean(log_error^2)),
      bias = mean(log_error),
      coverage_50 = mean(covered_50),
      coverage_95 = mean(covered_95),
      A = seesaw_A(est, phase, year),
      A_max = seesaw_A_max(est, phase, year, window),
      beta_error = mean(seesaw_beta(est, phase, year, window) - seesaw_beta(true_est, phase, year, window)),
      .groups = "drop"
    )

  orientation <- results$true_indexes |>
    group_by(species, rep) |>
    summarise(contrast_omega = mean(contrast_omega), .groups = "drop") |>
    mutate(orientation = if_else(contrast_omega < 0, -1, 1))
  metrics_rep <- metrics_rep |>
    left_join(orientation, by = c("species", "rep")) |>
    mutate(seesaw_error = 100 * (exp(orientation * beta_error) - 1))

  # Same replicates as the fits (see compare_to_truth())
  A_truth <- results$true_indexes |>
    semi_join(distinct(compared, species, rep), by = c("species", "rep")) |>
    left_join(results$occasions, by = c("species", "year")) |>
    arrange(species, rep, year) |>
    group_by(species, rep) |>
    summarise(
      A = seesaw_A(true_est, phase, year),
      A_max = seesaw_A_max(true_est, phase, year, window),
      .groups = "drop"
    )

  metrics <- metrics_rep |>
    group_by(species, model) |>
    summarise(
      n_reps = n(),
      rmse = mean(rmse),
      bias = mean(bias),
      coverage_50 = mean(coverage_50),
      coverage_95 = mean(coverage_95),
      A_median = median(A),
      A_max_median = median(A_max),
      seesaw_error_median = median(seesaw_error),
      .groups = "drop"
    ) |>
    arrange(species, rmse)

  list(
    compared = compared, metrics_rep = metrics_rep, metrics = metrics, A_truth = A_truth,
    failure_rates = results$failure_rates
  )
}

# The three self-test figures: <fig_prefix>.pdf (all metrics by species),
# -summary.pdf (across species), and -indexes.pdf (first three replicates);
# fig_prefix is a full path without the extension
plot_self_test <- function(m, full_model, window, fig_prefix, caption) {
  metric_levels <- c("rmse", "bias", "coverage_50", "coverage_95", "A", "A_max", "seesaw_error")
  metric_labels <- c(
    "RMSE (log)", "Bias (log)", "50% CI coverage", "95% CI coverage",
    "Seesaw A (%)", paste0("Max A, ", window, "-survey windows (%)"),
    "Seesaw error (%)"
  )
  # Nominal values; for A and max A, the median of the true index by species
  nominal <- bind_rows(
    tidyr::expand_grid(
      species = unique(m$metrics_rep$species),
      tibble::tibble(metric = c("bias", "coverage_50", "coverage_95", "seesaw_error"), value = c(0, 0.5, 0.95, 0))
    ),
    m$A_truth |>
      group_by(species) |>
      summarise(A = median(A), A_max = median(A_max)) |>
      tidyr::pivot_longer(c(A, A_max), names_to = "metric")
  ) |>
    mutate(metric = factor(metric, levels = metric_levels, labels = metric_labels))

  # Models ordered by mean RMSE across species
  model_order <- m$metrics |>
    group_by(model) |>
    summarise(rmse = mean(rmse)) |>
    arrange(rmse) |>
    pull(model)

  m$metrics_rep |>
    tidyr::pivot_longer(all_of(metric_levels), names_to = "metric") |>
    mutate(
      model = factor(model, levels = rev(model_order)),
      metric = factor(metric, levels = metric_levels, labels = metric_labels)
    ) |>
    ggplot(aes(value, model)) +
    geom_vline(aes(xintercept = value), data = nominal, linetype = 2, colour = "grey50") +
    geom_boxplot(outlier.shape = NA) +
    geom_point(position = position_jitter(height = 0.15), alpha = 0.4) +
    facet_grid(species ~ metric, scales = "free_x") +
    labs(x = NULL, y = NULL, caption = caption)
  ggsave(
    paste0(fig_prefix, ".pdf"),
    width = 18, height = 1 + 2.2 * length(unique(m$metrics_rep$species))
  )

  # Summary across species: one point per species, each a summary across
  # replicates, scaled so species are comparable; violins show all replicates
  # pooled across species. RMSE is relative to the species' mean full-data
  # operating model RMSE (how much is lost by sampling one half per occasion).
  # The seesaw error is signed and relative to the true index in the same
  # replicate (spurious seesaw from the model, not from the simulated
  # population); see self_test_metrics().
  summary_levels <- c("rmse_ratio", "coverage_50", "seesaw_error")
  summary_labels <- c("RMSE / full-data RMSE", "50% CI coverage", "Seesaw error (%),\nricher-region occasions")

  rep_summary <- m$metrics_rep |>
    group_by(species) |>
    mutate(rmse_ratio = rmse / mean(rmse[model == full_model])) |>
    ungroup() |>
    tidyr::pivot_longer(all_of(summary_levels), names_to = "metric") |>
    mutate(
      model = factor(model, levels = rev(model_order)),
      metric = factor(metric, levels = summary_levels, labels = summary_labels)
    )

  species_summary <- rep_summary |>
    group_by(species, model, metric) |>
    summarise(
      value = if (first(metric) == summary_labels[3]) median(value) else mean(value),
      .groups = "drop"
    )

  summary_nominal <- tibble::tibble(
    metric = factor(summary_labels, levels = summary_labels),
    value = c(1, 0.5, 0)
  )

  # Zoom each panel to the widest 10-90% quantile range of the replicates
  # across models, extended to keep the species points and nominal values in
  # view. oob_keep so the violins are still computed from all replicates (as
  # with coord_cartesian(), which can't be set per panel).
  summary_limits <- bind_rows(
    rep_summary |>
      group_by(metric, model) |>
      summarise(lo = quantile(value, 0.1), hi = quantile(value, 0.9), .groups = "drop"),
    rename(species_summary, lo = value) |> mutate(hi = lo),
    rename(summary_nominal, lo = value) |> mutate(hi = lo)
  ) |>
    group_by(metric) |>
    summarise(lo = min(lo), hi = max(hi)) |>
    arrange(metric)
  summary_x_scales <- purrr::map2(
    summary_limits$lo, summary_limits$hi,
    \(lo, hi) scale_x_continuous(limits = c(lo, hi), oob = scales::oob_keep)
  )

  ggplot(species_summary, aes(value, model)) +
    geom_violin(data = rep_summary, fill = "grey90", colour = NA, scale = "width") +
    geom_vline(aes(xintercept = value), data = summary_nominal, linetype = 2, colour = "grey50") +
    stat_summary(fun = median, geom = "crossbar", fun.min = median, fun.max = median, width = 0.6) +
    geom_point(aes(colour = species), position = position_jitter(height = 0.15, seed = 1)) +
    facet_wrap(~metric, nrow = 1, scales = "free_x") +
    ggh4x::facetted_pos_scales(x = summary_x_scales) +
    # Paired has 12 colours (Dark2 has only 8)
    scale_colour_brewer(palette = "Paired", labels = tools::toTitleCase) +
    labs(x = NULL, y = NULL, colour = NULL) +
    theme(legend.position = "bottom", plot.margin = margin(5.5, 12, 5.5, 5.5)) +
    guides(colour = guide_legend(nrow = 2))
  ggsave(paste0(fig_prefix, "-summary.pdf"), width = 9, height = 4)

  m$compared |>
    filter(rep <= 3) |>
    ggplot(aes(cal_year, est, colour = model)) +
    geom_line(aes(y = true_est), colour = "black") +
    geom_line() +
    facet_grid(species ~ rep, scales = "free_y", labeller = labeller(rep = label_both)) +
    labs(x = "Year", y = "Biomass index", colour = "Model", caption = "Black: true index")
  ggsave(
    paste0(fig_prefix, "-indexes.pdf"),
    width = 11, height = 1 + 2 * length(unique(m$compared$species))
  )
  invisible(NULL)
}
