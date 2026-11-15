# Why does the alternating-survey estimator keep an HS-low seesaw when the
# generating omega is flattened? Stages 0, 1, and 2A/2B of
# analysis/omega-swap-experiment-plan.md (cod only, few replicates; later
# stages only if these implicate them).
#
# Generating model: as 13's conditional_flat (omega shifted so the HS and QCS
# area-weighted means are equal, per delta component) with epsilon exactly
# zero, the real year effects, and the observation model. One simulation per
# replicate; every arm below is fit to the same catches:
#
#   A alternating, omega estimated, IID epsilon estimated (reference)
#   B alternating, omega estimated, epsilon off
#   C alternating, omega fixed at the generating field, epsilon off
#   D full data,   omega estimated, epsilon off
#   E full data,   omega fixed at the generating field, epsilon off
#
# Stage 2 (omega estimated, epsilon off, same catches; tow subsets chosen
# independently of the catches, seeded by replicate):
#
#   F overlap10: each occasion keeps its tow total, but 10% of the primary
#     region's tows are replaced by tows from the other region
#   G overlap25: as F, 25%
#   H overlap_area: each occasion's tow total split between regions in
#     proportion to area
#   I reversed: the complementary alternating subset (QCS on odd occasions,
#     HS on even); phase follows the region actually sampled
#
# Constant-omega control (epsilon off): the same replicate simulated with each
# component's omega a constant (its coast-wide area-weighted mean), omega
# estimated except in M:
#
#   J alternating, K reversed, L 25% overlap, N full data (estimated omega)
#   M alternating, spatial field off (the known constant is absorbed by the
#     year effects)
#
# Other species: B, C, G, J only.
#
# The primary outcome is raw (not oriented): the estimated minus true phase
# coefficient of log index, HS = +0.5 and QCS = -0.5, averaged over the
# 8-occasion windows. Below 0 (ratio below 1): HS occasions read low.
#
# Fixed omega: the generating nodes are the starting values of omega_s, which
# is mapped off (so it is not estimated and receives no shrinkage). The year
# effects are still estimated. Checked below against the simulation's field at
# the grid cells.

source(here::here("analysis/self-test-funcs.R"))
source(here::here("analysis/omega-swap-funcs.R"))
theme_set(ggsidekick::theme_sleek())

surveyjoin::cache_data()
surveyjoin::load_sql_data()

species_to_fit <- c("pacific cod", "spotted ratfish", "pacific halibut")
n_reps <- c("pacific cod" = 30L, "spotted ratfish" = 20L, "pacific halibut" = 20L)
n_workers <- max(1L, min(8L, future::availableCores() %/% 2L))
seesaw_window <- 8L
tolerance <- log(1.05)
est_model <- "IID RF, factor(year)"

# Cache roots. Setups and simulations are keyed by the generating
# configuration and the code that makes it; fits by the code that fits them
# (hashes of deparsed functions, so comments and formatting don't matter).
# Both are computed after the functions are defined (below).
config <- list(
  generating = "conditional_flat, epsilon exactly 0", mesh_cutoff = 10,
  est_model = est_model, seeds = c(seed_pars, seed_sim, seed_extra),
  sdmTMB = as.character(packageVersion("sdmTMB"))
)

arms <- tibble::tribble(
  ~arm, ~generating, ~design,      ~omega,      ~epsilon,
  "A",  "pattern",   "alternating", "estimated", "iid",
  "B",  "pattern",   "alternating", "estimated", "off",
  "C",  "pattern",   "alternating", "fixed",     "off",
  "D",  "pattern",   "full",        "estimated", "off",
  "E",  "pattern",   "full",        "fixed",     "off",
  "F",  "pattern",   "overlap10",   "estimated", "off",
  "G",  "pattern",   "overlap25",   "estimated", "off",
  "H",  "pattern",   "overlap_area", "estimated", "off",
  "I",  "pattern",   "reversed",    "estimated", "off",
  # spatially constant generating omega
  "J",  "constant",  "alternating", "estimated", "off",
  "K",  "constant",  "reversed",    "estimated", "off",
  "L",  "constant",  "overlap25",   "estimated", "off",
  "M",  "constant",  "alternating", "none",      "off",
  "N",  "constant",  "full",        "estimated", "off"
)
# Arms fit in each species (cod: all; others: the decisive subset)
species_arms <- \(species) if (species == "pacific cod") arms$arm else c("B", "C", "G", "J")
# Stage 1 arms are the first five
arm_levels <- arms$arm

# Rows of the tow data used by each design for replicate `rep`. Overlap
# designs keep each occasion's number of tows (the primary region's count) and
# take tows within depth-sorted systematic samples, independent of catches.
design_rows <- function(setup, design, rep) {
  d <- setup$dat
  if (design == "alternating") return(setup$is_sampled)
  if (design == "reversed") return(!setup$is_sampled)
  if (design == "full") return(rep(TRUE, nrow(d)))
  set.seed(seed_extra + 1e5L + rep)
  # systematic sample of n rows over depth, random start
  pick <- function(rows, n) {
    if (n <= 0L) return(integer(0))
    rows <- rows[order(d$depth_m[rows])]
    stopifnot(n <= length(rows))
    rows[floor((seq_len(n) - 1 + runif(1)) * length(rows) / n) + 1L]
  }
  area <- tapply(setup$grid$area[setup$grid$year == 1L], setup$grid$hs[setup$grid$year == 1L], sum)
  keep <- integer(0)
  for (y in sort(unique(d$year))) {
    primary <- if (d$primary_survey[d$year == y][1] == "SYN HS") 1L else 0L
    n <- sum(d$year == y & d$hs == primary)
    n_other <- switch(design,
      overlap10 = round(0.1 * n), overlap25 = round(0.25 * n),
      overlap_area = round(n * area[[as.character(1L - primary)]] / sum(area))
    )
    keep <- c(
      keep,
      pick(which(d$year == y & d$hs == primary), n - n_other),
      pick(which(d$year == y & d$hs != primary), n_other)
    )
  }
  seq_len(nrow(d)) %in% keep
}
# Phase (+0.5 when HS is the primary region) of each occasion under a design
design_phase <- function(design, phase) if (design == "reversed") -phase else phase

# Setup and simulation -------------------------------------------------------

setup_species <- function(species, dat, grid) {
  setup <- qcs_hs_setup(species, dat, grid)
  om <- fit_operating_model(est_model, dat, setup$mesh)
  if (isTRUE(om$om_failed)) return(c(list(species = species), om))
  nm <- names(om$fit$tmb_obj$env$last.par.best)
  idx <- which(nm == "omega_s")
  n_nodes <- nrow(setup$mesh$mesh$loc)
  omega_map <- om$fit$tmb_obj$env$map$omega_s
  kept <- if (is.null(omega_map)) seq_along(idx) else which(!is.na(omega_map))
  component <- (kept - 1L) %/% n_nodes + 1L
  if (n_distinct(component) < 2L) {
    return(list(species = species, om_failed = TRUE, error = "a spatial field collapsed"))
  }
  stopifnot(length(kept) == length(idx), length(idx) == n_nodes * 2L)
  c(setup, list(
    fit_om = om$fit, ran_pars = om$ran_pars,
    omega_idx = split(idx, component),
    eps_idx = which(nm == "epsilon_st"),
    b_idx = list(which(nm == "b_j"), which(nm == "b_j2")),
    shifter = make_region_shifter(setup$mesh, grid)
  ))
}

# One replicate with epsilon exactly zero, simulated twice from the same
# parameter draw and seed (so the same year effects, tow locations and
# efforts; observation errors share the random stream but not the realization):
#   pattern:  the flattened fitted omega (equal HS and QCS area-weighted means
#     per delta component, within-region pattern kept)
#   constant: each component's omega a constant, its coast-wide area-weighted
#     mean (the same constant the flattened field is centred on)
# simulate_re = character(0) keeps every supplied random effect, so the field
# is the one in `pars`.
simulate_rep <- function(setup, sampler, rep) {
  pars0 <- draw_pars(sampler, rep)
  pars <- list(pattern = pars0, constant = pars0)
  for (j in setup$omega_idx) {
    centre <- setup$shifter$overall(pars0[j, 1L])
    pars$pattern[j, ] <- setup$shifter$shift(pars0[j, 1L], rep(centre, 2L))
    pars$constant[j, ] <- centre
  }
  lapply(pars, \(p) simulate_one(setup, p, rep))
}

simulate_one <- function(setup, pars, rep) {
  pars[setup$eps_idx, ] <- 0
  report <- simulate_pars(setup, pars, rep, simulate_re = character(0))
  grid_rows <- !setup$is_tow
  g <- setup$grid
  stopifnot(max(abs(report$epsilon_st)) == 0)
  # eta = year effect + omega (nothing else enters the generating predictor)
  b <- lapply(setup$b_idx, \(i) pars[i, 1L])
  for (k in 1:2) {
    stopifnot(isTRUE(all.equal(
      report$eta_i[grid_rows, k], report$omega_s_A[grid_rows, k] + b[[k]][g$year]
    )))
  }
  truth <- g |>
    mutate(density = exp(report$eta_i[grid_rows, 1L] + report$eta_i[grid_rows, 2L]) * area)
  catch <- rep_catch(setup, report)
  list(
    catch = catch,
    truth_total = summarise(truth, true_est = sum(density), .by = year),
    truth_region = summarise(truth, true_est = sum(density), .by = c(year, hs)),
    # generating fields: nodes by delta component (for the fixed-omega fits),
    # and omega at tows and grid cells (each component)
    omega_nodes = do.call(cbind, lapply(setup$omega_idx, \(j) pars[j, 1L])),
    omega_rows = report$omega_s_A,
    year_effects = tibble::tibble(year = rep(seq_along(b[[1]]), 2L), component = rep(1:2, each = length(b[[1]])), true = unlist(b)),
    eps_max_abs = max(abs(report$epsilon_st)),
    true_region_omega = purrr::map_dbl(1:2, \(k) {
      w <- g$area[g$year == 1L]
      o <- report$omega_s_A[grid_rows, k][g$year == 1L]
      h <- g$hs[g$year == 1L]
      weighted.mean(o[h == 1L], w[h == 1L]) - weighted.mean(o[h == 0L], w[h == 0L])
    }),
    # mean true density per area and encounter rate of the simulated tows, by
    # region (removing the pattern changes both even at the same mean log omega)
    region_summary = truth |>
      summarise(density_per_area = sum(density) / sum(area), .by = c(year, hs)) |>
      summarise(density_per_area = mean(density_per_area), .by = hs) |>
      left_join(
        tibble::tibble(hs = setup$dat$hs, encounter = catch > 0) |> summarise(encounter_rate = mean(encounter), .by = hs),
        by = "hs"
      )
  )
}

# Fits -----------------------------------------------------------------------

fit_arm <- function(setup, sim, arm_row, rep, fit_file = NULL) {
  d <- setup$dat
  d$catch_weight <- sim$catch
  rows <- design_rows(setup, arm_row$design, rep)
  tow_rows <- which(rows)
  d <- d[tow_rows, ]
  # same mesh nodes for every design
  mesh <- make_mesh(d, c("X", "Y"), mesh = setup$mesh$mesh)
  control <- st_control
  spatial <- if (arm_row$omega == "none") "off" else "on"
  if (arm_row$omega == "fixed") {
    control <- sdmTMBcontrol(
      start = list(omega_s = sim$omega_nodes),
      map = list(omega_s = factor(rep(NA, length(sim$omega_nodes)))),
      collapse_spatial_variance = FALSE
    )
  }
  fit <- tryCatch(
    sdmTMB(
      catch_weight ~ 0 + factor(year), data = d, mesh = mesh, offset = log(d$effort),
      family = st_family, time = "year", spatial = spatial, spatiotemporal = arm_row$epsilon,
      share_range = spatial == "on", anisotropy = spatial == "on", control = control, silent = TRUE
    ),
    error = \(e) e
  )
  if (!inherits(fit, "sdmTMB")) return(list(status = "fit error", message = conditionMessage(fit)))
  checks <- tryCatch(sanity(fit, gradient_thresh = 0.01, silent = TRUE), error = \(e) list(all_ok = FALSE))
  ran_pars <- bind_rows(
    tryCatch(mutate(tidy(fit, "ran_pars", model = 1), model = 1), error = \(e) NULL),
    tryCatch(mutate(tidy(fit, "ran_pars", model = 2), model = 2), error = \(e) NULL)
  )
  out <- list(status = if (isTRUE(checks$all_ok)) "ok" else "failed sanity", sanity = unlist(checks), ran_pars = ran_pars)
  # With a constant generating field there is nothing for the spatial field to
  # fit (its variance and range are unidentified), so sanity() flags nearly
  # every fit. These are kept, labelled "degenerate", if the index can be
  # computed; summaries report sanity-ok fits separately
  if (out$status == "failed sanity" && arm_row$generating == "constant") out$status <- "degenerate"
  if (!out$status %in% c("ok", "degenerate")) return(out)
  if (!is.null(fit_file)) saveRDS(fit, fit_file)

  g <- setup$grid
  index <- get_index_ok(fit, g, nsplit = 2L)
  if (is.null(index)) return(c(list(status = "index error"), out[-1L]))
  # the seesaw uses only the point estimates; a degenerate field can leave the SEs NaN
  if (!all(is.finite(index$est))) return(c(list(status = "index error"), out[-1L]))
  # Uncorrected regional (and total) biomass from the predictions
  p <- predict(fit, newdata = g, offset = rep(0, nrow(g)))
  p$density <- exp(p$est1 + p$est2) * p$area
  # Fitted fields at the grid and at the fitted tows
  p_tow <- predict(fit)
  a1 <- g$year == min(g$year)
  # a collapsed (mapped-off) field has no omega column
  col0 <- \(x, nm) if (nm %in% names(x)) x[[nm]] else rep(0, nrow(x))
  grid_rows <- seq_len(nrow(sim$omega_rows)) > nrow(setup$dat)
  region_means <- \(x, h, w) c(weighted.mean(x[h == 1L], w[h == 1L]) - weighted.mean(x[h == 0L], w[h == 0L]))
  omega_contrast <- tibble::tibble(
    component = 1:2,
    grid = c(
      region_means(col0(p, "omega_s1")[a1], g$hs[a1], g$area[a1]),
      region_means(col0(p, "omega_s2")[a1], g$hs[a1], g$area[a1])
    ),
    tow = c(
      region_means(col0(p_tow, "omega_s1"), d$hs, rep(1, nrow(d))),
      region_means(col0(p_tow, "omega_s2"), d$hs, rep(1, nrow(d)))
    ),
    # true, same grid and tows
    grid_true = c(
      region_means(sim$omega_rows[grid_rows, 1L][a1], g$hs[a1], g$area[a1]),
      region_means(sim$omega_rows[grid_rows, 2L][a1], g$hs[a1], g$area[a1])
    ),
    tow_true = c(
      region_means(sim$omega_rows[tow_rows, 1L], d$hs, rep(1, nrow(d))),
      region_means(sim$omega_rows[tow_rows, 2L], d$hs, rep(1, nrow(d)))
    )
  )
  year_effects <- bind_rows(lapply(1:2, \(k) {
    tibble::tibble(
      year = sort(unique(d$year)), component = k,
      est = as.numeric(fit$tmb_obj$env$last.par.best[names(fit$tmb_obj$env$last.par.best) == c("b_j", "b_j2")[k]])
    )
  }))
  c(out, list(
    index = index,
    region = summarise(p, est = sum(density), .by = c(year, hs)),
    total_uncorrected = summarise(p, est = sum(density), .by = year),
    omega_contrast = omega_contrast,
    year_effects = year_effects,
    n_tows = count(d, year, hs),
    n_tows_total = nrow(d)
  ))
}

setup_worker <- function(species, dat, grid, file, log_file) {
  single_thread()
  run_logged(log_file, paste0("setup: ", species), {
    saveRDS(setup_species(species, dat, grid), file)
    NULL
  })
}
sim_worker <- function(setup, rep, file, log_file) {
  single_thread()
  run_logged(log_file, paste0("sim: ", setup$species, " rep ", rep), {
    setup$fit_om <- sdmTMB:::reload_model(setup$fit_om)
    saveRDS(simulate_rep(setup, sdmTMB:::.posterior_re_sampler(setup$fit_om), rep), file)
    NULL
  })
}
fit_worker <- function(setup, sim, arm_row, rep, file, fit_file, log_file) {
  single_thread()
  run_logged(log_file, paste0("fit: ", setup$species, " ", arm_row$arm), {
    res <- fit_arm(setup, sim, arm_row, rep, if (rep <= 2L) fit_file)
    saveRDS(res, file)
    if (res$status != "ok") paste0(" ", toupper(res$status))
  })
}

# Run ------------------------------------------------------------------------

grid <- qcs_hs_grid()
code_hash <- \(...) substr(rlang::hash(lapply(list(...), deparse)), 1, 8)
gen_id <- substr(rlang::hash(list(config, code_hash(setup_species, simulate_rep))), 1, 8)
fit_id <- code_hash(fit_arm, design_rows, st_control, st_family)
root <- here::here("data-generated/qcs-hs-residual-seesaw")
dir <- file.path(root, paste0("gen-", gen_id))
fit_dir <- file.path(dir, paste0("fits-", fit_id))
for (d in c(file.path(dir, c("setups", "sims")), fit_dir, file.path(fit_dir, "fitted-objects"))) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}
saveRDS(config, file.path(dir, "config.rds"))
log_file <- file.path(dir, "progress.log")
future::plan(future::multisession, workers = n_workers)
opts <- furrr::furrr_options(seed = TRUE, scheduling = Inf)
log_progress(log_file, "=== run started: ", sum(n_reps), " reps, ", n_workers, " workers ===")

setup_file <- \(s) file.path(dir, "setups", paste0(slug(s), ".rds"))
todo <- species_to_fit[!file.exists(setup_file(species_to_fit))]
if (length(todo) > 0L) {
  furrr::future_pwalk(
    list(species = todo, dat = lapply(todo, qcs_hs_data), file = setup_file(todo)),
    setup_worker, grid = grid, log_file = log_file, .options = opts
  )
}
setups <- setNames(lapply(species_to_fit, \(s) readRDS(setup_file(s))), species_to_fit)
stopifnot(!purrr::some(setups, \(x) isTRUE(x$om_failed)))

sims <- purrr::map_dfr(species_to_fit, \(s) tibble::tibble(species = s, rep = seq_len(n_reps[[s]]))) |>
  mutate(file = file.path(dir, "sims", sprintf("%s-rep%03d.rds", slug(species), rep)))
todo <- filter(sims, !file.exists(file))
if (nrow(todo) > 0L) {
  furrr::future_pwalk(
    list(setup = setups[todo$species], rep = todo$rep, file = todo$file),
    sim_worker, log_file = log_file, .options = opts
  )
}
sims <- filter(sims, file.exists(file))
sim_list <- setNames(lapply(sims$file, readRDS), paste(sims$species, sims$rep))

tasks <- tidyr::expand_grid(sims[c("species", "rep")], arms) |>
  filter(purrr::map2_lgl(species, arm, \(s, a) a %in% species_arms(s))) |>
  mutate(file = file.path(fit_dir, sprintf("%s-rep%03d-%s.rds", slug(species), rep, arm)),
    fit_file = file.path(fit_dir, "fitted-objects", sprintf("%s-rep%03d-%s.rds", slug(species), rep, arm)))
todo <- filter(tasks, !file.exists(file))
if (nrow(todo) > 0L) {
  lite <- lapply(setups, \(x) x[c("species", "dat", "is_sampled", "grid", "mesh", "mesh_alt")])
  furrr::future_pwalk(
    list(
      setup = lite[todo$species],
      sim = purrr::pmap(todo[c("species", "rep", "generating")], \(species, rep, generating) sim_list[[paste(species, rep)]][[generating]]),
      arm_row = purrr::map(seq_len(nrow(todo)), \(i) todo[i, names(arms)]),
      rep = todo$rep,
      fit_file = todo$fit_file,
      file = todo$file
    ),
    fit_worker, log_file = log_file, .options = opts
  )
}
log_progress(log_file, "=== run finished ===")
future::plan(future::sequential)

fits <- setNames(lapply(tasks$file, \(f) if (file.exists(f)) readRDS(f)), paste(tasks$species, tasks$rep, tasks$arm))
tasks$status <- vapply(fits, \(x) if (is.null(x)) "missing" else x$status, character(1L))

# Metrics --------------------------------------------------------------------

occasions <- purrr::map_dfr(setups, \(s) select(s$occasions, year, phase), .id = "species")
fit_status <- tasks |>
  count(species, arm, status) |>
  tidyr::pivot_wider(names_from = status, values_from = n, values_fill = 0L)
print(fit_status)

# Phase coefficient of log(est / true) in each window, averaged; also the
# whole-series detrended coefficient and the equal 4-versus-4 windows without
# detrending (sensitivity, Stage 0 item 4)
seesaw_metrics <- function(est, true, years, phase) {
  ratio <- est / true
  beta_windows <- seesaw_beta(ratio, phase, years, seesaw_window)
  # subtraction used in 13, to confirm agreement
  beta_sub <- seesaw_beta(est, phase, years, seesaw_window) - seesaw_beta(true, phase, years, seesaw_window)
  fourfour <- vapply(seq(min(years), max(years) - seesaw_window + 1L), \(y0) {
    i <- years >= y0 & years < y0 + seesaw_window
    mean(log(ratio[i & phase > 0])) - mean(log(ratio[i & phase < 0]))
  }, numeric(1L))
  c(
    seesaw = mean(beta_windows), seesaw_sub = mean(beta_sub),
    seesaw_whole = period2_metric(ratio, year = years, phase = phase)[["beta"]],
    seesaw_4v4 = mean(fourfour)
  )
}

per_fit <- purrr::pmap_dfr(tasks[c("species", "rep", "arm", "design", "generating")], \(species, rep, arm, design, generating) {
  f <- fits[[paste(species, rep, arm)]]
  if (is.null(f) || !f$status %in% c("ok", "degenerate")) return(NULL)
  sim <- sim_list[[paste(species, rep)]][[generating]]
  occ <- occasions |>
    filter(species == !!species) |>
    mutate(phase = design_phase(design, phase)) |>
    select(year, phase)
  tot <- f$index |>
    select(year, est, log_est, se) |>
    left_join(sim$truth_total, by = "year") |>
    left_join(occ, by = "year") |>
    arrange(year)
  m <- seesaw_metrics(tot$est, tot$true_est, tot$year, tot$phase)
  # the same metric from the uncorrected (summed predictions) totals: how much
  # of the seesaw the bias correction adds
  unc <- f$total_uncorrected |> arrange(year)
  m_unc <- seesaw_metrics(unc$est, tot$true_est, tot$year, tot$phase)[["seesaw"]]
  tot <- mutate(tot, log_error = log(est / true_est))
  # regional errors: the sampled region in each occasion is HS when phase > 0
  reg <- f$region |>
    left_join(sim$truth_region, by = c("year", "hs")) |>
    left_join(occ, by = "year") |>
    # the primary region of each occasion (not meaningful for full data)
    mutate(sampled = (hs == 1L) == (phase > 0), log_error = log(est / true_est))
  year_err <- f$year_effects |>
    left_join(sim$year_effects, by = c("year", "component")) |>
    # with the field off, the known constant is part of the year effect
    mutate(true = if (arms$omega[arms$arm == arm] == "none") true + sim$omega_nodes[1L, component] else true, err = est - true) |>
    left_join(occ, by = "year")
  ye_phase <- year_err |>
    summarise(phase_err = coef(lm(err ~ year + phase))[["phase"]], .by = component)
  oc <- f$omega_contrast
  tibble::tibble(
    species = species, rep = rep, arm = arm, generating = generating, sanity_ok = f$status == "ok",
    !!!as.list(m), seesaw_uncorrected = m_unc, correction_effect = m[["seesaw"]] - m_unc,
    # mean log index error by occasion phase
    err_hs_occ = mean(tot$log_error[tot$phase > 0]), err_qcs_occ = mean(tot$log_error[tot$phase < 0]),
    # sampled and unsampled region biomass errors, averaged over occasions
    err_sampled = mean(reg$log_error[reg$sampled]), err_unsampled = mean(reg$log_error[!reg$sampled]),
    err_hs_when_unsampled = mean(reg$log_error[reg$hs == 1L & !reg$sampled]),
    err_qcs_when_unsampled = mean(reg$log_error[reg$hs == 0L & !reg$sampled]),
    # unsampled-region error when QCS is unsampled (HS occasions) minus when HS
    # is (QCS occasions): the direction of the compensation
    unsampled_contrast = mean(reg$log_error[reg$hs == 0L & !reg$sampled]) - mean(reg$log_error[reg$hs == 1L & !reg$sampled]),
    err_hs_region = mean(reg$log_error[reg$hs == 1L]), err_qcs_region = mean(reg$log_error[reg$hs == 0L]),
    # omega contrast errors (estimated minus true HS - QCS; the true grid contrast is 0)
    omega_err_grid1 = oc$grid[1] - oc$grid_true[1], omega_err_grid2 = oc$grid[2] - oc$grid_true[2],
    omega_err_grid = sum(oc$grid - oc$grid_true),
    omega_err_tow = sum(oc$tow - oc$tow_true),
    ye_phase_err1 = ye_phase$phase_err[1], ye_phase_err2 = ye_phase$phase_err[2],
    ye_phase_err = sum(ye_phase$phase_err),
    # Stage 0: log(est) vs the returned log_est (bias-corrected), and
    # regions summed vs the total
    log_est_gap = mean(log(tot$est) - tot$log_est),
    region_sum_gap = mean(log(f$total_uncorrected$est[order(f$total_uncorrected$year)]) - log(tot$est)),
    true_region_sum_gap = max(abs(log(summarise(sim$truth_region, s = sum(true_est), .by = year)$s) - log(sim$truth_total$true_est))),
    n_occasions = nrow(tot), n_tows_total = f$n_tows_total, phases_ok = identical(tot$phase, design_phase(design, rep(c(0.5, -0.5), length.out = nrow(tot)))),
    sigma_O1 = if ("sigma_O" %in% f$ran_pars$term) f$ran_pars$estimate[f$ran_pars$term == "sigma_O" & f$ran_pars$model == 1][1] else NA_real_,
    sigma_O2 = if ("sigma_O" %in% f$ran_pars$term) f$ran_pars$estimate[f$ran_pars$term == "sigma_O" & f$ran_pars$model == 2][1] else NA_real_,
    sigma_E1 = if ("sigma_E" %in% f$ran_pars$term) f$ran_pars$estimate[f$ran_pars$term == "sigma_E" & f$ran_pars$model == 1][1] else NA_real_,
    sigma_E2 = if ("sigma_E" %in% f$ran_pars$term) f$ran_pars$estimate[f$ran_pars$term == "sigma_E" & f$ran_pars$model == 2][1] else NA_real_
  )
})

# Stage 0 checks
cat("\n== Stage 0 checks ==\n")
all_sims <- unlist(sim_list, recursive = FALSE)
cat("max |simulated epsilon|:", max(vapply(all_sims, \(s) s$eps_max_abs, numeric(1L))), "\n")
cat("max |true HS - QCS omega contrast| (grid, per component):",
    max(abs(unlist(lapply(all_sims, `[[`, "true_region_omega")))), "\n")
# Constant vs pattern: density and encounter rates (mean over replicates)
region_docs <- purrr::imap_dfr(sim_list, \(s, nm) {
  purrr::imap_dfr(s, \(x, generating) mutate(x$region_summary, generating = generating, key = nm))
}) |>
  tidyr::separate(key, c("species", "rep"), sep = " (?=[0-9]+$)") |>
  summarise(density_per_area = mean(density_per_area), encounter_rate = mean(encounter_rate), .by = c(species, generating, hs)) |>
  mutate(region = if_else(hs == 1L, "HS", "QCS"), hs = NULL) |>
  arrange(species, region, generating)
cat("\nGenerating processes: mean true density per area and tow encounter rate\n")
print(as.data.frame(region_docs))
cat("max |seesaw (log ratio) - seesaw (subtracted betas)|:", max(abs(per_fit$seesaw - per_fit$seesaw_sub)), "\n")
cat("all fits have 11 occasions in the design's phase order:", all(per_fit$n_occasions == 11L & per_fit$phases_ok), "\n")
cat("max |true regions summed - true total| (log):", max(per_fit$true_region_sum_gap), "\n")
per_fit |>
  group_by(species, arm) |>
  summarise(
    across(c(log_est_gap, region_sum_gap), \(x) mean(x)),
    max_abs_region_sum_gap = max(abs(region_sum_gap)),
    n_tows_total = mean(n_tows_total),
    .groups = "drop"
  ) |>
  print()
cat("Tow totals per replicate by arm (alternating, overlap, and reversed designs have different totals by occasion parity; equal budgets are within each occasion)\n")

# Summaries --------------------------------------------------------------------

# Each arm uses its own successful replicates; each paired contrast uses the
# replicates where both of its arms succeeded (so a failure-prone arm doesn't
# thin the others)
pf <- per_fit

mc <- \(x) c(mean = mean(x), se = sd(x) / sqrt(length(x)))
arm_summary <- pf |>
  group_by(species, arm) |>
  summarise(
    n = n(), n_sanity_ok = sum(sanity_ok),
    mean_log = mean(seesaw), se_log = sd(seesaw) / sqrt(n),
    gm_ratio = exp(mean_log), lwr = exp(mean_log - 1.96 * se_log), upr = exp(mean_log + 1.96 * se_log),
    median_ratio = exp(median(seesaw)), q25 = exp(quantile(seesaw, 0.25)), q75 = exp(quantile(seesaw, 0.75)),
    prop_below_1 = mean(seesaw < 0),
    generating = first(generating),
    sigma_O1 = median(sigma_O1, na.rm = TRUE), sigma_O2 = median(sigma_O2, na.rm = TRUE),
    # sensitivity
    whole_series = mean(seesaw_whole), four_v_four = mean(seesaw_4v4),
    err_hs_occ = mean(err_hs_occ), err_qcs_occ = mean(err_qcs_occ),
    err_sampled = mean(err_sampled), err_unsampled = mean(err_unsampled),
    unsampled_contrast = mean(unsampled_contrast),
    seesaw_uncorrected = mean(seesaw_uncorrected), correction_effect = mean(correction_effect),
    omega_err_grid = mean(omega_err_grid), ye_phase_err = mean(ye_phase_err),
    err_hs_region = mean(err_hs_region), err_qcs_region = mean(err_qcs_region),
    verdict = case_when(
      mean_log + 1.96 * se_log < tolerance & mean_log - 1.96 * se_log > -tolerance ~ "negligible (within +/- log 1.05)",
      mean_log - 1.96 * se_log > 0 | mean_log + 1.96 * se_log < 0 ~ "materially nonzero / unresolved",
      .default = "unresolved"
    ),
    .groups = "drop"
  )
cat("\n== Seesaw error by arm (log phase coefficient; ratio below 1: HS occasions read low) ==\n")
print(arm_summary |> mutate(across(where(is.numeric), \(x) round(x, 3))), n = Inf, width = Inf)

contrasts <- list(
  "A - B (epsilon estimated vs off)" = c("A", "B"),
  "B - C (omega estimated vs fixed)" = c("B", "C"),
  "B - D (alternating vs full, omega estimated)" = c("B", "D"),
  "C - E (alternating vs full, omega fixed)" = c("C", "E"),
  "D - E (full: omega estimated vs fixed)" = c("D", "E"),
  "F - B (10% overlap vs alternating)" = c("F", "B"),
  "G - B (25% overlap vs alternating)" = c("G", "B"),
  "H - B (area-proportional vs alternating)" = c("H", "B"),
  "I - B (reversed vs alternating)" = c("I", "B"),
  "G - D (25% overlap vs full)" = c("G", "D"),
  # constant omega
  "J - B (constant vs patterned omega, alternating)" = c("J", "B"),
  "K - I (constant vs patterned omega, reversed)" = c("K", "I"),
  "L - G (constant vs patterned omega, 25% overlap)" = c("L", "G"),
  "M - J (constant omega: fixed vs estimated)" = c("M", "J"),
  "J - N (constant omega: alternating vs full)" = c("J", "N"),
  "L - J (constant omega: 25% overlap vs alternating)" = c("L", "J"),
  "K - J (constant omega: reversed vs alternating)" = c("K", "J")
)
wide <- pf |> select(species, rep, arm, seesaw) |> tidyr::pivot_wider(names_from = arm, values_from = seesaw)
# a contrast of arms absent from the data is skipped
contrasts <- purrr::keep(contrasts, \(x) all(x %in% names(wide)))
paired <- purrr::imap_dfr(contrasts, \(x, name) {
  wide |>
    mutate(contrast = name, effect = .data[[x[1]]] - .data[[x[2]]]) |>
    filter(!is.na(effect)) |>
    group_by(species, contrast) |>
    summarise(
      n = n(), mean_effect = mean(effect), se = sd(effect) / sqrt(n),
      lwr = mean_effect - 1.96 * se, upr = mean_effect + 1.96 * se,
      ratio = exp(mean_effect), prop_positive = mean(effect > 0),
      .groups = "drop"
    )
})
cat("\n== Paired differences in log seesaw error (mean, 95% MC interval; tolerance +/- 0.049) ==\n")
print(paired |> mutate(across(where(is.numeric), \(x) round(x, 3))), n = Inf, width = Inf)

# Mechanism: do omega-contrast errors trade off against year-effect phase
# errors, and do they predict unsampled-region biomass error? (alternating
# arms with estimated omega: A, B)
mech <- pf |>
  filter(arm %in% c("A", "B", "D", "F", "G", "H", "I", "J", "K", "L", "N")) |>
  group_by(species, arm) |>
  summarise(
    mean_omega_err_grid = mean(omega_err_grid), mean_omega_err_tow = mean(omega_err_tow),
    mean_ye_phase_err = mean(ye_phase_err),
    mean_sigma_O = mean(c(sigma_O1, sigma_O2), na.rm = TRUE),
    cor_omega_ye = cor(omega_err_grid, ye_phase_err),
    cor_omega_unsampled = cor(omega_err_grid, err_unsampled),
    cor_omega_seesaw = cor(omega_err_grid, seesaw),
    cor_ye_seesaw = cor(ye_phase_err, seesaw),
    cor_unsampled_seesaw = cor(err_unsampled, seesaw),
    cor_unsampled_contrast_seesaw = cor(unsampled_contrast, seesaw),
    cor_omega_unsampled_contrast = cor(omega_err_grid, unsampled_contrast),
    .groups = "drop"
  )
cat("\n== Mechanism: estimated omega contrast error (HS - QCS; true 0), year-effect phase error ==\n")
print(mech |> mutate(across(where(is.numeric), \(x) round(x, 3))), width = Inf)
est_sigma <- pf |> filter(arm == "A") |> summarise(across(c(sigma_E1, sigma_E2), \(x) median(x, na.rm = TRUE)), .by = species)
cat("Median fitted sigma_E in arm A (components 1, 2):\n"); print(est_sigma)

# Outputs ----------------------------------------------------------------------

dir.create(here::here("data-outputs"), showWarnings = FALSE)
write.csv(arm_summary, here::here("data-outputs/residual-seesaw-arms.csv"), row.names = FALSE)
write.csv(paired, here::here("data-outputs/residual-seesaw-paired.csv"), row.names = FALSE)
write.csv(mech, here::here("data-outputs/residual-seesaw-mechanism.csv"), row.names = FALSE)
saveRDS(
  list(config = config, per_fit = per_fit, fit_status = fit_status, arm_summary = arm_summary, paired = paired, mech = mech),
  here::here("data-outputs/residual-seesaw.rds")
)

# Figure 1: raw signed seesaw errors, paired across arms; the mean and its 95%
# Monte Carlo interval (black) are separate from the replicate IQR (grey bars)
ratio_breaks <- c(0.5, 0.7, 1, 1.4, 2)
arm_labels <- c(
  A = "A. Alt.,\nomega est.,\nepsilon est.", B = "B. Alt.,\nomega est.,\nno epsilon",
  C = "C. Alt.,\nomega fixed", D = "D. Full data,\nomega est.", E = "E. Full data,\nomega fixed",
  F = "F. Alt. +\n10% overlap", G = "G. Alt. +\n25% overlap", H = "H. Area-\nproportional", I = "I. Alt.,\nreversed",
  J = "J. Const.,\nalt.", K = "K. Const.,\nreversed", L = "L. Const.,\n25% overlap", M = "M. Const.,\nalt., omega\nfixed", N = "N. Const.,\nfull"
)
pf |>
  mutate(arm = factor(arm, levels = names(arm_labels), labels = arm_labels)) |>
  ggplot(aes(arm, seesaw)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey60") +
  geom_line(aes(group = rep), colour = "grey80") +
  geom_point(colour = "grey60", alpha = 0.5) +
  stat_summary(fun.data = \(x) data.frame(y = mean(x), ymin = mean(x) - 1.96 * sd(x) / sqrt(length(x)), ymax = mean(x) + 1.96 * sd(x) / sqrt(length(x)))) +
  facet_wrap(~species, ncol = 1) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  scale_y_continuous(breaks = log(ratio_breaks), labels = paste0("×", ratio_breaks)) +
  labs(
    x = NULL, y = "Seesaw error (ratio of index error, HS vs. QCS occasions)",
    caption = "Flattened omega, epsilon exactly 0. Grey: paired replicates; black: mean and 95% Monte Carlo interval."
  )
ggsave(here::here("figs/residual-seesaw.pdf"), width = 14, height = 4.5 * n_distinct(pf$species), device = cairo_pdf)

# Figure 2: omega-contrast error, year-effect phase error, unsampled-region
# biomass error, and the seesaw (alternating arms, estimated omega)
pf |>
  filter(arm %in% c("A", "B", "G", "I")) |>
  select(species, arm, rep, seesaw, omega_err_grid, ye_phase_err, unsampled_contrast) |>
  tidyr::pivot_longer(c(ye_phase_err, unsampled_contrast, seesaw), names_to = "y", values_to = "value") |>
  mutate(y = factor(y,
    levels = c("ye_phase_err", "unsampled_contrast", "seesaw"),
    labels = c("Year-effect phase error (HS - QCS, sum of components)", "Unsampled-region log error,\nQCS unsampled minus HS unsampled", "Seesaw error (log)")
  )) |>
  ggplot(aes(omega_err_grid, value, colour = arm)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey60") +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey60") +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  facet_grid(species ~ y, scales = "free_y") +
  labs(x = "Estimated omega contrast error (HS - QCS, grid-weighted, sum of components)", y = NULL, colour = "Arm")
ggsave(here::here("figs/residual-seesaw-mechanism.pdf"), width = 11, height = 4, device = cairo_pdf)
