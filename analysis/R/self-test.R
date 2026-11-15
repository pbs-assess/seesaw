# Simulation self-tests with a known truth:
#
#   1. Fit the operating model to the full survey data.
#   2. For each replicate, draw the spatial fields from their approximate
#      posterior and entirely new spatiotemporal fields and observation error,
#      at the observed tows and, on the same field draws, on the prediction
#      grid (the true index).
#   3. Keep the tows of each sampling design and fit the estimation models.
#
# Every random draw for a replicate is seeded by its replicate number and
# every stage is cached, so raising n_reps later only runs the new replicates.

source(here::here("analysis/R/models.R"))

self_test_family <- delta_gamma(type = "poisson-link")
om_model <- "IID RF, factor(year)"

fit_self_test_model <- function(model, d, mesh) {
  fit_model(
    model,
    data = d, mesh = mesh, response = "catch_weight",
    family = self_test_family, offset = log(d$effort)
  )
}

# On failure, keep diagnostics so the cause can be seen without refitting
fit_operating_model <- function(d, mesh) {
  fit <- tryCatch(fit_self_test_model(om_model, d, mesh), error = \(e) e)
  if (!inherits(fit, "sdmTMB")) {
    return(list(om_failed = TRUE, error = conditionMessage(fit)))
  }
  checks <- tryCatch(
    unlist(sanity(fit, gradient_thresh = 0.01, silent = TRUE)),
    error = \(e) c(sanity_error = FALSE)
  )
  if (all(checks)) {
    return(list(fit = fit))
  }
  list(
    om_failed = TRUE,
    sanity = checks,
    pd_hessian = fit$pos_def_hessian,
    convergence = fit$model$message
  )
}

# Sampling designs ---------------------------------------------------------------

# Separate seed ranges so the parameter, simulation, and design random streams
# don't overlap; each is offset by the replicate number
seed_pars <- 1e6L
seed_sim <- 2e6L
seed_design <- 3e6L

# Same number of tows per occasion as the alternating design, spread over both
# regions
sparse_keep <- function(setup, rep) {
  set.seed(seed_design + rep)
  keep <- logical(nrow(setup$dat))
  for (y in unique(setup$dat$year)) {
    in_year <- which(setup$dat$year == y)
    n_alt <- sum(setup$is_sampled[in_year])
    keep[in_year[sample.int(length(in_year), n_alt)]] <- TRUE
  }
  keep
}

# Adds the proportion p of the off region's tows closest to the sampled region
# and drops as many at random from the sampled region, so effort matches the
# alternating design
overlap_keep <- function(setup, rep, p) {
  dat <- setup$dat
  grid_xy <- distinct(setup$grid, X, Y, survey)
  dist_to_other <- numeric(nrow(dat))
  for (s in unique(dat$survey_name)) {
    i <- dat$survey_name == s
    other <- as.matrix(grid_xy[grid_xy$survey != s, c("X", "Y")])
    dist_to_other[i] <- FNN::get.knnx(other, as.matrix(dat[i, c("X", "Y")]), k = 1L)$nn.dist[, 1L]
  }
  set.seed(seed_design + rep)
  keep <- setup$is_sampled
  for (y in unique(dat$year)) {
    primary <- which(dat$year == y & setup$is_sampled)
    off <- which(dat$year == y & !setup$is_sampled)
    k <- round(p * length(off))
    keep[off[order(dist_to_other[off])[seq_len(k)]]] <- TRUE
    keep[primary[sample.int(length(primary), k)]] <- FALSE
  }
  keep
}

# design: "full", "alternating", "sparse", or "overlapNN" (NN = percent of the
# off region's tows). Returns an empty tibble if the fit fails (see
# index_ok()).
fit_design <- function(setup, catch, design, model, rep, index_nsplit = 2L) {
  d <- setup$dat
  d$catch_weight <- catch
  mesh <- setup$mesh
  if (design == "alternating") {
    d <- d[setup$is_sampled, ]
    mesh <- setup$mesh_alt
  } else if (design != "full") {
    keep <- if (design == "sparse") {
      sparse_keep(setup, rep)
    } else {
      overlap_keep(setup, rep, as.numeric(sub("overlap", "", design)) / 100)
    }
    d <- d[keep, ]
    # same knots as the other designs
    mesh <- make_mesh(d, c("X", "Y"), mesh = setup$mesh$mesh)
  }
  index <- timed_index(fit_self_test_model(model, d, mesh), model, setup$grid, nsplit = index_nsplit)
  if (is.null(index)) tibble() else index
}

# Simulation -----------------------------------------------------------------------

# Tows and grid cells go in one newdata so both share the same field draws.
# eta_i excludes the offset, so grid rows give expected density per unit area.
sim_setup <- function(species, dat, is_sampled, grid, mesh, mesh_alt, occasions) {
  cols <- c("X", "Y", "year")
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

# Fixed effects at their MLEs and random effects drawn from their approximate
# posterior (as mle_mvn_samples = "multiple"), then the spatiotemporal fields
# replaced by new draws from the GMRF prior. `sampler` is
# sdmTMB:::.posterior_re_sampler(fit), built once per fit because it's slow.
simulate_rep <- function(setup, sampler, rep) {
  set.seed(seed_pars + rep)
  pars <- sampler(1L)
  report <- simulate(
    setup$fit_om,
    nsim = 1L,
    mcmc_samples = pars,
    simulate_re = "spatiotemporal",
    newdata = setup$sim_newdata,
    offset = log(setup$sim_newdata$effort),
    return_tmb_report = TRUE,
    seed = seed_sim + rep,
    silent = TRUE
  )[[1L]]
  tows <- setup$is_tow
  grid_rows <- !tows
  truth <- setup$grid |>
    mutate(density = exp(report$eta_i[grid_rows, 1L] + report$eta_i[grid_rows, 2L]) * area) |>
    summarise(true_est = sum(density), .by = year)
  list(
    catch = report$y_i[tows, 1L] * report$y_i[tows, 2L],
    truth = truth
  )
}

# Pipeline -------------------------------------------------------------------------

minutes_since <- \(t0) round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)

log_progress <- function(log_file, ...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S "), ..., "\n", sep = "", file = log_file, append = TRUE)
}

# expr returns a note for the log or NULL. Errors are logged, not raised, so
# the other tasks continue.
run_logged <- function(log_file, label, expr) {
  t0 <- Sys.time()
  tryCatch(
    {
      note <- expr
      log_progress(log_file, "done ", label, note, " [", minutes_since(t0), " min]")
    },
    error = \(e) log_progress(log_file, "ERROR ", label, ": ", conditionMessage(e))
  )
  invisible(NULL)
}

# The workers are top-level functions that take everything as arguments: a
# closure sent to a worker carries its enclosing environment, which would ship
# every species' setup to every task.
setup_worker <- function(species, dat, file, setup_species, log_file) {
  single_thread()
  run_logged(log_file, paste0("setup: ", species), {
    setup <- setup_species(species, dat)
    if (!is.null(setup$fit_om)) setup$fit_om <- strip_sdmTMB(setup$fit_om)
    saveRDS(setup, file)
    if (isTRUE(setup$om_failed)) {
      why <- if (!is.null(setup$error)) setup$error else paste(names(setup$sanity)[!setup$sanity], collapse = ", ")
      paste0(" (OM FAILED: ", why, ")")
    }
  })
}

# Workers read the setup (with its large operating-model fit) from disk rather
# than being sent it: serializing a copy per task made the main session the
# bottleneck and the simulations ran one at a time. Each task runs a chunk of
# one species' replicates so the setup is read and the sampler built once per
# chunk; a cache kept on the worker wouldn't help, as each future gets a fresh
# copy of its globals.
sim_worker <- function(species, setup_file, reps, files, log_file) {
  single_thread()
  loaded <- tryCatch(
    {
      setup <- readRDS(setup_file)
      setup$fit_om <- sdmTMB:::reload_model(setup$fit_om)
      list(setup = setup, sampler = sdmTMB:::.posterior_re_sampler(setup$fit_om))
    },
    error = \(e) {
      log_progress(log_file, "ERROR sim setup: ", species, ": ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(loaded)) return(invisible(NULL))
  for (i in seq_along(reps)) {
    run_logged(log_file, paste0("sim: ", species, " rep ", reps[i]), {
      saveRDS(simulate_rep(loaded$setup, loaded$sampler, reps[i]), files[i])
      NULL
    })
  }
}

fit_worker <- function(setup, catch, rep, design, model, file, index_nsplit, log_file) {
  single_thread()
  run_logged(log_file, paste0("fit: ", setup$species, " rep ", rep, " ", design, " ", model), {
    index <- fit_design(setup, catch, design, model, rep, index_nsplit)
    saveRDS(index, file)
    if (nrow(index) == 0L) " FAILED"
  })
}

# Three cached, parallel stages (setup per species, simulation per replicate,
# fit per replicate x row of `fits`); only what's missing is run. Errors are
# logged and not cached, so a rerun retries them; a failed operating model or
# sanity check is cached (delete its file to redo it). dat is prepared in the
# main session because surveyjoin's database isn't safe to load from parallel
# workers.
#
# setup_species(species, dat) returns a sim_setup() list plus fit_om, or the
# fit_operating_model() failure. `fits` has columns design and model.
#
# Returns list(indexes, fit_status); fit_status is "ok", "failed" (sanity
# checks, index, or a collapsed "RW RF"; see index_ok()), "fit error", "sim
# error", or "om failed" for every requested fit.
run_self_test <- function(dir, species, n_reps, prepare_data, setup_species, fits,
                          n_workers = 1L, index_nsplit = 2L) {
  for (d in file.path(dir, c("setups", "sims", "fits"))) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  setup_file <- \(s) file.path(dir, "setups", paste0(slug(s), ".rds"))
  sim_file <- \(s, rep) file.path(dir, "sims", sprintf("%s-rep%03d.rds", slug(s), rep))
  fit_file <- \(s, rep, design, model) {
    file.path(dir, "fits", sprintf("%s-rep%03d-%s-%s.rds", slug(s), rep, design, slug(model)))
  }
  log_file <- file.path(dir, "progress.log")

  future::plan(future::multisession, workers = n_workers)
  on.exit(future::plan(future::sequential), add = TRUE)
  opts <- furrr_opts(scheduling = Inf)
  log_progress(log_file, "=== run started: ", length(species), " species, ", n_reps, " reps, ", n_workers, " workers ===")

  # Setup
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
  setups <- lapply(set_names(species), \(s) if (file.exists(setup_file(s))) readRDS(setup_file(s)))
  failed <- names(keep(setups, \(x) is.null(x) || isTRUE(x$om_failed)))
  if (length(failed) > 0L) cli::cli_warn("No operating model for: {failed}; see {log_file}")
  setups <- setups[setdiff(species, failed)]

  # Simulation
  sims <- expand_grid(species = names(setups), rep = seq_len(n_reps)) |>
    mutate(file = sim_file(species, rep))
  todo <- filter(sims, !file.exists(file))
  if (nrow(todo) > 0L) {
    log_progress(log_file, nrow(todo), " simulations to run")
    # enough chunks to keep every worker busy
    n_chunks <- max(1L, ceiling(n_workers / n_distinct(todo$species)))
    chunks <- todo |>
      mutate(chunk = ceiling(row_number() * n_chunks / n()), .by = species) |>
      summarise(reps = list(rep), files = list(file), .by = c(species, chunk))
    furrr::future_pwalk(
      list(species = chunks$species, setup_file = setup_file(chunks$species), reps = chunks$reps, files = chunks$files),
      sim_worker,
      log_file = log_file,
      .options = opts
    )
  }
  sims <- filter(sims, file.exists(file))
  sim_out <- lapply(sims$file, readRDS)

  # Fits; each task gets only its species' data, grid, and meshes plus its
  # replicate's catches
  fit_tasks <- expand_grid(i_sim = seq_along(sim_out), fits) |>
    mutate(
      species = sims$species[i_sim],
      rep = sims$rep[i_sim],
      file = fit_file(species, rep, design, model)
    )
  todo <- filter(fit_tasks, !file.exists(file))
  if (nrow(todo) > 0L) {
    log_progress(log_file, nrow(todo), " of ", nrow(fit_tasks), " fits to run")
    setups_lite <- lapply(setups, \(x) x[c("species", "dat", "is_sampled", "grid", "mesh", "mesh_alt")])
    furrr::future_pwalk(
      list(
        setup = setups_lite[todo$species],
        catch = lapply(todo$i_sim, \(i) sim_out[[i]]$catch),
        rep = todo$rep, design = todo$design, model = todo$model, file = todo$file
      ),
      fit_worker,
      index_nsplit = index_nsplit, log_file = log_file,
      .options = opts
    )
  }
  log_progress(log_file, "=== run finished ===")

  done <- filter(fit_tasks, file.exists(file))
  done_indexes <- lapply(done$file, readRDS)
  done$status <- if_else(map_int(done_indexes, nrow) > 0L, "ok", "failed")
  fit_status <- expand_grid(species = species, rep = seq_len(n_reps), fits) |>
    left_join(select(fit_tasks, species, rep, design, model, file), by = c("species", "rep", "design", "model")) |>
    left_join(select(done, file, done_status = status), by = "file") |>
    mutate(status = case_when(
      species %in% failed ~ "om failed",
      is.na(file) ~ "sim error",
      !is.na(done_status) ~ done_status,
      .default = "fit error"
    )) |>
    select(species, rep, design, model, status)
  n_bad <- sum(fit_status$status != "ok")
  if (n_bad > 0L) {
    cli::cli_warn("{n_bad} of {nrow(fit_status)} requested fit(s) failed; see fit_status and {log_file}")
    print(count(filter(fit_status, status != "ok"), species, design, model, status), n = Inf)
  }

  ok <- done$status == "ok"
  truth <- imap_dfr(sim_out, \(sim, i) mutate(sim$truth, species = sims$species[i], rep = sims$rep[i]))
  occasions <- map_dfr(setups, \(x) select(x$occasions, year, cal_year, phase), .id = "species")
  indexes <- pmap_dfr(
    list(done_indexes[ok], done$species[ok], done$rep[ok], done$design[ok], done$model[ok]),
    \(index, s, rep, design, model) mutate(index, species = s, rep = rep, design = design, model = model, .before = 1L)
  ) |>
    select(species, rep, design, model, year, est, lwr, upr, log_est, se, any_of(c("st_collapsed", "fit_time", "index_time"))) |>
    left_join(truth, by = c("species", "rep", "year")) |>
    left_join(occasions, by = c("species", "year"))

  list(indexes = indexes, fit_status = fit_status)
}
