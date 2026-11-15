library(dplyr, warn.conflicts = FALSE)
library(tidyr)
library(purrr)
library(sdmTMB)

# Parallel workers need these attached too
analysis_packages <- c("dplyr", "tidyr", "purrr", "sdmTMB")
furrr_opts <- function(...) furrr::furrr_options(seed = TRUE, packages = analysis_packages, ...)

# year = "none": any year-to-year change is carried by the spatiotemporal field
paper_models <- list(
  "IID RF, factor(year)" = list(spatiotemporal = "iid", year = "factor"),
  "IID RF, RW year"      = list(spatiotemporal = "iid", year = "rw"),
  "RW RF, RW year"       = list(spatiotemporal = "rw",  year = "rw"),
  "AR1 RF, RW year"      = list(spatiotemporal = "ar1", year = "rw"),
  "RW RF"                = list(spatiotemporal = "rw",  year = "none")
)

# Lognormal prior on sigma_V with its mode at target_sd. The lognormal mode is
# exp(meanlog - sdlog^2), so with a fixed meanlog a wider prior would also pull
# sigma_V towards zero. Shifting meanlog by sdlog^2 makes this a normal penalty
# on log(sigma_V) centred at log(target_sd) (the scale sdmTMB optimizes on):
# target_sd sets where the estimate is pulled and sdlog only how strongly.
sigma_V_meanlog <- function(target_sd, sdlog) log(target_sd) + sdlog^2
sigma_V_prior <- function(target_sd, sdlog) {
  lognormal_prior(sigma_V_meanlog(target_sd, sdlog), sdlog)
}

# Fields whose SD collapses to zero (e.g., the encounter spatiotemporal field
# with alternating data) would otherwise fail to converge
default_control <- sdmTMBcontrol(collapse_spatial_variance = TRUE)
default_sigma_V_prior <- sigma_V_prior(0.3, 0.5)

# offset is already logged; covariates are extra RHS terms, e.g. "factor(gear)"
fit_model <- function(model, data, mesh, response, family, offset = NULL,
                      covariates = NULL, anisotropy = TRUE, control = default_control,
                      sigma_V_prior = default_sigma_V_prior) {
  spec <- paper_models[[model]]
  if (is.null(spec)) cli::cli_abort("Unknown model: {model}")
  rw_year <- spec$year == "rw"
  year_term <- if (spec$year == "factor") "0 + factor(year)" else "1"
  sdmTMB(
    reformulate(c(year_term, covariates), response),
    data = data, mesh = mesh, offset = offset, family = family,
    time = "year", spatial = "on", spatiotemporal = spec$spatiotemporal,
    share_range = TRUE, anisotropy = anisotropy,
    time_varying = if (rw_year) ~1,
    time_varying_type = "rw0",
    priors = if (rw_year && !is.null(sigma_V_prior)) sdmTMBpriors(sigma_V = sigma_V_prior) else sdmTMBpriors(),
    # models that carry information between years can predict missing years
    extra_time = if (spec$year != "factor") seq(min(data$year), max(data$year)),
    control = control, silent = TRUE
  )
}

fit_ok <- function(expr) {
  fit <- tryCatch(expr, error = \(e) NULL)
  ok <- inherits(fit, "sdmTMB") && isTRUE(tryCatch(
    all(unlist(sanity(fit, gradient_thresh = 0.01, silent = TRUE))),
    error = \(e) FALSE
  ))
  if (ok) fit else NULL
}

# nsplit > 1 predicts in chunks of years to reduce peak memory on large grids
get_index_ok <- function(fit, grid, nsplit = 2L) {
  tryCatch(
    get_index_split(
      fit,
      newdata = grid, area = grid$area, offset = rep(0, nrow(grid)), nsplit = nsplit,
      bias_correct = TRUE, silent = TRUE
    ),
    error = \(e) NULL
  )
}

# With collapse_spatial_variance, sdmTMB refits with a near-zero
# spatiotemporal field turned off. In "RW RF" (no year term) a collapse in
# every component leaves no time variation: a flat index that would score as
# no seesaw, so it counts as a failed fit. In the other models st_collapsed
# flags that the IID, RW, and AR1 RF variants have become the same model.
index_ok <- function(fit, model, grid, nsplit = 2L) {
  if (is.null(fit)) return(NULL)
  collapsed <- all(fit$spatiotemporal == "off")
  if (collapsed && paper_models[[model]]$year == "none") return(NULL)
  index <- get_index_ok(fit, grid, nsplit = nsplit)
  if (!is.null(index)) index$st_collapsed <- collapsed
  index
}

fit_index_models <- function(data, grid, mesh, models = names(paper_models), nsplit = 2L, ...) {
  map(set_names(models), \(m) index_ok(fit_ok(fit_model(m, data = data, mesh = mesh, ...)), m, grid, nsplit)) |>
    discard(is.null) |>
    bind_rows(.id = "model") |>
    # typed, so callers' mutate() and joins still work if every model fails
    bind_rows(tibble(model = character(), year = integer()))
}

# Skipped if `file` exists (delete it to refit). fit_species() returns NULL to
# leave a species out (e.g., too few positive sets). A species that errors is
# reported with status "error" rather than stopping the run and losing the
# other species' fits. `designs` names the designs fit_species() returns, if
# it returns a design column, so a design that fails for every species is
# still counted.
#
# Each species is also cached in a folder named after `file`, so if the run
# dies (e.g., a worker runs out of memory) a rerun fits only the species not
# yet done. Errors aren't cached, so a rerun retries them.
run_case_study <- function(species, fit_species, file, n_workers, designs = NULL) {
  if (file.exists(file)) return(invisible(NULL))
  species_dir <- sub("\\.rds$", "", file)
  dir.create(species_dir, showWarnings = FALSE, recursive = TRUE)
  species_file <- \(s) file.path(species_dir, paste0(slug(s), ".rds"))
  stopifnot(!anyDuplicated(species_file(species)))
  single_thread()
  future::plan(future::multisession, workers = n_workers)
  on.exit(future::plan(future::sequential), add = TRUE)
  res <- furrr::future_map(set_names(species), \(s) {
    f <- species_file(s)
    if (file.exists(f)) return(readRDS(f))
    out <- tryCatch(fit_species(s), error = \(e) structure(conditionMessage(e), class = "species_error"))
    if (!inherits(out, "species_error")) saveRDS(out, f)
    out
  }, .options = furrr_opts())
  errored <- keep(res, \(x) inherits(x, "species_error"))
  if (length(errored) > 0L) {
    cli::cli_warn("Fitting errored for {length(errored)} species")
    iwalk(errored, \(e, s) message(s, ": ", e))
  }
  res <- discard(res, \(x) is.null(x) || inherits(x, "species_error"))
  indexes <- bind_rows(res, .id = "species") |>
    bind_rows(tibble(species = character(), model = character(), year = integer())) |>
    select(any_of(c(
      "species", "design", "model", "family", "year", "est", "lwr", "upr", "log_est", "se", "phase", "st_collapsed"
    )))
  if (is.null(designs) && "design" %in% names(indexes)) designs <- unique(indexes$design)
  fitted <- species[species %in% c(names(res), names(errored))]
  fit_status <- expand_grid(species = fitted, design = designs, model = names(paper_models))
  keys <- names(fit_status)
  fit_status <- fit_status |>
    left_join(mutate(distinct(indexes, across(all_of(keys))), status = "ok"), by = keys) |>
    mutate(status = case_when(
      !is.na(status) ~ status,
      species %in% names(errored) ~ "error",
      .default = "failed"
    ))
  failed <- filter(fit_status, status != "ok")
  if (nrow(failed) > 0L) {
    cli::cli_warn("{nrow(failed)} fit(s) failed sanity checks or errored")
    print(failed, n = Inf)
  }
  dir.create(dirname(file), showWarnings = FALSE, recursive = TRUE)
  saveRDS(list(indexes = indexes, fit_status = fit_status), file)
}

slug <- \(x) gsub("^-|-$", "", gsub("[^a-z0-9]+", "-", tolower(x)))

single_thread <- function() {
  RhpcBLASctl::blas_set_num_threads(1L)
  RhpcBLASctl::omp_set_num_threads(1L)
}
