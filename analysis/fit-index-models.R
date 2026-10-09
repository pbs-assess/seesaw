# Shared set of sdmTMB index models for comparing seesaw behaviour across surveys.
# Survey-specific scripts prepare `dat`, `grid`, and `mesh` and call
# `fit_index_models()`.

# Returns the fit if it ran and passed sanity checks, otherwise NULL
fit_ok <- function(expr) {
  fit <- tryCatch(expr, error = \(e) NULL)
  ok <- inherits(fit, "sdmTMB") && isTRUE(tryCatch(
    all(unlist(sdmTMB::sanity(fit, gradient_thresh = 0.01, silent = TRUE))),
    error = \(e) FALSE
  ))
  if (ok) fit else NULL
}

get_index_ok <- function(fit, grid) {
  tryCatch(
    sdmTMB::get_index(fit, newdata = grid, area = grid$area, offset = rep(0, nrow(grid))),
    error = \(e) NULL
  )
}

# Clamp grid depths to the sampled range to avoid extrapolating the depth polynomial
clamp_depth <- function(grid, dat, depth_col = "depth_m") {
  depth_range <- range(dat[[depth_col]])
  grid[[depth_col]] <- pmin(pmax(grid[[depth_col]], depth_range[1]), depth_range[2])
  grid
}

#' @param dat Data with columns `year`, `X`, `Y`, `depth_m`, and the response.
#' @param grid Prediction grid with columns `year`, `X`, `Y`, `depth_m`, `area`.
#' @param mesh Mesh built from `dat`.
#' @param response Name of the response column.
#' @param family sdmTMB family, or a list of families. With a list, the base
#'   model is fit with each and the family of the converged fit with the
#'   lowest AIC (or lowest AIC overall if none converge) is used for all models.
#' @param offset Offset vector (already logged) for `dat`.
#' @param covariates Optional character vector of extra RHS terms (e.g.
#'   `"factor(survey_group)"`) kept in every model. Columns must be in `dat`.
#' @param all_data Optional list(data, mesh, offset) for the "all data" models,
#'   which also need a `survey_name` column.
#' @param control [sdmTMB::sdmTMBcontrol()] list used for all models.
fit_index_models <- function(dat, grid, mesh, response, family, offset, covariates = NULL, all_data = NULL,
                             control = sdmTMB::sdmTMBcontrol()) {
  all_yrs <- seq(min(dat$year), max(dat$year))
  base_formula <- stats::as.formula(paste(response, "~", paste(c("0 + factor(year)", covariates), collapse = " + ")))
  # RHS for models where the year effect is not a fixed factor
  no_year_formula <- stats::as.formula(paste(". ~", paste(c("1", covariates), collapse = " + ")))
  depth_formula <- stats::as.formula(paste(
    ". ~", paste(c("0 + factor(year)", covariates, "poly(log(depth_m), 2)"), collapse = " + ")
  ))
  rw_prior <- sdmTMB::sdmTMBpriors(sigma_V = sdmTMB::gamma_cv(0.3, 0.5))

  # The base fit is the update() template for the other models, so keep it
  # even if it fails sanity checks; it is only reported if it passes.
  # bquote() embeds the formula and family in the call so update() can
  # modify the formula and re-evaluate it outside fit_base().
  fit_base <- function(family) {
    tryCatch(eval(bquote(sdmTMB::sdmTMB(
      .(base_formula),
      data = dat,
      mesh = mesh,
      offset = .(offset),
      family = .(family),
      time = "year",
      spatial = "on",
      spatiotemporal = "iid",
      share_range = TRUE,
      anisotropy = TRUE,
      control = .(control),
      silent = FALSE
    ))), error = \(e) NULL)
  }

  if (inherits(family, "family")) family <- list(family)
  candidates <- purrr::compact(lapply(family, fit_base))
  if (length(candidates) == 0L) {
    return(dplyr::tibble())
  }
  converged <- purrr::compact(lapply(candidates, fit_ok))
  pool <- if (length(converged) > 0L) converged else candidates
  template <- pool[[which.min(vapply(pool, stats::AIC, numeric(1)))]]

  fits <- list()
  fits[["IID RF, factor(year)"]] <- fit_ok(template)

  if (!is.null(all_data)) {
    fits[["IID RF, factor(year), all data"]] <- fit_ok(update(
      template,
      data = all_data$data,
      mesh = all_data$mesh,
      offset = all_data$offset
    ))

    fits[["IID RF, factor(year), factor(survey), all data"]] <- fit_ok(update(
      template,
      formula. = . ~ factor(year) + factor(survey_name),
      data = all_data$data,
      mesh = all_data$mesh,
      offset = all_data$offset
    ))
  }

  fits[["RW RF"]] <- fit_ok(update(
    template,
    formula. = no_year_formula,
    spatiotemporal = "rw",
    extra_time = all_yrs
  ))

  fits[["AR1 RF"]] <- fit_ok(update(
    template,
    formula. = no_year_formula,
    spatiotemporal = "ar1",
    extra_time = all_yrs
  ))

  fits[["RW RF, RW year"]] <- fit_ok(update(
    template,
    formula. = no_year_formula,
    spatiotemporal = "rw",
    time_varying = ~1,
    time_varying_type = "rw0",
    priors = rw_prior,
    extra_time = all_yrs
  ))

  fits[["AR1 RF, RW year"]] <- fit_ok(update(
    template,
    formula. = no_year_formula,
    spatiotemporal = "ar1",
    time_varying = ~1,
    time_varying_type = "rw0",
    priors = rw_prior,
    extra_time = all_yrs
  ))

  fits[["IID RF, RW year"]] <- fit_ok(update(
    template,
    formula. = no_year_formula,
    spatiotemporal = "iid",
    time_varying = ~1,
    time_varying_type = "rw0",
    priors = rw_prior,
    extra_time = all_yrs
  ))

  fits[["Spatial only, RW year"]] <- fit_ok(update(
    template,
    formula. = no_year_formula,
    spatiotemporal = "off",
    time_varying = ~1,
    time_varying_type = "rw0",
    priors = rw_prior,
    extra_time = all_yrs
  ))

  fits[["IID RF, factor(year), depth"]] <- fit_ok(update(
    template,
    formula. = depth_formula
  ))

  purrr::compact(fits) |>
    purrr::map(\(f) get_index_ok(f, grid)) |>
    dplyr::bind_rows(.id = "model") |>
    dplyr::mutate(family = paste(template$family$family, collapse = "/"))
}
