# Fit-level summaries to relate to seesaw metrics across species (see
# analysis/18-seesaw-correlates.R): random-effect parameters, the persistent
# spatial structure (omega and the density it implies) by region, a north-south
# gradient, and basic data summaries. Computed from one fitted model, so they
# can be cached next to the indexes.

#' @param fit A fitted sdmTMB model (the IID RF, factor(year) model).
#' @param grid Prediction grid with `area`, `X`, `Y`, the `year` column the
#'   model was fit with, and `grid_group` (if not NULL). Only the first year
#'   is predicted: omega and the year-free terms are the same in every year.
#' @param dat The data the model was fit to.
#' @param grid_group,dat_group Column names defining regions (e.g. survey) in
#'   `grid` and `dat`, or NULL for none.
#' @param response Name of the response column in `dat` (for prevalence).
#' @return list(ran_pars, regions, gradient, data): tibbles.
fit_covariates <- function(fit, grid, dat, response, grid_group = NULL, dat_group = NULL) {
  tidy_ok <- \(m) tryCatch(
    dplyr::mutate(tidy(fit, "ran_pars", model = m), model = m),
    error = \(e) NULL
  )
  n_comp <- length(fit$family$family)
  ran_pars <- dplyr::bind_rows(lapply(seq_len(n_comp), tidy_ok))

  g <- grid[grid$year == min(grid$year), , drop = FALSE]
  p <- predict(fit, newdata = g, return_tmb_object = FALSE)
  # Persistent log density: everything but the spatiotemporal fields. With a
  # delta poisson-link model, log density is the sum of the two linear
  # predictors; otherwise the single one. Without a year effect in the group
  # contrast it cancels, so only omega and the covariates contribute.
  comps <- if (n_comp > 1L) 1:2 else ""
  omega <- Reduce(`+`, lapply(comps, \(i) p[[paste0("omega_s", i)]]))
  non_rf <- Reduce(`+`, lapply(comps, \(i) p[[paste0("est_non_rf", i)]]))
  p$omega <- omega
  p$log_dens <- omega + non_rf
  for (i in seq_len(n_comp)) p[[paste0("omega_c", i)]] <- p[[paste0("omega_s", comps[i])]]
  p$region <- if (is.null(grid_group)) "all" else as.character(g[[grid_group]])

  regions <- p |>
    dplyr::summarise(
      omega = weighted.mean(omega, area),
      mean_log_density = weighted.mean(log_dens, area),
      log_mean_density = log(weighted.mean(exp(log_dens), area)),
      dplyr::across(dplyr::starts_with("omega_c"), \(x) weighted.mean(x, area)),
      area = sum(area),
      .by = region
    )

  # Area-weighted slope of persistent log density along the northing (per 100 km)
  w <- p$area / sum(p$area)
  gradient <- tibble::tibble(
    northing_slope = unname(coef(lm(p$log_dens ~ p$Y, weights = w))[2]) * 100,
    omega_northing_slope = unname(coef(lm(p$omega ~ p$Y, weights = w))[2]) * 100
  )

  dat$.region <- if (is.null(dat_group)) "all" else as.character(dat[[dat_group]])
  data <- dat |>
    dplyr::summarise(
      n_tows = dplyr::n(),
      prevalence = mean(.data[[response]] > 0),
      mean_depth = mean(depth_m),
      .by = .region
    ) |>
    dplyr::rename(region = .region)

  list(ran_pars = ran_pars, regions = regions, gradient = gradient, data = data)
}

# fit_covariates() from the operating models cached by a self-test (their
# setups hold fit_om), for runs made before the summaries were stored
cached_fit_covariates <- function(dir, species, response = "catch_weight", grid_group = NULL, dat_group = NULL) {
  slug <- \(x) gsub("^-|-$", "", gsub("[^a-z0-9]+", "-", tolower(x)))
  out <- lapply(setNames(species, species), \(sp) {
    f <- file.path(dir, "setups", paste0(slug(sp), ".rds"))
    s <- if (file.exists(f)) readRDS(f)
    if (is.null(s) || is.null(s$fit_om)) return(NULL)
    tryCatch(
      fit_covariates(s$fit_om, s$grid, s$dat, response, grid_group, dat_group),
      error = \(e) NULL
    )
  })
  purrr::compact(out)
}
