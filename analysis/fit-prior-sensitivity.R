# Sensitivity of IID RF, RW year to the sigma_V prior: Pacific cod under
# pretend biennial QCS / HS sampling, against IID RF, factor(year) fit to the
# full data. Time is survey occasion, so sigma_V is per occasion (2 calendar
# years). Fits that fail sanity checks are kept and flagged, not dropped: the
# paper cites the no-prior fit failing.

library(dplyr)
source(here::here("analysis/R/models.R"))
source(here::here("analysis/R/survey-data.R"))
load_surveyjoin()

out_file <- here::here("data-generated/prior-sensitivity.rds")
dir.create(dirname(out_file), showWarnings = FALSE)
# Delete out_file to refit
if (!file.exists(out_file)) {
  prior_grid <- crossing(target_sd = c(0.1, 0.3, 0.5), sdlog = c(0.1, 0.3, 0.5, 0.8)) |>
    mutate(prior = paste0("target = ", target_sd, ", sdlog = ", sdlog))

  dat <- qcs_hs_data("pacific cod")
  dat_biennial <- filter(dat, survey_name == primary_survey)
  grid <- replicate_df(qcs_hs_grid(), "year", sort(unique(dat$year)))
  mesh <- make_mesh(dat, c("X", "Y"), cutoff = 10)
  mesh_biennial <- make_mesh(dat_biennial, c("X", "Y"), mesh = mesh$mesh)

  fit_cod <- function(model, d, mesh, sigma_V_prior = NULL) {
    fit_model(
      model,
      data = d, mesh = mesh, response = "catch_weight",
      family = delta_gamma(type = "poisson-link"), offset = log(d$effort),
      sigma_V_prior = sigma_V_prior
    )
  }

  summarise_fit <- function(fit) {
    if (inherits(fit, "error")) {
      return(list(sanity_ok = FALSE, index = NULL, sigma_V = NULL))
    }
    # tidy() omits sigma_V
    est <- as.list(fit$sd_report, "Estimate", report = TRUE)$sigma_V
    sigma_V <- if (length(est) > 0L) {
      tibble(
        component = c("encounter", "positive"),
        estimate = as.numeric(est),
        std.error = as.numeric(as.list(fit$sd_report, "Std. Error", report = TRUE)$sigma_V)
      )
    }
    list(
      sanity_ok = isTRUE(all(unlist(sanity(fit, gradient_thresh = 0.01, silent = TRUE)))),
      index = get_index_ok(fit, grid),
      sigma_V = sigma_V
    )
  }

  jobs <- c(
    list(reference = "reference", "no prior" = NULL),
    set_names(map2(prior_grid$target_sd, prior_grid$sdlog, sigma_V_prior), prior_grid$prior)
  )

  single_thread()
  future::plan(future::multisession, workers = min(length(jobs), parallelly::availableCores() / 2))
  res <- furrr::future_imap(jobs, \(prior, name) {
    single_thread()
    fit <- tryCatch(
      if (name == "reference") {
        fit_cod("IID RF, factor(year)", dat, mesh)
      } else {
        fit_cod("IID RF, RW year", dat_biennial, mesh_biennial, sigma_V_prior = prior)
      },
      error = \(e) e
    )
    summarise_fit(fit)
  }, .options = furrr_opts(scheduling = Inf))
  future::plan(future::sequential)

  stopifnot(res$reference$sanity_ok)
  fit_info <- tibble(prior = names(res)[-1], sanity_ok = map_lgl(res[-1], "sanity_ok")) |>
    left_join(prior_grid, by = "prior")
  saveRDS(
    list(
      prior_grid = prior_grid,
      occasions = distinct(dat, year, cal_year, primary_survey),
      reference = res$reference$index,
      fit_info = fit_info,
      indexes = imap_dfr(res[-1], \(r, name) if (!is.null(r$index)) mutate(r$index, prior = name)),
      sigma_V = imap_dfr(res[-1], \(r, name) if (!is.null(r$sigma_V)) mutate(r$sigma_V, prior = name))
    ),
    out_file
  )
}
