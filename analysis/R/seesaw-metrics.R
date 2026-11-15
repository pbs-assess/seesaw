library(dplyr, warn.conflicts = FALSE)

# Phase is +0.5 in the occasions one region is sampled and -0.5 in the other
# (each fit script sets which), so beta is the full difference between them in
# the detrended log index, and A = 100 * (exp(|beta|) - 1) is the biennial
# amplitude in percent.

period2_metric <- function(x, year, phase, conf = 0.95) {
  stopifnot(
    length(x) == length(year),
    length(phase) == length(x),
    all(x > 0),
    all(phase %in% c(-0.5, 0.5))
  )
  y <- log(x)
  fit <- lm(y ~ year + phase)

  beta <- coef(fit)[["phase"]]
  beta_se <- coef(summary(fit))["phase", "Std. Error"]
  ci <- confint(fit, "phase", level = conf)
  A_fun <- \(b) 100 * (exp(abs(b)) - 1)
  A_ci <- A_fun(ci)
  A_bounds <- if (ci[1L] <= 0 && ci[2L] >= 0) c(0, max(A_ci)) else range(A_ci)

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

# Windows span a fixed number of steps of `year`, so a window with a missing
# step has fewer observations. `year` is calendar year in the case studies and
# survey occasion in the self-tests.
moving_window <- function(x, year, phase, window) {
  starts <- seq(min(year), max(year) - window + 1L)
  out <- lapply(starts, function(y0) {
    ind <- which(year >= y0 & year < y0 + window)
    res <- period2_metric(x = x[ind], year = year[ind], phase = phase[ind])
    data.frame(start = y0, end = y0 + window - 1L, n = length(ind), t(res))
  })
  do.call(rbind, out)
}

mean_window_beta <- function(x, year, phase, window) {
  mean(moving_window(x, year = year, phase = phase, window = window)$beta)
}

seesaw_abs_error <- function(est, true, year, phase, window) {
  b_est <- mean_window_beta(est, year, phase, window)
  b_true <- mean_window_beta(true, year, phase, window)
  100 * (exp(abs(b_est - b_true)) - 1)
}

# Positive: the index reads too high in phase +0.5 occasions relative to
# phase -0.5 ones
seesaw_signed_error <- function(est, true, year, phase, window) {
  b_est <- mean_window_beta(est, year, phase, window)
  b_true <- mean_window_beta(true, year, phase, window)
  100 * (exp(b_est - b_true) - 1)
}

# Above 1, the index is wigglier than the truth (differencing weights the
# period-2 seesaw most); below 1, it is oversmoothed
year_to_year_sd_ratio <- function(est, true) {
  sd(diff(log(est))) / sd(diff(log(true)))
}

# Self-tests -----------------------------------------------------------------

# Compare fits on the same simulated populations. Filter `fit_status` to the
# fits a figure uses first. Replicates are dropped per species, so one species'
# failures cost the others nothing. A species with fewer than `min_complete` of
# its replicates complete is dropped entirely: its few survivors are likely the
# easy draws and would misrepresent it.
common_replicates <- function(indexes, fit_status, min_complete = 0.5) {
  by_rep <- fit_status |>
    summarise(complete = all(status == "ok"), .by = c(species, rep))
  by_species <- summarise(by_rep, frac = mean(complete), .by = species)
  dropped_species <- filter(by_species, frac < min_complete)
  if (nrow(dropped_species) > 0L) {
    cli::cli_inform(c(
      "Dropped {nrow(dropped_species)} species with < {100 * min_complete}% complete replicates:",
      set_names(sprintf("%s (%.0f%%)", dropped_species$species, 100 * dropped_species$frac), "*")
    ))
  }
  by_rep <- filter(by_rep, !species %in% dropped_species$species)
  keep <- filter(by_rep, complete)
  n_dropped <- sum(!by_rep$complete)
  if (n_dropped > 0L) cli::cli_inform("Dropped {n_dropped} species x replicate(s) with a failed fit")
  semi_join(indexes, keep, by = c("species", "rep"))
}

# rmse_ratio is relative to the species' mean RMSE of the full-data fit, so
# species are comparable
self_test_rep_metrics <- function(indexes, window) {
  q50 <- qnorm(0.75)
  indexes |>
    arrange(species, rep, design, model, year) |>
    summarise(
      rmse = sqrt(mean((log_est - log(true_est))^2)),
      coverage_50 = mean(abs(log_est - log(true_est)) <= q50 * se),
      seesaw_abs_error = seesaw_abs_error(est, true_est, year, phase, window),
      variability = year_to_year_sd_ratio(est, true_est),
      .by = c(species, rep, design, model)
    ) |>
    mutate(rmse_ratio = rmse / mean(rmse[design == "full"]), .by = species)
}

# The seesaw error of a model that estimates no seesaw at all: the reference
# for one that smooths the seesaw away
no_seesaw_error <- function(indexes, window) {
  indexes |>
    distinct(species, rep, year, true_est, phase) |>
    arrange(species, rep, year) |>
    summarise(beta_true = mean_window_beta(true_est, year, phase, window), .by = c(species, rep)) |>
    summarise(value = 100 * (exp(median(abs(beta_true))) - 1)) |>
    pull(value)
}

# Gradient simulations ---------------------------------------------------------

# Only seeds where every model succeeded
gradient_grid_metrics <- function(grid_out, window) {
  indexes <- grid_out$indexes |>
    semi_join(
      filter(grid_out$fit_status, all(status == "ok"), .by = c(delta, sigma_O, seed)),
      by = c("delta", "sigma_O", "seed")
    )
  by_seed <- indexes |>
    arrange(delta, sigma_O, seed, model, year) |>
    summarise(
      seesaw_error = seesaw_signed_error(est, true_est, year, phase, window),
      rmse = sqrt(mean((log_est - log(true_est))^2)),
      .by = c(delta, sigma_O, seed, model)
    )
  by_cell <- by_seed |>
    summarise(
      n = n(),
      seesaw_med = median(seesaw_error),
      seesaw_lo = quantile(seesaw_error, 0.25),
      seesaw_hi = quantile(seesaw_error, 0.75),
      rmse = mean(rmse),
      .by = c(model, delta, sigma_O)
    )
  list(by_seed = as_tibble(by_seed), by_cell = as_tibble(by_cell))
}
