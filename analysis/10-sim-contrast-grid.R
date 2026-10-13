# Simulation grid: is the seesaw governed by the regional contrast in the
# spatial field relative to its SD (Delta / sigma_O)?

library(dplyr)
library(ggplot2)
library(sdmTMB)

source(here::here("analysis/funcs.R"))
source(here::here("analysis/metric-functions.R"))
source(here::here("analysis/fit-index-models.R"))

theme_set(ggsidekick::theme_sleek())

deltas <- c(0.5, 1, 2)
sigma_Os <- c(1, 2, 4)
seeds <- seq_len(50L)

# Held fixed
n_year <- 12L
matern_range <- 0.1
sigma_E <- 1
phi <- 10
year_arima.sim <- list(ar = 0.5)
year_marginal_sd <- 0.3
year_mean <- 1
obs_n <- 200L
region_cutoff <- 0.5
gap <- 0 # contiguous regions, as QCS / HS
north_yrs <- seq(1, n_year - 1, 2)
south_yrs <- seq(2, n_year, 2)

models <- c("IID RF, factor(year)", "IID RF, RW year")

grid <- expand.grid(
  X = seq(0, 1, length.out = 100), Y = seq(0, 1, length.out = 100),
  year = seq_len(n_year)
)
is_north <- grid$Y >= region_cutoff
mean_Y_diff <- mean(grid$Y[is_north]) - mean(grid$Y[!is_north])

# Simulate ---------------------------------------------------------------------

simulate_cell <- function(delta, sigma_O, seed, mesh) {
  sim_one <- function(slope) {
    sim(
      grid, mesh,
      seed = seed, phi = phi, range = matern_range, region_cutoff = region_cutoff,
      sigma_O = sigma_O, sigma_E = sigma_E, coefs = c(slope, 0),
      # centre the gradient so it doesn't change the mean density
      year_mean = year_mean - slope * mean(grid$Y),
      year_arima.sim = year_arima.sim, year_marginal_sd = year_marginal_sd
    )$sim_dat
  }
  pilot <- sim_one(0)
  y1 <- pilot$year == 1L
  omega_diff <- mean(pilot$omega_s[y1 & is_north]) - mean(pilot$omega_s[y1 & !is_north])
  slope <- (delta - omega_diff) / mean_Y_diff
  s <- sim_one(slope)
  stopifnot(isTRUE(all.equal(s$omega_s, pilot$omega_s)))

  spatial <- s$omega_s + slope * s$Y
  realized <- mean(spatial[y1 & is_north]) - mean(spatial[y1 & !is_north])
  stopifnot(abs(realized - delta) < 1e-6)
  s
}

# Fit --------------------------------------------------------------------------

fit_model <- function(model, d, mesh) {
  args <- list(
    data = d, mesh = mesh, family = tweedie(), time = "year", spatial = "on",
    spatiotemporal = "iid", silent = TRUE,
    control = sdmTMBcontrol(nlminb_loops = 1L, newton_loops = 1L)
  )
  extra <- switch(model,
    "IID RF, factor(year)" = list(formula = observed ~ 0 + factor(year)),
    "IID RF, RW year" = list(
      formula = observed ~ 0, time_varying = ~1, time_varying_type = "rw",
      priors = sdmTMBpriors(sigma_V = lognormal_prior(log(0.3), 0.5))
    )
  )
  fit_ok(do.call(sdmTMB, c(args, extra)))
}

# Index plus the estimated contrast in omega and sigma_O (the diagnostic)
index_and_contrast <- function(fit) {
  if (is.null(fit)) {
    return(NULL)
  }
  y1 <- grid$year == 1L
  om <- predict(fit, newdata = grid[y1, ])$omega_s
  rp <- tidy(fit, "ran_pars")
  get_index(fit, newdata = grid, bias_correct = TRUE) |>
    mutate(
      delta_hat = mean(om[is_north[y1]]) - mean(om[!is_north[y1]]),
      sigma_O_hat = rp$estimate[rp$term == "sigma_O"]
    )
}

run_task <- function(delta, sigma_O, seed) {
  RhpcBLASctl::blas_set_num_threads(1L)
  RhpcBLASctl::omp_set_num_threads(1L)

  mesh_sim <- make_square_mesh(grid, max_edge = 0.02)
  s <- simulate_cell(delta, sigma_O, seed, mesh_sim)
  truth <- s |>
    group_by(year) |>
    summarise(true = sum(mu))

  d <- observe(
    s,
    sample_n = obs_n, seed = seed, gap = gap, region_cutoff = region_cutoff,
    north_yrs = north_yrs, south_yrs = south_yrs
  )
  mesh <- make_square_mesh(d, max_edge = 0.05)

  purrr::map(models, \(m) index_and_contrast(fit_model(m, d, mesh))) |>
    setNames(models) |>
    bind_rows(.id = "model") |>
    left_join(truth, by = "year") |>
    mutate(delta = delta, sigma_O = sigma_O, seed = seed, .before = 1L)
}

if (FALSE) {
  tictoc::tic()
  x <- run_task(delta = 1, sigma_O = 2, seed = 1)
  tictoc::toc()
}

# Run (one cache file per task; delete files to refit) -------------------------

cache_dir <- here::here("data-generated/contrast-grid-cache")
dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)

tasks <- tidyr::expand_grid(delta = deltas, sigma_O = sigma_Os, seed = seeds) |>
  mutate(cache_file = file.path(
    cache_dir, sprintf("delta-%s__sigma-O-%s__seed-%03d.rds", delta, sigma_O, seed)
  ))
todo <- filter(tasks, !file.exists(cache_file))

if (nrow(todo) > 0L) {
  future::plan(future::multicore, workers = max(1L, min(25L, future::availableCores() - 2L)))
  furrr::future_pwalk(
    todo,
    \(delta, sigma_O, seed, cache_file) saveRDS(run_task(delta, sigma_O, seed), cache_file),
    .options = furrr::furrr_options(seed = TRUE, scheduling = Inf)
  )
  future::plan(future::sequential)
}

out <- purrr::map_dfr(tasks$cache_file[file.exists(tasks$cache_file)], readRDS)
saveRDS(out, here::here("data-generated/contrast-grid.rds"))

# Metrics ----------------------------------------------------------------------

# Alternation follows year parity (north odd, south even), the default phase
A_of <- \(x, year) period2_metric(x, year = year)[["A"]]

metrics <- out |>
  arrange(delta, sigma_O, seed, model, year) |>
  group_by(delta, sigma_O, seed, model) |>
  summarise(
    A = A_of(est, year),
    A_true = A_of(true, year),
    rmse = sqrt(mean((log_est - log(true))^2)),
    ratio_hat = abs(first(delta_hat)) / first(sigma_O_hat),
    .groups = "drop"
  ) |>
  mutate(ratio = delta / sigma_O, excess_A = A - A_true)

cell_summary <- metrics |>
  group_by(model, delta, sigma_O, ratio) |>
  summarise(
    n = n(),
    excess_med = median(excess_A),
    excess_lo = quantile(excess_A, 0.25),
    excess_hi = quantile(excess_A, 0.75),
    rmse = mean(rmse),
    ratio_hat = median(ratio_hat),
    .groups = "drop"
  )
print(cell_summary, n = Inf)

# Plots ------------------------------------------------------------------------

# Excess seesaw against Delta / sigma_O; collapse of the sigma_O lines onto
# one curve means the ratio, not Delta or sigma_O alone, drives the seesaw
ggplot(cell_summary, aes(ratio, excess_med, colour = factor(sigma_O))) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey60") +
  geom_linerange(aes(ymin = excess_lo, ymax = excess_hi)) +
  geom_line() +
  geom_point() +
  scale_x_continuous(trans = "log2") +
  scale_colour_viridis_d(end = 0.85) +
  facet_wrap(~model) +
  labs(
    x = expression(Delta / sigma[O]),
    y = "Excess seesaw A over truth (%)",
    colour = expression(sigma[O]),
    caption = "Median and IQR across seeds"
  )
ggsave(here::here("figs/contrast-grid-excess-A.pdf"), width = 7, height = 3.5)

# Diagnostic: estimated ratio from each biennial fit against the true ratio
ggplot(metrics, aes(ratio, ratio_hat, colour = factor(sigma_O))) +
  geom_abline(linetype = 2, colour = "grey60") +
  geom_point(position = position_jitter(width = 0.03), alpha = 0.3) +
  scale_x_continuous(trans = "log2") +
  scale_y_continuous(trans = "log2") +
  scale_colour_viridis_d(end = 0.85) +
  facet_wrap(~model) +
  labs(
    x = expression(True ~ Delta / sigma[O]),
    y = expression(Estimated ~ Delta / sigma[O]),
    colour = expression(sigma[O])
  )
ggsave(here::here("figs/contrast-grid-ratio-diagnostic.pdf"), width = 7, height = 3.5)
