# A population on the unit square with a north-south difference Delta in its
# time-invariant spatial effect, sampled in one region per year (north in odd
# years, south in even)

source(here::here("analysis/R/models.R"))

gradient_deltas <- c(0, 0.5, 1, 2) # 0 = no regional contrast (control)
gradient_sigma_Os <- c(1, 2, 4)
gradient_seeds <- seq_len(100L)
gradient_models <- c("IID RF, factor(year)", "IID RF, RW year")

gradient_n_year <- 12L
gradient_obs_n <- 200L # samples per year, from the sampled half
region_cutoff <- 0.5 # contiguous regions, as QCS / HS
north_years <- seq(1L, gradient_n_year - 1L, 2L)

gradient_grid <- expand.grid(
  X = seq(0, 1, length.out = 100), Y = seq(0, 1, length.out = 100),
  year = seq_len(gradient_n_year)
)
gradient_grid$region <- ifelse(gradient_grid$Y >= region_cutoff, "north", "south")

# Outer buffer to limit SPDE boundary effects
make_square_mesh <- function(dat, max_edge, outer_max_edge = 0.15) {
  m <- fmesher::fm_mesh_2d_inla(
    loc.domain = cbind(c(0, 1, 1, 0), c(0, 0, 1, 1)),
    max.edge = c(max_edge, outer_max_edge),
    offset = c(0.03, 0.25)
  )
  suppressMessages(make_mesh(dat, xy_cols = c("X", "Y"), mesh = m))
}

# The year effects are shifted so the gradient doesn't change the mean
# density. Delta is the realized north-minus-south difference in the mean
# spatial effect (omega_s + slope * Y), set exactly by solving for the slope
# after a pilot simulation with the same field draws.
simulate_gradient <- function(delta, sigma_O, seed, mesh,
                              matern_range = 0.1, sigma_E = 1, phi = 10, tweedie_p = 1.5) {
  grid <- gradient_grid
  set.seed(seed * 927849)
  year_ar1 <- as.numeric(arima.sim(model = list(ar = 0.5), n = gradient_n_year, sd = sqrt(1 - 0.5^2))) * 0.3
  # the gradient covariate is a copy of Y (simulate_new() adds a column
  # per coefficient, which would duplicate Y)
  grid$gradient <- grid$Y
  sim_one <- function(slope) {
    set.seed(seed * 1029)
    s <- simulate_new(
      formula = ~ 0 + factor(year) + gradient,
      data = grid, time = "year", mesh = mesh, family = tweedie(),
      range = matern_range, sigma_O = sigma_O, sigma_E = sigma_E, rho = 0,
      phi = phi, tweedie_p = tweedie_p,
      B = c(1 - slope * mean(grid$Y) + year_ar1, slope)
    )
    tibble(s[, c("year", "X", "Y", "omega_s", "epsilon_st", "eta", "mu", "observed")], region = grid$region)
  }
  y1 <- grid$year == 1L
  north <- grid$region == "north"
  mean_Y_diff <- mean(grid$Y[y1 & north]) - mean(grid$Y[y1 & !north])
  pilot <- sim_one(0)
  omega_diff <- mean(pilot$omega_s[y1 & north]) - mean(pilot$omega_s[y1 & !north])
  slope <- (delta - omega_diff) / mean_Y_diff
  s <- sim_one(slope)
  stopifnot(isTRUE(all.equal(s$omega_s, pilot$omega_s)))
  spatial <- s$omega_s + slope * s$Y
  stopifnot(abs(mean(spatial[y1 & north]) - mean(spatial[y1 & !north]) - delta) < 1e-6)
  s
}

observe_alternating <- function(s, seed) {
  set.seed(seed * 7)
  s |>
    filter((year %in% north_years) == (region == "north")) |>
    slice_sample(n = gradient_obs_n, by = year) |>
    arrange(year, X, Y)
}

gradient_phase <- function(year) ifelse(year %in% north_years, 0.5, -0.5)

fit_gradient_model <- function(model, d, mesh) {
  fit_model(model, data = d, mesh = mesh, response = "observed", family = tweedie(), anisotropy = FALSE)
}

run_gradient_task <- function(delta, sigma_O, seed, return_fits = FALSE) {
  single_thread()
  # The grid runs in furrr workers, which use L'Ecuyer-CMRG; pinning it here
  # gives the same draws from the main session (fit-gradient-example.R)
  old_kind <- RNGkind("L'Ecuyer-CMRG")
  on.exit(RNGkind(old_kind[1L]), add = TRUE)
  s <- simulate_gradient(delta, sigma_O, seed, make_square_mesh(gradient_grid, max_edge = 0.02))
  truth <- summarise(s, true_est = sum(mu), .by = year)
  d <- observe_alternating(s, seed)
  mesh <- make_square_mesh(d, max_edge = 0.05)
  fits <- map(set_names(gradient_models), \(m) fit_ok(fit_gradient_model(m, d, mesh)))
  indexes <- compact(fits) |>
    map(\(f) tryCatch(get_index(f, newdata = gradient_grid, bias_correct = TRUE), error = \(e) NULL)) |>
    compact() |>
    bind_rows(.id = "model") |>
    bind_rows(tibble(model = character(), year = integer())) |>
    select(any_of(c("model", "year", "est", "lwr", "upr", "log_est", "se"))) |>
    left_join(truth, by = "year") |>
    mutate(phase = gradient_phase(year), delta = delta, sigma_O = sigma_O, seed = seed, .before = 1L)
  if (return_fits) list(s = s, fits = fits, indexes = indexes) else indexes
}
