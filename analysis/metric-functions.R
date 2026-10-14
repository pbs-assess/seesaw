period2_metric <- function(x, year = seq_along(x), conf = 0.95, phase = NULL) {
  stopifnot(
    length(x) == length(year),
    all(x > 0)
  )

  y <- log(x)

  # Alternating survey phases, coded so beta is the phase difference.
  # Defaults to calendar-year parity; pass `phase` (+/- 0.5) when the
  # alternation doesn't follow parity (e.g., HBLL OUT after the 2013 gap).
  if (is.null(phase)) {
    phase <- ifelse(year %% 2 == 0, 0.5, -0.5)
  }
  stopifnot(length(phase) == length(x), all(phase %in% c(-0.5, 0.5)))
  # period2 <- (-1)^year

  fit <- lm(y ~ year + phase)

  beta <- coef(fit)[["phase"]]
  beta_se <- coef(summary(fit))["phase", "Std. Error"]
  ci <- confint(fit, "phase", level = conf)
  A_fun <- \(b) 100 * (exp(abs(b)) - 1)
  A_ci <- A_fun(ci)
  A_bounds <- if (ci[1L] <= 0 && ci[2L] >= 0) {
    c(0, max(A_ci))
  } else {
    range(A_ci)
  }

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

# Windows span a fixed number of calendar years, so a window with a missing
# year contains fewer observations
moving_window <- function(x, year, window = 10L, phase = NULL, ...) {
  starts <- seq(min(year), max(year) - window + 1L)

  out <- lapply(starts, function(y0) {
    ind <- which(year >= y0 & year < y0 + window)

    res <- period2_metric(
      x = x[ind],
      year = year[ind],
      phase = if (is.null(phase)) NULL else phase[ind],
      ...
    )

    data.frame(
      start = y0,
      end = y0 + window - 1L,
      n = length(ind),
      t(res)
    )
  })

  do.call(rbind, out)
}
