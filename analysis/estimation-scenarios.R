build_model_specs <- function() {
  list(
    list(
      name = "RW RF",
      fit_args = list(
        formula = observed ~ 1,
        spatiotemporal = "rw"
      )
    ),
    list(
      name = "AR(1) RF",
      fit_args = list(
        formula = observed ~ 1,
        spatiotemporal = "ar1"
      )
    ),
    list(
      name = "AR(1) RF, factor(year)",
      fit_args = list(
        formula = observed ~ 0 + as.factor(year),
        spatiotemporal = "ar1"
      )
    ),
    list(
      name = "AR(1) RF, RW year",
      fit_args = list(
        formula = observed ~ 0,
        time_varying = ~1,
        time_varying_type = "rw",
        spatiotemporal = "ar1"
      )
    ),
    list(
      name = "RW RF, RW year",
      fit_args = list(
        formula = observed ~ 0,
        time_varying = ~1,
        time_varying_type = "rw",
        spatiotemporal = "rw"
      )
    ),
    list(
      name = "RW RF, fixed 0.3 RW year",
      fit_args = list(
        formula = observed ~ 0,
        time_varying = ~1,
        time_varying_type = "rw",
        spatiotemporal = "rw"
      ),
      control = sdmTMB::sdmTMBcontrol(
        start = list(ln_tau_V = matrix(log(0.3), nrow = 1, ncol = 1L)),
        map = list(ln_tau_V = rep(factor(NA), 1L))
      )
    ),
    list(
      name = "RW RF, fixed 0.1 RW year",
      fit_args = list(
        formula = observed ~ 0,
        time_varying = ~1,
        time_varying_type = "rw",
        spatiotemporal = "rw"
      ),
      control = sdmTMB::sdmTMBcontrol(
        start = list(ln_tau_V = matrix(log(0.1), nrow = 1, ncol = 1L)),
        map = list(ln_tau_V = rep(factor(NA), 1L))
      )
    ),
    list(
      name = "RW RF, factor(year_pairs)",
      data_key = "pairs",
      fit_args = list(
        formula = observed ~ 0 + as.factor(year_pairs),
        spatiotemporal = "rw"
      )
    ),
    list(
      name = "IID RF, N-S gradient estimated",
      fit_args = list(
        formula = observed ~ 0 + as.factor(year) + depth_cov,
        spatiotemporal = "iid"
      )
    ),
    list(
      name = "IID RF, s(year)",
      fit_args = list(
        formula = observed ~ s(year),
        spatiotemporal = "iid"
      )
    ),
    list(
      name = "IID RF, factor(year)",
      fit_args = list(
        formula = observed ~ 0 + factor(year),
        spatiotemporal = "iid"
      )
    ),
    list(
      name = "IID RF, factor(year) + factor(region)",
      fit_args = list(
        formula = observed ~ 0 + factor(year) + factor(region),
        spatiotemporal = "iid"
      )
    ),
    list(
      name = "IID RF, RW year",
      fit_args = list(
        formula = observed ~ 0,
        time_varying = ~1,
        time_varying_type = "rw",
        spatiotemporal = "iid"
      )
    ),
    list(
      name = "IID RF, RW year + gamma(0.3, 0.5) prior",
      fit_args = list(
        formula = observed ~ 0,
        time_varying = ~1,
        time_varying_type = "rw",
        spatiotemporal = "iid",
        priors = sdmTMB::sdmTMBpriors(sigma_V = gamma_cv(0.3, 0.5))
      )
    ),
    list(
      name = "IID RF, AR1 year",
      fit_args = list(
        formula = observed ~ 1,
        time_varying = ~1,
        time_varying_type = "ar1",
        spatiotemporal = "iid"
      )
    ),
    list(
      name = "IID RF, (1|year)",
      fit_args = list(
        formula = observed ~ 1 + (1 | fyear),
        spatiotemporal = "iid"
      )
    ),
    list(
      name = "Spatial only, factor(year)",
      fit_args = list(
        formula = observed ~ 0 + as.factor(year),
        spatiotemporal = "off"
      )
    ),
    list(
      name = "Spatial only, RW year",
      fit_args = list(
        formula = observed ~ 0,
        time_varying = ~1,
        time_varying_type = "rw",
        spatiotemporal = "off"
      )
    ),
    list(
      name = "Spatial only, AR(1) year",
      fit_args = list(
        formula = observed ~ 1,
        time_varying = ~1,
        time_varying_type = "ar1",
        spatiotemporal = "off"
      )
    )
  )
}

if (FALSE) {
  dgamma_cv <- function(x, mean, cv, log = FALSE) {
    shape <- 1 / cv^2
    rate  <- 1 / (mean * cv^2)   # equivalently: shape / mean
    dgamma(x, shape = shape, rate = rate, log = log)
  }
  rgamma_cv <- function(n, mean, cv) {
    rgamma(n, shape = 1/cv^2, rate = 1/(mean*cv^2))
  }
  x <- seq(0.001, 1, length.out = 300)
  plot(x, dgamma_cv(x, 0.3, 0.5), type = "l")
  abline(v = 0.3)
}
