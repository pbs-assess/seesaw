build_model_specs <- function() {
  list(
    # Independent years ---------------------------------------------------
    list(
      name = "IID RF, factor(year)",
      category = "Independent years",
      order = 1,
      fit_args = list(
        formula = observed ~ 0 + factor(year),
        spatiotemporal = "iid"
      )
    ),
    list(
      name = "IID RF, factor(year) + factor(region)",
      category = "Independent years",
      order = 2,
      fit_args = list(
        formula = observed ~ 0 + factor(year) + factor(region),
        spatiotemporal = "iid"
      )
    ),
    list(
      name = "IID RF, N-S gradient estimated",
      category = "Independent years",
      order = 3,
      fit_args = list(
        formula = observed ~ 0 + as.factor(year) + depth_cov,
        spatiotemporal = "iid"
      )
    ),
    list(
      name = "Spatial only, factor(year)",
      category = "Independent years",
      order = 4,
      fit_args = list(
        formula = observed ~ 0 + as.factor(year),
        spatiotemporal = "off"
      )
    ),
    # Smoothed years ------------------------------------------------------
    list(
      name = "IID RF, RW year",
      category = "Smoothed years",
      order = 5,
      fit_args = list(
        formula = observed ~ 0,
        time_varying = ~1,
        time_varying_type = "rw",
        spatiotemporal = "iid",
        priors = sdmTMB::sdmTMBpriors(sigma_V = sdmTMB::gamma_cv(0.3, 0.5))
      )
    ),
    list(
      name = "AR(1) RF, RW year",
      category = "Smoothed years",
      order = 6,
      fit_args = list(
        formula = observed ~ 0,
        time_varying = ~1,
        time_varying_type = "rw",
        spatiotemporal = "ar1",
        priors = sdmTMB::sdmTMBpriors(sigma_V = sdmTMB::gamma_cv(0.3, 0.5))
      )
    ),
    list(
      name = "RW RF, RW year",
      category = "Smoothed years",
      order = 7,
      fit_args = list(
        formula = observed ~ 0,
        time_varying = ~1,
        time_varying_type = "rw",
        spatiotemporal = "rw",
        priors = sdmTMB::sdmTMBpriors(sigma_V = sdmTMB::gamma_cv(0.3, 0.5))
      )
    ),
    list(
      name = "Spatial only, RW year",
      category = "Smoothed years",
      order = 8,
      fit_args = list(
        formula = observed ~ 0,
        time_varying = ~1,
        time_varying_type = "rw",
        spatiotemporal = "off"
      )
    ),
    # No year effects -----------------------------------------------------
    list(
      name = "AR(1) RF",
      category = "No year effects",
      order = 9,
      fit_args = list(
        formula = observed ~ 1,
        spatiotemporal = "ar1"
      )
    ),
    list(
      name = "RW RF",
      category = "No year effects",
      order = 10,
      fit_args = list(
        formula = observed ~ 1,
        spatiotemporal = "rw"
      )
    )
    # list(
    #   name = "AR(1) RF, factor(year)",
    #   fit_args = list(
    #     formula = observed ~ 0 + as.factor(year),
    #     spatiotemporal = "ar1"
    #   )
    # ),
    # list(
    #   name = "IID RF, s(year, k = 5)",
    #   fit_args = list(
    #     formula = observed ~ s(year, k = 5),
    #     spatiotemporal = "iid"
    #   )
    # ),
    # list(
    #   name = "IID RF, RW year",  # without prior
    #   fit_args = list(
    #     formula = observed ~ 0,
    #     time_varying = ~1,
    #     time_varying_type = "rw",
    #     spatiotemporal = "iid"
    #   )
    # ),
    # list(
    #   name = "IID RF, AR1 year",
    #   fit_args = list(
    #     formula = observed ~ 1,
    #     time_varying = ~1,
    #     time_varying_type = "ar1",
    #     spatiotemporal = "iid"
    #   )
    # ),
    # list(
    #   name = "IID RF, (1|year)",
    #   fit_args = list(
    #     formula = observed ~ 1 + (1 | fyear),
    #     spatiotemporal = "iid"
    #   )
    # ),
    # list(
    #   name = "RW RF, fixed 0.3 RW year",
    #   fit_args = list(
    #     formula = observed ~ 0,
    #     time_varying = ~1,
    #     time_varying_type = "rw",
    #     spatiotemporal = "rw"
    #   ),
    #   control = sdmTMB::sdmTMBcontrol(
    #     start = list(ln_tau_V = matrix(log(0.3), nrow = 1, ncol = 1L)),
    #     map = list(ln_tau_V = rep(factor(NA), 1L))
    #   )
    # ),
    # list(
    #   name = "RW RF, fixed 0.1 RW year",
    #   fit_args = list(
    #     formula = observed ~ 0,
    #     time_varying = ~1,
    #     time_varying_type = "rw",
    #     spatiotemporal = "rw"
    #   ),
    #   control = sdmTMB::sdmTMBcontrol(
    #     start = list(ln_tau_V = matrix(log(0.1), nrow = 1, ncol = 1L)),
    #     map = list(ln_tau_V = rep(factor(NA), 1L))
    #   )
    # ),
    # list(
    #   name = "RW RF, factor(year_pairs)",
    #   data_key = "pairs",
    #   fit_args = list(
    #     formula = observed ~ 0 + as.factor(year_pairs),
    #     spatiotemporal = "rw"
    #   )
    # ),
    # list(
    #   name = "Spatial only, AR(1) year",
    #   fit_args = list(
    #     formula = observed ~ 1,
    #     time_varying = ~1,
    #     time_varying_type = "ar1",
    #     spatiotemporal = "off"
    #   )
    # )
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
