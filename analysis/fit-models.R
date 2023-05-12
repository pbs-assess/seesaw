fit_models <- function(
    dat, catch, data_subset = NULL, mesh = NULL, cutoff = 20, family = tweedie(),
    offset = NULL, silent = TRUE,
    ctrl = sdmTMBcontrol(nlminb_loops = 1L, newton_loops = 1L)) {
  if (is.null(data_subset)) {
    data_subset <- unique(dat$species_common_name)
  } else {
    data_subset <- unique(dat[[data_subset]])
  }

  if (!is.null(offset)) offset <- dat[[offset]]

  message(cat("\n\tFitting models for data subset:", data_subset, "\n"))

  if (is.null(mesh)) {
    message(cat("\t- No mesh provided, making mesh with cutoff:", cutoff))
    mesh <- make_mesh(dat, c("X", "Y"), cutoff = 20)
  }

  missing_years <- sdmTMB:::find_missing_time(dat$year)

  if (length(missing_years) < 1L) {
    missing_years <- NULL
    message(cat("\t- No missing time to be filled in."))
  } else {
    message(cat("\t- Filling in extra_time with:", missing_years, "\n"))
  }

  dat <- droplevels(dat) # drop extra factor levels before running models
  fits <- list()
  model_ids <- c()
  i <- 1

  #cli::cli_inform("\tFitting 7: spatial only")
  message("\tFitting", i, ": spatial only")
  fits[[i]] <- try(
    sdmTMB(
      eval(parse(text = catch)) ~ 0 + as.factor(year),
      family = family,
      data = dat, time = "year", spatiotemporal = "off", spatial = "on",
      mesh = mesh,
      silent = silent,
      offset = offset,
      control = ctrl
    )
  )
  model_ids <- c(model_ids, "spatial only")
  i <- i + 1

  cli::cli_inform("\tFitting 1: st = 'rw'")
  fits[[i]] <- try(
    sdmTMB(
      eval(parse(text = catch)) ~ 1,
      family = family,
      data = dat, time = "year", spatiotemporal = "rw", spatial = "on",
      silent = silent, mesh = mesh,
      offset = offset,
      extra_time = missing_years,
      control = ctrl
    )
  )
  model_ids <- c(model_ids, "st = 'rw'")
  i <- i + 1

  cli::cli_inform("\tFitting 1: st = 'ar1'")
  fits[[i]] <- try(
    sdmTMB(
      eval(parse(text = catch)) ~ 1,
      family = family,
      data = dat, time = "year", spatiotemporal = "ar1", spatial = "on",
      silent = silent, mesh = mesh,
      offset = offset,
      extra_time = missing_years,
      control = ctrl
    )
  )
  model_ids <- c(model_ids, "st = 'ar1'")
  i <- i + 1

  cli::cli_inform("\tFitting 2: st IID covariate")
  fits[[i]] <- try(
    update(fits[[1]],
      formula = eval(parse(text = catch)) ~ 0 + as.factor(year) + poly(log_depth, 2),
      spatiotemporal = "iid", spatial = "on"
    )
  )
  model_ids <- c(model_ids, "st IID covariate")
  i <- i + 1

  cli::cli_inform("\tFitting 3: st IID s(year)")
  fits[[i]] <- try(
    sdmTMB(
      eval(parse(text = catch)) ~ s(year, k = 5),
      family = family,
      data = dat, time = "year", spatiotemporal = "iid", spatial = "on",
      silent = silent, mesh = mesh,
      offset = offset,
      control = ctrl
    )
  )
  model_ids <- c(model_ids, "st IID s(year)")
  i <- i + 1

  cli::cli_inform("\tFitting 4: st IID no covariate as.factor year")
  fits[[i]] <- try(
    sdmTMB(
      eval(parse(text = catch)) ~ 0 + as.factor(year),
      family = family,
      data = dat, time = "year", spatiotemporal = "iid", spatial = "on",
      mesh = mesh,
      silent = silent,
      offset = offset,
      control = ctrl
    )
  )
  model_ids <- c(model_ids, "st IID no covariate as.factor year")
  i <- i + 1

  cli::cli_inform("\tFitting 5: st time_varying RW")
  fits[[i]] <- try(
    sdmTMB(
      eval(parse(text = catch)) ~ 0,
      family = family,
      time_varying = ~1, time_varying_type = "rw",
      data = dat, time = "year", spatiotemporal = "iid", spatial = "on",
      mesh = mesh,
      offset = offset,
      silent = silent,
      extra_time = missing_years,
      control = ctrl
    )
  )
  model_ids <- c(model_ids, "st time-varying RW")
  i <- i + 1

  cli::cli_inform("\tst time_varying RW, fixed SD")

  .dim <- if (isTRUE(family$delta)) 2 else 1
  fits[[i]] <- try(
    sdmTMB(
      eval(parse(text = catch)) ~ 0,
      family = family,
      time_varying = ~1, time_varying_type = "rw",
      data = dat, time = "year", spatiotemporal = "iid", spatial = "on",
      mesh = mesh,
      offset = offset,
      silent = silent,
      extra_time = missing_years,
      control = sdmTMBcontrol(
        start = list(ln_tau_V = matrix(log(0.1), nrow = 1, ncol = .dim)),
        map = list(ln_tau_V = rep(factor(NA), .dim))
      )
    )
  )
  model_ids <- c(model_ids, "st time-varying RW; fixed 0.1 SD")
  i <- i + 1

  cli::cli_inform("\tspatial time-varying RW")
  fits[[i]] <- try(
    sdmTMB(
      eval(parse(text = catch)) ~ 0,
      family = family,
      time_varying = ~1, time_varying_type = "rw",
      data = dat, time = "year", spatiotemporal = "iid", spatial = "on",
      mesh = mesh,
      offset = offset,
      silent = silent,
      extra_time = missing_years,
      control = ctrl
    )
  )
  model_ids <- c(model_ids, "spatial time-varying RW")
  i <- i + 1

  cli::cli_inform("\tspatial time-varying AR(1)")
  fits[[i]] <- try(
    sdmTMB(
      eval(parse(text = catch)) ~ 0,
      family = family,
      time_varying = ~1, time_varying_type = "ar1",
      data = dat, time = "year", spatiotemporal = "iid", spatial = "on",
      mesh = mesh,
      offset = offset,
      silent = silent,
      extra_time = missing_years,
      control = ctrl
    )
  )
  model_ids <- c(model_ids, "spatial time-varying AR(1)")
  i <- i + 1

  cli::cli_inform("\tFitting : st (1|year)")
  fits[[i]] <- try(
    sdmTMB(
      eval(parse(text = catch)) ~ 1 + (1 | fyear),
      family = family,
      data = dat, time = "year", spatiotemporal = "iid", spatial = "on",
      mesh = mesh,
      silent = silent,
      offset = offset,
      control = ctrl
    )
  )
  model_ids <- c(model_ids, "st (1|year)")
  i <- i + 1

  cli::cli_inform("\tFitting 8: st (1 | region)")
  fits[[i]] <- try(
    sdmTMB(
      eval(parse(text = catch)) ~ 0 + fyear + region,
      family = family,
      data = dat, time = "year", spatiotemporal = "iid", spatial = "on",
      mesh = mesh,
      silent = silent,
      offset = offset,
      control = ctrl
    )
  )
  model_ids <- c(model_ids, "st (1 | region)")

  cli::cli_inform("\tFitting: st (1 | year-pairs) AR1")
  fits[[i]] <- try(
    sdmTMB(
      eval(parse(text = catch)) ~ 0 + as.factor(year_pair),
      family = family,
      data = dat, time = "year", spatiotemporal = "AR1", spatial = "on",
      mesh = mesh,
      silent = silent,
      offset = offset,
      extra_time = missing_years,
      control = ctrl
    )
  )
  model_ids <- c(model_ids, "st (1 | year-pairs) RW")

  cli::cli_inform("\tFitting: st (1 | year-pairs) RW")
  fits[[i]] <- try(
    sdmTMB(
      eval(parse(text = catch)) ~ 0 + as.factor(year_pair),
      family = family,
      data = dat, time = "year", spatiotemporal = "RW", spatial = "on",
      mesh = mesh,
      silent = silent,
      offset = offset,
      extra_time = missing_years,
      control = ctrl
    )
  )
  model_ids <- c(model_ids, "st (1 | year-pairs) RW")

  names(fits) <- paste(model_ids, data_subset, sep = ":")
  fits
}
