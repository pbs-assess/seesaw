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

  message("\tFitting", i, ": spatial only, as.factor(year)")
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
  model_ids <- c(model_ids, "spatial only, as.factor(year)")
  i <- i + 1

  cli::cli_inform("\tFitting: st RW")
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
  model_ids <- c(model_ids, "st RW")
  i <- i + 1

  cli::cli_inform("\tFitting: st AR1")
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
  model_ids <- c(model_ids, "st AR1")
  i <- i + 1

  cli::cli_inform("\tFitting: st IID, as.factor(year)")
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
  model_ids <- c(model_ids, "st IID, as.factor(year)")
  i <- i + 1

  cli::cli_inform("\tFitting: st IID, as.factor(year) w/ depth")
  fits[[i]] <- try(
    sdmTMB(
      formula = eval(parse(text = catch)) ~ 0 + as.factor(year) + poly(log_depth, 2),
      family = family,
      data = dat, time = "year", spatiotemporal = "iid", spatial = "on",
      mesh = mesh,
      silent = silent,
      offset = offset,
      control = ctrl
    )
  )
  model_ids <- c(model_ids, "st IID, as.factor(year) w/ depth")
  i <- i + 1

  cli::cli_inform("\tFitting 3: st IID, s(year)")
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
  model_ids <- c(model_ids, "st IID, s(year)")
  i <- i + 1

  cli::cli_inform("\tFitting: st IID, time_varying RW")
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
  model_ids <- c(model_ids, "st IID, time-varying RW")
  i <- i + 1

  cli::cli_inform("\tFitting: st RW, time_varying RW")
  fits[[i]] <- try(
    sdmTMB(
      eval(parse(text = catch)) ~ 0,
      family = family,
      time_varying = ~1, time_varying_type = "rw",
      data = dat, time = "year", spatiotemporal = "rw", spatial = "on",
      mesh = mesh,
      offset = offset,
      silent = silent,
      extra_time = missing_years,
      control = ctrl
    )
  )
  model_ids <- c(model_ids, "st RW, time-varying RW")
  i <- i + 1

  cli::cli_inform("\tst IID, time_varying RW, fixed SD")
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
  model_ids <- c(model_ids, "st IID, time-varying RW; fixed 0.1 SD")
  i <- i + 1

  cli::cli_inform("\tst IID, time-varying RW")
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
  model_ids <- c(model_ids, "st IID, time-varying RW")
  i <- i + 1

  cli::cli_inform("\tst IID, time-varying AR(1)")
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
  model_ids <- c(model_ids, "st IID, time-varying AR(1)")
  i <- i + 1

  cli::cli_inform("\tFitting : st IID, (1|year)")
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
  model_ids <- c(model_ids, "st IID, (1|year)")
  i <- i + 1

  cli::cli_inform("\tFitting 1: st RW, (1|year)")
  fits[[i]] <- try(
    sdmTMB(
      eval(parse(text = catch)) ~ 1 + (1 | fyear),
      family = family,
      data = dat, time = "year", spatiotemporal = "rw", spatial = "on",
      silent = silent, mesh = mesh,
      offset = offset,
      extra_time = missing_years,
      control = ctrl
    )
  )
  model_ids <- c(model_ids, "st RW, (1|year)")
  i <- i + 1

  cli::cli_inform("\tFitting 1: st AR1, (1|year)")
  fits[[i]] <- try(
    sdmTMB(
      eval(parse(text = catch)) ~ 1 + (1 | fyear),
      family = family,
      data = dat, time = "year", spatiotemporal = "ar1", spatial = "on",
      silent = silent, mesh = mesh,
      offset = offset,
      extra_time = missing_years,
      control = ctrl
    )
  )
  model_ids <- c(model_ids, "st AR1, (1|year)")
  i <- i + 1

  cli::cli_inform("\tFitting 8: st IID, (1|region)")
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
  model_ids <- c(model_ids, "st IID, (1|region)")
  i <- i + 1

  cli::cli_inform("\tFitting: st = AR1, as.factor(year_pair)")
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
  model_ids <- c(model_ids, "st = AR1, as.factor(year_pair)")
  i <- i + 1

  cli::cli_inform("\tFitting: st = RW, as.factor(year_pair)")
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
  model_ids <- c(model_ids, "st = RW, as.factor(year_pair)")

  names(fits) <- paste(model_ids, data_subset, sep = ":")
  fits
}
