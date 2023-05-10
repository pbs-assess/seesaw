# Test model on survey data
check_sanity <- function(x) {
  if (!all(unlist(sanity(x)))) {
    return(NA)
  } else {
    return(x)
  }
}

fit_models <- function(
    dat, catch, data_subset = NULL, mesh = NULL, cutoff = 20, family = tweedie(),
    offset = NULL, use_extra_time = TRUE, silent = TRUE,
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

  missing_years <- NULL
  if (use_extra_time) {
    missing_years <- sdmTMB:::find_missing_time(dat$year)
    #dat_years <- unique(dat$year)
    #all_years <- min(dat$year):max(dat$year)
    #missing_years <- all_years[!(all_years%in% dat_years)]
    message(cat("\t- Filling in extra_time with:", missing_years, "\n"))
  }

  dat <- droplevels(dat)  # drop extra factor levels before running models
  fits <- list()
  model_ids <- c()
  i <- 1

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

  cli::cli_inform("\tFitting 2: st IID covariate")
  fits[[i]] <- try(
    sdmTMB(
      eval(parse(text = catch)) ~ 0 + as.factor(year) + log_depth + I(log_depth^2),
      family = family,
      data = dat, time = "year", spatiotemporal = "iid", spatial = "on",
      silent = silent, mesh = mesh,
      offset = offset,
      control = ctrl
    )
  )
  model_ids <- c(model_ids, "st IID covariate")
  i <- i + 1

  cli::cli_inform("\tFitting 3: st IID s(year)")
  fits[[i]] <- try(
    sdmTMB(
      eval(parse(text = catch)) ~ s(year),
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
  model_ids <- c(model_ids, "st time_varying RW")
  i <- i + 1

  cli::cli_inform("\tFitting 6: st (1|year)")
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

  cli::cli_inform("\tFitting 7: spatial only")
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

  names(fits) <- paste(model_ids, data_subset, sep = ":")
  fits
}

# is_even <- function(column) {
#   ifelse({{column}} %% 2 == 0, TRUE, FALSE)
# }

get_pred_list <- function(fit_list, newdata) {
  fit_list %>%
  purrr::map(., function(.x) {
    if (inherits(.x, "sdmTMB")) {
      newdata <- newdata %>%
        filter(survey %in% unique(.x$data$survey_abbrev),
               year %in% unique(.x$data$year)
        ) %>%
        droplevels()
      out <- predict(.x, newdata = newdata, return_tmb_object = TRUE, extra_time = .x$extra_time)
      out$newdata_input <- newdata
    } else {
      out <- NA
    }
    out
  })
}

get_index_list <- function(pred_list) {
  purrr::map(pred_list, function(.x) {
    if (length(.x) > 1) {
      out <- get_index(.x, bias_correct = TRUE, area = .x$newdata_input$area)
    } else {
      out <- NA  # keep empty fits as visual cue that these did not fit when plotting
    }
  })
}

mk_index_df <- function(index_list) {
  enframe(index_list) %>%
    unnest(col = "value") %>%
    separate(col = 'name', into = c('id', 'group'), sep = ":") %>%
    mutate(id = as.numeric(id)) %>%
    right_join(., model_lookup)
}
