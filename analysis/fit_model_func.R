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
    offset = NULL, use_extra_time = TRUE,
    ctrl = sdmTMBcontrol(nlminb_loops = 1L, newton_loops = 1L)) {

  if (is.null(data_subset)) {
    data_subset <- unique(dat$species_common_name)
  } else {
    data_subset <- unique(dat[[data_subset]])
  }

  message(cat("\n\tFitting models for data subset:", data_subset, "\n"))

  if (is.null(mesh)) {
    message(cat("\t- No mesh provided, making mesh with cutoff:", cutoff))
    mesh <- make_mesh(dat, c("X", "Y"), cutoff = 20)
  }

  missing_years <- NULL
  if (use_extra_time) {
    #missing_years <- sdmTMB:::find_missing_time(dat$year)
    dat_years <- unique(dat$year)
    all_years <- min(dat$year):max(dat$year)
    missing_years <- all_years[!(all_years%in% dat_years)]
    message(cat("\t- Filling in extra_time with:", missing_years, "\n"))
  }

  dat <- droplevels(dat)  # drop extra factor levels before running models
  fits <- list()
  i <- 1

  cli::cli_inform("\tFitting 1: st = 'rw'")
  fits[[i]] <- try(
    sdmTMB(
      eval(parse(text = catch)) ~ 1,
      family = family,
      data = dat, time = "year", spatiotemporal = "rw", spatial = "on",
      silent = TRUE, mesh = mesh,
      offset = offset,
      extra_time = missing_years,
      control = ctrl
    )
  )
  model_ids <- i
  i <- i + 1

  cli::cli_inform("\tFitting 2: st IID covariate")
  fits[[i]] <- try(
    sdmTMB(
      eval(parse(text = catch)) ~ 0 + as.factor(year) + log_depth + I(log_depth^2),
      family = family,
      data = dat, time = "year", spatiotemporal = "iid", spatial = "on",
      silent = TRUE, mesh = mesh,
      offset = offset,
      control = ctrl
    )
  )
  model_ids <- c(model_ids, i)
  i <- i + 1

  cli::cli_inform("\tFitting 3: st IID s(year)")
  fits[[i]] <- try(
    sdmTMB(
      eval(parse(text = catch)) ~ s(year),
      family = family,
      data = dat, time = "year", spatiotemporal = "iid", spatial = "on",
      silent = TRUE, mesh = mesh,
      offset = offset,
      control = ctrl
    )
  )
  model_ids <- c(model_ids, i)
  i <- i + 1

  cli::cli_inform("\tFitting 4: st IID no covariate as.factor year")
  fits[[i]] <- try(
    sdmTMB(
      eval(parse(text = catch)) ~ 0 + as.factor(year),
      family = family,
      data = dat, time = "year", spatiotemporal = "iid", spatial = "on",
      mesh = mesh,
      offset = offset,
      control = ctrl
    )
  )
  model_ids <- c(model_ids, i)
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
      extra_time = missing_years,
      control = ctrl
    )
  )
  model_ids <- c(model_ids, i)
  i <- i + 1

  cli::cli_inform("\tFitting 6: st (1|year)")
  fits[[i]] <- try(
    sdmTMB(
      eval(parse(text = catch)) ~ 1 + (1 | fyear),
      family = family,
      data = dat, time = "year", spatiotemporal = "iid", spatial = "on",
      mesh = mesh,
      offset = offset,
      control = ctrl
    )
  )
  model_ids <- c(model_ids, i)
  i <- i + 1

  cli::cli_inform("\tFitting 7: spatial only")
  fits[[i]] <- try(
    sdmTMB(
      eval(parse(text = catch)) ~ 0 + as.factor(year),
      family = family,
      data = dat, time = "year", spatiotemporal = "off", spatial = "on",
      mesh = mesh,
      offset = offset,
      control = ctrl
    )
  )
  model_ids <- c(model_ids, i)
  i <- i + 1

  cli::cli_inform("\tFitting 8: st (1 | region)")
  fits[[i]] <- try(
    sdmTMB(
      eval(parse(text = catch)) ~ 0 + fyear + region,
      family = family,
      data = dat, time = "year", spatiotemporal = "iid", spatial = "on",
      mesh = mesh,
      offset = offset,
      control = ctrl
    )
  )
  model_ids <- c(model_ids, i)

  names(fits) <- paste(model_ids, data_subset, sep = ":")
  fits 
}

# is_even <- function(column) {
#   ifelse({{column}} %% 2 == 0, TRUE, FALSE)
# }

# -------
get_pred_list <- function(fit_list, newdata, newdata_extra_time = NULL) {
  fit_list %>% 
  #purrr::map(., function(.x) {
  furrr::future_map(., function(.x) {
    newdata <- newdata
    if (inherits(.x, "sdmTMB")) {
      if(!is.null(.x$extra_time)) newdata <- newdata_extra_time
      out <- predict(.x, newdata = newdata, return_tmb_object = TRUE)
    } else {
      out <- NA
    }
    out
  })
}

get_index_list <- function(pred_list) {
  #purrr::map(pred_list, function(.x) {
  furrr::future_map(pred_list, function(.x) {
    if (length(.x) > 1) {
      get_index(.x, bias_correct = TRUE, area = .x$data$area)
    } else {
      out <- NA  # keep empty fits as visual cue that these did not fit when plotting
    }
  })
}

mk_index_df <- function(index_list) {
  enframe(index_list) %>% 
  unnest(col = "value") %>% 
  separate(col = 'name', into = c('id', 'species'), sep = ":") %>% 
  mutate(id = as.numeric(id)) %>% 
  right_join(., model_lookup)
}