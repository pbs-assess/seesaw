sim_ar1_heavy <- function(n, marginal_sd, rho, mu = 0,
  heavy_sd_mult = 1, heavy_sd_frac = 0, s = 42) {
  sigma <- sqrt(1 - rho^2)
  set.seed(s)
  devs <- rnorm(n) * marginal_sd
  set.seed(s)
  devs_heavy <- rnorm(n) * marginal_sd * heavy_sd_mult
  hd <- rbinom(n, size = 1, prob = heavy_sd_frac)
  this_dev <- if (hd[1] == 0) devs[1] else devs_heavy[1]
  x0 <- this_dev * marginal_sd
  B <- numeric(n)
  B[1] <- rho * x0 + sigma * devs[1]
  for (i in seq(2, length(B))) {
    this_dev <- if (hd[i] == 0) devs[i] else devs_heavy[i]
    B[i] <- rho * B[i-1] + sigma * this_dev
  }
  B + mu
}

sim <- function(predictor_grid, mesh, phi = 8, seed = 123,
                region_cutoff = 0.5, rho = 0, sigma_E = 0.4, range = 0.8, tweedie_p = 1.6,
                sigma_O = 1.6, coefs = c(2.2, 3.8), year_mean = 1, north_effect = 0,
                year_arima.sim = list(ar = 0.8), year_marginal_sd = 0.2,
                svc_trend = 0,
                heavy_sd_mult = 1,
                heavy_sd_frac = 0) {

  predictor_grid$region <- NA
  predictor_grid$region[predictor_grid$Y < region_cutoff] <- "south"
  predictor_grid$region[predictor_grid$Y >= region_cutoff] <- "north"
  predictor_grid$region <- factor(predictor_grid$region, levels = c("south", "north"))

  # yrs <- as.numeric(stats::arima.sim(
  #   n = max(predictor_grid$year),
  #   model = year_arima.sim
  # ), sd = sqrt(1 - year_arima.sim$ar1^2))
  # yrs <- yrs * year_marginal_sd + year_mean

  yrs <- sim_ar1_heavy(
    n = max(predictor_grid$year),
    rho = year_arima.sim$ar,
    mu = year_mean,
    marginal_sd = year_marginal_sd,
    heavy_sd_mult = heavy_sd_mult,
    heavy_sd_frac = heavy_sd_frac,
    s = seed * 927849
  )
  # plot(yrs, type = "o", ylim = c(-5, 5))

  predictor_grid$depth_cov <- predictor_grid$Y

  # visualize SVC trend:
  # predictor_grid |>
  #   mutate(year_cent = year - mean(year)) |>
  #   mutate(Y_cent = Y - mean(Y)) |>
  #   mutate(mu = Y_cent * year_cent * svc_trend) |>
  #   ggplot(aes(X, Y, fill = exp(mu))) +
  #   geom_raster() +
  #   facet_wrap(~year) +
  #   scale_fill_viridis_c(trans = "log10") +
  #   coord_fixed()

  predictor_grid <- predictor_grid |>
    mutate(year_cent = year - mean(year)) |>
    mutate(Y_cent = Y - mean(Y))

  if (north_effect != 0) {
    formula <- ~ 0 + as.factor(year) + region + Y_cent * year_cent
    B <- c(
      yrs,
      north_effect,
      0, # Y_cent
      0, # year_cent
      svc_trend # Y_cent:year_cent
    )
  } else {
    formula <- ~ 0 + as.factor(year) + depth_cov + I(depth_cov^2) + Y_cent * year_cent
    B <- c(
      yrs,
      coefs,
      0, # Y_cent
      0, # year_cent
      svc_trend # Y_cent:year_cent
    )
  }

  sim_dat <- sdmTMB_simulate(
    formula = formula,
    data = predictor_grid,
    time = "year",
    mesh = mesh,
    family = tweedie(),
    range = range,
    sigma_E = sigma_E,
    rho = rho,
    phi = phi,
    tweedie_p = tweedie_p,
    sigma_O = sigma_O,
    seed = seed * 1029,
    B = B
  )
  sim_dat$region <- predictor_grid$region
  sim_dat$year <- predictor_grid$year
  sim_dat$depth_cov <- predictor_grid$depth_cov
  list(sim_dat = sim_dat, predictor_dat = predictor_grid, year_effects = yrs, coefs = coefs, B = B)
}

# sample_before_split = TRUE means apply sample_n first, then discard any gap
# smaple_before_split = FALSE applies sample_n *after* cutting gap
observe <- function(sim_dat, sample_n = 300L, seed = 10282, gap = 0, region_cutoff = 0.5,
                    north_yrs = seq(1, 9, 2), south_yrs = seq(2, 10, 2),
                    sample_before_split = FALSE) {

  assertthat::assert_that(gap >= 0 && gap <= 0.99)
  assertthat::assert_that(sample_n > 1L)
  assertthat::assert_that(all(north_yrs %in% sim_dat$year) &&
    all(south_yrs %in% sim_dat$year))
  assertthat::assert_that(region_cutoff > 0.05 && region_cutoff < 0.95)

  both_yrs <- intersect(north_yrs, south_yrs)

  set.seed(seed * 7)

  if (sample_before_split) {
    d <- sim_dat %>%
      group_by(year) %>%
      sample_n(sample_n)
    # Lose ~half the survey most years -----------------------------------------
    north <- d[(d$year %in% union(north_yrs, both_yrs) & d$region == "north"), ]
    south <- d[(d$year %in% union(south_yrs, both_yrs) & d$region == "south"), ]
    d <- bind_rows(north, south) |>
      dplyr::arrange(year, X, Y)
    # Remove strip in middle to increase gap?
    d <- d[!(d$Y > (region_cutoff - gap / 2) & d$Y < (region_cutoff + gap / 2)), ]
  } else { # cut gap first, then assign samples
    d <- sim_dat
    north <- d[(d$year %in% union(north_yrs, both_yrs) & d$region == "north"), ]
    south <- d[(d$year %in% union(south_yrs, both_yrs) & d$region == "south"), ]
    d <- bind_rows(north, south) |>
      dplyr::arrange(year, X, Y)
    # Remove strip in middle to increase gap?
    d <- d[!(d$Y > (region_cutoff - gap / 2) & d$Y < (region_cutoff + gap / 2)), ]
    d <- d %>%
      group_by(year) %>%
      sample_n(sample_n)
  }
  d$sampled_region <- as.character(d$region)
  d$sampled_region[d$year %in% both_yrs] <- "both"
  d
}

sim_fit_and_index <- function(n_year,
                              .seed,
                              gap_size = 0.3,
                              obs_sampled_size = 400L,
                              year_marginal_sd = 0.5,
                              obs_yrs = list(
                                north_yrs = seq(1, n_year - 1, 2),
                                south_yrs = seq(2, n_year, 2)
                              ),
                              phi = 8,
                              region_cutoff = 0.5,
                              range = 0.5,
                              sigma_O = 1,
                              sigma_E = 0.5,
                              svc_trend = 0,
                              heavy_sd_mult = 1,
                              heavy_sd_frac = 0,
                              sample_before_split = FALSE,
                              year_arima.sim = list(ar = 0.5),
                              make_plots = FALSE,
                              sim_coefs = c(2, 5)) {
  is_even <- function(x) x %% 2 == 0
  if (!is_even(n_year)) cli::cli_abort("Number of years must be even.")

  cli::cli_alert_info(glue::glue("Running simulation for seed ", .seed))

  cli::cli_alert_success("Creating mesh...")
  predictor_grid <- expand.grid(
    X = seq(0, 1, length.out = 100), Y = seq(0, 1, length.out = 100),
    year = seq_len(n_year)
  )
  mesh_sim <- make_mesh(predictor_grid, xy_cols = c("X", "Y"), cutoff = 0.1)

  cli::cli_alert_success("Simulating...")
  x <- sim(predictor_grid, mesh_sim,
    seed = .seed, phi = phi, range = range,
    region_cutoff = region_cutoff,
    year_arima.sim = year_arima.sim, year_marginal_sd = year_marginal_sd,
    coefs = sim_coefs, sigma_E = sigma_E,
    heavy_sd_mult = heavy_sd_mult,
    heavy_sd_frac = heavy_sd_frac,
    svc_trend = svc_trend,
    north_effect = 0, sigma_O = sigma_O
  )
  sim_dat <- x$sim_dat
  sim_dat$fyear <- as.factor(sim_dat$year)
  predictor_dat <- x$predictor_dat
  predictor_dat$fyear <- as.factor(predictor_dat$year)

  # Visualize what we just did ----------------------------------------------

  # # Year effects:
  if (make_plots) {
    g <- ggplot(data.frame(x = seq_len(n_year), y = x$year_effects), aes(x, y)) +
      geom_line()
    print(g)
    g <- ggplot(data.frame(x = seq_len(n_year), y = x$year_effects), aes(x, exp(y))) +
      geom_line()
    print(g)

    blank_theme_elements <- theme(panel.grid.major = element_line(colour = "grey90"),
      # panel.spacing.x = unit(20, "pt"),
      axis.text = element_blank(), axis.ticks = element_blank(),
      axis.title = element_blank(), legend.position = "right")

    g <- ggplot(sim_dat, aes(X, Y, fill = eta)) +
      geom_raster() +
      facet_wrap(vars(year)) +
      scale_fill_viridis_c() +
      ggsidekick::theme_sleek() +
      coord_equal(expand = FALSE) +
      blank_theme_elements +
      labs(fill = "True\nsimulated\nlog abundance")
    print(g)
    # ggsave("stitch/figs/spatio-temporal-truth.pdf", width = 8, height = 5)
  }

  # Sample N per year -----------------------------------------------------

  cli::cli_alert_success("Observing...")
  d <- observe(
    sim_dat,
    sample_n = obs_sampled_size,
    region_cutoff = region_cutoff,
    seed = .seed,
    north_yrs = obs_yrs$north_yrs,
    sample_before_split = sample_before_split,
    south_yrs = obs_yrs$south_yrs, gap = gap_size
  )

  # Visualize it ------------------------------------------------------------

  if (make_plots) {
    g <- ggplot(d, aes(X, Y, colour = log(observed))) +
      geom_point() +
      facet_wrap(vars(year)) +
      scale_colour_viridis_c() +
      ggsidekick::theme_sleek() +
      blank_theme_elements +
      coord_equal(expand = FALSE) +
      theme(panel.grid.major = element_line(colour = "grey90"))
    print(g)
  }
  # # panel.spacing.x = unit(20, "pt")
  # # ggsave("stitch/figs/spatio-temporal-observed.pdf", width = 8, height = 5)

  # Calculate known true biomass/abundance ----------------------------------

  regions_samp <- select(d, year, sampled_region) %>% distinct()

  actual <- group_by(sim_dat, year) %>%
    summarise(total = sum(mu)) |>
    left_join(regions_samp, by = "year")

  if (make_plots) {
    g <- ggplot(actual, aes(year, total)) +
      geom_line()
    print(g)
  }

  if (make_plots) return(NULL)

  # Fit models --------------------------------------------------------------

  cli::cli_alert_success("Fitting models...")
  mesh <- make_mesh(d, c("X", "Y"), cutoff = 0.1)
  priors <- sdmTMBpriors(
    matern_s = pc_matern(range_gt = 0.2, sigma_lt = 1.5),
    matern_st = pc_matern(range_gt = 0.2, sigma_lt = 0.5)
  )
  # priors <- sdmTMBpriors()

  fits <- list()
  nms <- c()
  i <- 1
  ctl <- sdmTMBcontrol(nlminb_loops = 1L, newton_loops = 1L)

  check_sanity <- function(x) {
    if (!all(unlist(sanity(x, gradient_thresh = 0.01)))) {
      return(NA)
    } else {
      return(x)
    }
  }

  cli::cli_inform("Fitting st = 'rw'")
  fits[[i]] <- sdmTMB(
    observed ~ 1,
    family = tweedie(),
    data = d, time = "year", spatiotemporal = "rw", spatial = "on",
    silent = TRUE, mesh = mesh,
    priors = priors,
    control = ctl
  )
  fits[[i]] <- check_sanity(fits[[i]])
  nms <- c(nms, "RW")
  i <- i + 1

  cli::cli_inform("Fitting st IID covariate")
  fits[[i]] <- sdmTMB(
    observed ~ 0 + as.factor(year) + depth_cov + I(depth_cov^2),
    family = tweedie(),
    data = d, time = "year", spatiotemporal = "iid", spatial = "on",
    silent = TRUE, mesh = mesh,
    priors = priors,
    control = ctl
  )
  fits[[i]] <- check_sanity(fits[[i]])
  nms <- c(nms, "IID covariate")
  i <- i + 1

  cli::cli_inform("Fitting st IID s(year)")
  fits[[i]] <- sdmTMB(
    observed ~ s(year),
    family = tweedie(),
    data = d, time = "year", spatiotemporal = "iid", spatial = "on",
    silent = TRUE, mesh = mesh,
    priors = priors,
    control = ctl
  )
  fits[[i]] <- check_sanity(fits[[i]])
  nms <- c(nms, "IID s(year)")
  i <- i + 1

  cli::cli_inform("Fitting st IID no covariate as.factor year")
  fits[[i]] <- sdmTMB(
    observed ~ 0 + as.factor(year),
    family = tweedie(),
    data = d, time = "year", spatiotemporal = "iid", spatial = "on",
    mesh = mesh,
    priors = priors,
    control = ctl
  )
  fits[[i]] <- check_sanity(fits[[i]])
  i <- i + 1
  nms <- c(nms, "IID")

  cli::cli_inform("Fitting st time_varying RW")
  fits[[i]] <- sdmTMB(
    observed ~ 0,
    family = tweedie(),
    time_varying = ~1,
    data = d, time = "year", spatiotemporal = "iid", spatial = "on",
    mesh = mesh,
    priors = priors,
    control = ctl
  )
  fits[[i]] <- check_sanity(fits[[i]])
  i <- i + 1
  nms <- c(nms, "IID RW year")

  cli::cli_inform("Fitting st (1|year)")
  fits[[i]] <- sdmTMB(
    observed ~ 1 + (1 | fyear),
    family = tweedie(),
    data = d, time = "year", spatiotemporal = "iid", spatial = "on",
    mesh = mesh,
    priors = priors,
    control = ctl
  )
  fits[[i]] <- check_sanity(fits[[i]])
  i <- i + 1
  nms <- c(nms, "IID (1|year)")

  # cli::cli_inform("Fitting st time_varying AR1")
  # fits[[i]] <- sdmTMB(
  #   observed ~ 0,
  #   family = tweedie(),
  #   time_varying = ~1, time_varying_type = "ar1",
  #   data = d, time = "year", spatiotemporal = "iid", spatial = "on",
  #   mesh = mesh,
  #   priors = priors,
  #   control = ctl
  # )
  # fits[[i]] <- check_sanity(fits[[i]])
  # i <- i + 1
  # nms <- c(nms, "IID AR1 year")

  cli::cli_inform("Fitting spatial only")
  fits[[i]] <- sdmTMB(
    observed ~ 0 + as.factor(year),
    family = tweedie(),
    data = d, time = "year", spatiotemporal = "off", spatial = "on",
    mesh = mesh,
    priors = priors,
    control = ctl
  )
  fits[[i]] <- check_sanity(fits[[i]])
  i <- i + 1
  nms <- c(nms, "Spatial only")

  cli::cli_inform("Fitting SVC trend model")

  mean_year <- mean(d$year)
  d$year_cent <- d$year - mean_year
  fits[[i]] <- sdmTMB(
    observed ~ 0 + as.factor(year),
    family = tweedie(),
    data = d, time = "year",
    spatiotemporal = "iid",
    spatial = "on",
    spatial_varying = ~ 0 + year_cent,
    mesh = mesh,
    priors = priors,
    control = ctl
  )
  fits[[i]] <- check_sanity(fits[[i]])
  i <- i + 1
  nms <- c(nms, "SVC trend, IID fields")

  cli::cli_inform("Fitting SVC trend model without ST fields")

  fits[[i]] <- sdmTMB(
    observed ~ 0 + as.factor(year),
    family = tweedie(),
    data = d, time = "year",
    spatiotemporal = "off",
    spatial = "on",
    spatial_varying = ~ 0 + year_cent,
    mesh = mesh,
    priors = priors,
    control = ctl
  )
  fits[[i]] <- check_sanity(fits[[i]])
  i <- i + 1
  nms <- c(nms, "SVC trend, spatial only")

  names(fits) <- nms

  # Predict on grid and calculate indexes -----------------------------------

  cli::cli_alert_success("Calculating indices...")
  nd <- select(predictor_dat, X, Y, year, fyear, region)
  nd$depth_cov <- nd$Y
  nd$year_cent <- nd$year - mean_year

  nd

  preds <- purrr::map(fits, function(.x) {
    if (inherits(.x, "sdmTMB")) {
      out <- predict(.x, newdata = nd, return_tmb_object = TRUE)
    } else {
      out <- NA
    }
    out
  })

  indexes <- purrr::map(preds, function(.x) {
    if (length(.x) > 1) {
      get_index(.x, bias_correct = TRUE)
    }
  })

  indexes <- indexes[!sapply(indexes, is.null)]

  indexes_df <- dplyr::bind_rows(indexes, .id = "model") |>
    mutate(with_depth = paste0("covariate = ", grepl("covariate", model))) |>
    mutate(type = gsub(" covariate", "", model))

  indexes_df <- left_join(indexes_df, actual, by = "year") |>
    mutate(seed = .seed)
  indexes_df
}

