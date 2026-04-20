sim_ar1_heavy <- function(
    n, marginal_sd, rho, mu = 0,
    heavy_sd_mult = 1, heavy_sd_frac = 0, s = 42) {
  set.seed(s)

  # Exact non-heavy path to match stats::arima.sim() behavior.
  if (heavy_sd_frac == 0) {
    x <- stats::arima.sim(
      model = list(ar = rho),
      n = n,
      sd = sqrt(1 - rho^2)
    )
    return(as.numeric(x) * marginal_sd + mu)
  }

  sigma <- sqrt(1 - rho^2)
  devs <- rnorm(n) * marginal_sd
  set.seed(s)
  devs_heavy <- rnorm(n) * marginal_sd * heavy_sd_mult
  hd <- rbinom(n, size = 1, prob = heavy_sd_frac)
  this_dev <- if (hd[1] == 0) devs[1] else devs_heavy[1]
  x0 <- this_dev
  B <- numeric(n)
  B[1] <- rho * x0 + sigma * this_dev
  for (i in seq(2, length(B))) {
    this_dev <- if (hd[i] == 0) devs[i] else devs_heavy[i]
    B[i] <- rho * B[i - 1] + sigma * this_dev
  }
  B + mu
}

sim <- function(
  predictor_grid, 
  mesh, 
  phi = 8,
  seed = 123,
  region_cutoff = 0.5,
  rho = 0,
  sigma_E = 0.4,
  range = 0.8,
  tweedie_p = 1.6,
  sigma_O = 1.6,
  coefs = c(2.2, 3.8),
  year_mean = 1,
  north_effect = 0,
  year_arima.sim = list(ar = 0.8),
  year_marginal_sd = 0.2,
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

  set.seed(seed * 1029)

  # x <- predictor_grid$depth_cov
  # plot(x, x*2 + (x^2)*5)
  #
  # x <- predictor_grid$depth_cov
  # lines(x, x*7, lty = 2)

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
    B = B
  )
  sim_dat$region <- predictor_grid$region
  sim_dat$year <- predictor_grid$year
  sim_dat$depth_cov <- predictor_grid$depth_cov
  list(sim_dat = sim_dat, predictor_dat = predictor_grid, year_effects = yrs, coefs = coefs, B = B)
}

resolve_sample_n <- function(sample_n, years) {
  years_chr <- as.character(sort(unique(years)))

  if (length(sample_n) == 1L) {
    vals <- rep(as.integer(sample_n), length(years_chr))
    names(vals) <- years_chr
    return(vals)
  }

  if (is.null(names(sample_n))) {
    if (length(sample_n) != length(years_chr)) {
      cli::cli_abort(c(
        "When `sample_n` has length > 1 and no names, it must match the number of years.",
        "x" = "{length(sample_n)} provided for {length(years_chr)} years."
      ))
    }
    vals <- as.integer(sample_n)
    names(vals) <- years_chr
    return(vals)
  }

  missing_years <- setdiff(years_chr, names(sample_n))
  if (length(missing_years) > 0L) {
    cli::cli_abort(c(
      "Named `sample_n` is missing one or more years.",
      "x" = "Missing year names: {paste(missing_years, collapse = ', ')}"
    ))
  }

  vals <- as.integer(sample_n[years_chr])
  names(vals) <- years_chr
  vals
}

sample_by_year <- function(d, sample_n, context = "sampling") {
  yrs <- sort(unique(d$year))
  n_by_year <- resolve_sample_n(sample_n, yrs)

  dplyr::group_by(d, year) |>
    dplyr::group_modify(\(.x, .y) {
      this_year <- as.character(.y$year[[1]])
      n_this <- n_by_year[[this_year]]
      if (nrow(.x) < n_this) {
        cli::cli_abort(c(
          "Not enough rows to sample in {.val {context}}.",
          "x" = "Year {.val {this_year}} has {nrow(.x)} rows but requested {n_this}."
        ))
      }
      dplyr::slice_sample(.x, n = n_this)
    }) |>
    dplyr::ungroup()
}

# sample_before_split = TRUE means apply sample_n first, then discard any gap
# smaple_before_split = FALSE applies sample_n *after* cutting gap
observe <- function(
    sim_dat, sample_n = 300L, seed = 10282, gap = 0, region_cutoff = 0.5,
    north_yrs = seq(1, 9, 2), south_yrs = seq(2, 10, 2),
    sample_before_split = FALSE) {
  assertthat::assert_that(gap >= 0 && gap <= 0.99)
  assertthat::assert_that(all(as.integer(sample_n) > 1L))
  assertthat::assert_that(all(north_yrs %in% sim_dat$year) &&
    all(south_yrs %in% sim_dat$year))
  assertthat::assert_that(region_cutoff > 0.05 && region_cutoff < 0.95)

  both_yrs <- intersect(north_yrs, south_yrs)

  set.seed(seed * 7)

  if (sample_before_split) {
    d <- sample_by_year(sim_dat, sample_n = sample_n, context = "pre-split sampling")
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
    d <- sample_by_year(d, sample_n = sample_n, context = "post-split sampling")
  }
  d$sampled_region <- as.character(d$region)
  d$sampled_region[d$year %in% both_yrs] <- "both"
  d
}

sim_fit_and_index <- function(
    n_year,
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
    save_plots = FALSE,
    return_fits = FALSE, 
    return_preds = FALSE, 
    return_sim_dat = FALSE,
    return_obs_dat = FALSE,
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

    blank_theme_elements <- theme(
      panel.grid.major = element_line(colour = "grey90"),
      # panel.spacing.x = unit(20, "pt"),
      axis.text = element_blank(), axis.ticks = element_blank(),
      axis.title = element_blank(), legend.position = "right"
    )

    g <- ggplot(sim_dat, aes(X, Y, fill = eta)) +
      geom_raster() +
      facet_wrap(vars(year)) +
      scale_fill_viridis_c() +
      ggsidekick::theme_sleek() +
      coord_equal(expand = FALSE) +
      blank_theme_elements +
      labs(fill = "True\nsimulated\nlog abundance")
    print(g)
    if (save_plots) ggsave("figs/spatio-temporal-truth.pdf", width = 8, height = 5)
    if (save_plots) ggsave("figs/spatio-temporal-truth.png", width = 8, height = 5)
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
  if (save_plots) ggsave("figs/spatio-temporal-observed.pdf", width = 8, height = 5)
  if (save_plots) ggsave("figs/spatio-temporal-observed.png", width = 8, height = 5)

  # Calculate known true biomass/abundance ----------------------------------

  regions_samp <- select(d, year, sampled_region) %>% distinct()

  if (return_sim_dat) {
    return(sim_dat)
  }
  if (return_obs_dat) {
    return(d)
  }
  actual <- group_by(sim_dat, year) %>%
    summarise(total = sum(mu)) |>
    left_join(regions_samp, by = "year")

  if (make_plots) {
    g <- ggplot(actual, aes(year, total)) +
      geom_line()
    print(g)
  }

  if (make_plots) {
    return(NULL)
  }

  # Fit models --------------------------------------------------------------

  cli::cli_alert_success("Fitting models...")
  mesh <- make_mesh(d, c("X", "Y"), cutoff = 0.1)
  priors <- sdmTMBpriors(
    matern_s = pc_matern(range_gt = 0.2, sigma_lt = 1.5),
    matern_st = pc_matern(range_gt = 0.2, sigma_lt = 0.5)
  )
  # priors <- sdmTMBpriors()

  ctl <- sdmTMBcontrol(nlminb_loops = 1L, newton_loops = 1L)

  check_sanity <- function(x) {
    if (!inherits(x, "sdmTMB")) {
      return(NA)
    }
    if (!all(unlist(sanity(x, gradient_thresh = 0.01)))) {
      return(NA)
    } else {
      return(x)
    }
  }
  lu <- data.frame(year = sort(unique(d$year)), year_pairs = as.factor(rep(seq(1, max(d$year)), each = 2)[seq_along(unique(d$year))]))
  d_pairs <- left_join(d, lu, by = join_by(year))

  source(here::here("analysis", "estimation-scenarios.R"))
  model_specs <- build_model_specs()

  data_lookup <- list(base = d, pairs = d_pairs)
  fit_one <- function(spec) {
    cli::cli_inform("Fitting {.val {spec$name}}")
    control_this <- if (is.null(spec$control)) ctl else spec$control
    data_key <- if (is.null(spec$data_key)) "base" else spec$data_key
    if (!data_key %in% names(data_lookup)) {
      cli::cli_abort(c(
        "Unknown data key in model spec {.val {spec$name}}.",
        "x" = "Expected one of: {paste(names(data_lookup), collapse = ', ')}.",
        "x" = "Got: {.val {data_key}}."
      ))
    }

    fit_args <- c(
      spec$fit_args,
      list(
        family = tweedie(),
        data = data_lookup[[data_key]],
        time = "year",
        spatial = "on",
        silent = TRUE,
        mesh = mesh,
        priors = priors,
        control = control_this
      )
    )
    fit <- tryCatch(do.call(sdmTMB, fit_args), error = function(e) NA)
    check_sanity(fit)
  }
  model_names <- purrr::map_chr(model_specs, "name")

  if (return_fits) {
    fits <- purrr::map(model_specs, fit_one)
    names(fits) <- model_names
    return(fits)
  }

  # Predict on grid and calculate indexes -----------------------------------

  cli::cli_alert_success("Calculating indices...")
  nd <- select(predictor_dat, X, Y, year, fyear, region)
  nd$depth_cov <- nd$Y
  # nd$year_cent <- nd$year - mean_year

  nd <- left_join(nd, lu, by = join_by(year))

  predict_one <- function(fit) {
    if (inherits(fit, "sdmTMB")) {
      predict(fit, newdata = nd, return_tmb_object = TRUE)
    } else {
      NA
    }
  }
  index_one <- function(pred) {
    if (length(pred) > 1) {
      get_index(pred, bias_correct = TRUE)
    }
  }

  if (return_preds) {
    preds <- setNames(vector("list", length(model_specs)), model_names)
    for (i in seq_along(model_specs)) {
      fit <- fit_one(model_specs[[i]])
      preds[[i]] <- predict_one(fit)
      rm(fit)
      if (i %% 4L == 0L) gc(FALSE)
    }
    return(preds)
  }

  # Default path: stream fit -> predict -> index per model to reduce memory.
  indexes <- list()
  for (i in seq_along(model_specs)) {
    fit <- fit_one(model_specs[[i]])
    pred <- predict_one(fit)
    idx <- index_one(pred)
    if (!is.null(idx)) {
      indexes[[model_names[[i]]]] <- idx
    }
    rm(fit, pred, idx)
    if (i %% 4L == 0L) gc(FALSE)
  }

  indexes_df <- dplyr::bind_rows(indexes, .id = "model") |>
    mutate(with_depth = paste0("covariate = ", grepl("covariate", model))) |>
    mutate(type = gsub(" covariate", "", model))

  indexes_df <- left_join(indexes_df, actual, by = "year") |>
    mutate(seed = .seed)
  indexes_df
}
