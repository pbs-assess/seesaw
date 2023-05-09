prep_data <- function(region = c("SYN QCS", "SYN HS", "SYN WCVI", "SYN WCHG")) {
  dat <- readRDS(here("data/all-survey-sets-2021-select.rds")) %>%
    dplyr::filter(survey_abbrev %in% region)
  # change from per m2 to per km2:
  dat$density <- dat[["density_kgpm2"]] * 1000000
  dat$log_depth <- log(dat$depth_m)
  dat$area_swept1 <- dat$doorspread_m * (dat$speed_mpm * dat$duration_min)
  dat$area_swept2 <- dat$tow_length_m * dat$doorspread_m
  dat$area_swept <- ifelse(!is.na(dat$area_swept2), dat$area_swept2, dat$area_swept1)
  dat <- dplyr::filter(dat, !is.na(area_swept))
  dat <- sdmTMB::add_utm_columns(dat, c("longitude", "latitude"), utm_crs = 32609)
  dat
}

read_species_list <- function(f) {
  x <- readr::read_delim(f, comment = "#", delim = "\n", col_types = "c")
  x$species
}

fit_index <- function(dat,
                      species = "arrowtooth flounder",
                      formula = catch_weight ~ s(log_depth, k = 5),
                      family = tweedie(),
                      cutoff = 20,
                      bias_correct = TRUE,
                      anisotropy = FALSE,
                      spatiotemporal = "rw",
                      silent = FALSE,
                      skip_sanity = FALSE,
                      priors = sdmTMBpriors(
                        matern_s = pc_matern(range_gt = 10, sigma_lt = 4),
                        matern_st = pc_matern(range_gt = 10, sigma_lt = 3)
                      ),
                      ...) {
  dat <- dplyr::filter(dat, species_common_name == tolower(species))
  region <- unique(dat$survey_abbrev)

  # bnd <- INLA::inla.nonconvex.hull(cbind(dat$X, dat$Y), convex = -0.04)
  # inla_mesh <- INLA::inla.mesh.2d(
  #   boundary = bnd,
  #   max.edge = c(50, 120),
  #   offset = c(-0.01, -0.02),
  #   cutoff = c(20, 100),
  #   min.angle = c(21, 21)
  # )
  # mesh <- make_mesh(dat, c("X", "Y"), mesh = inla_mesh)
  # plot(mesh$mesh, asp = 1)
  # points(dat$X, dat$Y, pch = ".")
  # mesh$mesh$n

  ctrl <- sdmTMBcontrol(newton_loops = 1L)

  mesh <- make_mesh(dat, c("X", "Y"), cutoff = 20)
  # plot(mesh$mesh, asp = 1)
  # points(dat$X, dat$Y, pch = ".")
  # mesh$mesh$n

  # yrs <- unique(dat$year)
  # all_yrs <- seq(min(dat$year), max(dat$year))
  # base::setdiff(yrs, all_yrs)

  nd <- readRDS(here("grids/synoptic_grid.rds")) %>%
    dplyr::filter(survey %in% region)
  fitted_yrs <- sort(unique(dat$year))
  nd <- make_grid(nd, years = fitted_yrs)
  nd <- na.omit(nd)
  nd$year <- as.integer(nd$year)
  nd$log_depth <- log(nd$depth)
  nd$cell_area <- 4e+6

  TMB::openmp(n = 2L, DLL = "sdmTMB")
  fit <- try(
    sdmTMB(
      formula,
      dat,
      mesh = mesh,
      time = "year",
      family = family,
      spatial = "on",
      spatiotemporal = spatiotemporal,
      offset = log(dat$area_swept / 100000),
      share_range = TRUE,
      anisotropy = anisotropy,
      silent = silent,
      control = ctrl,
      do_index = TRUE,
      priors = priors,
      index_args = list(area = nd$cell_area / 100000),
      predict_args = list(newdata = nd, return_tmb_object = TRUE),
      ...
    )
  )
  s <- sanity(fit)
  if ((s$all_ok || skip_sanity) && class(fit) != "try-error") {
    if (bias_correct) TMB::openmp(n = 1L, DLL = "sdmTMB")
    ind <- get_index(fit, bias_correct = bias_correct, area = nd$cell_area / 100000)
  } else {
    ind <- NA
  }

  list(sanity = s, fit = fit, index = ind)
}

make_grid <- function(.x, years) {
  years <- sort(unique(years))
  .nd <- do.call(
    "rbind",
    replicate(length(years), .x, simplify = FALSE)
  )
  .nd$year <- rep(years, each = nrow(.x))
  .nd
}
