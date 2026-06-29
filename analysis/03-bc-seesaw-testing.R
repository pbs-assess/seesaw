# 2026-06-26
# Running across 21 synoptic groundfish
# Could do for HBLL OUT too
# And maybe HBLL inside??
# Assess the seesawness for the various models
# Could probably also look at correlates in this case (e.g., estimate things like north-south gradient or difference in mean density by subregion)
#
# Then do this for the Norwegian survey?
# Can summarize the overall findings and show some case study panels in paper
#
# This could be more informative than the simulation study on which models work... although only that can get at the coverage and RMSE

library(sdmTMB)
library(ggplot2)
theme_set(theme_light())
library(dplyr)

surveyjoin::cache_data()
surveyjoin::load_sql_data()

# dat_test <- surveyjoin::get_data(regions = "pbs")
# spp_to_fit <- dat_test |>
#   tidyr::drop_na(effort, catch_weight, depth_m) |>
#   group_by(common_name) |>
#   summarise(prop_positive = mean(catch_weight > 0)) |>
#   filter(prop_positive > 0.25) |>
#   pull(common_name) |> sort()
# spp_to_fit
# dput(spp_to_fit)

do_fit <- function(.sp, .survey) {
  RhpcBLASctl::blas_set_num_threads(1L)
  RhpcBLASctl::omp_set_num_threads(1L)

  safe_sanity <- function(fit) {
    if (!inherits(fit, "sdmTMB")) {
      return(FALSE)
    }
    ok <- tryCatch(
      all(unlist(sanity(fit, gradient_thresh = 0.01))),
      error = function(e) {
        message(.sp, " - sanity check errored: ", conditionMessage(e))
        FALSE
      }
    )
    if (!ok) {
      message(.sp, " - failed sanity checks.")
    }
    ok
  }

  safe_fit <- function(expr) {
    fit <- tryCatch(
      expr,
      error = function(e) {
        message(.sp, " - fit errored: ", conditionMessage(e))
        NULL
      }
    )
    if (safe_sanity(fit)) {
      fit
    } else {
      NULL
    }
  }

  safe_index <- function(fit) {
    tryCatch(
      {
        pred <- predict(fit, newdata = grid, return_tmb_obj = TRUE)
        get_index(pred)
      },
      error = function(e) {
        message(.sp, " - index errored: ", conditionMessage(e))
        NULL
      }
    )
  }

  if (.survey == "synoptic") {
    dat0 <- surveyjoin::get_data(.sp, regions = "pbs") |>
      mutate(year = lubridate::year(lubridate::ymd(date))) |>
      select(survey_name, year, lon_start, lat_start, depth_m, effort, catch_weight, common_name)

    dat0 <- sdmTMB::add_utm_columns(dat0, ll_names = c("lon_start", "lat_start"), utm_crs = 3156)
    dat <- dat0 |>
      # Use only complete N/S sampling years
      filter(!(year %in% c(2003, 2004, 2020))) |>
      tidyr::drop_na(effort, catch_weight, depth_m) |>
      # Drop these surveys to be perfectly bienniel
      filter(!(year == 2007 & survey_name == "SYN WCHG")) |>
      filter(!(year == 2021 & survey_name == "SYN WCVI"))

    dat_all <- dat0 |>
      tidyr::drop_na(effort, catch_weight, depth_m)

    grid <- surveyjoin::dfo_synoptic_grid |>
      sdmTMB::replicate_df("year", unique(dat$year))
    grid <- sdmTMB::add_utm_columns(grid, c("lon", "lat"), utm_crs = 3156)

    mesh <- make_mesh(dat, c("X", "Y"), cutoff = 10)
    mesh_all <- make_mesh(dat_all, c("X", "Y"), mesh = mesh$mesh)
  }
  if (.survey == "hbll") {
    stop("Not implemented.")
  }

  all_yrs <- seq(min(dat$year), max(dat$year))

  fits <- list()

  base_model <- "IID RF, factor(year)"
  fits[[base_model]] <- safe_fit(sdmTMB(
    catch_weight ~ 0 + factor(year),
    data = dat,
    mesh = mesh,
    offset = log(dat$effort),
    family = delta_gamma(type = "poisson-link"),
    time = "year",
    spatial = "on",
    spatiotemporal = "off",
    share_range = TRUE,
    anisotropy = TRUE,
    silent = FALSE
  ))

  if (is.null(fits[[base_model]])) {
    return(dplyr::tibble())
  }

  fits[["IID RF, factor(year), all data"]] <- safe_fit(update(
    fits[[base_model]],
    data = dat_all,
    mesh = mesh_all,
    offset = log(dat_all$effort)
  ))

  fits[["IID RF, factor(year), factor(survey), all data"]] <- safe_fit(update(
    fits[[base_model]],
    formula. = . ~ factor(year) + factor(survey_name),
    data = dat_all,
    mesh = mesh_all,
    offset = log(dat_all$effort)
  ))

  # fits[["RW RF"]] <- safe_fit(update(
  #   fits[[base_model]],
  #   formula. = . ~ 1,
  #   spatiotemporal = "rw",
  #   extra_time = all_yrs
  # ))
  #
  # fits[["AR1 RF"]] <- safe_fit(update(
  #   fits[[base_model]],
  #   formula. = . ~ 1,
  #   spatiotemporal = "ar1",
  #   extra_time = all_yrs
  # ))
  #
  # fits[["RW RF, RW year"]] <- safe_fit(update(
  #   fits[[base_model]],
  #   formula. = . ~ 1,
  #   spatiotemporal = "rw",
  #   time_varying = ~1,
  #   time_varying_type = "rw0",
  #   priors = sdmTMB::sdmTMBpriors(sigma_V = sdmTMB::gamma_cv(0.3, 0.5)),
  #   extra_time = all_yrs
  # ))
  #
  # fits[["AR1 RF, RW year"]] <- safe_fit(update(
  #   fits[[base_model]],
  #   formula. = . ~ 1,
  #   spatiotemporal = "ar1",
  #   time_varying = ~1,
  #   time_varying_type = "rw0",
  #   priors = sdmTMB::sdmTMBpriors(sigma_V = sdmTMB::gamma_cv(0.3, 0.5)),
  #   extra_time = all_yrs
  # ))
  #
  # fits[["IID RF, RW year"]] <- safe_fit(update(
  #   fits[[base_model]],
  #   formula. = . ~ 1,
  #   spatiotemporal = "iid",
  #   time_varying = ~1,
  #   time_varying_type = "rw0",
  #   priors = sdmTMB::sdmTMBpriors(sigma_V = sdmTMB::gamma_cv(0.3, 0.5)),
  #   extra_time = all_yrs
  # ))
  #
  # fits[["Spatial only, RW year"]] <- safe_fit(update(
  #   fits[[base_model]],
  #   formula. = . ~ 1,
  #   spatiotemporal = "off",
  #   time_varying = ~1,
  #   time_varying_type = "rw0",
  #   priors = sdmTMB::sdmTMBpriors(sigma_V = sdmTMB::gamma_cv(0.3, 0.5)),
  #   extra_time = all_yrs
  # ))
  #
  # fits[["IID RF, factor(year), depth"]] <- safe_fit(update(
  #   fits[[base_model]],
  #   formula. = . ~ 1,
  #   spatiotemporal = "off",
  #   time_varying = ~1,
  #   time_varying_type = "rw0"
  # ))

  indexes <- purrr::map(fits, safe_index) |>
    dplyr::bind_rows(.id = "model")

  indexes$species <- .sp
  indexes
}

spp_to_fit_syn <- c(
  "arrowtooth flounder",
  "dover sole",
  "english sole",
  "flathead sole",
  # "greenstriped rockfish",
  "lingcod",
  "longnose skate",
  # "north pacific hake",
  "pacific cod",
  "pacific halibut",
  "pacific ocean perch",
  "pacific spiny dogfish",
  "petrale sole",
  "redbanded rockfish",
  "rex sole",
  "sablefish",
  # "sharpchin rockfish",
  "shortspine thornyhead",
  "slender sole",
  "spotted ratfish",
  "walleye pollock"
)

# from ICES JMS MPA paper
spp_to_fit_hbll <- c(
  "rougheye/blackspotted rockfish",
  "china rockfish",
  "copper rockfish",
  "redbanded rockfish",
  "north pacific spiny dogfish",
  "tiger rockfish",
  "lingcod",
  "canary rockfish",
  "quillback rockfish",
  "shortspine thornyhead",
  "yelloweye rockfish",
  "silvergray rockfish",
  "spotted rockfish",
  "big skate",
  "rosethorn rockfish",
  "souther rock sole",
  "longnose skate",
  "pacific cod",
  "arrowtooth flounder"
)

RhpcBLASctl::blas_set_num_threads(1L)
RhpcBLASctl::omp_set_num_threads(1L)

future::plan(future::multicore, workers = min(c(length(spp_to_fit), future::availableCores())))
out <- purrr::map_dfr(spp_to_fit, do_fit, .survey = "synoptic")
out <- furrr::future_map_dfr(spp_to_fit_syn, do_fit, .survey = "synoptic")
# out <- furrr::future_map_dfr(spp_to_fit_syn, do_fit, .survey = "hbll")
dir.create("data-generated")
saveRDS(out, file = "data-generated/bc-indexes2.rds")
#
# out <- readRDS("data-generated/bc-indexes.rds")
# out2 <- readRDS("data-generated/bc-indexes2.rds")
#
# lu <- data.frame(year = sort(unique(dat$year)))
# lu$even <- lu$year %% 2 == 0
# lu$survey_group <- ifelse(lu$even, "WCHG + WCVI", "QCS + HS")
#
# moving_window_acf <- function(x, window = 8L) {
#   n <- length(x)
#   if (n < window) {
#     return(NA_real_)
#   }
#   acf_vals <- vapply(seq_len(n - window + 1L), function(i) {
#     acf(x[i:(i + window - 1L)], plot = FALSE)$acf[2L]
#   }, numeric(1L))
#   min(acf_vals, na.rm = TRUE)
# }
#
# x <- out |>
#   group_by(species, model) |>
#   summarise(seesaw_index = acf(log_est, plot = FALSE)$acf[2])
#
# out |>
#   left_join(lu) |>
#   left_join(x) |>
#   group_by(species, model) |>
#   mutate(geomean = exp(mean(log(est))), est = est / geomean, lwr = lwr / geomean, upr = upr / geomean) |>
#   ggplot(aes(year, log(est), ymin = log(lwr), ymax = log(upr))) +
#   geom_ribbon(fill = "grey90") +
#   geom_linerange(aes(colour = survey_group)) +
#   geom_point(aes(colour = survey_group), size = 2) +
#   scale_colour_brewer(palette = "Dark2") +
#   facet_grid(forcats::fct_reorder(model, seesaw_index) ~ species) +
#   ylab("Biomass index") +
#   xlab("Year") +
#   labs(colour = "Survey\ngrouping") +
#   ggsidekick::theme_sleek()
# ggsave("figs/bc-testing.pdf", width = 30, height = 15)
#
# out |>
#   left_join(lu) |>
#   group_by(model) |>
#   summarise(n = n())
#
# out |>
#   left_join(lu) |>
#   group_by(species, model) |>
#   # summarise(seesaw_index = abs(mean(log_est[which(even)]) - mean(log_est[which(!even)]))) |>
#   summarise(seesaw_index = acf(log_est, plot = FALSE)$acf[2]) |>
#   group_by(model) |>
#   mutate(mean_acf = mean(seesaw_index)) |>
#   ggplot(aes(forcats::fct_reorder(model, mean_acf), seesaw_index)) +
#   geom_point(position = position_jitter(width = 0)) +
#   coord_flip() +
#   ylab("First-order index autocorrelation") +
#   scale_colour_gradient2(mid = "grey80") +
#   ggsidekick::theme_sleek() +
#   theme(axis.title.y = element_blank(), panel.grid.major = element_line(colour = "grey90", linewidth = 0.3), panel.grid.minor = element_line(colour = "grey90", linewidth = 0.3))
#
# ggsave("figs/bc-trawl-acf.pdf", width = 4.5, height = 3.5)
#
# out |>
#   left_join(lu) |>
#   group_by(species, model) |>
#   summarise(seesaw_index = moving_window_acf(log_est, window = 10)) |>
#   group_by(model) |>
#   mutate(mean_acf = mean(seesaw_index)) |>
#   ggplot(aes(forcats::fct_reorder(model, mean_acf), seesaw_index)) +
#   geom_point(position = position_jitter(width = 0)) +
#   coord_flip() +
#   ylab("Most negative lag-1 autocorrelation\nacross 10-year windows") +
#   scale_colour_gradient2(mid = "grey80") +
#   ggsidekick::theme_sleek() +
#   theme(axis.title.y = element_blank(), panel.grid.major = element_line(colour = "grey90", linewidth = 0.3), panel.grid.minor = element_line(colour = "grey90", linewidth = 0.3))
#
# ggsave("figs/bc-trawl-acf-moving-window.pdf", width = 4.5, height = 3.5)
