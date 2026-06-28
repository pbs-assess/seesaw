# 2026-06-26
# Running across 21 synoptic groundfish
# Could do for HBLL OUT too
# And maybe HBLL inside??
# Assess the seesawness for the various models
# Maybe add back a version that has depth included
# And add back a version that has all the bits and pieces (not perfectly biennial) added back too?
# And try a version with a factor for region in that case?
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

do_fit <- function(.sp) {
  RhpcBLASctl::blas_set_num_threads(1L)
  RhpcBLASctl::omp_set_num_threads(1L)

  safe_sanity <- function(fit, model_name) {
    if (!inherits(fit, "sdmTMB")) {
      return(FALSE)
    }
    ok <- tryCatch(
      all(unlist(sanity(fit, gradient_thresh = 0.01))),
      error = function(e) {
        message(.sp, " - ", model_name, " sanity check errored: ", conditionMessage(e))
        FALSE
      }
    )
    if (!ok) {
      message(.sp, " - ", model_name, " failed sanity checks.")
    }
    ok
  }

  safe_fit <- function(model_name, expr) {
    message(.sp, " - fitting ", model_name)
    fit <- tryCatch(
      expr,
      error = function(e) {
        message(.sp, " - ", model_name, " fit errored: ", conditionMessage(e))
        NULL
      }
    )
    if (safe_sanity(fit, model_name)) {
      fit
    } else {
      NULL
    }
  }

  safe_index <- function(fit, model_name) {
    tryCatch({
      pred <- predict(fit, newdata = grid, return_tmb_obj = TRUE)
      get_index(pred)
    }, error = function(e) {
      message(.sp, " - ", model_name, " index errored: ", conditionMessage(e))
      NULL
    })
  }

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

  fits <- list()

  base_model <- "IID RF, factor(year)"
  fits[[base_model]] <- safe_fit(base_model, sdmTMB(
    catch_weight ~ 0 + factor(year),
    data = dat,
    mesh = mesh,
    offset = log(dat$effort),
    family = delta_gamma(type = "poisson-link"),
    time = "year",
    spatial = "on",
    spatiotemporal = "iid",
    share_range = TRUE,
    anisotropy = TRUE,
    silent = FALSE
  ))

  if (is.null(fits[[base_model]])) {
    message(.sp, " - skipping species because the base model did not fit.")
    return(dplyr::tibble())
  }

  fits[["IID RF, factor(year), all data"]] <- safe_fit("IID RF, factor(year), all data", update(
    fits[[base_model]],
    data = data_all,
    mesh = mesh_all
  ))

  fits[["IID RF, factor(year), factor(survey), all data"]] <- safe_fit("IID RF, factor(year), factor(survey), all data", update(
    fits[[base_model]],
    formula. = . ~ factor(year) + factor(survey_name),
    data = data_all,
    mesh = mesh_all
  ))

  fits[["RW RF"]] <- safe_fit("RW RF", update(
    fits[[base_model]],
    formula. = . ~ 1,
    spatiotemporal = "rw"
  ))

  fits[["AR1 RF"]] <- safe_fit("AR1 RF", update(
    fits[[base_model]],
    formula. = . ~ 1,
    spatiotemporal = "ar1"
  ))

  fits[["RW RF, RW year"]] <- safe_fit("RW RF, RW year", update(
    fits[[base_model]],
    formula. = . ~ 1,
    spatiotemporal = "rw",
    time_varying = ~ 1,
    time_varying_type = "rw0"
  ))

  fits[["AR1 RF, RW year"]] <- safe_fit("AR1 RF, RW year", update(
    fits[[base_model]],
    formula. = . ~ 1,
    spatiotemporal = "ar1",
    time_varying = ~ 1,
    time_varying_type = "rw0"
  ))

  fits[["IID RF, RW year"]] <- safe_fit("IID RF, RW year", update(
    fits[[base_model]],
    formula. = . ~ 1,
    spatiotemporal = "iid",
    time_varying = ~ 1,
    time_varying_type = "rw0"
  ))

  fits[["Spatial only, RW year"]] <- safe_fit("Spatial only, RW year", update(
    fits[[base_model]],
    formula. = . ~ 1,
    spatiotemporal = "off",
    time_varying = ~ 1,
    time_varying_type = "rw0"
  ))

  fits[["IID RF, factor(year), depth"]] <- safe_fit("Spatial only, RW year", update(
    fits[[base_model]],
    formula. = . ~ 1,
    spatiotemporal = "off",
    time_varying = ~ 1,
    time_varying_type = "rw0"
  ))

  indexes <- purrr::imap(fits, \(x, nm) safe_index(x, nm)) |>
    dplyr::bind_rows(.id = "model")

  indexes$species <- .sp
  indexes
}


spp_to_fit <- c(
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

RhpcBLASctl::blas_set_num_threads(1L)
RhpcBLASctl::omp_set_num_threads(1L)

future::plan(future::multicore, workers = min(c(length(spp_to_fit), future::availableCores())))
out <- furrr::future_map_dfr(spp_to_fit, do_fit)
dir.create("data-generated")
saveRDS(out, file = "data-generated/bc-indexes.rds")

out <- readRDS("data-generated/bc-indexes.rds")

# lu <- data.frame(year = sort(unique(dat$year)))
# lu$even <- lu$year %% 2 == 0
# lu$survey_group <- ifelse(lu$even, "WCHG + WCVI", "QCS + HS")

out |>
  left_join(lu) |>
  filter(species != "sharpchin rockfish") |>
  group_by(species, model) |>
  mutate(geomean = exp(mean(log(est))), est = est / geomean, lwr = lwr / geomean, upr = upr/ geomean) |>
  ggplot(aes(year, log(est), ymin = log(lwr), ymax = log(upr))) +
  geom_ribbon(fill = "grey90") +
  geom_linerange(aes(colour = survey_group)) +
  geom_point(aes(colour = survey_group), size = 2) +
  scale_colour_brewer(palette = "Dark2") +
  facet_grid(model~species) +
  ylab("Biomass index") + xlab("Year") + labs(colour = "Survey\ngrouping") +
  ggsidekick::theme_sleek()

out |>
  left_join(lu) |>
  group_by(model) |>
  summarise(n = n())

out |>
  left_join(lu) |>
  group_by(species, model) |>
  summarise(seesaw_index = abs(mean(log_est[which(even)]) - mean(log_est[which(!even)]))) |>
  ggplot(aes(model, seesaw_index)) + geom_violin() + coord_flip() + scale_y_log10()

