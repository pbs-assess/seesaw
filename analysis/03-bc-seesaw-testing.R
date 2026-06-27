library(sdmTMB)
library(ggplot2)
theme_set(theme_light())
library(dplyr)

# surveyjoin::cache_data()
surveyjoin::load_sql_data()

Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1"
)

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
  dat <- surveyjoin::get_data(.sp, regions = "pbs") |>
    mutate(year = lubridate::year(lubridate::ymd(date))) |>
    select(survey_name, year, lon_start, lat_start, depth_m, effort, catch_weight, common_name)

  dat <- sdmTMB::add_utm_columns(dat, ll_names = c("lon_start", "lat_start"), utm_crs = 3156)
  dat <- dat |>
    # Use only complete N/S sampling years
    filter(!(year %in% c(2003, 2004, 2020))) |>
    tidyr::drop_na(effort, catch_weight, depth_m) |>
    # Drop these surveys to be perfectly bienniel
    filter(!(year == 2007 & survey_name == "SYN WCHG")) |>
    filter(!(year == 2021 & survey_name == "SYN WCVI"))
  table(dat$survey_name, dat$year)

  grid <- surveyjoin::dfo_synoptic_grid |>
    sdmTMB::replicate_df("year", unique(dat$year))
  grid <- sdmTMB::add_utm_columns(grid, c("lon", "lat"), utm_crs = 3156)

  map_data <- rnaturalearth::ne_countries(
    scale = "large",
    returnclass = "sf", country = "canada")
  bc_coast <- suppressWarnings(suppressMessages(
    sf::st_crop(map_data,
      c(xmin = -134, ymin = 46, xmax = -120, ymax = 57))))
  bc_coast <- sf::st_transform(bc_coast, crs = 3156)

  mesh <- make_mesh(dat, c("X", "Y"), cutoff = 15)

  fits <- list()

  fits[[1]] <- sdmTMB(
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
  )
  names(fits)[[1]] <- "IID RF, factor(year)"

  fits[[2]] <- update(
    fits[[1]],
    formula. = . ~ 1,
    spatiotemporal = "rw"
  )
  names(fits)[[2]] <- "RW RF"

  fits[[3]] <- update(
    fits[[1]],
    formula. = . ~ 1,
    spatiotemporal = "ar1"
  )
  names(fits)[[3]] <- "AR1 RF"

  fits[[4]] <- update(
    fits[[1]],
    formula. = . ~ 1,
    spatiotemporal = "rw",
    time_varying = ~ 1,
    time_varying_type = "rw0"
  )
  names(fits)[[4]] <- "RW RF, RW year"

  fits[[5]] <- update(
    fits[[1]],
    formula. = . ~ 1,
    spatiotemporal = "ar1",
    time_varying = ~ 1,
    time_varying_type = "rw0"
  )
  names(fits)[[5]] <- "AR1 RF, RW year"

  fits[[6]] <- update(
    fits[[1]],
    formula. = . ~ 1,
    spatiotemporal = "iid",
    time_varying = ~ 1,
    time_varying_type = "rw0"
  )
  names(fits)[[6]] <- "IID RF, RW year"

  fits[[7]] <- update(
    fits[[1]],
    formula. = . ~ 1,
    spatiotemporal = "off",
    time_varying = ~ 1,
    time_varying_type = "rw0"
  )
  names(fits)[[7]] <- "Spatial only, RW year"

  indexes <- purrr::map_dfr(fits, \(x) {
    pred <- predict(x, newdata = grid, return_tmb_obj = TRUE)
    get_index(pred)
  }, .id = "model")

  indexes$species <- .sp
}


spp_to_fit <- c("arrowtooth flounder", "dover sole", "english sole", "flathead sole",
  "greenstriped rockfish", "lingcod", "longnose skate", "north pacific hake",
  "pacific cod", "pacific halibut", "pacific ocean perch", "pacific spiny dogfish",
  "petrale sole", "redbanded rockfish", "rex sole", "sablefish",
  "sharpchin rockfish", "shortspine thornyhead", "slender sole",
  "spotted ratfish", "walleye pollock")

future::plan(future::multicore, workers = min(c(length(spp_to_fit), future::availableCores())))
out <- furrr::future_map_dfr(spp_to_fit, do_fit)
dir.create("data-generated")
saveRDS(out, file = "data-generated/bc-indexes.rds")

# lu <- data.frame(year = sort(unique(dat$year)))
# lu$even <- lu$year %% 2 == 0
# lu$survey_group <- ifelse(lu$even, "WCHG + WCVI", "QCS + HS")



# indexes |>
#   left_join(lu) |>
#   ggplot(aes(year, est, ymin = lwr, ymax = upr)) +
#   geom_ribbon(fill = "grey90") +
#   geom_linerange(aes(colour = survey_group)) +
#   geom_point(aes(colour = survey_group), size = 2) +
#   scale_colour_brewer(palette = "Dark2") +
#   facet_wrap(~model) +
#   ylab("Biomass index") + xlab("Year") + labs(colour = "Survey\ngrouping")
