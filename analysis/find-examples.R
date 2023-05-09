library(dplyr)
library(ggplot2)
library(sdmTMB)
source(here::here("stitch", "utils.R"))

mytheme <- function() ggsidekick::theme_sleek() +
  theme(panel.grid = element_line(colour = "grey85"),
  panel.grid.major = element_line(linewidth = rel(0.5)),
  panel.grid.minor = element_line(linewidth = rel(0.25))
)
theme_set(mytheme())

dat <-
  readRDS("data/all_surv_catch.rds") |>
  mutate(
    density_kgkm2 = density_kgpm2 * 1000000,
    log_depth = log(depth_m),
    area_swept1 = doorspread_m * (speed_mpm * duration_min),
    area_swept2 = tow_length_m * doorspread_m,
    area_swept = ifelse(!is.na(area_swept2), area_swept2, area_swept1)
  ) |>
  mutate(area_swept_km2 = area_swept / 1e6) |> # This should be used for offset???
  mutate(log_area_km2 = log(area_swept_km2)) |> # Value used for offset
  sdmTMB::add_utm_columns(c("longitude", "latitude"), utm_crs = 32609)

sp_dat <-
  filter(dat, stringr::str_detect(survey_abbrev, "SYN")) |>
  filter(species_common_name == "petrale sole") |>
  # simplify df columns
  select(
    survey_id, trip_id, fishing_event_id,
    year, month, day, latitude, longitude, X, Y,
    depth_m, log_depth,
    species_code, species_common_name,
    catch_weight, catch_count,
    density_kgkm2, density_pcpm2, density_ppkm2, survey_abbrev,
    # area_swept1, area_swept2, doorspread_m, speed_mpm, duration_min,
    # tow_length_m,
    area_swept, area_swept_km2, log_area_km2, hook_count, time_deployed
  ) |>
  # specify factor variables that will be used in models
  mutate(
    fyear = as.factor(year),
    region = as.factor(survey_abbrev)
  ) |>
  # use complete datasets
  tidyr::drop_na(area_swept) |> # drop empty area swept (no doorspread given)
  tidyr::drop_na(depth_m) # drop rows without depths

mesh <- make_mesh(sp_dat, c("X", "Y"), cutoff = 25)

# Inputs for predictions and index calcluations --------------------------------
syn_grid <-
  gfplot::synoptic_grid |>
  tibble() |> # because I accidentally print the full df too often
  dplyr::select(-survey_series_name, -utm_zone, -survey_domain_year)


SILENT <- FALSE

cli::cli_inform("Fitting st = 'rw'")
fit1 <-
  sdmTMB(
    catch_weight ~ 1,
    family = tweedie(),
    data = sp_dat, time = "year", spatiotemporal = "rw", spatial = "on",
    mesh = mesh,
    offset = "log_area_km2",
    silent = SILENT
  )
sanity(fit1)

cli::cli_inform("Fitting st IID covariate")
fit2 <- try(sdmTMB(
  catch_weight ~ 0 + as.factor(year) + log_depth + I(log_depth^2),
  family = tweedie(),
  data = sp_dat, time = "year", spatiotemporal = "iid", spatial = "on",
  silent = SILENT, mesh = mesh,
  offset = "log_area_km2"
))
sanity(fit2)

cli::cli_inform("Fitting st IID no covariate as.factor year")
fit4 <- sdmTMB(
  catch_weight ~ 0 + as.factor(year),
  family = tweedie(),
  data = sp_dat, time = "year", spatiotemporal = "iid", spatial = "on",
  mesh = mesh,
  offset = "log_area_km2",
  silent = SILENT
)
sanity(fit4)

cli::cli_inform("Fitting st time_varying RW")
fit5 <- sdmTMB(
  catch_weight ~ 0,
  family = tweedie(),
  time_varying = ~ 1, time_varying_type = "rw",
  data = sp_dat, time = "year", spatiotemporal = "iid", spatial = "on",
  mesh = mesh,
  offset = "log_area_km2",
  silent = SILENT
)
sanity(fit5)

cli::cli_inform("Fitting RW fields, as.factor(year)")
fit9 <- sdmTMB(
  catch_weight ~ 0 + as.factor(year),
  family = tweedie(),
  data = sp_dat,
  time = "year",
  spatial = "on",
  spatiotemporal = "rw",
  mesh = mesh,
  offset = "log_area_km2",
  silent = SILENT
)
sanity(fit9)

cli::cli_inform("Fitting spatial only, SVC time spline")

sp_dat$present <- as.integer(sp_dat$catch_weight > 0)

m <- mgcv::gam(present ~ s(year, k = 4), data = sp_dat)
mm <- model.matrix(m)
head(mm)
plot(mm[,1]) # intercept
plot(mm[,2])
plot(mm[,3])
plot(mm[,4])
sp_dat$spline1 <- as.numeric(mm[,2])
sp_dat$spline2 <- as.numeric(mm[,3])
sp_dat$spline3 <- as.numeric(mm[,4])
spline_lu <- sp_dat |>
  select(year, spline1, spline2, spline3) |>
  distinct()

fit10 <- sdmTMB(
  catch_weight ~ 1 + spline1 + spline2 + spline3,
  family = tweedie(),
  data = sp_dat,
  time = "year",
  spatial = "on",
  spatiotemporal = "off",
  spatial_varying = ~ 0 + spline1 + spline2 + spline3,
  time_varying = ~ 1,
  time_varying_type = "rw0",
  mesh = mesh,
  offset = "log_area_km2",
  silent = SILENT,
  do_fit = TRUE
)
sanity(fit10)
fit10

cli::cli_inform("Fitting IID fields, SVC time spline, RW year")

fit11 <- sdmTMB(
  catch_weight ~ 1 + spline1 + spline2 + spline3,
  family = tweedie(),
  data = sp_dat,
  time = "year",
  spatial = "on",
  spatiotemporal = "iid",
  spatial_varying = ~ 0 + spline1 + spline2 + spline3,
  time_varying = ~ 1,
  time_varying_type = "rw0",
  mesh = mesh,
  offset = "log_area_km2",
  silent = SILENT
)
sanity(fit11)
fit11

fit12 <- sdmTMB(
  catch_weight ~ as.factor(year) + as.factor(survey_abbrev),
  family = tweedie(),
  data = sp_dat,
  time = "year",
  spatial = "on",
  spatiotemporal = "iid",
  mesh = mesh,
  offset = "log_area_km2",
  silent = SILENT
)
sanity(fit12)
fit12

fits <- list(fit1, fit2, fit4, fit5, fit9, fit10, fit11, fit12)

fitted_yrs <- sort(unique(sp_dat$year))
nd <- make_grid(syn_grid, years = fitted_yrs) |>
  mutate(
    log_depth = log(depth),
    fyear = as.factor(year)
  ) |>
  rename(survey_abbrev = survey)
nd <- left_join(nd, spline_lu, by = join_by(year))

future::plan(future::multisession, workers = length(fits))
preds <-
  # purrr::map(fits, function(.x) {
  furrr::future_map(fits, function(.x) {
    if (inherits(.x, "sdmTMB")) {
      out <- predict(.x, newdata = nd, return_tmb_object = TRUE)
    } else {
      out <- NA
    }
    out
  })

indices <-
  # purrr::map(preds, function(.x) {
  furrr::future_map(preds, function(.x) {
    if (length(.x) > 1) {
      get_index(.x, bias_correct = F) # FIXME
    }
  })
future::plan(future::sequential)

model_lookup <-
  tibble(
    id = 1:7,
    desc = c(
      "st RW", # 6
      "st IID, depth covariate", # 3
      "st IID, s(year)", # 5
      "st IID, no covariate, as.factor(year)", # 2????
      "st IID, time-varying RW", # 4
      "st IID (1|year)", # 1
      "spatial only"
    ), # 7
    order = c(6, 3, 5, 2, 4, 1, 7)
  ) |>  # try matching order of the simulated plots
  filter(id %in% c(1, 2, 4, 5)) |>
  mutate(id = as.numeric(as.factor(id)))

model_lookup <- bind_rows(
  model_lookup,
  tibble(id = 5, desc = "st RW, as.factor(year)", order = 99),
  tibble(id = 6, desc = "st off, spline SVC", order = 99),
  tibble(id = 7, desc = "st IID spline SVC", order = 99),
  tibble(id = 8, desc = "st IID, as.factor(year) + as.factor(survey)", order = 99)
) |>
  mutate(order = seq_len(n()))


reg_lu <-
  sp_dat |>
  select(year, survey_abbrev) |>
  group_by(year) |>
  summarize(concat_regions = paste(sort(unique(survey_abbrev)), collapse = ", ")) |>
  mutate(concat_regions = gsub("SYN ", "", concat_regions))

index_df <-
  bind_rows(indices, .id = "id") |>
  mutate(id = as.numeric(id)) |>
  as_tibble() |>
  left_join(model_lookup) |>
  left_join(reg_lu) |>
  # Add this for quick and dirty north/south, but note that WCVI data from 2021
  # was in this dataset used.
  mutate(sampled_region = if_else(year %% 2 == 1, "north", "south"))

# Plot index over time ---------------------------------------------------------
ggplot(data = index_df, aes(x = year, y = est, ymin = lwr, ymax = upr)) +
  # geom_pointrange(aes(colour = sampled_region)) +
  geom_pointrange(aes(colour = concat_regions)) +
  geom_ribbon(alpha = 0.15, colour = NA) +
  # scale_colour_manual(values = c("#66C2A5", "#FC8D62")) +
  scale_colour_brewer(palette = "Dark2") +
  # ggthemes::scale_colour_colorblind() +
  labs(colour = "Sampled region", x = "Year", y = "Index (log distributed)") +
  facet_wrap(~ forcats::fct_reorder(desc, order)) +
  scale_y_log10() +
  ggtitle(stringr::str_to_title(unique(sp_dat$species_common_name)))

# does using year factor coefs fix it!?
b <- tidy(fit12, conf.int = TRUE) |> filter(grepl("year", term)) |>
  mutate(year = as.numeric(gsub("as\\.factor\\(year\\)", "", term))) |>
  mutate(sampled_region = if_else(year %% 2 == 1, "north", "south"))

ggplot(b, aes(x = year, y = exp(estimate), ymin = exp(conf.low), ymax = exp(conf.high))) +
  geom_pointrange(aes(colour = sampled_region)) +
  geom_ribbon(alpha = 0.20, colour = NA) +
  scale_colour_manual(values = c("#66C2A5", "#FC8D62")) +
  labs(colour = "Sampled region") +
  scale_y_log10()
# no


