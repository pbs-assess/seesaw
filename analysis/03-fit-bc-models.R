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

source(here::here("analysis/fit-index-models.R"))

do_fit_syn <- function(.sp) {
  RhpcBLASctl::blas_set_num_threads(1L)
  RhpcBLASctl::omp_set_num_threads(1L)

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
  grid$survey_name <- "SYN WCVI"
  grid <- clamp_depth(grid, dat)

  mesh <- make_mesh(dat, c("X", "Y"), cutoff = 10)
  mesh_all <- make_mesh(dat_all, c("X", "Y"), mesh = mesh$mesh)

  fit_index_models(
    dat = dat,
    grid = grid,
    mesh = mesh,
    response = "catch_weight",
    family = delta_gamma(type = "poisson-link"),
    offset = log(dat$effort),
    all_data = list(data = dat_all, mesh = mesh_all, offset = log(dat_all$effort))
  ) |>
    mutate(species = .sp)
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


RhpcBLASctl::blas_set_num_threads(1L)
RhpcBLASctl::omp_set_num_threads(1L)

future::plan(future::multisession, workers = future::availableCores()/2)
# out <- purrr::map_dfr(spp_to_fit_syn, do_fit_syn)
out <- furrr::future_map_dfr(spp_to_fit_syn, do_fit_syn)
dir.create("data-generated", showWarnings = FALSE)
saveRDS(out, file = "data-generated/bc-indexes4.rds")

out <- readRDS("data-generated/bc-indexes4.rds") |>
  filter(!species %in% c("redbanded rockfish", "shortspine thornyhead"))

dat <- surveyjoin::get_data("pacific cod", regions = "pbs") |>
  mutate(year = lubridate::year(lubridate::ymd(date))) |>
  filter(!(year %in% c(2003, 2004, 2020))) |>
  tidyr::drop_na(effort, catch_weight, depth_m) |>
  filter(!(year == 2007 & survey_name == "SYN WCHG")) |>
  filter(!(year == 2021 & survey_name == "SYN WCVI"))
lu <- data.frame(year = sort(unique(dat$year)))
lu$even <- lu$year %% 2 == 0
lu$survey_group <- ifelse(lu$even, "WCHG + WCVI", "QCS + HS")

source(here::here("analysis/metric-functions.R"))

# Phase is coded by calendar year: even = WCHG + WCVI, odd = QCS + HS
seesaw_window <- 10L
seesaw_mw <- out |>
  arrange(species, model, year) |>
  group_by(species, model) |>
  group_modify(\(.x, .y) moving_window(.x$est, year = .x$year, window = seesaw_window)) |>
  ungroup()

seesaw_summary <- seesaw_mw |>
  summarise(mean_A = mean(A), max_A = max(A), .by = c(species, model))

out |>
  left_join(lu) |>
  left_join(seesaw_summary) |>
  group_by(species, model) |>
  mutate(geomean = exp(mean(log(est))), est = est / geomean, lwr = lwr / geomean, upr = upr / geomean) |>
  ggplot(aes(year, log(est), ymin = log(lwr), ymax = log(upr))) +
  geom_ribbon(fill = "grey90") +
  geom_linerange(aes(colour = survey_group)) +
  geom_point(aes(colour = survey_group), size = 2) +
  scale_colour_brewer(palette = "Dark2") +
  facet_grid(forcats::fct_reorder(model, mean_A) ~ species) +
  ylab("Biomass index") +
  xlab("Year") +
  labs(colour = "Survey\ngrouping") +
  ggsidekick::theme_sleek()
ggsave("figs/bc-testing2.pdf", width = 30, height = 15)

out |>
  left_join(lu) |>
  group_by(model) |>
  summarise(n = n())

make_fig <- function(.ylab = "", include_all_data = FALSE) {
  if (!include_all_data) {
    dat_mw <- seesaw_mw |>
      filter(!grepl("all data", model)) |>
      mutate(all_data = FALSE)
  } else {
    dat_mw <- seesaw_mw |> mutate(all_data = grepl("all data", model))
  }

  dat_mw <- dat_mw |>
    mutate(model = reorder(model, A, FUN = mean))

  g <- dat_mw |>
    ggplot(aes(model, A)) +
    coord_flip(ylim = c(0, 200)) +
    ylab(.ylab) +
    ggsidekick::theme_sleek() +
    theme(axis.title.y = element_blank(), panel.grid.major = element_line(colour = "grey90", linewidth = 0.3), panel.grid.minor = element_line(colour = "grey90", linewidth = 0.3))

  if (include_all_data) {
    blue <- RColorBrewer::brewer.pal(8, "Blues")[3]
    orange <- RColorBrewer::brewer.pal(8, "Oranges")[3]
    g <- g + geom_violin(scale = "width", mapping = aes(colour = all_data, fill = all_data)) +
      scale_colour_manual(values = c(blue, orange)) +
      scale_fill_manual(values = c(blue, orange)) +
      guides(colour = "none", fill = "none")
  } else {
    blue <- RColorBrewer::brewer.pal(8, "Blues")[3]
    g <- g + geom_violin(scale = "width", colour = blue, fill = blue)
  }

  g +
    geom_point(position = position_jitter(width = 0.1), colour = "grey25", alpha = 0.3) +
    geom_point(stat = "summary", fun = mean, colour = "black") +
    scale_y_sqrt(limits = c(0, NA), expand = expansion(mult = c(0, 0.05)))
}

a_lab <- paste0("Estimated biennial amplitude (%)\nacross ", seesaw_window, "-survey windows")
make_fig(a_lab)
ggsave("figs/bc-trawl-A-moving-window.pdf", width = 5, height = 3.5)

make_fig(a_lab, include_all_data = TRUE)
ggsave("figs/bc-trawl-A-moving-window2.pdf", width = 5.5, height = 3.5)
