# just work with QCS and HS and pretend those were biennially sampled

library(ggplot2)
library(dplyr)
library(sdmTMB)
library(furrr)
theme_set(ggsidekick::theme_sleek())
Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1"
)

d <- readRDS("~/src/gfsynopsis-2024/report/data-cache-2025-03/pacific-cod.rds")$survey_sets
d <- filter(d, grepl("^SYN", survey_abbrev))
d <- mutate(d,
  density_kgkm2 = density_kgpm2 * 1e6,
  log_depth = log(depth_m),
  area_swept1 = doorspread_m * (speed_mpm * duration_min),
  area_swept2 = tow_length_m * doorspread_m,
  area_swept = ifelse(!is.na(area_swept2), area_swept2, area_swept1)
) |>
  filter(!is.na(area_swept))
# make it cleanly biennial for visualization:
d <- filter(d, !year %in% c(2003, 2004, 2020))
# d <- filter(d, !(year %in% c(2021) & survey_abbrev == "SYN WCVI"))
d <- filter(d, !(year %in% c(2007) & survey_abbrev == "SYN WCHG"))
d <- filter(d, survey_abbrev != "SYN WCHG")
d <- filter(d, survey_abbrev != "SYN WCVI")

d_raw <- d

target_surveys <- c("SYN QCS", "SYN HS")
d_raw <- d_raw |>
  filter(survey_abbrev %in% target_surveys)

analysis_years <- sort(unique(d_raw$year))
truth_label <- "Truth: both regions sampled all years"
dropped_survey <- "SYN HS"
scenario_drop_years <- list(
  "Scenario 1: both 2005-2015, then HS dropped (2016+)" = analysis_years[analysis_years >= 2016],
  "Scenario 2: HS dropped 2005-2015, then both (2016+)" = analysis_years[analysis_years >= 2005 & analysis_years <= 2015],
  "Scenario 3: both all years except HS dropped in 2015" = 2015,
  "Truth: both regions sampled all years" = integer(0)
)

build_scenario_data <- function(dat, drop_years, drop_survey) {
  dat |>
    filter(!(year %in% drop_years & survey_abbrev == drop_survey)) |>
    mutate(
      observed = density_kgkm2,
      log_effort = log(area_swept)
    )
}

family <- tweedie()

grid <- gfplot::synoptic_grid |> select(-survey_domain_year, -utm_zone) |>
  filter(survey %in% target_surveys)

run_scenario <- function(label, drop_years) {
  cli::cli_alert_info("Running scenario: {label}")
  d <- build_scenario_data(d_raw, drop_years = drop_years, drop_survey = dropped_survey)
  print(table(d$year, d$survey_abbrev))

  d <- add_utm_columns(d)
  mesh <- make_mesh(d, c("X", "Y"), cutoff = 8)

  fit <- sdmTMB(
    observed ~ factor(year),
    offset = "log_effort",
    data = d,
    time = "year",
    spatial = "on",
    spatiotemporal = "iid",
    mesh = mesh,
    silent = FALSE,
    family = family
  )

  nd <- replicate_df(grid, "year", sort(unique(d$year)))
  p <- predict(fit, newdata = nd, return_tmb_object = TRUE)
  ind <- get_index(p, bias_correct = TRUE)

  coverage_by_year <- d |>
    distinct(year, survey_abbrev) |>
    summarise(
      n_regions = n(),
      coverage = ifelse(n_regions > 1L, "Both regions", paste0(sub("^SYN ", "", survey_abbrev[1]), " only")),
      .by = year
    ) |>
    select(year, coverage)

  ind |>
    left_join(coverage_by_year, by = "year") |>
    mutate(
      scenario = factor(label, levels = names(scenario_drop_years))
    )
}

scenario_labels <- names(scenario_drop_years)
workers <- max(1L, future::availableCores() - 1L)
future::plan(future::multisession, workers = workers)
index_results <- furrr::future_map_dfr(
  scenario_labels,
  function(lbl) run_scenario(lbl, scenario_drop_years[[lbl]]),
  .options = furrr::furrr_options(seed = TRUE)
)
future::plan(future::sequential)

index_truth <- index_results |>
  filter(scenario == truth_label)
index_panel <- index_results |>
  filter(scenario != truth_label) |>
  mutate(coverage = factor(coverage, levels = c("Both regions", "QCS only", "HS only")))
index_truth_panel <- tidyr::crossing(
  scenario = unique(index_panel$scenario),
  index_truth |>
    select(year, est, lwr, upr)
)

ggplot() +
  geom_ribbon(
    data = index_truth_panel,
    aes(year, ymin = lwr, ymax = upr),
    fill = "grey70",
    alpha = 0.30
  ) +
  geom_line(
    data = index_truth_panel,
    aes(year, est),
    colour = "black"
  ) +
  geom_ribbon(
    data = index_panel,
    aes(year, ymin = lwr, ymax = upr),
    fill = "grey50",
    alpha = 0.12
  ) +
  geom_line(
    data = index_panel,
    aes(year, est),
    colour = "grey50"
  ) +
  geom_point(
    data = index_panel,
    aes(year, est, colour = coverage)
  ) +
  scale_colour_manual(
    values = c(
      "Both regions" = "black",
      "QCS only" = "#D95F02",
      "HS only" = "#1B9E77"
    ),
    name = "Observed coverage"
  ) +
  facet_wrap(~scenario) +
  labs(
    x = "Year",
    y = "Estimated index",
    title = "Pacific cod SYN QCS + HS sampling scenarios",
    subtitle = "Black line and grey band show truth (both regions sampled all years)"
  )
ggsave("figs/pcod-qcs-hs-scenario-comparison.pdf", width = 10, height = 7)
