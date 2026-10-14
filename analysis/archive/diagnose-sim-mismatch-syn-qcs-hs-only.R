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

base_qcs_drop_yrs <- seq(2001, 2025, 4)
base_hs_drop_yrs <- seq(2003, 2025, 4)
all_drop_years <- sort(unique(c(base_qcs_drop_yrs, base_hs_drop_yrs)))

progressive_overlap_years <- sort(intersect(seq(2005, 2023, 2), unique(d_raw$year)), decreasing = TRUE)
progressive_overlap_scenarios <- lapply(seq_along(progressive_overlap_years), function(i) {
  progressive_overlap_years[seq_len(i)]
})
progressive_overlap_labels <- vapply(progressive_overlap_scenarios, function(x) {
  x_sorted <- sort(x)
  if (length(x_sorted) == 1L) {
    paste0("Full coverage: ", x_sorted[1])
  } else {
    paste0("Full coverage: ", min(x_sorted), "-", max(x_sorted))
  }
}, character(1))
names(progressive_overlap_scenarios) <- progressive_overlap_labels

overlap_scenarios <- c(
  list("No full-coverage years" = integer(0)),
  progressive_overlap_scenarios,
  list("Full coverage: all years" = all_drop_years)
)

build_scenario_data <- function(dat, overlap_years) {
  qcs_drop_yrs <- setdiff(base_qcs_drop_yrs, overlap_years)
  hs_drop_yrs <- setdiff(base_hs_drop_yrs, overlap_years)

  d_scen <- dat |>
    filter(!(year %in% qcs_drop_yrs & survey_abbrev == "SYN QCS")) |>
    filter(!(year %in% hs_drop_yrs & survey_abbrev == "SYN HS"))

  years_with_both <- d_scen |>
    distinct(year, survey_abbrev) |>
    count(year, name = "n_regions") |>
    filter(n_regions > 1L) |>
    pull(year)

  d_scen |>
    mutate(
      biennial_region = ifelse(survey_abbrev == "SYN QCS", 1L, 2L),
      biennial_region = ifelse(year %in% years_with_both, 3L, biennial_region),
      observed = density_kgkm2,
      log_effort = log(area_swept)
    )
}

family <- tweedie()

grid <- gfplot::synoptic_grid |> select(-survey_domain_year, -utm_zone) |>
  filter(survey != "SYN WCHG") |> filter(survey != "SYN WCVI")
# max(grid$Y[grid$survey == "SYN WCVI"]) - min(grid$Y[grid$survey == "SYN QCS"])
ggplot(grid, aes(X, Y, colour = survey)) + geom_point()
table(grid$survey)

run_scenario <- function(label, overlap_years) {
  cli::cli_alert_info("Running scenario: {label}")
  d <- build_scenario_data(d_raw, overlap_years)
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

  ind |>
    mutate(
      scenario = factor(label, levels = names(overlap_scenarios)),
      overlap_years = if (length(overlap_years) == 0L) {
        "none"
      } else {
        paste(sort(overlap_years), collapse = ", ")
      }
    )
}

scenario_labels <- names(overlap_scenarios)
workers <- max(1L, future::availableCores() - 1L)
future::plan(future::multisession, workers = workers)
index_results <- furrr::future_map_dfr(
  scenario_labels,
  function(lbl) run_scenario(lbl, overlap_scenarios[[lbl]]),
  .options = furrr::furrr_options(seed = TRUE)
)
future::plan(future::sequential)

ggplot(index_results, aes(year, est, ymin = lwr, ymax = upr, colour = scenario, fill = scenario)) +
  geom_ribbon(alpha = 0.08, colour = NA) +
  geom_line() +
  geom_point() +
  labs(
    x = "Year",
    y = "Estimated index",
    colour = "Scenario",
    fill = "Scenario",
    title = "Effect of adding overlap years between QCS and HS"
  )

all_data_label <- "Full coverage: all years"

coverage_lu <- purrr::imap_dfr(overlap_scenarios, \(overlap_yrs, scen_lbl) {
  tibble::tibble(
    scenario = scen_lbl,
    year = sort(unique(index_results$year))
  ) |>
    mutate(
      coverage_type = case_when(
        scen_lbl == all_data_label ~ "Full coverage",
        year %in% overlap_yrs ~ "Full coverage",
        year %in% base_qcs_drop_yrs ~ "HS only",
        year %in% base_hs_drop_yrs ~ "QCS only",
        TRUE ~ "Full coverage"
      )
    )
})

plot_dat <- index_results |>
  mutate(scenario = as.character(scenario)) |>
  left_join(coverage_lu, by = c("scenario", "year")) |>
  mutate(coverage_type = factor(coverage_type, levels = c("Full coverage", "QCS only", "HS only")))

index_all_data <- plot_dat |>
  filter(scenario == all_data_label)
index_panel <- plot_dat |>
  filter(scenario != all_data_label) |>
  mutate(panel_scenario = factor(scenario, levels = names(overlap_scenarios)[names(overlap_scenarios) != all_data_label]))
panel_scenarios <- unique(index_panel$panel_scenario)
index_all_data_panel <- tidyr::crossing(
  panel_scenario = panel_scenarios,
  index_all_data |>
    select(year, est, lwr, upr)
)

ggplot() +
  geom_ribbon(
    data = index_all_data_panel,
    aes(year, ymin = lwr, ymax = upr),
    fill = "grey70",
    alpha = 0.30
  ) +
  geom_line(
    data = index_all_data_panel,
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
    aes(year, est, colour = coverage_type)
  ) +
  scale_colour_manual(
    values = c(
      "Full coverage" = "black",
      "QCS only" = "#1B9E77",
      "HS only" = "#D95F02"
    ),
    name = "Observed coverage"
  ) +
  facet_wrap(~panel_scenario) +
  labs(
    x = "Year",
    y = "Estimated index",
    title = "Pacific Cod QCS + HS",
    subtitle = "Grey band/black line: full coverage in all years"
  )
ggsave("figs/pcod-qcs-hs-retro-add.pdf", width = 9, height = 6)
