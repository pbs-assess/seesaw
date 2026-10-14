# Lognormal priors on sigma_V for the IID RF, RW year model, Pacific cod under
# pretend biennial QCS / HS sampling: the priors, the estimated sigma_V, the
# indexes, and their performance against the full-data reference.
#
# Plots only: the fits (sigma_V_prior() in analysis/priors.R: lognormal with
# its mode at the target SD, as throughout the project) and the full-data
# reference come from analysis/qcs-hs-biennial-prior-sensitivity.R (run that
# first). Time is survey occasion, so sigma_V is per occasion (2 calendar
# years).

library(dplyr)
library(ggplot2)

theme_set(ggsidekick::theme_sleek())

source(here::here("analysis/metric-functions.R"))

fits_file <- here::here("data-generated/qcs-hs-biennial-prior-sensitivity-lognormal-target.rds")
stopifnot(file.exists(fits_file))

surveyjoin::cache_data()
surveyjoin::load_sql_data()

species <- "pacific cod"
survey_names <- c("SYN QCS", "SYN HS")
ci_level <- 0.5
seesaw_window <- 8L

# As in analysis/qcs-hs-biennial-prior-sensitivity.R
prior_grid <- tidyr::crossing(target_sd = c(0.1, 0.3, 0.5), sdlog = c(0.1, 0.3, 0.5, 0.8)) |>
  mutate(prior = paste0("target = ", target_sd, ", sdlog = ", sdlog))

# Survey occasions (as in analysis/qcs-hs-biennial-prior-sensitivity.R) ------

occasions <- surveyjoin::get_data(species, regions = "pbs") |>
  mutate(year = lubridate::year(lubridate::ymd(date))) |>
  filter(survey_name %in% survey_names, effort > 0) |>
  distinct(year, survey_name) |>
  # QCS was sampled alone in 2003 and 2004; keep years with both surveys
  count(year) |>
  filter(n == length(survey_names)) |>
  transmute(
    cal_year = year,
    year = match(year, sort(unique(year))),
    primary_survey = if_else(year %% 2L == 1L, "SYN HS", "SYN QCS")
  )

# Fits -----------------------------------------------------------------------

res <- readRDS(fits_file)
stopifnot(res$truth$sanity_ok, all(prior_grid$prior %in% names(res)))
truth <- res$truth$index |>
  select(year, true_est = est, true_lwr = lwr, true_upr = upr) |>
  left_join(occasions, by = "year")
res <- res[prior_grid$prior]

fit_info <- mutate(prior_grid, sanity_ok = purrr::map_lgl(res[prior], "sanity_ok"))
if (!all(fit_info$sanity_ok)) {
  cli::cli_warn("Fits failing sanity checks (kept, flagged): {fit_info$prior[!fit_info$sanity_ok]}")
}

indexes <- purrr::imap_dfr(res, \(r, name) {
  if (!is.null(r$index)) mutate(r$index, prior = name)
}) |>
  left_join(truth, by = "year") |>
  left_join(fit_info, by = "prior") |>
  mutate(log_error = log_est - log(true_est))

sigma_V <- purrr::imap_dfr(res, \(r, name) {
  if (!is.null(r$sigma_V)) mutate(r$sigma_V, prior = name)
}) |>
  left_join(fit_info, by = "prior")

# Metrics ----------------------------------------------------------------------

phase_of <- \(p) if_else(p == "SYN HS", 0.5, -0.5)
seesaw <- function(est, primary_survey, year) {
  mw <- moving_window(est, year = year, window = seesaw_window, phase = phase_of(primary_survey))
  tibble::tibble(
    A = period2_metric(est, year = year, phase = phase_of(primary_survey))[["A"]],
    A_max = max(mw$A)
  )
}

q <- qnorm(1 - (1 - ci_level) / 2)
metrics <- indexes |>
  arrange(prior, year) |>
  group_by(target_sd, sdlog, sanity_ok) |>
  summarise(
    rmse = sqrt(mean(log_error^2)),
    bias = mean(log_error),
    coverage = mean(abs(log_error) <= q * se),
    seesaw(est, primary_survey, year),
    .groups = "drop"
  ) |>
  left_join(
    sigma_V |>
      select(target_sd, sdlog, component, estimate) |>
      tidyr::pivot_wider(names_from = component, values_from = estimate, names_prefix = "sigma_V_"),
    by = c("target_sd", "sdlog")
  ) |>
  arrange(target_sd, sdlog)

truth_seesaw <- seesaw(truth$true_est, truth$primary_survey, truth$year)
cat("\nTruth A:", round(truth_seesaw$A, 1), " A_max:", round(truth_seesaw$A_max, 1), "\n\n")
metrics |>
  mutate(across(where(is.numeric), \(x) round(x, 3))) |>
  print(n = Inf, width = Inf)

# Plots ----------------------------------------------------------------------

metric_levels <- c("rmse", "bias", "coverage", "A", "A_max", "sigma_V_encounter", "sigma_V_positive")
metrics |>
  tidyr::pivot_longer(any_of(metric_levels), names_to = "metric") |>
  mutate(metric = factor(metric, levels = metric_levels)) |>
  ggplot(aes(factor(sdlog), value, colour = factor(target_sd), group = target_sd)) +
  geom_line() +
  geom_point() +
  facet_wrap(~metric, scales = "free_y", nrow = 1) +
  scale_colour_viridis_d(end = 0.85) +
  labs(x = "Prior sdlog", y = NULL, colour = "Prior target SD")
ggsave(here::here("figs/qcs-hs-biennial-sigma-v-prior-family-metrics.pdf"), width = 14, height = 3)

# Prior densities with the estimated sigma_V for each component
densities <- prior_grid |>
  tidyr::crossing(x = seq(0.001, 1, length.out = 400)) |>
  mutate(density = dlnorm(x, meanlog = log(target_sd) + sdlog^2, sdlog = sdlog)) |>
  # Scale so each panel peaks at 1
  group_by(target_sd, sdlog) |>
  mutate(density = density / max(density)) |>
  ungroup()

panel_labs <- \(d) mutate(d,
  target_lab = factor(paste("Target =", target_sd)),
  sdlog_lab = factor(paste("sdlog =", sdlog), levels = paste("sdlog =", sort(unique(prior_grid$sdlog))))
)

ggplot(panel_labs(densities), aes(x, density)) +
  geom_line(colour = "grey40") +
  geom_vline(aes(xintercept = estimate, linetype = component), data = panel_labs(sigma_V), colour = "#1B9E77") +
  facet_grid(target_lab ~ sdlog_lab) +
  scale_x_continuous(breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1"), limits = c(0, 1), expand = expansion(mult = 0.02)) +
  scale_y_continuous(breaks = NULL, expand = expansion(mult = c(0, 0.05))) +
  scale_linetype_manual(values = c(encounter = "solid", positive = "22")) +
  labs(
    x = expression(sigma[V] ~ "(year random walk SD)"), y = "Prior density (scaled)",
    linetype = expression("Estimated" ~ sigma[V])
  ) +
  theme(legend.position = "bottom")
ggsave(here::here("figs/qcs-hs-biennial-sigma-v-prior-family-densities.pdf"), width = 7, height = 5)

# Indexes for each prior, same layout as the densities; each series scaled to
# a geometric mean of 1
geomean <- \(x) exp(mean(log(x)))
truth_panels <- truth |>
  mutate(across(c(true_est, true_lwr, true_upr), \(x) x / geomean(true_est))) |>
  tidyr::crossing(distinct(prior_grid, target_sd, sdlog)) |>
  panel_labs()
indexes_scaled <- indexes |>
  group_by(prior) |>
  mutate(est = est / geomean(est)) |>
  ungroup()

ggplot(panel_labs(indexes_scaled), aes(cal_year, est)) +
  geom_ribbon(aes(cal_year, ymin = true_lwr, ymax = true_upr), data = truth_panels, fill = "grey85", inherit.aes = FALSE) +
  geom_line(aes(y = true_est), data = truth_panels, colour = "black") +
  geom_line(colour = "#1B9E77") +
  geom_point(aes(shape = primary_survey), colour = "#1B9E77") +
  facet_grid(target_lab ~ sdlog_lab) +
  scale_x_continuous(breaks = c(2010, 2020)) +
  scale_shape_manual(values = c(`SYN HS` = 19, `SYN QCS` = 21)) +
  labs(
    x = "Year", y = "Relative biomass index", shape = "Region sampled",
    caption = "Black and grey: full-data reference index and its 95% CI."
  ) +
  theme(legend.position = "bottom")
ggsave(here::here("figs/qcs-hs-biennial-sigma-v-prior-family-indexes.pdf"), width = 7, height = 5)
