# Correlates of the spurious seesaw across species, from the fitted operating
# models cached by analysis/13-qcs-hs-omega-swap.R (full run, not the test run).
#
# Per species: random-effect parameters (sigma_O, sigma_E, range) of each delta
# component, the fitted HS - QCS difference in mean density from the
# operating model, prevalence, and mean depth, against the median seesaw error
# of the biennial fits (conditional arm: the real geography).

library(dplyr)
library(ggplot2)
library(sdmTMB)
theme_set(ggsidekick::theme_sleek())

dir <- here::here("data-generated/qcs-hs-omega-swap")
m <- readRDS(here::here("data-generated/qcs-hs-omega-swap-metrics.rds"))
species <- names(m$metrics_rep |> distinct(species) |> tidyr::pivot_wider(names_from = species, values_from = species))
slug <- \(s) gsub("[^a-z0-9]+", "-", s)

# Covariates from the operating model fit ----------------------------------------

covariates_one <- function(sp) {
  s <- readRDS(file.path(dir, "setups", paste0(slug(sp), ".rds")))
  rp <- s$ran_pars |>
    filter(term %in% c("sigma_O", "sigma_E", "range_min")) |>
    mutate(name = paste0(term, "_", if_else(model == 1, "enc", "pos"))) |>
    select(name, estimate) |>
    tidyr::pivot_wider(names_from = name, values_from = estimate)

  # Expected density on the prediction grid, by region (HS vs. QCS), averaged
  # over years
  p <- predict(s$fit_om, newdata = s$grid, return_tmb_object = FALSE)
  p$dens <- exp(p$est1 + p$est2)
  reg <- p |>
    summarise(
      log_dens = log(weighted.mean(dens, area)),
      mean_log_dens = weighted.mean(est1 + est2, area),
      .by = c(year, hs)
    ) |>
    summarise(across(c(log_dens, mean_log_dens), mean), .by = hs) |>
    arrange(hs)
  tibble::tibble(
    species = sp,
    density_diff = diff(reg$log_dens), # HS - QCS, log of mean density
    mean_log_density_diff = diff(reg$mean_log_dens),
    prevalence = mean(s$dat$catch_weight > 0),
    prevalence_diff = mean(s$dat$catch_weight[s$dat$hs == 1] > 0) -
      mean(s$dat$catch_weight[s$dat$hs == 0] > 0),
    mean_depth = mean(s$dat$depth_m)
  ) |>
    bind_cols(rp)
}
covs <- purrr::map_dfr(species, covariates_one)

# Seesaw outcomes -------------------------------------------------------------------

outcome <- m$metrics_rep |>
  filter(omega == "conditional") |>
  group_by(species) |>
  summarise(
    beta_error = median(beta_error, na.rm = TRUE),
    A_max_excess = median(A_max_excess, na.rm = TRUE),
    contrast_omega = mean(contrast_omega),
    .groups = "drop"
  )
d <- left_join(outcome, covs, by = "species") |>
  mutate(
    sigma_O_total = sqrt(sigma_O_enc^2 + sigma_O_pos^2),
    sigma_E_total = sqrt(sigma_E_enc^2 + sigma_E_pos^2),
    abs_beta_error = abs(beta_error),
    abs_contrast_omega = abs(contrast_omega),
    abs_density_diff = abs(density_diff),
    sigma_ratio = sigma_E_total / sigma_O_total
  )
print(d |> mutate(across(where(is.numeric), \(x) round(x, 2))), n = Inf, width = Inf)

# Plots ------------------------------------------------------------------------------

predictors <- c(
  abs_density_diff = "|HS - QCS| log mean density (fitted)",
  abs_contrast_omega = "|HS - QCS| omega contrast (log)",
  sigma_O_total = "sigma_O (both components combined)",
  sigma_E_total = "sigma_E (both components combined)",
  sigma_ratio = "sigma_E / sigma_O",
  range_min_enc = "Range, encounter (km)",
  range_min_pos = "Range, positive (km)",
  prevalence = "Prevalence",
  mean_depth = "Mean depth (m)"
)
long <- d |>
  tidyr::pivot_longer(all_of(names(predictors)), names_to = "predictor", values_to = "x") |>
  mutate(predictor = factor(predictor, levels = names(predictors), labels = predictors)) |>
  group_by(predictor) |>
  mutate(rho = cor(x, abs_beta_error, method = "spearman", use = "complete.obs")) |>
  ungroup() |>
  mutate(panel = factor(
    paste0(predictor, "\nSpearman ", sprintf("%.2f", rho)),
    levels = unique(paste0(predictor, "\nSpearman ", sprintf("%.2f", rho)))
  ))

ggplot(long, aes(x, abs_beta_error)) +
  geom_point() +
  ggrepel::geom_text_repel(aes(label = tools::toTitleCase(species)), size = 2.5, max.overlaps = Inf) +
  facet_wrap(~panel, scales = "free_x") +
  labs(
    x = NULL,
    y = "|Seesaw error| (log ratio, median over replicates)",
    caption = "One point per species; QCS/HS biennial self-test, conditional omega (analysis/13)."
  )
ggsave(here::here("figs/seesaw-correlates.pdf"), width = 11, height = 8, device = cairo_pdf)

saveRDS(d, here::here("data-generated/seesaw-correlates.rds"))

# Per replicate ------------------------------------------------------------------------
# Each replicate's realized contrast (omega, density) against its own seesaw
# error; sigma_O is a species-level value, so those points form vertical strips

rep_d <- m$metrics_rep |>
  filter(omega == "conditional", !is.na(beta_error)) |>
  left_join(select(d, species, sigma_O_total), by = "species") |>
  mutate(
    abs_beta_error = abs(beta_error),
    abs_contrast_omega = abs(contrast_omega),
    abs_contrast_density = abs(contrast_density)
  )

rep_predictors <- c(
  abs_contrast_omega = "|HS - QCS| omega contrast (log), replicate",
  abs_contrast_density = "|HS - QCS| log mean density contrast, replicate",
  sigma_O_total = "sigma_O (both components combined), species"
)
rep_long <- rep_d |>
  tidyr::pivot_longer(all_of(names(rep_predictors)), names_to = "predictor", values_to = "x") |>
  mutate(predictor = factor(predictor, levels = names(rep_predictors), labels = rep_predictors)) |>
  group_by(predictor) |>
  mutate(rho = cor(x, abs_beta_error, method = "spearman")) |>
  ungroup()
rep_long <- mutate(rep_long, panel = factor(
  paste0(predictor, "\nSpearman ", sprintf("%.2f", rho)),
  levels = unique(paste0(predictor, "\nSpearman ", sprintf("%.2f", rho)))
))
print(distinct(rep_long, predictor, rho))
# Within species (removes species differences)
rep_d |>
  group_by(species) |>
  summarise(
    rho_omega = cor(abs_contrast_omega, abs_beta_error, method = "spearman"),
    rho_density = cor(abs_contrast_density, abs_beta_error, method = "spearman")
  ) |>
  print(n = Inf)

ggplot(rep_long, aes(x, abs_beta_error, colour = tools::toTitleCase(species))) +
  geom_point(alpha = 0.5) +
  facet_wrap(~panel, scales = "free_x") +
  labs(
    x = NULL, y = "|Seesaw error| (log ratio), replicate", colour = NULL,
    caption = "One point per species and replicate; conditional omega (analysis/13)."
  ) +
  theme(legend.position = "bottom")
ggsave(here::here("figs/seesaw-correlates-replicate.pdf"), width = 11, height = 4.5, device = cairo_pdf)
