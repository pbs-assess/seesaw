library(dplyr)
source(here::here("analysis/R/models.R")) # sigma_V_meanlog()
source(here::here("analysis/R/plotting.R"))

ps <- readRDS(here::here("data-generated/prior-sensitivity.rds"))
prior_grid <- ps$prior_grid
# The no-prior fit is not shown
fit_info <- filter(ps$fit_info, prior %in% prior_grid$prior)
if (!all(fit_info$sanity_ok)) {
  cli::cli_warn("Fits failing sanity checks (kept, flagged): {fit_info$prior[!fit_info$sanity_ok]}")
}
reference <- ps$reference |>
  select(year, ref_est = est, ref_lwr = lwr, ref_upr = upr) |>
  left_join(ps$occasions, by = "year")
indexes <- ps$indexes |>
  inner_join(fit_info, by = "prior") |>
  left_join(ps$occasions, by = "year")
sigma_V <- inner_join(ps$sigma_V, fit_info, by = "prior")

panel_labs <- \(d) mutate(d,
  target_lab = factor(paste("Mode =", target_sd)),
  sdlog_lab = factor(paste("sdlog =", sdlog), levels = paste("sdlog =", sort(unique(prior_grid$sdlog))))
)

# each panel scaled to peak at 1
densities <- prior_grid |>
  crossing(x = seq(0.001, 1, length.out = 400)) |>
  mutate(density = dlnorm(x, meanlog = sigma_V_meanlog(target_sd, sdlog), sdlog = sdlog)) |>
  mutate(density = density / max(density), .by = c(target_sd, sdlog))

g_densities <- ggplot(panel_labs(densities), aes(x, density)) +
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
ggsave(here::here("figs/prior-sensitivity-densities.pdf"), g_densities, width = 7, height = 5)

geomean <- \(x) exp(mean(log(x)))
reference_panels <- reference |>
  mutate(across(c(ref_est, ref_lwr, ref_upr), \(x) x / geomean(ref_est))) |>
  crossing(distinct(prior_grid, target_sd, sdlog)) |>
  panel_labs()
indexes_scaled <- mutate(indexes, est = est / geomean(est), .by = prior)

g_indexes <- ggplot(panel_labs(indexes_scaled), aes(cal_year, est)) +
  geom_ribbon(aes(cal_year, ymin = ref_lwr, ymax = ref_upr), data = reference_panels, fill = "grey85", inherit.aes = FALSE) +
  geom_line(aes(y = ref_est), data = reference_panels, colour = "black") +
  geom_line(colour = "#1B9E77") +
  geom_point(aes(shape = primary_survey), colour = "#1B9E77") +
  facet_grid(target_lab ~ sdlog_lab) +
  scale_x_continuous(breaks = c(2010, 2020)) +
  scale_shape_manual(values = c(`SYN HS` = 19, `SYN QCS` = 21)) +
  labs(x = "Year", y = "Relative biomass index", shape = "Region sampled") +
  theme(legend.position = "bottom")
ggsave(here::here("figs/prior-sensitivity-indexes.pdf"), g_indexes, width = 7, height = 5)
