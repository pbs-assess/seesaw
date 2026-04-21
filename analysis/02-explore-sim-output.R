out_df <- readRDS("data-generated/sawtooth-sim-apr20.rds")

# out_df$label <- forcats::fct_inorder(out_df$label)

cols <- RColorBrewer::brewer.pal(3L, "Set2")
names(cols) <- c("north", "south", "both")

seed_to_plot <- 1
actual <- select(out_df, label, year, total, seed, sampled_region) |>
  filter(seed == seed_to_plot) |>
  distinct()
actual2 <- mutate(actual, label = gsub("obs", "\\\nobs", label)) |>
  mutate(label = gsub("black", "\\\nblack", label)) |>
  mutate(label = gsub("mixture", "\\\nmixture", label)) |>
  mutate(label = gsub("\\(", "\\\n\\(", label))

g <- out_df |>
  filter(seed == seed_to_plot) |>
  mutate(with_depth = gsub("covariate =", "cov =", with_depth)) |>
  mutate(type_facet = gsub(",\\s*", ",\n", type)) |>
  mutate(label = gsub("obs", "\\\nobs", label)) |>
  mutate(label = gsub("black", "\\\nblack", label)) |>
  mutate(label = gsub("mixture", "\\\nmixture", label)) |>
  mutate(label = gsub("\\(", "\\\n\\(", label)) |>
  ggplot(aes(year, est, ymin = lwr, ymax = upr)) +
  ggsidekick::theme_sleek() +
  geom_pointrange(aes(colour = sampled_region)) +
  # geom_line(colour = "grey50") +
  geom_ribbon(alpha = 0.20, colour = NA) +
  geom_line(
    data = actual2, mapping = aes(year, total),
    inherit.aes = FALSE, lty = 2
  ) +
  facet_grid(forcats::fct_inorder(label) ~ type_facet,
    scales = "free_y"
  ) +
  ylab("Abundance estimate") +
  xlab("Year") +
  labs(colour = "Sampled region") +
  scale_colour_manual(values = cols[c(2, 1, 3)]) +
  scale_y_log10() +
  scale_x_continuous(breaks = function(x) seq(ceiling(x[1]), floor(x[2]), by = 2))
# print(g)
ggsave("figs/saw-tooth-scenarios-2026-04-20.pdf", width = 28, height = 23)

# a minimal version for the main text:
seed_to_plot <- 3
actual <- select(out_df, label, year, total, seed, sampled_region) |>
  filter(seed == seed_to_plot) |>
  distinct()
actual2 <- mutate(actual, label = gsub("obs", "\\\nobs", label)) |>
  mutate(label = gsub("black", "\\\nblack", label)) |>
  mutate(label = gsub("mixture", "\\\nmixture", label)) |>
  mutate(label = gsub("\\(", "\\\n\\(", label))

labels <- c("Base", "No gap", "Large gap", "High range", "Low range")

.actual2 <- dplyr::filter(actual2, label %in% labels)
SCALER <- 1e4
out_df |>
  filter(seed == seed_to_plot) |>
  filter(label %in% labels, with_depth == "covariate = FALSE", model %in% c("RW RF", "IID RF, factor(year)", "IID RF, RW year")) |>
  mutate(with_depth = gsub("covariate =", "cov =", with_depth)) |>
  mutate(label = gsub("obs", "\\\nobs", label)) |>
  mutate(label = gsub("black", "\\\nblack", label)) |>
  mutate(label = gsub("mixture", "\\\nmixture", label)) |>
  mutate(label = gsub("\\(", "\\\n\\(", label)) |>
  ggplot(aes(year, est/SCALER, ymin = lwr/SCALER, ymax = upr/SCALER)) +
  ggsidekick::theme_sleek() +
  geom_pointrange(aes(colour = sampled_region)) +
  # geom_line(colour = "grey50") +
  geom_ribbon(alpha = 0.20, colour = NA) +
  geom_line(
    data = .actual2, mapping = aes(year, total/SCALER),
    inherit.aes = FALSE, lty = 2
  ) +
  facet_grid(forcats::fct_inorder(label) ~ type,
    scales = "free_y"
  ) +
  ylab("Abundance estimate") +
  xlab("Year") +
  labs(colour = "Sampled region") +
  scale_colour_manual(values = cols[c(2, 1, 3)]) +
  scale_y_log10() +
  coord_cartesian(ylim = c(10, 1000)) +
  scale_x_continuous(breaks = function(x) seq(ceiling(x[1]), floor(x[2]), by = 2))

ggsave("figs/example-indices-simulated.pdf", width = 7.5, height = 6.5)


# Look at one point in space... -------------------------------------------

# get_eg_cell <- function(obj, x, y) {
#   obj$data[
#     round(obj$data$X, 3) == round(x, 3) &
#       round(obj$data$Y, 3) == round(y, 3),
#   ]
# }
#
# p1 <- purrr::map_dfr(preds, get_eg_cell,
#   x = 0.50505051, y = 0.81818182,
#   .id = "model"
# ) |>
#   mutate(with_depth = grepl("covariate", model)) |>
#   mutate(type = gsub(" covariate", "", model))
#
# p1 |>
#   left_join(select(d, year, sampled_region) %>% distinct()) |>
#   ggplot(aes(year, est)) +
#   geom_line() +
#   facet_grid(with_depth ~ type) +
#   geom_point(aes(colour = sampled_region)) +
#   ggsidekick::theme_sleek()

# What about MRE, RMSE, see-saw, coverage etc. ? --------------------------

out_df |>
  group_by(model) |>
  mutate(log_residual = log(total) - log(est)) |>
  summarise(
    seesaw_index = abs(mean(log_residual[sampled_region == "north"]) -
      mean(log_residual[sampled_region == "south"])),
    mre = mean(log_residual),
    rmse = sqrt(mean(log_residual^2)),
    mean_se = mean(se),
    coverage = mean(total < upr & total > lwr)
  ) |>
  arrange(seesaw_index, rmse) |>
  mutate(model = forcats::fct_reorder(model, rev(seesaw_index))) |>
  tidyr::pivot_longer(cols = -model, names_to = "metric") |>
  mutate(metric = factor(metric,
    levels = c("seesaw_index", "rmse", "mre", "mean_se", "coverage")
  )) |>
  ggplot(aes(value, model)) +
  geom_point(pch = 21) +
  facet_wrap(~metric, scales = "free_x", nrow = 1L) +
  ggsidekick::theme_sleek() +
  theme(panel.grid.major.y = element_line(colour = "grey90"), axis.title.y.left = element_blank()) +
  xlab("Metric value")
ggsave("figs/dot-plot-metrics.pdf", width = 10, height = 4)

# ---------------------
# get at distribution

# assign an indicator for every 2nd year:
lu <- data.frame(year = sort(unique(out_df$year)))
is_odd <- function(x) x %% 2 != 0
lu$odd <- is_odd(lu$year)

temp <- out_df |>
  left_join(lu) |> 
  group_by(label, seed, model) |>
  mutate(log_residual = log(total) - log(est)) |>
  summarise(
    seesaw_index = abs(mean(log_residual[odd]) - mean(log_residual[!odd])),
    mre = mean(log_residual),
    rmse = sqrt(mean(log_residual^2)),
    mean_se = mean(se),
    coverage = mean(total < upr & total > lwr)
  ) |>
  tidyr::pivot_longer(cols = -c(model, seed, label), names_to = "metric") |>
  mutate(metric = factor(metric,
    levels = c("seesaw_index", "rmse", "mre", "mean_se", "coverage")
  )) |>
  group_by(label, model, metric) |>
  summarise(lwr = quantile(value, 0.2, na.rm = TRUE), upr = quantile(value, 0.8, na.rm = TRUE), med = quantile(value, 0.5, na.rm = TRUE))

saw_tooth_ind <- temp |> filter(metric == "seesaw_index") |>
  group_by(model) |>
  summarise(med_st_index = -median(med))

g <- temp |>
  left_join(saw_tooth_ind) |>
  filter(metric != "mre") |>
  mutate(label = gsub("obs", "\\\nobs", label)) |>
  mutate(label = gsub("black", "\\\nblack", label)) |>
  mutate(label = gsub("mixture", "\\\nmixture", label)) |>
  mutate(label = gsub("\\(", "\\\n\\(", label)) |>
  ggplot(aes(med, forcats::fct_reorder(model, med_st_index))) +
  geom_point(pch = 21) +
  geom_linerange(aes(xmin = lwr, xmax = upr)) +
  facet_grid(forcats::fct_inorder(label)~metric, scales = "free_x") +
  ggsidekick::theme_sleek() +
  theme(panel.grid.major.y = element_line(colour = "grey90"), axis.title.y.left = element_blank()) +
  xlab("Metric value")
g
ggsave("figs/saw-tooth-metrics-all-apr20.pdf", width = 7, height = 45)

# ---------------------------------------------------------------------
# together in one set of panels?

g <- temp |>
  left_join(saw_tooth_ind) |>
  filter(metric != "mre") |>
  mutate(metric = as.character(metric)) |> 
  mutate(metric = replace_values(metric,
    "seesaw_index" ~ "Seesaw metric",
    "coverage" ~ "CI coverage",
    "mean_se" ~ "Mean SE",
    "rmse" ~ "RMSE"
  )) |> 
  mutate(metric = factor(metric,
    levels = c("Seesaw metric", "RMSE", "Mean SE", "CI coverage")
  )) |>
  mutate(label = gsub("obs", "\\\nobs", label)) |>
  filter(!model %in% c("SVC trend, spatial only", "SVC trend, IID fields")) |> 
  ggplot(aes(x = forcats::fct_reorder(model, med_st_index), y = med, group = label)) +
  geom_point(pch = 21, position = position_dodge(width = 0.2), alpha = 0.5) +
  # geom_linerange(aes(ymin = lwr, ymax = upr), position = position_dodge(width = 0.5)) +
  facet_wrap(~metric, scales = "free_x", nrow = 1) +
  ggsidekick::theme_sleek() +
  theme(panel.grid.major.y = element_line(colour = "grey90"), axis.title.y = element_blank()) +
  ylab("Metric value") +
  coord_flip()
g
ggsave("figs/saw-tooth-metrics-condensed-apr20.pdf", width = 7, height = 4)

# ----------------------

temp <- out_df |>
  left_join(lu) |> 
  group_by(label, model, seed) |>
  mutate(log_residual = log(total) - log(est)) |>
  summarise(
    seesaw_index = abs(mean(log_residual[odd]) - mean(log_residual[!odd]))
  ) |>
  group_by(label, model) |>
  summarize(
    lwr = quantile(seesaw_index, 0.2, na.rm=TRUE),
    upr = quantile(seesaw_index, 0.8, na.rm=TRUE),
    med = median(seesaw_index, na.rm=TRUE)
  ) |>
  ungroup() |>
  filter(model == "IID RF, factor(year)")

temp |>
  ggplot(aes(med, forcats::fct_reorder(label, med))) +
  geom_vline(xintercept = temp$med[temp$model == "IID" & temp$label == "Base"], lty = 2, colour = "grey70") +
  geom_linerange(aes(xmin = lwr, xmax = upr)) +
  geom_point(pch = 21, size = 1.8) +
  # facet_wrap(~forcats::fct_inorder(label), scales = "free_x") +
  ggsidekick::theme_sleek() +
  theme(panel.grid.major.y = element_line(colour = "grey95"), axis.title.y.left = element_blank()) +
  xlab("Seesaw metric")
ggsave("figs/saw-tooth-bad-iid-apr20.pdf", width = 4.2, height = 4.5)

# Figures
# Spatial setup example
# time series example(s) from one seed?
# summary of when bad (dot-line plot)
# metrics across a variety of scenarios
# maybe a figure digging into an example where it goes wrong - show year effects, show spatial omegas and or fixed effect covariate predictions



# observation: even with covariate correctly fixed, seems to want to put
# variance into year factors and shrink random fields
# nothing stopping it from see-sawing fixed effects and shrinking random field
# variance

# Question: can an extra latent year at the beginning fix the 'burnin' problem!?

# Observations ------------------------------------------------------------

# - Try lower or increasing observation error ('phi') (e.g., 2 vs. 7)
# - Try lowering or increasing numbers of samples (e.g., 250 vs. 400)

# - Is it important to have at least one year with full coverage?
# - If so, does it matter if it's near the beginning or end?
# - Is it/when is it important to have some kind of fixed effect covariate (region/depth)?
# - What form (IID, RW, AR1) works best in what conditions?
# - Does an AR1 see-saw towards mean with uneven coverage?
# - When is RW at risk of over smoothing process?
