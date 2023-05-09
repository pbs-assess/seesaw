# question: same thing happen with missing chunks from one index!?

# To look at:
# - gap matters! how much; related to range I assume
# - **why does RW 'fix' things!?**
# - getting the depth covariate right seems to fix seesaw; *but* bigger CIs on
#   poorly sampled region
# - is coverage of RW still OK!?
# - happens with spatial-only model I think (from Quang)
# - dig in... need to articulate simply what's going wrong!
# - how do I know this is happening in reality (besides the seesaw)?
# - does including a region/north/south covariate fix it: No! build in to
#   routine code?

# What is an index of seesawness?
#
# Things to show:
# - spatial random field messing up
# - RW deviations adapting to fit?
# - fix depth covariates correctly - works?

# quadratic, linear, or breakpoint covariate effect?

library(sdmTMB)
library(ggplot2)
library(dplyr)
source("stitch/funcs.R")
dir.create("stitch/figs", showWarnings = FALSE)
# Simulation testing survey stitching with various models -----------------

source(here::here("stitch/scenarios.R"))
if (any(grepl("empty", purrr::map_chr(sc, "label")))) stop("Too many slots")
labels <- unname(purrr::map_chr(sc, "label"))
names(sc) <- labels
categories <- unname(purrr::map_chr(sc, "category"))
lu <- data.frame(label = labels, category = categories, stringsAsFactors = FALSE)
sc <- purrr::map(sc, ~ {
  .x$label <- NULL
  .x$category <- NULL
  .x
})

# tictoc::tic()
# out20 <- do.call(sim_fit_and_index, c(sc[[11]], .seed = 1, make_plots = FALSE))
# tictoc::toc()


# testing first:
out <- do.call(sim_fit_and_index, c(sc[[1]], svc_trend = -1, .seed = 1, make_plots = T))
out <- do.call(sim_fit_and_index, c(sc[[1]], svc_trend = -0.6, .seed = 2, make_plots = F))

actual <- select(out, year, total, seed, sampled_region) |>
  distinct()

ggplot(out, aes(year, est, ymin = lwr, ymax = upr)) +
  ggsidekick::theme_sleek() +
  geom_pointrange(aes(colour = sampled_region)) +
  # geom_line(colour = "grey50") +
  geom_ribbon(alpha = 0.20, colour = NA) +
  geom_line(
    data = actual, mapping = aes(year, total),
    inherit.aes = FALSE, lty = 2
  ) +
  facet_wrap( ~ paste(type, with_depth),
    scales = "free_y"
  )

future::plan(future::multisession, workers = 8L)
tictoc::tic()
seeds <- seq_len(16L)
# out_df <- purrr:::map_dfr(seeds[1], function(seed_i) {
out_df <- furrr::future_map_dfr(seeds, function(seed_i) {
  x <- list()
  for (i in seq_along(sc)) {
    sc[[i]]$.seed <- seed_i
    x[[i]] <- do.call(sim_fit_and_index, sc[[i]])
  }
  names(x) <- names(sc)
  bind_rows(x, .id = "label")
# })
}, .options = furrr::furrr_options(seed = TRUE))
tictoc::toc()
out_df2 <- left_join(out_df, lu, by = "label")
saveRDS(out_df2, "stitch/sawtooth-sim-may3.rds")
future::plan(future::sequential)

out_df <- readRDS("stitch/sawtooth-sim-may3.rds")

# out_df$label <- forcats::fct_inorder(out_df$label)

cols <- RColorBrewer::brewer.pal(3L, "Set2")
names(cols) <- c("north", "south", "both")

seed_to_plot <- 4
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
  facet_grid(forcats::fct_inorder(label) ~ paste(type, with_depth),
    scales = "free_y"
  ) +
  ylab("Abundance estimate") +
  xlab("Year") +
  labs(colour = "Sampled region") +
  scale_colour_manual(values = cols[c(2, 1, 3)]) +
  scale_y_log10() +
  scale_x_continuous(breaks = function(x) seq(ceiling(x[1]), floor(x[2]), by = 2))
# print(g)
ggsave("stitch/figs/saw-tooth-scenarios-may3-4.pdf", width = 15, height = 24)

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
  geom_point(pch = 21, size = 1.6) +
  facet_wrap(~metric, scales = "free_x", nrow = 1L) +
  ggsidekick::theme_sleek() +
  theme(panel.grid.major.y = element_line(colour = "grey90"), axis.title.y.left = element_blank()) +
  xlab("Metric value")

# ---------------------
# get at distribution

temp <- out_df |>
  group_by(label, seed, model) |>
  mutate(log_residual = log(total) - log(est)) |>
  summarise(
    seesaw_index1 = abs(mean(log_residual[sampled_region == "north" & year < 7]) -
        mean(log_residual[sampled_region == "south" & year < 7])),
    seesaw_index2 = abs(mean(log_residual[sampled_region == "north" & year >= 7]) -
        mean(log_residual[sampled_region == "south" & year >= 7])),
    mre = mean(log_residual),
    rmse = sqrt(mean(log_residual^2)),
    mean_se = mean(se),
    coverage = mean(total < upr & total > lwr)
  ) |>
  mutate(seesaw_index = (seesaw_index1 + seesaw_index2) / 2) |>
  tidyr::pivot_longer(cols = -c(model, seed, label), names_to = "metric") |>
  mutate(metric = factor(metric,
    levels = c("seesaw_index", "rmse", "mre", "mean_se", "coverage")
  )) |>
  group_by(label, model, metric) |>
  summarise(lwr = quantile(value, 0.2), upr = quantile(value, 0.8), med = quantile(value, 0.5))

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
ggsave("stitch/figs/saw-tooth-metrics-all-May3.pdf", width = 7, height = 25)

# ---------------------------------------------------------------------
# together in one set of panels?

g <- temp |>
  left_join(saw_tooth_ind) |>
  filter(metric != "mre") |>
  mutate(label = gsub("obs", "\\\nobs", label)) |>
  # ggplot(aes(x = forcats::fct_reorder(model, med_st_index), y = med, colour = label)) +
  ggplot(aes(x = forcats::fct_reorder(model, med_st_index), y = med, group = label)) +
  geom_point(pch = 21, position = position_dodge(width = 0.2), alpha = 0.8) +
  # geom_linerange(aes(ymin = lwr, ymax = upr), position = position_dodge(width = 0.5)) +
  facet_wrap(~metric, scales = "free_x", nrow = 1) +
  ggsidekick::theme_sleek() +
  theme(panel.grid.major.y = element_line(colour = "grey90"), axis.title.y = element_blank()) +
  ylab("Metric value") +
  coord_flip()
g
ggsave("stitch/figs/saw-tooth-metrics-condensed-dec14.pdf", width = 6.8, height = 3)


# ----------------------

temp <- out_df |>
  group_by(label, model, seed) |>
  mutate(log_residual = log(total) - log(est)) |>
  summarise(
    seesaw_index = abs(mean(log_residual[sampled_region == "north"]) -
      mean(log_residual[sampled_region == "south"]))
  ) |>
  group_by(label, model) |>
  summarize(
    lwr = quantile(seesaw_index, 0.2),
    upr = quantile(seesaw_index, 0.8),
    med = median(seesaw_index)
  ) |>
  ungroup() |>
  filter(model == "IID")

temp |>
  ggplot(aes(med, forcats::fct_reorder(label, med))) +
  geom_vline(xintercept = temp$med[temp$model == "IID" & temp$label == "Base"], lty = 2, colour = "grey70") +
  geom_linerange(aes(xmin = lwr, xmax = upr)) +
  geom_point(pch = 21, size = 1.8) +
  # facet_wrap(~forcats::fct_inorder(label), scales = "free_x") +
  ggsidekick::theme_sleek() +
  theme(panel.grid.major.y = element_line(colour = "grey95"), axis.title.y.left = element_blank()) +
  xlab("Seesaw metric")
ggsave("stitch/figs/saw-tooth-bad-iid-dec14.pdf", width = 4.2, height = 4.5)

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
