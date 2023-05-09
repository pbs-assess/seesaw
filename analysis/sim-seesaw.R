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

# study design:
# - no gap, small gap, large gap
# - no covariate effect; strong covariate effect
# - lower or increasing observation error ('phi') (e.g., 2 vs. 7)
# - try lowering or increasing numbers of samples (e.g., 250 vs. 400)
# - one year of complete coverage at beginning or end or none
# - impact of uneven north south vs. even?

# - what about ar1 on year effect? add to sdmTMB?
# - RW on year effect?

# quadratic, linear, or breakpoint covariate effect?

library(sdmTMB)
library(ggplot2)
library(dplyr)
source("stitch/funcs.R")
dir.create("stitch/figs")
# Simulation testing survey stitching with various models -----------------

N_YEAR <- 12 # even number
SEED <- 28817
SEED <- 288171

predictor_grid <- expand.grid(
  X = seq(0, 1, length.out = 100), Y = seq(0, 1, length.out = 100),
  year = seq_len(N_YEAR)
)
mesh_sim <- make_mesh(predictor_grid, xy_cols = c("X", "Y"), cutoff = 0.1)

xx <- seq(0, 1, length.out = 100)
plot(xx, xx * 2 + 5 * xx^2, ylab = "Depth effect", xlab = "Depth value")

x <- sim(predictor_grid, mesh_sim, seed = SEED,
  year_arima.sim = list(ar = 0.6), year_marginal_sd = 0.3,
  coefs = c(2, 5),
  north_effect = 0)
sim_dat <- x$sim_dat
predictor_dat <- x$predictor_dat

# Visualize what we just did ----------------------------------------------

# Year effects:
ggplot(data.frame(x = seq_len(N_YEAR), y = x$year_effects), aes(x, y)) +
  geom_line()
ggplot(data.frame(x = seq_len(N_YEAR), y = x$year_effects), aes(x, exp(y))) +
  geom_line()

blank_theme_elements <- theme(panel.grid.major = element_line(colour = "grey90"),
  # panel.spacing.x = unit(20, "pt"),
  axis.text = element_blank(), axis.ticks = element_blank(),
  axis.title = element_blank(), legend.position = "right")

# Spatio-temporal truth
ggplot(sim_dat, aes(X, Y, fill = eta)) +
  geom_raster() +
  facet_wrap(vars(year)) +
  scale_fill_viridis_c() +
  ggsidekick::theme_sleek() +
  coord_equal(expand = FALSE) +
  blank_theme_elements +
  labs(fill = "True\nsimulated\nlog abundance")
# ggsave("stitch/figs/spatio-temporal-truth.pdf", width = 8, height = 5)

# Sample N per year -----------------------------------------------------

d <- observe(sim_dat, sample_n = 400, seed = SEED,
  north_yrs = seq(1, N_YEAR - 1, 2), south_yrs = seq(2, N_YEAR, 2), gap = 0.3)

# Visualize it ------------------------------------------------------------

ggplot(d, aes(X, Y, colour = log(observed))) +
  geom_point() +
  facet_wrap(vars(year)) +
  scale_colour_viridis_c() +
  ggsidekick::theme_sleek() +
  blank_theme_elements +
  coord_equal(expand = FALSE) +
  theme(panel.grid.major = element_line(colour = "grey90"))
# panel.spacing.x = unit(20, "pt")
# ggsave("stitch/figs/spatio-temporal-observed.pdf", width = 8, height = 5)

# Calculate known true biomass/abundance ----------------------------------

actual <- group_by(sim_dat, year) %>%
  summarise(total = sum(mu))
ggplot(actual, aes(year, total)) +
  geom_line()

# Fit models --------------------------------------------------------------

mesh <- make_mesh(d, c("X", "Y"), cutoff = 0.1)
priors <- sdmTMBpriors(
  matern_s = pc_matern(range_gt = 0.3, sigma_lt = 0.4),
  matern_st = pc_matern(range_gt = 0.3, sigma_lt = 0.3)
)

fits <- list()
nms <- c()
i <- 1

fits[[i]] <- sdmTMB(
  observed ~ 1 + depth_cov + I(depth_cov^2), family = tweedie(),
  # observed ~ 1 + depth_cov, family = tweedie(),
  data = d, time = "year", spatiotemporal = "rw", spatial = "on",
  silent = TRUE, mesh = mesh,
  priors = priors
)
nms <- c(nms, "RW covariate")
i <- i + 1

fits[[i]] <- sdmTMB(
  observed ~ 1, family = tweedie(),
  data = d, time = "year", spatiotemporal = "rw", spatial = "on",
  silent = TRUE, mesh = mesh,
  # extra_time = 0L,
  priors = priors
)
nms <- c(nms, "RW")
i <- i + 1

fits[[i]] <- sdmTMB(
  # observed ~ 0 + as.factor(year) + depth_cov + I(depth_cov^2), family = tweedie(),
  observed ~ 0 + as.factor(year) + depth_cov, family = tweedie(),
  data = d, time = "year", spatiotemporal = "iid", spatial = "on",
  silent = TRUE, mesh = mesh,
  priors = priors
)
nms <- c(nms, "IID covariate")
i <- i + 1

# fits[[i]] <- sdmTMB(
#   observed ~ 0 + as.factor(year) + depth_cov + I(depth_cov^2), family = tweedie(),
#   data = d, time = "year", spatiotemporal = "iid", spatial = "on",
#   silent = TRUE, mesh = mesh,
#   priors = priors,
#   control = sdmTMBcontrol(
#     start = list(b_j = c(rep(0, 10), 2.2, 3.8)),
#     map = list(b_j = factor(c(seq(1, 10), NA, NA)))
#   )
# )
# nms <- c(nms, "IID covariate fixed")
# i <- i + 1

fits[[i]] <- sdmTMB(
  observed ~ s(year), family = tweedie(),
  data = d, time = "year", spatiotemporal = "iid", spatial = "on",
  silent = TRUE, mesh = mesh,
  priors = priors
)
nms <- c(nms, "IID s(year)")
i <- i + 1

fits[[i]] <- sdmTMB(
  # observed ~ s(year) + depth_cov + I(depth_cov^2), family = tweedie(),
  observed ~ s(year) + depth_cov, family = tweedie(),
  data = d, time = "year", spatiotemporal = "iid", spatial = "on",
  silent = TRUE, mesh = mesh,
  priors = priors
)
nms <- c(nms, "IID covariate s(year)")
i <- i + 1

fits[[i]] <- sdmTMB(
  observed ~ 0 + as.factor(year), family = tweedie(),
  data = d, time = "year", spatiotemporal = "iid", spatial = "on",
  mesh = mesh,
  priors = priors
)
i <- i + 1
nms <- c(nms, "IID")

fits[[i]] <- sdmTMB(
  observed ~ 0, family = tweedie(),
  time_varying = ~ 1,
  data = d, time = "year", spatiotemporal = "iid", spatial = "on",
  mesh = mesh,
  priors = priors
)
i <- i + 1
nms <- c(nms, "IID RW year")

# fits[[i]] <- sdmTMB(
#   observed ~ 0, family = tweedie(),
#   time_varying = ~ 1, time_varying_type = "ar1",
#   data = d, time = "year", spatiotemporal = "iid", spatial = "on",
#   mesh = mesh,
#   priors = priors
# )
# i <- i + 1
# nms <- c(nms, "IID AR1 year")

fits[[i]] <- sdmTMB(
  # observed ~ 0 + depth_cov + I(depth_cov^2), family = tweedie(),
  observed ~ 0 + depth_cov, family = tweedie(),
  time_varying = ~ 1,
  data = d, time = "year", spatiotemporal = "iid", spatial = "on",
  mesh = mesh,
  priors = priors
)
i <- i + 1
nms <- c(nms, "IID RW covariate year")

fits[[i]] <- sdmTMB(
  # observed ~ 0 + as.factor(year) + depth_cov + I(depth_cov^2),
  observed ~ 0 + as.factor(year) + depth_cov,
  family = tweedie(),
  data = d, time = "year", spatiotemporal = "off", spatial = "on",
  mesh = mesh,
  priors = priors
)
i <- i + 1
nms <- c(nms, "Spatial only covariate")

fits[[i]] <- sdmTMB(
  observed ~ 0 + as.factor(year),
  family = tweedie(),
  data = d, time = "year", spatiotemporal = "off", spatial = "on",
  mesh = mesh,
  priors = priors
)
i <- i + 1
nms <- c(nms, "Spatial only")

names(fits) <- nms

# Predict on grid and calculate indexes -----------------------------------

nd <- select(predictor_dat, X, Y, year, region)
nd$depth_cov <- nd$Y
preds <- purrr::map(fits, predict, newdata = nd, return_tmb_object = TRUE)
indexes <- purrr::map(preds, get_index, bias_correct = TRUE)

indexes_df <- dplyr::bind_rows(indexes, .id = "model") |>
  mutate(with_depth = paste0("covariate = ", grepl("covariate", model))) |>
  mutate(type = gsub(" covariate", "", model))

# Plot it -----------------------------------------------------------------

ylims <- range(indexes_df$est) * c(0.8, 1.2)
mult <- 1
g <- indexes_df |>
  left_join(select(d, year, sampled_region) %>% distinct()) %>%
  ggplot(aes(year, est / mult,
  ymin = lwr / mult, ymax = upr / mult)) +
  ggsidekick::theme_sleek() +
  geom_pointrange(aes(colour = sampled_region)) +
  geom_line() +
  geom_ribbon(alpha = 0.2, colour = NA) +
  geom_line(data = actual, mapping = aes(year, total),
    inherit.aes = FALSE, lty = 2) +
  facet_grid(with_depth~type) +
  # ggtitle("IID vs. RW vs. spatial model; with and without covariate",
  # subtitle = paste0("Dashed = true; dots/lines = estimated\n", "phi = ", PHI, ", N = ", SAMPLE_N)) +
  # subtitle = paste0("Dashed = true; dots/lines = estimated\n")) +
  ylab("Abundance estimate") + xlab("Year") +
  labs(colour = "Sampled region") +
  coord_cartesian(ylim = ylims) +
  scale_y_log10() +
  scale_x_continuous(breaks = function(x) seq(ceiling(x[1]), floor(x[2]), by = 2))
print(g)
# ggsave("stitch/figs/saw-tooth-example.pdf", width = 8, height = 5)

# Look at one point in space... -------------------------------------------

get_eg_cell <- function(obj, x, y) {
  obj$data[
    round(obj$data$X, 3) == round(x, 3) &
    round(obj$data$Y, 3) == round(y, 3), ]
}

p1 <- purrr::map_dfr(preds, get_eg_cell, x = 0.50505051, y = 0.81818182,
  .id = "model") |>
  mutate(with_depth = grepl("covariate", model)) |>
  mutate(type = gsub(" covariate", "", model))

p1 |> left_join(select(d, year, sampled_region) %>% distinct()) |>
  ggplot(aes(year, est)) + geom_line() +
  facet_grid(with_depth~type) +
  geom_point(aes(colour = sampled_region)) +
  ggsidekick::theme_sleek()

# What about MRE, RMSE, see-saw, coverage etc. ? --------------------------

indexes_df |>
  left_join(actual) |>
  left_join(select(d, year, sampled_region) %>% distinct()) |>
  group_by(model) |>
  mutate(log_residual = log(total) - log(est)) |>
  summarise(
    seesaw_index = abs(mean(log_residual[sampled_region == "north"]) - mean(log_residual[sampled_region == "south"])),
    mre = mean(log_residual),
    rmse = sqrt(mean(log_residual^2)),
    mean_se = mean(se),
    coverage = mean(total < upr & total > lwr)
  ) |>
  arrange(seesaw_index, rmse) |>
  mutate(model = forcats::fct_reorder(model, rev(seesaw_index))) |>
  tidyr::pivot_longer(cols = -model, names_to = "metric") |>
  mutate(metric = factor(metric,
    levels = c("seesaw_index", "rmse", "mre", "mean_se", "coverage"))) |>
  ggplot(aes(value, model)) +
  geom_point(pch = 21, size = 1.6) +
  facet_wrap(~metric, scales = "free_x", nrow = 1L) +
  ggsidekick::theme_sleek() +
  theme(panel.grid.major.y = element_line(colour = "grey90"), axis.title.y.left = element_blank()) +
  xlab("Metric value")

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
