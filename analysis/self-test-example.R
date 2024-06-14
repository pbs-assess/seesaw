# - simulate with new random fields with full vs. the actual survey coverage
# - do they substantially differ? (maybe chi-squared test idea from Rufener et al.)
# - do with underlying IID and with underlying RW... must be self consistent

library(dplyr)
library(ggplot2)
theme_set(ggsidekick::theme_sleek() + theme(tagger.panel.tag.text = element_text(color = "grey30", size = 9)))
library(sdmTMB)

# bring in set and grid data ------------------------------------

# CASE <- "quillback inside"
# SURVEY <- "inside"

CASE <- "lingcod synoptic"
SURVEY <- "synoptic"

if (CASE == "quillback inside") {
  d <- readRDS("~/src/gfsynopsis-2022/report/data-cache-nov-2023/quillback-rockfish.rds")$survey_sets
  d <- filter(d, survey_abbrev %in% c("HBLL INS N", "HBLL INS S"))
  # simplify, remove N in 2021 (was N + S)
  d <- filter(d, !(year == 2021 & survey_abbrev %in% "HBLL INS N"))
  d <- mutate(d, response = catch_count, log_effort = log(hook_count))
  d$biennial_region <- as.integer(as.factor(d$survey_abbrev))
  family <- "nb"
}

if (SURVEY == "inside") {
  grid1 <- gfplot::hbll_inside_n_grid$grid
  grid2 <- gfplot::hbll_inside_s_grid$grid
  grid <- bind_rows(grid1, grid2)
  grid <- rename(grid, lon = X, lat = Y) |>
    add_utm_columns(ll_names = c("lon", "lat"))
}

if (CASE == "lingcod synoptic") {
  d <- readRDS("~/src/gfsynopsis-2023/report/data-cache-2024-05/lingcod.rds")$survey_sets
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
  table(d$year, d$survey_abbrev)
  d <- filter(d, !year %in% c(2003, 2004, 2020))
  d <- filter(d, !(year %in% c(2021) & survey_abbrev == "SYN WCVI"))
  d <- filter(d, !(year %in% c(2007) & survey_abbrev == "SYN WCHG"))
  d$biennial_region <- ifelse(d$survey_abbrev %in% c("SYN WCVI", "SYN WCHG"), 1L, 2L)
  table(d$year, d$survey_abbrev)
  table(d$year, d$biennial_region)

  d <- mutate(d, observed = density_kgkm2, log_effort = log(area_swept))
  family <- tweedie()
}

if (SURVEY == "synoptic") {
  grid <- gfplot::synoptic_grid |> select(-survey_domain_year, -utm_zone)
}

# fit initial models --------------------------------------------

d <- add_utm_columns(d)
mesh <- make_mesh(d, c("X", "Y"), cutoff = 15)
plot(mesh)

fit <- sdmTMB(
  observed ~1,
  offset = "log_effort",
  data = d,
  time = "year",
  spatial = "on",
  time_varying = ~1,
  time_varying_type = "rw0",
  priors = sdmTMBpriors(sigma_V = gamma_cv(0.2, 0.5)),
  control = sdmTMBcontrol(multiphase = FALSE),
  spatiotemporal = "rw",
  mesh = mesh,
  silent = FALSE,
  family = family
)
fit
sanity(fit)

# fitii <- update(
#   fit,
#   spatiotemporal = "iid",
#   formula. = observed ~ 0 + factor(year)
# )

fitii <- sdmTMB(
  observed ~ 0 + factor(year),
  offset = "log_effort",
  data = d,
  time = "year",
  spatial = "on",
  spatiotemporal = "iid",
  mesh = mesh,
  silent = FALSE,
  family = family
)

dy <- group_by(d, year) |>
  summarise(biennial_region = factor(biennial_region[1]))

nd <- replicate_df(grid, "year", unique(d$year))

p <- predict(fit, newdata = nd, return_tmb_object = TRUE)
pii <- predict(fitii, newdata = nd, return_tmb_object = TRUE)
ind <- get_index(p, bias_correct = TRUE)
indii <- get_index(pii, bias_correct = TRUE)
ind <- left_join(ind, dy)
indii <- left_join(indii, dy)

ggplot(ind, aes(year, est, ymin = lwr, ymax = upr, colour = factor(biennial_region))) +
  geom_pointrange()
ggplot(indii, aes(year, est, ymin = lwr, ymax = upr, colour = factor(biennial_region))) +
  geom_pointrange()

b1 <- tidy(fit)
b2 <- tidy(fit, "ran_pars")

omega <- get_pars(fit)$omega_s
eps <- get_pars(fit)$epsilon_st

# want to create fake sampling in the unsampled region in a given year
# want to see if IID and RW give plausibly similar indexes in partial vs. full sampling
# to create these fake sampling locations, quick hack is to take sampling
# locations randomly from another sampled year

group_by(d, year) |>
  summarise(n = n())

get_random_year_samples <-
  function(data, region = c(1, 2)) {
    dd <- filter(data, biennial_region %in% region)
    yrs <- unique(dd$year)
    yr <- sample(yrs, 1L)
    filter(dd, year %in% yr) |>
      select(X, Y)
  }

set.seed(1)
locs <- group_by(d, year) |>
  group_split() |>
  purrr::map(\(.x) {
    real_loc <- select(.x, X, Y, year) |> mutate(type = "real")
    region <- unique(.x$biennial_region)
    # other_region <- if (region == "HBLL INS N") "HBLL INS S" else "HBLL INS N"
    other_region <- if (region == 1) 2 else 1
    other_loc <- get_random_year_samples(d, other_region) |>
      mutate(type = "fake") |>
      mutate(year = unique(.x$year))
    bind_rows(real_loc, other_loc)
  }) |>
  bind_rows()

ggplot(d, aes(X, Y)) +
  geom_point() +
  facet_wrap(~year)

ggplot(locs, aes(X, Y, colour = type)) +
  geom_point(alpha = 0.8, pch = ".") +
  facet_wrap(~year)

locs_combined <- bind_rows(locs, select(nd, X, Y, year) |> mutate(type = "grid"))

meshs <- make_mesh(locs_combined, c("X", "Y"), mesh = mesh$mesh)
plot(meshs$mesh)
meshs$mesh$n

# set.seed(123)
# year_effects <- numeric(length(unique(d$year)))
# (sigma_V <- as.numeric(fit$tmb_obj$report()$sigma_V))
# year_effects[1] <- unname(coef(fit)) + rnorm(1, 0, sigma_V)
# for (i in 2:length(year_effects)) {
#   year_effects[i] <- year_effects[i-1] + rnorm(1, 0, sigma_V)
# }
# plot(year_effects)

e <- as.list(fit$sd_report, "Estimate")
year_effects <- e$b_rw_t[,,1]
plot(exp(year_effects))

simf <- sdmTMB_simulate(
  ~ 0 + as.factor(year),
  data = locs_combined,
  sigma_O = b2$estimate[b2$term == "sigma_O"],
  sigma_E = b2$estimate[b2$term == "sigma_E"],
  phi = b2$estimate[b2$term == "phi"],
  tweedie_p = b2$estimate[b2$term == "tweedie_p"],
  range = b2$estimate[b2$term == "range"],
  fixed_re = list(omega_s = omega),
  spatiotemporal = "ar1",
  rho = 0.9,
  time = "year",
  seed = 123,
  family = family,
  mesh = meshs,
  # offset = 1,
  B = year_effects
)

simf$type <- locs_combined$type
simb <- simf
simb <- filter(simb, type == "real")
simg <- filter(simf, type == "grid")
simf <- filter(simf, type %in% c("real", "fake"))
nrow(simf)
nrow(simb)
nrow(simg)

ggplot(simf, aes(X, Y, colour = log(observed + 1))) +
  geom_point() +
  facet_wrap(~year) +
  scale_colour_viridis_c()

ggplot(simg, aes(X, Y, fill = log(observed + 1))) +
  geom_tile(width = 2, height = 2) +
  facet_wrap(~year) +
  scale_fill_viridis_c()

true_index <- simg |> group_by(year) |>
  summarise(est = sum(mu)) |>
  mutate(type = "True")

true_index |>
  ggplot(aes(year, est)) + geom_line() +
  labs(y = "True population index")

ggplot(simb, aes(X, Y, colour = log(observed + 1))) +
  geom_point() +
  facet_wrap(~year) +
  scale_colour_viridis_c()

mean(simf$observed)
mean(simf$observed == 0)
mean(d$observed)
mean(d$observed == 0)

# cross-fit models and simulations and test -----------------------------

sigma_V <- seq(0, 1, length.out = 500L)
p <- gamma_cv(0.2, 0.5)
plot(sigma_V, dgamma(sigma_V, shape = p[1], scale = p[2]), type = "l")
sdmTMBpriors(sigma_V = gamma_cv(0.2, 0.5))

meshs <- make_mesh(simf, c("X", "Y"), mesh = mesh$mesh)
mf <- sdmTMB(
  observed ~ 1,
  data = simf,
  time = "year",
  spatial = "on",
  time_varying = ~1,
  time_varying_type = "rw0",
  priors = sdmTMBpriors(sigma_V = gamma_cv(0.2, 0.5)),
  spatiotemporal = "ar1",
  mesh = meshs,
  silent = FALSE,
  family = family
)
meshb <- make_mesh(simb, c("X", "Y"), mesh = mesh$mesh)
mb <- sdmTMB(
  observed ~ 1,
  data = simb,
  time = "year",
  spatial = "on",
  time_varying = ~1,
  time_varying_type = "rw0",
  priors = sdmTMBpriors(sigma_V = gamma_cv(0.2, 0.5)),
  # control = sdmTMBcontrol(start = list(ln_tau_V = get_pars(mf)$ln_tau_V), map = list(ln_tau_V = factor(NA))),
  spatiotemporal = "ar1",
  mesh = meshb,
  silent = FALSE,
  family = family
)

mfii <- update(mf, spatiotemporal = "iid", formula. = observed ~ 0 + factor(year), time_varying = NULL, priors = sdmTMBpriors())
mbii <- update(mb, spatiotemporal = "iid", formula. = observed ~ 0 + factor(year), time_varying = NULL, control = sdmTMBcontrol(), priors = sdmTMBpriors())

pf <- predict(mf, newdata = nd, return_tmb_object = TRUE)
pb <- predict(mb, newdata = nd, return_tmb_object = TRUE)
pfii <- predict(mfii, newdata = nd, return_tmb_object = TRUE)
pbii <- predict(mbii, newdata = nd, return_tmb_object = TRUE)

indf <- get_index(pf, bias_correct = TRUE) |> left_join(dy, by = join_by(year))
indb <- get_index(pb, bias_correct = TRUE) |> left_join(dy, by = join_by(year))
indfii <- get_index(pfii, bias_correct = TRUE) |> left_join(dy, by = join_by(year))
indbii <- get_index(pbii, bias_correct = TRUE) |> left_join(dy, by = join_by(year))

indexes <- bind_rows(
  indf |> mutate(sampling = "full", model = "RW"),
  indb |> mutate(sampling = "biennial", model = "RW"),
  indfii |> mutate(sampling = "full", model = "IID"),
  indbii |> mutate(sampling = "biennial", model = "IID")
)

ggplot(indexes, aes(year, est, ymin = lwr, ymax = upr, colour = biennial_region)) +
  geom_pointrange() +
  facet_grid(sampling~model)

chi_square <- function(model_full, model_biennial) {
  ## as in Rufener, M.-C., Kristensen, K., Nielsen, J.R., and Bastardie, F.
  ## 2021. Bridging the gap between commercial fisheries and survey data to
  ## model the spatiotemporal dynamics of marine species. Ecological
  ## Applications 31(8): e02453. doi:10.1002/eap.2453.
  ## https://github.com/mcruf/LGNB/blob/master/R/Validation_and_Residuals.R
  par1 <- model_biennial$model$par
  par2 <- model_full$model$par
  obj1 <- model_biennial$tmb_obj
  obj2 <- model_full$tmb_obj

  ## Evaluate biennial likelihood from biennial model parameters
  obj1$env$beSilent()
  f1 <- as.numeric(obj1$fn(par1))
  ## Evaluate biennial likelihood from full model parameters
  f2 <- as.numeric(obj1$fn(par2))

  df <- attr(logLik(model_full), "df")
  fixed <- 1 - pchisq(2 * (f2 - f1), df = df)
  cat("Fixed theta p-value =", fixed, "\n")

  ## Similar test now including random effects
  f1.all <- obj1$env$f(obj1$env$last.par.best)
  f2.all <- obj1$env$f(obj2$env$last.par.best)
  df.all <- df + length(obj1$env$random)
  random <- 1 - pchisq(2 * (f2.all - f1.all), df = df.all)

  cat("Random theta p-value =", random, "\n")
}

chi_square(mfii, mbii)
chi_square(mf, mb)

# missing coverage?
mean(indbii$lwr > indfii$est | indbii$upr < indfii$est)
mean(indb$lwr > indfii$est | indb$upr < indfii$est)
mean(indb$lwr > indf$est | indb$upr < indf$est)
missing_covr <- data.frame(year = indbii$year[(indbii$lwr > indfii$est | indbii$upr < indfii$est)], model = "Model: IID")

ind <- left_join(ind, dy)
indii <- left_join(indii, dy)

ind_orig <- bind_rows(
  mutate(ind, model = "Model: RW"),
  mutate(indii, model = "Model: IID")
) |>
  mutate(sampling = "Sampling: biennial")

real <- ind_orig |>
  mutate(geo = exp(mean(log(est)))) |>
  mutate(lwr = lwr / geo) |>
  mutate(upr = upr / geo) |>
  mutate(est = est / geo) |>
  ggplot(aes(year, est, ymin = lwr, ymax = upr, colour = factor(biennial_region))) +
  geom_ribbon(fill = "grey90", colour = NA) +
  geom_pointrange() +
  theme(legend.position = "inside", legend.position.inside = c(0.7, 0.7)) + # axis.text.y = element_blank(), axis.ticks.y = element_blank(),
  scale_colour_brewer(palette = "Set2") +
  facet_grid(sampling~model) +
  coord_cartesian(ylim = c(0, 2.6)) +
  tagger::tag_facets(tag = "panel", tag_prefix = "(", position = "tl", tag_pool = letters) +
  labs(colour = "Biennial sampling region", shape = "Sampling type", x = "Year", y = "Relative biomass") +
  ggtitle("Real data")
real

true_index_scaled <- true_index |>
  mutate(geo = exp(mean(log(est)))) |>
  mutate(est = est / geo)

simulated <- indexes |>
  mutate(geo = exp(mean(log(est)))) |>
  mutate(lwr = lwr / geo) |>
  mutate(upr = upr / geo) |>
  mutate(est = est / geo) |>
  mutate(sampling = paste0("Sampling: ", sampling)) |>
  mutate(model = paste0("Model: ", model)) |>
  ggplot(aes(year, est, ymin = lwr, ymax = upr, colour = biennial_region)) +
  geom_ribbon(fill = "grey90", colour = NA) +
  geom_pointrange(position = position_dodge(width = 0.5)) +
  scale_shape_manual(values = c(21, 19)) +
  scale_colour_brewer(palette = "Set2") +
  facet_grid(sampling~model) +
  labs(colour = "Biennial sampling region", shape = "Sampling type", x = "Year", y = "Relative biomass") +
  theme(legend.position = "top") + # , axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  # coord_cartesian(ylim = c(2e5, 1e6)) +
  coord_cartesian(ylim = c(0, 2.6)) +
  ggtitle("Simulated data") +
  tagger::tag_facets(tag = "panel", tag_prefix = "(", position = "tl", tag_pool = letters[-c(1:2)]) +
  geom_line(data = true_index_scaled, mapping = aes(year, est), inherit.aes = FALSE)
  # geom_vline(data = missing_covr, mapping = aes(xintercept = year), lwd = 1.5, alpha = 0.3, colour = "blue")
simulated

library(patchwork)
design <- "
  1
  2
  2
"

real +
  simulated + theme(legend.position = "none") +
  plot_layout(design = design)

# ggsave("figs/self-cross-test-example-lingcod.pdf", width = 7, height = 7)

indexes |>
  mutate(sampling = paste0("Sampling: ", sampling)) |>
  mutate(model = paste0("Model: ", model)) |>
  mutate(type = paste(sampling, model, sep = "; ")) |>
  ggplot(aes(year, est, ymin = lwr, ymax = upr, colour = type, group = type)) +
  geom_ribbon(aes(fill = type), colour = NA, alpha = 0.3) +
  # geom_pointrange(position = position_dodge(width = 0.2), alpha = 0.8) +
  geom_line(alpha = 0.8) +
  scale_shape_manual(values = c(21, 19)) +
  scale_colour_brewer(palette = "Set2") +
  scale_fill_brewer(palette = "Set2") +
  labs(fill = "Type", colour = "Type", x = "Year", y = "Index") +
  # facet_grid(sampling~model) +
  # labs(colour = "Biennial sampling region", shape = "Sampling type", x = "Year", y = "Relative biomass") +
  theme(legend.position = "top", axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  coord_cartesian(ylim = c(2e5, 1e6))
  # geom_vline(xintercept = missing_covr$year, lwd = 1.5, alpha = 0.3, lty = 1)

# crank down on sequence of rw intercept SDs!? -----------------------------

sds <- c(1e-3, 1e-2, seq(0.1, 1, 0.2), 1e1, 1e2, 1e3)
sds
length(sds)

get_ind <- function(s) {
  mb2 <- sdmTMB(
    observed ~ 1,
    data = simb,
    time = "year",
    spatial = "on",
    time_varying = ~ 1,
    time_varying_type = "rw0",
    control = sdmTMBcontrol(start = list(ln_tau_V = matrix(log(s), nrow = 1, ncol = 1)), map = list(ln_tau_V = factor(NA))),
    spatiotemporal = "iid",
    mesh = meshb,
    silent = FALSE,
    family = family
  )
  pb2 <- predict(mb2, newdata = nd, return_tmb_object = TRUE)
  indb2 <- get_index(pb2, bias_correct = TRUE) |> left_join(dy, by = join_by(year))
  indb2$rw_sd <- s
  indb2
}

library(future)
future::plan(future::multisession)
out <- furrr::future_map_dfr(sds, get_ind)

out |>
  mutate(geo = exp(mean(log(est)))) |>
  mutate(lwr = lwr / geo) |>
  mutate(upr = upr / geo) |>
  mutate(est = est / geo) |>
  mutate(rw_sd = paste0("RW SD: ", rw_sd)) |>
  ggplot(aes(year, est, ymin = lwr, ymax = upr, colour = biennial_region)) +
  geom_ribbon(fill = "grey90", colour = NA) +
  geom_pointrange(position = position_dodge(width = 0.5)) +
  facet_wrap(~rw_sd) +
  scale_colour_brewer(palette = "Set2") +
  geom_line(data = true_index_scaled, mapping = aes(year, est), inherit.aes = FALSE) +
  coord_cartesian(ylim = c(0, 2.6))

# with a prior!?


sigma_V <- seq(0, 1, length.out = 1000)
p <- gamma_cv(0.2, 0.4)
plot(sigma_V, dgamma(x, shape = p[1], scale = p[2]), type = "l")
sdmTMBpriors(sigma_V = gamma_cv(0.2, 0.2))
mb2 <- sdmTMB(
  observed ~ 1,
  data = simb,
  time = "year",
  spatial = "on",
  time_varying = ~ 1,
  time_varying_type = "rw0",
  priors = sdmTMBpriors(sigma_V = gamma_cv(0.2, 0.2)),
  spatiotemporal = "iid",
  mesh = meshb,
  silent = FALSE,
  family = family
)
mb2
pb2 <- predict(mb2, newdata = nd, return_tmb_object = TRUE)
indb2 <- get_index(pb2, bias_correct = TRUE) |> left_join(dy, by = join_by(year))
indb2 |>
  mutate(geo = exp(mean(log(est)))) |>
  mutate(lwr = lwr / geo) |>
  mutate(upr = upr / geo) |>
  mutate(est = est / geo) |>
  ggplot(aes(year, est, ymin = lwr, ymax = upr, colour = biennial_region)) +
  geom_ribbon(fill = "grey90", colour = NA) +
  geom_pointrange(position = position_dodge(width = 0.5)) +
  scale_colour_brewer(palette = "Set2") +
  geom_line(data = true_index_scaled, mapping = aes(year, est), inherit.aes = FALSE) +
  coord_cartesian(ylim = c(0, 2.6))
