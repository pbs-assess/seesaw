# - simulate with new random fields with full vs. the actual survey coverage
# - do they substantially differ? (maybe chi-squared test idea from Rufener et al.)
# - do with underlying IID and with underlying RW... must be self consistent

library(dplyr)
library(ggplot2)
theme_set(ggsidekick::theme_sleek())
library(sdmTMB)

# bring in set and grid data ------------------------------------

CASE <- "quillback inside"
SURVEY <- "inside"

if (CASE == "quillback inside") {
  d <- readRDS("~/src/gfsynopsis-2022/report/data-cache-nov-2023/quillback-rockfish.rds")$survey_sets
  d <- filter(d, survey_abbrev %in% c("HBLL INS N", "HBLL INS S"))
  # simplify, remove N in 2021 (was N + S)
  d <- filter(d, !(year == 2021 & survey_abbrev %in% "HBLL INS N"))
  d <- mutate(d, response = catch_count, log_effort = log(hook_count))
  family <- nbinom2()
}

if (SURVEY == "inside") {
  grid1 <- gfplot::hbll_inside_n_grid$grid
  grid2 <- gfplot::hbll_inside_s_grid$grid
  grid <- bind_rows(grid1, grid2)
  grid <- rename(grid, lon = X, lat = Y) |>
    add_utm_columns(ll_names = c("lon", "lat"))
}

# fit initial models --------------------------------------------

d <- add_utm_columns(d)
mesh <- make_mesh(d, c("X", "Y"), cutoff = 10)
plot(mesh)

fit <- sdmTMB(
  response ~ 1,
  offset = "log_effort",
  data = d,
  time = "year",
  spatial = "on",
  spatiotemporal = "rw",
  mesh = mesh,
  silent = FALSE,
  family = family
)

fit

fitii <- update(
  fit,
  spatiotemporal = "iid",
  formula. = response ~ 0 + factor(year)
)
fitii2 <- update(fit,
  control = sdmTMBcontrol(profile = c("ln_phi")),
  spatiotemporal = "iid",
  formula. = response ~ 0 + factor(year)
)

dy <- group_by(d, year) |>
  summarise(north = grepl("N$", survey_abbrev)[1])

nd <- replicate_df(grid, "year", unique(d$year))

p <- predict(fit, newdata = nd, return_tmb_object = TRUE)
pii <- predict(fitii, newdata = nd, return_tmb_object = TRUE)
ind <- get_index(p, bias_correct = TRUE)

pii2 <- predict(fitii2, newdata = nd, return_tmb_object = TRUE)
ind2 <- get_index(pii2, bias_correct = TRUE)

indii <- get_index(pii, bias_correct = TRUE)

ind <- left_join(ind, dy)
indii <- left_join(indii, dy)

ggplot(ind, aes(year, est, ymin = lwr, ymax = upr, colour = north)) +
  geom_pointrange()

ggplot(indii, aes(year, est, ymin = lwr, ymax = upr, colour = north)) +
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
  function(data, .survey_abbrev = c("HBLL INS N", "HBLL INS S")) {
    .survey_abbrev <- match.arg(.survey_abbrev)
    dd <- filter(data, survey_abbrev %in% .survey_abbrev)
    yrs <- unique(dd$year)
    yr <- sample(yrs, 1)
    filter(dd, year %in% yr) |>
      select(X, Y)
  }

set.seed(1)
locs <- group_by(d, year) |>
  group_split() |>
  purrr::map(\(.x) {
    real_loc <- select(.x, X, Y, year) |> mutate(type = "real")
    region <- unique(.x$survey_abbrev)
    other_region <- if (region == "HBLL INS N") "HBLL INS S" else "HBLL INS N"
    other_loc <- get_random_year_samples(d, other_region) |>
      mutate(type = "fake") |>
      mutate(year = unique(.x$year))
    bind_rows(real_loc, other_loc)
  }) |>
  bind_rows()

meshs <- make_mesh(locs, c("X", "Y"), mesh = mesh$mesh)
plot(meshs)
meshs$mesh$n

ggplot(d, aes(X, Y)) +
  geom_point() +
  facet_wrap(~year)

ggplot(locs, aes(X, Y)) +
  geom_point() +
  facet_wrap(~year)

simf <- sdmTMB_simulate(
  ~1,
  data = locs,
  sigma_O = b2$estimate[b2$term == "sigma_O"],
  sigma_E = b2$estimate[b2$term == "sigma_E"],
  phi = b2$estimate[b2$term == "phi"],
  range = b2$estimate[b2$term == "range"],
  fixed_re = list(omega_s = omega, epsilon_st = eps),
  time = "year",
  seed = 12345,
  family = family,
  mesh = meshs,
  offset = rep(mean(d$log_effort), nrow(locs)),
  B = unname(coef(fit)) # intercept
)

simb <- simf
simb$type <- locs$type
simb <- filter(simb, type == "real")
nrow(simb)
nrow(simf)

ggplot(simf, aes(X, Y, colour = log(observed))) +
  geom_point() +
  facet_wrap(~year) +
  scale_colour_viridis_c()

mean(simf$observed)
mean(simf$observed == 0)
mean(d$response)
mean(d$response == 0)
head(simf)

# cross-fit models and simulations and test -----------------------------

mf <- sdmTMB(
  observed ~ 1,
  data = simf,
  time = "year",
  spatial = "on",
  time_varying = ~1,
  time_varying_type = "ar1",
  spatiotemporal = "rw",
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
  time_varying_type = "ar1",
  spatiotemporal = "rw",
  mesh = meshb,
  silent = FALSE,
  family = family
)

mfii <- update(mf, spatiotemporal = "iid", formula. = observed ~ 0 + factor(year), time_varying = NULL)
mbii <- update(mb, spatiotemporal = "iid", formula. = observed ~ 0 + factor(year), time_varying = NULL)

pf <- predict(mf, newdata = nd, return_tmb_object = TRUE)
pb <- predict(mb, newdata = nd, return_tmb_object = TRUE)
pfii <- predict(mfii, newdata = nd, return_tmb_object = TRUE)
pbii <- predict(mbii, newdata = nd, return_tmb_object = TRUE)

indf <- get_index(pf, bias_correct = TRUE) |> left_join(dy)
indb <- get_index(pb, bias_correct = TRUE) |> left_join(dy)
indfii <- get_index(pfii, bias_correct = TRUE) |> left_join(dy)
indbii <- get_index(pbii, bias_correct = TRUE) |> left_join(dy)

ii <- bind_rows(
  indf |> mutate(type = "RW full"),
  indb |> mutate(type = "RW biennial"),
  indfii |> mutate(type = "IID full"),
  indbii |> mutate(type = "IID biennial")
)

ggplot(ii, aes(year, est, ymin = lwr, ymax = upr, colour = north)) +
  geom_pointrange() +
  facet_wrap(~type)

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

# conclusions:
# - random walk random field is consistent with biennial or full sampling
# - IID does not (when including random effects)

ii |>
  mutate(biennial = grepl("biennial", type)) |>
  mutate(iid = paste("Random walk random field: ", !grepl("IID", type))) |>
  ggplot(aes(year, est, ymin = lwr, ymax = upr, shape = north, colour = biennial)) +
  geom_pointrange(position = position_dodge(width = 0.5)) +
  scale_shape_manual(values = c("FALSE" = 21, "TRUE" = 19)) +
  scale_colour_brewer(palette = "Set2") +
  facet_wrap(~iid) +
  labs(colour = "Biennial sampling", shape = "Northern region sampled\n(if biennial)", x = "Year", y = "Relative abundance") +
  theme(legend.position = "bottom")

ggsave("figs/self-cross-test-example.pdf", width = 7, height = 4)

# RW may be over smoothed
# consider relaxing the smoothing here?
