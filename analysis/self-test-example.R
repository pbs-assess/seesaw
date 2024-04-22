# - simulate with new random fields?? and with full vs. the actual survey coverage... do they substantially differ? (maybe chi-squared test idea from Rufener et al.)
# - do with underlying IID and with underlying RW... must be self consistent to pass

library(dplyr)
library(ggplot2)

d <- readRDS("~/src/gfsynopsis-2021/report/data-cache-nov-2023/quillback-rockfish.rds")$survey_sets
d <- filter(d, survey_abbrev %in% c("HBLL OUT N", "HBLL OUT S"))

glimpse(d)

library(sdmTMB)

d <- add_utm_columns(d)
mesh <- make_mesh(d, c("X", 'Y'), cutoff = 15)
plot(mesh)

fit <- sdmTMB(
  # catch_count ~ 0 + factor(year),
  catch_count ~ 1,
  offset = log(d$hook_count),
  data = d,
  time = "year",
  spatial = "on",
  spatiotemporal = 'rw',
  mesh = mesh,
  silent = FALSE,
  family = nbinom2()
)

fit

dy <- group_by(d, year) |>
  summarise(north = grepl("N$", survey_abbrev)[1])

grid <- gfplot::hbll_grid$grid
grid <- rename(grid, lon = X, lat = Y) |>
  add_utm_columns(ll_names = c("lon", "lat"))

nd <- replicate_df(grid, "year", unique(d$year))

p <- predict(fit, newdata = nd, return_tmb_object = TRUE)
ind <- get_index(p, bias_correct = TRUE)

ind <- left_join(ind, dy)

ggplot(ind, aes(year, est, ymin = lwr, ymax = upr, colour = north)) +
  geom_pointrange()

b1 <- tidy(fit)
b2 <- tidy(fit, 'ran_pars')

omega <- get_pars(fit)$omega_s
eps <- get_pars(fit)$epsilon_st

# want to create fake sampling in the unsampled region in a given year
# want to see if IID and RW give plausibly similar indexes in partial vs. full sampling
# to create these fake sampling locations, quick hack is to take sampling locations randomly from another sampled year

group_by(d, year) |>
  summarise(n = n())

get_random_year_samples <-
  function(data, .survey_abbrev = c("HBLL OUT N", "HBLL OUT S")) {
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
    other_region <- if (region == "HBLL OUT N") "HBLL OUT S" else "HBLL OUT N"
    other_loc <- get_random_year_samples(d, other_region) |>
      mutate(type = "fake") |>
      mutate(year = unique(.x$year))
    bind_rows(real_loc, other_loc)
  }) |>
  bind_rows()

meshs <- make_mesh(locs, c("X", "Y"), mesh = mesh$mesh)
plot(meshs)
meshs$mesh$n

ggplot(d, aes(X, Y)) + geom_point() +
  facet_wrap(~year)

ggplot(locs, aes(X, Y)) + geom_point() +
  facet_wrap(~year)

simf <- sdmTMB_simulate(
  ~ 1,
  data = locs,
  sigma_O = b2$estimate[b2$term == "sigma_O"],
  sigma_E = b2$estimate[b2$term == "sigma_E"],
  phi = b2$estimate[b2$term == "phi"],
  range = b2$estimate[b2$term == "range"],
  fixed_re = list(omega_s = omega, epsilon_st = eps),
  time = "year",
  seed = 12345,
  family = nbinom2(),
  mesh = meshs,
  offset = rep(log(450), nrow(locs)),
  B = unname(coef(fit)) # intercept
)

simb <- simf
simb$type <- locs$type
simb <- filter(simb, type == "real")
nrow(simb)
nrow(simf)

ggplot(simf, aes(X, Y, colour = log(observed))) + geom_point() +
  facet_wrap(~year) + scale_colour_viridis_c()

mean(simf$observed)
mean(simf$observed == 0)
mean(d$catch_count)
mean(d$catch_count == 0)
head(simf)

mf <- sdmTMB(
  observed ~ 1,
  data = simf,
  time = "year",
  spatial = "on",
  spatiotemporal = 'rw',
  mesh = meshs,
  silent = FALSE,
  family = nbinom2()
)
meshb <- make_mesh(simb, c("X", "Y"), mesh = mesh$mesh)
mb <- sdmTMB(
  observed ~ 1,
  data = simb,
  time = "year",
  spatial = "on",
  spatiotemporal = 'rw',
  mesh = meshb,
  silent = FALSE,
  family = nbinom2()
)

mfii <- update(mf, spatiotemporal = "iid", formula. = observed ~ 0 + factor(year))
mbii <- update(mb, spatiotemporal = "iid", formula. = observed ~ 0 + factor(year))

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
