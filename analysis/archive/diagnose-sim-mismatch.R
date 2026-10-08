library(ggplot2)
library(dplyr)
library(sdmTMB)
theme_set(ggsidekick::theme_sleek() + theme(tagger.panel.tag.text = element_text(color = "grey30", size = 9)))

d <- readRDS("~/src/gfsynopsis-2024/report/data-cache-2025-03/quillback-rockfish.rds")$survey_sets
d <- filter(d, survey_abbrev %in% c("HBLL INS N", "HBLL INS S"))
# simplify, remove N in 2021 (was N + S)
# d <- filter(d, !(year == 2021 & survey_abbrev %in% "HBLL INS N"))
d <- mutate(d, observed = catch_count, log_effort = log(hook_count))
d$biennial_region <- as.integer(as.factor(d$survey_abbrev))
table(d$biennial_region, d$year)
d$biennial_region[d$year == 2021] <- 3L
family <- sdmTMB::nbinom2()

grid1 <- gfplot::hbll_inside_n_grid$grid
grid2 <- gfplot::hbll_inside_s_grid$grid
grid <- bind_rows(grid1, grid2)
grid <- rename(grid, lon = X, lat = Y) |>
  add_utm_columns(ll_names = c("lon", "lat"))

# fit initial models --------------------------------------------

d <- add_utm_columns(d)
mesh <- make_mesh(d, c("X", "Y"), cutoff = 10)
# plot(mesh)

fit <- sdmTMB(
  observed ~ factor(year),
  offset = "log_effort",
  data = d,
  time = "year",
  spatial = "on",
  #time_varying = ~1,
  #time_varying_type = "rw0",
  # priors = sdmTMBpriors(sigma_V = gamma_cv(0.2, 0.4)),
  control = sdmTMBcontrol(multiphase = FALSE),
  spatiotemporal = "iid",
  #  extra_time = seq(min(d$year), max(d$year)),
  mesh = mesh,
  silent = FALSE,
  family = family
)
fit
sanity(fit)

dy <- group_by(d, year) |>
  summarise(biennial_region = factor(biennial_region[1]))

nd <- replicate_df(grid, "year", unique(d$year))

p <- predict(fit, newdata = nd, return_tmb_object = TRUE)
ind <- get_index(p, bias_correct = TRUE)

left_join(ind, dy) |> 
  ggplot(aes(year, est, ymin = lwr, ymax = upr, colour = factor(biennial_region))) +
  geom_pointrange()

