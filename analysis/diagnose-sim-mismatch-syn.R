library(ggplot2)
library(dplyr)
library(sdmTMB)
theme_set(ggsidekick::theme_sleek() + theme(tagger.panel.tag.text = element_text(color = "grey30", size = 9)))

d <- readRDS("~/src/gfsynopsis-2024/report/data-cache-2025-03/lingcod.rds")$survey_sets
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
# d <- filter(d, !(year %in% c(2021) & survey_abbrev == "SYN WCVI"))
d <- filter(d, !(year %in% c(2007) & survey_abbrev == "SYN WCHG"))
d <- filter(d, survey_abbrev != "SYN WCHG")
d$biennial_region <- ifelse(d$survey_abbrev %in% c("SYN WCVI", "SYN WCHG"), 1L, 2L)
d$biennial_region[d$year == 2021] <- 3L

d <- filter(d, survey_abbrev %in% c("SYN QCS", "SYN HS"))
table(d$year, d$survey_abbrev)
table(d$year, d$biennial_region)

d <- mutate(d, observed = density_kgkm2, log_effort = log(area_swept))
family <- tweedie()

grid <- gfplot::synoptic_grid |> select(-survey_domain_year, -utm_zone) |>
  filter(survey != "SYN WCHG")  |> filter(survey != "SYN WCVI")
# max(grid$Y[grid$survey == "SYN WCVI"]) - min(grid$Y[grid$survey == "SYN QCS"])
ggplot(grid, aes(X, Y, colour = survey)) + geom_point()
table(grid$survey)

# fit initial models --------------------------------------------

d <- add_utm_columns(d)
mesh <- make_mesh(d, c("X", "Y"), cutoff = 6)
plot(mesh)

fit <- sdmTMB(
  observed ~ factor(year),
  offset = "log_effort",
  data = d,
  time = "year",
  spatial = "on",
  #time_varying = ~1,
  #time_varying_type = "rw0",
  # priors = sdmTMBpriors(sigma_V = gamma_cv(0.2, 0.4)),
  # control = sdmTMBcontrol(multiphase = FALSE),
  spatiotemporal = "iid",
  #  extra_time = seq(min(d$year), max(d$year)),
  mesh = mesh,
  silent = FALSE,
  family = family
)
fit
b <- tidy(fit, "ran_pars")
sanity(fit)

.range <- b$estimate[b$term == "range"]
diff(range(grid$Y)) / .range

dy <- group_by(d, year) |>
  summarise(biennial_region = factor(biennial_region[1]))

nd <- replicate_df(grid, "year", unique(d$year))

p <- predict(fit, newdata = nd, return_tmb_object = TRUE)
ind <- get_index(p, bias_correct = TRUE)

left_join(ind, dy) |> 
  ggplot(aes(year, est, ymin = lwr, ymax = upr, colour = factor(biennial_region))) +
  geom_pointrange()

