library(tidyverse)
library(gfplot)
library(sdmTMB)
library(beepr)
library(patchwork)

source(here::here("analysis", "00-prep-example-data.R"))
source(here::here("analysis", "utils.R"))
source(here::here("analysis", "fit-funcs.R"))
source(here::here("analysis", "fit-models.R"))

mytheme <- function() ggsidekick::theme_sleek() # sometimes I add more layers to themes
theme_set(mytheme())

# family <- nbinom2(link = 'log') # could use this instead, but would need to change units and add offset

# Setup inside data and look at survey coverage over time
# ------------------------------------------------------------------------------
inside_survey_dat <- dat |> filter(str_detect(survey_abbrev, "HBLL INS"))

inside_region_colours <- tibble(
  region = c("HBLL INS N", "HBLL INS S", "Both", "No data"),
  colours = c(RColorBrewer::brewer.pal(3L, "Set2"), "grey70")
)

ggplot(data = inside_survey_dat) +
  geom_point(aes(x = X, y = Y, colour = survey_abbrev)) +
  scale_colour_manual(
    values = inside_region_colours$colours,
    breaks = inside_region_colours$region, na.translate = FALSE
  ) +
  facet_wrap(~fyear, drop = FALSE)

inside_survey_yrs <-
  dat |>
  filter(str_detect(survey_abbrev, "HBLL INS")) |>
  distinct(year, survey_abbrev) |>
  arrange(year) |>
  rename(region = "survey_abbrev") |>
  complete(year = 2003:2022, fill = list(region = "None", surveyed = FALSE)) |>
  mutate(surveyed = 1) |>
  pivot_wider(id_cols = year, names_from = region, values_from = surveyed) |>
  mutate(region = case_when(
    `HBLL INS N` == 1 & is.na(`HBLL INS S`) & is.na(None) ~ "HBLL INS N",
    is.na(`HBLL INS N`) & `HBLL INS S` == 1 & is.na(None) ~ "HBLL INS S",
    `HBLL INS N` == 1 & `HBLL INS S` == 1 & is.na(None) ~ "Both",
    is.na(`HBLL INS N`) & is.na(`HBLL INS S`) & None == 1 ~ "No data"
  )) |>
  select(year, region)

fitted_yrs_extra <- min(inside_survey_dat$year):max(inside_survey_dat$year)

inside_nd <-
  hbll_inside |>
  make_grid(years = fitted_yrs_extra) |>
  mutate(fyear = as.factor(year))


# QB and YE inside all years
# ------------------------------------------------------------------------------
inside_dat <-
  inside_survey_dat |>
  filter(species_common_name %in% c("quillback rockfish")) |>
  mutate(data_subset = paste(species_common_name, "Stitched N/S", sep = "-")) |>
  group_by(species_common_name) |>
  group_split()

future::plan(future::multisession, workers = 5) # or whatever number
fits1 <- inside_dat |>
  furrr::future_map(fit_models,
    catch = "catch_count", family = nbinom2(),
    offset = "hook_offset",
    data_subset = "data_subset"
  ) |>
  list_flatten(name_spec = "{inner}")
future::plan(future::sequential)

fits_cleaned1 <- fits1 |>
  map(check_sanity) # omit plots made from models that did not pass sanity check

preds1 <- get_pred_list(fits_cleaned1, newdata = inside_nd)
indices1 <- get_index_list(pred_list = preds1)
# beep()

index_df1 <-
  mk_index_df(indices1) |>
  left_join(inside_survey_yrs) |>
  separate(group, into = c("species", "group"), sep = "-")

p1 <-
  ggplot(data = index_df1, aes(x = year, y = est, ymin = lwr, ymax = upr)) +
  geom_pointrange(aes(colour = region)) +
  geom_ribbon(alpha = 0.20, colour = NA) +
  scale_colour_manual(
    values = inside_region_colours$colours,
    breaks = inside_region_colours$region, na.translate = FALSE
  ) +
  labs(colour = "Sampled region") +
  facet_wrap(species ~ desc, nrow = 2L, scales = "free_y") +
  ggtitle("Stitched N/S - All Years")

# QB and YE inside all years
# ------------------------------------------------------------------------------
inside_dat_no2021 <-
  inside_survey_dat |>
  filter(year != 2021) |>
  filter(species_common_name %in% c("quillback rockfish", "yelloweye rockfish")) |>
  mutate(data_subset = paste(species_common_name, "Stitched N/S no 2021", sep = "-")) |>
  group_by(species_common_name) |>
  group_split()

future::plan(future::multisession, workers = 5) # or whatever number
fits2 <- inside_dat_no2021 |>
  furrr::future_map(fit_models,
    catch = "catch_count", family = nbinom2(),
    offset = "hook_offset",
    data_subset = "data_subset"
  ) |>
  list_flatten(name_spec = "{inner}")
future::plan(future::sequential)

fits_cleaned2 <- fits2 |>
  map(check_sanity) # omit plots made from models that did not pass sanity check

preds2 <- get_pred_list(fits_cleaned2, newdata = inside_nd)
indices2 <- get_index_list(pred_list = preds2)
# beep()

index_df2 <-
  mk_index_df(indices2) |>
  left_join(inside_survey_yrs) |>
  separate(group, into = c("species", "group"), sep = "-") |>
  mutate(region = ifelse(region %in% c("Both", "None"), "No data", region))

p2 <-
  ggplot(data = index_df2, aes(x = year, y = est, ymin = lwr, ymax = upr)) +
  geom_pointrange(aes(colour = region)) +
  geom_ribbon(alpha = 0.20, colour = NA) +
  scale_colour_manual(
    values = inside_region_colours$colours,
    breaks = inside_region_colours$region, na.translate = FALSE
  ) +
  labs(colour = "Sampled region") +
  facet_wrap(species ~ desc, nrow = 2L, scales = "free_y") +
  ggtitle("Stitched N/S - No 2021")

p1 / p2
