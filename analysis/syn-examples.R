# General workflow for testing with data

# 1) Wrangle the survey data for trawl to start
# 2a) Fit test models to single species
# 2b) Convert to function using try() to fit all models
# 3a) Get index from each fitted model;
# 3b) Convert to function using try() to apply to multiple species
# 4)

# 00) Describe the differences in the species, particularly with respect to
#     how characteristics of the different species compare or relate to the
#     simulated data examples.
# E.g., + range of autocorrelation
#       + gap between surveys
#       + expectations in variability from year to year

# Questions
# - "Survey domain year" --> this has to do with when the grid was adjusted?
#   I.e., the last time it was adjusted to account for dropped cells/areas?
# - Can we walk through the different model scenarios to make sure I understand?
#   it might help me if I could make a table describing what each is doing
# - Why is offset area divided by 100 000 instead of 1 000 000? m2 to km2 should
#   be 1e-6?
# - Should priors be used? Should they be standardised like they were in the
# - sdmTMBcontrol(nlminb_loops = 1L, newton_loops = 1L); default? nlminb_loops?
# - index_args; an argument in sdmTMB; these should be a default in fitting the
#   model fitting?
# - How to decide what k should be when using s(year)
# - I need to remind myself what this bias_correct is about for the index fitting
# - If an offset is included in the original model; does this follow through in
#   when we use predict.sdmTMB

# Should 2021 WCVI be left out? Should the analysis be compared with and without
# this year? To see if it has any affect?
# distinct(dat, year, survey_abbrev) %>% arrange(survey_abbrev, year)
# dat %>% filter(year == 2021, survey_abbrev == "SYN WCVI")


library(tidyverse)
library(gfplot)
library(sdmTMB)
library(beepr)
library(patchwork)

source(here::here("analysis", "00_prep-example-data.R"))
source(here::here("analysis", "utils.R"))
source(here::here("analysis", "fit-models.R"))
source(here::here("analysis", "fit-funcs.R"))

mytheme <- function() ggsidekick::theme_sleek() # sometimes I add more layers to themes
theme_set(mytheme())

plot_index <- function(df) {
  df |>
    ggplot(mapping = aes(x = year, y = est, ymin = lwr, ymax = upr, colour = region)) +
    geom_pointrange() +
    geom_ribbon(alpha = 0.20, colour = NA) +
    scale_colour_manual(values = syn_region_colours$colours, breaks = syn_region_colours$region) +
    labs(colour = "Sampled region")
}


# Setup inside data and look at survey coverage over time
# ------------------------------------------------------------------------------
syn_survey_dat <- dat |>
  filter(str_detect(survey_abbrev, "SYN")) |>
  filter(!(year %in% c(2003, 2004, 2020))) |> # Use only complete N/S sampling years
  drop_na(trawl_offset) # There are some NA values in the offset

# Consider excluding 2007; or exclude the WCHG from 2007 to make the data more comparable?

# Come back to this because I don't know that the categories are straightforward
syn_region_colours <- tibble(
  region = c("QCS + HS", "WCHG + WCVI", "QCS + HS + WCVI", "QCS + HS + WCHG", "No data"),
  colours = c(RColorBrewer::brewer.pal(4L, "Set2"), "grey70")
)

# ggplot(data = syn_survey_dat) +
#   geom_point(aes(x = X, y = Y, colour = survey_abbrev)) +
#   scale_colour_manual(
#     values = RColorBrewer::brewer.pal(4L, "Set2"),
#     breaks = c("SYN QCS", "SYN HS", "SYN WCVI", "SYN WCHG"), na.translate = FALSE
#   ) +
#   facet_wrap(~fyear, drop = FALSE)

# dir.create("figs")
# ggsave(here::here("figs", "synoptic-survey-coverage.png"), width = 7.5, height = 6.5)

syn_survey_yrs <-
  syn_survey_dat |>
  filter(str_detect(survey_abbrev, "SYN")) |>
  distinct(year, survey_abbrev) |>
  arrange(year) |>
  rename(region = "survey_abbrev") |>
  mutate(surveyed = 1) |>
  pivot_wider(id_cols = year, names_from = region, values_from = surveyed) |>
  mutate(region = ifelse(year %% 2 == 0, "QCS + HS", "WCHG + WCVI")) |>
  mutate(
    region = replace(region, year == 2007, "QCS + HS + WCHG"),
    region = replace(region, year == 2021, "QCS + HS + WCVI")
  ) |>
  select(year, region) |>
  bind_rows(tibble(year = 2020, region = "No data"))

fitted_yrs_extra <- min(syn_survey_dat$year):max(syn_survey_dat$year)

syn_nd <-
  syn_grid |>
  rename(area = "cell_area") |>
  make_grid(years = fitted_yrs_extra) |>
  mutate(fyear = as.factor(year))


# Arrowtooth & Bocaccio examples
# ------------------------------------------------------------------------------
get_syn_index <- function(sp_dat) {
  future::plan(future::multisession, workers = 5) # or whatever number

  sp_dat <-
    sp_dat |>
    group_by(species_common_name) |>
    group_split()

  map_func <- furrr::future_map
  # map_func <- purrr::map

  fits1 <- sp_dat |>
    map_func(
      fit_models,
      catch = "catch_weight", silent = FALSE,
      family = tweedie(), offset = "trawl_offset", data_subset = "data_subset"
    ) |>
    list_flatten(name_spec = "{inner}")
  future::plan(future::sequential)

  fits_cleaned1 <- fits1 |>
    map(check_sanity) # omit plots made from models that did not pass sanity check

  message("Getting predictions")
  preds1 <- get_pred_list(fits_cleaned1, newdata = syn_nd)
  message("Getting index")
  indices1 <- get_index_list(pred_list = preds1)

  indices1
}

# ------------------------------------------------------------------------------
# Synoptic trawl example species
spp <- c("arrowtooth flounder", "bocaccio")

# Strict alternating
# -----------------------
sp_dat <- syn_survey_dat |>
  filter(species_common_name %in% spp) |>
  filter(
    !(year == 2007 & region == "SYN WCHG"),
    !(year == 2020),
    !(year == 2021 & region == "SYN WCVI")
  ) |>
  mutate(data_subset = paste(species_common_name, "Strict alternating N/S", sep = "-"))

# alternating_ind_list <- get_syn_index(sp_dat)
# saveRDS(alternating_ind_list, here::here('data-outputs', "syn_strict-alternating_ind-list.RDS"))
alternating_ind_list <- readRDS(here::here("data-outputs", "syn_strict-alternating_ind-list.RDS"))

alternating_ind_df <- mk_index_df(alternating_ind_list) |>
  separate(group, into = c("species", "group"), sep = "-") |>
  left_join(syn_survey_yrs, by = "year") |>
  mutate(region = ifelse(year %% 2 == 0, "QCS + HS", "WCHG + WCVI")) |>
  mutate(region = replace(region, year == 2020, "No data")) |>
  left_join(syn_region_colours, by = "region")

# Including broader spatial coverage: 2007 WCHG + QCS/HS & 2021 WCVI + QCS/HS &
# ------------------------------------------------------------------------------
sp_dat <- syn_survey_dat |>
  filter(species_common_name %in% spp) |>
  mutate(data_subset = paste(species_common_name, "Additional region sampling", sep = "-"))

# both_regions_ind_list <- get_syn_index(sp_dat)
# beep()
# saveRDS(both_regions_ind_list, here::here('data-outputs', "syn_extended-sp-coverage_ind-list.RDS"))
both_regions_ind_list <- readRDS(here::here("data-outputs", "syn_extended-sp-coverage_ind-list.RDS"))

both_regions_df <- mk_index_df(both_regions_ind_list) |>
  separate(group, into = c("species", "group"), sep = "-") |>
  left_join(syn_survey_yrs, by = "year") |>
  left_join(syn_region_colours, by = "region")

combined_ind_df <- bind_rows(alternating_ind_df, both_regions_df)

p1 <-
  combined_ind_df |>
  filter(species == "arrowtooth flounder", group == "Strict alternating N/S") |>
  plot_index() +
  facet_wrap(~ fct_reorder(desc, id), nrow = 2L, scales = "free_y", as.table = TRUE) +
  ggtitle("Arrowtooth Flounder - Alternating N/S sampling")

p2 <-
  combined_ind_df |>
  filter(species == "arrowtooth flounder", group == "Additional region sampling") |>
  plot_index() +
  facet_wrap(~ fct_reorder(desc, id), nrow = 2L, scales = "free_y", as.table = TRUE) +
  ggtitle("Arrowtooth Flounder - Additional region coverage")

p3 <-
  combined_ind_df |>
  filter(species == "bocaccio", group == "Strict alternating N/S") |>
  plot_index() +
  facet_wrap(~ fct_reorder(desc, id), nrow = 2L, scales = "free_y", as.table = TRUE) +
  ggtitle("Bocaccio - Alternating N/S sampling")

p4 <-
  combined_ind_df |>
  filter(species == "bocaccio", group == "Additional region sampling") |>
  plot_index() +
  facet_wrap(~ fct_reorder(desc, id), nrow = 2L, scales = "free_y", as.table = TRUE) +
  ggtitle("Bocaccio - Additional region coverage")

p1 / p2

p3 / p4
