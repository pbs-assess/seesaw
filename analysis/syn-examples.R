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

source(here::here("analysis", "00_prep-example-data.R"))
source(here::here("analysis", "utils.R"))
source(here::here("analysis", "fit-models.R"))
source(here::here("analysis", "fit-funcs.R"))

mytheme <- function() ggsidekick::theme_sleek() # sometimes I add more layers to themes
theme_set(mytheme())

# Setup inside data and look at survey coverage over time
# ------------------------------------------------------------------------------
syn_survey_dat <- dat |>
  filter(str_detect(survey_abbrev, "SYN")) |>
  drop_na(trawl_offset) # There are some NA values in the offset

# Come back to this because I don't know that the categories are straightforward
syn_region_colours <- tibble(
  region = c("SYN QCS", "SYN HS", "SYN WCVI", "SYN WCHG"),
  colours = c(RColorBrewer::brewer.pal(4L, "Set2"))
)

ggplot(data = syn_survey_dat) +
  geom_point(aes(x = X, y = Y, colour = survey_abbrev)) +
  scale_colour_manual(
    values = syn_region_colours$colours,
    breaks = syn_region_colours$region, na.translate = FALSE
  ) +
  facet_wrap(~fyear, drop = FALSE)

# Come back to this because I don't know that the categories are straightforward
syn_survey_yrs <-
  dat |>
  filter(str_detect(survey_abbrev, "SYN")) |>
  distinct(year, survey_abbrev) #|>
# arrange(year) |>
# rename(region = "survey_abbrev") |>
# complete(year = 2003:2022, fill = list(region = "None", surveyed = FALSE)) |>
# mutate(surveyed = 1) |>
# pivot_wider(id_cols = year, names_from = region, values_from = surveyed) |>
# mutate(region = case_when(
#   `HBLL INS N` == 1 & is.na(`HBLL INS S`) & is.na(None) ~ "HBLL INS N",
#   is.na(`HBLL INS N`) & `HBLL INS S` == 1 & is.na(None) ~ "HBLL INS S",
#   `HBLL INS N` == 1 & `HBLL INS S` == 1 & is.na(None) ~ "Both",
#   is.na(`HBLL INS N`) & is.na(`HBLL INS S`) & None == 1 ~ "No data"
# )) |>
# select(year, region)

fitted_yrs_extra <- min(syn_survey_dat$year):max(syn_survey_dat$year)

syn_nd <-
  syn_grid |>
  rename(area = "cell_area") |>
  make_grid(years = fitted_yrs_extra) |>
  mutate(fyear = as.factor(year))


# Arrowtooth & Bocaccio examples
# ------------------------------------------------------------------------------
spp <- c("arrowtooth flounder", "bocaccio")
sp_dat <- syn_survey_dat |>
  filter(species_common_name %in% spp) |>
  mutate(data_subset = paste(species_common_name, "Stitched N/S", sep = "-")) |>
  group_by(species_common_name) |>
  group_split()

future::plan(future::multisession, workers = 5) # or whatever number
map_func <- furrr::future_map
# map_func <- purrr::map

fits1 <- sp_dat %>%
  # furrr::future_map(fit_models, catch = "density_ppkm2", silent = FALSE) %>%
  # map(fit_models, catch = "density_ppkm2", silent = FALSE, family = sdmTMB::delta_gamma())
  # map(fit_models, catch = "density_ppkm2", silent = FALSE, family = sdmTMB::tweedie())
  map_func(
    fit_models,
    catch = "catch_weight", silent = FALSE,
    family = tweedie(), offset = "trawl_offset", data_subset = "data_subset"
  ) |>
  list_flatten(name_spec = "{inner}")
beep()
future::plan(future::sequential)

fits_cleaned1 <- fits1 |>
  map(check_sanity) # omit plots made from models that did not pass sanity check

saveRDS(fits_cleaned1, "syn-fits.RDS")
fits_cleaned1 <- readRDS("syn-fits.RDS")

preds1 <- get_pred_list(fits_cleaned1, newdata = syn_nd)
indices1 <- get_index_list(pred_list = preds1)
beep()

ids <- names(indices1) %>%
  strsplit(":") %>%
  purrr::map(~ .x[[1]]) %>%
  unlist() %>%
  unique()
lu <- tibble(id = seq_along(ids), desc = ids, order = as.integer(id))

index_df1 <- mk_index_df(indices1) |>
  left_join(lu, by = join_by(desc)) |>
  rename(species = "group")

p1 <-
  ggplot(data = index_df1, aes(x = year, y = est, ymin = lwr, ymax = upr)) +
  geom_pointrange() +
  geom_ribbon(alpha = 0.20, colour = NA) +
  labs(colour = "Sampled region") +
  facet_wrap(species ~ fct_reorder(desc, order), nrow = 4L, scales = "free_y") +
  ggtitle("Stitched N/S - All Years")
p1
