library(tidyverse)
library(sdmTMB)
library(patchwork)

source(here::here("analysis", "00-prep-example-data.R"))
source(here::here("analysis", "utils.R"))
source(here::here("analysis", "fit-funcs.R"))
source(here::here("analysis", "fit-models.R"))

mytheme <- function() ggsidekick::theme_sleek() # sometimes I add more layers to themes
theme_set(mytheme())

outside_survey_dat <- dat %>% filter(str_detect(survey_abbrev, "HBLL OUT"))

outside_survey_yrs <- dat %>%
  filter(str_detect(survey_abbrev, "HBLL OUT")) %>%
  distinct(year, survey_abbrev) %>%
  arrange(year) %>%
  rename(region = "survey_abbrev")

fitted_yrs_extra <- seq(min(outside_survey_yrs$year), max(outside_survey_yrs$year))

outside_nd <-
  hbll_outside %>%
  make_grid(years = fitted_yrs_extra) %>%
  mutate(fyear = as.factor(year))

# Yelloweye outside all
# ------------------------------------------------------------------------------

# spp <- c('yelloweye rockfish', 'north pacific spiny dogfish')
spp <- c("yelloweye rockfish")
sp_dat <- dat %>%
  filter(str_detect(survey_abbrev, "HBLL OUT")) %>%
  filter(species_common_name %in% spp) %>%
  mutate(data_subset = paste(species_common_name, "Stitched N/S", sep = "-")) %>%
  group_by(species_common_name) %>%
  group_split()

# future::plan(future::multisession, workers = 5) # or whatever number
# map_func <- furrr::future_map
map_func <- purrr::map

fits1 <- sp_dat %>%
  # furrr::future_map(fit_models, catch = "density_ppkm2", silent = FALSE) %>%
  # map(fit_models, catch = "density_ppkm2", silent = FALSE, family = sdmTMB::delta_gamma())
  # map(fit_models, catch = "density_ppkm2", silent = FALSE, family = sdmTMB::tweedie())
  map_func(
    fit_models,
    catch = "catch_count", silent = FALSE,
    family = sdmTMB::nbinom2(), offset = "hook_offset"
  )

fits_cleaned1 <- list_flatten(fits1, name_spec = "{inner}") %>%
  map(check_sanity)
ids <- names(fits_cleaned1) %>% strsplit(":") %>% purrr::map(~.x[[1]]) %>%
  unlist() %>% unique()

preds1 <- get_pred_list(fits_cleaned1, newdata = outside_nd)
indices1 <- get_index_list(pred_list = preds1) %>%
  setNames(names(fits_cleaned1))

future::plan(future::sequential)
# beep()

# temp:
lu <- tibble(id = seq_along(ids), desc = ids, order = as.integer(id))

index_df <- mk_index_df(indices1) %>%
  left_join(outside_survey_yrs, by = join_by(year)) %>%
  rename(species = "group") %>%
  left_join(lu, by = join_by(desc))

p1 <-
  ggplot(data = index_df, aes(x = year, y = est, ymin = lwr, ymax = upr)) +
  geom_pointrange(aes(colour = region)) +
  geom_ribbon(alpha = 0.20, colour = NA) +
  scale_colour_manual(values = c("#66C2A5", "#FC8D62"), na.translate = FALSE) +
  labs(colour = "Sampled region") +
  facet_wrap(species ~ fct_reorder(desc, order), scales = "free_y") +
  ggtitle("Stitched N/S") +
  scale_y_log10()

p1

# Yelloweye outside split N/S
# ------------------------------------------------------------------------------
if (FALSE) {
  outside_test_dat <- dat %>%
    filter(str_detect(survey_abbrev, "HBLL OUT")) %>%
    filter(species_common_name %in% c("yelloweye rockfish")) %>%
    split(f = .$survey_abbrev) %>%
    map(~ .x %>% mutate(data_subsett = paste(species_common_name, survey_abbrev, sep = "-")))

  future::plan(future::multisession, workers = 5) # or whatever number
  fits2 <- outside_test_dat %>%
    furrr::future_map(fit_models, catch = "density_ppkm2", data_subset = "data_subset") %>%
    list_flatten(name_spec = "{inner}")
  future::plan(future::sequential)

  fits_cleaned2 <- fits2 %>%
    map(., check_sanity) # omit plots made from models that did not pass sanity check

  preds2 <- get_pred_list(fits_cleaned2, newdata = outside_nd)
  indices2 <- get_index_list(pred_list = preds2)
  beep()

  index_df2 <-
    mk_index_df(indices2) %>%
    left_join(outside_survey_yrs) %>%
    separate(group, into = c("species", "region2"), sep = "-")

  p2 <-
    ggplot(data = index_df2, aes(x = year, y = est, ymin = lwr, ymax = upr)) +
    geom_pointrange(aes(colour = region)) +
    geom_ribbon(alpha = 0.20, colour = NA) +
    scale_colour_manual(values = c("#66C2A5", "#FC8D62"), na.translate = FALSE) +
    labs(colour = "Sampled region") +
    facet_wrap(region2 ~ fct_reorder(desc, order), nrow = 2L, scales = "free_y")

  (p1 / p2) + plot_layout(heights = c(1, 2))
}
