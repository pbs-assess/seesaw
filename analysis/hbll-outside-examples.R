library(tidyverse)
library(gfplot)
library(sdmTMB)
library(beepr)
library(patchwork)

source(here::here('analysis', '00_prep-example-data.R'))
source(here::here('analysis', 'utils.R'))
source(here::here('analysis', 'fit_model_func.R'))

mytheme <- function() ggsidekick::theme_sleek()  # sometimes I add more layers to themes
theme_set(mytheme())

outside_survey_dat <- dat %>% filter(str_detect(survey_abbrev, "HBLL OUT"))

outside_survey_yrs <- dat %>%
  filter(str_detect(survey_abbrev, "HBLL OUT")) %>%
  distinct(year, survey_abbrev) %>%
  arrange(year) %>%
  rename(region = 'survey_abbrev')

fitted_yrs_extra <- min(outside_survey_yrs$year):max(outside_survey_yrs$year)

outside_nd <-
  hbll_outside %>%
  make_grid(years = fitted_yrs_extra) %>%
  mutate(fyear = as.factor(year))


# Yelloweye outside all
# ------------------------------------------------------------------------------
outside_ye_dat <- dat %>%
  filter(str_detect(survey_abbrev, "HBLL OUT")) %>% 
  filter(species_common_name %in% c('yelloweye rockfish')) %>% 
  split(f = .$species_common_name) %>%
  map(~.x %>% mutate(data_subset = paste(species_common_name, 'Stitched N/S', sep = "-")))

future::plan(future::multisession, workers = 5) # or whatever number
fits1 <- outside_ye_dat %>%
  furrr::future_map(fit_models, catch = "density_ppkm2", data_subset = NULL) %>%
  list_flatten(name_spec = "{inner}")
future::plan(future::sequential)

fits_cleaned1 <- fits1 %>%
  map(check_sanity)  # omit plots made from models that did not pass sanity check

preds1 <- get_pred_list(fits_cleaned1, newdata = outside_nd)
indices1 <- get_index_list(pred_list = preds1)
beep()

index_df <-
  mk_index_df(indices1) %>%
  left_join(outside_survey_yrs) %>%
  rename(species = 'group')

p1 <-
ggplot(data = index_df, aes(x = year, y = est, ymin = lwr, ymax = upr)) +
  geom_pointrange(aes(colour = region)) +
  geom_ribbon(alpha = 0.20, colour = NA) +
  scale_colour_manual(values = c("#66C2A5", "#FC8D62"), na.translate = FALSE) +
  labs(colour = "Sampled region") +
  facet_wrap(species ~ fct_reorder(desc, order), nrow = 1L, scales = "free_y") +
  ggtitle("Stitched N/S")

# Yelloweye outside split N/S
# ------------------------------------------------------------------------------
outside_test_dat <- dat %>%
  filter(str_detect(survey_abbrev, "HBLL OUT")) %>%
  filter(species_common_name %in% c('yelloweye rockfish')) %>%
  split(f = .$survey_abbrev) %>%
  map(~.x %>% mutate(data_subset = paste(species_common_name, survey_abbrev, sep = "-")))

future::plan(future::multisession, workers = 5) # or whatever number
fits2 <- outside_test_dat %>%
  furrr::future_map(fit_models, catch = "density_ppkm2", data_subset = "data_subset") %>%
  list_flatten(name_spec = "{inner}")
future::plan(future::sequential)

fits_cleaned2 <- fits2 %>%
  map(check_sanity)  # omit plots made from models that did not pass sanity check

preds2 <- get_pred_list(fits_cleaned2, newdata = outside_nd)
indices2 <- get_index_list(pred_list = preds2)
beep()

index_df2 <-
  mk_index_df(indices2) %>%
  left_join(outside_survey_yrs) %>%
  separate(group, into = c('species', 'region2'), sep = "-")

p2 <-
ggplot(data = index_df2, aes(x = year, y = est, ymin = lwr, ymax = upr)) +
  geom_pointrange(aes(colour = region)) +
  geom_ribbon(alpha = 0.20, colour = NA) +
  scale_colour_manual(values = c("#66C2A5", "#FC8D62"), na.translate = FALSE) +
  labs(colour = "Sampled region") +
  facet_wrap(region2 ~ fct_reorder(desc, order), nrow = 2L, scales = "free_y") +
  ggtitle("Separate N/S")

(p1 / p2) + plot_layout(heights = c(1, 2))
