library(tidyverse)
# library(gfplot)
library(sdmTMB)
library(beepr)
library(patchwork)

source(here::here('analysis', '00_prep-example-data.R'))
source(here::here('analysis', 'utils.R'))
source(here::here('analysis', 'fit_model_func.R'))

mytheme <- function() ggsidekick::theme_sleek()  # sometimes I add more layers to themes
theme_set(mytheme())

outside_survey_yrs <- dat %>%
  filter(str_detect(survey_abbrev, "HBLL OUT")) %>%
  distinct(year, survey_abbrev) %>%
  arrange(year) %>%
  rename(region = 'survey_abbrev')

# Yelloweye outside all
# ------------------------------------------------------------------------------
outside_ye_dat <- dat %>%
  filter(str_detect(survey_abbrev, "HBLL OUT")) %>%
  filter(species_common_name %in% c('yelloweye rockfish')) %>%
  mutate(data_subset = paste(species_common_name, 'Stitched N/S', sep = "-")) %>%
  split(f = .$species_common_name)

# future::plan(future::multisession, workers = 5) # or whatever number
fits <- outside_ye_dat %>%
  # furrr::future_map(fit_models, catch = "density_ppkm2", silent = FALSE) %>%
  purrr::map(fit_models, catch = "density_ppkm2", silent = FALSE) %>%
  list_flatten(name_spec = "{inner}")
# future::plan(future::sequential)

fits_cleaned <- fits %>%
  map(check_sanity)  # omit plots made from models that did not pass sanity check

fitted_yrs <- sort(unique(outside_ye_dat[[1]]$year))
fitted_yrs_extra <- min(outside_ye_dat[[1]]$year):max(outside_ye_dat[[1]]$year)

nd <-
  hbll_outside %>% #filter(survey == "HBLL OUT") %>%
  make_grid(years = fitted_yrs) %>%
  mutate(fyear = as.factor(year))
nd_extra_time <-
  hbll_outside %>% #filter(survey == "HBLL OUT") %>%
  make_grid(years = fitted_yrs_extra) %>%
  mutate(fyear = as.factor(year))

preds <- get_pred_list(fits_cleaned, newdata = nd, newdata_extra_time = nd_extra_time)
indices <- get_index_list(pred_list = preds)
beep()

index_df <-
  mk_index_df(indices) %>%
  left_join(outside_survey_yrs)

p1 <-
ggplot(data = index_df, aes(x = year, y = est, ymin = lwr, ymax = upr)) +
  geom_pointrange(aes(colour = region)) +
  geom_ribbon(alpha = 0.20, colour = NA) +
  scale_colour_manual(values = c("#66C2A5", "#FC8D62"), na.translate = FALSE) +
  labs(colour = "Sampled region") +
  facet_wrap(species ~ fct_reorder(desc, order), nrow = 2L, scales = "free_y") +
  ggtitle("Stitched N/S")

p1
# Yelloweye outside split N/S
# ------------------------------------------------------------------------------
if (FALSE) {
outside_test_dat <- dat %>%
  filter(str_detect(survey_abbrev, "HBLL OUT")) %>%
  filter(species_common_name %in% c('yelloweye rockfish')) %>%
  split(f = .$survey_abbrev) %>%
  map(~.x %>% mutate(data_subset = paste(species_common_name, survey_abbrev, sep = "-")))

future::plan(future::multisession, workers = 5) # or whatever number
fits <- outside_test_dat %>%
  furrr::future_map(fit_models, catch = "density_ppkm2", data_subset = "data_subset") %>%
  list_flatten(name_spec = "{inner}")
future::plan(future::sequential)

fits_cleaned <- fits %>%
  map(., check_sanity)  # omit plots made from models that did not pass sanity check

fitted_yrs_extra <- min(outside_survey_yrs$year):max(outside_survey_yrs$year)

nd_extra_time <-
  hbll_outside %>%
  make_grid(years = fitted_yrs_extra) %>%
  mutate(fyear = as.factor(year))

pred_list <-
  fits_cleaned %>%
  map(., function(.x) {
    if (inherits(.x, "sdmTMB")) {
      newdata <- nd_extra_time %>%
        filter(survey %in% unique(.x$data$survey_abbrev),
               year %in% unique(.x$data$year)
        ) %>%
        droplevels()
      out <- predict(.x, newdata = newdata, return_tmb_object = TRUE, extra_time = .x$extra_time)
      out$newdata_input <- newdata
    } else {
      out <- NA
    }
    out
  })

ind_list <- pred_list %>%
  purrr::map(function(.x) {
    if (length(.x) > 1) {
      out <- get_index(.x, bias_correct = TRUE, area = .x$newdata_input$area)
    } else {
      out <- NA  # keep empty fits as visual cue that these did not fit when plotting
    }
    out
  })

index_df2 <-
  mk_index_df(ind_list) %>%
  left_join(outside_survey_yrs) %>%
  separate(species, into = c('species', 'region2'), sep = "-")

p2 <-
ggplot(data = index_df2, aes(x = year, y = est, ymin = lwr, ymax = upr)) +
  geom_pointrange(aes(colour = region)) +
  geom_ribbon(alpha = 0.20, colour = NA) +
  scale_colour_manual(values = c("#66C2A5", "#FC8D62"), na.translate = FALSE) +
  labs(colour = "Sampled region") +
  facet_wrap(region2 ~ fct_reorder(desc, order), nrow = 2L, scales = "free_y")

(p1 / p2) + plot_layout(heights = c(1, 2))
}
