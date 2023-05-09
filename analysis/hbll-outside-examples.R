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

# Try looking at quillback and yelloweye:
# ------------------------------------------------------------------------------
outside_test_dat <- dat %>% 
  filter(str_detect(survey_abbrev, "HBLL OUT")) %>% 
  filter(species_common_name %in% c('yelloweye rockfish')) %>% 
  #split(f = .$survey_abbrev) %>%
  split(f = .$species_common_name) #%>%
  #map(~.x %>% mutate(data_subset = paste(species_common_name, survey_abbrev, sep = "-")))

#future::plan(future::multisession, workers = 5) # or whatever number
fits <- outside_test_dat %>% 
  #map(fit_models, catch = "density_ppkm2", data_subset = "data_subset") %>%
  map(fit_models, catch = "density_ppkm2") %>%
  list_flatten(name_spec = "{inner}")
#future::plan(future::sequential)

# saveRDS(fits, 'fits.RDS')
# fits <- readRDS('fits.RDS')

fits_cleaned <- fits %>%
  map(., check_sanity)  # omit plots made from models that did not pass sanity check

fitted_yrs <- sort(unique(outside_test_dat[[1]]$year))
fitted_yrs_extra <- min(outside_test_dat[[1]]$year):max(outside_test_dat[[1]]$year)
#fitted_yrs_extra <- sort(c(fitted_yrs, sdmTMB:::find_missing_time(outside_test_dat[[1]]$year)))
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
  left_join(distinct(outside_test_dat[[1]], year, region))

ggplot(data = index_df, aes(x = year, y = est, ymin = lwr, ymax = upr)) +
  geom_pointrange(aes(colour = region)) +
  geom_ribbon(alpha = 0.20, colour = NA) +
  scale_colour_manual(values = c("#66C2A5", "#FC8D62"), na.translate = FALSE) +
  labs(colour = "Sampled region") +
  facet_wrap(species ~ fct_reorder(desc, order), nrow = 1L, scales = "free_y")


# Below is broken...

fitted_yrs_extra <- sort(min(dat$year):max(dat$year))
nd_extra_time <-
  hbll_outside %>% #filter(survey == "HBLL OUT") %>% 
  make_grid(years = fitted_yrs_extra) %>% 
  mutate(fyear = as.factor(year))

pred_list <- fits_cleaned %>%
  purrr::map(function(.x, newdata = nd) {
    message(inherits(.x, "sdmTMB"))
    if (inherits(.x, "sdmTMB")) {
      if (!is.null)
      nd <- newdata |> filter(year %in% unique(.x$data$year), survey %in% unique(.x$data$survey_abbrev))
      message(unique(.x$data$survey_abbrev))
      out <- try(predict(.x, newdata = nd, return_tmb_object = TRUE))
    } else {
      out <- NA
    }
    out
  })

ind_list <- pred_list %>%
  purrr::map(function(.x) {
    if (length(.x) > 1) {
      get_index(.x, bias_correct = TRUE, area = .x$data$area)
    } else {
      out <- NA  # keep empty fits as visual cue that these did not fit when plotting
    }
  })

index_df <- ind_list |>
  enframe() |>
  unnest(col = "value") |>
  separate(col = 'name', into = c('id', 'group'), sep = ":") |>
  mutate(id = as.numeric(id)) |>
  right_join(model_lookup)

ggplot(data = index_df, aes(x = year, y = est, ymin = lwr, ymax = upr)) + 
  geom_pointrange(aes(colour = group)) +
  geom_ribbon(alpha = 0.20, colour = NA) +
  scale_colour_manual(values = c("#66C2A5", "#FC8D62"), na.translate = FALSE) + 
  labs(colour = "Sampled region") + 
  facet_wrap(group ~ fct_reorder(desc, order), nrow = 2L, scales = "free_y")
