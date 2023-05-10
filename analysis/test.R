library(tidyverse)
library(gfplot)
library(sdmTMB)
library(beepr)
library(patchwork)

source(here::here('analysis', '00_prep-example-data.R'))
source(here::here('analysis', 'utils.R'))
source(here::here('analysis', 'fit_model_func.R'))

ggplot(data = inside_survey_dat) +
geom_point(aes(x = X, y = Y, colour = survey_abbrev)) +
scale_colour_manual(values = inside_region_colours$colours,
  breaks = inside_region_colours$region, na.translate = FALSE) +
facet_wrap(~ fyear, drop = FALSE