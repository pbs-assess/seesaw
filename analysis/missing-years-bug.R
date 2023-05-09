library(tidyverse)
source(here::here('analysis', '00_prep-example-data.R'))

qb <- dat |>
  filter(str_detect(survey_abbrev, "HBLL OUT N")) %>% 
  filter(species_common_name == "quillback rockfish")
qb_mesh <- make_mesh(qb, xy_cols = c("X", "Y"), cutoff = 20)


distinct(qb, year)  # actual missing years: 2007, 2009, 2011, 2013, 2014, 2016, 2018, 2020
missing_years <- sdmTMB:::find_missing_time(qb$year)
# Missing years are identified as: 
# [1] 2014 2016 2018 2020

fit <- sdmTMB(
    density_ppkm2 ~ 0,
    family = tweedie(),
    time_varying = ~1, time_varying_type = "rw",
    data = qb, time = "year", spatiotemporal = "iid", spatial = "on",
    mesh = qb_mesh,
    extra_time = missing_years
  )
