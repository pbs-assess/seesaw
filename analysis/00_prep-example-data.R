# Data preparations
# ------------------------------------------------------------------------------
# Specify model order and labels for plotting ----------------------------------
model_lookup <- 
  tibble(id = 1:8, 
         desc = c("st = 'rw'", # 6
                  "st IID covariate", # 3
                  "st IID s(year)", # 5
                  "st IID no covariate as.factor year", # 2????
                  "st time_varying RW", # 4
                  "st (1|year)", # 1
                  "spatial only",  # 7
                  "st (1|region)"),  # 8
         order = c(6, 3, 5, 2, 4, 1, 7, 8))  # try matching order of the simulated plots

# Prepare grids ----------------------------------------------------------------
syn_grid <- 
  gfplot::synoptic_grid %>%
  tibble() %>%  # because I accidentally print the full df too often
  dplyr::select(-survey_series_name, -utm_zone, -survey_domain_year) %>% 
  mutate(log_depth = log(depth), region = as.factor(survey))

# This inside grid is corrected for area overlapping with land
# From https://github.com/Blue-Matter/quillback-rockfish/blob/master/data-generated/hbll-inside-grid.rds
  # QUESTION: Is there source code available for generating this grid? 
  #  Should it be included in gfplot?
hbll_inside <- 
  readRDS(here::here('data-raw', 'hbll-inside-grid.rds')) %>%
  sdmTMB::add_utm_columns(c("longitude", "latitude"), utm_crs = 32609) %>% 
  select(survey, depth, area, X, Y) %>%  # block looks like a unique spatial id
  mutate(log_depth = log(depth), region = as.factor(survey))

hbll_n_grid  <- 
  gfplot::hbll_n_grid$grid %>% 
  mutate(survey = "HBLL OUT N")
hbll_s_grid  <- 
  gfplot::hbll_s_grid$grid %>%
  mutate(survey = "HBLL OUT S")

hbll_outside <- bind_rows(hbll_n_grid, hbll_s_grid) %>% 
  as_tibble() %>% 
  rename(longitude = "X", latitude = "Y") %>% 
  add_utm_columns(c("longitude", "latitude"), utm_crs = 32609) %>% 
  mutate(log_depth = log(depth), region = as.factor(survey)) %>% 
  mutate(area = 4)  # CHECK ME: area of hbll outside should be 4?

# Clean survey data ------------------------------------------------------------
dat <- 
  readRDS(here::here('data-raw', 'all_surv_catch.rds')) %>%  # What was the code that actually made this (from SOPO)
  mutate(density_kgkm2 = density_kgpm2 * 1e6, 
         log_depth = log(depth_m), 
         area_swept1 = doorspread_m * (speed_mpm * duration_min), 
         area_swept2 = tow_length_m * doorspread_m, 
         area_swept = ifelse(!is.na(area_swept2), area_swept2, area_swept1)) %>% 
  mutate(trawl_offset = log(area_swept / 1e5)) %>%  # Value used for offset
  #filter(!(year == 2021 & survey_abbrev == "SYN WCVI")) %>%  # this region not usually surveyed in odd years
  sdmTMB::add_utm_columns(c("longitude", "latitude"), utm_crs = 32609) %>% 
  # simplify df columns
  select(#survey_id, trip_id, fishing_event_id, 
         survey_abbrev,
         year, month, day, latitude, longitude, X, Y,
         depth_m, log_depth,
         #species_code, 
         species_common_name, 
         catch_weight, catch_count,
         density_kgkm2, #density_pcpm2, 
         density_ppkm2, 
         area_swept, trawl_offset, hook_count, time_deployed) %>%
  # specify factor variables that will be used in models
  mutate(fyear = as.factor(year), 
         region = as.factor(survey_abbrev)) %>%
  # use complete datasets
  drop_na(depth_m)          # drop rows without depths