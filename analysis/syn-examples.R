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

library(tidyverse)
library(gfplot)
library(sdmTMB)
library(beepr)
library(patchwork)

source(here::here("analysis", "00_prep-example-data.R"))
source(here::here("analysis", "utils.R"))
source(here::here("analysis", "fit-models.R"))
source(here::here("analysis", "fit-funcs.R"))

theme_set(ggsidekick::theme_sleek())

get_positive_sets <- function(df, response_column) {
  positive_sets <- filter(df, {{ response_column }} != 0)
  prop_pos <- round(nrow(positive_sets) / nrow(df), digits = 3)
  prop_pos
}

pal <- c(
  "QCS + HS" = "#D55E00",
  "QCS + HS + WCVI" = "#CC79A7",
  "QCS + HS + WCHG" = "#009E73",
  "WCHG + WCVI" = "#56B4E9",
  "No data" = "grey50"
)

plot_index <- function(df) {
  df |>
    ggplot(mapping = aes(x = year, y = est, ymin = lwr, ymax = upr, colour = region)) +
    geom_pointrange() +
    geom_ribbon(alpha = 0.20, colour = NA) +
    scale_colour_manual(values = pal) +
    labs(colour = "Sampled region")
}

get_syn_index <- function(sp_dat, fit_file, index_file) {
  # future::plan(future::multisession, workers = 10) # or whatever number
  # map_func <- furrr::future_map
  map_func <- purrr::map
  message("Fitting models")
  fits1 <- sp_dat |>
    map_func(
      fit_models,
      catch = "catch_weight", silent = FALSE,
      family = tweedie(), offset = "trawl_offset", data_subset = "data_subset"
    ) |>
    list_flatten(name_spec = "{inner}")
  # future::plan(future::sequential)

  fits_cleaned <- fits1 |>
    map(check_sanity) # omit plots made from models that did not pass sanity check

  if (!is.null(fit_file)) saveRDS(fits_cleaned, fit_file)

  message("Getting predictions")
  preds1 <- get_pred_list(fits_cleaned, newdata = syn_nd)
  message("Getting index")
  indices1 <- get_index_list(pred_list = preds1)

  if (!is.null(index_file)) saveRDS(fits_cleaned, index_file)

  indices1
}

# Setup inside data and look at survey coverage over time
# ------------------------------------------------------------------------------
syn_survey_dat <- dat |>
  filter(str_detect(survey_abbrev, "SYN")) |>
  filter(!(year %in% c(2003, 2004, 2020))) |> # Use only complete N/S sampling years
  drop_na(trawl_offset) |> # There are some NA values in the offset
  # SOPO only species also tended to be the species that had best trawl coverage
  drop_na(type) # include only SOPO species for a start

# syn_survey_dat |> distinct(species_common_name, type) %>% print(n = 40)

positive_sets <-
  syn_survey_dat %>%
  split(f = list(.$species_common_name, .$year)) %>%
  map(\(x) get_positive_sets(df = x, response_column = density_kgkm2)) %>%
  enframe(name = "species.year", "positive_sets") %>%
  unnest(col = "positive_sets") %>%
  separate(species.year, into = c("species", "year"), sep = "\\.") %>%
  group_by(species) %>%
  summarise(mean_pos_sets = round(mean(positive_sets), digits = 2))

# arrange(positive_sets, mean_pos_sets) %>% print(n = 40)

# syn_region_colours <- tibble(
#   region = c("QCS + HS", "WCHG + WCVI", "QCS + HS + WCVI", "QCS + HS + WCHG", "No data"),
#   colours = c(RColorBrewer::brewer.pal(4L, "Set2"), "grey70")
# )

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
  mutate(pre_f_year = year) |>
  mutate(year_pair = cut(year, seq(min(year), max(year) + 2, 2), right = FALSE)) |>
  mutate(pre_f_year = ifelse(pre_f_year == 2020, 2021, pre_f_year)) |>
  mutate(fyear = as.factor(pre_f_year))

# distinct(syn_nd, year, year_pair)

# What species to include or exclude
# spp <- c("shortraker rockfish") # Not enough data to estimate factor year
# spp <- c("lingcod", "dover sole", "pacific cod", "rex sole", "walleye pollock")
# spp <- c("petrale sole")
spp <- c("walleye pollock")
spp <- c("lingcod")
# spp <- c("lingcod", "dover sole", "pacific cod", "rex sole")
# spp <- c("arrowtooth flounder", "bocaccio", "silvergray rockfish", "yellowtail rockfish")
fit_dir <- here::here("data-outputs", "syn", "fits")
ind_dir <- here::here("data-outputs", "syn", "inds")
fig_dir <- here::here("figs", "syn")
# ------------------------------------------------------------------------------
# Strict alternating
# -----------------------
sp_dat_alt <- syn_survey_dat |>
  filter((species_common_name %in% spp)) |>
  filter(
    !(year == 2007 & region == "SYN WCHG"),
    !(year == 2020),
    !(year == 2021 & region == "SYN WCVI")
  ) |>
  mutate(year_pair = cut(year, seq(min(year), max(year) + 2, 2), right = FALSE)) |>
  mutate(data_subset = paste(species_common_name, "Strict alternating", sep = "-"))

# distinct(sp_dat_alt, year, year_pair)

# sp_dat_both <- syn_survey_dat |>
#   filter((species_common_name %in% spp)) |>
#   mutate(data_subset = paste(species_common_name, "Additional region sampling", sep = "-"))

sp_dat <- sp_dat_alt |>
  group_by(data_subset) |>
  group_split()

system(paste0("mkdir -p ", ind_dir))
system(paste0("mkdir -p ", fit_dir))
ind_list <- get_syn_index(sp_dat,
  fit_file = here::here(fit_dir, paste0(paste(spp, collapse = "-"), "_all-mods", ".RDS")),
  index_file = here::here(ind_dir, paste0(paste(spp, collapse = "-"), "_all-mods", ".RDS"))
)
beep()

ling_fits <- readRDS(here::here(fit_dir, paste0(paste(spp, collapse = "-"), "_all-mods", ".RDS")))

ling_preds <- get_pred_list(ling_fits, newdata = syn_nd)
ling_inds <- get_index_list(ling_preds)
beep()

#ind_list <- readRDS(here::here(fit_dir, paste0(paste(spp, collapse = "-"), "_all-mods", ".RDS")))

ind_df2 <- mk_index_df(ling_inds) |>
  separate(group, into = c("species", "group"), sep = "-") |>
  left_join(syn_survey_yrs, by = "year") |>
  mutate(region = ifelse(year %% 2 == 1, "QCS + HS", "WCHG + WCVI")) |>
  mutate(region = replace(region, year == 2020, "No data")) #|>
  #filter(!(desc %in% c("st IID, (1|region)"))) # This doesn't make sense for strict alternating

saveRDS(ind_df2, here::here(fit_dir, paste0(paste(spp, collapse = "-"), "_all-mods_df", ".RDS")))

ind_df2 <- readRDS(here::here(ind_dir, paste0(paste(spp, collapse = "-"), "_all-mods_df", ".RDS")))

seesaw_order <- ind_df2 |>
  group_by(desc) |>
  mutate(seesaw = mean(log_est[region == "QCS + HS"]) - mean(log_est[region == "WCHG + WCVI"])) |>
  summarise(seesaw = seesaw[1]) |>
  arrange(seesaw) |>
  select(desc) |>
  mutate(order = row_number())

ylim_val <- ind_df2 |>
  filter(desc == "st RW, time-varying RW") |>
  slice(which.max(upr)) |>
  pluck('upr')

left_join(ind_df2, seesaw_order) |>
  plot_index() +
  facet_wrap(~ forcats::fct_reorder(desc, order), scales = "fixed", nrow = 3L) +
  coord_cartesian(ylim = c(0, ylim_val * 1.5)) +
  ggtitle(spp)
# left_join(positive_sets, by = "species") |>
# filter(mean_pos_sets > 0.06)

plot_index(test2) +
  facet_wrap(~desc, scales = "free_y", nrow = 2L) +
  ggtitle("Walleye Pollock")
# facet_wrap(paste0(species, "(", mean_pos_sets, ")") ~ desc, scales = "free_y")


# alternating_ind_list <- get_syn_index(sp_dat,
#   fit_file = here::here(fit_dir, paste0("alternating_", paste(spp, collapse = "-"), "_all-mods", ".RDS")),
#   index_file = here::here(ind_dir, paste0("alternating_", paste(spp, collapse = "-"), "_all-mods", ".RDS")))
alternating_ind_list <- readRDS(here::here(ind_dir, "alternating_all-mods.RDS"))

# alternating_ind_df <- mk_index_df(alternating_ind_list) |>
#   separate(group, into = c("species", "group"), sep = "-") |>
#   left_join(syn_survey_yrs, by = "year") |>
#   mutate(region = ifelse(year %% 2 == 0, "QCS + HS", "WCHG + WCVI")) |>
#   mutate(region = replace(region, year == 2020, "No data")) |>
#   left_join(syn_region_colours, by = "region") |>
#   left_join(positive_sets, by = "species") |>
#   filter(mean_pos_sets > 0.06)

# plot_index(alternating_ind_df) +
#   facet_wrap(paste0(species, "(", mean_pos_sets, ")") ~ desc, scales = "free_y")

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
