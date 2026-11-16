# BC hard-bottom longline (HBLL outside). North and south alternate years, but
# the alternation flips calendar-year parity after the 2013 gap (N: 2006-2012
# even, 2015-2025 odd), so phase is coded by the region sampled rather than by
# year parity.

library(dplyr)
source(here::here("analysis/R/models.R"))

# From the gfsynopsis-2025 data pull
# (~/src/gfsynopsis-2025/report/data-cache-2026-07/survey-sets.rds, filtered to
# HBLL OUT N / S). This tracked file is the source of truth.
sets <- readRDS(here::here("data-raw/hbll-out-sets.rds")) |>
  drop_na(lglsp_hook_count, count_bait_only, catch_count) |>
  # Hook competition adjustment as in gfsynopsis-2025 (prep_stitch_dat()):
  # offset = log(hooks / adjust), adjust = -log(p) / (1 - p), p = prop. baited hooks
  mutate(
    count_bait_only = replace(count_bait_only, count_bait_only == 0, 1),
    prop_bait_hooks = count_bait_only / lglsp_hook_count,
    hook_adjust_factor = -log(prop_bait_hooks) / (1 - prop_bait_hooks),
    offset = log(lglsp_hook_count / hook_adjust_factor)
  ) |>
  filter(is.finite(offset)) |>
  mutate(species_common_name = replace(
    species_common_name, species_common_name == "north pacific spiny dogfish", "pacific spiny dogfish"
  ))

hbll_phase <- sets |>
  distinct(year, survey_abbrev) |>
  mutate(phase = ifelse(survey_abbrev == "HBLL OUT N", 0.5, -0.5)) |>
  select(year, phase)
stopifnot(!any(duplicated(hbll_phase$year)))

grid0 <- bind_rows(
  gfplot::hbll_n_grid$grid |> mutate(survey_abbrev = "HBLL OUT N"),
  gfplot::hbll_s_grid$grid |> mutate(survey_abbrev = "HBLL OUT S")
) |>
  rename(lon = X, lat = Y) |>
  mutate(area = gfplot::hbll_n_grid$cell_area) |>
  add_utm_columns(c("lon", "lat"), utm_crs = 3156)

longline_species <- c(
  "rougheye/blackspotted rockfish complex",
  "china rockfish",
  "copper rockfish",
  "pacific spiny dogfish",
  "tiger rockfish",
  "lingcod",
  "canary rockfish",
  "quillback rockfish",
  "yelloweye rockfish",
  "silvergray rockfish",
  "spotted ratfish",
  "big skate",
  "rosethorn rockfish",
  "southern rock sole",
  "longnose skate",
  "pacific cod",
  "arrowtooth flounder",
  "pacific halibut"
)
stopifnot(all(longline_species %in% sets$species_common_name))

longline_fit_args <- function(species) {
  dat <- sets |>
    filter(species_common_name == species) |>
    add_utm_columns(c("longitude", "latitude"), utm_crs = 3156)
  mesh <- make_mesh(dat, c("X", "Y"), cutoff = 10)
  list(data = dat, mesh = mesh, response = "catch_count", offset = dat$offset)
}

families <- list(nbinom1 = nbinom1(), nbinom2 = nbinom2())

# NB1 or NB2 by AIC on IID RF, factor(year), preferring fits that pass sanity
# checks. Chosen once per species, up front, so the species x model tasks agree.
choose_family <- function(species) {
  single_thread()
  fit_args <- longline_fit_args(species)
  candidates <- families |>
    map(\(f) tryCatch(do.call(fit_model, c(list("IID RF, factor(year)", family = f), fit_args)), error = \(e) NULL)) |>
    compact()
  if (length(candidates) == 0L) return(NA_character_)
  converged <- compact(lapply(candidates, \(f) fit_ok(f)))
  pool <- if (length(converged) > 0L) converged else candidates
  names(pool)[which.min(map_dbl(pool, AIC))]
}

family_file <- here::here("data-generated/case-study-bc-longline-family.rds")
if (!file.exists(family_file) && !file.exists(here::here("data-generated/case-study-bc-longline.rds"))) {
  future::plan(future::multisession, workers = min(70L, length(longline_species), parallelly::availableCores()))
  family_choice <- furrr::future_map_chr(set_names(longline_species), choose_family, .options = furrr_opts(scheduling = Inf))
  future::plan(future::sequential)
  saveRDS(family_choice, family_file)
}

fit_bc_longline <- function(species, model) {
  single_thread()
  family_name <- readRDS(family_file)[[species]]
  # no fit in either family
  if (is.na(family_name)) return(NULL)
  fit_args <- longline_fit_args(species)
  grid <- replicate_df(grid0, "year", sort(unique(fit_args$data$year)))

  do.call(fit_index_models, c(list(models = model, grid = grid, family = families[[family_name]]), fit_args)) |>
    left_join(hbll_phase, by = "year") |>
    mutate(family = family_name)
}

run_case_study(
  longline_species, fit_bc_longline,
  file = here::here("data-generated/case-study-bc-longline.rds"),
  n_workers = min(70L, parallelly::availableCores())
)
