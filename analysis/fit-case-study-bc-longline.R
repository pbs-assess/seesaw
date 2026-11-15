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

fit_bc_longline <- function(species) {
  single_thread()
  dat <- sets |>
    filter(species_common_name == species) |>
    add_utm_columns(c("longitude", "latitude"), utm_crs = 3156)
  grid <- replicate_df(grid0, "year", sort(unique(dat$year)))
  mesh <- make_mesh(dat, c("X", "Y"), cutoff = 10)
  fit_args <- list(data = dat, mesh = mesh, response = "catch_count", offset = dat$offset)

  # NB1 or NB2 by AIC on IID RF, factor(year), preferring fits that pass
  # sanity checks
  families <- list(nbinom1 = nbinom1(), nbinom2 = nbinom2())
  candidates <- families |>
    map(\(f) tryCatch(do.call(fit_model, c(list("IID RF, factor(year)", family = f), fit_args)), error = \(e) NULL)) |>
    compact()
  if (length(candidates) == 0L) return(NULL)
  converged <- compact(lapply(candidates, \(f) fit_ok(f)))
  pool <- if (length(converged) > 0L) converged else candidates
  family_name <- names(pool)[which.min(map_dbl(pool, AIC))]

  do.call(fit_index_models, c(list(grid = grid, family = families[[family_name]]), fit_args)) |>
    left_join(hbll_phase, by = "year") |>
    mutate(family = family_name)
}

run_case_study(
  longline_species, fit_bc_longline,
  file = here::here("data-generated/case-study-bc-longline.rds"),
  n_workers = min(length(longline_species), parallelly::availableCores())
)
