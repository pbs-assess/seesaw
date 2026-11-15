# QCS + HS self-test. Every design but "full" has the same number of tows per
# occasion.
#
# Monitor with: tail -f data-generated/self-test-qcs-hs/progress.log

library(dplyr)
source(here::here("analysis/R/self-test.R"))
source(here::here("analysis/R/survey-data.R"))
load_surveyjoin()

n_reps <- 30L
n_workers <- max(1L, min(60L, parallelly::availableCores() - 2L))

# The operating model has IID spatiotemporal fields: with AR1, rho can go to
# ~1, so each replicate's new epsilon_st is in effect a second, persistent
# spatial field drawn from the prior, swamping the posterior omega_s. With
# IID, persistent structure stays in omega_s and the new epsilon_st are
# year-to-year noise.
design_models <- c("IID RF, factor(year)", "IID RF, RW year")
fits <- bind_rows(
  tibble(design = "full", model = om_model),
  tibble(design = "alternating", model = names(paper_models)),
  expand_grid(design = c("sparse", "overlap10", "overlap25"), model = design_models)
)

base_grid <- qcs_hs_grid()

setup_species <- function(species, dat) {
  is_sampled <- dat$survey_name == dat$primary_survey
  mesh <- make_mesh(dat, c("X", "Y"), cutoff = 10)
  setup <- sim_setup(
    species = species,
    dat = dat,
    is_sampled = is_sampled,
    grid = replicate_df(base_grid, "year", sort(unique(dat$year))),
    mesh = mesh,
    mesh_alt = make_mesh(dat[is_sampled, ], c("X", "Y"), mesh = mesh$mesh),
    occasions = distinct(dat, year, cal_year, phase)
  )
  om <- fit_operating_model(dat, mesh)
  if (isTRUE(om$om_failed)) return(c(list(species = species), om))
  c(setup, list(fit_om = om$fit))
}

results <- run_self_test(
  dir = here::here("data-generated/self-test-qcs-hs"),
  species = qcs_hs_species,
  n_reps = n_reps,
  prepare_data = qcs_hs_data,
  setup_species = setup_species,
  fits = fits,
  n_workers = n_workers
)
saveRDS(results, here::here("data-generated/self-test-qcs-hs.rds"))
