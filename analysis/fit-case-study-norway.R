# Norwegian coastal gillnet-fyke net survey (data prep in R/norway-data.R).
# The "alternating" design continues the north/south alternation through
# 2024-2026; "all" uses all the data. Phase follows the alternating design in
# both.

source(here::here("analysis/R/models.R"))
source(here::here("analysis/R/norway-data.R"))

fit_norway <- function(species, model) {
  single_thread()
  dat_all <- filter(sets, englishname == species)
  dat_alt <- filter(dat_all, is_alternating(dat_all))
  mesh_alt <- make_mesh(dat_alt, c("X", "Y"), cutoff = 3)
  mesh_all <- make_mesh(dat_all, c("X", "Y"), mesh = mesh_alt$mesh)
  designs <- list(
    alternating = list(dat = dat_alt, mesh = mesh_alt),
    all = list(dat = dat_all, mesh = mesh_all)
  )
  imap(designs, \(x, design) {
    set_worker_log(label = paste0(species, " [", design, "]"))
    fit_index_models(
      models = model,
      data = x$dat, mesh = x$mesh,
      grid = replicate_df(grid0, "year", sort(unique(x$dat$year))),
      response = "catch_weight",
      family = delta_gamma(type = "poisson-link"),
      covariates = "factor(gear)",
      anisotropy = FALSE
    )
  }) |>
    bind_rows(.id = "design") |>
    left_join(norway_phase, by = "year")
}

# Species with >= 10% positive sets and at least one positive set every year
# in the alternating design
norway_species <- sets |>
  filter(is_alternating(sets)) |>
  summarise(n_pos = sum(catch_weight > 0), n = n(), .by = c(englishname, year)) |>
  summarise(prop_pos = sum(n_pos) / sum(n), min_pos = min(n_pos), .by = englishname) |>
  filter(prop_pos >= 0.1, min_pos >= 1) |>
  pull(englishname) |>
  sort()

run_case_study(
  norway_species, fit_norway,
  file = here::here("data-generated/case-study-norway.rds"),
  n_workers = min(70L, parallelly::availableCores()),
  designs = c("alternating", "all")
)
