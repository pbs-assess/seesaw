# Biomass index by subarea (North, South) for the Norwegian survey, from all
# the data (including the full coverage in 2024-2026), to show whether the two
# subareas trend differently.

source(here::here("analysis/R/models.R"))
source(here::here("analysis/R/norway-data.R"))

out_file <- here::here("data-generated/norway-subarea.rds")
subarea_species <- c("Pollack", "Atlantic cod", "Saithe")
subarea_model <- "IID RF, factor(year)"

fit_subarea <- function(species) {
  single_thread()
  dat <- filter(sets, englishname == species)
  # the alternating data's mesh, as in fit-case-study-norway.R, so this is
  # the same fit as its "all" design
  mesh_alt <- make_mesh(filter(dat, is_alternating(dat)), c("X", "Y"), cutoff = 3)
  mesh <- make_mesh(dat, c("X", "Y"), mesh = mesh_alt$mesh)
  fit <- fit_ok(fit_norway_model(subarea_model, dat, mesh))
  if (is.null(fit)) cli::cli_abort("{species}: fit failed sanity checks")
  years <- sort(unique(dat$year))
  sampled <- distinct(dat, year, subarea)
  map(c(North = "North", South = "South"), \(s) {
    grid <- replicate_df(filter(grid0, subarea == s), "year", years)
    get_index_ok(fit, grid)
  }) |>
    bind_rows(.id = "subarea") |>
    # a subarea's index in a year it wasn't sampled is extrapolated, not observed
    semi_join(sampled, by = c("year", "subarea"))
}

if (!file.exists(out_file)) {
  single_thread()
  future::plan(future::multisession, workers = length(subarea_species))
  indexes <- furrr::future_map(set_names(subarea_species), fit_subarea, .options = furrr_opts()) |>
    bind_rows(.id = "species") |>
    mutate(model = subarea_model)
  future::plan(future::sequential)
  saveRDS(indexes, out_file)
}
