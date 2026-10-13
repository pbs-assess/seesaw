# Omega swap for the QCS + HS biennial self-test: does the real spatial
# structure (the fitted omega, including the HS / QCS contrast) cause the
# spurious seesaw?
#
# As analysis/08-qcs-hs-biennial-self-test.R (operating model without an HS
# term), but each replicate is simulated twice from the same parameter draws
# and seed:
#
#   conditional: omega from its approximate posterior (the real geography)
#   new: omega replaced by an unconditional draw with the same range, sigma_O,
#     and anisotropy (no real geography)
#
# The new spatiotemporal fields and observation draws are identical between
# the pair, so any difference in the excess seesaw is caused by omega. Each
# replicate's realized HS / QCS contrast is kept so the excess seesaw can be
# compared with contrasts that arise by chance in the new fields.

library(dplyr)
library(ggplot2)
library(sdmTMB)

theme_set(ggsidekick::theme_sleek())

source(here::here("analysis/fit-index-models.R"))
source(here::here("analysis/metric-functions.R"))

surveyjoin::cache_data()
surveyjoin::load_sql_data()

survey_names <- c("SYN QCS", "SYN HS")
species_to_fit <- c(
  "pacific cod",
  "spotted ratfish",
  "pacific ocean perch",
  "pacific halibut",
  "pacific spiny dogfish",
  "lingcod",
  "sablefish",
  "dover sole",
  # control: almost no seesaw in analysis/qcs-hs-biennial-screen.R
  "arrowtooth flounder"
  # shortspine thornyhead dropped: its OM encounter AR1 field sits at rho = 1,
  # so each replicate's new epsilon is a near-static field
)
n_reps <- 20L
omega_types <- c("conditional", "new")
ci_levels <- c(0.5, 0.95)
# As in 08: with 11 occasions, 8 gives 4 windows
seesaw_window <- 8L

om_model <- "AR1 RF, factor(year)"
# The only biennial model with a substantial excess seesaw in 08 / 12
est_model <- "IID RF, factor(year)"

family <- delta_gamma(type = "poisson-link")
control <- sdmTMBcontrol(collapse_spatial_variance = TRUE)

# Data -----------------------------------------------------------------------

# Model time as survey occasion (see 08)
prepare_species_data <- function(species) {
  dat <- surveyjoin::get_data(species, regions = "pbs") |>
    mutate(year = lubridate::year(lubridate::ymd(date))) |>
    filter(survey_name %in% survey_names, effort > 0) |>
    select(survey_name, year, lon_start, lat_start, depth_m, effort, catch_weight) |>
    tidyr::drop_na()

  complete_years <- dat |>
    distinct(year, survey_name) |>
    count(year) |>
    filter(n == length(survey_names)) |>
    pull(year)

  dat |>
    filter(year %in% complete_years) |>
    mutate(
      occasion = match(year, sort(unique(year))),
      primary_survey = if_else(occasion %% 2L == 1L, "SYN HS", "SYN QCS"),
      cal_year = year,
      year = occasion
    ) |>
    sdmTMB::add_utm_columns(ll_names = c("lon_start", "lat_start"), utm_crs = 3156)
}

base_grid <- surveyjoin::dfo_synoptic_grid |>
  filter(survey %in% survey_names) |>
  mutate(hs = as.integer(survey == "SYN HS")) |>
  sdmTMB::add_utm_columns(c("lon", "lat"), utm_crs = 3156)

# Models ---------------------------------------------------------------------

# As in 08
fit_model <- function(model, d, mesh) {
  args <- list(
    formula = catch_weight ~ 0 + factor(year),
    data = d,
    mesh = mesh,
    offset = log(d$effort),
    family = family,
    time = "year",
    spatial = "on",
    spatiotemporal = switch(model,
      "IID RF, factor(year)" = "iid",
      "AR1 RF, factor(year)" = "ar1"
    ),
    share_range = TRUE,
    anisotropy = TRUE,
    control = control,
    silent = TRUE
  )
  do.call(sdmTMB, args)
}

# Operating model and paired simulation ----------------------------------------

# dat is prepared in the main session (surveyjoin's database is not safe to
# load from parallel workers)
setup_species <- function(species, dat) {
  RhpcBLASctl::blas_set_num_threads(1L)
  RhpcBLASctl::omp_set_num_threads(1L)

  is_biennial <- dat$survey_name == dat$primary_survey
  grid <- sdmTMB::replicate_df(base_grid, "year", sort(unique(dat$year)))
  mesh <- make_mesh(dat, c("X", "Y"), cutoff = 10)
  mesh_biennial <- make_mesh(dat[is_biennial, ], c("X", "Y"), mesh = mesh$mesh)

  fit_om <- fit_ok(fit_model(om_model, dat, mesh))
  if (is.null(fit_om)) {
    cli::cli_warn("Operating model did not pass sanity checks for {species}")
    return(NULL)
  }

  # One set of parameter draws shared by both omega types: fixed effects at
  # their MLEs, random effects from their approximate posterior (as
  # mle_mvn_samples = "multiple")
  set.seed(42)
  par_draws <- sdmTMB:::.posterior_re_sampler(fit_om)(n_reps)
  # (the draws are unnamed; rows follow the joint parameter vector)
  idx_omega <- which(names(fit_om$tmb_obj$env$last.par.best) == "omega_s")
  stopifnot(length(idx_omega) > 0L)

  # Unconditional omega for each draw (its own range, sigma_O, anisotropy).
  # omega_s is nodes x delta components, column-major as in the par vector;
  # a component collapsed by collapse_spatial_variance is absent from the par
  # vector, so keep only the components that are present
  new_omega <- simulate(
    fit_om,
    nsim = n_reps,
    mcmc_samples = par_draws,
    simulate_re = "spatial",
    return_tmb_report = TRUE,
    seed = 99,
    silent = TRUE
  ) |>
    lapply(\(r) r$omega_s)
  n_nodes <- nrow(new_omega[[1]])
  keep <- seq_len(length(idx_omega) / n_nodes)
  stopifnot(length(idx_omega) %% n_nodes == 0L)
  par_new <- par_draws
  for (i in seq_len(n_reps)) par_new[idx_omega, i] <- as.vector(new_omega[[i]][, keep])
  stopifnot(!isTRUE(all.equal(par_new[idx_omega, 1L], par_draws[idx_omega, 1L])))

  sim_newdata <- bind_rows(
    select(dat, X, Y, year, effort),
    mutate(select(grid, X, Y, year), effort = 1)
  )
  is_tow <- seq_len(nrow(sim_newdata)) <= nrow(dat)

  draws <- list(conditional = par_draws, new = par_new)
  sims <- purrr::map(draws[omega_types], \(p) {
    # One call (and seed) per replicate: the number of gamma draws depends on
    # the encounters, so a single stream would desynchronize the pair after
    # the first replicate
    reports <- lapply(seq_len(n_reps), \(i) {
      simulate(
        fit_om,
        nsim = 1L,
        mcmc_samples = p[, i, drop = FALSE],
        simulate_re = "spatiotemporal",
        newdata = sim_newdata,
        offset = log(sim_newdata$effort),
        return_tmb_report = TRUE,
        seed = 42 + i,
        silent = TRUE
      )[[1L]]
    })
    truth <- purrr::imap_dfr(reports, \(r, rep) {
      eta <- r$eta_i[!is_tow, 1L] + r$eta_i[!is_tow, 2L]
      omega <- r$omega_s_A[!is_tow, 1L] + r$omega_s_A[!is_tow, 2L]
      grid |>
        mutate(density = exp(eta) * area, eta = eta, omega = omega) |>
        group_by(year) |>
        summarise(
          true_est = sum(density),
          # realized log density contrast (HS minus QCS, area-weighted)
          contrast = weighted.mean(eta[hs == 1L], area[hs == 1L]) -
            weighted.mean(eta[hs == 0L], area[hs == 0L]),
          # the part of it from omega
          contrast_omega = weighted.mean(omega[hs == 1L], area[hs == 1L]) -
            weighted.mean(omega[hs == 0L], area[hs == 0L])
        ) |>
        mutate(rep = rep, .before = 1L)
    })
    list(
      catches = lapply(reports, \(r) r$y_i[is_tow, 1L] * r$y_i[is_tow, 2L]),
      epsilon = lapply(reports, \(r) r$epsilon_st),
      truth = truth
    )
  })
  # Pairing check: only omega differs between the two
  stopifnot(identical(sims$conditional$epsilon, sims$new$epsilon))

  list(
    species = species,
    dat = dat,
    is_biennial = is_biennial,
    grid = grid,
    mesh_biennial = mesh_biennial,
    om_pars = bind_rows(tidy(fit_om, "ran_pars", model = 1), tidy(fit_om, "ran_pars", model = 2)),
    occasions = distinct(dat, year, cal_year, primary_survey),
    catches = purrr::map(sims, "catches"),
    true_indexes = purrr::map_dfr(sims, "truth", .id = "omega")
  )
}

# Fit replicates --------------------------------------------------------------

fit_rep <- function(setup, omega, rep) {
  RhpcBLASctl::blas_set_num_threads(1L)
  RhpcBLASctl::omp_set_num_threads(1L)

  d <- setup$dat[setup$is_biennial, ]
  d$catch_weight <- setup$catches[[omega]][[rep]][setup$is_biennial]

  fit <- fit_ok(fit_model(est_model, d, setup$mesh_biennial))
  if (is.null(fit)) return(tibble::tibble())
  get_index_ok(fit, setup$grid) |>
    mutate(species = setup$species, omega = omega, rep = rep, .before = 1L)
}

if (FALSE) {
  tictoc::tic()
  s <- setup_species("spotted ratfish", prepare_species_data("spotted ratfish"))
  tictoc::toc()
  s$true_indexes |> group_by(omega) |> summarise(across(c(contrast, contrast_omega), mean))
  x <- fit_rep(s, "new", 1L)
}

# One cache file per species; delete a file to refit that species
fits_dir <- here::here("data-generated/qcs-hs-omega-swap")
dir.create(fits_dir, showWarnings = FALSE)
fits_file <- \(species) file.path(fits_dir, paste0(gsub(" ", "-", species), ".rds"))

to_fit <- species_to_fit[!file.exists(fits_file(species_to_fit))]
if (length(to_fit) > 0L) {
  dats <- lapply(to_fit, prepare_species_data)
  future::plan(future::multisession, workers = max(1L, min(20L, future::availableCores() - 2L)))
  setups <- furrr::future_map2(
    to_fit, dats, setup_species,
    .options = furrr::furrr_options(seed = TRUE, scheduling = Inf)
  ) |>
    purrr::compact()

  # One species x omega type x replicate per future; each gets only its own
  # species' setup
  tasks <- tidyr::expand_grid(i = seq_along(setups), omega = omega_types, rep = seq_len(n_reps))
  rep_indexes <- furrr::future_pmap(
    list(setup = setups[tasks$i], omega = tasks$omega, rep = tasks$rep), fit_rep,
    .options = furrr::furrr_options(seed = TRUE, scheduling = Inf)
  )
  future::plan(future::sequential)

  for (i in seq_along(setups)) {
    saveRDS(
      list(
        indexes = bind_rows(rep_indexes[tasks$i == i]),
        true_indexes = mutate(setups[[i]]$true_indexes, species = setups[[i]]$species, .before = 1L),
        occasions = mutate(setups[[i]]$occasions, species = setups[[i]]$species, .before = 1L),
        om_pars = mutate(setups[[i]]$om_pars, species = setups[[i]]$species, .before = 1L)
      ),
      fits_file(setups[[i]]$species)
    )
  }
}

swap <- lapply(fits_file(species_to_fit), \(f) if (file.exists(f)) readRDS(f))
indexes <- purrr::map_dfr(swap, "indexes")
true_indexes <- purrr::map_dfr(swap, "true_indexes")
occasions <- purrr::map_dfr(swap, "occasions")

# Metrics ----------------------------------------------------------------------

compared <- indexes |>
  left_join(true_indexes, by = c("species", "omega", "rep", "year")) |>
  left_join(occasions, by = c("species", "year")) |>
  mutate(log_error = log_est - log(true_est))

for (level in ci_levels) {
  q <- qnorm(1 - (1 - level) / 2)
  compared[[paste0("covered_", level * 100)]] <- abs(compared$log_error) <= q * compared$se
}

phase_of <- \(primary_survey) if_else(primary_survey == "SYN HS", 0.5, -0.5)
seesaw_A_max <- function(est, primary_survey, year) {
  max(moving_window(est, year = year, window = seesaw_window, phase = phase_of(primary_survey))$A)
}

metrics_rep <- compared |>
  arrange(species, omega, rep, year) |>
  group_by(species, omega, rep) |>
  summarise(
    rmse = sqrt(mean(log_error^2)),
    coverage_50 = mean(covered_50),
    coverage_95 = mean(covered_95),
    A_max = seesaw_A_max(est, primary_survey, year),
    .groups = "drop"
  )

A_truth <- true_indexes |>
  left_join(occasions, by = c("species", "year")) |>
  arrange(species, omega, rep, year) |>
  group_by(species, omega, rep) |>
  summarise(
    A_max_true = seesaw_A_max(true_est, primary_survey, year),
    contrast = mean(contrast),
    contrast_omega = mean(contrast_omega),
    .groups = "drop"
  )

metrics_rep <- metrics_rep |>
  left_join(A_truth, by = c("species", "omega", "rep")) |>
  mutate(A_max_excess = A_max - A_max_true)

metrics <- metrics_rep |>
  group_by(species, omega) |>
  summarise(
    n_reps = n(),
    abs_contrast = mean(abs(contrast)),
    abs_contrast_omega = mean(abs(contrast_omega)),
    rmse = mean(rmse),
    coverage_50 = mean(coverage_50),
    A_max_excess = median(A_max_excess),
    .groups = "drop"
  )

# Paired effect of the real omega on the excess seesaw (conditional minus new,
# same replicate)
paired <- metrics_rep |>
  select(species, rep, omega, A_max_excess) |>
  tidyr::pivot_wider(names_from = omega, values_from = A_max_excess) |>
  mutate(effect = conditional - new) |>
  group_by(species) |>
  summarise(
    n_pairs = sum(!is.na(effect)),
    effect_median = median(effect, na.rm = TRUE),
    prop_positive = mean(effect > 0, na.rm = TRUE),
    .groups = "drop"
  )

saveRDS(
  list(compared = compared, metrics_rep = metrics_rep, metrics = metrics, paired = paired),
  here::here("data-generated/qcs-hs-omega-swap-metrics.rds")
)

purrr::map_dfr(swap, "om_pars") |>
  select(species, model, term, estimate, std.error) |>
  print(n = Inf)
metrics |>
  mutate(across(where(is.numeric), \(x) round(x, 2))) |>
  print(n = Inf, width = Inf)
print(paired, n = Inf)

# Plots ----------------------------------------------------------------------

omega_labels <- c(conditional = "Fitted omega", new = "New omega")
species_labeller <- labeller(species = \(x) tools::toTitleCase(x))

metrics_rep |>
  mutate(omega = factor(omega, levels = names(omega_labels), labels = omega_labels)) |>
  ggplot(aes(omega, A_max_excess)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey60") +
  geom_line(aes(group = rep), colour = "grey80") +
  geom_point(colour = "grey60") +
  stat_summary(aes(group = 1), fun = median, geom = "line", colour = "black") +
  stat_summary(fun = median, geom = "point", colour = "black") +
  facet_wrap(~species, nrow = 1, labeller = species_labeller) +
  labs(
    x = NULL,
    y = paste0("Max A - true max A (%), ", seesaw_window, "-survey windows"),
    caption = paste0(
      "QCS/HS biennial self-test (", est_model, ") with omega replaced by new draws. ",
      "Grey: paired replicates (same epsilon and seed); black: median."
    )
  ) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(here::here("figs/qcs-hs-omega-swap.pdf"), width = 13, height = 4)

# Excess seesaw against the realized contrast: do chance contrasts in the new
# fields produce a seesaw like the real one?
metrics_rep |>
  mutate(omega = factor(omega, levels = names(omega_labels), labels = omega_labels)) |>
  ggplot(aes(abs(contrast), A_max_excess, colour = omega)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey60") +
  geom_point() +
  scale_colour_manual(values = c("black", "grey60")) +
  scale_y_continuous(trans = "pseudo_log", breaks = c(-100, -10, 0, 10, 30, 100, 300, 1000)) +
  facet_wrap(~species, nrow = 1, labeller = species_labeller) +
  labs(
    x = "|Realized HS - QCS log-density contrast| (replicate mean)",
    y = paste0("Max A - true max A (%), ", seesaw_window, "-survey windows"),
    colour = NULL
  ) +
  theme(legend.position = "bottom")
ggsave(here::here("figs/qcs-hs-omega-swap-contrast.pdf"), width = 13, height = 4.5)
