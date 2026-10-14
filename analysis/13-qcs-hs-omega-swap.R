# Omega swap for the QCS + HS biennial self-test: does the real spatial
# structure (the fitted omega, including the HS / QCS contrast) cause the
# spurious seesaw?
#
# As analysis/08-qcs-hs-biennial-self-test.R (operating model without an HS
# term), but each replicate is simulated four times from the same parameter
# draws and seed, changing only omega:
#
#   conditional: omega from its approximate posterior (the real geography)
#   conditional_flat: conditional with each region's (HS, QCS) mean shifted to
#     the overall mean: the real within-region pattern without the HS / QCS
#     contrast
#   new_contrast: new with each region's mean shifted to that of the paired
#     conditional draw: the real HS / QCS contrast with a random within-region
#     pattern
#   new: omega replaced by an unconditional draw with the same range, sigma_O,
#     and anisotropy (no real geography)
#
# The new spatiotemporal fields and observation draws are identical across
# the four, so any difference in the excess seesaw is caused by omega. Each
# replicate's realized HS / QCS contrast is kept so the excess seesaw can be
# compared with contrasts that arise by chance in the new fields.
#
# The pipeline is in analysis/self-test-funcs.R. Everything is cached per
# replicate, so n_reps can be raised later and only the new replicates are
# run. All four omega types are always simulated (it's cheap), so adding one
# to omega_types later only runs its fits.

source(here::here("analysis/self-test-funcs.R"))
theme_set(ggsidekick::theme_sleek())

surveyjoin::cache_data()
surveyjoin::load_sql_data()

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
  "arrowtooth flounder",
  "shortspine thornyhead"
)
n_reps <- 50L
# Quick check: fewer species and replicates, cached and plotted separately
test_run <- FALSE
if (test_run) {
  species_to_fit <- c("pacific ocean perch", "pacific halibut", "pacific spiny dogfish", "dover sole")
  n_reps <- 10L
}
out_suffix <- if (test_run) "-test" else ""
# Species for the example omega maps (must be in species_to_fit)
map_species <- "pacific halibut"
n_workers <- max(1L, min(60L, future::availableCores() - 2L))
# Omega types to fit: any subset of "conditional", "conditional_flat",
# "new_contrast", "new" (the random-pattern pair is the 2 x 2 check; see the
# header)
omega_types <- c("conditional", "conditional_flat")
# As in 08: with 11 occasions, 8 gives 4 windows
seesaw_window <- 8L

# IID spatiotemporal fields, as in 08: with AR1, rho can go to ~1, so each
# replicate's new epsilon_st is a second, near-static spatial field that the
# omega swap leaves untouched. With IID, persistent structure stays in omega.
om_model <- "IID RF, factor(year)"
# The only biennial model with a substantial excess seesaw in 08 / 12
est_model <- "IID RF, factor(year)"

# Omega types -------------------------------------------------------------------

# Shift omega (one delta component, at the mesh nodes) so each region's
# area-weighted mean on the grid moves exactly to a target. Nodes take the
# region of their nearest grid cell and each region's nodes get one shift; the
# projection mixes the two shifts across the HS / QCS boundary, so they're
# solved from the 2 x 2 system of how each moves both region means.
#
# This equalizes mean log density (omega), not mean density: with different
# within-region variances, a density contrast can remain after flattening
# (rep_truth()'s contrast_density)
make_region_shifter <- function(mesh, grid) {
  xy <- as.matrix(grid[, c("X", "Y")])
  A <- fmesher::fm_basis(mesh$mesh, loc = xy)
  node_hs <- grid$hs[RANN::nn2(xy, mesh$mesh$loc[, 1:2], k = 1)$nn.idx[, 1]]
  region_area <- c(sum(grid$area[grid$hs == 0L]), sum(grid$area[grid$hs == 1L]))
  # rows: QCS, HS
  W <- rbind(
    (grid$hs == 0L) * grid$area / region_area[1],
    (grid$hs == 1L) * grid$area / region_area[2]
  )
  WA <- as.matrix(W %*% A)
  means <- \(omega) as.vector(WA %*% omega)
  # nodes x regions (QCS, HS) indicator; region means move by M %*% shifts
  D <- cbind(node_hs == 0L, node_hs == 1L) * 1
  M <- WA %*% D
  list(
    project = \(omega) as.vector(A %*% omega),
    means = means,
    overall = \(omega) sum(means(omega) * region_area) / sum(region_area),
    shift = \(omega, target) omega + as.vector(D %*% solve(M, target - means(omega)))
  )
}

# The four omega types' parameter vectors for one replicate, all from the same
# draw (see the header)
omega_draws <- function(setup, sampler, rep) {
  pars <- draw_pars(sampler, rep)
  # Unconditional omega with the draw's range, sigma_O, and anisotropy
  # (nodes x delta components)
  new_omega <- simulate(
    setup$fit_om,
    nsim = 1L,
    mcmc_samples = pars,
    simulate_re = "spatial",
    return_tmb_report = TRUE,
    seed = seed_extra + rep,
    silent = TRUE
  )[[1L]]$omega_s
  pars_new <- pars
  for (k in names(setup$omega_idx)) pars_new[setup$omega_idx[[k]], ] <- new_omega[, as.integer(k)]

  pars_flat <- pars
  pars_new_contrast <- pars_new
  for (j in setup$omega_idx) {
    omega_j <- pars[j, 1L]
    pars_flat[j, ] <- setup$shifter$shift(omega_j, rep(setup$shifter$overall(omega_j), 2L))
    pars_new_contrast[j, ] <- setup$shifter$shift(pars_new[j, 1L], setup$shifter$means(omega_j))
  }
  list(conditional = pars, conditional_flat = pars_flat, new_contrast = pars_new_contrast, new = pars_new)
}

# Run ------------------------------------------------------------------------

qcs_hs_base_grid <- qcs_hs_grid()

setup_species <- function(species, dat) {
  setup <- qcs_hs_setup(species, dat, qcs_hs_base_grid)
  om <- fit_operating_model(om_model, dat, setup$mesh)
  if (isTRUE(om$om_failed)) return(c(list(species = species), om))
  setup <- c(setup, list(fit_om = om$fit, ran_pars = om$ran_pars))

  # Positions of omega_s in the joint parameter vector, named by delta
  # component (nodes in order). A component collapsed by
  # collapse_spatial_variance is mapped off and absent from the parameter
  # vector, so each surviving parameter's component comes from the map
  # (nodes x components, column-major; NULL if both are on)
  idx <- which(names(om$fit$tmb_obj$env$last.par.best) == "omega_s")
  if (length(idx) == 0L) {
    return(list(species = species, om_failed = TRUE, error = "both spatial fields collapsed: no omega to swap"))
  }
  n_nodes <- nrow(setup$mesh$mesh$loc)
  omega_map <- om$fit$tmb_obj$env$map$omega_s
  kept <- if (is.null(omega_map)) seq_along(idx) else which(!is.na(omega_map))
  component <- (kept - 1L) %/% n_nodes + 1L
  stopifnot(length(kept) == length(idx), length(idx) == n_nodes * n_distinct(component))
  setup$omega_idx <- split(idx, component)
  setup$shifter <- make_region_shifter(setup$mesh, qcs_hs_base_grid)

  # Example maps: replicate 1's omega types on the grid, summed over the delta
  # components (the log density contribution)
  draws <- omega_draws(setup, sdmTMB:::.posterior_re_sampler(om$fit), 1L)
  stopifnot(!isTRUE(all.equal(draws$new[idx, 1L], draws$conditional[idx, 1L])))
  setup$omega_maps <- purrr::map_dfr(draws, \(p) {
    value <- Reduce(`+`, lapply(setup$omega_idx, \(j) setup$shifter$project(p[j, 1L])))
    mutate(select(qcs_hs_base_grid, X, Y, hs, area), value = value)
  }, .id = "omega")
  setup
}

simulate_rep <- function(setup, sampler, rep) {
  reports <- lapply(omega_draws(setup, sampler, rep), \(p) simulate_pars(setup, p, rep))
  # Pairing check: only omega differs between the omega types
  for (r in reports[-1L]) stopifnot(identical(r$epsilon_st, reports[[1L]]$epsilon_st))
  list(
    catches = lapply(reports, \(r) rep_catch(setup, r)),
    truth = purrr::map_dfr(reports, \(r) rep_truth(setup, r), .id = "arm")
  )
}

if (FALSE) {
  s <- setup_species("spotted ratfish", qcs_hs_data("spotted ratfish"))
  sim <- simulate_rep(s, sdmTMB:::.posterior_re_sampler(s$fit_om), 1L)
  sim$truth |>
    group_by(arm) |>
    summarise(across(c(contrast, contrast_omega), mean))
}

results <- run_self_test(
  dir = here::here(paste0("data-generated/qcs-hs-omega-swap", out_suffix)),
  species = species_to_fit,
  n_reps = n_reps,
  prepare_data = qcs_hs_data,
  setup_species = setup_species,
  simulate_rep = simulate_rep,
  models = tibble::tibble(design = "alternating", model = est_model),
  arms = omega_types,
  n_workers = n_workers
)

omega_maps <- results$setups[[map_species]]$omega_maps
# Arms are the omega types; the true indexes hold all four
compared <- mutate(compare_to_truth(results), omega = arm)
true_indexes <- results$true_indexes |>
  mutate(omega = arm) |>
  filter(omega %in% omega_types)

# Metrics ----------------------------------------------------------------------

# Signed seesaw error: beta (HS-year minus QCS-year log difference, detrended)
# of the estimated index minus that of the true index, averaged over windows.
# A uses |beta| and a max over windows, so noise alone inflates it; noise
# averages to zero here

metrics_rep <- compared |>
  arrange(species, omega, rep, year) |>
  group_by(species, omega, rep) |>
  summarise(
    rmse = sqrt(mean(log_error^2)),
    coverage_50 = mean(covered_50),
    coverage_95 = mean(covered_95),
    A_max = seesaw_A_max(est, phase, year, seesaw_window),
    beta_error = mean(
      seesaw_beta(est, phase, year, seesaw_window) - seesaw_beta(true_est, phase, year, seesaw_window)
    ),
    .groups = "drop"
  )

A_truth <- true_indexes |>
  left_join(results$occasions, by = c("species", "year")) |>
  arrange(species, omega, rep, year) |>
  group_by(species, omega, rep) |>
  summarise(
    A_max_true = seesaw_A_max(true_est, phase, year, seesaw_window),
    contrast = mean(contrast),
    contrast_omega = mean(contrast_omega),
    contrast_density = mean(contrast_density),
    .groups = "drop"
  )

metrics_rep <- metrics_rep |>
  left_join(A_truth, by = c("species", "omega", "rep")) |>
  mutate(
    A_max_excess = A_max - A_max_true,
    # as a percentage: how much higher HS years are than QCS years in the
    # estimated index, relative to the truth
    seesaw_error = 100 * (exp(beta_error) - 1)
  )

metrics <- metrics_rep |>
  group_by(species, omega) |>
  summarise(
    n_reps = n(),
    abs_contrast = mean(abs(contrast)),
    abs_contrast_omega = mean(abs(contrast_omega)),
    abs_contrast_density = mean(abs(contrast_density)),
    rmse = mean(rmse),
    coverage_50 = mean(coverage_50),
    A_max_excess = median(A_max_excess),
    seesaw_error = median(seesaw_error),
    .groups = "drop"
  )

# Paired effects on the excess seesaw (same replicate): of the real contrast
# given the real within-region pattern (fitted minus flat), of the real omega
# (fitted minus new), and of the real contrast given a random within-region
# pattern (new + contrast minus new); only those whose omega types were run
# (on the log scale for the signed error, so effects are symmetric)
comparisons <- list(
  fitted_minus_flat = c("conditional", "conditional_flat"),
  fitted_minus_new = c("conditional", "new"),
  new_contrast_minus_new = c("new_contrast", "new")
) |>
  purrr::keep(\(x) all(x %in% omega_types))
paired_wide <- metrics_rep |>
  select(species, rep, omega, A_max_excess, beta_error) |>
  tidyr::pivot_longer(c(A_max_excess, beta_error), names_to = "metric") |>
  tidyr::pivot_wider(names_from = omega, values_from = value)
paired <- purrr::imap_dfr(comparisons, \(x, name) {
  mutate(paired_wide, comparison = name, effect = .data[[x[1]]] - .data[[x[2]]])
}) |>
  mutate(comparison = factor(comparison, levels = names(comparisons))) |>
  group_by(metric, species, comparison) |>
  summarise(
    n_pairs = sum(!is.na(effect)),
    effect_median = median(effect, na.rm = TRUE),
    prop_positive = mean(effect > 0, na.rm = TRUE),
    .groups = "drop"
  )

saveRDS(
  list(
    compared = compared, metrics_rep = metrics_rep, metrics = metrics, paired = paired,
    failure_rates = results$failure_rates
  ),
  here::here(paste0("data-generated/qcs-hs-omega-swap", out_suffix, "-metrics.rds"))
)

purrr::map_dfr(results$setups, "ran_pars", .id = "species") |>
  select(species, model, term, estimate, std.error) |>
  print(n = Inf)
metrics |>
  mutate(across(where(is.numeric), \(x) round(x, 2))) |>
  print(n = Inf, width = Inf)
print(paired, n = Inf)

# Plots ----------------------------------------------------------------------

omega_labels <- c(
  conditional = "Fitted",
  conditional_flat = "Fitted, contrast removed",
  new_contrast = "New, fitted contrast added",
  new = "New"
)
species_labeller <- labeller(species = \(x) stringr::str_wrap(tools::toTitleCase(x), 12))
seesaw_label <- "Seesaw error (ratio of index error, HS vs. QCS survey years)"
# beta_error is a log ratio; label it as a ratio
ratio_breaks <- c(0.125, 0.25, 0.5, 1, 2, 4, 8)
ratio_labels <- paste0("×", ratio_breaks)

# Orient each species by the sign of its fitted HS - QCS omega contrast (from
# the true indexes, not the fits), so positive is an error in the direction of
# the contrast; the same flip applies to every omega type. Orienting by the
# measured error would make even a null effect look positive.
orientation <- A_truth |>
  filter(omega == "conditional") |>
  group_by(species) |>
  summarise(fitted_contrast = mean(contrast_omega)) |>
  mutate(orientation = if_else(fitted_contrast < 0, -1, 1))
metrics_oriented <- metrics_rep |>
  left_join(orientation, by = "species") |>
  mutate(beta_error_oriented = orientation * beta_error)

# All four omega types (supplement)
metrics_oriented |>
  mutate(omega = factor(omega, levels = names(omega_labels), labels = omega_labels)) |>
  ggplot(aes(omega, beta_error_oriented)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey60") +
  geom_line(aes(group = rep), colour = "grey80") +
  geom_point(colour = "grey60", alpha = 0.5) +
  stat_summary(aes(group = 1), fun = median, geom = "line", colour = "black") +
  stat_summary(fun = median, geom = "point", colour = "black") +
  scale_y_continuous(breaks = log(ratio_breaks), labels = ratio_labels) +
  facet_wrap(~species, nrow = 1, labeller = species_labeller, scales = "free_y") +
  labs(
    x = NULL,
    y = sub(" (", "\n(", seesaw_label, fixed = TRUE),
    caption = paste0(
      "QCS/HS biennial self-test (", est_model, ") with omega modified or replaced. ",
      "Grey: paired replicates (same epsilon and seed); black: median. ",
      "Oriented so above 1 is an error in the direction of each species' fitted HS - QCS omega contrast."
    )
  ) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(here::here(paste0("figs/qcs-hs-omega-swap", out_suffix, ".pdf")), width = 13, height = 4.5, device = cairo_pdf)

# Seesaw error against the realized contrast (not oriented): its sign follows
# the contrast's
metrics_rep |>
  mutate(omega = factor(omega, levels = names(omega_labels), labels = omega_labels)) |>
  ggplot(aes(contrast, beta_error, colour = omega)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey60") +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey60") +
  geom_point() +
  scale_colour_manual(values = setNames(c("black", "#E69F00", "#56B4E9", "grey60"), omega_labels)) +
  scale_y_continuous(breaks = log(ratio_breaks), labels = ratio_labels) +
  facet_wrap(~species, nrow = 1, labeller = species_labeller, scales = "free_y") +
  labs(
    x = "Realized HS - QCS log-density contrast (replicate mean)",
    y = "Seesaw error (ratio of index error, HS / QCS years)",
    colour = NULL
  ) +
  theme(legend.position = "bottom")
ggsave(here::here(paste0("figs/qcs-hs-omega-swap", out_suffix, "-contrast.pdf")), width = 13, height = 4.5, device = cairo_pdf)

# Main figure: (a) example omega with and without the HS / QCS contrast, (b)
# seesaw error for the same two

main_arms <- c(conditional = "Fitted omega", conditional_flat = "HS-QCS contrast removed")
main_colours <- c("#1F6FB4", "grey55")

maps <- omega_maps |>
  filter(omega %in% names(main_arms)) |>
  select(omega, X, Y, hs, area, value) |>
  # centred on the coast-wide mean (the same in both by construction)
  mutate(value = value - weighted.mean(value, area), .by = omega)
# Fitted = contrast removed + difference: the difference is the removed
# regional contrast (two near-flat blocks, smoothed across the boundary)
map_labels <- c(conditional = "Fitted omega", conditional_flat = "Contrast removed", difference = "Difference")
maps <- maps |>
  tidyr::pivot_wider(names_from = omega, values_from = value) |>
  mutate(difference = conditional - conditional_flat) |>
  tidyr::pivot_longer(names(map_labels), names_to = "omega") |>
  mutate(omega = factor(omega, levels = names(map_labels), labels = map_labels))

# Each region's mean, labelled at the region's centre
region_means <- maps |>
  group_by(omega, hs) |>
  summarise(X = mean(X), Y = mean(Y), value = weighted.mean(value, area), .groups = "drop") |>
  mutate(
    value = round(value, 1) + 0, # no "-0.0"
    label = paste0(if_else(hs == 1L, "HS", "QCS"), ": ", if_else(value == 0, "0.0", sprintf("%+.1f", value)))
  )

fill_limit <- quantile(abs(maps$value), 0.99)
p_maps <- ggplot(maps, aes(X, Y)) +
  geom_tile(aes(fill = value), width = 2, height = 2) +
  geom_label(aes(label = label), data = region_means, size = 3) +
  scale_fill_gradient2(
    low = "#2166AC", mid = "grey95", high = "#B2182B",
    limits = c(-1, 1) * fill_limit, oob = scales::squish
  ) +
  facet_wrap(~omega) +
  coord_fixed() +
  labs(
    title = paste0(tools::toTitleCase(map_species), ", one posterior draw"),
    fill = "Spatial\neffect\n(log)", x = NULL, y = NULL
  ) +
  theme(axis.text = element_blank(), axis.ticks = element_blank())

dots <- metrics_oriented |>
  filter(omega %in% names(main_arms), !is.na(beta_error)) |>
  group_by(species, omega, fitted_contrast) |>
  summarise(
    lwr = quantile(beta_error_oriented, 0.25),
    med = median(beta_error_oriented),
    upr = quantile(beta_error_oriented, 0.75),
    .groups = "drop"
  ) |>
  mutate(
    species_label = paste0(tools::toTitleCase(species), " (", sprintf("%+.1f", fitted_contrast), ")"),
    species_label = forcats::fct_reorder(species_label, abs(fitted_contrast)),
    y = as.integer(species_label) + if_else(omega == "conditional", 0.13, -0.13),
    omega = factor(omega, levels = names(main_arms), labels = main_arms)
  )

p_dots <- ggplot(dots, aes(med, y, colour = omega)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey60") +
  geom_linerange(aes(xmin = lwr, xmax = upr)) +
  geom_point() +
  expand_limits(x = log(0.5)) +
  scale_colour_manual(values = main_colours) +
  scale_x_continuous(breaks = log(ratio_breaks), labels = ratio_labels) +
  scale_y_continuous(
    breaks = seq_along(levels(dots$species_label)),
    labels = levels(dots$species_label),
    expand = expansion(add = 0.5)
  ) +
  labs(x = paste0(seesaw_label, ",\noriented by the sign of the fitted omega contrast (shown with each species)"), y = NULL, colour = NULL) +
  theme(legend.position = "bottom")

patchwork::wrap_plots(p_maps, p_dots, ncol = 1, heights = c(2.4, 0.6 + 0.4 * n_distinct(dots$species))) +
  patchwork::plot_annotation(tag_levels = "a")
ggsave(
  here::here(paste0("figs/qcs-hs-omega-swap", out_suffix, "-main.pdf")),
  width = 7, height = 4.2 + 0.4 * n_distinct(dots$species), device = cairo_pdf
)
