# Imposed north-south gradient x survey coverage; not a claim that overlap fixes it.
# Run: Rscript analysis/17-minimal-gradient-simulation.R [n_reps] [output_dir] [workers]
# Replicates run in parallel; each worker fits its paired scenarios sequentially.
# Example smoke run: Rscript analysis/17-minimal-gradient-simulation.R 1 /tmp/seesaw-gradient
# Same population/catches across designs; same latent-field seed across gradients.
# "Both regions" means spatially distributed sampling throughout each region,
# not a small shared boundary strip. Full coverage means extent, not a census.
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
library(sdmTMB)
library(dplyr)
library(ggplot2)
source(here::here("analysis/funcs.R"))
source(here::here("analysis/metric-functions.R"))

args <- commandArgs(trailingOnly = TRUE)
n_reps <- if (length(args)) as.integer(args[1]) else 50L
stopifnot(length(n_reps) == 1L, !is.na(n_reps), n_reps > 0L)
out_dir <- if (length(args) >= 2L) args[2] else here::here("data-generated/minimal-gradient")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
n_workers <- if (length(args) >= 3L) as.integer(args[3]) else
  max(1L, min(5L, future::availableCores()/2, n_reps))
stopifnot(length(n_workers) == 1L, !is.na(n_workers), n_workers > 0L)
n_workers <- min(n_workers, n_reps)
rep_dir <- file.path(out_dir, "replicates")
dir.create(rep_dir, recursive = TRUE, showWarnings = FALSE)

# Keep the original process settings; vary only gradient and coverage.
n_year <- 12L
n_tows <- 200L
gradients <- c(-2, 0, 2)
sigma_e <- 1.5 # Set to 0 AND st_type to "off" for a separate static-only check.
st_type <- "iid"
years <- seq_len(n_year)
phase <- ifelse(years %% 2L == 1L, 0.5, -0.5)
allocations <- list(
  alternating = ifelse(phase > 0, n_tows, 0L),
  one_full_year = replace(ifelse(phase > 0, n_tows, 0L), 6L, n_tows / 2L),
  partial_each_year = ifelse(phase > 0, 0.75, 0.25) * n_tows,
  full_each_year = rep(n_tows / 2L, n_year)
)
grid <- expand.grid(
  X = (seq_len(100) - 0.5) / 100,
  Y = (seq_len(100) - 0.5) / 100,
  year = years
)
grid$area <- 1 / 10000
mesh <- make_square_mesh(grid, max_edge = 0.05)
saveRDS(list(
  n_reps = n_reps, n_workers = n_workers, gradients = gradients, allocations = allocations,
  sigma_e = sigma_e, st_type = st_type, range = 0.1, sigma_o = 2,
  phi = 10, n_tows = n_tows, session = sessionInfo()
), file.path(out_dir, "settings.rds"))

run_rep <- function(rep_id) {
  # Keep the original RNG algorithm despite future's parallel RNG streams.
  # sim() and sampling below each set their own replicate-specific seeds.
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  results <- list()
  reference_eta <- NULL
  for (gradient in gradients) {
    # sim() uses gradient * Y. Adjusting the annual mean gives gradient * (Y - .5).
    population <- sim(
      grid, mesh, seed = rep_id, coefs = c(gradient, 0),
      year_mean = 1 - gradient / 2, range = 0.1, sigma_O = 2,
      sigma_E = sigma_e, phi = 10, rho = 0,
      year_arima.sim = list(ar = 0.5), year_marginal_sd = 0.3
    )$sim_dat
    base_eta <- population$eta - gradient * (grid$Y - 0.5)
    if (is.null(reference_eta)) reference_eta <- base_eta
    stopifnot(isTRUE(all.equal(base_eta, reference_eta, tolerance = 1e-8)))
    population$area <- grid$area
    truth <- population |>
      summarise(true_est = sum(mu * area), .by = year) |>
      arrange(year)
    # This includes random-field contrasts; gradient = 0 is not zero contrast.
    realized_contrast <- population |>
      summarise(contrast = mean(eta[region == "north"]) -
        mean(eta[region == "south"]), .by = year) |>
      summarise(contrast = mean(contrast)) |>
      pull(contrast)

    # Shared random order gives nested locations wherever allocations allow.
    set.seed(100000L + rep_id)
    population$order <- runif(nrow(population))
    for (design in names(allocations)) {
      message("Rep ", rep_id, ", gradient ", gradient, ", ", design)
      north_n <- allocations[[design]]
      dat <- population |>
        group_by(year, region) |>
        group_modify(\(d, key) {
          n <- if (key$region == "north") north_n[key$year] else n_tows - north_n[key$year]
          d[order(d$order)[seq_len(n)], ]
        }) |>
        ungroup()
      stopifnot(all(count(dat, year)$n == n_tows))
      # Fixed mesh geometry; only the observation projection changes by design.
      fit_mesh <- make_mesh(dat, c("X", "Y"), mesh = mesh$mesh)
      result <- tryCatch({
        fit <- sdmTMB(
          observed ~ 0 + factor(year), data = dat, mesh = fit_mesh,
          family = tweedie(), time = "year", spatial = "on",
          spatiotemporal = st_type, silent = TRUE,
          control = sdmTMBcontrol(collapse_spatial_variance = TRUE)
        )
        checks <- sanity(fit, gradient_thresh = 0.01, silent = TRUE)
        if (!isTRUE(checks$all_ok)) {
          list(status = "failed_sanity", checks = checks)
        } else {
          pred <- predict(fit, newdata = grid, return_tmb_object = TRUE)
          index <- get_index(pred, area = grid$area, bias_correct = TRUE) |>
            left_join(truth, by = "year") |>
            arrange(year)
          stopifnot(nrow(index) == n_year, identical(as.integer(index$year), years))
          error <- log(index$est / index$true_est)
          beta <- mean(moving_window(
            exp(error), year = years, window = 8L, phase = phase
          )$beta)
          list(status = "ok", index = index, beta = beta,
            rmse = sqrt(mean(error^2)))
        }
      }, error = \(e) list(status = "error", message = conditionMessage(e)))
      result <- c(list(rep = rep_id, gradient = gradient, design = design,
        realized_contrast = realized_contrast), result)
      results[[length(results) + 1L]] <- result
    }
  }
  # Unique checkpoint per worker task; rerunning recomputes requested replicates.
  file <- file.path(rep_dir, sprintf("rep-%03d.rds", rep_id))
  saveRDS(results, file)
  file
}

# Return only checkpoint paths to keep worker transfers small. Restore the
# previous future plan even if a replicate errors; completed checkpoints remain.
run_replicates <- function() {
  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)
  if (n_workers == 1L) {
    future::plan(future::sequential)
  } else {
    future::plan(future::multisession, workers = n_workers)
  }
  message("Running ", n_reps, " replicates with ", n_workers, " worker(s)")
  furrr::future_map(
    seq_len(n_reps), run_rep,
    .options = furrr::furrr_options(seed = TRUE, scheduling = Inf,
      packages = c("sdmTMB", "dplyr"))
  )
}
rep_files <- run_replicates()
results <- do.call(c, lapply(rep_files, readRDS))
saveRDS(results, file.path(out_dir, "results.rds"))

metrics <- bind_rows(lapply(results, \(x) tibble(
  rep = x$rep, gradient = x$gradient, design = x$design, status = x$status,
  realized_contrast = x$realized_contrast,
  beta = if (is.null(x$beta)) NA_real_ else x$beta,
  rmse = if (is.null(x$rmse)) NA_real_ else x$rmse
)))
write.csv(metrics, file.path(out_dir, "metrics.csv"), row.names = FALSE)
summary <- metrics |>
  summarise(
    requested = n(), succeeded = sum(status == "ok"),
    mean_log = mean(beta, na.rm = TRUE),
    mc_se = sd(beta, na.rm = TRUE) / sqrt(succeeded),
    ratio = exp(mean_log), lower = exp(mean_log - 1.96 * mc_se),
    upper = exp(mean_log + 1.96 * mc_se), .by = c(gradient, design)
  )
write.csv(summary, file.path(out_dir, "summary.csv"), row.names = FALSE)
print(summary, n = Inf)

# Pair-specific successful matches; failures in unrelated arms do not filter pairs.
paired <- metrics |>
  filter(status == "ok", design != "alternating") |>
  inner_join(
    metrics |> filter(status == "ok", design == "alternating") |>
      select(rep, gradient, baseline = beta), by = c("rep", "gradient")
  ) |>
  mutate(change = beta - baseline) |>
  summarise(n = n(), mean_change = mean(change), mc_se = sd(change) / sqrt(n),
    .by = c(gradient, design))
write.csv(paired, file.path(out_dir, "paired.csv"), row.names = FALSE)

p <- ggplot(summary, aes(gradient, mean_log, colour = design, group = design)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey60") +
  geom_line() + geom_point() +
  geom_errorbar(aes(ymin = mean_log - 1.96 * mc_se, ymax = mean_log + 1.96 * mc_se)) +
  scale_x_continuous(breaks = gradients) +
  labs(x = "Imposed north-south gradient coefficient",
    y = "Mean log seesaw error (north-primary minus south-primary)", colour = NULL,
    caption = "Same tow budget. Bars: 95% Monte Carlo intervals among successful fits.") +
  theme_bw() + theme(legend.position = "bottom")
ggsave(file.path(out_dir, "seesaw.pdf"), p, width = 8, height = 5)
