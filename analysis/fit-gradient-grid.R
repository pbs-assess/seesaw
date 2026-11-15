# How does the seesaw depend on the north-south contrast in the spatial field
# (Delta) and the spatial SD (sigma_O)?
#
# One cache file per task (delete files to refit). A task that errors isn't
# cached, so a rerun retries it; failed fits are cached as absent models.

library(dplyr)
source(here::here("analysis/R/gradient-sim.R"))

cache_dir <- here::here("data-generated/gradient-grid")
dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)

tasks <- expand_grid(delta = gradient_deltas, sigma_O = gradient_sigma_Os, seed = gradient_seeds) |>
  mutate(cache_file = file.path(
    cache_dir, sprintf("delta-%s__sigma-O-%s__seed-%03d.rds", delta, sigma_O, seed)
  ))
todo <- filter(tasks, !file.exists(cache_file))

if (nrow(todo) > 0L) {
  future::plan(future::multisession, workers = max(1L, min(18L, parallelly::availableCores() - 2L)))
  furrr::future_pwalk(
    todo,
    \(delta, sigma_O, seed, cache_file) {
      out <- tryCatch(run_gradient_task(delta, sigma_O, seed), error = \(e) {
        message("Task delta = ", delta, ", sigma_O = ", sigma_O, ", seed = ", seed, " errored: ", conditionMessage(e))
        NULL
      })
      if (!is.null(out)) saveRDS(out, cache_file)
    },
    .options = furrr_opts(scheduling = Inf)
  )
  future::plan(future::sequential)
}

done <- file.exists(tasks$cache_file)
indexes <- map_dfr(tasks$cache_file[done], readRDS)

fit_status <- tasks |>
  expand_grid(model = gradient_models) |>
  left_join(mutate(distinct(indexes, delta, sigma_O, seed, model), ok = TRUE), by = c("delta", "sigma_O", "seed", "model")) |>
  mutate(status = case_when(
    !file.exists(cache_file) ~ "task error",
    !is.na(ok) ~ "ok",
    .default = "fit failed"
  )) |>
  select(delta, sigma_O, seed, model, status)
n_bad <- sum(fit_status$status != "ok")
if (n_bad > 0L) cli::cli_warn("{n_bad} of {nrow(fit_status)} requested fit(s) failed")

saveRDS(list(indexes = indexes, fit_status = fit_status), here::here("data-generated/gradient-grid.rds"))
