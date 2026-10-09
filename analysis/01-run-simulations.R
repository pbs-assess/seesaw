library(sdmTMB)
library(ggplot2)
library(dplyr)
source("analysis/funcs.R")
dir.create("figs", showWarnings = FALSE)

# Simulation testing survey stitching with various models -----------------

source(here::here("analysis/simulation-scenarios.R"))
if (any(grepl("empty", purrr::map_chr(sc, "label")))) stop("Too many slots")
labels <- unname(purrr::map_chr(sc, "label"))
names(sc) <- labels
categories <- unname(purrr::map_chr(sc, "category"))
lu <- data.frame(label = labels, category = categories, stringsAsFactors = FALSE)
sc <- purrr::map(sc, ~ {
  .x$label <- NULL
  .x$category <- NULL
  .x
})

if (FALSE) {
  # testing first:
  tictoc::tic()
  out20 <- do.call(sim_fit_and_index, c(sc[[1]], .seed = 1, make_plots = FALSE))
  tictoc::toc()
  out <- do.call(sim_fit_and_index, c(sc[[1]], .seed = 1, make_plots = T, save_plots = T))
  out <- do.call(sim_fit_and_index, c(sc[[12]], .seed = 1))
  actual <- select(out20, year, total, seed, sampled_region) |>
    distinct()
  actual
  out1 <- do.call(sim_fit_and_index, c(sc[[1]], .seed = 1, make_plots = FALSE))
  ggplot(out1, aes(year, est, ymin = lwr, ymax = upr)) +
    ggsidekick::theme_sleek() +
    geom_pointrange(aes(colour = sampled_region)) +
    geom_ribbon(alpha = 0.20, colour = NA) +
    geom_line(
      data = actual, mapping = aes(year, total),
      inherit.aes = FALSE, lty = 2
    ) +
    facet_wrap( ~ model,
      scales = "free_y"
    )
}

Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1"
)

sanitize_scenario_name <- function(x) {
  x <- tolower(x)
  x <- gsub("[^a-z0-9]+", "-", x)
  x <- gsub("(^-+|-+$)", "", x)
  x <- gsub("-+", "-", x)
  ifelse(nchar(x) == 0L, "scenario", x)
}

# Set to TRUE to fit only the model that produces the seesaw (IID fields +
# factor(year)), e.g., to diagnose which scenarios drive it:
seesaw_only <- TRUE
models_to_fit <- if (seesaw_only) "IID RF, factor(year)" else NULL

run_name <- if (seesaw_only) "sawtooth-sim-oct09-seesaw-only" else "sawtooth-sim-oct09"
output_file <- file.path("data-generated", paste0(run_name, ".rds"))
cache_dir <- file.path("data-generated", paste0(run_name, "-cache"))
dir.create("data-generated", showWarnings = FALSE)
dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)

seeds <- seq_len(50L)
scenario_slugs <- make.unique(vapply(names(sc), sanitize_scenario_name, character(1)), sep = "-dup-")
tasks <- tidyr::crossing(
  seed = seeds,
  scen_i = seq_along(sc)
) |>
  mutate(
    scenario_slug = scenario_slugs[scen_i],
    scenario_label = names(sc)[scen_i],
    cache_file = file.path(cache_dir, sprintf("seed-%03d__scenario-%s.rds", seed, scenario_slug))
  )

todo <- tasks |>
  filter(!file.exists(cache_file))

if (nrow(todo) > 0L) {
  NCORES <- future::availableCores()
  workers <- max(1L, min(NCORES - 2L, nrow(todo)))

  future::plan(
    future::multisession,
    workers = workers
  )
  tictoc::tic()
  furrr::future_pwalk(
    todo,
    function(seed, scen_i, scenario_slug, scenario_label, cache_file) {
      out <- do.call(sim_fit_and_index, c(sc[[scen_i]], list(.seed = seed, models = models_to_fit)))
      out$label <- scenario_label
      saveRDS(out, cache_file)
    },
    .options = furrr::furrr_options(seed = TRUE, scheduling = 1)
  )
  tictoc::toc()
  future::plan(future::sequential)
}

out_df <- purrr::map_dfr(tasks$cache_file, readRDS)
out_df2 <- left_join(out_df, lu, by = "label")
saveRDS(out_df2, output_file)
