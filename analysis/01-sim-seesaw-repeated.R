# add full coverage sampling
# try with random walk field as truth - does is do OK there??
# try random walk time series - obs error - state space- fit and simulate over and over
# does it flatten??
# free up range
# fix SD on random walk mean?
# fitting a square dataset with a round model - title?
# sept 19 - share range = FALSE!!

# question: same thing happen with missing chunks from one index!?

# To look at:
# - gap matters! how much; related to range I assume
# - **why does RW 'fix' things!?**
# - getting the depth covariate right seems to fix seesaw; *but* bigger CIs on
#   poorly sampled region
# - is coverage of RW still OK!?
# - happens with spatial-only model I think (from Quang)
# - dig in... need to articulate simply what's going wrong!
# - how do I know this is happening in reality (besides the seesaw)?
# - does including a region/north/south covariate fix it: No! build in to
#   routine code?

# What is an index of seesawness?
#
# Things to show:
# - spatial random field messing up
# - RW deviations adapting to fit?
# - fix depth covariates correctly - works?

# quadratic, linear, or breakpoint covariate effect?

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
  tictoc::tic()
  out20 <- do.call(sim_fit_and_index, c(sc[[1]], .seed = 3, make_plots = FALSE))
  tictoc::toc()

  # testing first:
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
    facet_wrap( ~ type,
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

output_file <- "data-generated/sawtooth-sim-apr20.rds"
cache_dir <- "data-generated/sawtooth-sim-apr20-cache"
dir.create("data-generated", showWarnings = FALSE)
dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)

seeds <- seq_len(20)
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
      out <- do.call(sim_fit_and_index, c(sc[[scen_i]], .seed = seed))
      out$label <- scenario_label
      saveRDS(out, cache_file)
    },
    .options = furrr::furrr_options(seed = TRUE, scheduling = 1)
  )
  tictoc::toc()
  future::plan(future::sequential)
}

if (any(!file.exists(tasks$cache_file))) {
  stop("Some cache files were not created.")
}

out_df <- purrr::map_dfr(tasks$cache_file, readRDS)
out_df2 <- left_join(out_df, lu, by = "label")
saveRDS(out_df2, output_file)
