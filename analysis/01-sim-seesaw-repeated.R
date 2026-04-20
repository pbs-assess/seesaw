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
  out20 <- do.call(sim_fit_and_index, c(sc[[11]], .seed = 1, make_plots = FALSE))
  tictoc::toc()

  # testing first:
  out <- do.call(sim_fit_and_index, c(sc[[1]], .seed = 1, make_plots = T, save_plots = T))

  actual <- select(out20, year, total, seed, sampled_region) |>
    distinct()
  actual

  out1 <- do.call(sim_fit_and_index, c(sc[[1]], .seed = 1, make_plots = FALSE))

  ggplot(out20, aes(year, est, ymin = lwr, ymax = upr)) +
    ggsidekick::theme_sleek() +
    geom_pointrange(aes(colour = sampled_region)) +
    geom_ribbon(alpha = 0.20, colour = NA) +
    geom_line(
      data = actual, mapping = aes(year, total),
      inherit.aes = FALSE, lty = 2
    ) +
    facet_wrap( ~ paste(type, with_depth),
      scales = "free_y"
    )
}

future::plan(future::multisession, workers = 8L)
tictoc::tic()
seeds <- seq_len(16L)
# out_df <- purrr:::map_dfr(seeds[1], function(seed_i) {
out_df <- furrr::future_map_dfr(seeds, function(seed_i) {
  x <- list()
  for (i in seq_along(sc)) {
    sc[[i]]$.seed <- seed_i
    x[[i]] <- do.call(sim_fit_and_index, sc[[i]])
  }
  names(x) <- names(sc)
  bind_rows(x, .id = "label")
# })
}, .options = furrr::furrr_options(seed = TRUE))
tictoc::toc()
out_df2 <- left_join(out_df, lu, by = "label")
dir.create("data-generated", showWarnings = FALSE)
saveRDS(out_df2, "data-generated/sawtooth-sim-may3.rds")
future::plan(future::sequential)

