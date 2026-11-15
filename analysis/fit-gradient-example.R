# One cell and seed of the gradient grid, run with the same code so its index
# matches the grid's, but keeping the simulation and fits

source(here::here("analysis/R/gradient-sim.R"))

example <- run_gradient_task(delta = 2, sigma_O = 1, seed = 1, return_fits = TRUE)
stopifnot(!some(example$fits, is.null))
example$fits <- lapply(example$fits, strip_sdmTMB)
saveRDS(example, here::here("data-generated/gradient-example.rds"))
