n_year <- 12L

base <- list(
  n_year = n_year,
  label = "empty",
  gap_size = 0.25,
  sim_coefs = c(2, 5),
  phi = 7,
  obs_sampled_size = 200L, # or 200 if blocking first
  obs_yrs = list(
    north_yrs = seq(1, n_year - 1, 2), south_yrs = c(seq(2, n_year, 2))
  ),
  region_cutoff = 0.50,
  range = 0.5,
  sigma_O = 1.2,
  sigma_E = 0.4,
  year_arima.sim = list(ar = 0.6),
  year_marginal_sd = 0.3
)

# N here should match number of scenarios below:
sc <- purrr::map(seq_len(23L), ~base)

i <- 1

sc[[i]]$label <- "Base"
sc[[i]]$category <- "Base"
i <- i + 1

# Gap

sc[[i]]$label <- "No gap"
sc[[i]]$category <- "Gap"
sc[[i]]$gap_size <- 0
i <- i + 1

sc[[i]]$label <- "Large gap"
sc[[i]]$category <- "Gap"
sc[[i]]$gap_size <- 0.5
i <- i + 1

sc[[i]]$label <- "Small gap"
sc[[i]]$category <- "Gap"
sc[[i]]$gap_size <- 0.1
i <- i + 1

# Covariate

sc[[i]]$label <- "No covariate"
sc[[i]]$category <- "Abundance gradient"
sc[[i]]$sim_coefs <- c(0, 0)
i <- i + 1

# Observation error

sc[[i]]$label <- "Low observation error"
sc[[i]]$category <- "Observation error"
sc[[i]]$phi <- 3
i <- i + 1

sc[[i]]$label <- "High observation error"
sc[[i]]$category <- "Observation error"
sc[[i]]$phi <- 14
i <- i + 1

# Observation sample size

sc[[i]]$label <- "Low sample size"
sc[[i]]$category <- "Sample size"
sc[[i]]$obs_sampled_size <- 100L # FIXME
i <- i + 1

sc[[i]]$label <- "High sample size"
sc[[i]]$category <- "Sample size"
sc[[i]]$obs_sampled_size <- 400L # FIXME
i <- i + 1

# Year of overlap

sc[[i]]$label <- "Year of overlap (same effort)"
sc[[i]]$category <- "Year overlap"
sc[[i]]$obs_yrs <- list(north_yrs = c(seq(1, n_year - 1, 2), n_year), south_yrs = c(seq(2, n_year, 2)))
i <- i + 1

# Year of overlap (and double effort overall)

sc[[i]]$label <- "Year of overlap (double effort)"
sc[[i]]$category <- "Year overlap"
sc[[i]]$sample_before_split <- TRUE
sc[[i]]$obs_yrs <- list(north_yrs = c(seq(1, n_year - 1, 2), n_year), south_yrs = c(seq(2, n_year, 2)))
i <- i + 1

# Unequal regions

sc[[i]]$label <- "Unequal regions"
sc[[i]]$category <- "Region size"
sc[[i]]$region_cutoff <- 0.25
i <- i + 1

# Range size

sc[[i]]$label <- "Low range"
sc[[i]]$category <- "Range"
sc[[i]]$range <- 0.2  # .25, .5, 1 cod nfld: 180/
i <- i + 1

sc[[i]]$label <- "High range"
sc[[i]]$category <- "Range"
sc[[i]]$range <- 0.8  # .25, .5, 1 cod nfld: 180/
i <- i + 1

# Marginal RF SD

sc[[i]]$label <- "Low sigma O"
sc[[i]]$category <- "Spatial SD"
sc[[i]]$sigma_O <- 0.4
i <- i + 1

# Annual variability stuff N(0, SD)

sc[[i]]$label <- "High year SD"
sc[[i]]$category <- "Annual variability"
sc[[i]]$year_marginal_sd <- 1.0
i <- i + 1

sc[[i]]$label <- "No year SD"
sc[[i]]$category <- "Annual variability"
sc[[i]]$year_marginal_sd <- 0.0001
sc[[i]]$year_arima.sim <- list(ar = 0)
i <- i + 1

# Annual variability AR1 stuff

sc[[i]]$label <- "Year AR(1) rho = 1"
sc[[i]]$category <- "Annual correlation"
sc[[i]]$year_arima.sim <- list(ar = 0.9999)
i <- i + 1

sc[[i]]$label <- "Year AR(1) rho = 0"
sc[[i]]$category <- "Annual correlation"
sc[[i]]$year_arima.sim <- list(ar = 0)
i <- i + 1

# Annual extremes

sc[[i]]$label <- "Years black-swan like"
sc[[i]]$category <- "Annual extremes"
sc[[i]]$heavy_sd_mult <- 6
sc[[i]]$heavy_sd_frac <- 0.1
i <- i + 1

sc[[i]]$label <- "Years heavy-tailed mixture"
sc[[i]]$category <- "Annual extremes"
sc[[i]]$heavy_sd_mult <- 3
sc[[i]]$heavy_sd_frac <- 1/3
i <- i + 1

# SVC trend from north to south; moderate

sc[[i]]$label <- "SVC trend -0.2"
sc[[i]]$category <- "SVC"
sc[[i]]$svc_trend <- -0.2
i <- i + 1

# SVC trend from north to south; very strong

sc[[i]]$label <- "SVC trend -0.6"
sc[[i]]$category <- "SVC"
sc[[i]]$svc_trend <- -0.6
i <- i + 1
