# Lognormal prior on sigma_V (the year random walk SD) centred on target_sd.
#
# sdmTMB::lognormal_prior() adds the lognormal density of sigma_V, whose mode
# is exp(meanlog - sdlog^2), so with a fixed meanlog a wider prior would also
# pull sigma_V towards zero. Shifting meanlog by sdlog^2 puts the mode at
# target_sd for any sdlog; this is exactly (up to a constant) a normal penalty
# on log(sigma_V) centred at log(target_sd), the scale sdmTMB optimizes on.
# target_sd sets where the penalized estimate is pulled and sdlog only how
# strongly. (The prior's median is target_sd * exp(sdlog^2).)
sigma_V_prior <- function(target_sd, sdlog) {
  sdmTMB::lognormal_prior(log(target_sd) + sdlog^2, sdlog)
}
