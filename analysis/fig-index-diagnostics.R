# Not a paper figure: every case-study index, for checking a refit by eye
# (run-all.sh --diagnostics)

library(dplyr)
source(here::here("analysis/R/seesaw-metrics.R"))
source(here::here("analysis/R/plotting.R"))

dir.create(here::here("diagnostics"), showWarnings = FALSE)
seesaw_window <- 10L

plot_indexes <- function(out) {
  mean_A <- out |>
    arrange(species, model, year) |>
    summarise(mean_A = mean(moving_window(est, year, phase, window = seesaw_window)$A), .by = c(species, model))
  out |>
    left_join(mean_A, by = c("species", "model")) |>
    mutate(
      geomean = exp(mean(log(est))), est = est / geomean, lwr = lwr / geomean, upr = upr / geomean,
      .by = c(species, model)
    ) |>
    ggplot(aes(year, log(est), ymin = log(lwr), ymax = log(upr))) +
    geom_ribbon(fill = "grey90") +
    geom_linerange(aes(colour = factor(phase))) +
    geom_point(aes(colour = factor(phase))) +
    scale_colour_brewer(palette = "Dark2") +
    facet_grid(forcats::fct_reorder(model, mean_A) ~ species) +
    labs(x = "Year", y = "Log index (geometric mean 1)", colour = "Phase")
}

for (cs in c("bc-trawl", "goa-bc-trawl", "bc-longline", "norway")) {
  res <- readRDS(here::here("data-generated", paste0("case-study-", cs, ".rds")))
  print(filter(res$fit_status, status != "ok"), n = Inf)
  out <- res$indexes
  designs <- if ("design" %in% names(out)) unique(out$design) else NA
  for (d in designs) {
    out_d <- if (is.na(d)) out else filter(out, design == d)
    plot_indexes(out_d)
    ggsave(
      here::here("diagnostics", paste0(cs, if (!is.na(d)) paste0("-", d), ".pdf")),
      width = 30, height = 15
    )
  }
}
