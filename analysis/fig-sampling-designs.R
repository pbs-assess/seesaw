# Coverage is left out: with this many designs it is cluttered and
# uninformative.

library(dplyr)
source(here::here("analysis/R/seesaw-metrics.R"))
source(here::here("analysis/R/plotting.R"))

window <- 8L # survey occasions, as in the self-test figure
year_models <- c("factor(year)" = "IID RF, factor(year)", "RW year" = "IID RF, RW year")
design_labels <- c(
  full = "Full data", sparse = "Sparse annual", overlap25 = "Overlap 25%",
  overlap10 = "Overlap 10%", alternating = "Alternating"
)

results <- readRDS(here::here("data-generated/self-test-qcs-hs.rds"))
fit_status <- filter(results$fit_status, model %in% year_models)
indexes <- results$indexes |>
  semi_join(fit_status, by = c("design", "model")) |>
  common_replicates(fit_status)

rep_metrics <- self_test_rep_metrics(indexes, window)
year_model_rows <- \(model) factor(names(year_models)[match(model, year_models)], levels = names(year_models))
rep_metrics <- bind_rows(
  filter(rep_metrics, design != "full"),
  # the full-data fit (IID RF, factor(year)) in both rows
  expand_grid(filter(rep_metrics, design == "full") |> select(-model), model = year_models)
) |>
  mutate(
    year_model = year_model_rows(model),
    design = factor(design_labels[design], levels = design_labels)
  )

g <- plot_metric_summary(
  rep_metrics, c("rmse_ratio", "seesaw_abs_error", "variability"),
  y = "design", row = "year_model", connect = TRUE,
  floor = no_seesaw_error(indexes, window)
)
ggsave(here::here("figs/sampling-designs.pdf"), g, width = 8.1, height = 4.5)
