library(dplyr)
library(patchwork)
source(here::here("analysis/R/seesaw-metrics.R"))
source(here::here("analysis/R/plotting.R"))

metrics <- c("rmse_ratio", "coverage_50", "seesaw_abs_error", "variability")

# Moving windows in survey occasions: 8 for QCS + HS (only 11 occasions; 10
# would give 2 windows), 10 for the WCGBTS (21 occasions)
plot_self_test <- function(file, window, ...) {
  results <- readRDS(here::here("data-generated", file))
  fit_status <- filter(results$fit_status, design %in% c("full", "alternating"))
  indexes <- results$indexes |>
    semi_join(fit_status, by = c("design", "model")) |>
    common_replicates(fit_status)
  rep_metrics <- self_test_rep_metrics(indexes, window) |>
    mutate(model = if_else(design == "full", paste0(model, " (full data)"), model))
  model_order <- rep_metrics |>
    summarise(rmse = mean(rmse), .by = c(species, model)) |>
    summarise(rmse = mean(rmse), .by = model) |>
    arrange(rmse) |>
    pull(model)
  rep_metrics |>
    mutate(model = factor(model, levels = model_order)) |>
    plot_metric_summary(metrics, y = "model", floor = no_seesaw_error(indexes, window), seesaw_max = 400, ...)
}

row_theme <- list(
  guides(colour = guide_legend(ncol = 1)),
  theme(
    legend.position = "right",
    legend.background = element_blank(),
    legend.text = element_text(size = 7),
    legend.key.size = unit(0.3, "cm"),
    panel.spacing.x = unit(0.5, "lines"),
    plot.margin = margin(5.5, 12, 5.5, 5.5)
  )
)

p_qcs <- plot_self_test("self-test-qcs-hs.rds", window = 8L) + row_theme +
  ggtitle("BC QCS-HS split")
p_nwfsc <- plot_self_test("self-test-nwfsc.rds", window = 10L, rmse_pad = 0.05) + row_theme +
  ggtitle("NWFSC split") +
  theme(strip.text.x = element_blank(), strip.background = element_blank())

g <- (p_qcs / p_nwfsc) +
  plot_annotation(tag_levels = "a")
ggsave(here::here("figs/self-tests.pdf"), g, width = 9.9, height = 5.7)
