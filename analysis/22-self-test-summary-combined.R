# Combined self-test summary: QCS + HS biennial (a) above WCGBTS N / S split (b),
# from the cached metrics saved by 08-qcs-hs-biennial-self-test.R and
# 11-nwfsc-split-self-test.R

library(dplyr)
library(ggplot2)
library(patchwork)
source(here::here("analysis/self-test-funcs.R"))
theme_set(ggsidekick::theme_sleek())

# Full-data operating model name, as in the scripts above
om_model <- "IID RF, factor(year)"

m_qcs <- readRDS(here::here("data-generated/qcs-hs-biennial-self-test-metrics.rds"))
m_nwfsc <- readRDS(here::here("data-generated/nwfsc-split-self-test-metrics.rds"))

full_model <- paste0(om_model, " (full data)")

# Light grid lines; species legends to the right of each row
grid_lines <- list(
  theme(panel.grid.major = element_line(colour = "grey90", linewidth = rel(0.5))),
  guides(colour = guide_legend(ncol = 1)),
  theme(
    legend.position = "right",
    legend.background = element_blank(),
    legend.text = element_text(size = 7),
    legend.key.size = unit(0.3, "cm"),
    panel.spacing.x = unit(0.5, "lines")
  )
)

# Column (facet) names on the top row only
# Moving windows (survey occasions) as in 08 (8: QCS + HS has only 11
# occasions) and 11 (10)
p_qcs <- plot_self_test_summary(m_qcs, full_model, window = 8L, seesaw_metric = "error_variability", seesaw_max = 400) + grid_lines +
  ggtitle("BC QCS-HS split")
p_nwfsc <- plot_self_test_summary(m_nwfsc, full_model, window = 10L, rmse_pad = 0.05, seesaw_metric = "error_variability", seesaw_max = 400) + grid_lines +
  ggtitle("NWFSC split") +
  theme(strip.text.x = element_blank(), strip.background = element_blank())

(p_qcs / p_nwfsc) +
  plot_annotation(tag_levels = "a")

ggsave(here::here("figs/self-test-summary-combined.pdf"), width = 9.9, height = 5.7)
