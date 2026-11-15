library(dplyr)
source(here::here("analysis/R/seesaw-metrics.R"))
source(here::here("analysis/R/plotting.R"))

# 8-year windows: the series is only 12 years long
m <- gradient_grid_metrics(readRDS(here::here("data-generated/gradient-grid.rds")), window = 8L)
print(m$by_cell, n = Inf)

dodge <- position_dodge(width = 0.15)
g <- ggplot(m$by_cell, aes(delta, seesaw_med, colour = factor(sigma_O))) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey60") +
  geom_linerange(aes(ymin = seesaw_lo, ymax = seesaw_hi), position = dodge) +
  geom_line(position = dodge) +
  geom_point(position = dodge) +
  scale_x_continuous(
    breaks = sort(unique(m$by_cell$delta)),
    labels = \(x) sprintf("%g\n(%.1fx)", x, exp(x)),
    expand = expansion(mult = 0.08)
  ) +
  scale_colour_viridis_d(end = 0.85) +
  facet_wrap(~model) +
  labs(
    x = expression("North-south difference in log density" ~ (Delta)),
    y = "Signed seesaw error (%)\n(north years too high)",
    colour = expression(Spatial ~ SD ~ (sigma[O]))
  )
# cairo_pdf so the plotmath Delta renders (the default pdf device drops it)
ggsave(here::here("figs/gradient-grid.pdf"), g, width = 7, height = 3.5, device = cairo_pdf)
