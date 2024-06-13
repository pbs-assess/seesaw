library(sdmTMB)
library(ggplot2)
library(dplyr)
source("analysis/funcs.R")
dir.create("figs", showWarnings = FALSE)
source(here::here("analysis/scenarios.R"))
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

dir.create("data/generated/", showWarnings = FALSE, recursive = TRUE)

f <- "data/generated/illustrative_data.rda"
# big file, don't cache:
fits <- do.call(sim_fit_and_index, c(sc[[1]], .seed = 3, return_fits = TRUE))

if (!file.exists(f)) { # save a little time
  preds <- do.call(sim_fit_and_index, c(sc[[1]], .seed = 3, return_preds = TRUE))
  sim_dat <- do.call(sim_fit_and_index, c(sc[[1]], .seed = 3, return_sim_dat = TRUE))
  index <- do.call(sim_fit_and_index, c(sc[[1]], .seed = 123))
  save(preds, sim_dat, index, file = f)
} else {
  load(f)
}

filter(index, type == "IID", with_depth == "covariate = FALSE") |>
  ggplot(aes(year, est, ymin = lwr, ymax = upr)) +
  ggsidekick::theme_sleek() +
  geom_pointrange(aes(colour = sampled_region)) +
  geom_ribbon(alpha = 0.20, colour = NA)

dfit <- fits$IID$data
ns <- dfit |>
  select(year, region) |>
  distinct()

d <- preds$RW$data
d |> ggplot(aes(X, Y, fill = epsilon_st)) +
  geom_raster() +
  facet_wrap(~year) +
  scale_fill_gradient2(limits = c(-1.7, 1.7)) +
  coord_fixed(expand = FALSE) +
  ggsidekick::theme_sleek() +
  labs(fill = "Spatiotemporal\neffects") +
  theme(axis.title = element_blank(), axis.ticks = element_blank(), axis.text = element_blank())
ggsave("figs/epsilon_st_rw2.png", width = 8, height = 8)

eps <- list()
yrs <- sort(unique(d$year))
eps[[1]] <- d$epsilon_st[d$year == 1]
for (i in yrs[-1]) {
  eps[[i]] <- d$epsilon_st[d$year == i] - d$epsilon_st[d$year == (i - 1)]
}
d$eps_diff <- do.call(c, eps)

d |> ggplot(aes(X, Y, fill = eps_diff)) +
  geom_raster() +
  facet_wrap(~year) +
  scale_fill_gradient2(limits = c(-1.7, 1.7)) +
  coord_fixed(expand = FALSE) +
  ggsidekick::theme_sleek() +
  labs(fill = "Spatiotemporal\neffects") +
  theme(axis.title = element_blank(), axis.ticks = element_blank(), axis.text = element_blank())
ggsave("figs/epsilon_st_rw_diff2.png", width = 8, height = 8)

d |> ggplot(aes(X, Y, fill = omega_s + epsilon_st)) +
  geom_raster() +
  facet_wrap(~year) +
  scale_fill_gradient2() +
  coord_fixed(expand = FALSE) +
  ggsidekick::theme_sleek() +
  labs(fill = "Spatial +\nspatiotemporal\neffects") +
  theme(axis.title = element_blank(), axis.ticks = element_blank(), axis.text = element_blank())
ggsave("figs/omega_plus_epsilon_st_rw.png", width = 8, height = 8)

d <- preds$IID$data
d |> ggplot(aes(X, Y, fill = epsilon_st)) +
  geom_raster() +
  facet_wrap(~year) +
  scale_fill_gradient2() +
  coord_fixed(expand = FALSE) +
  ggsidekick::theme_sleek() +
  labs(fill = "Spatiotemporal\neffects") +
  theme(axis.title = element_blank(), axis.ticks = element_blank(), axis.text = element_blank())
ggsave("figs/epsilon_st_iid.png", width = 8, height = 8)

d |> ggplot(aes(X, Y, fill = omega_s + epsilon_st)) +
  geom_raster() +
  facet_wrap(~year) +
  scale_fill_gradient2() +
  coord_fixed(expand = FALSE) +
  ggsidekick::theme_sleek() +
  labs(fill = "Spatial +\nspatiotemporal\neffects") +
  theme(axis.title = element_blank(), axis.ticks = element_blank(), axis.text = element_blank())
ggsave("figs/omega_plus_epsilon_st_iid.png", width = 8, height = 8)

LIMS <- c(0.1, 8000)
ggplot2::theme_set(ggsidekick::theme_sleek())
g1 <- ggplot(filter(sim_dat, year == 2), aes(X, Y, fill = exp(eta - epsilon_st))) +
  geom_raster() +
  scale_fill_viridis_c(trans = "log10", limits = LIMS) +
  coord_fixed(expand = FALSE) +
  ggtitle("True spatial effect")
g1
ggsave("figs/omega-true.png", width = 6, height = 5)

d <- preds$IID$data
g2 <- ggplot(filter(d, year == 2), aes(X, Y, fill = exp(est - epsilon_st))) +
  geom_raster() +
  scale_fill_viridis_c(trans = "log10", limits = LIMS) +
  coord_fixed(expand = FALSE) +
  ggtitle("IID fields, as.factor(year)\nestimated spatial effect")
g2
ggsave("figs/omega-iid.png", width = 6, height = 5)

dd <- preds$RW$data
g3 <- ggplot(filter(dd, year == 2), aes(X, Y, fill = exp(est - epsilon_st))) +
  geom_raster() +
  scale_fill_viridis_c(trans = "log10", limits = LIMS) +
  coord_fixed(expand = FALSE) +
  ggtitle("RW fields\nestimated spatial effect")
g3
ggsave("figs/omega-rw.png", width = 6, height = 5)

cowplot::plot_grid(
  g1 + guides(fill = "none"),
  g2 + guides(fill = "none"),
  g3 + guides(fill = "none"),
  nrow = 1L, align = "h"
)
ggsave("figs/omegas-example.png", width = 10, height = 4)

d <- mutate(d, est_stuff = est - epsilon_st)

sim_dat2 <- left_join(sim_dat, select(d, X, Y, year, est_stuff)) |>
  mutate(omega_diff = est_stuff - (eta - epsilon_st))
ggplot(sim_dat2, aes(X, Y, fill = exp(omega_diff))) +
  geom_raster() +
  scale_fill_viridis_c(trans = "log10") +
  theme_light() +
  labs(fill = "Ratio of estimated to\ntrue spatial effects") +
  coord_fixed(expand = FALSE)
ggsave("figs/omega-ratio.png", width = 6, height = 5)

d <- preds$IID$data
ggplot(d, aes(X, Y, fill = exp(est_non_rf))) +
  geom_raster() +
  scale_fill_viridis_c(trans = "log10") +
  facet_wrap(~year) +
  coord_fixed(expand = FALSE) +
  labs(fill = "Year effect") +
  theme(axis.title = element_blank(), axis.ticks = element_blank(), axis.text = element_blank())
ggsave("figs/year-effects-facet.png", width = 6, height = 5)

d |>
  select(-region) |>
  left_join(ns) |>
  select(year, est_non_rf, region) |>
  distinct() |>
  ggplot(aes(year, est_non_rf)) +
  geom_line(colour = "grey70") +
  geom_point(aes(colour = region), size = 4, pch = 20) +
  theme_light() +
  ylab("Year effect") +
  labs(colour = "Region\nsampled")
ggsave("figs/seesaw-years.png", width = 5, height = 3.5)
