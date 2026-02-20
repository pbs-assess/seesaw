library(sdmTMB)
library(ggplot2)
theme_set(ggsidekick::theme_sleek())
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
# fits <- do.call(sim_fit_and_index, c(sc[[1]], .seed = 3, return_fits = TRUE))

if (!file.exists(f)) { # save a little time
  preds <- do.call(sim_fit_and_index, c(sc[[1]], .seed = 3, return_preds = TRUE))
  sim_dat <- do.call(sim_fit_and_index, c(sc[[1]], .seed = 3, return_sim_dat = TRUE))
  obs_dat <- do.call(sim_fit_and_index, c(sc[[1]], .seed = 3, return_obs_dat = TRUE))
  index <- do.call(sim_fit_and_index, c(sc[[1]], .seed = 123))
  save(preds, obs_dat, sim_dat, index, file = f)
} else {
  load(f)
}

filter(index, type == "IID", with_depth == "covariate = FALSE") |>
  ggplot(aes(year, est, ymin = lwr, ymax = upr)) +
  ggsidekick::theme_sleek() +
  geom_pointrange(aes(colour = sampled_region)) +
  geom_ribbon(alpha = 0.20, colour = NA)

# dfit <- fits$IID$data
# ns <- dfit |>
#   select(year, region) |>
#   distinct()

ns <- structure(list(year = 1:12, region = structure(c(2L, 1L, 2L,
  1L, 2L, 1L, 2L, 1L, 2L, 1L, 2L, 1L), levels = c("south", "north"
  ), class = "factor")), class = c("grouped_df", "tbl_df", "tbl",
    "data.frame"), row.names = c(NA, -12L), groups = structure(list(
      year = 1:12, .rows = structure(list(1L, 2L, 3L, 4L, 5L, 6L,
        7L, 8L, 9L, 10L, 11L, 12L), ptype = integer(0), class = c("vctrs_list_of",
          "vctrs_vctr", "list"))), row.names = c(NA, -12L), .drop = TRUE, class = c("tbl_df",
            "tbl", "data.frame")))

d <- preds$RW$data
d |> ggplot(aes(X, Y, fill = epsilon_st)) +
  geom_raster() +
  facet_wrap(~year) +
  # scale_fill_gradient2(limits = c(-1.7, 1.7)) +
  scale_fill_gradient2() +
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
  theme(axis.title = element_blank(), axis.ticks = element_blank())
ggsave("figs/omega_plus_epsilon_st_iid.png", width = 8, height = 8)

LIMS <- c(0.1, 8000)
ggplot2::theme_set(ggsidekick::theme_sleek())
g1 <- ggplot(filter(sim_dat, year == 2), aes(X, Y, fill = exp(eta - epsilon_st))) +
  geom_raster() +
  scale_fill_viridis_c(trans = "log10", limits = LIMS) +
  coord_fixed(expand = FALSE) +
  ggtitle("True latent spatial effect") +
  guides(fill = "none") +
  theme(axis.text = element_blank(), axis.ticks = element_blank()) +
  geom_hline(yintercept = 0.5)
g1
ggsave("figs/omega-true.png", width = 6, height = 5)

d <- preds$IID$data
g2 <- ggplot(filter(d, year == 2), aes(X, Y, fill = exp(est))) +
  geom_raster() +
  scale_fill_viridis_c(trans = "log10", limits = LIMS) +
  coord_fixed(expand = FALSE) +
  ggtitle("IID fields, factor(year)\nEstimated latent spatial effect") +
  guides(fill = "none") +
  theme(axis.text = element_blank(), axis.ticks = element_blank()) +
  geom_hline(yintercept = 0.5)
g2
ggsave("figs/omega-iid.png", width = 6, height = 5)

dd <- preds$RW$data
g3 <- ggplot(filter(dd, year == 2), aes(X, Y, fill = exp(est))) +
  geom_raster() +
  scale_fill_viridis_c(trans = "log10", limits = LIMS) +
  coord_fixed(expand = FALSE) +
  ggtitle("RW fields\nEstimated latent spatial effect") +
  guides(fill = "none") +
  theme(axis.text = element_blank(), axis.ticks = element_blank()) +
  geom_hline(yintercept = 0.5)
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
  theme(axis.title = element_blank(), axis.ticks = element_blank())
ggsave("figs/year-effects-facet.png", width = 6, height = 5)

g_ts <- d |>
  select(-region) |>
  left_join(ns) |>
  select(year, est_non_rf, region) |>
  distinct() |>
  mutate(region = stringr::str_to_title(region)) |>
  mutate(region = factor(region, levels = c("North", "South"))) |>
  ggplot(aes(year, exp(est_non_rf))) +
  geom_line(colour = "grey70") +
  geom_point(aes(colour = region), size = 4, pch = 20) +
  theme_light() +
  ylab("Year effect") +
  labs(colour = "Region\nsampled") +
  ggsidekick::theme_sleek() +
  scale_x_continuous(breaks = seq(0, 12, 2)) +
  xlab("Year") +
  scale_colour_brewer(palette = "Set2") +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.1))) +
  theme(legend.position = "right") +
  theme(axis.text = element_blank(), axis.ticks = element_blank()) +
  ggtitle("IID fields, factor(year)\nEstimated year effect")
g_ts
ggsave("figs/seesaw-years.png", width = 5, height = 3.5)

patchwork::wrap_plots(g1, g3, g2, g_ts, ncol = 4)
ggsave("text/figs/sim-illustration.pdf", width = 10, height = 3)

blank_theme_elements <- theme(
  panel.grid.major = element_line(colour = "grey90"),
  # panel.spacing.x = unit(20, "pt"),
  axis.text = element_blank(), axis.ticks = element_blank(),
  axis.title = element_blank(), legend.position = "right"
)


LIM <- range(c(sim_dat$eta), log(obs_dat$observed), na.rm = TRUE))
g_sim <-
  ggplot(filter(sim_dat, year %in% c(1:4)), aes(X, Y, fill = eta)) +
  geom_raster() +
  facet_wrap(vars(year), nrow = 1) +
  scale_fill_viridis_c(limits = c(-5, 8.5)) +
  ggsidekick::theme_sleek() +
  blank_theme_elements +
  coord_equal(expand = FALSE) +
  theme(panel.grid.major = element_line(colour = "grey90")) +
  labs(fill = "Log true\ndensity")
g_sim


g_obs <-
  ggplot(filter(obs_dat, year %in% c(1:4)), aes(X, Y, colour = log(observed))) +
  geom_point() +
  facet_wrap(vars(year), nrow = 1) +
  scale_colour_viridis_c() +
  ggsidekick::theme_sleek() +
  blank_theme_elements +
  coord_equal(expand = FALSE) +
  theme(panel.grid.major = element_line(colour = "grey90")) +
  labs(colour = "Log observed\ndensity")
g_obs

patchwork::wrap_plots(g_sim, g_obs, nrow = 2)
