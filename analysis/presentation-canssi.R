library(dplyr)
library(ggplot2)

d <- readRDS("~/Downloads/combined_ind_df.RDS")
d <- readRDS("~/Downloads/walleye pollock_all-mods_df.RDS")
# d <- readRDS("data-outputs/syn/inds/lingcod_all-mods_2.RDS")

d <- ind_df2
syn_region_colours <- tibble(
  region = c("QCS + HS", "WCHG + WCVI", "QCS + HS + WCVI", "QCS + HS + WCHG", "No data"),
  colours = c(RColorBrewer::brewer.pal(4L, "Set2"), "grey70"))

# pal <- unname(colorBlindness::availableColors()[-1])
# dput(pal)

# lu <- tribble(
#   ~desc, ~desc2,
#   "st = 'rw'",
#   "st IID covariate",
#   "st IID s(year)",
#   "st IID no covariate as.factor year",
#   "st time-varying RW",
#   "st time-varying RW; fixed 0.1 SD",
#   "spatial time-varying RW",
#   "spatial time-varying AR(1)",
#   "st (1|year)",
#   "spatial only",
#   "st (1 | region)")

pal <- c(
  "QCS + HS" = "#D55E00",
  "QCS + HS + WCVI" = "#CC79A7",
  "QCS + HS + WCHG" = "#009E73",
  "WCHG + WCVI" = "#56B4E9",
  "No data" = "grey50"
)

dd <- d |>
  # filter(species == "lingcod") |>
  # filter(group == "Strict alternating N/S") |>
  filter(!grepl("fixed 0.1", desc), !is.na(est)) |>
  filter(!grepl("only", desc), !is.na(est)) |>
  mutate(desc = gsub(" = ", " ", desc)) |>
  mutate(desc = gsub("st ", "ST ", desc)) |>
  mutate(desc = gsub("AR1", "AR(1)", desc)) |>
  group_by(desc) |>
  mutate(seesaw = abs(mean(log_est[region == "QCS + HS"]) - mean(log_est[region == "WCHG + WCVI"])))

dd <- dd |> group_by(desc) |>
  mutate(seesaw = mean(abs(diff(log_est))))

# because forcats::fct_reorder() isn't working!?
fct_order <- dd |> group_by(desc) |>
  summarise(seesaw = seesaw[1]) |> arrange(seesaw) |> pull(desc)
dd$desc <- factor(dd$desc, levels = fct_order)

dd |>
  ggplot(aes(x = year, y = est, ymin = lwr, ymax = upr, colour = region)) +
  geom_pointrange() +
  geom_ribbon(alpha = 0.20, colour = NA) +
  scale_colour_manual(values = pal) +
  labs(colour = "Sampled region") +
  facet_wrap(~desc, scales = "free_y", nrow=3) +
  # scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, .05))) +
  # scale_y_log10() +
  ylab("Index (relative biomass)") + xlab("Year") +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.x = element_blank(), legend.position = "top") +
  ggtitle(stringr::str_to_title(unique(dd$species)))
ggsave("figs/lingcod-eg.png", width = 10.5, height = 6)

a <- readRDS("~/src/gfsynopsis-2021/report/data-cache-feb-2023/arrowtooth-flounder.rds")
a <- a$survey_sets
wchg <- filter(a, survey_abbrev == "SYN WCHG")

ggplot(wchg, aes(longitude, latitude, size = density_kgpm2)) +
  geom_point(pch = 21) +
  facet_wrap(~year)+
  coord_fixed() +
  scale_size_area(max_size = 9) +
  guides(size = "none")

dd |>
  filter(desc == "ST IID, as.factor(year)") |>
  ggplot(aes(x = year, y = est, ymin = lwr, ymax = upr, colour = region)) +
  geom_pointrange() +
  geom_ribbon(alpha = 0.20, colour = NA) +
  scale_colour_manual(values = pal) +
  labs(colour = "Sampled region") +
  # facet_wrap(~desc, scales = "free_y", nrow=3) +
  # scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, .05))) +
  # scale_y_log10() +
  ylab("Index (relative biomass)") + xlab("Year") +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.x = element_blank(), legend.position = "top") +
  ggtitle(stringr::str_to_title(unique(dd$species)))
ggsave("figs/lingcod-canssi-bad.png", width = 6, height =  4)
ggsave("figs/pollock-canssi-bad.png", width = 6, height =  4)

qb <- readRDS("~/Downloads/qb.rds")

inside_region_colours <- tibble(
  region = c("HBLL INS N", "HBLL INS S", "Both", "No data"),
  colours = c(as.character(pal[c(1, 3, 2)]), "grey70")
)

filter(qb, desc == "st IID, as.factor(year)") |>
  ggplot(aes(x = year, y = est, ymin = lwr, ymax = upr, colour = region)) +
  geom_pointrange() +
  geom_ribbon(alpha = 0.20, colour = NA) +
  # scale_colour_brewer(palette = "Dark2") +
  scale_colour_manual(
    values = inside_region_colours$colours,
    breaks = inside_region_colours$region, na.translate = FALSE
  ) +
  labs(colour = "Sampled region") +
  # facet_wrap(~desc, scales = "free_y", nrow=3) +
  # scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, .05))) +
  # scale_y_log10() +
  ylab("Index (relative biomass)") + xlab("Year") +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.x = element_blank(), legend.position = "top") +
  ggtitle(stringr::str_to_title(unique(qb$species)))
ggsave("figs/qb-canssi.png", width = 6, height = 4)
