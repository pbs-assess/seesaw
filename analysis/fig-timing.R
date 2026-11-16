library(dplyr)
library(tidyr)
library(ggplot2)
source(here::here("analysis/R/plotting.R"))

# Self-test fits are the only caches that record wall times for every model
files <- c("QCS + HS" = "self-test-qcs-hs.rds", "NWFSC" = "self-test-nwfsc.rds")
files <- files[file.exists(here::here("data-generated", files))]

times <- purrr::imap(files, \(file, survey) {
  readRDS(here::here("data-generated", file))$indexes |>
    distinct(species, rep, design, model, fit_time, index_time) |>
    mutate(survey = survey)
}) |>
  bind_rows() |>
  pivot_longer(c(fit_time, index_time), names_to = "step", values_to = "minutes") |>
  filter(!is.na(minutes)) |>
  mutate(
    minutes = minutes / 60,
    step = factor(step, c("fit_time", "index_time"), c("Model fitting", "Index calculation")),
    survey = factor(survey, names(files)),
    species = tools::toTitleCase(species)
  )

model_order <- times |>
  filter(step == "Model fitting") |>
  summarise(minutes = median(minutes), .by = model) |>
  arrange(minutes) |>
  pull(model)
times <- mutate(times, model = factor(model, model_order))

g <- ggplot(times, aes(minutes, model)) +
  geom_jitter(aes(colour = species), height = 0.2, alpha = 0.4) +
  stat_summary(fun = median, geom = "point", shape = 23, fill = "black") +
  facet_grid(survey ~ step, scales = "free_x") +
  scale_colour_brewer(palette = "Paired") +
  scale_x_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  labs(x = "Time (minutes)", y = NULL, colour = NULL)

ggsave(here::here("figs/timing.pdf"), g, width = 9, height = 2.2 + 1.7 * length(files))
