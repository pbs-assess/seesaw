# Norway case study
# shallow net survey / garn-ruse tokt

# drop crabs for now bc only have counts, need diff family
# could drop 'extra coverage' in 2024-2026 so survey has no spatial overlap,
#    see how much it helped to add 1,2,3 years of overlapping covarage

# data processing:
#   C:\Users\a37789\Documents\coastal-cod\garn-ruse-tokt\code\0-data-garn-ruse.R
#   C:\Users\a37789\Documents\coastal-cod\garn-ruse-tokt\code\1-explore-catch.R

library(sdmTMB)
library(ggplot2)
theme_set(theme_light())
library(dplyr)
casestudy <- 'norway'
# datnames <- c('2013-2023, no spatial overlap',
#               'all years, no spatial overlap',
#               'all years, 2 years spatial overlap')

do_fit <- function(.sp, .survey, .out_dir) {
  safe_sanity <- function(fit) {
    if (!inherits(fit, "sdmTMB")) {
      return(FALSE)
    }
    ok <- tryCatch(
      all(unlist(sanity(fit, gradient_thresh = 0.01))),
      error = function(e) {
        message(.sp, " - sanity check errored: ", conditionMessage(e))
        FALSE
      }
    )
    if (!ok) {
      message(.sp, " - failed sanity checks.")
    }
    ok
  }

  safe_fit <- function(expr) {
    fit <- tryCatch(
      expr,
      error = function(e) {
        message(.sp, " - fit errored: ", conditionMessage(e))
        NULL
      }
    )
    if (safe_sanity(fit)) {
      fit
    } else {
      NULL
    }
  }

  safe_index <- function(fit,grid) {
    tryCatch(
      {
        get_index_split(fit, newdata = grid, nsplit = 3L, bias_correct = TRUE)
      },
      error = function(e) {
        message(.sp, " - index errored: ", conditionMessage(e))
        NULL
      }
    )
  }

  if (.survey == "norway") {
    dat0 <- readRDS(file.path('data-raw','garn_ruse_catch_21sp.rds'))

    datnames <- c('2013-2023, no spatial overlap',
                  'all years, no spatial overlap',
                  'all years, 2 years spatial overlap')
    dat <- vector("list", length(datnames))
    names(dat) <- datnames
    dat[[3]] <- dat0 |>
      filter(commonname == .sp) |>
      mutate(depth_m = depth,
             effort = 1,
             catch_weight = biomass_kg,
             catch_count = catchcount,
             common_name = commonname) |>
      select(haulid, year, lon, lat, depth_m, effort, catch_weight, catch_count, common_name, site, gear) |>
      tidyr::drop_na(lon, lat) |>
      tidyr::drop_na(effort, catch_weight)
    dat[[3]] <- sdmTMB::add_utm_columns(dat[[3]], ll_names = c("lon", "lat"), utm_crs=32632, units='km') # most in UTM zone 32N
    dat[[1]] <- dat[[3]] |> filter(year < 2024)
    dat[[2]] <- dat[[3]] |> filter(year < 2024)
    dat[[2]] <- rbind(dat[[2]],
                      filter(dat[[3]], year == 2024 & lat > 64.9),
                      filter(dat[[3]], year == 2025 & lat < 64.9))

    grid0 <- readRDS(file.path('data-raw','garn_ruse_predgrid.rds'))
    grid <- vector("list", length(datnames))
    names(grid) <- datnames
    grid[[1]] <- grid0 |> sdmTMB::replicate_df("year", unique(dat[[1]]$year))
    grid[[2]] <- grid0 |> sdmTMB::replicate_df("year", unique(dat[[2]]$year))
    grid[[3]] <- grid0 |> sdmTMB::replicate_df("year", unique(dat[[3]]$year))

    mesh <- vector("list", length(datnames))
    mesh[[3]] <- make_mesh(dat[[3]], c("X", "Y"), cutoff = 3)
    mesh[[1]] <- make_mesh(dat[[1]], c("X", "Y"), mesh=mesh[[3]]$mesh)
    mesh[[2]] <- make_mesh(dat[[2]], c("X", "Y"), mesh=mesh[[3]]$mesh)
    # mesh[[2]] <- make_mesh(dat[[2]], c("X", "Y"), cutoff = 3)
  }

  all_yrs <- seq(min(dat0$year), max(dat0$year))

  out_dir1 <- file.path(species_out_dir, .sp)
  dir.create(out_dir1, showWarnings = FALSE)
  tmp <- vector("list", length(datnames))
  for(i in 1:length(datnames)){
    fits <- list()

    base_model <- paste("IID RF, factor(year)",datnames[i],sep=' - ')
    fits[[base_model]] <- safe_fit(sdmTMB(
        catch_weight ~ 0 + factor(year) + factor(gear),
        data = dat[[i]],
        mesh = mesh[[i]],
        offset = log(dat[[i]]$effort),
        family = delta_gamma(type = "poisson-link"),
        time = "year",
        spatial = "on",
        spatiotemporal = "iid",
        share_range = TRUE,
        anisotropy = TRUE,
        silent = FALSE
    ))
    # if (is.null(fits[[base_model]])) {
    #   return(dplyr::tibble())
    # } else {
    #   saveRDS(fits[[base_model]], file.path(out_dir1, paste0(.sp,'-m1-d',i,'.rds')))
    # }
    if(!is.null(fits[[base_model]])) {
      saveRDS(fits[[base_model]], file.path(out_dir1, paste0(.sp,'-m1-d',i,'.rds')))
    }

    # fits[[paste("RW RF",datnames[i],sep=' - ')]] <- safe_fit(update(
    #   fits[[base_model]],
    #   formula. = . ~ factor(gear),
    #   spatiotemporal = "rw",
    #   extra_time = all_yrs
    # ))
    fits[[paste("RW RF",datnames[i],sep=' - ')]] <- safe_fit(sdmTMB(
      catch_weight ~ factor(gear),
      data = dat[[i]],
      mesh = mesh[[i]],
      offset = log(dat[[i]]$effort),
      family = delta_gamma(type = "poisson-link"),
      time = "year",
      spatial = "on",
      spatiotemporal = "rw",
      extra_time = all_yrs,
      share_range = TRUE,
      anisotropy = TRUE,
      silent = FALSE
    ))
    if(!is.null(fits[[paste("RW RF",datnames[i],sep=' - ')]])) {
      saveRDS(fits[[paste("RW RF",datnames[i],sep=' - ')]], file.path(out_dir1, paste0(.sp,'-m2-d',i,'.rds')))
    }

    # fits[[paste("AR1 RF",datnames[i],sep=' - ')]] <- safe_fit(update(
    #   fits[[base_model]],
    #   formula. = . ~ factor(gear),
    #   spatiotemporal = "ar1",
    #   extra_time = all_yrs
    # ))
    fits[[paste("AR1 RF",datnames[i],sep=' - ')]] <- safe_fit(sdmTMB(
      catch_weight ~ factor(gear),
      data = dat[[i]],
      mesh = mesh[[i]],
      offset = log(dat[[i]]$effort),
      family = delta_gamma(type = "poisson-link"),
      time = "year",
      spatial = "on",
      spatiotemporal = "ar1",
      extra_time = all_yrs,
      share_range = TRUE,
      anisotropy = TRUE,
      silent = FALSE
    ))
    if(!is.null(fits[[paste("AR1 RF",datnames[i],sep=' - ')]])) {
      saveRDS(fits[[paste("AR1 RF",datnames[i],sep=' - ')]], file.path(out_dir1, paste0(.sp,'-m3-d',i,'.rds')))
    }

    # fits[[paste("RW RF, RW year",datnames[i],sep=' - ')]] <- safe_fit(update(
    #   fits[[base_model]],
    #   formula. = . ~ factor(gear),
    #   spatiotemporal = "rw",
    #   time_varying = ~1,
    #   time_varying_type = "rw0",
    #   priors = sdmTMB::sdmTMBpriors(sigma_V = sdmTMB::gamma_cv(0.3, 0.5)),
    #   extra_time = all_yrs
    # ))
    fits[[paste("RW RF, RW year",datnames[i],sep=' - ')]] <- safe_fit(sdmTMB(
      catch_weight ~ factor(gear),
      data = dat[[i]],
      mesh = mesh[[i]],
      offset = log(dat[[i]]$effort),
      family = delta_gamma(type = "poisson-link"),
      time = "year",
      spatial = "on",
      spatiotemporal = "rw",
        time_varying = ~1,
        time_varying_type = "rw0",
        priors = sdmTMB::sdmTMBpriors(sigma_V = sdmTMB::gamma_cv(0.3, 0.5)),
      extra_time = all_yrs,
      share_range = TRUE,
      anisotropy = TRUE,
      silent = FALSE
    ))
    if(!is.null(fits[[paste("RW RF, RW year",datnames[i],sep=' - ')]])) {
      saveRDS(fits[[paste("RW RF, RW year",datnames[i],sep=' - ')]], file.path(out_dir1, paste0(.sp,'-m4-d',i,'.rds')))
    }

    # fits[[paste("AR1 RF, RW year",datnames[i],sep=' - ')]] <- safe_fit(update(
    #   fits[[base_model]],
    #   formula. = . ~ factor(gear),
    #   spatiotemporal = "ar1",
    #   time_varying = ~1,
    #   time_varying_type = "rw0",
    #   priors = sdmTMB::sdmTMBpriors(sigma_V = sdmTMB::gamma_cv(0.3, 0.5)),
    #   extra_time = all_yrs
    # ))
    fits[[paste("AR1 RF, RW year",datnames[i],sep=' - ')]] <- safe_fit(sdmTMB(
      catch_weight ~ factor(gear),
      data = dat[[i]],
      mesh = mesh[[i]],
      offset = log(dat[[i]]$effort),
      family = delta_gamma(type = "poisson-link"),
      time = "year",
      spatial = "on",
      spatiotemporal = "ar1",
      time_varying = ~1,
      time_varying_type = "rw0",
      priors = sdmTMB::sdmTMBpriors(sigma_V = sdmTMB::gamma_cv(0.3, 0.5)),
      extra_time = all_yrs,
      share_range = TRUE,
      anisotropy = TRUE,
      silent = FALSE
    ))
    if(!is.null(fits[[paste("AR1 RF, RW year",datnames[i],sep=' - ')]])) {
      saveRDS(fits[[paste("AR1 RF, RW year",datnames[i],sep=' - ')]], file.path(out_dir1, paste0(.sp,'-m5-d',i,'.rds')))
    }

    # fits[[paste("IID RF, RW year",datnames[i],sep=' - ')]] <- safe_fit(update(
    #   fits[[base_model]],
    #   formula. = . ~ factor(gear),
    #   spatiotemporal = "iid",
    #   time_varying = ~1,
    #   time_varying_type = "rw0",
    #   priors = sdmTMB::sdmTMBpriors(sigma_V = sdmTMB::gamma_cv(0.3, 0.5)),
    #   extra_time = all_yrs
    # ))
    fits[[paste("IID RF, RW year",datnames[i],sep=' - ')]] <- safe_fit(sdmTMB(
      catch_weight ~ factor(gear),
      data = dat[[i]],
      mesh = mesh[[i]],
      offset = log(dat[[i]]$effort),
      family = delta_gamma(type = "poisson-link"),
      time = "year",
      spatial = "on",
      spatiotemporal = "iid",
      time_varying = ~1,
      time_varying_type = "rw0",
      priors = sdmTMB::sdmTMBpriors(sigma_V = sdmTMB::gamma_cv(0.3, 0.5)),
      extra_time = all_yrs,
      share_range = TRUE,
      anisotropy = TRUE,
      silent = FALSE
    ))
    if(!is.null(fits[[paste("IID RF, RW year",datnames[i],sep=' - ')]])) {
      saveRDS(fits[[paste("IID RF, RW year",datnames[i],sep=' - ')]], file.path(out_dir1, paste0(.sp,'-m6-d',i,'.rds')))
    }

    # fits[[paste("Spatial only, RW year",datnames[i],sep=' - ')]] <- safe_fit(update(
    #   fits[[base_model]],
    #   formula. = . ~ factor(gear),
    #   spatiotemporal = "off",
    #   time_varying = ~1,
    #   time_varying_type = "rw0",
    #   priors = sdmTMB::sdmTMBpriors(sigma_V = sdmTMB::gamma_cv(0.3, 0.5)),
    #   extra_time = all_yrs
    # ))
    fits[[paste("Spatial only, RW year",datnames[i],sep=' - ')]] <- safe_fit(sdmTMB(
      catch_weight ~ factor(gear),
      data = dat[[i]],
      mesh = mesh[[i]],
      offset = log(dat[[i]]$effort),
      family = delta_gamma(type = "poisson-link"),
      time = "year",
      spatial = "on",
      spatiotemporal = "off",
      time_varying = ~1,
      time_varying_type = "rw0",
      priors = sdmTMB::sdmTMBpriors(sigma_V = sdmTMB::gamma_cv(0.3, 0.5)),
      extra_time = all_yrs,
      share_range = TRUE,
      anisotropy = TRUE,
      silent = FALSE
    ))
    if(!is.null(fits[[paste("Spatial only, RW year",datnames[i],sep=' - ')]])) {
      saveRDS(fits[[paste("Spatial only, RW year",datnames[i],sep=' - ')]], file.path(out_dir1, paste0(.sp,'-m7-d',i,'.rds')))
    }

    # fits[[paste("Spatial only, RW year noprior",datnames[i],sep=' - ')]] <- safe_fit(update(
    #   fits[[base_model]],
    #   formula. = . ~ factor(gear),
    #   spatiotemporal = "off",
    #   time_varying = ~1,
    #   time_varying_type = "rw0",
    #   extra_time = all_yrs # added
    # ))
    fits[[paste("Spatial only, RW year noprior",datnames[i],sep=' - ')]] <- safe_fit(sdmTMB(
      catch_weight ~ factor(gear),
      data = dat[[i]],
      mesh = mesh[[i]],
      offset = log(dat[[i]]$effort),
      family = delta_gamma(type = "poisson-link"),
      time = "year",
      spatial = "on",
      spatiotemporal = "off",
      time_varying = ~1,
      time_varying_type = "rw0",
      extra_time = all_yrs,
      share_range = TRUE,
      anisotropy = TRUE,
      silent = FALSE
    ))
    if(!is.null(fits[[paste("Spatial only, RW year noprior",datnames[i],sep=' - ')]])) {
      saveRDS(fits[[paste("Spatial only, RW year noprior",datnames[i],sep=' - ')]], file.path(out_dir1, paste0(.sp,'-m8-d',i,'.rds')))
    }

    # fits[[paste("Spatial only, factor(year)",datnames[i],sep=' - ')]] <- safe_fit(update(
    #   fits[[base_model]],
    #   spatiotemporal = "off"
    # ))

    tmp[[i]] <- fits |> purrr::map(\(x) safe_index(x, grid[[i]])) |>
        dplyr::bind_rows(.id = "model")
  }
  indexes <- dplyr::bind_rows(tmp)
  indexes$species <- .sp
  indexes
}

species_slug <- function(x) {
  gsub("(^-|-$)", "", gsub("[^[:alnum:]]+", "-", tolower(x)))
}

fit_species <- function(.sp, .survey, .out_dir, .overwrite = FALSE) {
  out_file <- file.path(.out_dir, paste0(species_slug(.sp), ".rds"))
  err_file <- file.path(.out_dir, paste0(species_slug(.sp), "-error.txt"))

  if (!.overwrite && file.exists(out_file)) {
    out <- readRDS(out_file)
    if (nrow(out) > 0L) {
      message(.sp, " - reading existing result.")
      return(out)
    }
  }

  message(.sp, " - starting.")
  out <- tryCatch(
    do_fit(.sp, .survey = .survey, .out_dir=.out_dir),
    error = function(e) {
      msg <- conditionMessage(e)
      message(.sp, " - species fit errored: ", msg)
      writeLines(msg, err_file)
      dplyr::tibble()
    }
  )
  if (nrow(out) > 0L) {
    saveRDS(out, out_file)
  }
  out
}

read_species_results <- function(.spp, .out_dir) {
  out_files <- file.path(.out_dir, paste0(species_slug(.spp), ".rds"))
  out_files <- out_files[file.exists(out_files)]
  purrr::map_dfr(out_files, \(x) {
    out <- readRDS(x)
    if (nrow(out) > 0L) out else NULL
  })
}

dat0 <- readRDS(file.path('data-raw','garn_ruse_catch_21sp.rds'))
spp <- unique(dat0$commonname)
dontfit <- c("taskekrabbe", "strandkrabbe")
spp_to_fit <- spp[!spp %in% dontfit]
# # spp_to_fit <- "torsk"
# # spp_to_fit <- "berggylt"
# # spp_to_fit <- "bergnebb"
# # spp_to_fit <- "blåstål"
# spp_to_fit <- "grasgylt"
# spp_to_fit <- "grønngylt"
# spp_to_fit <- "kveite"
# # spp_to_fit <- "lange"
# spp_to_fit <- "lomre"
# spp_to_fit <- "lusuer"
# # spp_to_fit <- "lyr"
# spp_to_fit <- "pigghå"
# spp_to_fit <- "piggskate"
# # spp_to_fit <- "rødspette"
# spp_to_fit <- "sandflyndre"
# spp_to_fit <- "sei"
# # spp_to_fit <- "strandkrabbe"
# spp_to_fit <- "sypike"
# # spp_to_fit <- "taskekrabbe"
# spp_to_fit <- "vanlig ulke"
# spp_to_fit <- "ål"

# spp_to_fit <- c("berggylt",
#                 "bergnebb",
#                 "blåstål",
#                 "grasgylt",
#                 "grønngylt",
#                 "kveite",
#                 "lange",
#                 "lomre",
#                 "lusuer",
#                 "lyr",
#                 "pigghå",
#                 "piggskate",
#                 "rødspette",
#                 "sandflyndre",
#                 "sei",
#                 # "strandkrabbe",
#                 "sypike",
#                 # "taskekrabbe",
#                 "torsk",
#                 "vanlig ulke",
#                 "ål")

# spp_to_fit <- c("torsk","lyr")

RhpcBLASctl::blas_set_num_threads(1L)
RhpcBLASctl::omp_set_num_threads(1L)

dir.create("data-generated", showWarnings = FALSE)
species_out_dir <- file.path("data-generated", paste0(casestudy,"-species"))
dir.create(species_out_dir, showWarnings = FALSE)

workers <- min(length(spp_to_fit), max(1L, future::availableCores() - 1L))
future::plan(future::multisession, workers = workers)
# out <- purrr::map_dfr(spp_to_fit_syn, fit_species, .survey = "synoptic", .out_dir = species_out_dir)

out <- tryCatch(
  furrr::future_map_dfr(
    spp_to_fit,
    fit_species,
    .survey = casestudy,
    .out_dir = species_out_dir,
    .options = furrr::furrr_options(seed = TRUE, scheduling = 1)
  ),
  error = function(e) {
    read_species_results(spp_to_fit, species_out_dir)
  }
)
# out <- furrr::future_map_dfr(spp_to_fit_syn, do_fit, .survey = "hbll")
future::plan(future::sequential)

# not strictly odd/even split
# out$even <- out$year %% 2 == 0
# out$survey_group <- out$even
out <- out |>
  mutate(survey_group = case_when(
    year %in% c(2013,2016,2018,2020,2022) ~ "north",
    year %in% c(2015,2017,2019,2021,2023) ~ "south",
    year %in% c(2024,2025) ~ "north + south"
  ))
out <- out |>
  mutate(dataset = stringr::str_split_i(model, " - ", 2),
         model = stringr::str_split_i(model, " - ", 1))
out$survey_group[out$dataset == 'all years, no spatial overlap' & out$year == 2024] = 'north'
out$survey_group[out$dataset == 'all years, no spatial overlap' & out$year == 2025] = 'south'
out$species <- factor(out$species, levels=c('lyr','torsk'), labels=c('Pollack','Cod'))

saveRDS(out, file = file.path("data-generated",paste0(casestudy,"-indexes.rds")))
# out <- readRDS(file.path("data-generated",paste0(casestudy,"-indexes.rds")))

moving_window_acf <- function(x, window = 10L) {
  n <- length(x)
  if (n < window) {
    return(NA_real_)
  }
  acf_vals <- vapply(seq_len(n - window + 1L), function(i) {
    acf(x[i:(i + window - 1L)], plot = FALSE)$acf[2L]
  }, numeric(1L))
  min(acf_vals, na.rm = TRUE)
}

moving_window_amp <- function(x, window = 10L) {
  n <- length(x)
  if (n < window) {
    return(NA_real_)
  }
  amp_log <- function(y) mean(abs(diff(y, differences = 2))) / 4

  amp_val <- vapply(seq_len(n - window + 1L), function(i) {
    amp_log(x[i:(i + window - 1L)])
  }, numeric(1L))
  max(amp_val, na.rm = TRUE)
}

x <- out |>
  group_by(species, model, dataset) |>
  summarise(seesaw_index = acf(log_est, plot = FALSE)$acf[2])
x

out1 <- out |> filter(species %in% c('Pollack','Cod'))
out1 |>
  left_join(x) |>
  group_by(species, model, dataset) |>
  mutate(geomean = exp(mean(log(est))), est = est / geomean, lwr = lwr / geomean, upr = upr / geomean) |>
  ggplot(aes(year, log(est), ymin = log(lwr), ymax = log(upr))) +
  # geom_ribbon(fill = "grey90") +
  geom_line(aes(linetype = dataset)) +
  geom_linerange(aes(colour = survey_group, group = interaction(survey_group, dataset)), position=position_dodge(width = 0.3)) +
  geom_point(aes(colour = survey_group, shape=dataset, group = interaction(survey_group, dataset)), size = 2, position=position_dodge(width = 0.3)) +
  scale_colour_brewer(palette = "Dark2") +
  ggh4x::facet_grid2(model ~ species, scales = "free_y", independent = "y") +
  # facet_grid(model ~ species, scales = "free_y") +
  # facet_wrap(model ~ species, scales = "free_y")
  ylab("Biomass index") +
  xlab("Year") +
  labs(colour = "Survey coverage") +
  scale_x_continuous(breaks = seq(min(out$year), max(out$year), by = 2)) +
  ggsidekick::theme_sleek() +
  theme(strip.text.y = element_text(angle = 0),
        legend.position = "top", legend.direction='vertical', legend.box='horizontal')
ggsave(file.path("figs",paste0(casestudy,"-testing.pdf")), width = 8.5, height = 9.7, units='in')

out1 |>
  filter(dataset %in% c('all years, 2 years spatial overlap','all years, no spatial overlap')) |>
  group_by(species, model, dataset) |>
  # summarise(seesaw_index = abs(mean(log_est[which(even)]) - mean(log_est[which(!even)]))) |>
  summarise(seesaw_index = acf(log_est, plot = FALSE)$acf[2]) |>
  ungroup() |>
  # group_by(model) |>
  # mutate(mean_acf = mean(seesaw_index)) |>
  # ggplot(aes(forcats::fct_reorder(model, mean_acf), seesaw_index)) +
  ggplot(aes(x=model, y=seesaw_index, colour=species, shape=dataset)) +
  geom_point(aes(group=species), position=position_dodge(width = 0.4)) +
  coord_flip() +
  ylab("First-order index autocorrelation") +
  ggsidekick::theme_sleek() +
  theme(axis.title.y = element_blank(), panel.grid.major = element_line(colour = "grey90", linewidth = 0.3), panel.grid.minor = element_line(colour = "grey90", linewidth = 0.3))
ggsave(file.path("figs",paste0(casestudy,"-acf.pdf")), width = 8, height = 3.5)

# make_fig <- function(what, .ylab = "", include_all_data = FALSE, exclude = NULL) {
#   if (!include_all_data) {
#     out1 <- out |>
#       filter(!grepl("all data", model)) |>
#       mutate(all_data = FALSE)
#   } else {
#     out1 <- out |> mutate(all_data = grepl("all data", model))
#   }
#
#   if (!is.null(exclude)) {
#     out1 <- filter(out1, !species %in% exclude)
#   }
#
#   log_se_to_cv <- function(sigma_log) sqrt(exp(sigma_log^2) - 1)
#
#   out1 <- out1 |>
#     summarise(
#       min_acf = moving_window_acf(log_est, window = 10),
#       max_amp = moving_window_amp(log_est, window = 10),
#       mean_se = mean(log_se_to_cv(se)),
#       .by = c(species, model, all_data)
#     ) |>
#     group_by(model) |>
#     mutate(mean_acf = mean(min_acf)) |>
#     ungroup() |>
#     mutate(model = forcats::fct_reorder(model, mean_acf))
#
#   g <- out1 |>
#     ggplot(aes(model, {{ what }})) +
#     coord_flip() +
#     ylab(.ylab) +
#     ggsidekick::theme_sleek() +
#     theme(axis.title.y = element_blank(), panel.grid.major = element_line(colour = "grey90", linewidth = 0.3), panel.grid.minor = element_line(colour = "grey90", linewidth = 0.3))
#
#   if (include_all_data) {
#     blue <- RColorBrewer::brewer.pal(8, "Blues")[3]
#     orange <- RColorBrewer::brewer.pal(8, "Oranges")[3]
#     g <- g + geom_violin(data = out1, scale = "width", mapping = aes(colour = all_data, fill = all_data)) +
#       scale_colour_manual(values = c(blue, orange)) +
#       scale_fill_manual(values = c(blue, orange)) +
#       guides(colour = "none", fill = "none")
#   } else {
#     blue <- RColorBrewer::brewer.pal(8, "Blues")[3]
#     g <- g + geom_violin(scale = "width", colour = blue, fill = blue)
#   }
#
#   g + geom_point(position = position_jitter(width = 0), colour = "grey25")
#   return(out1)
# }
#
# g1 <- make_fig(min_acf, "Most negative\nlag-1 autocorrelation\nacross 10-year windows")
# g2 <- make_fig(max_amp, "Maximum biennial\noscillation amplitude\nacross 10-year windows")
#   # coord_flip(ylim = c(0, 1.05), expand = FALSE)
# patchwork::wrap_plots(g1, g2, axes = "collect")
# ggsave(file.path("figs",paste0(casestudy,"-acf-moving-window.pdf")), width = 6, height = 3.5)

# -------------------------------------------------------
# catch maps by species, sum by site
touse <- c("tidyverse","sdmTMB","sf","viridis","ggOceanMaps","ggrepel")
lapply(touse, require, character.only=TRUE, quietly=TRUE)

df <- readRDS(file.path('data-raw','garn_ruse_catch_21sp.rds'))
dfsite <- df |>  filter(!commonname %in% c("taskekrabbe", "strandkrabbe")) |>
  group_by(commonname, year, site) |>
  summarize(totkg = sum(biomass_kg),
            lon = mean(lon),
            lat = mean(lat)) %>% as.data.frame()
spnames <- unique(dfsite$commonname)

plotdir <- file.path("figs", paste0(casestudy,"-maps"))
if(!dir.exists(plotdir)) dir.create(plotdir)
for(i in 1:length(spnames)){
  sp = spnames[i]
  dfc <- filter(dfsite, commonname == sp)
  maxswept <- max(dfc$totkg)
  g <- basemap(limits = c(5, 12, 61.9, 67.1), land.col = "grey95") +
    scale_y_continuous(breaks = seq(62,67,by=1), name='', expand=c(0,0)) +
    scale_x_continuous(breaks = seq(5,12,by=3), name='', expand=c(0,0)) +
    ggspatial::geom_spatial_point(data = dfc %>% filter(totkg == 0), aes(lon, lat), shape ="x", show.legend=F, size=3, color='red')+
    ggspatial::geom_spatial_point(data = dfc %>% filter(totkg > 0), aes(lon, lat, size = totkg, color=totkg), alpha = 0.6)+
    facet_wrap(vars(year), nrow=2) +
    scale_size_continuous(limits = c(0,maxswept)) +
    scale_color_viridis(limits = c(0,maxswept)) +
    guides(color= guide_legend(), size=guide_legend()) +
    labs(x = "", y = "", size = "Biomass (kg)",color = "Biomass (kg)") +
    theme_bw() +
    theme(axis.text=element_text(size = 10), legend.position = 'top')
  ggsave(plot = g, filename = file.path(plotdir, paste0(sp,'-map.png')), width = 45, height = 25, units="cm", dpi = 300)
}



