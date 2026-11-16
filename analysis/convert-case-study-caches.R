# One-off: split the old per-species case-study caches into the per-species x
# model layout that run_case_study() now uses, so finished fits aren't redone.
# Old files are left in place (delete them afterwards). A model absent from an
# old file had failed, so it is cached as an empty result rather than refit.
#
# Rscript analysis/convert-case-study-caches.R [data-generated dir]

source(here::here("analysis/R/models.R"))

args <- commandArgs(trailingOnly = TRUE)
data_dir <- if (length(args) > 0L) args[[1]] else here::here("data-generated")

task_dirs <- list.dirs(data_dir, recursive = FALSE) |>
  keep(\(d) grepl("^case-study-", basename(d)))

for (d in task_dirs) {
  old_files <- list.files(d, pattern = "\\.rds$", full.names = TRUE) |>
    discard(\(f) grepl("__", basename(f)))
  for (f in old_files) {
    old <- readRDS(f)
    for (model in names(paper_models)) {
      new_file <- file.path(d, paste0(tools::file_path_sans_ext(basename(f)), "__", slug(model), ".rds"))
      if (file.exists(new_file)) next
      new <- if (is.null(old)) NULL else if (model %in% old$model) filter(old, model == !!model) else old[0, ]
      saveRDS(new, new_file)
    }
  }
  cat(basename(d), ": ", length(old_files), " species converted\n", sep = "")
}
