# Install the pinned sdmTMB used for all analyses (1.2.0.9025).
# Run once per machine: Rscript install-sdmtmb.R
if (!requireNamespace("pak", quietly = TRUE)) install.packages("pak")
pak::pak("pbs-assess/sdmTMB@7a4a785")
