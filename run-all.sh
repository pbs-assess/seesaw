#!/usr/bin/env bash
# Run every fit script, then every paper figure script. Run from anywhere,
# e.g. overnight on a server:
#   nohup ./run-all.sh > /dev/null 2>&1 &
#
#   ./run-all.sh               fits (each skips what is already cached), then figures
#   ./run-all.sh --figs-only   just the figure scripts (seconds to minutes)
#   ./run-all.sh --diagnostics also the non-paper index diagnostics (diagnostics/)
#
# A script that fails doesn't stop the rest. Progress is in logs/run-all.log
# and each script's output in logs/<script>.log. Package versions are written
# to session-info.txt (tracked).
#
# Fits are slow (hours to days in total; the self-tests most of all). Each
# fit script caches its output in data-generated/; delete a cache to refit.

cd "$(dirname "$0")" || exit 1
mkdir -p logs

figs_only=false
diagnostics=false
for arg in "$@"; do
  case "$arg" in
    --figs-only) figs_only=true ;;
    --diagnostics) diagnostics=true ;;
    *) echo "Unknown option: $arg" >&2; exit 1 ;;
  esac
done

fits=(
  fit-case-study-bc-trawl
  fit-case-study-goa-bc-trawl
  fit-case-study-bc-longline
  fit-case-study-norway
  fit-gradient-grid
  fit-gradient-example
  fit-prior-sensitivity
  fit-self-test-qcs-hs
  fit-self-test-nwfsc # longest; its fits are cached as they finish, so it can be restarted
)
figs=(
  fig-case-study-maps
  fig-case-study-seesaw
  fig-self-tests
  fig-gradient-mechanism
  fig-norway-full-coverage
  fig-sampling-designs
  fig-prior-sensitivity
  fig-gradient-grid
)

run() {
  echo "$(date '+%F %T') start $1" >> logs/run-all.log
  Rscript "analysis/$1.R" > "logs/$1.log" 2>&1
  echo "$(date '+%F %T') done $1 (exit $?)" >> logs/run-all.log
}

Rscript -e 'pkgs <- c("sdmTMB", "TMB", "fmesher", "surveyjoin", "gfplot", "dplyr", "ggplot2", "patchwork", "ggh4x", "furrr", "future", "sf"); for (p in pkgs) library(p, character.only = TRUE); print(sessionInfo())' > session-info.txt 2>&1

if [ "$figs_only" = false ]; then
  for s in "${fits[@]}"; do run "$s"; done
fi
for s in "${figs[@]}"; do run "$s"; done
if [ "$diagnostics" = true ]; then run fig-index-diagnostics; fi

# figs/ should hold only the paper figures (the list in copy-figs-to-overleaf.sh)
for f in figs/*; do
  grep -q "^  $(basename "$f" .pdf)\$" copy-figs-to-overleaf.sh || echo "Not a paper figure: $f" | tee -a logs/run-all.log
done
