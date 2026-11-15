#!/usr/bin/env bash
# Copy the paper figures to the Overleaf clone. Run from anywhere; the
# destination can be given as the first argument. Only copies: committing and
# pushing the Overleaf repo is left to you. The list below is also the list of
# what belongs in figs/ (run-all.sh warns about anything else).
set -euo pipefail
cd "$(dirname "$0")"
dest="${1:-$HOME/src/seesaw-overleaf/figs}"
figs=(
  case-study-maps
  case-study-seesaw
  self-tests
  gradient-mechanism
  norway-full-coverage
  sampling-designs
  prior-sensitivity-densities
  prior-sensitivity-indexes
  gradient-grid
)
for f in "${figs[@]}"; do
  cp -v "figs/$f.pdf" "$dest/"
done
