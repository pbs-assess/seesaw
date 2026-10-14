#!/usr/bin/env bash
# Run the case-study index fits and the self-tests one after the other, e.g.
# overnight:
#   nohup bash analysis/run-analyses.sh > /dev/null 2>&1 &
# A script that fails doesn't stop the rest. Progress is in logs/run-analyses.log
# and each script's output in logs/<script>.log.
#
# The case-study scripts (03, 04, 06) skip fitting if their fits file in
# data-generated/ exists; delete it to refit. 11 is last because it's the
# longest; its fits are cached as they finish, so it can be restarted.

cd "$(dirname "$0")/.." || exit 1
mkdir -p logs

for s in \
  03-fit-syn-models \
  04-fit-afsc-pbs-models \
  06-fit-hbll-models \
  08-qcs-hs-biennial-self-test \
  13-qcs-hs-omega-swap \
  11-nwfsc-split-self-test; do
  echo "$(date '+%F %T') start $s" >> logs/run-analyses.log
  Rscript "analysis/$s.R" > "logs/$s.log" 2>&1
  status=$?
  echo "$(date '+%F %T') done $s (exit $status)" >> logs/run-analyses.log
done
