#!/usr/bin/env bash
# Run the self-test scripts one after the other, e.g. overnight:
#   nohup bash analysis/run-self-tests.sh > /dev/null 2>&1 &
# A script that fails doesn't stop the rest. 11 is last because it's the
# longest; its fits are cached as they finish, so it can be restarted.

cd "$(dirname "$0")/.." || exit 1
mkdir -p logs

for s in 08-qcs-hs-biennial-self-test 13-qcs-hs-omega-swap 11-nwfsc-split-self-test; do
  echo "$(date '+%F %T') start $s" >> logs/run-self-tests.log
  Rscript "analysis/$s.R" > "logs/$s.log" 2>&1
  echo "$(date '+%F %T') done $s (exit $?)" >> logs/run-self-tests.log
done
