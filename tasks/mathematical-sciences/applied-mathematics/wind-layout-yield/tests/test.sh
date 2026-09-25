#!/bin/bash
# Always leaves a reward behind, including when pytest crashes.
set -uo pipefail
mkdir -p /logs/verifier
echo 0 > /logs/verifier/reward.txt
if pytest --ctrf /logs/verifier/ctrf.json /tests/test_yield.py -rA; then
  echo 1 > /logs/verifier/reward.txt
else
  echo 0 > /logs/verifier/reward.txt
fi
