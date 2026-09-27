#!/bin/bash
set -euo pipefail
umask 077

# The verifier channel must remain root-only even when /logs is supplied by
# the evaluation harness as a writable mount.  Reassert ownership and mode at
# runtime before creating any reward or report files.
mkdir -p /logs/verifier
chown root:root /logs/verifier
chmod 0700 /logs/verifier

printf '0\n' > /logs/verifier/reward.txt
chmod 0600 /logs/verifier/reward.txt
pytest --ctrf /logs/verifier/ctrf.json /tests/test_verifier.py -rA
chmod 0600 /logs/verifier/ctrf.json 2>/dev/null || true
printf '1\n' > /logs/verifier/reward.txt
chmod 0600 /logs/verifier/reward.txt
