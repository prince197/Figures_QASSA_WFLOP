#!/bin/bash
# usage: run.sh <outdir-or-empty>
docker run --rm -v "$1":/app/output:ro wly-tests bash -c "bash /tests/test.sh >/tmp/log 2>&1; tail -3 /tmp/log; cat /logs/verifier/reward.txt; cat /logs/verifier/metrics.json 2>/dev/null | tr -d '\n ' ; echo"
