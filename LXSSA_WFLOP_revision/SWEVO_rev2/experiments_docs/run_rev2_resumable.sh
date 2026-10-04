#!/bin/bash
# Resumable driver: every experiment is split into small shards; a shard whose CSV exists is skipped.
cd "$(dirname "$0")"
run() { script=$1; exp=$2; n=$3; pfx=$4
  for ((i=0;i<n;i++)); do f=rev2_${exp}_s${i}of${n}.csv
    [ -s "$f" ] && continue
    echo "$(date -u +%H:%M:%S) start $exp $i/$n"
    python3 $script $exp $i $n ${pfx}2 >> run_rev2_${exp}.log 2>&1 || echo "FAILED $exp $i"
  done; }
run rev2_laplace.py laplace 12 "--procs "
run rev2_ga.py ga 10 "--procs "
run rev2_spacing.py spacing 40 "--procs "
run rev2_ga.py gahr 2 "--procs "
run rev2_ga.py gaiea 2 "--procs "
run rev2_site.py lg16 4 "--procs="
run rev2_gradient.py grad 40 "--procs="
run rev2_site.py lg16b 30 "--procs="
echo "$(date -u +%H:%M:%S) ALLDONE"
