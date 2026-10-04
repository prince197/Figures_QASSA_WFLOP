#!/bin/sh
# Revision-2 experiments on a many-core machine. Copy the rev2_*.py files into LXSSA_WFLOP_revision/analysis of
# branch claude/clever-mendel-hc0j56 (repo prince197/Figures_QASSA_WFLOP) and run:  sh run_rev2_hpc.sh 64
P=${1:-64}
cd "$(dirname "$0")"
python3 rev2_laplace.py laplace 0 1 --procs $P
python3 rev2_ga.py ga 0 1 --procs $P
python3 rev2_spacing.py spacing 0 1 --procs $P
python3 rev2_ga.py gahr 0 1 --procs $P
python3 rev2_ga.py gaiea 0 1 --procs $P
python3 rev2_site.py lg16 0 1 --procs=$P
python3 rev2_site.py lg16b 0 1 --procs=$P
python3 rev2_gradient.py grad 0 1 --procs=$P
python3 rev2_analysis.py
