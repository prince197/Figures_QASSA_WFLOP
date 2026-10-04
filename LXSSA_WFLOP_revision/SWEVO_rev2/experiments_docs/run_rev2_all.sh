#!/bin/sh
# Revision-2 experiments, sequential stages, 2 processes each (this machine: 2 cores)
cd "$(dirname "$0")"
date; python3 rev2_laplace.py laplace 0 1 --procs 2 > run_rev2_laplace.log 2>&1
date; python3 rev2_ga.py ga 0 1 --procs 2 > run_rev2_ga.log 2>&1
date; python3 rev2_spacing.py spacing 0 1 --procs 2 > run_rev2_spacing.log 2>&1
date; python3 rev2_ga.py gahr 0 1 --procs 2 > run_rev2_gahr.log 2>&1
date; python3 rev2_ga.py gaiea 0 1 --procs 2 > run_rev2_gaiea.log 2>&1
date; python3 rev2_site.py lg16 0 1 --procs=2 > run_rev2_lg16.log 2>&1
date; python3 rev2_gradient.py grad 0 1 --procs=2 > run_rev2_grad.log 2>&1
date; python3 rev2_site.py lg16b 0 1 --procs=2 > run_rev2_lg16b.log 2>&1
date; echo ALLDONE
