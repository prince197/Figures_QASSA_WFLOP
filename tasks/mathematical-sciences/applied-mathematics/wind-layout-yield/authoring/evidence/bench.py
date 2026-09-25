"""Truth-informed benchmark: optimise directly on the clean hub-height truth (1 deg bins).

This is not a solver the agent could run (it reads the sealed truth); it only sets the
best-known achievable mean power used as the denominator of the layout-quality gate.
"""
import json
import sys

import numpy as np

import opt
import truth_eval as te

d = sys.argv[1]
starts = int(sys.argv[2]) if len(sys.argv) > 2 else 4
t = np.load(f"{d}/truth.npz")
c = te.load_curve(f"{d}/env/power_curve.csv")
cen, w, G = opt.build_bins(t["U"], t["theta"], c, 1)
ev = opt.Evaluator(cen, w, G)
best = None
for s in range(starts):
    xy, f = opt.optimize(24, 1000.0, ev, seed=1000 + s, iters=100000)
    P = te.farm_power(xy, t["U"], t["theta"], c)
    if best is None or P > best[0]:
        best = (P, xy)
np.savetxt(f"{d}/bench_layout.csv", best[1], delimiter=",", fmt="%.3f")
print(json.dumps(dict(bench=best[0])))
