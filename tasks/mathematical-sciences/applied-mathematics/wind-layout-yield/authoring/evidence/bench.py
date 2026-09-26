"""Truth-informed benchmark: optimise directly on the sealed 20-year hub-height truth.
Only sets the best-known achievable power (P_best); an agent cannot run this."""
import json, os, sys
import numpy as np
d = sys.argv[1]
os.environ["WLY_DATA"] = f"{d}/env"
import solve as S
import truth_eval as te
t = np.load(f"{d}/truth.npz")
curve = S.load_curve()
model = S.Model(t["U"], t["theta"], curve, width=1.0)
rng = np.random.default_rng(4242)
best = None
for s in range(int(sys.argv[2]) if len(sys.argv) > 2 else 4):
    xy, f = S.search(model, S.random_layout(rng), rng, 60000)
    if best is None or f > best[1]:
        best = (xy, f)
xy, f = S.search(model, best[0], rng, 60000, step0=60.0, jump=0.0)
xy, f = S.anneal(model, xy, rng, 200000)
P = te.farm_power(xy, t["U"], t["theta"], te.load_curve(f"{d}/env/power_curve.csv"))
np.savetxt(f"{d}/bench_layout.csv", xy, delimiter=",", fmt="%.3f")
print(json.dumps(dict(bench=P)))
