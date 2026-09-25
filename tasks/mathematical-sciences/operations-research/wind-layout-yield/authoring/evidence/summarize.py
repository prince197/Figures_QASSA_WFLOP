"""Summarise calibration: shortfall vs best-known layout per seed, yield error, pass/fail at the gates."""
import json
import sys

import numpy as np

SHORT, YERR = float(sys.argv[2]), float(sys.argv[3])
rows = [json.loads(l) for l in open(sys.argv[1])]
routes = {}
for r in rows:
    best = max([r["bench"]] + [v["true_net"] for v in r.values() if isinstance(v, dict)])
    for k, v in r.items():
        if not isinstance(v, dict):
            continue
        s = 100 * (1 - v["true_net"] / best)
        y = v["pred_err_pct"]
        routes.setdefault(k, []).append((s, y, s <= SHORT and abs(y) <= YERR))

print(f"gates: shortfall <= {SHORT}%  |yield err| <= {YERR}%   seeds: {len(rows)}\n")
print("| route | n | shortfall % (median, worst) | yield error % (median, worst abs) | pass rate |")
print("|---|---|---|---|---|")
for k, v in routes.items():
    a = np.array(v)
    s, y, p = a[:, 0], a[:, 1], a[:, 2]
    wy = y[np.argmax(np.abs(y))]
    print(f"| {k} | {len(a)} | {np.median(s):.2f}, {s.max():.2f} | {np.median(y):+.2f}, {wy:+.2f} | {int(p.sum())}/{len(a)} |")
