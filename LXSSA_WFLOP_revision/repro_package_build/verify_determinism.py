#!/usr/bin/env python3
"""Determinism sample: rerun stored runs with the original drivers and compare them field by field with the records.

Usage:  python3 verify_determinism.py [--analysis-dir DIR] [--cases 1:500:6,2:1000:15] [--seeds 1] [--algs PSOBV,SSABV]

For every algorithm, case and seed, the run is recomputed with the driver that produced the stored record
(PSO-VNS: mpce_experiments.run_grid -> mpce_psobv_s*of2.csv; SSA-VNS: full_grid_experiments.run_grid ->
fresh_bgrid.csv), serialized with pandas exactly as the drivers do, and every CSV field except the wall-clock
`Seconds` is compared as text with the stored row (the objective is written as the shortest round-trip repr of
the double, so equal text means bit-identical doubles; coordinates and the convergence curve are compared as the
stored strings). Nothing is written except to stdout. Exit code 1 if any field differs.
"""
import argparse, csv, glob, io, os, sys, time

ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
ap.add_argument("--analysis-dir", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "analysis"))
ap.add_argument("--cases", default="1:500:6,2:1000:15")
ap.add_argument("--seeds", default="1")
ap.add_argument("--algs", default="PSOBV,SSABV")
a = ap.parse_args()
AD = os.path.abspath(a.analysis_dir)
sys.path.insert(0, AD)
os.chdir(AD)
import pandas as pd
import mpce_experiments as ME
import full_grid_experiments as FG

STORE = {"PSOBV": sorted(glob.glob(os.path.join(AD, "mpce_psobv_s*of2.csv"))),
         "SSABV": [os.path.join(AD, "fresh_bgrid.csv")]}


def stored_row(alg, ds, r, n, seed):
    for fn in STORE[alg]:
        with open(fn, newline="") as fh:
            for row in csv.DictReader(fh):
                if (row["Algorithm"], row["Dataset"], row["Radius"], row["Turbines"], row["Seed"]) == \
                        (alg, str(ds), str(r), str(n), str(seed)) and row.get("Init", "random") == "random" \
                        and row.get("Budget", "6030") == "6030":
                    return os.path.basename(fn), row
    raise KeyError((alg, ds, r, n, seed))


def rerun(alg, ds, r, n, seed):
    if alg == "SSABV":
        res = FG.run_grid((alg, ds, r, n, seed))
    else:
        res = ME.run_grid((alg, ds, r, n, seed, 6030, "random"))
    buf = io.StringIO()
    pd.DataFrame([res]).to_csv(buf, index=False)
    return next(csv.DictReader(io.StringIO(buf.getvalue())))


fails = 0
total = 0
for alg in a.algs.split(","):
    for c in a.cases.split(","):
        ds, r, n = map(int, c.split(":"))
        for seed in map(int, a.seeds.split(",")):
            src, old = stored_row(alg, ds, r, n, seed)
            t0 = time.time()
            new = rerun(alg, ds, r, n, seed)
            keys = [k for k in new if k != "Seconds"]
            diff = [k for k in keys if new[k] != old.get(k)]
            total += 1
            fails += bool(diff)
            print(f"{alg:6s} DS{ds} r={r} N={n:2d} seed {seed}: {'IDENTICAL' if not diff else 'DIFFERENT ' + ','.join(diff)}"
                  f" ({len(keys)} fields incl. Objective={new['Objective']}, Calls={new['Calls']}, {len(new['Curve'].split(';'))}"
                  f" curve points; stored in {src}; rerun {time.time() - t0:.1f} s, stored Seconds {float(old['Seconds']):.1f})",
                  flush=True)
print(f"{total - fails} of {total} reruns bit-identical to the stored records (all fields except Seconds)")
sys.exit(1 if fails else 0)
