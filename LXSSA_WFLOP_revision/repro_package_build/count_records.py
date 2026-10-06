#!/usr/bin/env python3
"""Count the per-run records of the study from the run files (the 36,190 / 11,640 accounting of Table tab:design).

Usage:  python3 count_records.py [--analysis-dir DIR]
Main study (36,190): the 11 methods of the benchmark (SSA, LX-SSA, DE, old-setting PSO from fresh_grid.csv; VNS from
fresh_vgrid.csv; SSA-VNS and LX-SSA-VNS from fresh_bgrid.csv; PSO, RS-VNS, MS-SLSQP, PSO-VNS from mpce_*) = 22,440,
RSD-VNS 2,040, split 1,080, PSO coefficient sweep 3,240, budget 3,600 and feasible initialization 1,620 (both without
the superseded Horns Rev rows, which hrfix replaces), Horns Rev 1 (hrfix) 970, IEA37 1,200.
Additional experiments (11,640): the rev2_<family>_s<i>of<k>.csv shards. Sensitivity studies: every analysis/rev3_*.csv that has the
run-record columns (Algorithm, Seed, Objective, Feasible) is counted and listed. Exit code 1 if a total differs.
"""
import argparse, glob, os, re, sys
import pandas as pd

ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
ap.add_argument("--analysis-dir", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "analysis"))
D = ap.parse_args().analysis_dir


def rd(pat):
    fs = sorted(glob.glob(os.path.join(D, pat)))
    if not fs:
        raise SystemExit(f"missing: {pat}")
    return pd.concat([pd.read_csv(f, usecols=lambda c: c not in ("Coordinates", "Curve")) for f in fs], ignore_index=True)


def nohr(df):
    return int((df.Dataset.astype(str) != "HR").sum())


c = {}
fg = rd("fresh_grid.csv")
c["benchmark: SSA, LX-SSA, DE, old-setting PSO (fresh_grid)"] = int(fg.Algorithm.isin(["SSA", "LXSSA", "DE", "PSO"]).sum())
c["benchmark: VNS (fresh_vgrid)"] = len(rd("fresh_vgrid.csv"))
c["benchmark: SSA-VNS, LX-SSA-VNS (fresh_bgrid)"] = len(rd("fresh_bgrid.csv"))
for e, lab in (("psoc", "PSO"), ("rsvns", "RS-VNS"), ("slsqp", "MS-SLSQP")):
    c[f"benchmark: {lab} (mpce_{e})"] = len(rd(f"mpce_{e}_s*of*.csv"))
c["benchmark: PSO-VNS (mpce_psobv, without superseded HR rows)"] = nohr(rd("mpce_psobv_s*of*.csv"))
c["RSD-VNS (mpce_rsdisc)"] = len(rd("mpce_rsdisc_s*of*.csv"))
c["split (mpce_psosplit + mpce_omega90)"] = len(rd("mpce_psosplit_s*of*.csv")) + len(rd("mpce_omega90_s*of*.csv"))
c["PSO coefficient sweep (mpce_csweep)"] = len(rd("mpce_csweep_s*of*.csv"))
c["budget 30,030 / 120,030 (mpce_b30k/b120k(+p), without HR)"] = nohr(pd.concat([rd(f"mpce_{e}_s*of*.csv") for e in ("b30k", "b30kp", "b120k", "b120kp")]))
c["feasible initialization (mpce_feas + feasp, without HR)"] = nohr(pd.concat([rd("mpce_feas_s*of*.csv"), rd("mpce_feasp_s*of*.csv")]))
c["Horns Rev 1 (mpce_hrfix, 24 shards)"] = len(rd("mpce_hrfix_s*of24.csv"))
c["IEA37 CS1 (mpce_iea16/iea36(+p))"] = sum(len(rd(f"mpce_{e}_s*of*.csv")) for e in ("iea16", "iea36", "iea16p", "iea36p"))
orig = sum(c.values())
for k, v in c.items():
    print(f"  {k:62s} {v:6d}")
print(f"main-study records: {orig} (expected 36190)")
r2 = {}
for f in sorted(glob.glob(os.path.join(D, "rev2_*_s*of*.csv"))):
    fam = re.sub(r"_s\d+of\d+\.csv$", "", os.path.basename(f))
    r2[fam] = r2.get(fam, 0) + len(pd.read_csv(f, usecols=["Algorithm"]))
print("additional-experiment records:", sum(r2.values()), "(expected 11640):", ", ".join(f"{k} {v}" for k, v in r2.items()))
r3 = {}
for f in sorted(glob.glob(os.path.join(D, "rev3_*.csv"))):
    try:
        cols = pd.read_csv(f, nrows=0).columns
    except Exception:
        continue
    if {"Algorithm", "Seed", "Objective", "Feasible"} <= set(cols):
        r3[os.path.basename(f)] = len(pd.read_csv(f, usecols=["Algorithm"]))
print("sensitivity-study files in run-record format (runs or re-evaluations):", sum(r3.values()), "rows" if r3 else "(none present)")
for k, v in r3.items():
    print(f"  {k:62s} {v:6d}")
sys.exit(0 if orig == 36190 and sum(r2.values()) == 11640 else 1)
