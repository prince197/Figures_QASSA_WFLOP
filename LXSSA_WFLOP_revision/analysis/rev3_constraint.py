"""Revision 3, item C1 (reviewer concern R5): alternative constraint handling.

Question: do the method conclusions survive a constraint handling that does not mix m^2 and m violations and does
not depend on mu = 1e10?  Design, seeds, budget, outcomes and comparison family: rev3_constraint_manifest.md
(written before any run).

Usage (from analysis/):
    python3 rev3_constraint.py validate [--procs 2]        bit-for-bit check of variant pen on seeds 1-2 against
                                                           the stored runs + order check of the Deb surrogate
    python3 rev3_constraint.py run VARIANT [--procs 2]     VARIANT in pen | deb | proj; seeds 31-60
Output: rev3_constraint_<VARIANT>.csv (columns of mpce_experiments.run_grid, coordinates with 17 significant digits
        via record_io.encode_coordinates, extra columns at the end), rev3_constraint_validate.csv / .log.

Variants (everything else -- optimizers, seeds, initial populations, budget accounting, Tracker -- unchanged):
  pen   (i)   the implemented handling: penalized wake loss F_p of Eq. (6) (authors_objective.make_objective,
              mu = 1e10, g_b in m^2, g_s in m) and box clipping to [-r, r] inside the operators (velocity kept in PSO).
  deb   (ii)  Deb's feasibility rules (Deb 2000, CMAME 186:311-338) on DIMENSIONLESS violations
              v = sum_i max(0, g_b_i) / r^2 + sum_{i<j} max(0, g_s_ij) / l_min, realized by the scalar surrogate
                  F_deb = L                 if the layout is feasible (study rule: l_ij >= l_min - 1e-6 m and
                                            sqrt(x_i^2 + y_i^2) <= r + 1e-6 m, identical to run_grid's check)
                        = BIG + v           otherwise,          BIG = 2^18 = 262,144 > every ideal value (max 2.107e5)
              Box clipping unchanged.
  proj  (iii) as deb, but every turbine outside the circle is projected radially onto it (p <- r p / |p|) inside the
              objective wrapper, IN PLACE: all five optimizers pass the stored array (or a row view of it) to the
              objective, so the repaired position replaces the candidate (Lamarckian repair; checked in validate).
              The operators' box clipping still acts first, so the repair is "clip to the box, then project".

Order equivalence of F_deb (stated in the manifest, checked in `validate`): every optimizer used here (PSO, RCGA, SSA,
DE, basic VNS and the PSO/SSA phase of the hybrids) uses objective values only through comparisons (<, argmin, argmax,
argsort); none uses them arithmetically (no fitness-proportional selection, no averaging, no differences of F).
L lies in [0, ideal] subset [0, BIG) and v > 0 for an infeasible layout, so feasible < infeasible always; two feasible
layouts compare by L exactly; two infeasible ones compare by BIG + v, a monotone (non-decreasing) rounding of v, so
the order is never reversed and two violations closer than ulp(BIG)/2 = 2^-35 ~ 2.9e-11 (dimensionless) become a tie
(the incumbent is kept, the tie rule of Deb's comparison) -- far below the 1e-6 m feasibility tolerance
(1e-6 / l_min = 3.2e-9; boundary 1e-6 m at r = 500 m -> g_b / r^2 ~ 4e-9).
"""
import os, sys, time, math, json
import numpy as np, pandas as pd
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import mpce_experiments as mx
import rev2_ga as rg
from authors_objective import make_objective as authors_make_objective
from wflop_model import farm_objective, min_spacing, R
from record_io import encode_coordinates

METHODS = ["PSOBV", "PSOC", "GA", "SSABV", "DE"]                 # PSO-VNS, PSO, GA (RCGA), SSA-VNS, DE
CASES = [(1, 500, 10), (1, 750, 6), (1, 1000, 15), (2, 500, 10), (2, 750, 6), (2, 1000, 15)]
SEEDS = range(31, 61)
BUDGET = 6030
SMIN = 8 * R
TOL = 1e-6
BIG = 2.0 ** 18
VARIANTS = ("pen", "deb", "proj")

_STATS = {}            # per-process, per-task evaluation statistics (one task at a time per process)
_POS = {}


def layout_stats(xy, radius, smin=SMIN):
    """(feasible flag of the study, v dimensionless, max boundary excess m, max spacing deficit m)."""
    n = len(xy)
    rr2 = (xy ** 2).sum(1)
    dd = np.sqrt(((xy[:, None] - xy[None]) ** 2).sum(-1))[np.triu_indices(n, 1)]
    gb = rr2 - radius ** 2
    gs = smin - dd
    v = np.maximum(0.0, gb).sum() / radius ** 2 + np.maximum(0.0, gs).sum() / smin
    feas = (n == 1 or min_spacing(xy) >= smin - TOL) and np.sqrt(rr2).max() <= radius + TOL   # = run_grid's rule
    return bool(feas), float(v), float(max(0.0, np.sqrt(rr2).max() - radius)), float(max(0.0, gs.max()) if n > 1 else 0.0)


def _note(feas):
    _STATS["n"] += 1
    if feas:
        _STATS["nf"] += 1
        if _STATS["first"] is None:
            _STATS["first"] = _STATS["n"]


def objective_factory(variant):
    """Replacement for mpce_experiments.make_objective(ds, radius) for the given variant."""
    def factory(dataset, radius, smin=SMIN, **kw):
        if variant == "pen":
            base = authors_make_objective(dataset, radius, smin, **kw)

            def f(x):
                val = base(x)                               # the stored objective, unchanged
                _note(layout_stats(x.reshape(-1, 2), radius, smin)[0])
                return val
            return f

        def f(x):
            if variant == "proj":                           # radial projection, written back into x
                p = x.reshape(-1, 2).copy()
                rn = np.sqrt((p ** 2).sum(1))
                out = rn > radius
                if out.any():
                    p[out] *= (radius / rn[out])[:, None]
                    x[:] = p.ravel()
                    _STATS["proj"] += int(out.sum())
            xy = x.reshape(-1, 2)
            obj, ideal = farm_objective(xy, dataset)
            if not ideal < BIG:
                raise AssertionError("BIG must exceed the ideal value")
            feas, v, _, _ = layout_stats(xy, radius, smin)
            _note(feas)
            if feas:
                return ideal - obj                          # wake loss L (as the first term of F_p)
            if not v > 0:
                raise AssertionError("infeasible layout with v = 0")
            return BIG + v
        return f
    return factory


def run_method_c(alg, seed, budget, f, wake, feasible, dim, lb, ub, radius, smin, bcons=None):
    """rev2_ga.run_method_ga (GA label; every other label -> mpce_experiments.run_method unchanged); keeps pos."""
    pos, tr = rg.run_method_ga(alg, seed, budget, f, wake, feasible, dim, lb, ub, radius, smin, bcons)
    _POS["pos"] = np.array(pos, dtype=float).copy()
    return pos, tr


def run_task(task):
    variant, alg, ds, rad, n, seed = task
    _STATS.clear(); _STATS.update(n=0, nf=0, first=None, proj=0); _POS.clear()
    orig_mo, orig_rm = mx.make_objective, mx.run_method
    mx.make_objective = objective_factory(variant)
    mx.run_method = run_method_c
    try:
        row = mx.run_grid((alg, ds, rad, n, seed, BUDGET, "random"))
    finally:
        mx.make_objective, mx.run_method = orig_mo, orig_rm
    xy = _POS["pos"].reshape(-1, 2)
    row["Coordinates3"] = row["Coordinates"]                  # run_grid's 3-decimal string (for validation only)
    row["Coordinates"] = encode_coordinates(xy)
    feas, v, bex, sdef = layout_stats(xy, rad)
    if feas != row["Feasible"]:
        raise AssertionError("feasibility rule mismatch")
    if _STATS["n"] != row["Calls"]:
        raise AssertionError("evaluation count mismatch")
    row.update(Variant=variant, Violation=v, MaxBoundExcess=bex, MaxSpacingDeficit=sdef,
               FeasEvalPct=100.0 * _STATS["nf"] / max(1, _STATS["n"]),
               FirstFeasCall=_STATS["first"] if _STATS["first"] is not None else -1,
               Projections=_STATS["proj"])
    return row


COLS = ["Algorithm", "Dataset", "Radius", "Turbines", "Seed", "Budget", "Init", "Calls", "Objective", "Ideal",
        "WakeLoss", "Feasible", "MinSpacing", "Seconds", "Coordinates", "Curve",
        "Variant", "Violation", "MaxBoundExcess", "MaxSpacingDeficit", "FeasEvalPct", "FirstFeasCall", "Projections"]


def tasks(variant, seeds=SEEDS):
    # heaviest first (N = 15, then 10, then 6)
    cases = sorted(CASES, key=lambda c: (-c[2], c[0]))
    return [(variant, a, *c, s) for c in cases for a in METHODS for s in seeds]


# ------------------------------------------------------------------------------------------------ validation
def stored_runs():
    def rd(fns, alg):
        d = pd.concat([pd.read_csv(os.path.join(HERE, f), float_precision="round_trip") for f in fns], ignore_index=True)
        d = d[d.Algorithm == alg].copy()
        d["Dataset"] = d.Dataset.astype(str)
        return d
    parts = [rd(["mpce_psobv_s0of2.csv", "mpce_psobv_s1of2.csv"], "PSOBV"), rd(["mpce_psoc_s0of1.csv"], "PSOC"),
             rd([f"rev2_ga_s{i}of10.csv" for i in range(10)], "GA"), rd(["fresh_bgrid.csv"], "SSABV"),
             rd(["fresh_grid.csv"], "DE")]
    d = pd.concat(parts, ignore_index=True)
    if "Budget" in d:
        d = d[(d.Budget.isna()) | (d.Budget == BUDGET)]
    if "Init" in d:
        d = d[(d.Init.isna()) | (d.Init == "random")]
    return d


def _fmt_like(stored, xy):
    dec = len(stored.split(";")[0].split()[0].split(".")[1]) if "." in stored.split(";")[0].split()[0] else 0
    if dec >= 15:
        return encode_coordinates(xy)
    return ";".join(f"{a:.{dec}f} {b:.{dec}f}" for a, b in xy)


def deb_order_check(n_samples=4000, seed=12345):
    """Pairwise check that F_deb realizes Deb's lexicographic order on sampled layouts of the six cases: counts of
    reversed pairs (must be 0) and of pairs tied by F_deb but not by the exact key."""
    rng = np.random.default_rng(seed)
    res = []
    for ds, rad, n in CASES:
        _STATS.clear(); _STATS.update(n=0, nf=0, first=None, proj=0)
        f = objective_factory("deb")(ds, rad)
        keys, F = [], []
        for k in range(n_samples):
            # mixture: random box layouts (mostly infeasible) and spread layouts on a ring (some feasible)
            if k % 2 == 0:
                x = rng.uniform(-rad, rad, 2 * n)
            else:
                ang = 2 * np.pi * np.arange(n) / n + rng.uniform(0, 2 * np.pi) + rng.normal(0, 0.003, n)
                rr = rad * rng.uniform(0.995, 1.0 + 2e-9, n)
                x = np.c_[rr * np.cos(ang), rr * np.sin(ang)].ravel()
            xy = x.reshape(-1, 2)
            feas, v, _, _ = layout_stats(xy, rad)
            obj, ideal = farm_objective(xy, ds)
            keys.append((0, ideal - obj) if feas else (1, v)); F.append(f(x))
        F = np.array(F); kf = np.array([k[0] for k in keys]); kv = np.array([k[1] for k in keys])
        # exact lexicographic comparison vs scalar comparison over all pairs
        lt_key = (kf[:, None] < kf[None]) | ((kf[:, None] == kf[None]) & (kv[:, None] < kv[None]))
        eq_key = (kf[:, None] == kf[None]) & (kv[:, None] == kv[None])
        lt_F = F[:, None] < F[None]; eq_F = F[:, None] == F[None]
        reversed_ = int((lt_key & (F[:, None] > F[None])).sum())
        new_ties = int((eq_F & ~eq_key).sum() // 2)
        mism = int(((lt_key != lt_F) & ~(eq_F & ~eq_key)).sum())
        res.append(dict(Dataset=ds, Radius=rad, Turbines=n, Samples=n_samples, Feasible=int((kf == 0).sum()),
                        Pairs=n_samples * (n_samples - 1) // 2, Reversed=reversed_, TiesIntroduced=new_ties,
                        OtherMismatch=mism, MinV=float(kv[kf == 1].min()), MaxL=float(kv[kf == 0].max()) if (kf == 0).any() else float("nan")))
    return res


def inplace_check():
    """Verify that each optimizer hands the stored array (or a row view) to the objective, so that the in-place
    projection of variant proj is written back: run each method briefly on one case with a wrapper that marks the
    evaluated array and checks that the returned position is projected (inside the disc)."""
    out = []
    for alg in METHODS:
        row = run_task(("proj", alg, 1, 500, 10, 1))
        xy = np.array([[float(u) for u in p.split()] for p in row["Coordinates"].split(";")])
        out.append(dict(Algorithm=alg, Projections=row["Projections"],
                        MaxRadiusMinusR=float(np.sqrt((xy ** 2).sum(1)).max() - 500.0)))
    return out


def validate(procs):
    t0 = time.time()
    S = stored_runs()
    tl = [("pen", a, *c, s) for c in CASES for a in METHODS for s in (1, 2)]
    with Pool(procs) as pool:
        rows = pool.map(run_task, tl, chunksize=1)
    out = []
    for r in rows:
        m = S[(S.Algorithm == r["Algorithm"]) & (S.Dataset == str(r["Dataset"])) & (S.Radius == r["Radius"])
              & (S.Turbines == r["Turbines"]) & (S.Seed == r["Seed"])]
        if len(m) != 1:
            out.append(dict(Algorithm=r["Algorithm"], Dataset=r["Dataset"], Radius=r["Radius"], Turbines=r["Turbines"],
                            Seed=r["Seed"], Found=len(m)))
            continue
        m = m.iloc[0]
        xy = np.array([[float(u) for u in p.split()] for p in r["Coordinates"].split(";")])
        stored_c = m.Coordinates
        out.append(dict(Algorithm=r["Algorithm"], Dataset=r["Dataset"], Radius=r["Radius"], Turbines=r["Turbines"],
                        Seed=r["Seed"], Found=1, ObjectiveEqual=bool(float(m.Objective) == r["Objective"]),
                        WakeLossEqual=bool(float(m.WakeLoss) == r["WakeLoss"]),
                        FeasibleEqual=bool(bool(m.Feasible) == r["Feasible"]), CallsEqual=bool(int(m.Calls) == r["Calls"]),
                        CoordinatesEqual=bool(_fmt_like(stored_c, xy) == stored_c),
                        StoredDecimals=len(stored_c.split(";")[0].split()[0].split(".")[1]),
                        CurveEqual=bool(str(m.Curve) == r["Curve"]),
                        Objective=repr(float(r["Objective"])), StoredObjective=repr(float(m.Objective))))
    V = pd.DataFrame(out)
    V.to_csv(os.path.join(HERE, "rev3_constraint_validate.csv"), index=False)
    D = deb_order_check()
    P = inplace_check()
    rep = dict(bitwise=dict(runs=len(V), all_found=bool((V.Found == 1).all()),
                            objective_equal=int(V.get("ObjectiveEqual", pd.Series(dtype=bool)).sum()),
                            wakeloss_equal=int(V.get("WakeLossEqual", pd.Series(dtype=bool)).sum()),
                            coordinates_equal=int(V.get("CoordinatesEqual", pd.Series(dtype=bool)).sum()),
                            feasible_equal=int(V.get("FeasibleEqual", pd.Series(dtype=bool)).sum()),
                            calls_equal=int(V.get("CallsEqual", pd.Series(dtype=bool)).sum()),
                            curve_equal=int(V.get("CurveEqual", pd.Series(dtype=bool)).sum()),
                            per_method=V.groupby("Algorithm")[["ObjectiveEqual", "CoordinatesEqual", "CurveEqual"]].sum()
                            .astype(int).to_dict(orient="index")),
               deb_order=D, inplace_projection=P, seconds=round(time.time() - t0, 1))
    with open(os.path.join(HERE, "rev3_constraint_validate.json"), "w") as fh:
        json.dump(rep, fh, indent=1, default=float)
    print(json.dumps(rep, indent=1, default=float), flush=True)


def run(variant, procs):
    tl = tasks(variant)
    fn = os.path.join(HERE, f"rev3_constraint_{variant}.csv")
    part = fn + ".partial"
    done = set()
    if os.path.exists(part):
        P = pd.read_csv(part, float_precision="round_trip")
        done = set(zip(P.Algorithm, P.Dataset, P.Radius, P.Turbines, P.Seed))
    todo = [t for t in tl if (t[1], t[2], t[3], t[4], t[5]) not in done]
    t0 = time.time()
    print(f"{variant}: {len(tl)} runs, {len(todo)} to do, procs {procs}", flush=True)
    with Pool(procs) as pool:
        for k, row in enumerate(pool.imap_unordered(run_task, todo, chunksize=1), 1):
            row.pop("Coordinates3", None)
            pd.DataFrame([row])[COLS].to_csv(part, mode="a", header=not os.path.exists(part), index=False)
            if k % 25 == 0 or k == len(todo):
                el = time.time() - t0
                print(f"{variant}: {k}/{len(todo)} done, {el:.0f} s elapsed, ETA {el / k * (len(todo) - k):.0f} s",
                      flush=True)
    P = pd.read_csv(part, float_precision="round_trip", dtype={"Coordinates": str, "Curve": str})
    if len(P) != len(tl) or P.duplicated(["Algorithm", "Dataset", "Radius", "Turbines", "Seed"]).any():
        raise SystemExit(f"{variant}: incomplete ({len(P)} of {len(tl)})")
    order = {t[1:]: i for i, t in enumerate(tl)}
    P["_o"] = [order[(a, d, r, n, s)] for a, d, r, n, s in zip(P.Algorithm, P.Dataset, P.Radius, P.Turbines, P.Seed)]
    P = P.sort_values("_o").drop(columns="_o")
    # 17 significant digits are kept by re-reading as strings: rewrite from the partial file verbatim per row
    P.to_csv(fn, index=False, float_format=None)
    os.remove(part)
    print(variant, fn, len(P), "runs; calls", P.Calls.min(), P.Calls.max(), flush=True)


if __name__ == "__main__":
    argv = sys.argv[1:]
    procs = 2
    args = []
    i = 0
    while i < len(argv):
        a = argv[i]
        if a.startswith("--procs="):
            procs = int(a.split("=")[1])
        elif a == "--procs":
            procs = int(argv[i + 1]); i += 1
        else:
            args.append(a)
        i += 1
    procs = min(procs, 2)                       # brief: at most 2 worker processes for this topic
    if args[0] == "validate":
        validate(procs)
    elif args[0] == "run":
        if args[1] not in VARIANTS:
            raise SystemExit(f"variant must be one of {VARIANTS}")
        run(args[1], procs)
    else:
        raise SystemExit(__doc__)
