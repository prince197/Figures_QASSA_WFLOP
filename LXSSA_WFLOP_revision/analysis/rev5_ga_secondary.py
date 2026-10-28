"""GA (rev2_ga.RCGA) in the two secondary benchmark studies that so far lacked it: budget scaling /
feasibility-preserving initialization (mpce_experiments feas / b30k / b120k) and the 5D / 6D spacing
re-optimization (rev2_spacing.py spacing).

Usage:  python3 rev5_ga_secondary.py EXP [SHARD NSHARDS] [--procs=P]
        python3 rev5_ga_secondary.py verify          (seed-pairing and runner checks, no output file)
Output: rev5_ga_<EXP>_s<SHARD>of<NSHARDS>.csv with the columns of the records of the other methods of the study
        (mpce_<exp>_s*of*.csv: Algorithm ... Seconds, Coordinates, Curve; spacing: + Spacing, SMin).

Experiments (label GA, N_p = 30; the design of the other methods of each study, unchanged)
  feas     6 largest benchmark cases (mpce_experiments.LARGE), seeds 1..30, 6,030 calls, feasibility-preserving
           initialization (as mpce_experiments feas / feasx / feasp)                                  180 runs
  b30k     6 largest cases, seeds 1..30, 30,030 calls, random initialization (as b30k / b30kp)      180 runs
  b120k    6 largest cases, seeds 1..30, 120,030 calls, random initialization (as b120k / b120kp)   180 runs
  hr       the Horns Rev 1 16-turbine block of the same budget / initialization design (as hrfix): 120,030 calls
           random, seeds 1..10, and 6,030 calls feasible, seeds 1..30 (6,030 and 30,030 random exist:
           rev2_gahr_s*of2.csv)                                                                       40 runs
  spacing  the 12 split cases (mpce_experiments.SPLITCASES), minimum spacing 5D and 6D, seeds 1..30, 6,030
           calls, random initialization (as rev2_spacing.py spacing)                                 720 runs
GA at 6,030 calls with random initialization on the benchmark (the 6,030 column of the budget study and the 4D
column of the spacing study) is the stored rev2_ga_s*of10.csv (68 cases, seeds 1..30).

Nothing is re-implemented. feas / b30k / b120k / hr: each run is the unchanged mpce_experiments.run_grid /
run_hr executed through rev2_ga._with_ga, which swaps mpce_experiments.run_method for rev2_ga.run_method_ga (as
rev2_ga.py ga / gahr and rev4_galg.py). spacing: each run is the unchanged rev2_spacing.run_grid_sp, whose
dispatcher mpce_experiments.run_method_x is swapped for _run_method_x_ga (GA -> rev2_ga.run_method_ga with the
spacing-specific penalty / feasibility test passed in by run_grid_sp; every other label -> run_method_x
unchanged). The GA is rev2_ga.RCGA unchanged (binary tournament, SBX eta_c = 15 p_c = 0.9, polynomial mutation
eta_m = 20 p_m = 1/n, clipping, elitism); it seeds NumPy with the run seed and draws the common init_pop first
(init_hook.GEN = feasible_init generator for Init = "feasible"), so its seed-s initial population equals that of
the other methods (checked by `verify`).
"""
import sys, time
import numpy as np, pandas as pd
from multiprocessing import Pool
import rev2_ga
import mpce_experiments as mx
import rev2_spacing as rsp

S30 = range(1, 31)
S10 = range(1, 11)


# ------------------------------------------------------------------ spacing: GA through run_grid_sp
def _run_method_x_ga(alg, seed, budget, f, wake, feasible, dim, lb, ub, radius, smin, bcons=None):
    """Dispatcher for rev2_spacing.run_grid_sp: GA -> rev2_ga.run_method_ga, every other label -> the original
    mpce_experiments.run_method_x."""
    if alg == "GA":
        return rev2_ga.run_method_ga(alg, seed, budget, f, wake, feasible, dim, lb, ub, radius, smin, bcons)
    return _ORIG_RMX(alg, seed, budget, f, wake, feasible, dim, lb, ub, radius, smin, bcons)


_ORIG_RMX = mx.run_method_x


def run_grid_sp_ga(task):
    """rev2_spacing.run_grid_sp (unchanged) with _run_method_x_ga in place of mpce_experiments.run_method_x.
    run_grid_sp writes its Coordinates with record_io.encode_coordinates (17 significant digits) itself."""
    mx.run_method_x = _run_method_x_ga
    try:
        return rsp.run_grid_sp(task)
    finally:
        mx.run_method_x = _ORIG_RMX


def _ga_call(fn, task):
    return rev2_ga._with_ga(fn, task)


# ------------------------------------------------------------------ designs
def tasks(exp):
    """List of (wrapper, function, task); heaviest first."""
    if exp in ("feas", "b30k", "b120k"):
        budget = {"feas": 6030, "b30k": 30030, "b120k": 120030}[exp]
        init = "feasible" if exp == "feas" else "random"
        return [(_ga_call, mx.run_grid, ("GA", *c, s, budget, init)) for c in mx.LARGE[::-1] for s in S30]
    if exp == "hr":
        return ([(_ga_call, mx.run_hr, ("GA", 16, s, 120030, "random")) for s in S10]
                + [(_ga_call, mx.run_hr, ("GA", 16, s, 6030, "feasible")) for s in S30])
    if exp == "spacing":
        cases = sorted(mx.SPLITCASES, key=lambda c: (-c[2], -c[1], c[0]))          # largest N first (as rev2_spacing)
        return [(None, run_grid_sp_ga, ("GA", *c, s, 6030, "random", sp)) for c in cases for sp in ("6D", "5D")
                for s in S30]
    raise ValueError(exp)


def _call(item):
    wrap, fn, t = item
    return fn(t) if wrap is None else wrap(fn, t)


# ------------------------------------------------------------------ verification
class _Stop(Exception):
    pass


def _first_calls(kind, alg, task, ncalls=30):
    """First ncalls objective arguments of one run (the run is aborted after them)."""
    import authors_objective
    seen = []
    orig_mo = authors_objective.make_objective

    def mo(*a, **k):
        f = orig_mo(*a, **k)

        def g(x):
            seen.append(np.array(x, float).copy())
            if len(seen) >= ncalls:
                raise _Stop
            return f(x)
        return g
    mx.make_objective = mo; rsp.make_objective = mo
    try:
        t = (alg,) + tuple(task)
        if kind == "spacing":
            run_grid_sp_ga(t)
        else:
            rev2_ga._with_ga(mx.run_grid, t)
    except _Stop:
        pass
    finally:
        mx.make_objective = orig_mo; rsp.make_objective = orig_mo
        import init_hook
        init_hook.GEN = None
    return np.array(seen[:ncalls])


def verify():
    """(1) Seed pairing: the first 30 evaluations of GA equal those of PSO (PSOC) and PSO-VNS (PSOBV) for seeds 1, 2
    in each design; (2) runner check: one stored record of another method per design is reproduced bit for bit by
    the same code path (the wrappers pass every non-GA label through unchanged)."""
    ok = True
    designs = [("grid", (1, 500, 10), 30030, "random", None),
               ("grid", (2, 1000, 15), 120030, "random", None),
               ("grid", (1, 750, 12), 6030, "feasible", None),
               ("spacing", (2, 750, 8), 6030, "random", "5D"),
               ("spacing", (1, 1000, 10), 6030, "random", "6D")]
    for kind, c, b, init, sp in designs:
        for s in (1, 2):
            task = (*c, s, b, init) + ((sp,) if sp else ())
            ga = _first_calls(kind, "GA", task)
            for other in ("PSOC", "PSOBV"):
                o = _first_calls(kind, other, task)
                eq = ga.shape == o.shape == (30, 2 * c[2]) and np.array_equal(ga, o)
                ok &= eq
                print(f"pairing {kind:7s} case {c} B={b} {init:8s} {sp or '4D'} seed {s}: GA vs {other}: first 30 "
                      f"evaluations identical = {eq}", flush=True)
    # runner reproduction of stored records (cheap ones)
    chk = [("mpce_feas_s0of1.csv", mx.run_grid, ("PSOC", 1, 500, 10, 1, 6030, "feasible"), None),
           ("rev2_spacing", run_grid_sp_ga, ("PSOC", 2, 750, 8, 1, 6030, "random", "5D"), "5D")]
    import glob
    for src, fn, t, sp in chk:
        files = sorted(glob.glob(src if src.endswith(".csv") else f"{src}_s*of*.csv"))
        D = pd.concat([pd.read_csv(f, usecols=lambda c: c not in ("Coordinates", "Curve")) for f in files])
        m = ((D.Algorithm == t[0]) & (D.Dataset.astype(str) == str(t[1])) & (D.Radius == t[2]) & (D.Turbines == t[3])
             & (D.Seed == t[4]) & (D.Budget == t[5]) & (D.Init == t[6]))
        if sp:
            m &= D.Spacing == sp
        rec = D[m].iloc[0]
        row = fn(t) if fn is run_grid_sp_ga else rev2_ga._with_ga(fn, t)
        import init_hook
        init_hook.GEN = None
        same = (row["Objective"] == rec.Objective) and (row["Calls"] == rec.Calls) and (bool(row["Feasible"]) == bool(rec.Feasible))
        ok &= same
        print(f"runner check {src} {t}: Objective {row['Objective']!r} vs stored {rec.Objective!r}; Calls "
              f"{row['Calls']} vs {rec.Calls}; identical = {same}", flush=True)
    print("VERIFY", "PASS" if ok else "FAIL", flush=True)
    return ok


if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    procs = 4
    for a in sys.argv[1:]:
        if a.startswith("--procs="):
            procs = int(a.split("=")[1])
        elif a.startswith("--procs"):
            raise SystemExit("use --procs=P")
    exp = args[0]
    if exp == "verify":
        sys.exit(0 if verify() else 1)
    shard, nsh = (int(args[1]), int(args[2])) if len(args) >= 3 else (0, 1)
    tl = tasks(exp)[shard::nsh]
    print(exp, "shard", shard, "of", nsh, ":", len(tl), "runs,", procs, "processes", flush=True)
    t0 = time.time()
    with Pool(procs) as pool:
        rows = pool.map(_call, tl, chunksize=1)
    out = pd.DataFrame(rows)
    fn = f"rev5_ga_{exp}_s{shard}of{nsh}.csv"
    out.to_csv(fn, index=False)
    print(exp, fn, len(out), "runs in", round(time.time() - t0), "s; calls", out.Calls.min(), out.Calls.max(),
          "; sum of per-run Seconds", round(out.Seconds.sum()), flush=True)
