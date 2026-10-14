"""GA (rev2_ga.RCGA) on the Lillgrund 16-turbine block, completing the nine-method comparison of rev2_site.py.

Usage:  python3 rev4_galg.py EXP [SHARD NSHARDS] [--procs=P]
Output: rev4_<EXP>_s<SHARD>of<NSHARDS>.csv, the columns of rev2_lg16*.csv (Algorithm ... Seconds, Coordinates,
        Curve, Spacing).

Experiments (label GA, random initialization, seeds 1..30, N_p = 30)
  galg   30 seeds at  6,030 evaluations   (as rev2_site.py lg16)
  galgb  30 seeds at 30,030 evaluations   (as rev2_site.py lg16b)

Nothing is re-implemented: each run is rev2_site.run_lg (the unchanged Lillgrund runner used for the other nine
methods: same objective, penalty, l_min = 3 D, box, feasibility test, Tracker, record columns) executed through
rev2_ga._with_ga, which swaps mpce_experiments.run_method for rev2_ga.run_method_ga exactly as rev2_ga.py does for
the benchmark (ga), Horns Rev 1 (gahr) and IEA37 (gaiea) runs. The GA is rev2_ga.RCGA unchanged (SBX eta_c = 15,
p_c = 0.9, polynomial mutation eta_m = 20, p_m = 1/n, binary tournament, clipping, elitism); it seeds NumPy with the
run seed and draws the common init_pop first, so the seed-s initial population equals that of the other methods.
"""
import sys, time
import pandas as pd
from multiprocessing import Pool
import rev2_ga
import rev2_site

S30 = range(1, 31)


def tasks(exp):
    if exp == "galg":
        return [(rev2_site.run_lg, ("GA", 16, s, 6030, "random")) for s in S30]
    if exp == "galgb":
        return [(rev2_site.run_lg, ("GA", 16, s, 30030, "random")) for s in S30]
    raise ValueError(exp)


def _call(pair):
    fn, t = pair
    return rev2_ga._with_ga(fn, t)


if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    procs = 4
    for a in sys.argv[1:]:
        if a.startswith("--procs="):
            procs = int(a.split("=")[1])
        elif a.startswith("--procs"):
            raise SystemExit("use --procs=P")
    exp = args[0]
    shard, nsh = (int(args[1]), int(args[2])) if len(args) >= 3 else (0, 1)
    tl = tasks(exp)[shard::nsh]
    t0 = time.time()
    with Pool(procs) as pool:
        rows = pool.map(_call, tl, chunksize=1)
    out = pd.DataFrame(rows)
    fn = f"rev4_{exp}_s{shard}of{nsh}.csv"
    out.to_csv(fn, index=False)
    print(exp, fn, len(out), "runs in", round(time.time() - t0), "s; calls", out.Calls.min(), out.Calls.max(), flush=True)
