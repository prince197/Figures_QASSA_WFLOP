"""Calibration: run the authors' SSA and LX-SSA (phi=0, chi=1, pop 30, 100 iterations) with the
reconstructed objective and compare with the recorded follow-up runs."""
import sys, numpy as np, pandas as pd
from multiprocessing import Pool
from scipy.stats import mannwhitneyu
from authors_optimizers import SSA, LXSSA
from authors_objective import make_objective
from wflop_model import farm_objective

CASES = [(1, 500, 4), (1, 750, 8), (1, 1000, 8), (2, 500, 4), (2, 750, 8), (2, 1000, 8)]


def job(a):
    alg, pen, ds, rad, n, seed = a
    f = make_objective(ds, rad, penalty=pen)
    opt = SSA(30, 100, seed=seed) if alg == "SSA" else LXSSA(30, 100, 0.0, 1.0, seed=seed)
    pos, score, _ = opt.optimize(f, 2 * n, -rad, rad)
    obj, ideal = farm_objective(pos.reshape(-1, 2), ds)
    return dict(Algorithm=alg, Penalty=pen, Dataset=ds, Radius=rad, Turbines=n, Seed=seed,
                Score=score, Objective=obj, Feasible=score < 1e9)


if __name__ == "__main__":
    algs = sys.argv[1].split(",")
    pens = sys.argv[2].split(",")
    jobs = [(a, p, *c, s) for a in algs for p in pens for c in CASES for s in range(1, 31)]
    with Pool(4) as pool:
        out = pd.DataFrame(pool.map(job, jobs, chunksize=4))
    out.to_csv(f"calibration_{'_'.join(algs)}_{'_'.join(pens)}.csv", index=False)
    rec = pd.read_csv("../selected_30_run_data.csv")
    for (alg, pen, ds, rad, n), s in out.groupby(["Algorithm", "Penalty", "Dataset", "Radius", "Turbines"]):
        r = rec[(rec.Algorithm == alg) & (rec.Dataset == ds) & (rec.Radius == rad) & (rec.Turbines == n)].EnergyProduction
        sf = s[s.Feasible].Objective
        print(f"{alg:6s} {pen:6s} {ds} {rad:5d} {n}  rerun {sf.mean():10.2f} ± {sf.std():7.2f} feas {s.Feasible.mean():.2f} | "
              f"recorded {r.mean():10.2f} ± {r.std():7.2f} | MWU p={mannwhitneyu(sf, r).pvalue:.3g}")
