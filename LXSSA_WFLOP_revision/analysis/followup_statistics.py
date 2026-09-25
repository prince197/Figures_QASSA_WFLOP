"""Additional statistics for the DS1/DS2 follow-up cases: case-level Friedman test with Holm
post-hoc comparisons, bootstrap CIs, Vargha--Delaney A12, batch runtimes and best layouts."""
import numpy as np, pandas as pd
from scipy.stats import friedmanchisquare, rankdata, norm, mannwhitneyu

df = pd.read_csv("../selected_30_run_data.csv")
ALGS = ["LXSSA", "SSA", "PSO", "DE"]
cases = df.groupby(["Dataset", "Radius", "Turbines"])
means = cases.apply(lambda s: s.groupby("Algorithm").EnergyProduction.mean()).reindex(columns=ALGS)
print(means.round(2))
stat, p = friedmanchisquare(*[means[a] for a in ALGS])
ranks = np.vstack([rankdata(-row) for row in means.values])
avg = ranks.mean(0)
n, k = means.shape
ff = (n - 1) * stat / (n * (k - 1) - stat)   # Iman--Davenport
print("Friedman chi2=%.4f p=%.6g ; Iman-Davenport F=%.4f ; avg ranks" % (stat, p, ff), dict(zip(ALGS, avg)))
se = np.sqrt(k * (k + 1) / (6 * n))
post = []
for j, a in enumerate(ALGS[1:], 1):
    z = (avg[j] - avg[0]) / se
    post.append((a, z, 2 * norm.sf(abs(z))))
ps = np.array([q[2] for q in post]); order = np.argsort(ps); holm = np.empty(3)
m = 3; run = 0
for i, o in enumerate(order):
    run = max(run, min(1, (m - i) * ps[o])); holm[o] = run
for (a, z, praw), ph in zip(post, holm):
    print(f"  LXSSA vs {a}: z={z:.3f} p={praw:.4g} pHolm={ph:.4g}")
pd.DataFrame(dict(Algorithm=ALGS, AverageRank=avg)).assign(FriedmanChi2=stat, FriedmanP=p, ImanDavenportF=ff) \
    .to_csv("friedman_case_level.csv", index=False)

rng = np.random.default_rng(2026)
rows = []
for key, sub in cases:
    lx = sub[sub.Algorithm == "LXSSA"].EnergyProduction.values
    for b in ALGS[1:]:
        y = sub[sub.Algorithm == b].EnergyProduction.values
        boots = [rng.choice(lx, lx.size).mean() - rng.choice(y, y.size).mean() for _ in range(10000)]
        lo, hi = np.percentile(boots, [2.5, 97.5])
        u = mannwhitneyu(lx, y).statistic
        rows.append(dict(Dataset=key[0], Radius=key[1], Turbines=key[2], Baseline=b,
                         MeanDiff=lx.mean() - y.mean(), CI_low=lo, CI_high=hi, A12=u / (lx.size * y.size)))
ci = pd.DataFrame(rows); ci.to_csv("bootstrap_ci.csv", index=False); print(ci.round(3).to_string())

rt = cases.apply(lambda s: s.groupby("Algorithm").Runtime.first()).reindex(columns=ALGS)
print("batch seconds per run\n", rt.round(4))
rt.to_csv("batch_runtime_per_run.csv")

best = []
for key, sub in cases:
    for a in ALGS:
        r = sub[sub.Algorithm == a].sort_values("EnergyProduction").iloc[-1]
        best.append(dict(Dataset=key[0], Radius=key[1], Turbines=key[2], Algorithm=a, Seed=r.Seed,
                         Objective=r.EnergyProduction, WakeLoss=r.WakeLoss, Coordinates=r.Coordinates))
pd.DataFrame(best).to_csv("best_layouts_selected_cases.csv", index=False)
