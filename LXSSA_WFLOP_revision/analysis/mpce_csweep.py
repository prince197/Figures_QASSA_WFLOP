"""Sensitivity study: PSO coefficient sweep and bound handling (experiment csweep of mpce_experiments.py).

Usage (from analysis/):  python3 mpce_csweep.py
Inputs : mpce_csweep_s<i>of<k>.csv (all shards; labels PSOW07C12 ... PSOW07C20, PSOOLD_VZERO, PSOOLD_VMAX),
         mpce_psoc_s0of1.csv (constriction reference PSOC, restricted to the 12 split cases; same run_grid, same
         seeds), fresh_grid.csv (stored runs of the old setting, label PSO; used only to verify that PSOW07C20
         reproduces them).
Outputs: mpce_summary_csweep.json, mpce_numbers_csweep.tex (\\NS... macros), mpce_supp_csweep.tex (tab:S-csweep,
         fig:S-csweep), ../figures_mpce/csweep.pdf (+ .png). Checks: mpce_check_csweep.py (S01 ...).

Design: stand-alone PSO, Np = 30, 200 iterations (6,030 evaluations), random starts, the 12 split cases x 30 seeds
(psosplit design). Sweep: w = 0.7, c1 = c2 = c in {1.2, 1.4, 1.6, 1.7, 1.8, 1.9, 2.0} (c = 2 = old setting);
order-2 boundary (Proposition S1, c1 = c2, w = 0.7): 2c < 24(1 - w^2)/(7 - 5w), i.e. c < 12(1 - w^2)/(7 - 5w) =
1.7486. Bound handling of the old setting: clip position, keep velocity (as implemented; = PSOW07C20), clip and
zero the clipped velocity components (PSOOLD_VZERO), velocity clamping |v| <= 0.2 (ub - lb) plus clipping
(PSOOLD_VMAX). Every setting starts from the same seeded initial swarm (seed-paired).

Statistics per setting: feasible final layouts (% of 360 runs; cases with >= 15 of 30 feasible), mean wake loss
(% of ideal) of the feasible runs averaged over cases, average rank per case among the 10 settings (7 sweep points,
2 bound variants, PSOC) with the pipeline's feasibility-aware rule (mpce_results.rank_rule: < 15 feasible runs ->
ranked below all qualifying settings), dynamics logged by mpce_experiments.PSOBH (final swarm spread / r, clipped
coordinates, feasible particle evaluations in iterations 101-200). Paired tests: McNemar (exact binomial on the
discordant case-seed pairs) for feasibility, Wilcoxon signed-rank on case means for loss (cases where both have
>= 15 feasible runs).
"""
import os, sys, glob, json, re
import numpy as np, pandas as pd
from scipy.stats import binomtest, wilcoxon

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from mpce_results import rank_rule  # noqa: E402  (the pipeline's feasibility-aware ranking rule)

W = 0.7
C_SWEEP = (1.2, 1.4, 1.6, 1.7, 1.8, 1.9, 2.0)
SWEEP = [f"PSOW07C{int(round(10 * c))}" for c in C_SWEEP]
VAR = ["PSOOLD_VZERO", "PSOOLD_VMAX"]
REF = "PSOC"
ALL = SWEEP + VAR + [REF]
C_OF = dict(zip(SWEEP, C_SWEEP)); C_OF.update(PSOOLD_VZERO=2.0, PSOOLD_VMAX=2.0, PSOC=1.49618)
W_OF = {a: W for a in SWEEP + VAR}; W_OF[REF] = 0.7298
BOUND_OF = {a: "clip" for a in SWEEP + [REF]}; BOUND_OF.update(PSOOLD_VZERO="vzero", PSOOLD_VMAX="vmax")
C_STAR = 12 * (1 - W ** 2) / (7 - 5 * W)                 # order-2 boundary for c1 = c2 at w = 0.7 (1.748571...)
SPLITCASES = [(ds, r, n) for ds in (1, 2) for r, n in ((500, 6), (750, 8), (1000, 10), (500, 10), (750, 12), (1000, 15))]
CASE = ["Dataset", "Radius", "Turbines"]
DYN = ["Spread", "VelMean", "ClipPct", "ClipLatePct", "FeasEvalLatePct"]
DSR = {1: "I", 2: "II"}
WORD = dict(zip("0123456789", ["Zero", "One", "Two", "Three", "Four", "Five", "Six", "Seven", "Eight", "Nine"]))
MAC = {a: "C" + "".join(WORD[d] for d in a[-2:]) for a in SWEEP}
MAC.update(PSOOLD_VZERO="VZero", PSOOLD_VMAX="VMax", PSOC="Constr")
LABEL = {a: f"$c={C_OF[a]:.1f}$" for a in SWEEP}
LABEL.update(PSOW07C20="$c=2.0$ (old; clip, keep $v$)", PSOOLD_VZERO="old, $v\\gets0$ on clip",
             PSOOLD_VMAX="old, $|v|\\le v_{\\max}$", PSOC="constriction")


def load():
    fs = sorted(glob.glob(os.path.join(HERE, "mpce_csweep_s*of*.csv")))
    if not fs:
        raise SystemExit("no mpce_csweep_s*of*.csv")
    D = pd.concat([pd.read_csv(f) for f in fs], ignore_index=True)
    pc = pd.read_csv(os.path.join(HERE, "mpce_psoc_s0of1.csv"))
    pc = pc[pd.MultiIndex.from_frame(pc[CASE]).isin(SPLITCASES)].copy()
    for k in DYN:
        pc[k] = np.nan
    D = pd.concat([D, pc], ignore_index=True)
    D["LossPct"] = 100 * D.WakeLoss / D.Ideal
    D["Feasible"] = D.Feasible.astype(bool)
    D["SpreadRel"] = D.Spread / D.Radius
    return D, fs


def verify(D):
    """PSOW07C20 must reproduce the stored old-setting runs (fresh_grid.csv, label PSO) bit for bit."""
    fg = pd.read_csv(os.path.join(HERE, "fresh_grid.csv"))
    fg = fg[(fg.Algorithm == "PSO") & pd.MultiIndex.from_frame(fg[CASE]).isin(SPLITCASES)]
    a = D[D.Algorithm == "PSOW07C20"].set_index(CASE + ["Seed"]).sort_index()
    b = fg.set_index(CASE + ["Seed"]).sort_index()
    common = a.index.intersection(b.index)
    same = ((a.loc[common, "WakeLoss"] == b.loc[common, "WakeLoss"]) & (a.loc[common, "Coordinates"] == b.loc[common, "Coordinates"])
            & (a.loc[common, "Curve"] == b.loc[common, "Curve"]) & (a.loc[common, "Calls"] == b.loc[common, "Calls"]))
    return dict(old_runs_compared=int(len(common)), old_runs_identical=int(same.sum()), old_runs_stored=int(len(b)))


def feas_same_across_datasets(D):
    """Whether every run is feasible in data set I iff the seed-paired run of the same geometry is feasible in data
    set II (the penalty dominates the objective while a layout is infeasible, so feasibility depends on the geometry)."""
    X = D.pivot_table(index=["Algorithm", "Radius", "Turbines", "Seed"], columns="Dataset", values="Feasible").dropna()
    return bool((X[1] == X[2]).all()), int(len(X))


def mcnemar(D, a, b):
    # if feasibility is identical in the two data sets, each geometry-seed pair is counted once (data set I only)
    if FEAS_SAME[0]:
        D = D[D.Dataset == 1]
    X = D[D.Algorithm.isin([a, b])].pivot_table(index=CASE + ["Seed"], columns="Algorithm", values="Feasible").dropna()
    fa, fb = X[a].astype(bool), X[b].astype(bool)
    n10, n01 = int((fa & ~fb).sum()), int((~fa & fb).sum())
    p = binomtest(n10, n10 + n01, 0.5).pvalue if n10 + n01 else 1.0
    return dict(a_only=n10, b_only=n01, p=float(p), pairs=int(len(X)))


def loss_test(S, a, b):
    Q = S[S.Qualified].pivot_table(index=CASE, columns="Algorithm", values="Loss")
    if a not in Q or b not in Q:
        return dict(cases=0)
    Q = Q[[a, b]].dropna()
    d = (Q[a] - Q[b]).values
    nz = d[np.abs(d) > 1e-12]
    p = float(wilcoxon(nz).pvalue) if len(nz) >= 1 else 1.0
    return dict(cases=int(len(Q)), mean_diff_pp=float(d.mean()) if len(d) else None, a_better=int((d < 0).sum()),
                b_better=int((d > 0).sum()), p=p)


FEAS_SAME = [False, 0]


def analyse(D):
    FEAS_SAME[:] = feas_same_across_datasets(D)
    g = D.groupby(CASE + ["Algorithm"])
    S = pd.DataFrame({"N": g.size(), "NFeas": g.Feasible.sum()})
    S["Loss"] = D[D.Feasible].groupby(CASE + ["Algorithm"]).LossPct.mean()
    S["Mean"] = D[D.Feasible].groupby(CASE + ["Algorithm"]).Objective.mean()
    S = S.reset_index()
    S["NFeas"] = S.NFeas.astype(int)
    S["Qualified"] = S.NFeas >= np.ceil(S.N / 2)
    S["Rank"] = np.nan
    for _, idx in S.groupby(CASE).groups.items():
        s = S.loc[idx]
        S.loc[idx, "Rank"] = rank_rule(s.Mean.values, s.NFeas.values, s.N.values)
    # cases in which every setting has >= 15 feasible runs (common basis for the loss comparison)
    Qall = S.pivot_table(index=CASE, columns="Algorithm", values="Qualified").reindex(columns=ALL).fillna(False).astype(bool)
    common = Qall[Qall.all(axis=1)].index
    Lq = S.pivot_table(index=CASE, columns="Algorithm", values="Loss")
    out = dict(feas_identical_across_datasets=FEAS_SAME[0], feas_pairs_compared_across_datasets=FEAS_SAME[1],
               mcnemar_pairs_basis="data set I only (180 geometry-seed pairs)" if FEAS_SAME[0] else "all 360 case-seed pairs",
               c_star=C_STAR, c_star_order1=2 * (1 + W), n_cases=len(SPLITCASES), seeds=int(D.Seed.nunique()),
               common_loss_cases=[list(map(int, c)) for c in common], n_common_loss_cases=int(len(common)))
    st = {}
    for a in ALL:
        d = D[D.Algorithm == a]; s = S[S.Algorithm == a]
        st[a] = dict(label=a, w=W_OF[a], c=C_OF[a], bound=BOUND_OF[a], runs=int(len(d)), feasible=int(d.Feasible.sum()),
                     feas_pct=float(100 * d.Feasible.mean()), cases_qualified=int(s.Qualified.sum()),
                     cases_all_feasible=int((s.NFeas == s.N).sum()), cases_no_feasible=int((s.NFeas == 0).sum()),
                     loss_common=float(Lq.loc[common, a].mean()) if len(common) else None,
                     loss_qualified=float(s[s.Qualified].Loss.mean()) if s.Qualified.any() else None,
                     loss_any_feasible=float(s.Loss.mean()) if s.Loss.notna().any() else None,
                     loss_all12=float(s.Loss.mean()) if bool(s.Qualified.all()) and len(s) == len(SPLITCASES) else None,
                     avg_rank=float(s.Rank.mean()),
                     spread_rel_median=float(d.SpreadRel.median()) if d.SpreadRel.notna().any() else None,
                     clip_pct_mean=float(d.ClipPct.mean()) if d.ClipPct.notna().any() else None,
                     clip_late_pct_mean=float(d.ClipLatePct.mean()) if d.ClipLatePct.notna().any() else None,
                     feas_eval_late_pct_mean=float(d.FeasEvalLatePct.mean()) if d.FeasEvalLatePct.notna().any() else None,
                     cpu_hours=float(d.Seconds.sum() / 3600),
                     per_case={f"{ds}-{r}-{n}": dict(nfeas=int(x.NFeas), loss=None if pd.isna(x.Loss) else float(x.Loss),
                                                      rank=float(x.Rank))
                               for ds, r, n, x in ((x.Dataset, x.Radius, x.Turbines, x) for x in s.itertuples())})
        # feasibility by farm size (N >= 10 vs N < 10)
        st[a]["feas_pct_n_ge10"] = float(100 * d[d.Turbines >= 10].Feasible.mean())
        st[a]["feas_pct_n_lt10"] = float(100 * d[d.Turbines < 10].Feasible.mean())
    out["settings"] = st
    # transition along the sweep
    f = [st[a]["feas_pct"] for a in SWEEP]
    steps = []
    for i in range(len(SWEEP) - 1):
        dc = C_SWEEP[i + 1] - C_SWEEP[i]
        steps.append(dict(c_from=C_SWEEP[i], c_to=C_SWEEP[i + 1], drop_pp=f[i] - f[i + 1], drop_pp_per_0p1=(f[i] - f[i + 1]) / dc * 0.1,
                          contains_c_star=bool(C_SWEEP[i] < C_STAR < C_SWEEP[i + 1])))
    out["steps"] = steps
    jmax = int(np.argmax([s["drop_pp_per_0p1"] for s in steps]))
    out["largest_drop"] = steps[jmax]
    below = [a for a in SWEEP if C_OF[a] < C_STAR]; above = [a for a in SWEEP if C_OF[a] > C_STAR]
    out["feas_just_below"] = st[below[-1]]["feas_pct"]; out["feas_just_above"] = st[above[0]]["feas_pct"]
    out["c_just_below"] = C_OF[below[-1]]; out["c_just_above"] = C_OF[above[0]]
    # 50 % crossing of the pooled feasibility (linear interpolation between sweep points)
    c50 = None
    for i in range(len(f) - 1):
        if f[i] >= 50 > f[i + 1]:
            c50 = C_SWEEP[i] + (f[i] - 50) / (f[i] - f[i + 1]) * (C_SWEEP[i + 1] - C_SWEEP[i])
            break
    out["c50"] = c50
    fl = [st[a]["feas_pct_n_ge10"] for a in SWEEP]
    c50l = None
    for i in range(len(fl) - 1):
        if fl[i] >= 50 > fl[i + 1]:
            c50l = C_SWEEP[i] + (fl[i] - 50) / (fl[i] - fl[i + 1]) * (C_SWEEP[i + 1] - C_SWEEP[i])
            break
    out["c50_n_ge10"] = c50l
    # step sizes of the other responses along the sweep (per 0.1 in c): is anything discontinuous at c*?
    resp = {}
    for key, fac in (("loss_common", 1.0), ("spread_rel_median", 100.0), ("feas_eval_late_pct_mean", 1.0),
                     ("clip_pct_mean", 1.0), ("feas_pct_n_ge10", 1.0), ("avg_rank", 1.0)):
        v = [st[a][key] * fac for a in SWEEP]
        resp[key] = [dict(c_from=C_SWEEP[i], c_to=C_SWEEP[i + 1], change_per_0p1=(v[i + 1] - v[i]) / (C_SWEEP[i + 1] - C_SWEEP[i]) * 0.1,
                          contains_c_star=bool(C_SWEEP[i] < C_STAR < C_SWEEP[i + 1])) for i in range(len(v) - 1)]
    out["response_steps"] = resp
    out["spread_ratio_c17_over_c14"] = st["PSOW07C17"]["spread_rel_median"] / st["PSOW07C14"]["spread_rel_median"]
    out["spread_ratio_c20_over_c17"] = st["PSOW07C20"]["spread_rel_median"] / st["PSOW07C17"]["spread_rel_median"]
    out["feas_monotone_nonincreasing"] = bool(all(f[i] >= f[i + 1] for i in range(len(f) - 1)))
    out["feas_all_100_below"] = bool(all(st[a]["feas_pct"] == 100 for a in below))
    # loss minimum along the sweep (common cases)
    lc = [st[a]["loss_common"] for a in SWEEP]
    out["loss_min_c"] = C_SWEEP[int(np.nanargmin(lc))] if len(common) else None
    ar = {a: st[a]["avg_rank"] for a in ALL}
    out["best_avg_rank"] = min(ar, key=ar.get)
    out["rank_order"] = sorted(ar, key=ar.get)
    # paired comparisons
    pairs = [("PSOOLD_VZERO", "PSOW07C20"), ("PSOOLD_VMAX", "PSOW07C20"), ("PSOOLD_VMAX", "PSOOLD_VZERO"),
             ("PSOOLD_VMAX", "PSOC"), ("PSOOLD_VZERO", "PSOC"), ("PSOW07C17", "PSOW07C18"), ("PSOW07C16", "PSOW07C17"),
             ("PSOW07C18", "PSOW07C19"), ("PSOC", "PSOW07C17"), ("PSOW07C14", "PSOC"), ("PSOOLD_VMAX", "PSOW07C17"),
             ("PSOOLD_VMAX", "PSOW07C18")]
    out["paired"] = {f"{a}_vs_{b}": dict(feas=mcnemar(D, a, b), loss=loss_test(S, a, b)) for a, b in pairs}
    out["cpu_hours_csweep"] = float(D[D.Algorithm != REF].Seconds.sum() / 3600)
    return S, out


def fmt(x, d=1):
    return "--" if x is None or (isinstance(x, float) and not np.isfinite(x)) else f"{x:.{d}f}"


def macros(out, ver):
    st = out["settings"]; M = {}
    M["NSBoundC"] = f"{out['c_star']:.4f}"
    M["NSBoundPhi"] = f"{2 * out['c_star']:.2f}"
    M["NSCases"] = str(out["n_cases"]); M["NSSeeds"] = str(out["seeds"])
    M["NSRuns"] = str(out["n_cases"] * out["seeds"])
    M["NSVmaxFrac"] = "0.2"
    for a in ALL:
        k = MAC[a]; s = st[a]
        M[f"NSFeas{k}"] = fmt(s["feas_pct"], 1)
        M[f"NSFeasN{k}"] = str(s["feasible"])
        M[f"NSQual{k}"] = str(s["cases_qualified"])
        M[f"NSLoss{k}"] = fmt(s["loss_common"], 2)
        if s["loss_all12"] is not None:
            M[f"NSLossAll{k}"] = fmt(s["loss_all12"], 2)
        M[f"NSRank{k}"] = fmt(s["avg_rank"], 2)
        M[f"NSFeasLarge{k}"] = fmt(s["feas_pct_n_ge10"], 1)
        if s["spread_rel_median"] is not None:
            M[f"NSSpread{k}"] = fmt(100 * s["spread_rel_median"], 1)
            M[f"NSClip{k}"] = fmt(s["clip_pct_mean"], 1)
            M[f"NSFeasEval{k}"] = fmt(s["feas_eval_late_pct_mean"], 1)
    M["NSCommonCases"] = str(out["n_common_loss_cases"])
    M["NSCBelow"] = f"{out['c_just_below']:.1f}"; M["NSCAbove"] = f"{out['c_just_above']:.1f}"
    M["NSFeasBelow"] = fmt(out["feas_just_below"], 1); M["NSFeasAbove"] = fmt(out["feas_just_above"], 1)
    ld = out["largest_drop"]
    M["NSDropFrom"] = f"{ld['c_from']:.1f}"; M["NSDropTo"] = f"{ld['c_to']:.1f}"; M["NSDropPP"] = fmt(ld["drop_pp"], 1)
    M["NSCFifty"] = fmt(out["c50"], 2)
    M["NSCFiftyLarge"] = fmt(out["c50_n_ge10"], 2)
    M["NSLossMinC"] = fmt(out["loss_min_c"], 1)
    for key, name in (("PSOOLD_VZERO_vs_PSOW07C20", "VZeroOld"), ("PSOOLD_VMAX_vs_PSOW07C20", "VMaxOld"),
                      ("PSOOLD_VMAX_vs_PSOC", "VMaxConstr"), ("PSOW07C17_vs_PSOW07C18", "SevenEight"),
                      ("PSOW07C16_vs_PSOW07C17", "SixSeven"), ("PSOOLD_VMAX_vs_PSOW07C17", "VMaxSeven"),
                      ("PSOOLD_VMAX_vs_PSOW07C18", "VMaxEight"), ("PSOW07C14_vs_PSOC", "FourConstr"),
                      ("PSOW07C18_vs_PSOW07C19", "EightNine")):
        p = out["paired"][key]
        M[f"NSMcN{name}A"] = str(p["feas"]["a_only"]); M[f"NSMcN{name}B"] = str(p["feas"]["b_only"])
        M[f"NSMcP{name}"] = fmt(p["feas"]["p"], 3) if p["feas"]["p"] >= 0.001 else "\\ensuremath{<}0.001"
        if p["loss"].get("cases"):
            M[f"NSLossDiff{name}"] = fmt(p["loss"]["mean_diff_pp"], 2)
            M[f"NSLossCases{name}"] = str(p["loss"]["cases"])
            M[f"NSLossBetter{name}"] = str(p["loss"]["a_better"])
            M[f"NSLossP{name}"] = fmt(p["loss"]["p"], 3) if p["loss"]["p"] >= 0.001 else "\\ensuremath{<}0.001"
    M["NSCPUHours"] = fmt(out["cpu_hours_csweep"], 1)
    M["NSOldIdentical"] = str(ver["old_runs_identical"]); M["NSOldCompared"] = str(ver["old_runs_compared"])
    L = ["% generated by mpce_csweep.py from mpce_summary_csweep.json -- do not edit by hand"]
    L += [f"\\newcommand{{\\{k}}}{{{v}}}" for k, v in M.items()]
    open(os.path.join(HERE, "mpce_numbers_csweep.tex"), "w").write("\n".join(L) + "\n")
    return M


def supp_tex(out, S):
    st = out["settings"]
    cs = out["c_star"]
    L = ["%% generated by mpce_csweep.py -- do not edit by hand",
         "%% Supplementary: PSO coefficient sweep and bound handling; requires booktabs,",
         "%% graphicx and mpce_numbers_csweep.tex", "",
         "\\begin{figure*}[!t]", "\\centering",
         "\\includegraphics[width=\\textwidth]{figures_mpce/csweep.pdf}",
         "\\caption{Stand-alone PSO with $w=0.7$ and $c_1=c_2=c$ on the 12 split cases $\\times$ 30 seeds (6{,}030 "
         "evaluations, random starts; line with circles), the old setting ($c=2$) with two other bound-handling rules "
         "(clipped velocity components set to zero, squares; velocity clamping $|v|\\le v_{\\max}=0.2\\,(u-l)$, diamonds; drawn slightly to the right of $c=2$), and the "
         "constriction setting (triangles; $w=0.7298$, $c=1.49618$; stored runs of the main comparison, spread not logged). The dashed vertical line is the order-2 stability boundary "
         f"$c^\\ast=12(1-w^2)/(7-5w)={cs:.4f}$ of Proposition~\\ref{{prop:S-pso-stability}} (stagnation, no bounds). "
         "Left: feasible final layouts (\\% of 360 runs). Middle: mean wake loss of the feasible runs, averaged over the "
         f"{out['n_common_loss_cases']} cases in which every setting has at least 15 feasible runs. Right: median final "
         "swarm spread (mean over particles of the RMS turbine distance to the global best), relative to the farm radius "
         "(log scale). All settings start from the same seeded initial swarm. Values: Table~\\ref{tab:S-csweep}.}",
         "\\label{fig:S-csweep}", "\\end{figure*}", "",
         "\\begin{table*}[!t]", "\\centering",
         "\\caption{PSO coefficient sweep and bound handling (12 split cases $\\times$ 30 seeds, 6{,}030 evaluations, random "
         "starts). Feasible: final layouts feasible (\\% of 360 runs); $N\\ge10$: the same for the 180 runs with "
         "$N=10, 12, 15$; qualified: cases with at least 15 of 30 feasible runs; $L$: mean wake loss (\\%) of the feasible "
         f"runs, averaged over the {out['n_common_loss_cases']} cases in which every setting is qualified; $L_{{12}}$: the same over "
         "all 12 cases (only for settings qualified in all 12); rank: average "
         "rank over the 12 cases among the ten settings (feasibility-aware rule of the main comparison: settings with fewer "
         "than 15 feasible runs rank below all qualified ones); spread: median final swarm spread relative to the radius "
         "(\\%); clipped: coordinates leaving the box per iteration (\\%, iterations 1--200); feasible evals: particle "
         "evaluations with a feasible layout in iterations 101--200 (\\%). The row $c=2.0$ reproduces the stored runs of "
         "the old setting bit for bit; the constriction row is the stored PSO data of the main comparison (dynamics not "
         "logged). Order-2: whether $(w,c)$ satisfies the order-2 condition of Proposition~\\ref{prop:S-pso-stability} "
         "($c<\\NSBoundC$ at $w=0.7$), which assumes no velocity clamping. For every setting and seed, the final layout "
         "is feasible in Data Set~I if and only if it is feasible in Data Set~II (same geometry; the penalty dominates "
         "while a layout is infeasible), so the feasibility columns rest on 6 geometries $\\times$ 30 seeds.}",
         "\\label{tab:S-csweep}", "\\scriptsize\\setlength{\\tabcolsep}{3.5pt}",
         "\\resizebox{\\ifdim\\width>\\textwidth\\textwidth\\else\\width\\fi}{!}{%",
         "\\begin{tabular}{cclcccccccccc}", "\\toprule",
         "$w$ & $c$ & Bound handling & Order-2 & Feasible (\\%) & $N\\ge10$ (\\%) & Qualified & $L$ (\\%) & $L_{12}$ (\\%) & Rank "
         "& Spread (\\%) & Clipped (\\%) & Feas.\\ evals (\\%) \\\\", "\\midrule"]
    bh = dict(clip="clip, keep $v$", vzero="clip, $v\\gets0$", vmax="$|v|\\le v_{\\max}$, clip")
    for a in SWEEP + VAR + [REF]:
        s = st[a]
        if a == "PSOOLD_VZERO":
            L.append("\\midrule")
        name = (f"0.7 & {s['c']:.1f}" if a in SWEEP else "0.7 & 2.0" if a in VAR else "0.7298 & 1.496")
        stable = "yes" if (2 * s["c"] < 24 * (1 - s["w"] ** 2) / (7 - 5 * s["w"])) else "no"
        sp = fmt(100 * s["spread_rel_median"], 1) if s["spread_rel_median"] is not None else "--"
        L.append(f"{name} & {bh[s['bound']]} & {stable} & {fmt(s['feas_pct'])} & {fmt(s['feas_pct_n_ge10'])} & "
                 f"{s['cases_qualified']} & {fmt(s['loss_common'], 2)} & {fmt(s['loss_all12'], 2)} & {fmt(s['avg_rank'], 2)} & {sp} & "
                 f"{fmt(s['clip_pct_mean'])} & {fmt(s['feas_eval_late_pct_mean'])} \\\\")
    L += ["\\bottomrule", "\\end{tabular}}", "\\end{table*}", ""]
    open(os.path.join(HERE, "mpce_supp_csweep.tex"), "w").write("\n".join(L))


def figure(out, D):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.ticker
    from mpce_results import FIG_STYLE, TW_MAIN, MUTED   # the figure style shared by all figures of the paper
    # Okabe-Ito slots (dataviz validator, --pairs all, light: normal-vision dE >= 15.6); identity also carried by
    # marker shape; the constriction setting (= PSO of the main comparison) has the PSO colour of all figures
    COL = dict(sweep="#0072b2", vzero="#e69f00", ref="#009e73", vmax="#cc79a7")
    plt.rcParams.update(FIG_STYLE)
    st = out["settings"]
    # printed in the main text at 0.8 of the text width (\swResNarrowGraphics), in the supplement at full width
    fig, axes = plt.subplots(1, 3, figsize=(0.8 * TW_MAIN, 2.15))
    spec = [("feas_pct", "feasible final layouts (%)", False, 1.0),
            ("loss_common", f"mean wake loss (%)\n(feasible runs, {out['n_common_loss_cases']} common cases)", False, 1.0),
            ("spread_rel_median", "final swarm spread / $r$ (%)", True, 100.0)]
    xs = np.array(C_SWEEP)
    handles = {}
    for ax, (key, yl, logy, sc) in zip(axes, spec):
        ax.axvline(out["c_star"], color=MUTED, ls=(0, (4, 2)), lw=0.9, zorder=1)
        y = np.array([st[a][key] if st[a][key] is not None else np.nan for a in SWEEP], float) * sc
        h, = ax.plot(xs, y, color=COL["sweep"], marker="o", ms=4, lw=1.4, mec="white", mew=0.7, zorder=3,
                     label="$w=0.7$, $c_1=c_2=c$; clip, keep $v$ ($c=2$: old setting)")
        handles["sweep"] = h
        for a, col, mk, dx, lab in (("PSOOLD_VZERO", COL["vzero"], "s", 0.035, "$c=2$, clipped $v$ set to 0"),
                                    ("PSOOLD_VMAX", COL["vmax"], "D", 0.07, "$c=2$, $|v|\\leq 0.2\\,(u-l)$"),
                                    ("PSOC", COL["ref"], "^", 0.0, "constriction ($w=0.7298$, $c=1.496$)")):
            v = st[a][key]
            if v is None or not np.isfinite(v):
                continue
            h, = ax.plot([C_OF[a] + dx], [v * sc], ls="none", marker=mk, ms=5.5, color=col, mec="white", mew=0.7,
                         zorder=4, label=lab)
            handles[a] = h
        if logy:
            ax.set_yscale("log")
            ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}"))
        ax.set_xlabel("$c$ ($c_1=c_2$)"); ax.set_ylabel(yl)
        ax.set_xlim(1.12, 2.13)
        ax.set_xticks([1.2, 1.4, 1.6, 1.8, 2.0])
        ax.set_xticklabels(["1.2", "1.4", "1.6", "1.8", "2.0"])
    axes[0].set_ylim(-3, 103)
    axes[0].text(out["c_star"] - 0.02, 8, f"order-2\nboundary\n$c^*={out['c_star']:.4f}$", ha="right", va="bottom",
                 fontsize=7.5, color=MUTED, linespacing=1.1)
    fig.tight_layout(w_pad=0.8, rect=(0, 0, 1, 0.80))
    order = ["sweep", "PSOOLD_VZERO", "PSOOLD_VMAX", "PSOC"]
    fig.legend([handles[k] for k in order if k in handles], [handles[k].get_label() for k in order if k in handles],
               loc="upper center", ncol=2, bbox_to_anchor=(0.5, 1.0), handletextpad=0.3, columnspacing=1.0)
    d = os.path.join(HERE, "..", "figures_mpce")
    fig.savefig(os.path.join(d, "csweep.pdf"), bbox_inches="tight")
    fig.savefig(os.path.join(d, "csweep.png"), dpi=170, bbox_inches="tight")
    plt.close(fig)


def main():
    D, fs = load()
    ver = verify(D)
    S, out = analyse(D)
    out["verification"] = ver
    out["inputs"] = [os.path.basename(f) for f in fs] + ["mpce_psoc_s0of1.csv (PSOC, 12 split cases)"]
    json.dump(out, open(os.path.join(HERE, "mpce_summary_csweep.json"), "w"), indent=1)
    M = macros(out, ver)
    supp_tex(out, S)
    figure(out, D)
    st = out["settings"]
    print(f"order-2 boundary c* = {out['c_star']:.5f}; verification {ver}")
    print(f"{'setting':14s} {'feas%':>6s} {'N>=10':>6s} {'qual':>4s} {'Lcommon':>8s} {'L12':>7s} {'rank':>5s} {'spread%':>8s} {'clip%':>6s} {'feasEv%':>7s}")
    for a in ALL:
        s = st[a]
        print(f"{a:14s} {s['feas_pct']:6.1f} {s['feas_pct_n_ge10']:6.1f} {s['cases_qualified']:4d} {fmt(s['loss_common'], 3):>8s} "
              f"{fmt(s['loss_all12'], 3):>7s} {s['avg_rank']:5.2f} {fmt(None if s['spread_rel_median'] is None else 100 * s['spread_rel_median']):>8s} "
              f"{fmt(s['clip_pct_mean']):>6s} {fmt(s['feas_eval_late_pct_mean']):>7s}")
    for s in out["steps"]:
        print(f"  {s['c_from']:.1f}->{s['c_to']:.1f}: drop {s['drop_pp']:.1f} pp ({s['drop_pp_per_0p1']:.1f} per 0.1){'  <- c*' if s['contains_c_star'] else ''}")
    for k, v in out["response_steps"].items():
        print(" ", k, " ".join(f"{x['c_from']:.1f}-{x['c_to']:.1f}:{x['change_per_0p1']:+.3f}" for x in v))
    print("c50", out["c50"], "c50 N>=10", out["c50_n_ge10"], "loss min c", out["loss_min_c"], "common cases", out["n_common_loss_cases"], "rank order", out["rank_order"])
    for k, v in out["paired"].items():
        print(" ", k, v)
    print("CPU h", out["cpu_hours_csweep"], "; macros", len(M))


if __name__ == "__main__":
    main()
