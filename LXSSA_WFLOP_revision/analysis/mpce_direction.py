"""Direction-resolution robustness of the MPCE results (Phase 6, lead decision D16; answers R4 issues 1-3 and 5).

Usage (from the repository root or analysis/):  python3 analysis/mpce_direction.py [--procs 2]
Then:  python3 analysis/mpce_check_dir.py   (checks F01...)

Question: the benchmark wind rose has 24 direction bins of 15 deg and the Horns Rev model 72 bins of 5 deg. A Jensen
top-hat wake evaluated only at the bin centres leaves "blind" direction gaps between neighbouring bin centres, which
optimized layouts can exploit. Do the paper's conclusions survive a finer direction resolution?  All stored final
layouts are RE-EVALUATED (not re-optimized) with finer direction bins.

1. Benchmark (Kusiak-Song, 68 cases = data sets I/II x r = 500/750/1000 m x N = 2..10/12/15; 6,030 evaluations,
   random initialization, 30 seed-paired runs; the eight main methods and the component variants LX-SSA-VNS,
   RS-VNS and RSD-VNS). Each 15-deg bin j (centre theta_j, frequency w_j, Weibull scale psi_j and shape k_j) is split
   into S equal sub-bins of 15/S deg with centres theta_j - 7.5 + (m + 0.5) 15/S, m = 0..S-1. Every sub-bin inherits
   the SAME Weibull parameters (psi_j, k_j) and the frequency w_j / S (a piecewise-constant wind rose: the
   frequency of the sector is spread uniformly over its sub-bins; the speed distribution is unchanged). S = 1
   reproduces the benchmark exactly (check F01), S = 3 gives 5-deg bins and S = 15 gives 1-deg bins (primary).
   The wake-free (ideal) objective is identical for every S, so the wake loss changes only through the wakes.
   The Gaussian wake (Bastankhah-Porte-Agel, k* = 0.04, the "Gaussian wake" row of tab:robust-final) is also
   evaluated with S = 1 and S = 15.
   Statistics are those of mpce_results.py (imported read-only): the feasibility-aware ranking rule (a method is
   ranked by its mean feasible objective only if >= 15 of its 30 runs are feasible), Friedman / Holm post hoc,
   the case-mean Wilcoxon with imputation, and the TOST equivalence rule (90 % percentile bootstrap CI over the
   cases, 10,000 resamples, fixed seed, margin EQ_MARGIN = 0.05 pp). Feasibility does not change (same layouts).
2. Horns Rev 1, 16-turbine block (all 970 runs of experiment hrfix, 901 feasible, and the installed block).
   Paper model: 5-deg bins centred at 2.5, 7.5, ... (each 30-deg sector -> 6 bins with the sector's Weibull
   parameters and 1/6 of its frequency). Re-evaluated with 1-deg bins centred at 0.5, 1.5, ... (primary: each
   sector -> 30 bins with its Weibull A, k and 1/30 of its frequency), and, as phase checks, 1-deg bins at integer
   centres and 5-deg bins at 0, 5, ... (bins on a sector boundary go to the next sector clockwise, so every sector
   keeps exactly 30 / 6 bins and the frequencies sum to 1). Speed bins (1 m/s, 3..25 m/s) are unchanged.
3. PyWake check (pywake_check.py, run separately because PyWake is optional): its CSV output is read here.
4. IEA37 published layouts projected radially onto the boundary (iea37_projected.py, run from here).

Outputs (analysis/): mpce_summary_dir.json, mpce_numbers_dir.tex (\\NF... macros), mpce_supp_direction.tex
(supplementary tables, labels tab:F-...), mpce_direction_layouts.csv.gz (per-layout re-evaluations).
"""
import os, sys, json, math, time, argparse
import numpy as np, pandas as pd
from multiprocessing import Pool
from scipy.stats import kendalltau, wilcoxon

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import wflop_model as W
import hornsrev_model as H

FOCUS = "PSOBV"
MAIN8 = ["PSOBV", "PSOC", "SSABV", "SSA", "LXSSA", "DE", "BVNS", "SLSQP"]
ABL9 = ["PSOBV", "SSABV", "LXBV", "RSVNS", "RSDVNS", "BVNS", "PSOC", "SSA", "LXSSA"]
ALLM = MAIN8 + ["LXBV", "RSVNS", "RSDVNS"]
HR_ORDER = ["PSOBV", "SSABV", "LXBV", "RSVNS", "BVNS", "PSOC", "SSA", "LXSSA", "DE", "SLSQP"]
LAB = {"PSOBV": "PSO-VNS", "PSOC": "PSO", "SSABV": "SSA-VNS", "SSA": "SSA", "LXSSA": "LX-SSA", "DE": "DE",
       "BVNS": "VNS", "SLSQP": "MS-SLSQP", "LXBV": "LX-SSA-VNS", "RSVNS": "RS-VNS", "RSDVNS": "RSD-VNS"}
MAC = {"PSOBV": "PSOVNS", "PSOC": "PSO", "SSABV": "SSAVNS", "SSA": "SSA", "LXSSA": "LXSSA", "DE": "DE",
       "BVNS": "VNS", "SLSQP": "MSSLSQP", "LXBV": "LXSSAVNS", "RSVNS": "RSVNS", "RSDVNS": "RSDVNS"}
SUBS = (1, 3, 15)                       # sub-bins per 15-deg benchmark bin: 15, 5, 1 deg
RES = {1: "Fifteen", 3: "Five", 15: "One"}
RES_DEG = {1: "15", 3: "5", 15: "1"}
NONTRIVIAL_PP = 0.2                     # "non-trivial" cases: PSO-VNS mean wake loss >= 0.2 % (benchmark model)
ABL_CONTR = [("PSOBV", "PSOC"), ("SSABV", "RSDVNS"), ("LXBV", "RSDVNS"), ("PSOBV", "RSDVNS"), ("RSDVNS", "RSVNS"),
             ("SSABV", "RSVNS"), ("LXBV", "RSVNS"), ("PSOBV", "RSVNS"), ("PSOBV", "BVNS"), ("SSABV", "LXBV")]
# Horns Rev direction-bin settings: name -> (bin width deg, first centre deg)
HR_BINS = {"5deg_2.5": (5.0, 2.5), "1deg_0.5": (1.0, 0.5), "1deg_0": (1.0, 0.0), "5deg_0": (5.0, 0.0)}
HR_PRIMARY = "1deg_0.5"
HR_SET = [(6030, "random"), (6030, "feasible"), (30030, "random"), (120030, "random")]
HR_SET_NAME = {(6030, "random"): "SixKR", (6030, "feasible"): "SixKF", (30030, "random"): "ThirtyK",
               (120030, "random"): "OneTwentyK"}
HR_SET_LAB = {(6030, "random"): "6k R", (6030, "feasible"): "6k F", (30030, "random"): "30k R", (120030, "random"): "120k R"}


# ====================================================================== models
def bench_bins(sub):
    """Sub-bin centres (rad, mathematical angle as wflop_model.THETA) and parent 15-deg bin of each sub-bin."""
    j = np.repeat(np.arange(24), sub)
    m = np.tile(np.arange(sub), 24)
    deg = 15.0 * j + (m + 0.5) * 15.0 / sub          # sub = 1: 7.5, 22.5, ... = wflop_model.THETA
    return np.deg2rad(deg), j


_BB = {s: bench_bins(s) for s in SUBS}


def bench_objective(xy, ds, sub, wake="jensen"):
    """Benchmark objective (15 x expected farm power, kW) and wake-free value with `sub` sub-bins per 15-deg bin.
    Each sub-bin carries the parent bin's Weibull (psi_j, k_j) and frequency w_j / sub."""
    k, psi, w = W.DATASETS[int(ds)]
    th, j = _BB[sub]
    kk, pp, ww = k[j], psi[j], w[j] / sub
    delta = W.jensen_deficits(xy, th) if wake == "jensen" else W.gaussian_deficits(xy, th, kstar=0.04)
    psi_i = pp[:, None] * (1 - delta)
    ep = W.expected_power_linear(np.broadcast_to(kk[:, None], psi_i.shape), psi_i)
    obj = W.BIN_WIDTH * np.sum(ww[:, None] * ep)
    ideal = W.BIN_WIDTH * len(xy) * np.sum(ww * W.expected_power_linear(kk, pp))
    return float(obj), float(ideal)


def _bench_eval(args):
    ds, coords = args
    xy = W.parse_coords(coords)
    out = {}
    for s in SUBS:
        out[f"J{s}"], out[f"I{s}"] = bench_objective(xy, ds, s)
    for s in (1, 15):
        out[f"G{s}"], _ = bench_objective(xy, ds, s, wake="gauss")
    return out


def hr_bins(step, offset):
    """Horns Rev direction bins (meteorological 'from' directions, deg), their frequencies and speed-bin
    probabilities. Every bin carries the Weibull A, k of its 30-deg sector (floor((wd + 15) / 30) mod 12) and the
    sector frequency x step / 30; speed bins 3..25 m/s of width 1 m/s as hornsrev_model."""
    wd = np.arange(offset, 360.0, step)
    sec = (np.floor((wd + 15.0) / 30.0).astype(int)) % 12
    f = H.SEC_F[sec] * step / 30.0
    A, K = H.SEC_A[sec], H.SEC_K[sec]
    lo, hi = H.WS - 0.5, H.WS + 0.5
    pws = np.exp(-(lo[None] / A[:, None]) ** K[:, None]) - np.exp(-(hi[None] / A[:, None]) ** K[:, None])
    return dict(wd=wd, f=f, pws=pws, sec=sec, step=step, offset=offset)


def hr_aep(xy, b, with_wake=True, local_ct=False):
    """AEP (GWh/yr) of layout xy (N x 2, m) with bins b = hr_bins(...). Same equations as hornsrev_model.aep_gwh
    (Jensen top-hat, k = 0.04, hub-centre in-wake test, RSS). local_ct=True is a DIAGNOSTIC variant (used by
    pywake_check.py only): the thrust coefficient of each upstream turbine is taken at its own effective (waked)
    speed instead of the free-stream speed, turbines being processed in downstream order as in PyWake's
    PropagateDownwind."""
    xy = np.asarray(xy, float)
    n = len(xy)
    toward = np.deg2rad(270.0 - b["wd"])
    WS = H.WS
    if not with_wake:
        u = np.broadcast_to(WS[None, :, None], (len(toward), len(WS), n))
    else:
        dx = xy[:, 0][:, None] - xy[:, 0][None, :]
        dy = xy[:, 1][:, None] - xy[:, 1][None, :]
        c = np.cos(toward)[:, None, None]; s = np.sin(toward)[:, None, None]
        x = dx[None] * c + dy[None] * s
        lat = np.abs(-dx[None] * s + dy[None] * c)
        inw = (x > 0) & (lat < H.RR + H.KW * x)
        g2 = np.where(inw, (H.RR / (H.RR + H.KW * np.maximum(x, 0))) ** 4, 0.0)       # (n_wd, i, j) squared factor
        if not local_ct:
            G = np.sqrt(g2.sum(2))
            u = WS[None, :, None] * (1 - H.A_CT[None, :, None] * G[:, None, :])
        else:
            g = np.sqrt(g2)
            xd = xy[:, 0][None, :] * np.cos(toward)[:, None] + xy[:, 1][None, :] * np.sin(toward)[:, None]
            u = np.empty((len(toward), len(WS), n))
            for d in range(len(toward)):
                order = np.argsort(xd[d], kind="stable")
                a_up = np.zeros((len(WS), n))                     # 1 - sqrt(1 - CT(u_eff)) of processed turbines
                for i in order:
                    s2 = ((a_up * g[d, i][None, :]) ** 2).sum(1)
                    ui = WS * (1 - np.sqrt(s2))
                    u[d, :, i] = ui
                    a_up[:, i] = 1 - np.sqrt(1 - np.interp(ui, H.WS_TAB, H.CT_TAB))
    p = np.interp(u, H.WS_TAB, H.P_TAB)
    return float(np.einsum("dsn,ds,d->", p, b["pws"], b["f"]) * H.HOURS / 1e6)


# ====================================================================== formatting
def num(v, d=2, sign=False):
    if v is None or (isinstance(v, float) and not math.isfinite(v)):
        return "\\TBD{}"
    s = f"{v:+,.{d}f}" if sign else f"{v:,.{d}f}"
    s = s.replace(",", "{,}")
    return s.replace("-", "\\ensuremath{-}", 1) if s.startswith("-") else s


def pval(p):
    if p is None or not np.isfinite(p):
        return "--"
    if p < 1e-3:
        m, e = f"{p:.1e}".split("e")
        return f"\\ensuremath{{{m}\\times10^{{{int(e)}}}}}"
    if p >= 0.995:
        return "1.0"
    return f"{p:.2g}" if p < 0.1 else f"{p:.2f}"


def tp(p):
    if p is None or not np.isfinite(p):
        return "--"
    if p < 1e-3:
        m, e = f"{p:.1e}".split("e")
        return f"${m}\\times10^{{{int(e)}}}$"
    return f"{p:.3f}"


def tnum(v, d=2, sign=False):
    if v is None or not np.isfinite(v):
        return "--"
    s = f"{v:+,.{d}f}" if sign else f"{v:,.{d}f}"
    return s.replace(",", "{,}").replace("-", "$-$")


def ci_txt(ci, d=3):
    return "[" + ", ".join(num(v, d) for v in ci) + "]"


def ci_tab(ci, d=3):
    return "[" + ", ".join(tnum(v, d) for v in ci) + "]"


def listing(items):
    items = list(items)
    if not items:
        return "none"
    if len(items) == 1:
        return items[0]
    return ", ".join(items[:-1]) + " and " + items[-1]


def table(env, caption, label, spec, header, lines, sep="2.5pt", pos="!htb", resize=False, foot=None):
    body = "\n".join(lines)
    tab = f"\\begin{{tabular}}{{{spec}}}\n\\toprule\n{header} \\\\\n\\midrule\n{body}\n\\bottomrule\n\\end{{tabular}}"
    if resize:
        tab = "\\resizebox{\\textwidth}{!}{%\n" + tab + "}"
    f = f"\n\\par\\smallskip{{\\scriptsize {foot}}}" if foot else ""
    return (f"\\begin{{{env}}}[{pos}]\n\\centering\n\\caption{{{caption}}}\n\\label{{{label}}}\n"
            f"\\scriptsize\\setlength{{\\tabcolsep}}{{{sep}}}\n{tab}{f}\n\\end{{{env}}}\n")


def jclean(o):
    if isinstance(o, dict):
        return {str(k): jclean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [jclean(v) for v in o]
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating, float)):
        return None if not np.isfinite(o) else float(o)
    if isinstance(o, np.bool_):
        return bool(o)
    return o


# ====================================================================== benchmark
def bench_reevaluate(G, procs, log):
    F = G[G.Feasible].copy()
    todo = F.drop_duplicates(["Dataset", "Coordinates"])
    t = time.time()
    args = list(zip(todo.Dataset.astype(int), todo.Coordinates))
    if procs > 1:
        with Pool(procs) as pool:
            res = pool.map(_bench_eval, args, chunksize=50)
    else:
        res = [_bench_eval(a) for a in args]
    R = pd.DataFrame(res, index=todo.index)
    R["Dataset"], R["Coordinates"] = todo.Dataset.values, todo.Coordinates.values
    F = F.merge(R, on=["Dataset", "Coordinates"], how="left")
    log(f"  benchmark: {len(todo)} distinct feasible layouts re-evaluated in {time.time() - t:.0f} s")
    return F


def with_values(G, F, obj_col, ideal_col):
    """G (all runs) with Objective / LossPct replaced by a re-evaluation (feasible rows; infeasible rows unchanged,
    they never enter a mean)."""
    X = G.copy()
    key = ["Algorithm", "Dataset", "Radius", "Turbines", "Seed"]
    m = X[key].merge(F[key + [obj_col, ideal_col]], on=key, how="left")
    ok = X.Feasible.values
    X.loc[ok, "Objective"] = m.loc[ok, obj_col].values
    X.loc[ok, "Ideal"] = m.loc[ok, ideal_col].values
    X["WakeLoss"] = X.Ideal - X.Objective
    X["LossPct"] = 100 * X.WakeLoss / X.Ideal
    return X


def case_diff(S, a, b):
    """Per-case mean-loss difference L(a) - L(b) (pp), both methods qualified; Series indexed by case."""
    P = S.pivot_table(index=["Dataset", "Radius", "Turbines"], columns="Algorithm", values="Loss")
    Q = S.pivot_table(index=["Dataset", "Radius", "Turbines"], columns="Algorithm", values="Qualified").astype(bool)
    both = Q[a] & Q[b]
    return (P[a] - P[b])[both]


def within_case_tau(S0, S1, methods):
    """Mean Kendall tau between two per-case orderings (methods qualifying in both) and % cases with the same best
    method (convention of mpce_robustness.robust_table)."""
    A = S0.pivot_table(index=["Dataset", "Radius", "Turbines"], columns="Algorithm", values="Mean").reindex(columns=methods)
    B = S1.pivot_table(index=["Dataset", "Radius", "Turbines"], columns="Algorithm", values="Mean").reindex(columns=methods)
    Q = S0.pivot_table(index=["Dataset", "Radius", "Turbines"], columns="Algorithm", values="Qualified").reindex(columns=methods).astype(bool)
    taus, same = [], []
    for a, b, q in zip(A.values, B.values, Q.values):
        ok = q & np.isfinite(a) & np.isfinite(b)
        if ok.sum() >= 2:
            t = kendalltau(a[ok], b[ok])[0]
            if np.isfinite(t):
                taus.append(t)
        if ok.sum() >= 1:
            same.append(np.argmax(np.where(ok, a, -np.inf)) == np.argmax(np.where(ok, b, -np.inf)))
    return float(np.mean(taus)), float(100 * np.mean(same))


def bench_block(M, G, F, log):
    """All benchmark statistics per resolution."""
    out = dict(sub_bins=list(SUBS), subbin_rule="each 15-deg bin j -> S sub-bins of 15/S deg, centres theta_j - 7.5 + "
               "(m + 0.5) 15/S; every sub-bin keeps (psi_j, k_j) and gets frequency w_j / S", res={})
    Gm = G[G.Algorithm.isin(ALLM)]
    variants = {"rec": None}
    for s in SUBS:
        variants[f"J{s}"] = (f"J{s}", f"I{s}")
    for s in (1, 15):
        variants[f"G{s}"] = (f"G{s}", "I1")
    SS = {}
    for name, cols in variants.items():
        X = Gm if cols is None else with_values(Gm, F, *cols)
        SS[name] = (X, M.case_stats(X[X.Algorithm.isin(MAIN8)], MAIN8), M.case_stats(X[X.Algorithm.isin(ABL9)], ABL9))
    S_rec = SS["rec"][1]
    d_rec = case_diff(S_rec, FOCUS, "PSOC")
    cases = S_rec[["Dataset", "Radius", "Turbines"]].drop_duplicates()
    lossF = S_rec[S_rec.Algorithm == FOCUS].set_index(["Dataset", "Radius", "Turbines"]).Loss
    nontriv = lossF[lossF >= NONTRIVIAL_PP].index
    ds2n10 = [c for c in lossF.index if c[0] == "2" and c[2] >= 10]
    Fm = F[F.Algorithm.isin(ALLM)]
    for name, (X, S, SA) in SS.items():
        r = {}
        FR = M.friedman_block(M.rank_matrix(S, MAIN8), MAIN8, focus=FOCUS)
        FA = M.friedman_block(M.rank_matrix(SA, ABL9), ABL9, focus=FOCUS)
        CW = M.case_mean_wilcoxon(S, FOCUS, [a for a in MAIN8 if a != FOCUS])
        r["friedman"] = FR
        r["ablation_friedman"] = FA
        r["case_mean_wilcoxon"] = CW
        r["mean_loss_qualified"] = {a: float(S[(S.Algorithm == a) & S.Qualified].Loss.mean()) for a in MAIN8}
        r["mean_loss_qualified_abl"] = {a: float(SA[(SA.Algorithm == a) & SA.Qualified].Loss.mean()) for a in ABL9}
        r["mean_loss_by_ds"] = {a: {ds: float(S[(S.Algorithm == a) & S.Qualified & (S.Dataset == ds)].Loss.mean())
                                    for ds in ("1", "2")} for a in MAIN8}
        order = sorted(MAIN8, key=lambda a: FR["avg_rank"][a])
        r["rank_order"] = order
        # PSO-VNS - PSO
        d = case_diff(S, FOCUS, "PSOC")
        eq = M.tost(d.values)
        dn = d.reindex(nontriv).dropna()
        eqn = M.tost(dn.values)
        d2 = d.reindex(ds2n10).dropna()
        p2 = wilcoxon(d2.values).pvalue if (np.abs(d2.values) > 1e-9).any() else 1.0
        common = d_rec.index.intersection(d.index)
        a_, b_ = d_rec[common].values, d[common].values
        r["pso_pair"] = dict(
            all=eq, all_wilcoxon=CW["PSOC"], nontrivial=eqn, n_nontrivial=int(len(dn)),
            nontrivial_wilcoxon_p=float(wilcoxon(dn.values).pvalue) if (np.abs(dn.values) > 1e-9).any() else 1.0,
            ds2_n10=dict(n=int(len(d2)), mean=float(d2.mean()), focus_better=int((d2 < -1e-9).sum()),
                         other_better=int((d2 > 1e-9).sum()), wilcoxon_p=float(p2),
                         ci95=M.boot_ci(d2.values)),
            sign_changes_vs_recorded=int((np.sign(a_) * np.sign(b_) < 0).sum()),
            n_cases_nonzero_both=int(((np.abs(a_) > 1e-9) & (np.abs(b_) > 1e-9)).sum()),
            focus_better=int((d < -1e-9).sum()), other_better=int((d > 1e-9).sum()), ties=int((np.abs(d) <= 1e-9).sum()),
            abs_gt_margin=int((np.abs(d) > M.EQ_MARGIN).sum()), max_focus_better=float(-d.min()), max_other_better=float(d.max()),
            per_case={f"{c[0]}-{c[1]}-{c[2]}": float(v) for c, v in d.items()})
        # component-analysis contrasts
        CWA = {}
        for a, b in ABL_CONTR:
            w_ = M.case_mean_wilcoxon(SA, a, [b])[b]
            CWA[f"{a}-{b}"] = w_
        hp = M.holm([v["p"] for v in CWA.values()])
        for (k, v), h in zip(CWA.items(), hp):
            v["p_holm_family"] = float(h)
        r["ablation_contrasts"] = CWA
        if name != "rec":
            r["tau_within_case"], r["same_best_pct"] = within_case_tau(S_rec, S, MAIN8)
            r["tau_within_case_abl"], r["same_best_pct_abl"] = within_case_tau(SS["rec"][2], SA, ABL9)
            ar0 = SS["rec"][1]
            r0 = M.friedman_block(M.rank_matrix(ar0, MAIN8), MAIN8, focus=FOCUS)["avg_rank"]
            r["kendall_tau_avg_rank_order"] = float(kendalltau([r0[a] for a in MAIN8], [FR["avg_rank"][a] for a in MAIN8])[0])
            obj, idl = variants[name]
            base = "J1"
            rel = {ds: float(100 * (Fm[Fm.Dataset == ds][obj] / Fm[Fm.Dataset == ds][base] - 1).mean()) for ds in ("1", "2")}
            r["rel_change_objective_pct"] = rel
        out["res"][name] = r
        log(f"  [{name:4s}] best {LAB[FR['best_ranked']]}, ranks " +
            ", ".join(f"{LAB[a]} {FR['avg_rank'][a]:.2f}" for a in order[:4]) +
            f"; PSO-VNS-PSO {eq['mean_dloss_pp']:+.4f} CI90 [{eq['ci90_mean_dloss_pp'][0]:+.4f}, {eq['ci90_mean_dloss_pp'][1]:+.4f}]"
            f" eq={eq['equivalent']} p={CW['PSOC']['p']:.3g}; SLSQP gap {CW['SLSQP']['mean_dloss_pp']:+.3f};"
            f" DSII N>=10 {r['pso_pair']['ds2_n10']['mean']:+.3f} ({r['pso_pair']['ds2_n10']['focus_better']}/"
            f"{r['pso_pair']['ds2_n10']['n']}); sign changes {r['pso_pair']['sign_changes_vs_recorded']}")
    out["n_cases"] = int(len(cases))
    out["n_nontrivial"] = int(len(nontriv))
    out["n_layouts"] = int(len(Fm))
    out["n_layouts_by_method"] = {a: int((Fm.Algorithm == a).sum()) for a in ALLM}
    # reproduction of the benchmark: S = 1 re-evaluation vs recorded objective
    dev = (Fm.J1 - Fm.Objective).abs() / Fm.Objective.abs()
    out["reproduction"] = dict(median_rel_dev=float(dev.median()), max_rel_dev=float(dev.max()),
                               n_rel_dev_above_1e_4=int((dev > 1e-4).sum()),
                               ideal_equal_all_subs=bool(np.allclose(Fm.I1, Fm.I15, rtol=1e-12) and np.allclose(Fm.I1, Fm.I3, rtol=1e-12)),
                               ideal_matches_recorded=bool(np.allclose(Fm.I1, Fm.Ideal, rtol=1e-9)))
    # layout-level loss increase 15 -> 1 deg per method
    Fm = Fm.assign(L1=100 * (1 - Fm.J1 / Fm.I1), L15=100 * (1 - Fm.J15 / Fm.I15), L3=100 * (1 - Fm.J3 / Fm.I3))
    out["layout_loss_rise_pp"] = {a: float((Fm[Fm.Algorithm == a].L15 - Fm[Fm.Algorithm == a].L1).mean()) for a in ALLM}
    out["layout_loss_rise_rel_pct"] = {a: float(100 * ((Fm[Fm.Algorithm == a].L15 - Fm[Fm.Algorithm == a].L1).sum()
                                                        / Fm[Fm.Algorithm == a].L1.sum())) for a in ALLM}
    out["n_layouts_loss_decreases"] = int((Fm.L15 < Fm.L1 - 1e-9).sum())
    return out, Fm


# ====================================================================== Horns Rev
def hr_block(M, ALL, log):
    X = ALL[(ALL.Dataset == "HR") & (ALL.Turbines == 16)].copy().reset_index(drop=True)
    xyi, _ = H.site(16)
    B = {k: hr_bins(*v) for k, v in HR_BINS.items()}
    inst = {k: hr_aep(xyi, b) for k, b in B.items()}
    ideal = {k: hr_aep(xyi, b, with_wake=False) for k, b in B.items()}
    t = time.time()
    for k, b in B.items():
        X[k] = [hr_aep(W.parse_coords(c), b) for c in X.Coordinates]
    log(f"  Horns Rev: {len(X)} runs x {len(B)} bin settings in {time.time() - t:.0f} s; installed " +
        ", ".join(f"{k} {v:.3f}" for k, v in inst.items()))
    out = dict(bins={k: dict(width_deg=v[0], first_centre_deg=v[1], n_bins=int(len(B[k]["wd"])),
                             freq_sum=float(B[k]["f"].sum()),
                             bins_per_sector=sorted(set(np.bincount(B[k]["sec"], minlength=12).tolist())))
                     for k, v in HR_BINS.items()},
               primary=HR_PRIMARY, installed=inst, ideal=ideal,
               installed_loss_pct={k: 100 * (1 - inst[k] / ideal[k]) for k in B},
               n_runs=int(len(X)), n_feasible=int(X.Feasible.sum()))
    rec = X[X.Feasible]
    dev = (rec["5deg_2.5"] - rec.Objective).abs().max()
    out["reproduction_max_abs_dev_gwh"] = float(dev)
    out["installed_reproduces_hornsrev_model"] = bool(abs(inst["5deg_2.5"] - H.aep_gwh(xyi)) < 1e-9)
    # PyWake NOJ values of the same runs (optional; pywake_check.py --layouts)
    pwf = os.path.join(HERE, "pywake_check_hr16runs.csv")
    have_pw = os.path.exists(pwf)
    if have_pw:
        PW = pd.read_csv(pwf)
        X = X.merge(PW[["Algorithm", "Seed", "Budget", "Init", "PyWakeNOJ_1deg"]], on=["Algorithm", "Seed", "Budget", "Init"], how="left")
        meta = json.load(open(os.path.join(HERE, "pywake_check_hr16runs.json")))
        inst["pywake_noj_1deg"] = meta["installed_aep_gwh"]
    F = X[X.Feasible]
    cols = list(HR_BINS) + (["PyWakeNOJ_1deg"] if have_pw else [])
    ikey = {c: c for c in HR_BINS}
    ikey["PyWakeNOJ_1deg"] = "pywake_noj_1deg"
    meth = {}
    for (bud, ini) in HR_SET:
        for a in HR_ORDER:
            s = X[(X.Budget == bud) & (X.Init == ini) & (X.Algorithm == a)]
            if not len(s):
                continue
            f = s[s.Feasible]
            e = dict(runs=int(len(s)), feasible=int(len(f)))
            for c in cols:
                e[f"mean_{c}"] = float(f[c].mean()) if len(f) else None
                e[f"best_{c}"] = float(f[c].max()) if len(f) else None
                e[f"above_{c}"] = int((f[c] > inst[ikey[c]]).sum())
            e["mean_drop_gwh"] = float((f["5deg_2.5"] - f[HR_PRIMARY]).mean()) if len(f) else None
            meth[f"{bud}_{ini}_{a}"] = e
    out["methods"] = meth
    out["above_installed_total"] = {c: int((F[c] > inst[ikey[c]]).sum()) for c in cols}
    out["above_installed_by_method"] = {c: {a: int((F[F.Algorithm == a][c] > inst[ikey[c]]).sum()) for a in HR_ORDER} for c in cols}
    out["installed_drop_gwh"] = float(inst["5deg_2.5"] - inst[HR_PRIMARY])
    out["mean_drop_gwh_all"] = float((F["5deg_2.5"] - F[HR_PRIMARY]).mean())
    out["mean_drop_gwh_by_method"] = {a: float((F[F.Algorithm == a]["5deg_2.5"] - F[F.Algorithm == a][HR_PRIMARY]).mean()) for a in HR_ORDER}
    # PSO-VNS vs PSO (and every method) on jointly feasible seeds, per setting and model
    pairs = {}
    for (bud, ini) in HR_SET:
        s = X[(X.Budget == bud) & (X.Init == ini)]
        for b in HR_ORDER:
            if b == FOCUS:
                continue
            for c in cols:
                A_ = s[(s.Algorithm == FOCUS) & s.Feasible].set_index("Seed")[c]
                B_ = s[(s.Algorithm == b) & s.Feasible].set_index("Seed")[c]
                j = A_.index.intersection(B_.index)
                if len(j) < 2:
                    continue
                dd = (A_[j] - B_[j]).values
                p = float(wilcoxon(dd).pvalue) if (np.abs(dd) > 1e-12).any() else 1.0
                pairs[f"{bud}_{ini}_{b}_{c}"] = dict(n=int(len(j)), mean_diff_gwh=float(dd.mean()), p=p,
                                                      focus_better=int((dd > 0).sum()), other_better=int((dd < 0).sum()))
            # the paper's run-level test (infeasible runs ranked last; mpce_results.paired_vs), per model
    out["pairs_joint_feasible"] = pairs
    # paper's test (goodness with infeasible runs below every feasible run), Holm over the other 9 methods
    ptest = {}
    for (bud, ini) in HR_SET:
        s = X[(X.Budget == bud) & (X.Init == ini)]
        for c in cols:
            Y = s.copy()
            Y["Objective"] = Y[c]
            Y = Y[Y.Objective.notna() | ~Y.Feasible]
            rows = M.paired_vs(Y, FOCUS, [b for b in HR_ORDER if b != FOCUS and b in set(Y.Algorithm)])
            ptest[f"{bud}_{ini}_{c}"] = {r["Baseline"]: dict(p_holm=r["PHolm"], rb=r["RB"], outcome=r["Outcome"]) for r in rows}
    out["paper_test"] = ptest
    # highest mean AEP per setting and model (methods with >= half of the runs feasible, and all)
    hi = {}
    for (bud, ini) in HR_SET:
        for c in cols:
            ms = {a: meth[f"{bud}_{ini}_{a}"] for a in HR_ORDER if f"{bud}_{ini}_{a}" in meth}
            q = {a: e[f"mean_{c}"] for a, e in ms.items() if e["feasible"] >= math.ceil(e["runs"] / 2) and e[f"mean_{c}"] is not None}
            al = {a: e[f"mean_{c}"] for a, e in ms.items() if e[f"mean_{c}"] is not None}
            oq = sorted(q, key=lambda a: -q[a]); oa = sorted(al, key=lambda a: -al[a])
            hi[f"{bud}_{ini}_{c}"] = dict(best_qualified=oq[0], best_qualified_aep=q[oq[0]], second_qualified=oq[1],
                                          second_qualified_aep=q[oq[1]], best_any=oa[0], best_any_aep=al[oa[0]],
                                          order_qualified=oq)
    out["highest_mean"] = hi
    for (bud, ini) in HR_SET:
        h5, h1 = hi[f"{bud}_{ini}_5deg_2.5"], hi[f"{bud}_{ini}_{HR_PRIMARY}"]
        log(f"  HR {HR_SET_LAB[(bud, ini)]:7s}: highest mean 5deg {LAB[h5['best_qualified']]} {h5['best_qualified_aep']:.2f}, "
            f"1deg {LAB[h1['best_qualified']]} {h1['best_qualified_aep']:.2f} (2nd {LAB[h1['second_qualified']]} "
            f"{h1['second_qualified_aep']:.2f})")
    for bud in (6030, 30030):
        for c in cols:
            q = pairs.get(f"{bud}_random_PSOC_{c}")
            if q:
                log(f"  HR PSO-VNS - PSO {bud} {c}: {q['mean_diff_gwh']:+.3f} GWh (n={q['n']}, p={q['p']:.3g})")
    log(f"  HR runs above installed: " + ", ".join(f"{c} {v}" for c, v in out["above_installed_total"].items()))
    keep = ["Algorithm", "Seed", "Budget", "Init", "Feasible", "Objective"] + cols
    return out, X[keep]


# ====================================================================== macros / tables
class Mac:
    def __init__(self):
        self.lines = []

    def __call__(self, name, val, comment=""):
        assert name.isalpha(), name
        self.lines.append(f"\\newcommand{{\\{name}}}{{{val}}}" + (f"   % {comment}" if comment else ""))


ORD = {1: "first", 2: "second", 3: "third", 4: "fourth", 5: "fifth", 6: "sixth", 7: "seventh", 8: "eighth",
       9: "ninth", 10: "tenth", 11: "eleventh", 12: "twelfth", 13: "thirteenth", 14: "fourteenth"}


def write_macros(BM, HR, PW, IP, fn, stamp):
    m = Mac()
    R = BM["res"]
    rec, one, five = R["rec"], R["J15"], R["J3"]
    # ---- sub-bin rule, sizes
    m("NFSubBins", "15", "1-deg sub-bins per 15-deg benchmark bin")
    m("NFNLayouts", num(BM["n_layouts"], 0), "feasible final layouts re-evaluated (68 cases, 11 methods)")
    m("NFNCases", str(BM["n_cases"]))
    # ---- mean wake loss per method (qualified cases) at 15 (recorded), 5, 1 deg
    for a in MAIN8:
        m(f"NFLoss{MAC[a]}Fifteen", num(rec["mean_loss_qualified"][a], 3), "mean wake loss (%), 15-deg bins (recorded)")
        m(f"NFLoss{MAC[a]}Five", num(five["mean_loss_qualified"][a], 3))
        m(f"NFLoss{MAC[a]}One", num(one["mean_loss_qualified"][a], 3), "mean wake loss (%), 1-deg sub-bins")
        m(f"NFRank{MAC[a]}Fifteen", num(rec["friedman"]["avg_rank"][a], 2))
        m(f"NFRank{MAC[a]}One", num(one["friedman"]["avg_rank"][a], 2))
    for ds, nm in (("1", "DSI"), ("2", "DSII")):
        m(f"NFLossPSOVNS{nm}Fifteen", num(rec["mean_loss_by_ds"][FOCUS][ds], 2))
        m(f"NFLossPSOVNS{nm}One", num(one["mean_loss_by_ds"][FOCUS][ds], 2))
    rise = {a: one["mean_loss_qualified"][a] - rec["mean_loss_qualified"][a] for a in MAIN8}
    rel = {a: 100 * rise[a] / rec["mean_loss_qualified"][a] for a in MAIN8}
    m("NFLossRiseMin", num(min(rise.values()), 2), "pp, smallest increase of the mean wake loss over the 8 methods")
    m("NFLossRiseMax", num(max(rise.values()), 2))
    m("NFLossRiseMinBy", LAB[min(rise, key=rise.get)])
    m("NFLossRiseMaxBy", LAB[max(rise, key=rise.get)])
    m("NFLossRiseRelMin", num(min(rel.values()), 0), "% of the 15-deg wake loss")
    m("NFLossRiseRelMax", num(max(rel.values()), 0))
    m("NFLossRisePSOVNS", num(rise[FOCUS], 2))
    m("NFLossRiseRelPSOVNS", num(rel[FOCUS], 0))
    m("NFLossRiseMarginRatio", num(rise[FOCUS] / 0.05, 0), "increase of PSO-VNS / 0.05 pp margin")
    m("NFNLayoutsLossDecreases", str(BM["n_layouts_loss_decreases"]), "layouts whose wake loss is LOWER at 1 deg")
    # ---- ranking at 1 deg
    m("NFBestOne", LAB[one["friedman"]["best_ranked"]])
    m("NFBestFive", LAB[five["friedman"]["best_ranked"]])
    m("NFSecondOne", LAB[one["rank_order"][1]])
    m("NFRankOrderOne", listing(f"{LAB[a]} ({one['friedman']['avg_rank'][a]:.2f})" for a in one["rank_order"]))
    ns = [LAB[a] for a, p in one["friedman"]["p_holm_vs_focus"].items() if p >= 0.05]
    m("NFPostHocNotSigOne", listing(ns), "methods whose avg rank is not significantly different from PSO-VNS (Holm)")
    m("NFTauOne", num(one["tau_within_case"], 2), "mean within-case Kendall tau, 1 deg vs benchmark")
    m("NFSameBestOne", num(one["same_best_pct"], 0), "% cases with unchanged best method")
    m("NFTauRankOrderOne", num(one["kendall_tau_avg_rank_order"], 2), "Kendall tau of the 8 average ranks")
    m("NFDeltaDSIOne", num(one["rel_change_objective_pct"]["1"], 1, sign=True), "mean relative change of the objective (%)")
    m("NFDeltaDSIIOne", num(one["rel_change_objective_pct"]["2"], 1, sign=True))
    m("NFRobustRowOne", f"Direction bins 1$^\\circ$ & {tnum(one['rel_change_objective_pct']['1'], 1, True)} & "
      f"{tnum(one['rel_change_objective_pct']['2'], 1, True)} & {one['friedman']['avg_rank'][FOCUS]:.2f} & "
      f"{LAB[one['friedman']['best_ranked']]} & {one['tau_within_case']:.2f} & {one['same_best_pct']:.0f} \\\\",
      "row for tab:robust-final (Model & DS I & DS II & Rank & Best & tau & Same)")
    g1, g15 = R["G1"], R["G15"]
    m("NFGaussOneBest", LAB[g15["friedman"]["best_ranked"]], "best-ranked method, Gaussian wake with 1-deg bins")
    m("NFGaussOneRankPSOVNS", num(g15["friedman"]["avg_rank"][FOCUS], 2))
    m("NFGaussOneRankPSO", num(g15["friedman"]["avg_rank"]["PSOC"], 2))
    m("NFGaussFifteenRankPSOVNS", num(g1["friedman"]["avg_rank"][FOCUS], 2))
    m("NFGaussFifteenRankPSO", num(g1["friedman"]["avg_rank"]["PSOC"], 2))
    m("NFGaussOneRowOne", f"Gaussian, 1$^\\circ$ bins & {tnum(g15['rel_change_objective_pct']['1'], 1, True)} & "
      f"{tnum(g15['rel_change_objective_pct']['2'], 1, True)} & {g15['friedman']['avg_rank'][FOCUS]:.2f} & "
      f"{LAB[g15['friedman']['best_ranked']]} & {g15['tau_within_case']:.2f} & {g15['same_best_pct']:.0f} \\\\")
    # ---- PSO-VNS - PSO
    for key, r in (("Fifteen", rec), ("Five", five), ("One", one)):
        pp = r["pso_pair"]
        e = pp["all"]
        m(f"NFPairMean{key}", num(e["mean_dloss_pp"], 3), f"PSO-VNS - PSO case-mean wake loss (pp), {key}")
        m(f"NFPairCI{key}", ci_txt(e["ci90_mean_dloss_pp"]), "90% case-bootstrap CI (pp)")
        m(f"NFPairEq{key}", "yes" if e["equivalent"] else "no", "equivalent at +-0.05 pp (90% CI inside)")
        m(f"NFPairMinMargin{key}", num(e["min_margin_pp"], 3))
        m(f"NFPairP{key}", pval(pp["all_wilcoxon"]["p"]), "case-mean Wilcoxon p")
        m(f"NFPairPTost{key}", pval(e["p_tost"]))
        m(f"NFPairWins{key}", str(pp["focus_better"]), "cases with lower PSO-VNS mean loss")
        m(f"NFPairLosses{key}", str(pp["other_better"]))
        m(f"NFPairTies{key}", str(pp["ties"]))
        m(f"NFPairAbsGtMargin{key}", str(pp["abs_gt_margin"]), "cases with |difference| > 0.05 pp")
        m(f"NFPairMaxPSOVNSBetter{key}", num(pp["max_focus_better"], 2))
        m(f"NFPairMaxPSOBetter{key}", num(pp["max_other_better"], 2))
        en = pp["nontrivial"]
        m(f"NFNonTrivMean{key}", num(en["mean_dloss_pp"], 3), "non-trivial cases (PSO-VNS loss >= 0.2 %)")
        m(f"NFNonTrivCI{key}", ci_txt(en["ci90_mean_dloss_pp"]))
        m(f"NFNonTrivEq{key}", "yes" if en["equivalent"] else "no")
        m(f"NFNonTrivP{key}", pval(pp["nontrivial_wilcoxon_p"]))
        d2 = pp["ds2_n10"]
        m(f"NFDSIINTenMean{key}", num(d2["mean"], 3), "Data Set II, N >= 10: PSO-VNS - PSO (pp)")
        m(f"NFDSIINTenWins{key}", str(d2["focus_better"]))
        m(f"NFDSIINTenLosses{key}", str(d2["other_better"]))
        m(f"NFDSIINTenP{key}", pval(d2["wilcoxon_p"]))
        m(f"NFSLSQPGap{key}", num(-r["case_mean_wilcoxon"]["SLSQP"]["mean_dloss_pp"], 3),
          "MS-SLSQP - PSO-VNS case-mean wake loss (pp)")
        m(f"NFSLSQPGapP{key}", pval(r["case_mean_wilcoxon"]["SLSQP"]["p_holm"]))
    m("NFDSIINTenN", str(one["pso_pair"]["ds2_n10"]["n"]))
    m("NFNonTrivN", str(BM["n_nontrivial"]))
    m("NFNonTrivThr", num(NONTRIVIAL_PP, 1))
    m("NFSignChangesOne", str(one["pso_pair"]["sign_changes_vs_recorded"]),
      "cases in which PSO-VNS - PSO changes sign (15 deg recorded -> 1 deg; zero differences not counted)")
    m("NFSignChangesFive", str(five["pso_pair"]["sign_changes_vs_recorded"]))
    m("NFSignChangesBase", str(one["pso_pair"]["n_cases_nonzero_both"]), "cases with a nonzero difference under both")
    m("NFPairMarginRatioOne", num(abs(one["pso_pair"]["all"]["mean_dloss_pp"]) / 0.05, 1))
    # ---- component analysis at 1 deg
    for a in ABL9:
        m(f"NFAblRank{MAC[a]}One", num(one["ablation_friedman"]["avg_rank"][a], 2))
        m(f"NFAblRank{MAC[a]}Fifteen", num(rec["ablation_friedman"]["avg_rank"][a], 2))
    m("NFAblBestOne", LAB[one["ablation_friedman"]["best_ranked"]])
    for k in ABL_CONTR:
        kk = f"{k[0]}-{k[1]}"
        for key, r in (("Fifteen", rec), ("One", one)):
            c = r["ablation_contrasts"][kk]
            m(f"NFAbl{MAC[k[0]]}vs{MAC[k[1]]}{key}", num(c["mean_dloss_pp"], 3), f"{kk} case-mean (pp)")
            m(f"NFAbl{MAC[k[0]]}vs{MAC[k[1]]}{key}P", pval(c["p"]))
    # ---- Horns Rev
    h = HR
    P1 = HR_PRIMARY
    m("NFHRInstFive", num(h["installed"]["5deg_2.5"], 2), "installed 16-turbine block, 5-deg bins (paper)")
    m("NFHRInstOne", num(h["installed"][P1], 2), "installed block, 1-deg bins")
    m("NFHRInstLossFive", num(h["installed_loss_pct"]["5deg_2.5"], 2))
    m("NFHRInstLossOne", num(h["installed_loss_pct"][P1], 2))
    m("NFHRInstLossFivePhaseZero", num(h["installed_loss_pct"]["5deg_0"], 2), "installed, 5-deg bins centred at 0 deg")
    m("NFHRInstDrop", num(h["installed_drop_gwh"], 2), "GWh/yr lost by the installed block, 5 -> 1 deg")
    m("NFHRMeanDrop", num(h["mean_drop_gwh_all"], 2), "mean GWh/yr lost by the feasible optimized layouts")
    dm = {a: v for a, v in h["mean_drop_gwh_by_method"].items() if np.isfinite(v)}
    m("NFHRDropMax", num(max(dm.values()), 2)); m("NFHRDropMaxBy", LAB[max(dm, key=dm.get)])
    m("NFHRDropMin", num(min(dm.values()), 2)); m("NFHRDropMinBy", LAB[min(dm, key=dm.get)])
    for a in HR_ORDER:
        if a in dm:
            m(f"NFHRDrop{MAC[a]}", num(dm[a], 2))
    m("NFHRNRuns", str(h["n_runs"])); m("NFHRNFeas", str(h["n_feasible"]))
    m("NFHRAboveFiveTotal", str(h["above_installed_total"]["5deg_2.5"]))
    m("NFHRAboveOneTotal", str(h["above_installed_total"][P1]), "feasible runs above the installed block, 1 deg")
    m("NFHRAboveOneIntTotal", str(h["above_installed_total"]["1deg_0"]))
    if "PyWakeNOJ_1deg" in h["above_installed_total"]:
        m("NFHRAbovePyWakeTotal", str(h["above_installed_total"]["PyWakeNOJ_1deg"]), "PyWake NOJ, 1 deg")
        m("NFHRInstPyWakeOne", num(h["installed"]["pywake_noj_1deg"], 2))
    for (bud, ini), sn in HR_SET_NAME.items():
        for a in HR_ORDER:
            e = h["methods"].get(f"{bud}_{ini}_{a}")
            if e is None:
                continue
            m(f"NFHR{sn}{MAC[a]}Five", num(e["mean_5deg_2.5"], 2))
            m(f"NFHR{sn}{MAC[a]}One", num(e[f"mean_{P1}"], 2))
            m(f"NFHR{sn}{MAC[a]}AboveFive", str(e["above_5deg_2.5"]))
            m(f"NFHR{sn}{MAC[a]}AboveOne", str(e[f"above_{P1}"]))
            m(f"NFHR{sn}{MAC[a]}Feas", f"{e['feasible']}/{e['runs']}")
        for c, cn in (("5deg_2.5", "Five"), (P1, "One")):
            hi = h["highest_mean"][f"{bud}_{ini}_{c}"]
            m(f"NFHR{sn}Best{cn}", LAB[hi["best_qualified"]], "highest mean AEP (methods with >= half feasible runs)")
            m(f"NFHR{sn}Best{cn}AEP", num(hi["best_qualified_aep"], 2))
            m(f"NFHR{sn}Second{cn}", LAB[hi["second_qualified"]])
            m(f"NFHR{sn}Second{cn}AEP", num(hi["second_qualified_aep"], 2))
            m(f"NFHR{sn}BestAny{cn}", LAB[hi["best_any"]])
        if "PyWakeNOJ_1deg" in h["above_installed_total"]:
            hi = h["highest_mean"][f"{bud}_{ini}_PyWakeNOJ_1deg"]
            m(f"NFHR{sn}BestPyWake", LAB[hi["best_qualified"]]); m(f"NFHR{sn}BestPyWakeAEP", num(hi["best_qualified_aep"], 2))
    for bud, sn in ((6030, "SixK"), (30030, "ThirtyK")):
        for c, cn in (("5deg_2.5", "Five"), (P1, "One"), ("PyWakeNOJ_1deg", "PyWake")):
            q = h["pairs_joint_feasible"].get(f"{bud}_random_PSOC_{c}")
            if q is None:
                continue
            m(f"NFHRPair{sn}{cn}", num(q["mean_diff_gwh"], 2, sign=True), "PSO-VNS - PSO mean AEP (GWh/yr), jointly feasible seeds")
            m(f"NFHRPair{sn}{cn}P", pval(q["p"]))
            m(f"NFHRPair{sn}{cn}N", str(q["n"]))
            m(f"NFHRPair{sn}{cn}Wins", str(q["focus_better"]))
    # 30,030: is PSO-VNS significantly better than each other method (paper test / jointly feasible AEP) at 1 deg?
    pt = h["paper_test"][f"30030_random_{P1}"]
    worse = [LAB[b] for b, v in pt.items() if v["outcome"] != "W"]
    m("NFHRThirtyKNotSigOne", listing(worse), "methods PSO-VNS does NOT beat at 30,030 (paper's run-level test, Holm), 1 deg")
    m("NFHRThirtyKMaxPOne", pval(max(v["p_holm"] for v in pt.values())))
    pt5 = h["paper_test"]["30030_random_5deg_2.5"]
    m("NFHRThirtyKMaxPFive", pval(max(v["p_holm"] for v in pt5.values())))
    jf = {b: h["pairs_joint_feasible"].get(f"30030_random_{b}_{P1}") for b in HR_ORDER if b != FOCUS}
    m("NFHRThirtyKJointNotSigOne", listing(LAB[b] for b, v in jf.items() if v and not (v["p"] < 0.05 and v["mean_diff_gwh"] > 0)),
      "methods not significantly worse than PSO-VNS on the jointly feasible seeds at 30,030, 1 deg (unadjusted)")
    # ---- PyWake
    if PW is not None:
        def pw(farm, bins, model, src="PyWake"):
            r = PW[(PW.Farm == farm) & (PW.Bins == bins) & (PW.Model == model) & (PW.Source == src)]
            return r.iloc[0] if len(r) else None
        v = pw("HR80", "ours_5deg_2.5", "NOJ_k0.04")
        o = pw("HR80", "ours_5deg_2.5", "paper_model", "ours")
        m("NFPyWakeVersion", str(v.Version))
        m("NFPyWakeAEP", num(v.AEP_GWh, 2), "PyWake NOJ(k=0.04), 80 turbines, our 5-deg bins (GWh/yr)")
        m("NFPyWakeLossPct", num(v.WakeLossPct, 2))
        m("NFPyWakeOurAEP", num(o.AEP_GWh, 2)); m("NFPyWakeOurLossPct", num(o.WakeLossPct, 2))
        m("NFPyWakeDiffPct", num(100 * (o.AEP_GWh / v.AEP_GWh - 1), 2), "our AEP relative to PyWake (%)")
        m("NFPyWakeLossDiffPP", num(o.WakeLossPct - v.WakeLossPct, 2, sign=True), "our wake loss - PyWake wake loss (pp)")
        m("NFPyWakeLossDiffRel", num(100 * (o.WakeLossPct / v.WakeLossPct - 1), 0, sign=True), "% of PyWake's wake loss")
        m("NFPyWakeIdealAEP", num(v.IdealAEP_GWh, 2)); m("NFPyWakeIdealDiff", num(abs(v.IdealAEP_GWh - o.IdealAEP_GWh), 4))
        d = pw("HR80", "pywake_default_1deg_0", "NOJ_k0.04")
        m("NFPyWakeDefaultAEP", num(d.AEP_GWh, 2)); m("NFPyWakeDefaultLossPct", num(d.WakeLossPct, 2))
        rc = pw("HR80", "ours_5deg_2.5", "NOJ_k0.04_RotorCenter")
        m("NFPyWakeRotorCenterAEP", num(rc.AEP_GWh, 2)); m("NFPyWakeRotorCenterLossPct", num(rc.WakeLossPct, 2))
        rm = pw("HR80", "ours_5deg_2.5", "NOJ_k0.04_RotorCenter_momentum")
        if rm is not None:
            m("NFPyWakeRotorCenterMomAEP", num(rm.AEP_GWh, 2)); m("NFPyWakeRotorCenterMomLossPct", num(rm.WakeLossPct, 2))
        lc = pw("HR80", "ours_5deg_2.5", "paper_model_localCT", "ours")
        m("NFOurLocalCTAEP", num(lc.AEP_GWh, 2), "our model with C_T at the local (waked) speed")
        m("NFOurLocalCTLossPct", num(lc.WakeLossPct, 2))
        o1 = pw("HR80", "1deg_0.5", "paper_model", "ours"); v1 = pw("HR80", "1deg_0.5", "NOJ_k0.04")
        m("NFPyWakeOneAEP", num(v1.AEP_GWh, 2)); m("NFPyWakeOneLossPct", num(v1.WakeLossPct, 2))
        m("NFPyWakeOurOneAEP", num(o1.AEP_GWh, 2)); m("NFPyWakeOurOneLossPct", num(o1.WakeLossPct, 2))
        m("NFPyWakeOneDiffPct", num(100 * (o1.AEP_GWh / v1.AEP_GWh - 1), 2))
        s = pw("HR16", "ours_5deg_2.5", "NOJ_k0.04"); so = pw("HR16", "ours_5deg_2.5", "paper_model", "ours")
        m("NFPyWakeSixteenAEP", num(s.AEP_GWh, 2)); m("NFPyWakeSixteenLossPct", num(s.WakeLossPct, 2))
        m("NFPyWakeSixteenOurAEP", num(so.AEP_GWh, 2)); m("NFPyWakeSixteenOurLossPct", num(so.WakeLossPct, 2))
        m("NFPyWakeSixteenDiffPct", num(100 * (so.AEP_GWh / s.AEP_GWh - 1), 2))
        m("NFPyWakePMaxDiff", f"{float(PW[PW.Source == 'PyWake'].PMaxAbsDiff.max()):.0e}".replace("e-", "\\times10^{-") + "}",
          "max |P_PyWake - P_ours| over (wd, ws) bins")
    # ---- IEA37 projected
    if IP is not None:
        chg = IP["changed_status"]
        m("NFIEAChanged", listing(f"participant~{c['participant'].replace('par', '')} ({c['turbines']} turbines, "
                                  f"{c['max_excess_m'] * 1000:.1f}~mm)" if c['max_excess_m'] < 0.1 else
                                  f"participant~{c['participant'].replace('par', '')} ({c['turbines']} turbines, "
                                  f"{c['max_excess_m']:.1f}~m)" for c in chg), "participants feasible only after projection")
        m("NFIEANChanged", str(len(chg)))
        m("NFIEATolStrict", "1~mm")
        for n, nm in ((16, "Sixteen"), (36, "ThirtySix")):
            s = IP["scenarios"][str(n)]
            m(f"NFIEABestStrict{nm}", num(s["strict"]["best_aep"], 1)); m(f"NFIEABestStrictBy{nm}", s["strict"]["best_by"].replace("par", ""))
            m(f"NFIEABestProj{nm}", num(s["projected"]["best_aep"], 1)); m(f"NFIEABestProjBy{nm}", s["projected"]["best_by"].replace("par", ""))
            m(f"NFIEABestProjLoss{nm}", num(s["projected"]["best_loss_pct"], 2))
            m(f"NFIEABestStrictLoss{nm}", num(s["strict"]["best_loss_pct"], 2))
            m(f"NFIEAPubFeasProj{nm}", str(s["projected"]["n_feasible"])); m(f"NFIEAPubFeasStrict{nm}", str(s["strict"]["n_feasible"]))
            m(f"NFIEAProjExcess{nm}", num(s["projected"]["best_max_excess_m"], 3), "m, boundary excess of the projected best layout")
            m(f"NFIEAProjSpacing{nm}", num(s["projected"]["best_min_spacing_m"], 0))
            for b, bn in ((6030, "SixK"), (30030, "ThirtyK")):
                g = s["ours"].get(str(b))
                if not g:
                    continue
                for conv, cn in (("strict", "Strict"), ("projected", "Proj")):
                    x = g[conv]
                    m(f"NFIEAGap{nm}{bn}{cn}", num(x["best_gap_pct"], 2), f"best PSO-VNS vs best feasible published, AEP %, {conv}")
                    m(f"NFIEAMeanGap{nm}{bn}{cn}", num(x["mean_gap_pct"], 2))
                    m(f"NFIEAGapLoss{nm}{bn}{cn}", num(x["best_gap_loss_pp"], 2, sign=True), "wake loss PSO-VNS best - published best (pp)")
                    m(f"NFIEAMeanGapLoss{nm}{bn}{cn}", num(x["mean_gap_loss_pp"], 2, sign=True))
                    m(f"NFIEARatio{nm}{bn}{cn}", num(x["best_ratio_pct"], 1), "best PSO-VNS AEP as % of best feasible published")
                    m(f"NFIEARank{nm}{bn}{cn}", ORD[x["best_rank"]], f"rank of best PSO-VNS among {x['n_feasible_published']} feasible published + it")
                    m(f"NFIEAMeanRank{nm}{bn}{cn}", ORD[x["mean_rank"]])
                    m(f"NFIEAOurBestGap{nm}{bn}{cn}", num(x["our_best_gap_pct"], 2), "best of all our methods")
                    m(f"NFIEAOurBestBy{nm}{bn}", LAB.get(g["our_best_method"], g["our_best_method"]))
                m(f"NFIEAPSOVNSLoss{nm}{bn}", num(g["psovns_best_loss_pct"], 2))
                m(f"NFIEAPSOVNSMeanLoss{nm}{bn}", num(g["psovns_mean_loss_pct"], 2))
    hdr = [f"% generated by mpce_direction.py ({stamp}) -- do not edit by hand",
           "% Direction-resolution robustness (Phase 6, D16; R4 issues 1-3, 5). Benchmark: every 15-deg bin split into 15",
           "% 1-deg sub-bins with the same Weibull (psi_j, k_j) and frequency w_j/15 (re-evaluation, not re-optimization).",
           "% Horns Rev: 1-deg bins centred at 0.5, 1.5, ... (sector frequency /30, same sector Weibull). Pair macros: PSO-VNS",
           "% minus PSO (negative = PSO-VNS lower wake loss) for the benchmark; PSO-VNS minus PSO AEP (positive = PSO-VNS",
           "% higher AEP) for Horns Rev. Checks: mpce_check_dir.py (F01...)."]
    open(fn, "w").write("\n".join(hdr + m.lines) + "\n")
    return len(m.lines)


def write_supp(BM, HR, PW, IP, fn, stamp):
    R = BM["res"]
    rec, one, five, g1, g15 = R["rec"], R["J15"], R["J3"], R["G1"], R["G15"]
    T = [f"% generated by mpce_direction.py ({stamp}) -- do not edit by hand", "% requires booktabs", ""]
    # ---- tab:F-bench
    lines = []
    for a in MAIN8:
        cw = {k: r["case_mean_wilcoxon"].get(a) for k, r in (("rec", rec), ("one", one))}
        c = [LAB[a]] + [tnum(r["mean_loss_qualified"][a], 3) for r in (rec, five, one)] + \
            [tnum(r["friedman"]["avg_rank"][a], 2) for r in (rec, one, g1, g15)]
        if a == FOCUS:
            c += ["--", "--", "--", "--"]
        else:
            c += [tnum(cw["rec"]["mean_dloss_pp"], 3, True), tp(cw["rec"]["p_holm"]), tnum(cw["one"]["mean_dloss_pp"], 3, True), tp(cw["one"]["p_holm"])]
        lines.append(" & ".join(c) + " \\\\")
    lines.append("\\midrule")
    lines.append(f"Best-ranked & & & & {LAB[rec['friedman']['best_ranked']]} & {LAB[one['friedman']['best_ranked']]} & "
                 f"{LAB[g1['friedman']['best_ranked']]} & {LAB[g15['friedman']['best_ranked']]} & & & & \\\\")
    lines.append(f"$\\bar\\tau$ / same best (\\%) & & & & -- & {one['tau_within_case']:.2f} / {one['same_best_pct']:.0f} & "
                 f"{g1['tau_within_case']:.2f} / {g1['same_best_pct']:.0f} & {g15['tau_within_case']:.2f} / {g15['same_best_pct']:.0f} & & & & \\\\")
    T.append(table("table*", "Direction resolution of the benchmark: all feasible final layouts of the 68 cases (6,030 evaluations, random "
                   "initialization, 30 runs) re-evaluated, not re-optimized, with every 15$^\\circ$ direction bin split into 3 (5$^\\circ$) or 15 "
                   "(1$^\\circ$) sub-bins. Each sub-bin keeps the Weibull parameters of its bin and receives the bin frequency divided by the "
                   "number of sub-bins (piecewise-constant wind rose; wake-free power unchanged). Mean wake loss (\\%) over the cases in which "
                   "the method qualifies (at least 15 feasible runs); average rank over the 68 cases (ranking rule of Table~\\ref{M-tab:friedman68}) "
                   "for the Jensen benchmark model and the Gaussian wake ($k^*=0.04$) at 15$^\\circ$ and 1$^\\circ$; $\\Delta L$: mean over the cases "
                   "of the PSO-VNS minus method case-mean wake loss (pp; negative = PSO-VNS better) and Holm-adjusted case-mean Wilcoxon $p$. "
                   "$\\bar\\tau$: mean within-case Kendall correlation with the benchmark ordering; same best: cases with unchanged best method.",
                   "tab:F-bench", "lccccccccccc",
                   "& \\multicolumn{3}{c}{Mean wake loss (\\%)} & \\multicolumn{4}{c}{Average rank} & \\multicolumn{2}{c}{$\\Delta L$, 15$^\\circ$} & "
                   "\\multicolumn{2}{c}{$\\Delta L$, 1$^\\circ$} \\\\\n\\cmidrule(lr){2-4}\\cmidrule(lr){5-8}\\cmidrule(lr){9-10}\\cmidrule(lr){11-12}\n"
                   "Method & 15$^\\circ$ & 5$^\\circ$ & 1$^\\circ$ & Jensen 15$^\\circ$ & Jensen 1$^\\circ$ & Gauss 15$^\\circ$ & Gauss 1$^\\circ$ & "
                   "$\\Delta L$ & $p_{\\rm Holm}$ & $\\Delta L$ & $p_{\\rm Holm}$", lines, resize=True))
    # ---- tab:F-equiv
    lines = []
    for nm, r in (("15$^\\circ$ (benchmark)", rec), ("15$^\\circ$ (re-evaluated)", R["J1"]), ("5$^\\circ$", five), ("1$^\\circ$", one)):
        pp = r["pso_pair"]
        for sub, e, n, p in (("all", pp["all"], BM["n_cases"], pp["all_wilcoxon"]["p"]),
                             ("non-trivial", pp["nontrivial"], pp["n_nontrivial"], pp["nontrivial_wilcoxon_p"])):
            lines.append(f"{nm} & {sub} ({n}) & {tnum(e['mean_dloss_pp'], 3, True)} & {ci_tab(e['ci90_mean_dloss_pp'])} & "
                         f"{e['min_margin_pp']:.3f} & {'yes' if e['equivalent'] else 'no'} & {tp(p)} & "
                         + (f"{pp['focus_better']}/{pp['ties']}/{pp['other_better']} & {pp['abs_gt_margin']} & "
                            f"{pp['sign_changes_vs_recorded'] if r is not rec else '--'}" if sub == "all" else "& & ")
                         + " \\\\")
        d2 = pp["ds2_n10"]
        lines.append(f"{nm} & DS II, $N\\ge10$ ({d2['n']}) & {tnum(d2['mean'], 3, True)} & {ci_tab(d2['ci95'])}$^{{b}}$ & & & {tp(d2['wilcoxon_p'])} & "
                     f"{d2['focus_better']}/{d2['n'] - d2['focus_better'] - d2['other_better']}/{d2['other_better']} & & \\\\")
        if r is not one:
            lines.append("\\midrule")
    T.append(table("table*", f"PSO-VNS vs.\\ PSO under the benchmark direction bins and finer sub-bins (re-evaluation): mean over the cases of the "
                   f"per-case difference of the mean wake losses (PSO-VNS minus PSO, pp; negative = PSO-VNS better), 90\\% percentile case-bootstrap "
                   f"CI (10,000 resamples), smallest equivalence margin, equivalence at $\\pm${BM_MARGIN:.2f}~pp (90\\% CI inside the margin; TOST at "
                   f"$\\alpha=0.05$), case-mean Wilcoxon $p$, cases with lower/equal/higher PSO-VNS loss, cases with $|\\Delta L|>{BM_MARGIN:.2f}$~pp, and "
                   f"cases whose sign differs from the benchmark (zero differences not counted). Non-trivial: the {BM['n_nontrivial']} cases with a "
                   f"PSO-VNS wake loss of at least {NONTRIVIAL_PP:.1f}\\% in the benchmark model. $^{{b}}$95\\% CI.",
                   "tab:F-equiv", "llccccccccc",
                   "Direction bins & Cases & $\\Delta L$ (pp) & 90\\% CI & Min.\\ margin & Equiv. & $p$ & W/T/L & $>$margin & Sign changes",
                   lines, resize=True))
    # ---- tab:F-abl
    lines = []
    for a in ABL9:
        lines.append(f"{LAB[a]} & {rec['mean_loss_qualified_abl'][a]:.3f} & {one['mean_loss_qualified_abl'][a]:.3f} & "
                     f"{rec['ablation_friedman']['avg_rank'][a]:.2f} & {one['ablation_friedman']['avg_rank'][a]:.2f} \\\\")
    lines.append("\\midrule")
    for a, b in ABL_CONTR:
        k = f"{a}-{b}"
        c0, c1 = rec["ablation_contrasts"][k], one["ablation_contrasts"][k]
        lines.append(f"{LAB[a]} $-$ {LAB[b]} & \\multicolumn{{2}}{{c}}{{{tnum(c0['mean_dloss_pp'], 3, True)} ({tp(c0['p'])})}} & "
                     f"\\multicolumn{{2}}{{c}}{{{tnum(c1['mean_dloss_pp'], 3, True)} ({tp(c1['p'])})}} \\\\")
    T.append(table("table", "Component analysis under 15$^\\circ$ (benchmark) and 1$^\\circ$ direction bins (re-evaluation of the same final "
                   "layouts). Top: mean wake loss (\\%, qualified cases) and average rank among the nine variants. Bottom: case-mean "
                   "difference of the wake losses (first minus second method, pp) and unadjusted case-mean Wilcoxon $p$.",
                   "tab:F-abl", "lcccc", "& \\multicolumn{2}{c}{Loss (\\%)} & \\multicolumn{2}{c}{Avg.\\ rank} \\\\\n"
                   "\\cmidrule(lr){2-3}\\cmidrule(lr){4-5}\nVariant / contrast & 15$^\\circ$ & 1$^\\circ$ & 15$^\\circ$ & 1$^\\circ$", lines))
    # ---- tab:F-hr
    h = HR
    P1 = HR_PRIMARY
    has_pw = "PyWakeNOJ_1deg" in h["above_installed_total"]
    lines = [f"Installed block & -- & {h['installed']['5deg_2.5']:.2f} & {h['installed'][P1]:.2f}"
             + (f" & {h['installed']['pywake_noj_1deg']:.2f}" if has_pw else "") + " & & & " + ("& " if has_pw else "") + "\\\\", "\\midrule"]
    for (bud, ini) in HR_SET:
        for a in HR_ORDER:
            e = h["methods"].get(f"{bud}_{ini}_{a}")
            if e is None:
                continue
            c = [HR_SET_LAB[(bud, ini)] if a == HR_ORDER[0] else "", LAB[a], f"{e['feasible']}/{e['runs']}",
                 tnum(e["mean_5deg_2.5"] if e["feasible"] else np.nan, 2), tnum(e[f"mean_{P1}"] if e["feasible"] else np.nan, 2)]
            if has_pw:
                c.append(tnum(e["mean_PyWakeNOJ_1deg"] if e["feasible"] else np.nan, 2))
            c += [str(e["above_5deg_2.5"]), str(e[f"above_{P1}"])] + ([str(e["above_PyWakeNOJ_1deg"])] if has_pw else [])
            lines.append(" & ".join(c) + " \\\\")
        lines.append("\\midrule")
    lines = lines[:-1]
    pwh = " & PyWake NOJ 1$^\\circ$" if has_pw else ""
    T.append(table("table*", "Horns Rev~1, 16-turbine block: every run of the study re-evaluated (not re-optimized) with 1$^\\circ$ direction bins "
                   "(centres 0.5$^\\circ$, 1.5$^\\circ$, \\ldots; each 30$^\\circ$ sector keeps its Weibull $A$, $k$ and spreads its frequency "
                   "uniformly over its 30 bins; 1~m/s speed bins unchanged)" + (" and with PyWake's NOJ model ($k=0.04$) at the same 1$^\\circ$ bins" if has_pw else "")
                   + ". Mean AEP (GWh/yr) of the feasible runs and number of feasible runs whose AEP exceeds that of the installed block "
                   "under the same model. R/F: random / feasibility-preserving initialization; 6k, 30k, 120k: evaluations (10 seeds at 120k).",
                   "tab:F-hr", "llc" + "c" * (3 if has_pw else 2) + "c" * (3 if has_pw else 2),
                   "& & & \\multicolumn{%d}{c}{Mean AEP (GWh/yr)} & \\multicolumn{%d}{c}{Runs above installed} \\\\\n"
                   "\\cmidrule(lr){4-%d}\\cmidrule(lr){%d-%d}\nSetting & Method & Feas. & 5$^\\circ$ (paper) & 1$^\\circ$%s & 5$^\\circ$ & 1$^\\circ$%s"
                   % ((3, 3, 6, 7, 9, pwh, " & PyWake") if has_pw else (2, 2, 5, 6, 7, "", "")), lines, sep="2pt"))
    # ---- tab:F-hrpair
    lines = []
    for bud in (6030, 30030):
        for b in HR_ORDER:
            if b == FOCUS:
                continue
            cells = []
            for c in ["5deg_2.5", P1] + (["PyWakeNOJ_1deg"] if has_pw else []):
                q = h["pairs_joint_feasible"].get(f"{bud}_random_{b}_{c}")
                cells.append("--" if q is None else f"{tnum(q['mean_diff_gwh'], 2, True)} ({tp(q['p'])})")
            q = h["pairs_joint_feasible"].get(f"{bud}_random_{b}_{P1}")
            lines.append(f"{bud:,} & {LAB[b]} & {q['n'] if q else 0} & " + " & ".join(cells) + " \\\\")
        if bud == 6030:
            lines.append("\\midrule")
    T.append(table("table", "Horns Rev~1, 16 turbines, random initialization: PSO-VNS minus each method, mean AEP difference (GWh/yr; positive = "
                   "PSO-VNS better) on the seeds on which both runs are feasible, with the unadjusted Wilcoxon signed-rank $p$, under the paper's "
                   "5$^\\circ$ bins, 1$^\\circ$ bins" + (" and PyWake NOJ (1$^\\circ$)" if has_pw else "") + ". $n$: jointly feasible seeds.",
                   "tab:F-hrpair", "llc" + "c" * (3 if has_pw else 2),
                   "Budget & vs.\\ & $n$ & 5$^\\circ$ & 1$^\\circ$" + (" & PyWake" if has_pw else ""), lines, sep="2pt"))
    # ---- tab:F-pywake
    if PW is not None:
        lines = []
        blab = {"ours_5deg_2.5": "5$^\\circ$ at 2.5$^\\circ$ (paper)", "1deg_0.5": "1$^\\circ$ at 0.5$^\\circ$",
                "pywake_default_1deg_0": "1$^\\circ$ at 0$^\\circ$ (PyWake default)"}
        mlab = {"paper_model": "This paper (hornsrev\\_model)", "paper_model_localCT": "This paper, $C_T$ at local speed$^{a}$",
                "NOJ_k0.04": "PyWake NOJ", "NOJ_k0.04_RotorCenter": "PyWake NOJ, rotor centre",
                "NOJ_k0.04_RotorCenter_momentum": "PyWake NOJ, rotor centre, momentum $a(C_T)$"}
        for farm in ("HR80", "HR16"):
            for b in ("ours_5deg_2.5", "1deg_0.5", "pywake_default_1deg_0"):
                for _, r in PW[(PW.Farm == farm) & (PW.Bins == b)].iterrows():
                    lines.append(f"{'80' if farm == 'HR80' else '16'} & {blab[b]} & {mlab.get(r.Model, r.Model)} & {r.AEP_GWh:.2f} & "
                                 f"{r.IdealAEP_GWh:.2f} & {r.WakeLossPct:.2f} \\\\")
            if farm == "HR80":
                lines.append("\\midrule")
        v = PW.iloc[0]
        T.append(table("table", f"Horns Rev~1 wake model against PyWake {v.Version} (NOJ, $k=0.04$, V80 curves and Weibull climate of "
                       f"PyWake's Hornsrev1Site) at identical direction and speed bins (the bin probabilities agree to "
                       f"{float(PW[PW.Source == 'PyWake'].PMaxAbsDiff.max()):.0e}). AEP and wake-free AEP in GWh/yr, wake loss in \\%. PyWake's NOJ "
                       "uses rotor-area overlap, $C_T$ at the effective speed of the upstream turbine and the Madsen $a(C_T)$ polynomial by default. "
                       "$^{a}$Diagnostic variant of our model only.", "tab:F-pywake", "cllccc",
                       "$N$ & Direction bins & Model & AEP & Wake-free & Loss (\\%)", lines, sep="2pt"))
    # ---- tab:F-iea
    if IP is not None:
        lines = []
        for n in (16, 36):
            for p in IP["participants"]:
                if p["turbines"] != n:
                    continue
                lines.append(f"{n} & {p['participant'].replace('par', '')} & {p['aep_published']:,.1f} & {p['max_excess_m']:.4f} & "
                             f"{p['min_spacing_m']:.1f} & {'yes' if p['feasible_strict'] else 'no'} & "
                             f"{tnum(p['aep_projected'], 1) if p['projected_ok'] else '--'} & "
                             f"{p['min_spacing_projected_m']:.1f} & {'yes' if p['feasible_projected'] else 'no'} \\\\".replace(",", "{,}"))
            if n == 16:
                lines.append("\\midrule")
        T.append(table("table*", "IEA37 Case Study~1 published layouts: largest distance outside the boundary circle, minimum spacing and "
                       "feasibility at the published 1~mm tolerance (strict), and after projecting every turbine outside the circle radially "
                       "onto it (AEP recomputed with our calculator, which reproduces the official one; spacing re-checked against 2$D$ = "
                       "260~m with a $10^{-6}$~m tolerance). Layouts violating the spacing cannot be repaired by the projection.",
                       "tab:F-iea", "ccccccccc", "$N$ & Part. & AEP publ.\\ (MWh) & Excess (m) & Min.\\ spacing (m) & Feas.\\ strict & "
                       "AEP projected & Spacing proj.\\ (m) & Feas.\\ projected", lines, sep="2.5pt"))
        lines = []
        for n in (16, 36):
            s = IP["scenarios"][str(n)]
            for b in (6030, 30030):
                g = s["ours"].get(str(b))
                if not g:
                    continue
                for conv in ("strict", "projected"):
                    x = g[conv]
                    lines.append(f"{n} & {b:,} & {conv} & {x['ref_by'].replace('par', '')} & {x['ref_aep']:,.1f} & {x['ref_loss_pct']:.2f} & "
                                 f"{tnum(x['best_gap_pct'], 2, True)} & {tnum(x['best_gap_loss_pp'], 2, True)} & {ORD[x['best_rank']]} & "
                                 f"{tnum(x['mean_gap_pct'], 2, True)} & {tnum(x['mean_gap_loss_pp'], 2, True)} & "
                                 f"{tnum(x['our_best_gap_pct'], 2, True)} ({LAB.get(g['our_best_method'], g['our_best_method'])}) \\\\".replace(",", "{,}", 2))
            if n == 16:
                lines.append("\\midrule")
        T.append(table("table*", "IEA37 Case Study~1: gap of PSO-VNS (best and mean of 30 runs) and of the best run of all our methods to the "
                       "best feasible published layout, under the strict 1~mm tolerance and with boundary-projected published layouts; AEP gap "
                       "in \\% and wake-loss gap in percentage points (pp, positive = more wake loss than the reference). Rank: position of the "
                       "best PSO-VNS run among the feasible published layouts.", "tab:F-iea-gap", "cclcccccccccc",
                       "$N$ & Budget & Convention & Ref. & Ref.\\ AEP & Ref.\\ loss (\\%) & Best gap (\\%) & Best gap (pp) & Rank & "
                       "Mean gap (\\%) & Mean gap (pp) & Best of ours (\\%)", lines, sep="2pt", resize=True))
    open(fn, "w").write("\n".join(T))


BM_MARGIN = 0.05


def main(argv=None):
    global BM_MARGIN
    ap = argparse.ArgumentParser()
    ap.add_argument("--procs", type=int, default=2)
    ap.add_argument("--data-dir", default=HERE)
    ap.add_argument("--out-dir", default=HERE)
    args = ap.parse_args(argv)
    import mpce_results as M                     # read-only use of the pipeline's loader and statistics
    M.set_focus(FOCUS)
    BM_MARGIN = M.EQ_MARGIN
    LOG = []

    def log(s=""):
        print(s, flush=True); LOG.append(s)
    t0 = time.time()
    stamp = time.strftime("%Y-%m-%d %H:%M:%S")
    ALL, avail, fb, _ = M.load(args.data_dir, False)
    G = ALL[(ALL.Budget == 6030) & (ALL.Init == "random") & ALL.Dataset.isin(["1", "2"]) & ALL.Algorithm.isin(ALLM)].reset_index(drop=True)
    log(f"[1] Benchmark: {len(G)} runs ({G.Feasible.sum()} feasible), {G.groupby(['Dataset', 'Radius', 'Turbines']).ngroups} cases")
    F = bench_reevaluate(G, args.procs, log)
    BM, Fm = bench_block(M, G, F, log)
    log(f"  reproduction (S = 1 vs recorded objective): {BM['reproduction']}")
    log("  mean wake loss 15 -> 5 -> 1 deg: " + ", ".join(
        f"{LAB[a]} {BM['res']['rec']['mean_loss_qualified'][a]:.3f}/{BM['res']['J3']['mean_loss_qualified'][a]:.3f}/"
        f"{BM['res']['J15']['mean_loss_qualified'][a]:.3f}" for a in MAIN8))
    log("\n[2] Horns Rev 16")
    HR, HX = hr_block(M, ALL, log)
    log("\n[3] PyWake check (pywake_check.csv)")
    pwf = os.path.join(HERE, "pywake_check.csv")
    PW = pd.read_csv(pwf) if os.path.exists(pwf) else None
    log("  " + ("present" if PW is not None else "MISSING -> PyWake macros not written; the agreement claim must be dropped"))
    log("\n[4] IEA37 projected published layouts")
    import iea37_projected as IPm
    IP = IPm.run(M, ALL, log)
    summary = dict(generated=stamp, eq_margin_pp=M.EQ_MARGIN, benchmark=BM, horns_rev=HR,
                   pywake=None if PW is None else PW.to_dict(orient="records"), iea37_projected=IP,
                   data_availability={k: v.get("status") for k, v in avail.items()}, fallbacks=fb, log=LOG)
    json.dump(jclean(summary), open(os.path.join(args.out_dir, "mpce_summary_dir.json"), "w"), indent=1)
    keep = ["Algorithm", "Dataset", "Radius", "Turbines", "Seed", "Objective", "Ideal", "J1", "J3", "J15", "I1", "G1", "G15"]
    Fm[keep].to_csv(os.path.join(args.out_dir, "mpce_direction_layouts.csv.gz"), index=False, float_format="%.8g")
    HX.drop(columns=[]).to_csv(os.path.join(args.out_dir, "mpce_direction_hr16.csv"), index=False, float_format="%.8g")
    n = write_macros(BM, HR, PW, IP, os.path.join(args.out_dir, "mpce_numbers_dir.tex"), stamp)
    write_supp(BM, HR, PW, IP, os.path.join(args.out_dir, "mpce_supp_direction.tex"), stamp)
    log(f"\nwrote mpce_summary_dir.json, mpce_numbers_dir.tex ({n} macros), mpce_supp_direction.tex in {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
