"""Revision 3 (item C2 / R8): analysis of experiment rev3_fine (direct optimization with 1-deg direction bins).

Usage:  python3 rev3_fine_analysis.py [--data-dir D] [--out-dir O]
Inputs (--data-dir, default this folder): rev3_fine_s<i>of<k>.csv (primary arm, required) and, if present,
rev3_fine15_s<i>of<k>.csv (optional 15-deg control arm); a complete shard set is preferred, otherwise all shard files
are merged and the output is marked PARTIAL. Reference arm: the stored main-comparison runs of the same five methods,
cases and seeds 1-30 (mpce_psobv, mpce_psoc, mpce_slsqp, mpce_rsdisc, rev2_ga; always read from this folder), whose
feasible final layouts are re-evaluated here with mpce_direction.bench_objective(xy, ds, 15) (the evaluator of the
paper's 1-deg re-evaluation; infeasible rows keep their label and never enter a mean).
Outputs (--out-dir, default this folder): rev3_fine.json, rev3_fine_tables.tex (tab:S-r3-fine).

Statistics (prespecified in rev3_fine_manifest.md; functions of mpce_results (MR), mpce_inference_extra (MX) and
rev2_analysis (RA) reused unchanged): wake loss L = 100 (Ideal - Objective) / Ideal under the 1-deg objective;
feasibility-aware ranks (MR.case_stats / rank_rule), Friedman + Holm z (MR.friedman_block); PSO-VNS - PSO contrast
(RA.contrast: case-level bootstrap 90 % CI + TOST (MR.tost), seed-level within-case bootstrap 90 % CI + TOST
(MX.seed_level), case-mean Wilcoxon (MX.cm_test), Hodges-Lehmann); all-run paired scores (feasible beats infeasible,
two feasible runs by objective with ties |diff| <= 1e-6, two infeasible runs tie; (W + T/2) / n, 95 % case-bootstrap
interval with 20,000 resamples, exact two-sided sign test over discordant pairs, Holm over the 4 comparisons).
Arms: direct = rev3_fine (seeds 31-60, 1-deg objective); reeval = stored 15-deg runs re-evaluated at 1 deg (seeds
1-30); rec15 = the same stored runs with their recorded 15-deg values; optional ctrl1 / ctrl15 = rev3_fine15 runs
(seeds 31-60) evaluated at 1 deg / 15 deg.
"""
import os, sys, glob, re, json, argparse, warnings
import numpy as np, pandas as pd
from scipy.stats import kendalltau, binomtest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import mpce_results as MR                  # noqa: E402
import mpce_inference_extra as MX          # noqa: E402
import rev2_analysis as RA                 # noqa: E402
import wflop_model as W                    # noqa: E402
import mpce_direction as MD                # noqa: E402

warnings.filterwarnings("ignore", category=RuntimeWarning)
METHODS = ["PSOBV", "PSOC", "GA", "SLSQP", "RSDVNS"]
FOCUS = "PSOBV"
CASES = [("1", 500, 10), ("1", 750, 6), ("1", 750, 12), ("1", 1000, 15),
         ("2", 500, 10), ("2", 750, 6), ("2", 750, 12), ("2", 1000, 15)]
CASE = MR.CASE
SRC = {"PSOBV": "mpce_psobv", "PSOC": "mpce_psoc", "SLSQP": "mpce_slsqp", "RSDVNS": "mpce_rsdisc", "GA": "rev2_ga"}
LAB = {"PSOBV": "PSO-VNS", "PSOC": "PSO", "GA": "GA", "SLSQP": "MS-SLSQP", "RSDVNS": "RSD-VNS"}
EQ = MR.EQ_MARGIN
ALLRUN_B, ALLRUN_SEED, TIE_OBJ = 20000, 20261004, 1e-6
ARM_LAB = {"direct": "direct 1-deg optimization (seeds 31-60)", "reeval": "15-deg layouts re-evaluated at 1 deg (seeds 1-30)",
           "rec15": "15-deg layouts, recorded 15-deg objective (seeds 1-30)",
           "ctrl1": "15-deg runs, seeds 31-60, re-evaluated at 1 deg", "ctrl15": "15-deg runs, seeds 31-60, 15-deg objective"}


def log(m=""):
    print(m, flush=True)


# ------------------------------------------------------------------ data
def read_shards(prefix, data_dir):
    files = glob.glob(os.path.join(data_dir, f"{prefix}_s*of*.csv"))
    by_k = {}
    for fn in files:
        m = re.search(rf"{prefix}_s(\d+)of(\d+)\.csv$", os.path.basename(fn))
        if m:
            by_k.setdefault(int(m.group(2)), {})[int(m.group(1))] = fn
    if not by_k:
        return None, dict(status="missing", rows=0)
    complete = {k: v for k, v in by_k.items() if set(v) == set(range(k))}
    if complete:
        k = max(complete, key=lambda kk: sum(os.path.getsize(f) for f in complete[kk].values()))
        use, status = list(complete[k].values()), f"shard set complete ({k} shards)"
    else:
        use, status = [f for v in by_k.values() for f in v.values()], "PARTIAL (incomplete shard set)"
    parts = []
    for f in sorted(use):
        try:
            parts.append(pd.read_csv(f, float_precision="round_trip"))
        except pd.errors.EmptyDataError:
            pass
    if not parts:
        return None, dict(status="empty", rows=0)
    df = MR.std_cols(pd.concat(parts, ignore_index=True))
    df = df.drop_duplicates(MR.KEY + ["Rose"], keep="last").reset_index(drop=True)
    return df, dict(status=status, rows=int(len(df)), files=[os.path.basename(f) for f in sorted(use)])


def in_cases(D):
    return D[[c in CASES for c in zip(D.Dataset, D.Radius, D.Turbines)]].reset_index(drop=True)


def stored_runs():
    """Stored main-comparison runs (seeds 1-30, 6,030 evaluations, random starts) of the 5 methods on the 8 cases,
    re-evaluated at 1 deg (feasible rows) with mpce_direction.bench_objective(xy, ds, 15)."""
    parts = []
    for alg, src in SRC.items():
        for f in sorted(glob.glob(os.path.join(HERE, src + "_s*of*.csv"))):
            d = pd.read_csv(f, float_precision="round_trip", usecols=lambda c: c != "Curve")
            parts.append(d[d.Algorithm == alg])
    D = MR.std_cols(pd.concat(parts, ignore_index=True))
    D = D[(D.Budget == 6030) & (D.Init == "random") & D.Seed.between(1, 30)]
    D = in_cases(D.drop_duplicates(MR.KEY, keep="first"))
    # consistency with the loader of the paper (mpce_inference_extra.load) for the four methods it contains
    A = in_cases(RA.main_runs())
    A = A[A.Algorithm.isin(METHODS)]
    m = D.merge(A[MR.KEY + ["Objective", "Feasible"]], on=MR.KEY, suffixes=("", "_paper"))
    cons = dict(n_matched=int(len(m)), n_methods_in_paper_loader=int(A.Algorithm.nunique()),
                objective_max_rel_dev=float(((m.Objective - m.Objective_paper).abs() / m.Objective.abs()).max()),
                objective_equal_within_2ulp=bool(np.allclose(m.Objective, m.Objective_paper, rtol=4.5e-16, atol=0)),
                note_objective="the paper loader reads the CSVs with pandas' default float parser (can be off by 1 ulp)",
                feasible_equal=bool((m.Feasible == m.Feasible_paper).all()))
    j1, i1 = [], []
    for ds, c, feas in zip(D.Dataset, D.Coordinates, D.Feasible):
        if feas:
            o, i = MD.bench_objective(W.parse_coords(c), int(ds), 15)
        else:
            o, i = np.nan, np.nan
        j1.append(o); i1.append(i)
    D["J15"], D["I15"] = j1, i1
    return D.reset_index(drop=True), cons


def arm(D, obj, ideal):
    """Copy of D with Objective / Ideal from columns obj / ideal on the feasible rows (infeasible rows unchanged:
    they never enter a mean)."""
    X = D.copy()
    ok = X.Feasible.values
    if obj != "Objective":
        X.loc[ok, "Objective"] = X.loc[ok, obj].values
        X.loc[ok, "Ideal"] = X.loc[ok, ideal].values
    X["WakeLoss"] = X.Ideal - X.Objective
    X["LossPct"] = 100 * X.WakeLoss / X.Ideal
    return X


# ------------------------------------------------------------------ statistics
def allrun(G, a, b):
    """All-run paired outcome of a vs b over all seed pairs of all cases."""
    rows = []
    for c, g in G.groupby(CASE):
        A_ = g[g.Algorithm == a].set_index("Seed"); B_ = g[g.Algorithm == b].set_index("Seed")
        for s in A_.index.intersection(B_.index):
            fa, fb = bool(A_.Feasible[s]), bool(B_.Feasible[s])
            if fa and fb:
                d = A_.Objective[s] - B_.Objective[s]
                o = 0 if abs(d) <= TIE_OBJ else (1 if d > 0 else -1)
            else:
                o = 0 if fa == fb else (1 if fa else -1)
            rows.append((c, o))
    if not rows:
        return None
    cases = sorted({c for c, _ in rows})
    O = {c: np.array([o for cc, o in rows if cc == c]) for c in cases}
    allo = np.concatenate([O[c] for c in cases])
    W_, T_, L_ = int((allo == 1).sum()), int((allo == 0).sum()), int((allo == -1).sum())
    n = len(allo)
    score = (W_ + T_ / 2) / n
    rng = np.random.default_rng(ALLRUN_SEED)
    wc = np.array([(O[c] == 1).sum() + 0.5 * (O[c] == 0).sum() for c in cases]); nc = np.array([len(O[c]) for c in cases])
    I = rng.integers(0, len(cases), (ALLRUN_B, len(cases)))
    bs = wc[I].sum(1) / nc[I].sum(1)
    p = float(binomtest(W_, W_ + L_, 0.5).pvalue) if W_ + L_ > 0 else 1.0
    return dict(first=a, second=b, W=W_, T=T_, L=L_, n_pairs=n, score=float(score),
                ci95=[float(np.quantile(bs, 0.025)), float(np.quantile(bs, 0.975))], sign_test_p=p,
                n_cases=len(cases), resamples=ALLRUN_B, seed=ALLRUN_SEED)


def best_methods(S):
    """Per case: the set of methods with rank 1 (ties kept)."""
    out = {}
    for c, g in S.groupby(CASE):
        r = g.Rank.min()
        out[c] = set(g[g.Rank == r].Algorithm)
    return out


def arm_block(G, name, seed_level=True):
    S = MR.case_stats(G, METHODS)
    fr = RA.ranking(S, METHODS, FOCUS)
    comp = RA.completeness(G, METHODS, CASE)
    o = dict(label=ARM_LAB[name], n_runs=int(len(G)), completeness=comp,
             seeds=sorted(int(s) for s in G.Seed.unique()), friedman=fr,
             feasible_pct={a: float(100 * G[G.Algorithm == a].Feasible.mean()) for a in METHODS if (G.Algorithm == a).any()},
             n_qualified={a: int(S[(S.Algorithm == a)].Qualified.sum()) for a in METHODS},
             mean_loss_pct_qualified={a: float(S[(S.Algorithm == a) & S.Qualified].Loss.mean()) for a in METHODS},
             per_case={f"{c[0]}-{c[1]}-{c[2]}": {a: dict(loss=(float(v.Loss) if np.isfinite(v.Loss) else None),
                                                          nfeas=int(v.NFeas), rank=float(v.Rank))
                                                  for a, v in g.set_index("Algorithm").iterrows()}
                       for c, g in S.groupby(CASE)})
    if fr and "avg_rank" in fr:
        ar = fr["avg_rank"]; best = min(ar.values())
        o["leader"] = sorted(a for a in METHODS if abs(ar[a] - best) < 1e-12)
        o["rank_order"] = sorted(METHODS, key=lambda a: ar[a])
    o["best_per_case"] = {f"{c[0]}-{c[1]}-{c[2]}": sorted(v) for c, v in best_methods(S).items()}
    o["contrast_psovns_pso"] = RA.contrast(G, S, FOCUS, "PSOC", seed_level=seed_level)
    ar_ = {b: allrun(G, FOCUS, b) for b in METHODS if b != FOCUS}
    keys = [b for b in ar_ if ar_[b]]
    for b, h in zip(keys, MR.holm([ar_[b]["sign_test_p"] for b in keys]) if keys else []):
        ar_[b]["sign_test_p_holm4"] = float(h)
    o["allrun_vs_focus"] = ar_
    return o, S


def gain(Sa, Sb, a):
    """Per case mean loss of method a in arm Sa minus arm Sb (cases where a qualifies in both); case bootstrap."""
    pa = Sa[(Sa.Algorithm == a) & Sa.Qualified].set_index(CASE).Loss
    pb = Sb[(Sb.Algorithm == a) & Sb.Qualified].set_index(CASE).Loss
    j = pa.index.intersection(pb.index)
    d = (pa[j] - pb[j]).values
    if len(d) == 0:
        return None
    return dict(n_cases=int(len(d)), mean_pp=float(d.mean()), ci90=MR.boot_ci(d, level=0.90) if len(d) >= 2 else [None, None],
                lower_in=int((d < -1e-9).sum()), higher_in=int((d > 1e-9).sum()),
                per_case={f"{c[0]}-{c[1]}-{c[2]}": float(v) for c, v in zip(j, d)})


def compare(blocks, SS, x, y):
    """Leader / order comparison of arm x (direct) with arm y (reference)."""
    fx, fy = blocks[x]["friedman"], blocks[y]["friedman"]
    if not (fx and fy and "avg_rank" in fx and "avg_rank" in fy):
        return None
    tau = kendalltau([fx["avg_rank"][a] for a in METHODS], [fy["avg_rank"][a] for a in METHODS])[0]
    bx, by = best_methods(SS[x]), best_methods(SS[y])
    common = [c for c in bx if c in by]
    same = [bool(bx[c] & by[c]) for c in common]
    return dict(arms=[x, y], leader_x=blocks[x].get("leader"), leader_y=blocks[y].get("leader"),
                leader_changes=blocks[x].get("leader") != blocks[y].get("leader"),
                kendall_tau_avg_rank=float(tau), same_best_cases=int(sum(same)), n_cases=len(common),
                same_best_pct=float(100 * np.mean(same)) if same else None,
                gain_x_minus_y={a: gain(SS[x], SS[y], a) for a in METHODS})


def seed_paired_gain(D1, D2, a):
    """Seed-paired difference of 1-deg loss (direct - ctrl1) for method a over jointly feasible runs; per case
    mean, then mean over cases; case-bootstrap 90 % CI."""
    k = CASE + ["Seed"]
    A_ = D1[(D1.Algorithm == a) & D1.Feasible].set_index(k).LossPct
    B_ = D2[(D2.Algorithm == a) & D2.Feasible].set_index(k).LossPct
    j = A_.index.intersection(B_.index)
    if not len(j):
        return None
    dd = (A_[j] - B_[j]).groupby(level=[0, 1, 2]).mean()
    return dict(n_pairs=int(len(j)), n_cases=int(len(dd)), mean_pp=float(dd.mean()),
                ci90=MR.boot_ci(dd.values, level=0.90) if len(dd) >= 2 else [None, None])


# ------------------------------------------------------------------ LaTeX
def fnum(v, d=3, sign=True):
    return RA.fnum(v, d, sign)


def fci(c, d=3):
    return RA.fci(c, d)


def latex(summ, partial):
    B = summ["arms"]
    dr, rv = B.get("direct"), B["reeval"]
    pt = " [PARTIAL DATA -- preview only]" if partial else ""

    def rk(blk, a):
        f = blk["friedman"] if blk else None
        return fnum(f["avg_rank"][a], 2, False) if f and "avg_rank" in f else "--"

    def ml(blk, a):
        if not blk:
            return "--"
        v, nq = blk["mean_loss_pct_qualified"].get(a), blk["n_qualified"].get(a)
        s = fnum(v, 3, False)
        return s + (f"$^{{({nq})}}$" if nq is not None and nq < len(CASES) and s != "--" else "")

    lines = []
    for a in METHODS:
        g = summ["compare"]["direct_vs_reeval"]["gain_x_minus_y"].get(a) if summ.get("compare", {}).get("direct_vs_reeval") else None
        al = dr["allrun_vs_focus"].get(a) if dr and a != FOCUS else None
        cells = [LAB[a],
                 fnum(dr["feasible_pct"].get(a), 1, False) if dr else "--", rk(dr, a), ml(dr, a),
                 "--" if not al else f"{al['score']:.3f} {fci(al['ci95'], 2)}",
                 fnum(rv["feasible_pct"].get(a), 1, False), rk(rv, a), ml(rv, a),
                 "--" if not g else f"{fnum(g['mean_pp'])} {fci(g['ci90'], 2)}"]
        lines.append(" & ".join(cells) + " \\\\")
    lines.append("\\midrule")
    lines.append("\\emph{PSO-VNS $-$ PSO} & \\multicolumn{2}{c}{$n$} & $\\overline{\\Delta L}$ & 90\\% CI seed (Eq.) & "
                 "\\multicolumn{3}{c}{90\\% CI case (Eq.)} & $p_W$ \\\\")
    lines.append("\\midrule")
    for key, lab in (("direct", "Direct 1$^\\circ$"), ("reeval", "Re-eval.\\ 1$^\\circ$"), ("rec15", "Recorded 15$^\\circ$"),
                     ("ctrl1", "Control 1$^\\circ$")):
        blk = B.get(key)
        if not blk:
            continue
        r = blk.get("contrast_psovns_pso")
        if not r or not r["n_cases"]:
            lines.append(f"{lab} & \\multicolumn{{8}}{{l}}{{no case in which both qualify}} \\\\")
            continue
        sl, cl = r["seed_level"], r["case_level"]
        eqs = "" if not sl else (" (yes)" if sl["equivalent"] else " (no)")
        eqc = "" if not cl else (" (yes)" if cl["equivalent"] else " (no)")
        lines.append(f"{lab} & \\multicolumn{{2}}{{c}}{{{r['n_cases']}}} & {fnum(r['mean_dloss_pp'])} & "
                     f"{fci(sl['ci90'] if sl else None)}{eqs} & "
                     f"\\multicolumn{{3}}{{c}}{{{fci(cl['ci90_mean_dloss_pp'] if cl else None)}{eqc}}} & "
                     f"{RA.fp(r['wilcoxon']['p'])} \\\\")
    cmp_ = summ.get("compare", {}).get("direct_vs_reeval")
    fr_d = dr["friedman"] if dr else None
    lead = ("leader (best average rank) with direct optimization: " + ", ".join(LAB[a] for a in dr.get("leader", [])) +
            "; with re-evaluation: " + ", ".join(LAB[a] for a in rv.get("leader", [])) +
            (f"; Kendall $\\tau$ between the two average-rank orders {cmp_['kendall_tau_avg_rank']:.2f}; same best method "
             f"in {cmp_['same_best_cases']} of {cmp_['n_cases']} cases" if cmp_ else "")) if dr and fr_d else ""
    frt = ""
    if fr_d and "chi2" in fr_d:
        frt = f" Friedman (direct arm): $\\chi^2_F={fr_d['chi2']:.1f}$, $p={RA.fp(fr_d['p']).strip('$')}$."
    head = ("& \\multicolumn{4}{c}{Direct 1$^\\circ$ optimization (seeds 31--60)} & "
            "\\multicolumn{3}{c}{Re-eval.\\ 1$^\\circ$ (seeds 1--30)} & Direct $-$ re-eval. \\\\\n"
            "\\cmidrule(lr){2-5}\\cmidrule(lr){6-8}\\cmidrule(lr){9-9}\n"
            "Method & Feas. & Rank & $\\bar L$ & All-run [95\\% CI] & Feas. & Rank & $\\bar L$ & "
            "$\\Delta\\bar L$ [90\\% CI]")
    cap = ("Direct optimization with 1$^\\circ$ direction bins versus 1$^\\circ$ re-evaluation of layouts optimized with "
           "15$^\\circ$ bins: 8 cases (data sets I and II; $r=500$~m, $N=10$; $r=750$~m, $N=6$ and 12; $r=1000$~m, $N=15$), "
           "five methods (PSO-VNS, PSO, GA, MS-SLSQP, RSD-VNS), 30 seed-paired runs per method and case, 6{,}030 "
           "evaluations, random starts. Endpoint: wake loss $L$ (\\% of the wake-free objective) under the 1$^\\circ$ "
           "objective (15 sub-bins per benchmark bin); differences first minus second, negative = first better." + pt)
    note = ("$\\bar L$ in \\%, $\\Delta\\bar L$ in pp. Feas.: feasible runs (\\%). Rank: average feasibility-aware rank over the 8 cases (qualified = at least 15 "
            "of 30 runs feasible). $\\bar L$: mean over the cases in which the method qualifies (superscript: number of "
            "such cases if fewer than 8) of the mean wake loss of its feasible runs. All-run: paired score of PSO-VNS "
            "against the method over 240 seed pairs ($(W+T/2)/240$; feasible beats infeasible, two feasible runs "
            "compared by objective, two infeasible runs tie; above 0.5 = PSO-VNS better), with 95\\% case-bootstrap "
            "interval. Direct $-$ re-eval.: per-case mean $L$ of direct 1$^\\circ$ optimization minus that of the "
            "re-evaluated 15$^\\circ$ layouts (different seed sets), mean over the cases in which the method qualifies "
            "in both, 90\\% case-bootstrap interval. Lower panel: Direct 1$^\\circ$ = direct 1$^\\circ$ optimization (seeds 31--60); Re-eval.\\ 1$^\\circ$ = "
            "15$^\\circ$ layouts re-evaluated at 1$^\\circ$ (seeds 1--30); Recorded 15$^\\circ$ = the same layouts with "
            "their recorded 15$^\\circ$ objective; Control 1$^\\circ$ (if run) = 15$^\\circ$ runs with seeds 31--60 "
            "re-evaluated at 1$^\\circ$; $n$ cases in which both qualify; $\\overline{\\Delta L}$ in pp; seed: within-case "
            "seed-paired bootstrap; case: bootstrap over the cases (10{,}000 resamples, fixed seeds); Eq.: 90\\% "
            "interval inside $(-0.05, 0.05)$~pp; $p_W$: case-mean Wilcoxon (with 8 cases, $p\\ge0.0078$). Exploratory; "
            "" + lead + "." + frt)
    return MR.table("table*", cap, "tab:S-r3-fine", "lcccccccc", head, lines, sep="2pt", pos="!htb", note=note)


# ------------------------------------------------------------------ main
def clean(o):
    if isinstance(o, dict):
        return {str(k): clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple, set)):
        return [clean(v) for v in (sorted(o) if isinstance(o, set) else o)]
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, (np.floating, float)):
        return None if not np.isfinite(o) else float(o)
    if isinstance(o, np.bool_):
        return bool(o)
    return o


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default=HERE)
    ap.add_argument("--out-dir", default=HERE)
    args = ap.parse_args(argv)
    summ = dict(experiment="rev3_fine", manifest="rev3_fine_manifest.md", methods=METHODS,
                cases=[list(c) for c in CASES], margin_pp=EQ, data={}, arms={}, compare={},
                conventions=dict(endpoint="wake loss % of wake-free objective under the 1-deg objective (15 sub-bins)",
                                 sign="first minus second, pp; negative = first better",
                                 allrun=f"feasible beats infeasible; both feasible: higher objective wins, tie if |diff| <= {TIE_OBJ}; "
                                        f"both infeasible tie; 95% case bootstrap, {ALLRUN_B} resamples, seed {ALLRUN_SEED}",
                                 case_boot=dict(resamples=MR.BOOT_N, seed=MR.BOOT_SEED),
                                 seed_boot=dict(resamples=MX.BOOT_N, seed=MX.SEED_BOOT_SEED)))
    Dn, st = read_shards("rev3_fine", args.data_dir)
    summ["data"]["rev3_fine"] = st
    Dc, stc = read_shards("rev3_fine15", args.data_dir)
    summ["data"]["rev3_fine15"] = stc
    partial = False
    R, cons = stored_runs()
    summ["data"]["reference"] = dict(rows=int(len(R)), consistency_with_paper_loader=cons,
                                     feasible_reevaluated=int(R.Feasible.sum()))
    log(f"reference arm: {len(R)} stored runs; consistency with mpce_inference_extra.load: {cons}")
    arms = {"reeval": arm(R, "J15", "I15"), "rec15": arm(R, "Objective", "Ideal")}
    if Dn is not None:
        Dn = in_cases(Dn[Dn.Rose == "1deg"])
        bad = Dn[~Dn.Seed.between(31, 60) | (Dn.Budget != 6030)]
        assert not len(bad), "unexpected rows in rev3_fine"
        # Objective is the 1-deg objective; check against ObjFine (same evaluator)
        summ["data"]["rev3_fine"]["objective_equals_objfine"] = bool((Dn.Objective == Dn.ObjFine).all())
        summ["data"]["rev3_fine"]["calls_all_6030"] = bool((Dn.Calls == 6030).all())
        arms["direct"] = arm(Dn, "Objective", "Ideal")
    if Dc is not None:
        Dc = in_cases(Dc[Dc.Rose == "15deg"])
        arms["ctrl15"] = arm(Dc, "Objective", "Ideal")
        arms["ctrl1"] = arm(Dc, "ObjFine", "IdealFine")
    SS = {}
    for name, G in arms.items():
        try:
            blk, S = arm_block(G, name)
        except Exception as e:                       # partial previews: report and continue
            log(f"[{name}] statistics not computed: {e!r}")
            summ["arms"][name] = dict(error=repr(e), n_runs=int(len(G)))
            partial = True
            continue
        if not blk["completeness"]["complete"]:
            partial = partial or name == "direct"
        summ["arms"][name] = blk; SS[name] = S
        fr = blk["friedman"]
        c = blk["contrast_psovns_pso"]
        log(f"[{name}] runs {len(G)}; ranks " + (", ".join(f"{LAB[a]} {fr['avg_rank'][a]:.2f}" for a in blk.get('rank_order', []))
                                               if fr and 'avg_rank' in fr else "--") +
            f"; leader {blk.get('leader')}; feasible % " + ", ".join(f"{LAB[a]} {v:.1f}" for a, v in blk['feasible_pct'].items()))
        if c and c.get("n_cases"):
            log(f"   PSO-VNS - PSO: n {c['n_cases']} dL {c['mean_dloss_pp']:+.4f} case90 "
                f"{np.round(c['case_level']['ci90_mean_dloss_pp'], 4).tolist() if c['case_level'] else None} seed90 "
                f"{np.round(c['seed_level']['ci90'], 4).tolist() if c.get('seed_level') else c.get('seed_level_error')} "
                f"pW {c['wilcoxon']['p']:.3g}")
        for b, r in blk["allrun_vs_focus"].items():
            if r:
                log(f"   all-run PSO-VNS vs {LAB[b]}: {r['W']}/{r['T']}/{r['L']} score {r['score']:.3f} "
                    f"[{r['ci95'][0]:.3f}, {r['ci95'][1]:.3f}] p_holm {r.get('sign_test_p_holm4', float('nan')):.3g}")
    for x, y in (("direct", "reeval"), ("direct", "rec15"), ("reeval", "rec15"), ("ctrl1", "reeval"), ("direct", "ctrl1")):
        if x in SS and y in SS:
            summ["compare"][f"{x}_vs_{y}"] = compare(summ["arms"], SS, x, y)
            c = summ["compare"][f"{x}_vs_{y}"]
            if c:
                log(f"[compare {x} vs {y}] leaders {c['leader_x']} / {c['leader_y']} (changes: {c['leader_changes']}); "
                    f"tau {c['kendall_tau_avg_rank']:.2f}; same best {c['same_best_cases']}/{c['n_cases']}; gain " +
                    ", ".join(f"{LAB[a]} {g['mean_pp']:+.3f}" for a, g in c["gain_x_minus_y"].items() if g))
    if "direct" in arms and "ctrl1" in arms:
        summ["compare"]["direct_minus_ctrl1_seed_paired"] = {a: seed_paired_gain(arms["direct"], arms["ctrl1"], a) for a in METHODS}
    summ["partial"] = partial
    os.makedirs(args.out_dir, exist_ok=True)
    with open(os.path.join(args.out_dir, "rev3_fine.json"), "w") as f:
        json.dump(clean(summ), f, indent=1)
    if "direct" in summ["arms"] and "friedman" in summ["arms"]["direct"] and "reeval" in summ["arms"]:
        tex = latex(summ, partial)
    else:
        tex = "% rev3_fine_tables.tex: direct-arm statistics not available (no or too few rev3_fine runs)\n"
    with open(os.path.join(args.out_dir, "rev3_fine_tables.tex"), "w") as f:
        f.write("% generated by rev3_fine_analysis.py -- do not edit\n" + tex)
    log(f"wrote {os.path.join(args.out_dir, 'rev3_fine.json')} and rev3_fine_tables.tex (partial={partial})")


if __name__ == "__main__":
    main()
