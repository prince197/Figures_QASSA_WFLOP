"""Revision 3, item C1 (R5): analysis of the constraint-handling study (pre-specified in rev3_constraint_manifest.md).

Usage (from analysis/):  python3 rev3_constraint_analysis.py [--out-dir DIR]
Inputs:  rev3_constraint_{pen,deb,proj}.csv (seeds 31-60), stored seeds 1-30 of the five methods (replication check of
         variant pen; read through rev3_constraint.stored_runs), rev3_constraint_validate.json.
Outputs: rev3_constraint.json, rev3_constraint_tables.tex (tab:S-r3-constraint, tab:S-r3-constraint-cases).

Conventions: wake loss in % of the ideal; Delta L = PSO-VNS minus comparator (pp; negative = PSO-VNS better);
margin 0.05 pp; all-run paired score (W + T/2)/n; exact two-sided sign tests; Holm within each variant's family of
4 comparators x (6 cases + pooled) = 28 tests; two-stage (cases, then seeds) percentile bootstrap, 20,000 resamples,
RNG seed 20261004.
"""
import os, sys, json, argparse, warnings
import numpy as np, pandas as pd
from scipy.stats import binomtest, wilcoxon, rankdata, kendalltau, friedmanchisquare

HERE = os.path.dirname(os.path.abspath(__file__))
warnings.filterwarnings("ignore", message="Mean of empty slice")
sys.path.insert(0, HERE)
import mpce_results as MR
import rev3_constraint as RC

VARS = ["pen", "deb", "proj"]
VLAB = {"pen": "(i) penalty $F_p$, box clipping", "deb": "(ii) Deb rules, normalized $v$, box clipping",
        "proj": "(iii) Deb rules, radial projection"}
VSHORT = {"pen": "(i)", "deb": "(ii)", "proj": "(iii)"}
METH = RC.METHODS
FOCUS = "PSOBV"
OTHERS = [m for m in METH if m != FOCUS]
LAB = {"PSOBV": "PSO-VNS", "PSOC": "PSO", "GA": "GA", "SSABV": "SSA-VNS", "DE": "DE"}
CASES = RC.CASES
CLAB = {c: f"{'I' if c[0] == 1 else 'II'}/{c[1]}/{c[2]}" for c in CASES}
MARGIN = 0.05
B = 20000
SEED = 20261004
OBJ_TIE = 1e-9
ALPHA = 0.05


def load(v, data_dir=HERE):
    d = pd.read_csv(os.path.join(data_dir, f"rev3_constraint_{v}.csv"), float_precision="round_trip",
                    dtype={"Coordinates": str, "Curve": str})
    d["Feasible"] = d.Feasible.astype(bool)
    d["LossPct"] = 100.0 * d.WakeLoss / d.Ideal
    return d


def cube(d, col, seeds):
    """array [method, case, seed] of column col."""
    idx = d.set_index(["Algorithm", "Dataset", "Radius", "Turbines", "Seed"])[col]
    return np.array([[[idx[(m, *c, s)] for s in seeds] for c in CASES] for m in METH])


def outcome(fa, fb, oa, ob, la, lb, margin=None):
    """per-pair outcome of a vs b: 1 win, 0.5 tie, 0 loss (arrays)."""
    o = np.full(fa.shape, 0.5)
    o[fa & ~fb] = 1.0
    o[~fa & fb] = 0.0
    both = fa & fb
    diff = oa - ob                                   # objective, higher is better
    win = both & (diff > OBJ_TIE); loss = both & (diff < -OBJ_TIE)
    if margin is not None:
        dl = la - lb
        win = both & (dl < -margin); loss = both & (dl > margin)
    o[win] = 1.0; o[loss] = 0.0
    return o


class Boot:
    def __init__(self, ncase, nseed, b=B, seed=SEED):
        rng = np.random.default_rng(seed)
        self.ci = rng.integers(0, ncase, (b, ncase))
        self.si = rng.integers(0, nseed, (b, ncase, nseed))
        rng2 = np.random.default_rng(seed + 1)
        self.s1 = rng2.integers(0, nseed, (b, nseed))

    def two_stage(self, x):
        """x: [case, seed] (NaN = missing); percentile CI of the pooled nanmean."""
        xs = x[self.ci[:, :, None], self.si]
        with np.errstate(invalid="ignore"):
            st = np.nanmean(xs.reshape(len(xs), -1), axis=1)
        return st

    def seeds(self, x):
        xs = x[self.s1]
        with np.errstate(invalid="ignore"):
            return np.nanmean(xs, axis=1)


def pci(st, level=0.95):
    st = st[np.isfinite(st)]
    if len(st) == 0:
        return [None, None]
    a = (1 - level) / 2
    return [float(np.quantile(st, a)), float(np.quantile(st, 1 - a))]


def wilson(k, n, z=1.959963984540054):
    if n == 0:
        return [None, None]
    p = k / n; den = 1 + z * z / n
    c = (p + z * z / (2 * n)) / den; h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return [float(c - h), float(c + h)]


def sign_p(w, l):
    return 1.0 if w + l == 0 else float(binomtest(int(w), int(w + l), 0.5).pvalue)


def wil_p(d):
    d = d[np.isfinite(d)]
    d = d[d != 0]
    if len(d) < 1:
        return None
    try:
        return float(wilcoxon(d).pvalue)
    except ValueError:
        return None


def block_ranks(F, O):
    """F, O: [method, case, seed]; run-level feasibility-aware ranks [method, case, seed]."""
    val = np.where(F, -np.round(O, 6), np.inf)
    R = np.empty(val.shape)
    for c in range(val.shape[1]):
        for s in range(val.shape[2]):
            R[:, c, s] = rankdata(val[:, c, s])
    return R


def analyse_variant(d, seeds, boot):
    F = cube(d, "Feasible", seeds).astype(bool)
    O = cube(d, "Objective", seeds).astype(float)
    L = cube(d, "LossPct", seeds).astype(float)
    res = dict(methods={}, comparisons={}, ranks={})
    # 1-2 feasibility and conditional mean wake loss
    for i, m in enumerate(METH):
        nf_case = F[i].sum(1)
        cm = [float(np.mean(L[i, c][F[i, c]])) if F[i, c].any() else None for c in range(len(CASES))]
        valid = [x for x in cm if x is not None]
        dm = d[d.Algorithm == m]
        inf = dm[~dm.Feasible]
        bnd = inf.MaxBoundExcess > RC.TOL; spc = inf.MaxSpacingDeficit > RC.TOL
        res["methods"][m] = dict(
            feasible=int(F[i].sum()), runs=int(F[i].size), feas_pct=100.0 * F[i].mean(),
            feas_wilson95=[100 * x for x in wilson(int(F[i].sum()), F[i].size)],
            feas_case_pct={CLAB[c]: 100.0 * nf_case[k] / F.shape[2] for k, c in enumerate(CASES)},
            cond_loss_case={CLAB[c]: cm[k] for k, c in enumerate(CASES)},
            cond_loss_mean_of_cases=float(np.mean(valid)) if valid else None, cond_loss_cases_used=len(valid),
            infeasible_final=dict(boundary_only=int((bnd & ~spc).sum()), spacing_only=int((~bnd & spc).sum()),
                                  both=int((bnd & spc).sum())),
            feas_eval_pct_mean=float(dm.FeasEvalPct.mean()),
            first_feas_call_median=float(dm.FirstFeasCall[dm.FirstFeasCall > 0].median()) if (dm.FirstFeasCall > 0).any() else None,
            runs_without_feasible_eval=int((dm.FirstFeasCall < 0).sum()),
            projections_mean=float(dm.Projections.mean()), seconds_median=float(dm.Seconds.median()),
            # checkpoint 100 of the curve = evaluation 3,030 (the PSO-VNS / SSA-VNS switch): a feasible layout found by then
            feasible_by_3030_pct=float(100 * np.mean([np.isfinite(float(c.split(";")[100])) for c in dm.Curve])),
            feasible_by_3030_case_pct={CLAB[c]: float(100 * np.mean([np.isfinite(float(x.split(";")[100])) for x in dm[
                (dm.Dataset == c[0]) & (dm.Radius == c[1]) & (dm.Turbines == c[2])].Curve])) for c in CASES})
    # 3-4 paired comparisons
    fi = METH.index(FOCUS)
    tests_s, tests_w = [], []
    for m in OTHERS:
        j = METH.index(m)
        o = outcome(F[fi], F[j], O[fi], O[j], L[fi], L[j])
        om = outcome(F[fi], F[j], O[fi], O[j], L[fi], L[j], margin=MARGIN)
        dl = np.where(F[fi] & F[j], L[fi] - L[j], np.nan)
        comp = dict(cases={})
        for k, c in enumerate(CASES):
            w, t, l = int((o[k] == 1).sum()), int((o[k] == 0.5).sum()), int((o[k] == 0).sum())
            dk = dl[k]
            comp["cases"][CLAB[c]] = dict(
                W=w, T=t, L=l, score=float(o[k].mean()), score_ci95=pci(boot.seeds(o[k])), sign_p=sign_p(w, l),
                score_margin=float(om[k].mean()), both_feasible=int(np.isfinite(dk).sum()),
                dL_mean=float(np.nanmean(dk)) if np.isfinite(dk).any() else None,
                dL_ci95=pci(boot.seeds(dk)) if np.isfinite(dk).any() else [None, None],
                dL_ci90=pci(boot.seeds(dk), 0.90) if np.isfinite(dk).any() else [None, None],
                wilcoxon_p=wil_p(dk))
            tests_s.append((m, CLAB[c], comp["cases"][CLAB[c]]["sign_p"]))
            tests_w.append((m, CLAB[c], comp["cases"][CLAB[c]]["wilcoxon_p"]))
        w, t, l = int((o == 1).sum()), int((o == 0.5).sum()), int((o == 0).sum())
        st = boot.two_stage(o); stm = boot.two_stage(om); sd = boot.two_stage(dl)
        comp["pooled"] = dict(W=w, T=t, L=l, score=float(o.mean()), score_ci95=pci(st), sign_p=sign_p(w, l),
                              score_margin=float(om.mean()), score_margin_ci95=pci(stm),
                              Wm=int((om == 1).sum()), Tm=int((om == 0.5).sum()), Lm=int((om == 0).sum()),
                              both_feasible=int(np.isfinite(dl).sum()), dL_mean=float(np.nanmean(dl)),
                              dL_ci95=pci(sd), dL_ci90=pci(sd, 0.90), wilcoxon_p=wil_p(dl.ravel()))
        tests_s.append((m, "pooled", comp["pooled"]["sign_p"]))
        tests_w.append((m, "pooled", comp["pooled"]["wilcoxon_p"]))
        res["comparisons"][m] = comp
    for tests, key in ((tests_s, "sign_p_holm"), (tests_w, "wilcoxon_p_holm")):
        idx = [k for k, t in enumerate(tests) if t[2] is not None]
        adj = MR.holm([tests[k][2] for k in idx])
        for k, a in zip(idx, adj):
            m, cl, _ = tests[k]
            tgt = res["comparisons"][m]["pooled"] if cl == "pooled" else res["comparisons"][m]["cases"][cl]
            tgt[key] = float(a)
        res[f"family_size_{key}"] = len(idx)
    for m in OTHERS:
        for tgt in [res["comparisons"][m]["pooled"]] + list(res["comparisons"][m]["cases"].values()):
            sig = tgt["sign_p_holm"] < ALPHA
            tgt["category"] = ("PSO-VNS better" if tgt["score"] > 0.5 else "PSO-VNS worse") if sig else "ns"
            lo, hi = tgt["dL_ci90"]; lo95, hi95 = tgt["dL_ci95"]
            if lo is None:
                tgt["margin_verdict"] = "n/a"
            elif hi95 < -MARGIN:
                tgt["margin_verdict"] = "better beyond margin"
            elif lo >= -MARGIN and hi <= MARGIN:
                tgt["margin_verdict"] = "equivalent"
            elif lo95 > MARGIN:
                tgt["margin_verdict"] = "worse beyond margin"
            else:
                tgt["margin_verdict"] = "inconclusive"
    # 5 ranks
    R = block_ranks(F, O)
    rb = {}
    for i, m in enumerate(METH):
        rb[m] = boot.two_stage(R[i])
        res["ranks"][m] = dict(run_level=float(R[i].mean()), run_level_ci95=pci(rb[m]))
    case_r = np.zeros((len(METH), len(CASES)))
    for k, c in enumerate(CASES):
        fmean = [O[i, k][F[i, k]].mean() if F[i, k].any() else np.nan for i in range(len(METH))]
        case_r[:, k] = MR.rank_rule(fmean, F[:, k].sum(1), [F.shape[2]] * len(METH))
    for i, m in enumerate(METH):
        res["ranks"][m]["case_rule"] = float(case_r[i].mean())
        res["ranks"][m]["case_rule_by_case"] = {CLAB[c]: float(case_r[i, k]) for k, c in enumerate(CASES)}
    order = sorted(METH, key=lambda m: res["ranks"][m]["run_level"])
    res["order_run_level"] = order
    res["order_case_rule"] = sorted(METH, key=lambda m: (res["ranks"][m]["case_rule"], res["ranks"][m]["run_level"]))
    pair_diff = {}
    for a in range(len(order) - 1):
        m1, m2 = order[a], order[a + 1]
        pair_diff[f"{m1}-{m2}"] = dict(diff=float(R[METH.index(m1)].mean() - R[METH.index(m2)].mean()),
                                       ci95=pci(rb[m1] - rb[m2]))
    res["adjacent_rank_differences"] = pair_diff
    blocks = R.reshape(len(METH), -1)
    res["friedman_run_level_p"] = float(friedmanchisquare(*blocks).pvalue)
    return res, dict(F=F, O=O, L=L, R=R, rb=rb)


def compare_orders(A, Bv, ra, rbv):
    """Rank-order change from variant a to variant b."""
    oa, ob = A["order_run_level"], Bv["order_run_level"]
    va = [A["ranks"][m]["run_level"] for m in METH]; vb = [Bv["ranks"][m]["run_level"] for m in METH]
    swaps = []
    for x in range(len(METH)):
        for y in range(x + 1, len(METH)):
            m1, m2 = METH[x], METH[y]
            if (oa.index(m1) < oa.index(m2)) != (ob.index(m1) < ob.index(m2)):
                def ovl(R, mm1, mm2):
                    c1, c2 = R["ranks"][mm1]["run_level_ci95"], R["ranks"][mm2]["run_level_ci95"]
                    return bool(c1[0] <= c2[1] and c2[0] <= c1[1])
                swaps.append(dict(pair=[m1, m2], overlap_in_reference=ovl(A, m1, m2), overlap_in_variant=ovl(Bv, m1, m2)))
    ca, cb = A["order_case_rule"], Bv["order_case_rule"]
    return dict(identical=oa == ob, order_ref=oa, order_var=ob, kendall_tau=float(kendalltau(va, vb).statistic),
                swaps=swaps, case_rule_identical=ca == cb, case_rule_ref=ca, case_rule_var=cb)


def cross_variant(arr, va, vb, boot):
    out = {}
    for i, m in enumerate(METH):
        Fa, Fb = arr[va]["F"][i], arr[vb]["F"][i]
        b_ = int((Fa & ~Fb).sum()); c_ = int((~Fa & Fb).sum())
        dl = np.where(Fa & Fb, arr[vb]["L"][i] - arr[va]["L"][i], np.nan)
        out[m] = dict(feas_ref=int(Fa.sum()), feas_var=int(Fb.sum()), only_ref=b_, only_var=c_,
                      mcnemar_p=sign_p(b_, c_), both_feasible=int(np.isfinite(dl).sum()),
                      dL_var_minus_ref=float(np.nanmean(dl)) if np.isfinite(dl).any() else None,
                      dL_ci95=pci(boot.two_stage(dl)) if np.isfinite(dl).any() else [None, None])
    return out


def replication(dpen):
    S = RC.stored_runs()
    S = S[S.Seed.between(1, 30) & S.Dataset.isin(["1", "2"])].copy()
    S["Dataset"] = S.Dataset.astype(int)
    S["Feasible"] = S.Feasible.astype(str).str.lower().isin(["true", "1"])
    S["LossPct"] = 100.0 * S.WakeLoss / S.Ideal
    out = {}
    for m in METH:
        rows = {}
        for tag, d in (("stored_seeds_1_30", S), ("new_seeds_31_60", dpen)):
            dm = d[(d.Algorithm == m)]
            dm = dm[[ (int(a), int(b), int(c)) in CASES for a, b, c in zip(dm.Dataset, dm.Radius, dm.Turbines)]]
            cms = [dm[(dm.Dataset == c[0]) & (dm.Radius == c[1]) & (dm.Turbines == c[2]) & dm.Feasible].LossPct.mean()
                   for c in CASES]
            cms = [x for x in cms if np.isfinite(x)]
            rows[tag] = dict(runs=int(len(dm)), feas_pct=float(100 * dm.Feasible.mean()),
                             cond_loss_mean_of_cases=float(np.mean(cms)) if cms else None, cases_used=len(cms))
        out[m] = rows
    return out


# ------------------------------------------------------------------------------------------------ LaTeX
def f3(x):
    return "--" if x is None else f"{x:.3f}"


def fci(c, d=3):
    return "--" if c[0] is None else f"[{c[0]:.{d}f}, {c[1]:.{d}f}]"


def fp(p):
    if p is None:
        return "--"
    return "$<$0.001" if p < 0.001 else f"{p:.3f}"


def latex(R, comp):
    L = []
    L.append(r"% rev3_constraint_tables.tex -- written by analysis/rev3_constraint_analysis.py (do not edit by hand)")
    L.append(r"\begin{table}[!htbp]")
    L.append(r"\centering")
    L.append(r"\caption{Alternative constraint handling (item C1): five methods under (i) the implemented penalty $F_p$ "
             r"($\mu=10^{10}$, $g^{\rm b}$ in m$^2$, $g^{\rm s}$ in m) with box clipping, (ii) Deb's feasibility rules on "
             r"dimensionless violations ($g^{\rm b}/r^2$, $g^{\rm s}/\ell_{\min}$) with box clipping, and (iii) Deb's rules "
             r"with radial projection onto the circle; six cases (data sets I and II; $r=500$~m, $N=10$; $r=750$~m, $N=6$; "
             r"$r=1000$~m, $N=15$), new seeds 31--60 (identical in all variants), 6,030 evaluations, random starts, $N_p=30$. "
             r"Feas.: feasible final layouts (\%, 180 runs). Loss: mean wake loss of the feasible runs (\% of ideal; mean of "
             r"the case means). W/T/L and Score $=(W+T/2)/180$: all-run paired outcome of PSO-VNS against the method "
             r"(feasible beats infeasible; two feasible runs by wake loss; two infeasible runs tie), 95\% two-stage "
             r"bootstrap interval over cases and seeds; $p_{\rm H}$: exact sign test, Holm-adjusted over the 28 tests of "
             r"each variant (4 methods $\times$ 6 cases and pooled). $\Delta L$: PSO-VNS minus method (pp; negative = "
             r"PSO-VNS better) over seed pairs with both runs feasible, 95\% interval. Rank: feasibility-aware average "
             r"rank over the 180 (case, seed) blocks (1 = best; infeasible runs tied below feasible ones) and, in "
             r"parentheses, the paper's case-level rule.}")
    L.append(r"\label{tab:S-r3-constraint}")
    L.append(r"\scriptsize\setlength{\tabcolsep}{2.6pt}")
    L.append(r"\begin{tabular}{@{}llrrcccrc@{}}")
    L.append(r"\toprule")
    L.append(r"Variant & Method & Feas. & Loss & W/T/L & Score [95\% CI] & $p_{\rm H}$ & $\Delta L$ [95\% CI] & Rank \\")
    L.append(r"\midrule")
    for vi, v in enumerate(VARS):
        if v not in R:
            continue
        r = R[v]
        for k, m in enumerate(METH):
            mm = r["methods"][m]
            loss = f3(mm["cond_loss_mean_of_cases"])
            if mm["cond_loss_cases_used"] < len(CASES):
                loss += f"$^{{{mm['cond_loss_cases_used']}}}$"
            rk = f"{r['ranks'][m]['run_level']:.2f} ({r['ranks'][m]['case_rule']:.2f})"
            if m == FOCUS:
                cells = ["--", "--", "--", "--"]
            else:
                p = r["comparisons"][m]["pooled"]
                cells = [f"{p['W']}/{p['T']}/{p['L']}", f"{p['score']:.3f} {fci(p['score_ci95'])}", fp(p["sign_p_holm"]),
                         f"{p['dL_mean']:+.3f} {fci(p['dL_ci95'])}" if p["both_feasible"] else "--"]
            first = VSHORT[v] if k == 0 else ""
            L.append(f"{first} & {LAB[m]} & {mm['feas_pct']:.1f} & {loss} & " + " & ".join(cells) + f" & {rk} \\\\")
        if vi < len(VARS) - 1:
            L.append(r"\midrule")
    L.append(r"\bottomrule")
    L.append(r"\end{tabular}")
    L.append(r"\par\smallskip\parbox{\textwidth}{\scriptsize Superscript at Loss: number of cases (of 6) with at least "
             r"one feasible run. Run-level rank order: " + "; ".join(
                 f"{VSHORT[v]} " + " $<$ ".join(LAB[m] for m in R[v]["order_run_level"]) for v in VARS if v in R) + ".}")
    L.append(r"\end{table}")
    L.append("")
    # per-case table
    L.append(r"\begin{table}[!htbp]")
    L.append(r"\centering")
    L.append(r"\caption{Alternative constraint handling (item C1), per case: feasible runs (of 30) and all-run paired score "
             r"of PSO-VNS against each method, $(W+T/2)/30$, for the variants (i) penalty with box clipping, (ii) Deb's rules "
             r"with box clipping and (iii) Deb's rules with radial projection; seeds 31--60, 6,030 evaluations, random "
             r"starts. $^{*}$: Holm-adjusted exact sign test $p<0.05$ (family of 28 tests per variant). Case: data set / "
             r"$r$ (m) / $N$.}")
    L.append(r"\label{tab:S-r3-constraint-cases}")
    L.append(r"\scriptsize\setlength{\tabcolsep}{2.4pt}")
    L.append(r"\begin{tabular}{@{}l" + "ccc" * len(METH) + "@{}}")
    L.append(r"\toprule")
    L.append(" & " + " & ".join(rf"\multicolumn{{3}}{{c}}{{{LAB[m]}}}" for m in METH) + r" \\")
    L.append("".join(rf"\cmidrule(lr){{{2 + 3 * k}-{4 + 3 * k}}}" for k in range(len(METH))))
    L.append("Case & " + " & ".join("(i) & (ii) & (iii)" for _ in METH) + r" \\")
    L.append(r"\midrule")
    L.append(r"\multicolumn{" + str(1 + 3 * len(METH)) + r"}{@{}l}{\emph{Feasible runs (of 30)}} \\")
    for c in CASES:
        cells = []
        for m in METH:
            for v in VARS:
                cells.append(str(int(round(R[v]["methods"][m]["feas_case_pct"][CLAB[c]] * 30 / 100))) if v in R else "--")
        L.append(f"{CLAB[c]} & " + " & ".join(cells) + r" \\")
    L.append(r"\midrule")
    L.append(r"\multicolumn{" + str(1 + 3 * len(METH)) + r"}{@{}l}{\emph{Score of PSO-VNS against the method}} \\")
    for c in CASES:
        cells = []
        for m in METH:
            for v in VARS:
                if m == FOCUS or v not in R:
                    cells.append("--")
                else:
                    q = R[v]["comparisons"][m]["cases"][CLAB[c]]
                    cells.append(f"{q['score']:.2f}" + ("$^{*}$" if q["sign_p_holm"] < ALPHA else ""))
        L.append(f"{CLAB[c]} & " + " & ".join(cells) + r" \\")
    L.append(r"\bottomrule")
    L.append(r"\end{tabular}")
    L.append(r"\end{table}")
    return "\n".join(L) + "\n"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", default=HERE)
    ap.add_argument("--data-dir", default=HERE, help="where rev3_constraint_<variant>.csv are read (default analysis/)")
    a = ap.parse_args()
    seeds = list(RC.SEEDS)
    boot = Boot(len(CASES), len(seeds))
    data, R, arr = {}, {}, {}
    for v in VARS:
        fn = os.path.join(a.data_dir, f"rev3_constraint_{v}.csv")
        if not os.path.exists(fn):
            print("missing", fn); continue
        data[v] = load(v, a.data_dir)
        assert len(data[v]) == len(METH) * len(CASES) * len(seeds), (v, len(data[v]))
        assert (data[v].Calls == RC.BUDGET).all(), v
        R[v], arr[v] = analyse_variant(data[v], seeds, boot)
    out = dict(manifest="rev3_constraint_manifest.md", seeds=[seeds[0], seeds[-1]], budget=RC.BUDGET,
               cases=[CLAB[c] for c in CASES], methods=METH, margin_pp=MARGIN, bootstrap=dict(B=B, rng_seed=SEED),
               variants=R)
    try:
        with open(os.path.join(HERE, "rev3_constraint_validate.json")) as fh:
            out["validation"] = json.load(fh)
    except FileNotFoundError:
        out["validation"] = None
    out["order_change"] = {}
    out["cross_variant"] = {}
    out["survival"] = {}
    if "pen" in R:
        for v in VARS[1:]:
            if v not in R:
                continue
            oc = compare_orders(R["pen"], R[v], arr["pen"], arr[v])
            out["order_change"][f"pen->{v}"] = oc
            out["cross_variant"][f"{v}_vs_pen"] = cross_variant(arr, "pen", v, boot)
            cats = {m: (R["pen"]["comparisons"][m]["pooled"]["category"], R[v]["comparisons"][m]["pooled"]["category"])
                    for m in OTHERS}
            case_changes = [dict(method=m, case=cl, ref=R["pen"]["comparisons"][m]["cases"][cl]["category"],
                                 var=R[v]["comparisons"][m]["cases"][cl]["category"])
                            for m in OTHERS for cl in R[v]["comparisons"][m]["cases"]
                            if R["pen"]["comparisons"][m]["cases"][cl]["category"] != R[v]["comparisons"][m]["cases"][cl]["category"]]
            unresolved_only = all(s["overlap_in_reference"] or s["overlap_in_variant"] for s in oc["swaps"])
            out["survival"][v] = dict(pooled_categories=cats,
                                      pooled_categories_unchanged=all(x == y for x, y in cats.values()),
                                      rank_order_unchanged=oc["identical"],
                                      swaps_only_unresolved=unresolved_only if oc["swaps"] else True,
                                      case_category_changes=case_changes,
                                      survives=all(x == y for x, y in cats.values()) and (oc["identical"] or unresolved_only))
        out["replication_pen"] = replication(data["pen"])
    if "deb" in R and "proj" in R:
        out["cross_variant"]["proj_vs_deb"] = cross_variant(arr, "deb", "proj", boot)
        out["order_change"]["deb->proj"] = compare_orders(R["deb"], R["proj"], arr["deb"], arr["proj"])
    with open(os.path.join(a.out_dir, "rev3_constraint.json"), "w") as fh:
        json.dump(out, fh, indent=1, default=float)
    with open(os.path.join(a.out_dir, "rev3_constraint_tables.tex"), "w") as fh:
        fh.write(latex(R, out))
    # console summary
    for v in R:
        print(f"== {v}: order {R[v]['order_run_level']}  case-rule {R[v]['order_case_rule']}")
        for m in METH:
            mm = R[v]["methods"][m]; rk = R[v]["ranks"][m]
            line = f"  {m:6s} feas {mm['feas_pct']:5.1f}  loss {f3(mm['cond_loss_mean_of_cases'])} ({mm['cond_loss_cases_used']})  rank {rk['run_level']:.2f} {fci(rk['run_level_ci95'], 2)} case {rk['case_rule']:.2f}"
            if m != FOCUS:
                p = R[v]["comparisons"][m]["pooled"]
                line += f"  W/T/L {p['W']}/{p['T']}/{p['L']} score {p['score']:.3f} {fci(p['score_ci95'])} pH {fp(p['sign_p_holm'])} {p['category']}; dL {p['dL_mean']:+.3f} {fci(p['dL_ci95'])} ({p['margin_verdict']})"
            print(line)
    print(json.dumps(out.get("survival", {}), indent=1, default=float))


if __name__ == "__main__":
    main()
