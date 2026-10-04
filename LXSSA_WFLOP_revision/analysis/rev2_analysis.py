"""Review round 2: statistics and supplementary tables of the new experiments (rev2_*.csv).

Usage:  python3 rev2_analysis.py [--data-dir D] [--out-dir O] [--only BLOCK[,BLOCK...]]
        BLOCK in {spacing, laplace, ga, site, grad}; default: all blocks whose input CSVs exist.

Inputs (read from --data-dir, default this folder; every block is skipped with a note if its CSVs are missing;
shards rev2_<EXP>_s<i>of<k>.csv are merged, a complete shard set is preferred, an incomplete one is used and marked
PARTIAL in the summary):
  spacing   rev2_spacing  (9 methods, 12 split cases, 5D / 6D) + the stored 4D main runs of the same 9 methods
            (fresh_grid / fresh_vgrid / fresh_bgrid / mpce_psoc / mpce_psobv / mpce_slsqp / mpce_rsdisc, loaded with
            mpce_inference_extra.load, i.e. exactly the per-run data of the paper)
  laplace   rev2_laplace  (LXNR, LXREP, LXU, LXC025, LXC05, LXC2; 12 split cases) + stored LXSSA / SSA runs
  ga        rev2_ga (GA, 68 cases) + stored runs of the 8 main methods; rev2_gahr (GA, Horns Rev 16) + mpce_hrfix;
            rev2_gaiea (GA, IEA37 16 / 36) + mpce_iea16/36 (+ PSO-VNS arms mpce_iea16p/36p)
  site      rev2_lg16 (6,030) and rev2_lg16b (30,030): Lillgrund 16-turbine block, 9 methods
  grad      rev2_grad (SLSQPX, PSOSLSQPX on IEA37 16 / 36) + mpce_iea16/36(p)

Statistical protocol: the functions of mpce_results.py (MR) and mpce_inference_extra.py (MX) are imported and reused
unchanged (neither file is modified):
  * wake loss (%) = LossPct = 100 WakeLoss / Ideal (Ideal = wake-free benchmark objective; for the sites: AEP loss
    relative to the wake-free AEP); Delta L = difference of per-case mean wake losses of the feasible runs, in pp
    (MR.std_cols, MR.case_stats, MR.case_mean_diffs; sign: first minus second, negative = first better)
  * feasibility-aware ranking, qualification = at least half of the runs feasible (15 of 30) (MR.rank_rule via
    MR.case_stats), rank matrix, Friedman + Iman-Davenport + Holm z tests (MR.friedman_block)
  * case-level Wilcoxon on the case means with imputation of maximal differences where only one method qualifies
    (MR.case_mean_wilcoxon for Holm families, MX.cm_test for single contrasts)
  * run-level W/T/L: seed-paired Wilcoxon per case, infeasible runs below every feasible run, Holm over the
    comparisons of the case (MR.paired_vs)
  * 90 % CIs and TOST at m = MR.EQ_MARGIN = 0.05 pp: case level = bootstrap over cases (MR.tost), seed level =
    within-case seed-paired bootstrap (MX.seed_level, stratified version)
  * Hodges-Lehmann estimate with exact Wilcoxon-inversion CI (MX.hodges_lehmann)
Seed pairing: within every group (case x budget x spacing) only seeds present for every method of the pool are kept
(no-op on complete data; recorded as dropped_unpaired_runs).

Outputs (--out-dir): rev2_summary.json, rev2_tables.tex (tab:R2-spacing-rank, tab:R2-spacing-contrast,
tab:R2-laplace, tab:R2-ga, tab:R2-ga-hr, tab:R2-ga-iea, tab:R2-site, tab:R2-grad, tab:R2-grad-calls).
"""
import os, sys, glob, re, json, time, argparse, warnings
import numpy as np, pandas as pd
from scipy.stats import friedmanchisquare

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import mpce_results as MR               # noqa: E402  (unchanged; only functions and constants are used)
import mpce_inference_extra as MX       # noqa: E402

warnings.filterwarnings("ignore", category=RuntimeWarning)
EQ = MR.EQ_MARGIN
assert abs(EQ - MX.EQ_MARGIN) < 1e-12, (EQ, MX.EQ_MARGIN)
CASE = MR.CASE
SPLITCASES = MR.SPLITCASES
MAIN8 = list(MR.MAIN8)
SP9 = ["PSOBV", "PSOC", "SSABV", "SSA", "LXSSA", "DE", "BVNS", "SLSQP", "RSDVNS"]      # spacing / Lillgrund pool
LAPV = ["LXNR", "LXREP", "LXU", "LXC025", "LXC05", "LXC2"]
LAPPOOL = ["SSA", "LXSSA"] + LAPV
CALLS_PER_ITER = {"SSA": (30, 200), "LXSSA": (60, 100), "LXNR": (45, 134), "LXREP": (30, 200), "LXU": (60, 100),
                  "LXC025": (60, 100), "LXC05": (60, 100), "LXC2": (60, 100)}    # (calls / iteration, T) at 6,030 (rev2_laplace.py)
GRADM = ["SLSQPX", "PSOSLSQPX"]
HR10 = [a for a in MR.HR_ORDER]                                                   # the ten methods of tab:hr-site
EXCL = [("1", 500, 10), ("2", 500, 10)]                                           # no feasible layout at 5D / 6D
SPACINGS = ["4D", "5D", "6D"]
CONTRASTS_SP = [("PSOBV", "PSOC"), ("SSABV", "RSDVNS")]

LAB = dict(MR.LAB)
LAB.update(GA="GA", LXNR="LX-SSA, no re-evaluation", LXREP="LX-SSA, Laplace replaces midpoint",
           LXU="LX-SSA, $\\gamma\\sim U(-2,2)$", LXC025="LX-SSA, $\\chi=0.25$", LXC05="LX-SSA, $\\chi=0.5$",
           LXC2="LX-SSA, $\\chi=2$", SLSQPX="MS-SLSQP, exact gradient", PSOSLSQPX="PSO-SLSQP, exact gradient")
LAB_TXT = {k: v.replace("$\\gamma\\sim U(-2,2)$", "gamma~U(-2,2)").replace("$\\chi=", "chi=").replace("$", "")
           .replace("\\%", "%") for k, v in LAB.items()}
NOTES = []


def log(msg=""):
    print(msg, flush=True)


# ------------------------------------------------------------------ data
def read_rev2(exp, data_dir):
    """Merge rev2_<exp>_s<i>of<k>.csv (Coordinates / Curve dropped). Returns (DataFrame or None, status)."""
    files = glob.glob(os.path.join(data_dir, f"rev2_{exp}_s*of*.csv"))
    by_k = {}
    for fn in files:
        m = re.search(rf"rev2_{exp}_s(\d+)of(\d+)\.csv$", os.path.basename(fn))
        if m:
            by_k.setdefault(int(m.group(2)), {})[int(m.group(1))] = fn
    if not by_k:
        return None, dict(status="missing", rows=0)
    complete = {k: v for k, v in by_k.items() if set(v) == set(range(k))}
    if complete:
        k = max(complete, key=lambda kk: sum(os.path.getsize(f) for f in complete[kk].values()))
        use, status = list(complete[k].values()), "complete"
    else:
        use, status = [f for v in by_k.values() for f in v.values()], "PARTIAL (incomplete shard set)"
    parts = []
    for f in sorted(use):
        try:
            parts.append(pd.read_csv(f, usecols=lambda c: c not in ("Coordinates", "Curve")))
        except pd.errors.EmptyDataError:
            pass
    if not parts:
        return None, dict(status="empty", rows=0, files=[os.path.basename(f) for f in use])
    df = MR.std_cols(pd.concat(parts, ignore_index=True))
    key = MR.KEY + (["Spacing"] if "Spacing" in df else [])
    df = df.drop_duplicates(key, keep="last").reset_index(drop=True)
    df["Source"] = f"rev2_{exp}"
    return df, dict(status=status, rows=int(len(df)), files=[os.path.basename(f) for f in sorted(use)],
                    methods=sorted(df.Algorithm.unique()))


def read_mpce(exp):
    """Stored runs of the paper (mpce_<exp>_s*of*.csv in this folder) without Coordinates / Curve."""
    files = glob.glob(os.path.join(HERE, f"mpce_{exp}_s*of*.csv"))
    if not files:
        return None
    df = pd.concat([pd.read_csv(f, usecols=lambda c: c not in ("Coordinates", "Curve")) for f in sorted(files)],
                   ignore_index=True)
    return MR.std_cols(df).drop_duplicates(MR.KEY, keep="last").reset_index(drop=True)


_MAIN = {}


def main_runs():
    """Per-run data of the paper, 68 cases, 6,030 evaluations, random initialization (mpce_inference_extra.load)."""
    if "A" not in _MAIN:
        _MAIN["A"] = MX.load(HERE)
    return _MAIN["A"]


def pair_seeds(D, grp, methods):
    """Keep, within every group, only the seeds that every method of `methods` has (no-op on complete data)."""
    D = D[D.Algorithm.isin(methods)]
    keep = []
    for _, g in D.groupby(grp):
        sets = [set(g[g.Algorithm == a].Seed) for a in methods if (g.Algorithm == a).any()]
        common = set.intersection(*sets) if sets and len(sets) == len([a for a in methods if (D.Algorithm == a).any()]) else set()
        keep.append(g[g.Seed.isin(common)])
    out = pd.concat(keep, ignore_index=True) if keep else D.iloc[:0]
    return out, int(len(D) - len(out))


def completeness(D, methods, grp, runs=30):
    c = D.groupby(grp + ["Algorithm"]).size()
    short = {str(k): int(v) for k, v in c.items() if v < runs}
    miss = [a for a in methods if a not in set(D.Algorithm)]
    return dict(complete=not short and not miss, n_groups_with_fewer_runs=len(short),
                groups_with_fewer_runs_examples=dict(list(short.items())[:5]), methods_missing=miss,
                runs_per_group_min=int(c.min()) if len(c) else 0)


def ptag(*comp):
    """Caption suffix for incomplete data (preview runs)."""
    return "" if all(c["complete"] for c in comp) else " [PARTIAL DATA -- preview only]"


# ------------------------------------------------------------------ statistics (reused functions only)
def contrast(G, S, a, b, seed_level=True):
    """a - b (pp of wake loss; negative = a better) over the cases in which both qualify: case-level bootstrap 90 %
    CI + TOST (MR.tost), seed-level bootstrap 90 % CI + TOST (MX.seed_level), case-mean Wilcoxon with imputation
    (MX.cm_test), Hodges-Lehmann (MX.hodges_lehmann, 95 % and 90 % exact CIs)."""
    if a not in set(S.Algorithm) or b not in set(S.Algorithm):
        return None
    d, cnt = MR.case_mean_diffs(S, a, b)
    r = dict(first=a, second=b, sign="first minus second, pp; negative = first better", **cnt)
    r["n_cases"] = int(len(d))
    r["mean_dloss_pp"] = float(d.mean()) if len(d) else None
    r["case_level"] = MR.tost(d.values) if len(d) >= 2 else None
    cw = MX.cm_test(S, a, b)
    r["wilcoxon"] = dict(p=cw["p"], rb=cw["rb"], n_cases=cw["n_cases"], wins_first=cw["wins"], wins_second=cw["losses"],
                         imputed=cw["imputed"], dropped=cw["dropped"],
                         note="case-mean Wilcoxon incl. imputed cases (only one method qualifies); rb > 0 = first better")
    r["seed_level"] = None
    if seed_level and len(d) >= 1:
        try:
            sl, _ = MX.seed_level(G, a, b, d)
            r["seed_level"] = {k: v for k, v in sl.items()}
        except AssertionError as e:                    # fewer than 30 seed-paired runs in some case (partial data)
            r["seed_level_error"] = f"not computed (needs 30 seed-paired runs per case): {e}"
    if len(d) >= 1:
        r["hodges_lehmann"] = MX.hodges_lehmann(d.values)
        r["hodges_lehmann_90"] = MX.hodges_lehmann(d.values, alpha=0.10)
    return r


def wtl_block(G, focus, others):
    """Run-level W/T/L of focus vs each other method (MR.paired_vs per case; Holm over the comparisons of a case)."""
    rows = []
    for (ds, r, n), sub in G.groupby(CASE):
        for x in MR.paired_vs(sub, focus, [b for b in others if b in set(sub.Algorithm)]):
            rows.append(dict(Dataset=ds, Radius=r, Turbines=n, **x))
    C = pd.DataFrame(rows)
    out = {}
    for b in others:
        x = C[C.Baseline == b] if len(C) else C
        vc = x.Outcome.value_counts() if len(x) else pd.Series(dtype=int)
        out[b] = dict(W=int(vc.get("W", 0)), T=int(vc.get("T", 0)), L=int(vc.get("L", 0)),
                      median_rb=float(x.RB.median()) if len(x) else None)
    return out


def ranking(S, methods, focus, cases=None):
    R = MR.rank_matrix(S, methods)
    if cases is not None:
        R = R[[c in cases for c in R.index]]
    if len(R) < 2:
        return None
    try:
        return MR.friedman_block(R, methods, focus)
    except ValueError as e:                                    # e.g. every case completely tied
        return dict(error=str(e), avg_rank={a: float(R[a].mean()) for a in methods}, n_cases=int(len(R)))


def wtl_s(c):
    return "--" if c is None else f"{c['W']}/{c['T']}/{c['L']}"


# ------------------------------------------------------------------ LaTeX formatting
def fnum(v, d=3, sign=True):
    if v is None or not np.isfinite(v):
        return "--"
    s = f"{v:+.{d}f}" if sign else f"{v:.{d}f}"
    return s.replace("-", "$-$")


def fci(c, d=3):
    if c is None or c[0] is None:
        return "--"
    return f"[{fnum(c[0], d, False)}, {fnum(c[1], d, False)}]"


def fbig(v, d=1):
    return "--" if v is None or not np.isfinite(v) else f"{v:,.{d}f}".replace(",", "{,}")


def fp(p):
    return "--" if p is None else MR.fmt_p(p)


def eqv(r):
    return "--" if not r else ("yes" if r["equivalent"] else "no")


# ------------------------------------------------------------------ T-spacing
def block_spacing(data_dir, summ, tex):
    N, st = read_rev2("spacing", data_dir)
    summ["data"]["rev2_spacing"] = st
    if N is None:
        log("[spacing] rev2_spacing missing -> skipped"); return
    A = main_runs()
    M4 = A[A.Algorithm.isin(SP9) & pd.Series([c in SPLITCASES for c in zip(A.Dataset, A.Radius, A.Turbines)], index=A.index)].copy()
    M4["Spacing"] = "4D"
    D = pd.concat([M4, N[N.Algorithm.isin(SP9)]], ignore_index=True)
    D, drop = pair_seeds(D, CASE + ["Spacing"], SP9)
    out = dict(methods=SP9, cases=[list(c) for c in SPLITCASES], excluded_cells=[list(c) for c in EXCL],
               excluded_note="(DS 1 and 2, 500 m, N = 10) have no feasible layout at 5D and 6D (packing capacity 8 / 7); "
                             "'excl' = these two cases removed from the case set (10 cases), 'all' = all 12 split cases. "
                             "Contrasts never use them at 5D / 6D (neither method qualifies -> dropped); at 4D they are "
                             "feasible and enter the 'all' contrasts.",
               dropped_unpaired_runs=drop, spacing={})
    cases_all = set(SPLITCASES); cases_ex = cases_all - set(EXCL)
    for sp in SPACINGS:
        Gs = D[D.Spacing == sp].reset_index(drop=True)
        if not len(Gs):
            out["spacing"][sp] = None; continue
        S = MR.case_stats(Gs, SP9)
        o = dict(completeness=completeness(Gs, SP9, CASE), source="stored main runs (4D)" if sp == "4D" else "rev2_spacing")
        for tag, cs in (("excl", cases_ex), ("all", cases_all)):
            Gc = Gs[[c in cs for c in zip(Gs.Dataset, Gs.Radius, Gs.Turbines)]]
            Sc = S[[c in cs for c in zip(S.Dataset, S.Radius, S.Turbines)]]
            fr = ranking(Sc, SP9, "PSOBV")
            o[tag] = dict(n_cases=int(Sc[CASE].drop_duplicates().shape[0]), friedman=fr,
                          feasible_pct={a: float(100 * Gc[Gc.Algorithm == a].Feasible.mean()) for a in SP9 if (Gc.Algorithm == a).any()},
                          mean_loss_pct_qualified={a: float(Sc[(Sc.Algorithm == a) & Sc.Qualified].Loss.mean()) for a in SP9},
                          n_qualified={a: int(Sc[(Sc.Algorithm == a)].Qualified.sum()) for a in SP9},
                          contrasts={f"{a}-{b}": contrast(Gc, Sc, a, b) for a, b in CONTRASTS_SP})
        out["spacing"][sp] = o
        ex = o["excl"]
        if ex["friedman"]:
            log(f"[spacing] {sp}: avg ranks (10 cases) " + ", ".join(f"{a} {v:.2f}" for a, v in sorted(ex['friedman']['avg_rank'].items(), key=lambda t: t[1])))
        for k, r in ex["contrasts"].items():
            if r and r["case_level"]:
                log(f"[spacing] {sp} {k}: n={r['n_cases']} dL={r['mean_dloss_pp']:+.4f} case90={np.round(r['case_level']['ci90_mean_dloss_pp'], 4)} "
                    f"seed90={np.round(r['seed_level']['ci90'], 4) if r['seed_level'] else None} pW={r['wilcoxon']['p']:.3g}")
    # Holm over the 6 contrast tests (2 contrasts x 3 spacings, excl case set), for reference
    keys = [(sp, k) for sp in SPACINGS if out["spacing"].get(sp) for k, r in out["spacing"][sp]["excl"]["contrasts"].items() if r]
    if keys:
        ph = MR.holm([out["spacing"][sp]["excl"]["contrasts"][k]["wilcoxon"]["p"] for sp, k in keys])
        out["holm_over_contrasts_excl"] = {f"{sp}:{k}": float(h) for (sp, k), h in zip(keys, ph)}
    summ["spacing"] = out

    # --- table A: ranks and feasibility
    SPT = ptag(*[out["spacing"][sp]["completeness"] for sp in SPACINGS if out["spacing"].get(sp)])
    sps = [sp for sp in SPACINGS if out["spacing"].get(sp)]
    def rk(sp, tag, a):
        f = out["spacing"][sp][tag]["friedman"]
        return f["avg_rank"].get(a) if f else None
    lines = []
    for a in SP9:
        c = []
        for sp in sps:
            c += [fnum(rk(sp, "excl", a), 2, False) + " (" + fnum(rk(sp, "all", a), 2, False) + ")",
                  fnum(out["spacing"][sp]["excl"]["feasible_pct"].get(a), 1, False)]
        lines.append(f"{LAB[a]} & " + " & ".join(c) + " \\\\")
    fr_txt = []
    for sp in sps:
        f = out["spacing"][sp]["excl"]["friedman"]
        if f and "chi2" in f:
            fr_txt.append(f"{sp}: $\\chi^2_F={f['chi2']:.1f}$, $p={fp(f['p']).strip('$')}$, $F_{{ID}}={f['iman_davenport']:.1f}$")
    head = (" & " + " & ".join(f"\\multicolumn{{2}}{{c}}{{{sp} ({'stored runs' if sp == '4D' else 'rev.~2'})}}" for sp in sps) + " \\\\\n"
            + "".join(f"\\cmidrule(lr){{{2 + 2 * i}-{3 + 2 * i}}}" for i in range(len(sps))) + "\n"
            + "Method & " + " & ".join("Rank & Feas." for _ in sps))
    tex.append(MR.table(
        "table*", "Minimum spacing $4D$, $5D$, $6D$: average rank and feasibility of the nine methods on the split cases (6{,}030 evaluations, 30 seed-paired runs)." + SPT,
        "tab:R2-spacing-rank", "l" + "cc" * len(sps), head, lines, sep="3pt", pos="!htb",
        note="Rank: average feasibility-aware rank (Table~\\ref{tab:friedman68} rule: a method with fewer than 15 feasible runs "
             "of 30 in a case ranks below the qualified methods, by number of feasible runs) over the 10 split cases "
             "EXCLUDING the two cases $r=500$~m, $N=10$ (data sets I and II), which admit no feasible layout at $5D$ and $6D$ "
             "(packing capacity 8 and 7 turbines); in parentheses: over all 12 split cases (the two excluded cases then enter "
             "with all methods tied at $5D$/$6D$). Feas.: feasible runs (\\%) over the 10 cases. $4D$: stored runs of the main "
             "comparison restricted to the split cases. Friedman over the 10 cases: " + "; ".join(fr_txt) + "."))
    # --- table B: contrasts
    lines = []
    for a, b in CONTRASTS_SP:
        lines.append(f"\\multicolumn{{11}}{{l}}{{\\emph{{{LAB[a]} $-$ {LAB[b]}}}}} \\\\")
        for sp in sps:
            for tag in (("excl", "all") if sp == "4D" else ("excl",)):
                r = out["spacing"][sp][tag]["contrasts"].get(f"{a}-{b}")
                if not r or r["n_cases"] == 0:
                    lines.append(f"{sp} & -- & \\multicolumn{{9}}{{c}}{{no case in which both qualify}} \\\\"); continue
                cl, sl, hl = r["case_level"], r["seed_level"], r.get("hodges_lehmann")
                lab = sp + (" (12 cases)" if tag == "all" else "")
                lines.append(f"{lab} & {r['n_cases']} & {fnum(r['mean_dloss_pp'])} & {fci(sl['ci90'] if sl else None)} & {eqv(sl)} & "
                             f"{fci(cl['ci90_mean_dloss_pp'] if cl else None)} & {eqv(cl)} & "
                             f"{fnum(cl['min_margin_pp'], 3, False) if cl else '--'} & {fp(r['wilcoxon']['p'])} & "
                             f"{fnum(hl['hl']) if hl else '--'} & {fci(hl['ci95'] if hl else None)} \\\\")
    head = ("Spacing & $n$ & $\\overline{\\Delta L}$ (pp) & 90\\% CI seed & Eq. & 90\\% CI case & Eq. & $m_{\\min}$ & $p_W$ & HL & HL 95\\% CI")
    tex.append(MR.table(
        "table*", "Minimum spacing: component contrasts PSO-VNS $-$ PSO and SSA-VNS $-$ RSD-VNS at $4D$, $5D$ and $6D$ (split cases, 6{,}030 evaluations)." + SPT,
        "tab:R2-spacing-contrast", "lcccccccccc", head, lines, sep="2.5pt", pos="!htb",
        note="$\\overline{\\Delta L}$: mean difference of the per-case mean wake losses (pp of the wake-free benchmark objective; "
             "first minus second, negative = first better) over the $n$ cases in which both methods have at least 15 feasible "
             "runs of 30; the two cases $r=500$~m, $N=10$ are excluded (at $5D$/$6D$ no method has a feasible run; at $4D$ the "
             "row ``12 cases'' includes them). 90\\%% CI seed: percentile bootstrap of the seed-paired runs within the cases "
             "(fixed benchmark); 90\\%% CI case: percentile bootstrap over the cases (%s resamples, fixed seeds, as "
             "Table~\\ref{tab:equivalence}); Eq.: TOST equivalence at $m=%.2f$~pp (90\\%% CI inside $(-m,m)$); $m_{\\min}$: "
             "smallest equivalence margin (case level); $p_W$: two-sided Wilcoxon signed-rank test on the case means "
             "(unadjusted; cases in which only one method qualifies imputed as in Table~\\ref{tab:friedman68}); HL: "
             "Hodges--Lehmann estimate with exact 95\\%% CI." % (f"{MR.BOOT_N:,}".replace(",", "{,}"), EQ)))


# ------------------------------------------------------------------ T-laplace
def block_laplace(data_dir, summ, tex):
    N, st = read_rev2("laplace", data_dir)
    summ["data"]["rev2_laplace"] = st
    if N is None:
        log("[laplace] rev2_laplace missing -> skipped"); return
    A = main_runs()
    ref = A[A.Algorithm.isin(["SSA", "LXSSA"]) & pd.Series([c in SPLITCASES for c in zip(A.Dataset, A.Radius, A.Turbines)], index=A.index)]
    N = N[N.Algorithm.isin(LAPV)]
    N = N[[c in SPLITCASES for c in zip(N.Dataset, N.Radius, N.Turbines)]]
    pool = ["SSA", "LXSSA"] + [a for a in LAPV if a in set(N.Algorithm)]
    G, drop = pair_seeds(pd.concat([ref, N], ignore_index=True), CASE, pool)
    S = MR.case_stats(G, pool)
    fr = ranking(S, pool, "LXSSA")
    Rm = MR.rank_matrix(S, pool)
    allq = S.pivot_table(index=CASE, columns="Algorithm", values="Qualified").reindex(columns=pool).fillna(0).astype(bool).all(axis=1)
    out = dict(pool=pool, completeness=completeness(G, LAPPOOL, CASE), dropped_unpaired_runs=drop, friedman=fr,
               n_cases=int(S[CASE].drop_duplicates().shape[0]),
               calls_per_iteration={a: CALLS_PER_ITER[a][0] for a in pool}, iterations={a: CALLS_PER_ITER[a][1] for a in pool},
               feasible_pct={a: float(100 * G[G.Algorithm == a].Feasible.mean()) for a in pool},
               mean_loss_pct_qualified={a: float(S[(S.Algorithm == a) & S.Qualified].Loss.mean()) for a in pool},
               n_qualified={a: int(S[S.Algorithm == a].Qualified.sum()) for a in pool},
               n_cases_all_qualified=int(allq.sum()),
               mean_loss_pct_all_qualified={a: float(S.pivot_table(index=CASE, columns="Algorithm", values="Loss")[a][allq].mean()) for a in pool},
               avg_rank_among_pool={a: float(Rm[a].mean()) for a in pool}, vs={})
    for refm in ("LXSSA", "SSA"):
        oth = [a for a in pool if a != refm]
        cw = MR.case_mean_wilcoxon(S, refm, oth)                 # Holm over the 7 other methods of the pool
        wt = wtl_block(G, refm, oth)
        o = {}
        for a in oth:
            r = contrast(G, S, a, refm, seed_level=True)
            r["p_holm"] = cw[a]["p_holm"]
            r["holm_family"] = f"case-mean Wilcoxon, {refm} vs each of the other {len(oth)} methods of the pool"
            r["run_level_wtl_from_variant_side"] = dict(W=wt[a]["L"], T=wt[a]["T"], L=wt[a]["W"])
            o[a] = r
        out["vs"][refm] = o
    summ["laplace"] = out
    log("[laplace] avg ranks: " + ", ".join(f"{a} {out['avg_rank_among_pool'][a]:.2f}" for a in pool))
    for a in pool:
        for refm in ("LXSSA", "SSA"):
            r = out["vs"][refm].get(a)
            if r and r["case_level"]:
                log(f"[laplace] {a} - {refm}: dL={r['mean_dloss_pp']:+.4f} 90%CI={np.round(r['case_level']['ci90_mean_dloss_pp'], 4)} pHolm={r['p_holm']:.3g}")
    lines = []
    for a in pool:
        cpi, T = CALLS_PER_ITER[a]
        c = [f"{cpi} / {T}", fnum(out["mean_loss_pct_qualified"][a], 3, False) + (f"$^{{{out['n_qualified'][a]}}}$" if out["n_qualified"][a] < out["n_cases"] else ""),
             fnum(out["feasible_pct"][a], 1, False), fnum(out["avg_rank_among_pool"][a], 2, False)]
        for refm in ("LXSSA", "SSA"):
            r = out["vs"][refm].get(a)
            if r is None:
                c += ["--", "--", "--"]; continue
            cl = r["case_level"]
            c += [fnum(r["mean_dloss_pp"]), fci(cl["ci90_mean_dloss_pp"] if cl else None), fp(r["p_holm"])]
        lines.append(f"{LAB[a]} & " + " & ".join(c) + " \\\\")
    head = ("& & & & & \\multicolumn{3}{c}{vs.\\ LX-SSA} & \\multicolumn{3}{c}{vs.\\ SSA} \\\\\n\\cmidrule(lr){6-8}\\cmidrule(lr){9-11}\n"
            "Method & Calls/it. / $T$ & Loss (\\%) & Feas. & Rank & $\\overline{\\Delta L}$ & 90\\% CI & $p_{\\rm Holm}$ & $\\overline{\\Delta L}$ & 90\\% CI & $p_{\\rm Holm}$")
    ftxt = (f"Friedman $\\chi^2_F={fr['chi2']:.1f}$ ({len(pool) - 1} d.f.), $p={fp(fr['p']).strip('$')}$; Iman--Davenport $F_F={fr['iman_davenport']:.1f}$, "
            f"$p={fp(fr['iman_davenport_p']).strip('$')}$." if fr and "chi2" in fr else "")
    tex.append(MR.table(
        "table*", "Role and parameterization of the Laplace step: LX-SSA variants on the %d split cases (6{,}030 evaluations, 30 seed-paired runs)." % out["n_cases"] + ptag(out["completeness"]),
        "tab:R2-laplace", "lcccccccccc", head, lines, sep="2.5pt", pos="!htb",
        note="Calls/it.\\ / $T$: objective evaluations per iteration and number of iterations at 6{,}030 evaluations. Loss: mean "
             "wake loss (\\%%) of the feasible runs, averaged over the cases in which the method has at least 15 feasible runs of "
             "30 (superscript: number of such cases when fewer than %d). Feas.: feasible runs (\\%%). Rank: average "
             "feasibility-aware rank among the %d methods of this table. $\\overline{\\Delta L}$: mean difference of the "
             "per-case mean wake losses, method minus reference (pp; negative = method better) over the cases in which both "
             "qualify, with 90\\%% percentile bootstrap CI over the cases; $p_{\\rm Holm}$: case-mean Wilcoxon signed-rank test, "
             "Holm-adjusted over the %d comparisons with the same reference. LX-SSA and SSA: stored runs of the main "
             "comparison. %s" % (out["n_cases"], len(pool), len(pool) - 1, ftxt)))


# ------------------------------------------------------------------ T-ga
def hr_table(H, meth, inst, label, caption, summ_key, summ):
    """Horns Rev 16 table in the format of tab:hr-site (main_text_tables of mpce_results.py)."""
    ideal = float(H.Ideal.iloc[0])
    settings = [("6030R", 6030, "random"), ("6030F", 6030, "feasible"), ("30030R", 30030, "random"), ("120030R", 120030, "random")]
    lines = [f"Installed layout & {inst:.2f} & -- & {100 * (1 - inst / ideal):.2f} & -- & -- & -- \\\\", "\\midrule"]
    res = {}
    for a in meth:
        x6 = H[(H.Algorithm == a) & (H.Budget == 6030) & (H.Init == "random")]
        nf6 = int(x6.Feasible.sum())
        c = [f"{x6[x6.Feasible].Objective.mean():.2f}" if nf6 else "--", f"{nf6}/{len(x6)}" if len(x6) else "n/r"]
        res[a] = {}
        for tag, b, init in settings:
            y = H[(H.Algorithm == a) & (H.Budget == b) & (H.Init == init)]
            if not len(y):
                c.append("n/r"); res[a][tag] = None; continue
            nf = int(y.Feasible.sum())
            v = float(y[y.Feasible].LossPct.mean()) if nf else np.nan
            res[a][tag] = dict(loss_pct=v, feasible=nf, runs=int(len(y)), mean_aep=float(y[y.Feasible].Objective.mean()) if nf else None,
                               best_aep=float(y[y.Feasible].Objective.max()) if nf else None,
                               runs_above_installed=int((y.Feasible & (y.Objective > inst)).sum()))
            cell = "--" if not nf else f"{v:.2f}"
            if 0 < nf < len(y):
                cell += f"$^{{{nf}}}$"
            c.append(cell)
        lines.append(f"{LAB[a]} & " + " & ".join(c) + " \\\\")
    summ[summ_key] = dict(installed_aep=inst, ideal_aep=ideal, installed_loss_pct=100 * (1 - inst / ideal), methods=res)
    head = ("& \\multicolumn{2}{c}{6{,}030 evaluations, R} & \\multicolumn{4}{c}{Loss (\\%)} \\\\\n\\cmidrule(lr){2-3}\\cmidrule(lr){4-7}\n"
            "Layout / method & AEP & Feas. & 6k R & 6k F & 30k & 120k")
    return MR.table("table", caption, label, "lcccccc", head, lines, sep="3pt", pos="!htb",
                    note="Wake-free AEP %.2f GWh/yr. Format of Table~\\ref{tab:hr-site}. AEP: mean (GWh/yr) over the feasible runs "
                         "of 30 seeds at 6{,}030 evaluations, random initialization (R); Feas.: feasible runs. Loss: mean AEP loss "
                         "(\\%%, relative to the wake-free AEP) of the feasible runs at 6{,}030 evaluations with random (R) and "
                         "feasibility-preserving (F) initialization and at 30{,}030 and 120{,}030 evaluations (R; 10 seeds at "
                         "120{,}030); superscript: feasible runs when not all are feasible; n/r: not run." % ideal)


def iea_table(IE, meth, label, caption, summ, key, extra_note=""):
    """IEA37 CS1 in the format of tab:iea37 (best feasible AEP), extended by mean AEP and feasible runs."""
    pub = MR.load_published(HERE)
    pubv = {}
    for n in (16, 36):
        P = pub[pub.Turbines == n] if pub is not None else None
        base = float(P[P.Baseline].AEP.iloc[0]) if P is not None and P.Baseline.any() else MR.BASELINE_IEA[n]
        pf = P[(~P.Baseline) & P.Feasible] if P is not None else None
        pubv[n] = dict(base=base, best_feasible=float(pf.AEP.max()) if pf is not None and len(pf) else np.nan)
    try:
        J = json.load(open(os.path.join(HERE, "iea37_projected.json")))
        for n in (16, 36):
            pubv[n]["best_projected"] = float(J["scenarios"][str(n)]["projected"]["best_aep"])
    except (OSError, KeyError, ValueError):
        pass
    cols = [(n, b) for n in (16, 36) for b in (6030, 30030)]
    lines = []
    for lab_, k in (("Example layout", "base"), ("Best published, strict", "best_feasible"), ("Best published, projected", "best_projected")):
        if all(k in pubv[n] for n in (16, 36)):
            lines.append(f"{lab_} & " + " & ".join(f"\\multicolumn{{6}}{{c}}{{{fbig(pubv[n][k])}}}" for n in (16, 36)) + " \\\\")
    lines.append("\\midrule")
    res = {}
    for a in meth:
        c = []
        for n, b in cols:
            y = IE[(IE.Turbines == n) & (IE.Budget == b) & (IE.Algorithm == a)]
            if not len(y):
                c += ["n/r", "n/r", "n/r"]; continue
            yf = y[y.Feasible]
            best = float(yf.Objective.max()) if len(yf) else np.nan
            mean = float(yf.Objective.mean()) if len(yf) else np.nan
            res.setdefault(a, {})[f"{n}T_{b}"] = dict(best=best, mean=mean, sd=float(yf.Objective.std()) if len(yf) > 1 else None,
                                                     loss_pct_mean=float(yf.LossPct.mean()) if len(yf) else None,
                                                     feasible=int(len(yf)), runs=int(len(y)),
                                                     best_vs_best_feasible_published_pct=float(100 * (best / pubv[n]["best_feasible"] - 1)) if np.isfinite(best) and np.isfinite(pubv[n]["best_feasible"]) else None)
            c += [fbig(best), fbig(mean), f"{len(yf)}/{len(y)}"]
        lines.append(f"{LAB[a]} & " + " & ".join(c) + " \\\\")
    summ[key] = dict(published={str(n): v for n, v in pubv.items()}, methods=res)
    head = ("& \\multicolumn{6}{c}{16 turbines ($r=1300$~m)} & \\multicolumn{6}{c}{36 turbines ($r=2000$~m)} \\\\\n"
            "\\cmidrule(lr){2-7}\\cmidrule(lr){8-13}\n"
            "& \\multicolumn{3}{c}{6{,}030} & \\multicolumn{3}{c}{30{,}030} & \\multicolumn{3}{c}{6{,}030} & \\multicolumn{3}{c}{30{,}030} \\\\\n"
            "\\cmidrule(lr){2-4}\\cmidrule(lr){5-7}\\cmidrule(lr){8-10}\\cmidrule(lr){11-13}\n"
            "Method / source & Best & Mean & Feas. & Best & Mean & Feas. & Best & Mean & Feas. & Best & Mean & Feas.")
    return MR.table("table*", caption, label, "l" + "ccc" * 4, head, lines, sep="2pt", resize=True, pos="!htb",
                    note="AEP (MWh, official IEA37 model) of the best feasible run and mean over the feasible runs (30 seeds, random "
                         "initialization); Feas.: feasible runs; --: no feasible run; n/r: not run. Published rows as in "
                         "Table~\\ref{tab:iea37} (strict: every turbine within 1~mm of the boundary; projected: turbines outside "
                         "projected radially onto it). " + extra_note)


def block_ga(data_dir, summ, tex):
    N, st = read_rev2("ga", data_dir)
    summ["data"]["rev2_ga"] = st
    if N is not None:
        A = main_runs()
        pool = MAIN8 + ["GA"]
        G, drop = pair_seeds(pd.concat([A[A.Algorithm.isin(MAIN8)], N[N.Algorithm == "GA"]], ignore_index=True), CASE, pool)
        S = MR.case_stats(G, pool)
        oth = [a for a in pool if a != "GA"]
        fr_ga = ranking(S, pool, "GA")
        fr_pv = ranking(S, pool, "PSOBV")
        fr8 = ranking(MR.case_stats(G, MAIN8), MAIN8, "PSOBV")
        cw = MR.case_mean_wilcoxon(S, "GA", oth)
        wt = wtl_block(G, "GA", oth)
        out = dict(pool=pool, completeness=completeness(G, pool, CASE), dropped_unpaired_runs=drop,
                   n_cases=int(S[CASE].drop_duplicates().shape[0]), friedman=fr_ga, friedman_vs_psovns=fr_pv,
                   avg_rank_without_ga=fr8["avg_rank"] if fr8 else None,
                   feasible_pct={a: float(100 * G[G.Algorithm == a].Feasible.mean()) for a in pool},
                   mean_loss_pct_qualified={a: float(S[(S.Algorithm == a) & S.Qualified].Loss.mean()) for a in pool},
                   n_qualified={a: int(S[S.Algorithm == a].Qualified.sum()) for a in pool},
                   case_mean_wilcoxon_ga_vs=cw, wtl_ga_vs=wt, contrasts={})
        for b in oth:
            out["contrasts"][f"GA-{b}"] = contrast(G, S, "GA", b, seed_level=False)
        summ["ga68"] = out
        log("[ga] avg ranks (9 methods): " + ", ".join(f"{a} {v:.2f}" for a, v in sorted(fr_ga["avg_rank"].items(), key=lambda t: t[1])))
        order = sorted(pool, key=lambda a: fr_ga["avg_rank"][a])
        lines = []
        for a in order:
            if a == "GA":
                c = ["--"] * 5
            else:
                r = out["contrasts"][f"GA-{a}"]
                cl = r["case_level"] if r else None
                c = [fp(fr_ga["p_holm_vs_focus"][a]), fp(cw[a]["p_holm"]), fnum(r["mean_dloss_pp"] if r else None),
                     fci(cl["ci90_mean_dloss_pp"] if cl else None), wtl_s(wt[a])]
            lines.append(f"{LAB[a]} & {fr_ga['avg_rank'][a]:.2f} & {fr_ga['sole_best_count'][a]} & {out['feasible_pct'][a]:.1f} & "
                         + " & ".join(c) + " \\\\")
        head = ("Method & Avg.\\ rank & Best & Feas. & $p_z$ & $p_W$ & $\\overline{\\Delta L}$ (GA $-$ method) & 90\\% CI & W/T/L (GA)")
        tex.append(MR.table(
            "table*", "Genetic algorithm (GA) added to the method pool: case-level analysis over the %d benchmark cases (6{,}030 evaluations, 30 seed-paired runs)." % out["n_cases"] + ptag(out["completeness"]),
            "tab:R2-ga", "lcccccccc", head, lines, sep="2.5pt", pos="!htb",
            note="Avg.\\ rank: average feasibility-aware rank among the nine methods (rule of Table~\\ref{tab:friedman68}); Best: "
                 "cases in which the method alone ranks first; Feas.: feasible runs (\\%%); $p_z$: Holm-adjusted $p$ of the "
                 "average-rank test against GA; $p_W$: Wilcoxon signed-rank test on the per-case mean wake losses, GA vs.\\ the "
                 "method (Holm-adjusted over the 8 methods; cases in which only one qualifies imputed); $\\overline{\\Delta L}$: "
                 "mean difference of the case-mean wake losses, GA minus method (pp; positive = GA worse), over the cases in "
                 "which both have at least 15 feasible runs, 90\\%% bootstrap CI over the cases; W/T/L: cases in which GA is "
                 "significantly better / not different / worse (seed-paired Wilcoxon, Holm over the 8 comparisons of GA in "
                 "each case). Friedman $\\chi^2_F=%.1f$ (%d d.f.), $p=%s$; Iman--Davenport $F_F=%.1f$, $p=%s$. Without GA the "
                 "ranks of the eight methods are those of Table~\\ref{tab:friedman68}." %
                 (fr_ga["chi2"], len(pool) - 1, fp(fr_ga["p"]).strip("$"), fr_ga["iman_davenport"], fp(fr_ga["iman_davenport_p"]).strip("$"))))
    else:
        log("[ga] rev2_ga missing -> 68-case GA analysis skipped")
    # --- Horns Rev 16
    Nh, st = read_rev2("gahr", data_dir)
    summ["data"]["rev2_gahr"] = st
    if Nh is not None:
        H = read_mpce("hrfix")
        H = H[(H.Dataset == "HR") & (H.Turbines == 16)]
        Nh = Nh[(Nh.Algorithm == "GA") & (Nh.Dataset == "HR")]
        if abs(float(Nh.Ideal.iloc[0]) - float(H.Ideal.iloc[0])) > 1e-6:
            NOTES.append("GA Horns Rev wake-free AEP differs from mpce_hrfix: model mismatch")
        HH = pd.concat([H, Nh], ignore_index=True)
        import hornsrev_model as hr
        inst = float(hr.aep_gwh(hr.site(16)[0]))
        meth = HR10 + ["GA"]
        tex.append(hr_table(HH, meth, inst, "tab:R2-ga-hr",
                            "Horns Rev~1 16-turbine block with GA added (5$^\\circ$ direction bins): AEP, feasibility and AEP loss.",
                            "ga_hr16", summ))
        # run-level tests of GA vs each method (seed-paired, Holm over the comparisons) per budget
        tests = {}
        for b in (6030, 30030):
            y = HH[(HH.Budget == b) & (HH.Init == "random")]
            y, _ = pair_seeds(y, CASE, [a for a in meth if a in set(y.Algorithm)])
            tests[b] = {t["Baseline"]: dict(p_holm=t["PHolm"], rb=t["RB"], outcome=t["Outcome"], n_pairs=t["NPairs"])
                        for t in MR.paired_vs(y, "GA", [a for a in meth if a != "GA" and a in set(y.Algorithm)])}
        summ["ga_hr16"]["ga_run_level_tests"] = tests
        g = summ["ga_hr16"]["methods"]["GA"]
        log(f"[ga] HR16 GA: " + ", ".join(f"{k}: {v['feasible']}/{v['runs']} loss {v['loss_pct']:.2f}" for k, v in g.items() if v))
    else:
        log("[ga] rev2_gahr missing -> Horns Rev GA rows skipped")
    # --- IEA37
    Ni, st = read_rev2("gaiea", data_dir)
    summ["data"]["rev2_gaiea"] = st
    if Ni is not None:
        IE = pd.concat([d for d in (read_mpce(e) for e in ("iea16", "iea36", "iea16p", "iea36p")) if d is not None]
                       + [Ni[Ni.Algorithm == "GA"]], ignore_index=True)
        IE = IE[IE.Dataset.str.startswith("IEA37")]
        tex.append(iea_table(IE, MAIN8 + ["GA"], "tab:R2-ga-iea",
                             "IEA37 Case Study~1 with GA added: AEP (MWh) of the best run, mean AEP and feasible runs.",
                             summ, "ga_iea37"))
        log("[ga] IEA37 GA: " + json.dumps({k: (round(v["best"], 1), v["feasible"]) for k, v in summ["ga_iea37"]["methods"].get("GA", {}).items()}))
    else:
        log("[ga] rev2_gaiea missing -> IEA37 GA rows skipped")


# ------------------------------------------------------------------ T-site (Lillgrund)
def block_site(data_dir, summ, tex):
    parts = []
    for exp in ("lg16", "lg16b", "lg16all"):
        d, st = read_rev2(exp, data_dir)
        summ["data"][f"rev2_{exp}"] = st
        if d is not None:
            parts.append(d)
    if not parts:
        log("[site] rev2_lg16 / rev2_lg16b missing -> skipped"); return
    L = pd.concat(parts, ignore_index=True).drop_duplicates(MR.KEY, keep="last")
    L = L[(L.Dataset == "LG") & L.Algorithm.isin(SP9)]
    import rev2_site_model as lg
    xy0, _ = lg.site(16)
    inst = float(lg.aep_gwh(xy0)); ideal_m = float(lg.aep_gwh(xy0, with_wake=False))
    ideal = float(L.Ideal.iloc[0])
    if abs(ideal - ideal_m) > 1e-6:
        NOTES.append(f"Lillgrund wake-free AEP of the runs ({ideal:.4f}) differs from rev2_site_model ({ideal_m:.4f})")
    out = dict(installed_aep=inst, ideal_aep=ideal, installed_loss_pct=100 * (1 - inst / ideal), smin_m=float(lg.SMIN), budgets={})
    buds = [b for b in (6030, 30030) if (L.Budget == b).any()]
    for b in buds:
        y = L[(L.Budget == b) & (L.Init == "random")]
        mets = [a for a in SP9 if a in set(y.Algorithm)]
        y, drop = pair_seeds(y, CASE, mets)
        Sy = MR.case_stats(y, mets).set_index("Algorithm")
        tt = {t["Baseline"]: t for t in MR.paired_vs(y, "PSOBV", [a for a in mets if a != "PSOBV"])} if "PSOBV" in mets else {}
        piv = y.assign(S=MR.goodness(y)).pivot_table(index="Seed", columns="Algorithm", values="S").dropna()
        try:
            pf = float(friedmanchisquare(*[piv[a] for a in mets])[1]) if len(piv) > 1 else None
        except ValueError:
            pf = None
        out["budgets"][b] = dict(dropped_unpaired_runs=drop, completeness=completeness(y, mets, CASE), run_level_friedman_p=pf, methods={
            a: dict(mean=float(Sy.Mean[a]) if Sy.NFeas[a] else None, sd=float(Sy.SD[a]) if Sy.NFeas[a] > 1 else None,
                    best=float(Sy.Best[a]) if Sy.NFeas[a] else None, loss_pct=float(Sy.Loss[a]) if Sy.NFeas[a] else None,
                    feasible=int(Sy.NFeas[a]), runs=int(Sy.N[a]), rank=float(Sy.Rank[a]),
                    mean_vs_installed_pct=float(100 * (Sy.Mean[a] / inst - 1)) if Sy.NFeas[a] else None,
                    runs_above_installed=int(((y.Algorithm == a) & y.Feasible & (y.Objective > inst)).sum()),
                    p_holm_psovns_vs=tt.get(a, {}).get("PHolm"), rb_psovns_vs=tt.get(a, {}).get("RB"))
            for a in mets})
        log(f"[site] LG16 {b}: " + ", ".join(f"{a} {v['mean'] if v['mean'] is None else round(v['mean'], 2)}({v['feasible']}) r{v['rank']:g}"
                                            for a, v in out["budgets"][b]["methods"].items()) + f"; installed {inst:.2f}")
    summ["lillgrund"] = out
    lines = ["Installed block & " + " & ".join([f"{inst:.2f}", "--", f"{100 * (1 - inst / ideal):.2f}", "0.00", "--"] * len(buds)) + " \\\\", "\\midrule"]
    order = sorted(SP9, key=lambda a: out["budgets"][buds[0]]["methods"].get(a, {}).get("rank", 99))
    for a in order:
        c = []
        for b in buds:
            v = out["budgets"][b]["methods"].get(a)
            if v is None:
                c += ["n/r"] * 5; continue
            c += [f"{v['mean']:.2f}" if v["mean"] is not None else "--", f"{v['feasible']}/{v['runs']}",
                  f"{v['loss_pct']:.2f}" if v["loss_pct"] is not None else "--",
                  fnum(v["mean_vs_installed_pct"], 2), f"{v['rank']:g}"]
        lines.append(f"{LAB[a]} & " + " & ".join(c) + " \\\\")
    head = (" & " + " & ".join(f"\\multicolumn{{5}}{{c}}{{{b:,} evaluations}}".replace(",", "{,}") for b in buds) + " \\\\\n"
            + "".join(f"\\cmidrule(lr){{{2 + 5 * i}-{6 + 5 * i}}}" for i in range(len(buds))) + "\n"
            + "Layout / method & " + " & ".join("AEP & Feas. & Loss & $\\Delta_{\\rm inst}$ & Rank" for _ in buds))
    fr_txt = "; ".join(f"{b:,}".replace(",", "{,}") + f": $p={fp(out['budgets'][b]['run_level_friedman_p']).strip('$')}$"
                       for b in buds if out["budgets"][b]["run_level_friedman_p"] is not None)
    tex.append(MR.table(
        "table*", "Lillgrund 16-turbine block (second measured-wind site; 5$^\\circ$ direction bins, minimum spacing $3D$): AEP, feasibility, AEP loss and rank." + ptag(*[out["budgets"][b]["completeness"] for b in buds]),
        "tab:R2-site", "l" + "ccccc" * len(buds), head, lines, sep="2.5pt", pos="!htb",
        note="Siemens SWT-2.3-93 ($D=93$~m), 12-sector Weibull climate of the Lillgrund met mast, Jensen wake $k=0.04$, same "
             "model and bins as Table~\\ref{tab:hr-site}; wake-free AEP %.2f GWh/yr; minimum spacing $3D=279$~m ($4D$ is "
             "geometrically infeasible for 16 turbines in this block). AEP: mean (GWh/yr) over the feasible runs of 30 "
             "seed-paired runs, random initialization; Feas.: feasible runs; Loss: mean AEP loss (\\%%) relative to the "
             "wake-free AEP; $\\Delta_{\\rm inst}$: mean AEP relative to the installed block (\\%%); Rank: feasibility-aware "
             "rank of the mean AEP (fewer than 15 feasible runs = ranked below, by feasible runs). Run-level Friedman test "
             "(infeasible runs below every feasible run): %s." % (ideal, fr_txt)))


# ------------------------------------------------------------------ T-grad
def block_grad(data_dir, summ, tex):
    N, st = read_rev2("grad", data_dir)
    summ["data"]["rev2_grad"] = st
    if N is None:
        log("[grad] rev2_grad missing -> skipped"); return
    extra = [c for c in ("FunCalls", "GradCalls", "CG", "Unused") if c in N]
    IE = pd.concat([d for d in (read_mpce(e) for e in ("iea16", "iea36", "iea16p", "iea36p")) if d is not None]
                   + [N[N.Algorithm.isin(GRADM)]], ignore_index=True)
    IE = IE[IE.Dataset.str.startswith("IEA37")]
    gm = [a for a in GRADM if a in set(N.Algorithm)]
    tex.append(iea_table(IE, MAIN8 + gm, "tab:R2-grad",
                         "IEA37 Case Study~1 with exact-gradient methods added: AEP (MWh) of the best run, mean AEP and feasible runs.",
                         summ, "grad_iea37",
                         extra_note="MS-SLSQP: forward-difference gradients (paper); exact gradient: analytic AEP gradient, each "
                                    "gradient charged $c_g$ objective evaluations (Table~\\ref{tab:R2-grad-calls})."))
    calls, lines = {}, []
    for n in (16, 36):
        for b in (6030, 30030):
            y = IE[(IE.Turbines == n) & (IE.Budget == b) & (IE.Init == "random")]
            if not len(y):
                continue
            mets = [a for a in MAIN8 + gm if a in set(y.Algorithm)]
            y, _ = pair_seeds(y, CASE, mets)
            tests = {}
            for a, refs in (("SLSQPX", ["SLSQP"]), ("PSOSLSQPX", ["PSOBV", "PSOC", "SLSQP"])):
                if a in mets:
                    tests[a] = {t["Baseline"]: dict(p_holm=t["PHolm"], rb=t["RB"], outcome=t["Outcome"], n_pairs=t["NPairs"])
                                for t in MR.paired_vs(y, a, [r for r in refs if r in mets])}
            for a in gm:
                x = N[(N.Algorithm == a) & (N.Turbines == n) & (N.Budget == b)]
                if not len(x):
                    continue
                v = dict(runs=int(len(x)), calls_mean=float(x.Calls.mean()))
                for c in extra:
                    v[c.lower() + "_mean"] = float(x[c].mean())
                if "CG" in x:
                    v["cg"] = sorted(set(int(u) for u in x.CG.dropna()))
                v["tests_vs"] = tests.get(a, {})
                calls[f"{a}:{n}T_{b}"] = v
                t = tests.get(a, {})
                lines.append(f"{LAB[a]} & {n} & {fbig(b, 0)} & {fbig(v.get('funcalls_mean', np.nan), 0)} & "
                             f"{fbig(v.get('gradcalls_mean', np.nan), 1)} & " +
                             f"{', '.join(map(str, v.get('cg', []))) or '--'} & {v.get('unused_mean', np.nan):.2f} & " +
                             " & ".join(f"{fp(t[r]['p_holm'])} ({t[r]['outcome']})" if r in t else "--" for r in ("SLSQP", "PSOBV", "PSOC")) + " \\\\")
    summ["grad_iea37"]["calls"] = calls
    summ["grad_iea37"]["tests_note"] = ("seed-paired Wilcoxon of the method vs each reference (infeasible runs below every feasible "
                                        "run), Holm over the references of the method; outcome W/T/L from the method's side")
    tex.append(MR.table(
        "table*", "Exact-gradient methods on IEA37 Case Study~1: evaluation accounting and run-level tests.",
        "tab:R2-grad-calls", "lccccccccc",
        "Method & $N$ & Budget & FunCalls & GradCalls & $c_g$ & Unused & vs.\\ MS-SLSQP & vs.\\ PSO-VNS & vs.\\ PSO",
        lines, sep="2.5pt", pos="!htb",
        note="Means over the 30 runs. FunCalls: objective evaluations; GradCalls: exact gradient evaluations, each charged "
             "$c_g$ evaluations (measured cost ratio, rounded up); Unused: evaluations left when a gradient no longer fit "
             "into the budget, so FunCalls + $c_g\\,$GradCalls + Unused = budget. Last columns: seed-paired Wilcoxon $p$ "
             "(Holm over the references of the method; W/T/L from the method's side) against MS-SLSQP with "
             "forward-difference gradients and, for PSO-SLSQP, against PSO-VNS and PSO (same PSO phase 1)."))
    for k, v in calls.items():
        log(f"[grad] {k}: " + ", ".join(f"{kk}={vv}" for kk, vv in v.items() if kk != "tests_vs"))


# ------------------------------------------------------------------ reproduction of the paper's numbers
def reproduction_check(summ):
    """The reused pipeline, run on the paper's 68-case data, must reproduce mpce_summary.json (8-method average ranks,
    PSO-VNS - PSO case-level TOST) and mpce_summary_extra.json (seed-level 90 % CI, Hodges-Lehmann, case-mean Wilcoxon)."""
    A = main_runs()
    G = A[A.Algorithm.isin(MAIN8)].reset_index(drop=True)
    S = MR.case_stats(G, MAIN8)
    fr = ranking(S, MAIN8, "PSOBV")
    r = contrast(G, S, "PSOBV", "PSOC")
    chk = {}
    try:
        ms = json.load(open(os.path.join(HERE, "mpce_summary.json")))
        chk["avg_rank_max_abs_diff"] = max(abs(fr["avg_rank"][a] - ms["main"]["friedman"]["avg_rank"][a]) for a in MAIN8)
        e = ms["equivalence"]["pairs"]["PSOBV-PSOC"]
        chk["case_ci90_max_abs_diff"] = max(abs(u - v) for u, v in zip(r["case_level"]["ci90_mean_dloss_pp"], e["ci90_mean_dloss_pp"]))
        chk["case_p_tost_diff"] = abs(r["case_level"]["p_tost"] - e["p_tost"])
    except (OSError, KeyError) as ex:
        chk["mpce_summary_error"] = str(ex)
    try:
        mx = json.load(open(os.path.join(HERE, "mpce_summary_extra.json")))
        sd = mx["equivalence_levels"]["pairs"]["PSOBV-PSOC"]["seed"]
        chk["seed_ci90_max_abs_diff"] = max(abs(u - v) for u, v in zip(r["seed_level"]["ci90"], sd["ci90"]))
        h = mx["hodges_lehmann"]["PSOBV-PSOC"]
        chk["hl_diff"] = abs(r["hodges_lehmann"]["hl"] - h["hl"]) + sum(abs(u - v) for u, v in zip(r["hodges_lehmann"]["ci95"], h["ci95"]))
        chk["wilcoxon_p_diff"] = abs(r["wilcoxon"]["p"] - mx["case_level"]["PSOBV-PSOC"]["p"])
    except (OSError, KeyError) as ex:
        chk["mpce_summary_extra_error"] = str(ex)
    chk["pass"] = all(v < 1e-9 for k, v in chk.items() if k.endswith("diff"))
    summ["reproduction_check"] = chk
    log("[check] reproduction of the paper (68 cases, PSO-VNS - PSO, 8-method ranks): "
        + ", ".join(f"{k}={v:.2e}" if isinstance(v, float) else f"{k}={v}" for k, v in chk.items()))


# ------------------------------------------------------------------ main
def clean(o):
    if isinstance(o, dict):
        return {str(k): clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [clean(v) for v in o]
    if isinstance(o, (np.bool_,)):
        return bool(o)
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating, float)):
        return None if not np.isfinite(o) else float(o)
    return o


BLOCKS = dict(spacing=block_spacing, laplace=block_laplace, ga=block_ga, site=block_site, grad=block_grad)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default=HERE, help="folder of the rev2_*.csv files")
    ap.add_argument("--out-dir", default=HERE)
    ap.add_argument("--only", default=",".join(BLOCKS))
    args = ap.parse_args(argv)
    t0 = time.time()
    os.makedirs(args.out_dir, exist_ok=True)
    summ = dict(generated=time.strftime("%Y-%m-%d %H:%M:%S"), data_dir=os.path.abspath(args.data_dir),
                protocol=dict(loss="wake loss (%) = 100 WakeLoss / Ideal (wake-free benchmark objective; sites: AEP loss "
                                    "relative to the wake-free AEP); Delta L in pp = difference of per-case means of the "
                                    "feasible runs, first minus second (negative = first better)",
                              ranking=MR.RANK_RULE, equivalence_margin_pp=EQ,
                              case_bootstrap=dict(resamples=MR.BOOT_N, seed=MR.BOOT_SEED),
                              seed_bootstrap=dict(resamples=MX.BOOT_N, seed=MX.SEED_BOOT_SEED, version="stratified (within-case)"),
                              functions="mpce_results: std_cols, case_stats, rank_rule, rank_matrix, friedman_block, "
                                        "case_mean_wilcoxon, case_mean_diffs, tost, paired_vs, goodness, holm, fmt_p, table, "
                                        "load_published; mpce_inference_extra: load, cm_test, seed_level, hodges_lehmann"),
                data={}, notes=NOTES)
    tex = ["%% generated by rev2_analysis.py -- do not edit by hand\n"]
    try:
        reproduction_check(summ)
    except Exception as e:                                 # pragma: no cover
        summ["reproduction_check"] = dict(error=f"{type(e).__name__}: {e}", **{"pass": False})
    for name in args.only.split(","):
        name = name.strip()
        if not name:
            continue
        log(f"\n==== {name}")
        try:
            BLOCKS[name](args.data_dir, summ, tex)
        except Exception as e:                             # one failing block must not stop the others
            import traceback
            traceback.print_exc()
            summ.setdefault("errors", {})[name] = f"{type(e).__name__}: {e}"
    with open(os.path.join(args.out_dir, "rev2_summary.json"), "w") as fh:
        json.dump(clean(summ), fh, indent=1)
    with open(os.path.join(args.out_dir, "rev2_tables.tex"), "w") as fh:
        fh.write("\n".join(tex))
    log(f"\nwrote {os.path.join(args.out_dir, 'rev2_summary.json')} and rev2_tables.tex ({len(tex) - 1} tables) in {time.time() - t0:.0f} s")
    if summ.get("errors"):
        log(f"ERRORS: {summ['errors']}")


if __name__ == "__main__":
    main()
