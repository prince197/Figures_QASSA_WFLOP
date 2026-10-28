"""GA in the remaining supplement tables whose main-comparison rows still use the eight-method pool:
S-tab:X-cluster, S-tab:X-loco, S-tab:X-hl (rows of ALL_PAIRS), S-tab:X-equiv-levels, S-tab:X-bayes (EQ_PAIRS; main-
comparison rows PSO-VNS vs each method) and S-tab:cost (cost per evaluation by method).

Usage (from any directory):  python3 analysis/rev5_ga_infer2.py
Outputs (analysis/; no existing file is modified, nothing under SWEVO_rev2/ is written, no optimizer is run):
  rev5_ga_inference_extra2_tables.tex  S-tab:X-cluster, S-tab:X-loco, S-tab:X-equiv-levels, S-tab:X-bayes, S-tab:X-hl
                                       with the pair PSO-VNS vs GA (complete replacements, same labels, wording of the
                                       latex_source copies of mpce_supp_inference.tex)
  rev5_ga_inference_extra2_macros.tex  \\renewcommand of every \\NX... macro of mpce_numbers_extra.tex whose value changes
                                       with the pair PSO-VNS vs GA in these analyses, and \\newcommand of the new
                                       PSO-VNS-vs-GA macros (input after mpce_numbers_extra.tex and
                                       rev5_ga_inference_extra_macros.tex)
  rev5_ga_cost_tables.tex              S-tab:cost with the GA row (stored rev2_ga runs; separate batch)
  rev5_ga_infer2.json, rev5_ga_infer2.log

Procedure (the generating functions are imported and called):
  X-cluster / X-loco / X-hl / X-equiv-levels / X-bayes:
      mpce_inference_extra.cluster_block / loco_block / hl_block / equivalence_levels_block on the per-run data of
      mpce_inference_extra.load + the GA records rev2_ga_s*of10.csv (68 cases, seeds 1-30, 6,030 evaluations, random;
      std_cols as all other records). The pair lists are extended by PSO-VNS vs GA (ALL_PAIRS: after PSO-VNS vs DE, as
      in tab:friedman68-full; EQ_PAIRS: main-comparison rows, after PSO-VNS vs DE); the Friedman ranking of the LOCO
      runs uses the nine-method pool (module constant MAIN8 set to the nine methods for that call). Table text by
      mpce_inference_extra.write_tables, macros by write_macros. None of these tables has a multiplicity adjustment
      (every p is unadjusted; the TOST / Bayesian / bootstrap results of each pair use their own fixed seeds), so the
      rows of the existing pairs cannot change; this is asserted.
  tab:cost: the statements of mpce_results.main (section "computational cost") copied into cost_rows, first checked
      to reproduce the stored table of mpce_supplementary.tex exactly, then applied to the GA runs
      (mpce_results.load_ga_secondary, source rev2_ga, 2,040 runs).
Step 0 (reproduction): every block is recomputed from the run records and asserted to equal mpce_summary_extra.json;
write_tables / write_macros reproduce mpce_supp_inference.tex and mpce_numbers_extra.tex; after the wording edits of
the latex_source copies, the five tables equal those of latex_source/analysis/mpce_supp_inference.tex; tab:cost equals
analysis/ and latex_source/analysis/mpce_supplementary.tex.
One process; OMP / BLAS threads = 1.
"""
import os, sys, re, glob, json, time, copy, difflib, tempfile
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import mpce_inference_extra as X            # noqa: E402
import mpce_results as MR                   # noqa: E402

LS = os.path.join(ROOT, "SWEVO_rev2", "latex_source")
SCR = tempfile.mkdtemp(prefix="rev5_ga_infer2_")
LOG = []
MAIN8 = list(X.MAIN8)
MAIN9 = MAIN8[:6] + ["GA"] + MAIN8[6:]          # PSOBV PSOC SSABV SSA LXSSA DE GA BVNS SLSQP (tab:friedman68-full)
assert MAIN8[5] == "DE"
ALL8 = list(X.ALL_PAIRS)
MAIN_PAIRS9 = [(X.FOCUS, b) for b in MAIN9 if b != X.FOCUS]
ALL9 = list(dict.fromkeys(MAIN_PAIRS9 + X.ABL_CONTR))
assert ALL9 == ALL8[:5] + [("PSOBV", "GA")] + ALL8[5:], ALL9
EQ8 = list(X.EQ_PAIRS)
EQ9 = EQ8[:EQ8.index("PSOBV-DE") + 1] + ["PSOBV-GA"] + EQ8[EQ8.index("PSOBV-DE") + 1:]
LABELS = ("tab:X-cluster", "tab:X-loco", "tab:X-equiv-levels", "tab:X-bayes", "tab:X-hl")


def log(msg=""):
    print(msg, flush=True)
    LOG.append(str(msg))


def read(p):
    return open(p, encoding="utf-8").read()


def jround(o):
    return json.loads(json.dumps(X.clean(o)))


def same(x, y, path=""):
    """exact equality of nested JSON-like structures (floats compared exactly; all keys of both)."""
    if isinstance(x, dict):
        assert isinstance(y, dict) and set(x) == set(y), (path, sorted(set(x) ^ set(y)))
        for k in x:
            same(x[k], y[k], path + "/" + str(k))
        return
    if isinstance(x, (list, tuple)):
        assert len(x) == len(y), path
        for i, (u, v) in enumerate(zip(x, y)):
            same(u, v, f"{path}[{i}]")
        return
    assert x == y, (path, x, y)


def must_equal(a, b, what):
    if a != b:
        d = "\n".join(list(difflib.unified_diff(b.splitlines(), a.splitlines(), "expected", "regenerated", lineterm=""))[:40])
        raise AssertionError(f"{what}: regenerated text differs\n{d}")
    log(f"  reproduced exactly: {what}")


def tables_of(text):
    out = {}
    for m in re.finditer(r"\\begin\{(table\*?)\}.*?\\end\{\1\}", text, re.S):
        lab = re.search(r"\\label\{([^}]*)\}", m.group(0))
        out[lab.group(1) if lab else f"nolabel{len(out)}"] = m.group(0)
    return out


def rep1(t, a, b):
    assert t.count(a) == 1, (t.count(a), a)
    return t.replace(a, b)


def rows_of(tab):
    out = []
    for line in tab.splitlines():
        s = line.strip()
        if s.endswith("\\\\") and "&" in s and not s.startswith(("\\multicolumn", "&")):
            out.append([c.strip() for c in s[:-2].split(" & ")])
    return out


def cell_changes(old, new, key_cols=1):
    def keyed(rows):
        out, seen = {}, {}
        for r in rows:
            k = tuple(r[:key_cols]); seen[k] = seen.get(k, 0) + 1
            out[k if seen[k] == 1 else k + (f"(#{seen[k]})",)] = r
        return out
    O, N = keyed(rows_of(old)), keyed(rows_of(new))
    ch, added = [], []
    for k, r in N.items():
        if k not in O:
            added.append(" & ".join(r)); continue
        o = O[k]
        if len(o) != len(r):
            ch.append(dict(row=" / ".join(k), col="(row)", old=" & ".join(o), new=" & ".join(r))); continue
        for i, (u, v) in enumerate(zip(o, r)):
            if u != v:
                ch.append(dict(row=" / ".join(k), col=i, old=u, new=v))
    removed = [" & ".join(r) for k, r in O.items() if k not in N]
    return dict(changed=ch, added_rows=added, removed_rows=removed)


def text_changes(old, new):
    strip = lambda t: [w for line in t.splitlines() if not (line.strip().endswith("\\\\") and "&" in line) for w in line.split()]
    a, b = strip(old), strip(new)
    sm = difflib.SequenceMatcher(a=a, b=b, autojunk=False)
    return [dict(old=" ".join(a[i1:i2]), new=" ".join(b[j1:j2])) for op, i1, i2, j1, j2 in sm.get_opcodes() if op != "equal"]


class Pairs9:
    """ALL_PAIRS / EQ_PAIRS (and, for rank=True, MAIN8 / MAIN_PAIRS) of mpce_inference_extra extended by GA for the
    duration of a call."""
    def __init__(self, pool=False):
        self.pool = pool

    def __enter__(self):
        self.old = (X.MAIN8, X.MAIN_PAIRS, X.ALL_PAIRS, X.EQ_PAIRS, dict(X.LAB), dict(X.MAC))
        X.LAB["GA"] = "GA"; X.MAC["GA"] = "GA"
        X.ALL_PAIRS = list(ALL9); X.EQ_PAIRS = list(EQ9)
        if self.pool:
            X.MAIN8 = list(MAIN9); X.MAIN_PAIRS = list(MAIN_PAIRS9)
        return self

    def __exit__(self, *a):
        X.MAIN8, X.MAIN_PAIRS, X.ALL_PAIRS, X.EQ_PAIRS = self.old[:4]
        X.LAB.clear(); X.LAB.update(self.old[4]); X.MAC.clear(); X.MAC.update(self.old[5])


# ====================================================================== data
def load_ga():
    files = sorted(glob.glob(os.path.join(HERE, "rev2_ga_s*of10.csv")))
    assert len(files) == 10
    ga = pd.concat([pd.read_csv(f, usecols=lambda c: c not in ("Coordinates", "Curve")) for f in files], ignore_index=True)
    ga = X.std_cols(ga)
    assert not ga.duplicated(X.KEY).any()
    ga = ga[(ga.Budget == 6030) & (ga.Init == "random") & ga.Dataset.isin(["1", "2"])].reset_index(drop=True)
    assert len(ga) == 2040 and set(ga.Algorithm) == {"GA"} and sorted(ga.Seed.unique()) == list(range(1, 31))
    return ga


def load_all():
    A = X.load(HERE)
    G8 = A[A.Algorithm.isin(set(MAIN8) | set(X.ABL))]
    ga = load_ga()
    G9 = pd.concat([G8, ga[[c for c in G8.columns if c in ga.columns]]], ignore_index=True)
    assert set(G8.columns) <= set(ga.columns) | {"LossPct"}
    S8 = X.case_stats(G8, sorted(set(MAIN8) | set(X.ABL)), X.BASE_THR)
    S9 = X.case_stats(G9, sorted(set(MAIN9) | set(X.ABL)), X.BASE_THR)
    I = S9.pivot_table(index=X.CASE, columns="Algorithm", values="Ideal")
    assert float(((I["GA"] - I["PSOBV"]).abs() / I["PSOBV"]).max()) < 1e-9
    cases = S8[X.CASE].drop_duplicates()
    clusters = sorted({(d, int(r)) for d, r in zip(cases.Dataset, cases.Radius)})
    return G8, G9, S8, S9, clusters, ga


# ====================================================================== blocks
def blocks(G, S, clusters, pairs):
    full = {}
    for a, b in pairs:
        x = X.cm_test(S, a, b); x["verdict"] = X.verdict(x); full[X.pk(a, b)] = x
    CL = X.cluster_block(S, clusters, pairs)
    LO = X.loco_block(S, clusters, full, pairs)
    EL = X.equivalence_levels_block(G, S, clusters)
    HL = X.hl_block(S, full)
    return dict(case_level=full, cluster=CL, loco=LO, equivalence_levels=EL, hodges_lehmann=HL)


def sub(d, keys):
    return {k: d[k] for k in keys}


# latex_source wording of the five tables (edits found between analysis/ and latex_source/analysis/mpce_supp_inference.tex)
LS_EDITS = {
    "tab:X-equiv-levels": [("(generalization to exchangeable cases, as", "(sensitivity to exchangeable-case assumptions, as"),
                           ("(generalization to new farms; ", "(dependence within the six designed benchmark groups; ")],
    "tab:X-bayes": [("signed-rank test \\cite{Benavoli2017}", "signed-rank test~\\cite{Benavoli2017}")],
}


def ls_edit(lab, t):
    for a, b in LS_EDITS.get(lab, []):
        t = rep1(t, a, b)
    return t


def x_tables(J, G8, G9, S8, S9, clusters):
    stored = json.load(open(os.path.join(HERE, "mpce_summary_extra.json")))
    KEYS = ("case_level", "cluster", "loco", "equivalence_levels", "hodges_lehmann")
    # ---- step 0: recompute the eight-method blocks from the data and compare with the stored JSON
    t0 = time.time()
    B8 = jround(blocks(G8, S8, clusters, ALL8))
    st = copy.deepcopy(stored)
    rep = st["equivalence_levels"].pop("reproduction")
    for k in KEYS:
        same(B8[k], st[k], k)
        log(f"  reproduced exactly: mpce_summary_extra.json {k} (recomputed from the run records)")
    log(f"  ({time.time() - t0:.0f} s)")
    fn = os.path.join(SCR, "supp0.tex"); X.write_tables(stored, fn)
    must_equal(read(fn), read(os.path.join(HERE, "mpce_supp_inference.tex")), "analysis/mpce_supp_inference.tex (all 15 tables)")
    fm = os.path.join(SCR, "macros0.tex"); X.write_macros(stored, fm)
    must_equal(read(fm), read(os.path.join(HERE, "mpce_numbers_extra.tex")), "analysis/mpce_numbers_extra.tex")
    t_old = tables_of(read(fn))
    lsx = tables_of(read(os.path.join(LS, "analysis", "mpce_supp_inference.tex")))
    old = {}
    for lab in LABELS:
        old[lab] = ls_edit(lab, t_old[lab])
        must_equal(old[lab], lsx[lab], f"S-{lab} of latex_source/analysis/mpce_supp_inference.tex (after its wording edits)")
    # ---- GA pair added (LOCO Friedman ranking in the nine-method pool)
    t0 = time.time()
    with Pairs9(pool=True):
        B9 = jround(blocks(G9, S9, clusters, ALL9))
    log(f"  nine-method blocks computed ({time.time() - t0:.0f} s)")
    # existing pairs unchanged (no multiplicity adjustment in these blocks)
    k8 = [X.pk(a, b) for a, b in ALL8]
    same({k: B9["case_level"][k] for k in k8}, st["case_level"], "case_level (8 pairs) with GA")
    same({k: B9["cluster"][k] for k in k8}, st["cluster"], "cluster (old pairs) with GA")
    same({k: B9["loco"]["pairs"][k] for k in k8}, st["loco"]["pairs"], "loco pairs (old pairs) with GA")
    for r9, r8 in zip(B9["loco"]["per_cluster"], st["loco"]["per_cluster"]):
        same({k: r9["tests"][k] for k in k8}, r8["tests"], "loco per-cluster tests")
        assert r9["dropped"] == r8["dropped"] and r9["n_cases"] == r8["n_cases"]
    same({k: B9["hodges_lehmann"][k] for k in k8}, st["hodges_lehmann"], "hodges_lehmann (old pairs) with GA")
    same({k: B9["equivalence_levels"]["pairs"][k] for k in EQ8}, st["equivalence_levels"]["pairs"], "equivalence pairs (old) with GA")
    same({k: v for k, v in B9["equivalence_levels"].items() if k != "pairs"},
         {k: v for k, v in st["equivalence_levels"].items() if k != "pairs"}, "equivalence curve / levels with GA")
    log("  asserted: every value of the existing pairs is unchanged with GA in the pool (only the GA rows are new)")
    # LOCO ranking with GA (caption of tab:X-loco)
    lo = B9["loco"]
    ns_sets = [r["posthoc_nonsig"] for r in lo["per_cluster"]]
    lo["posthoc_only_PSOC_GA_always"] = bool(all(s == ["GA", "PSOC"] for s in ns_sets))
    J["loco_ranking_with_ga"] = [dict(dropped=r["dropped"], best_ranked=r["best_ranked"], posthoc_nonsig=r["posthoc_nonsig"],
                                      avg_rank_PSOBV=r["avg_rank_PSOBV"], avg_rank_PSOC=r["avg_rank_PSOC"]) for r in lo["per_cluster"]]
    J["loco_ranking_8"] = [dict(dropped=r["dropped"], best_ranked=r["best_ranked"], posthoc_nonsig=r["posthoc_nonsig"],
                                avg_rank_PSOBV=r["avg_rank_PSOBV"], avg_rank_PSOC=r["avg_rank_PSOC"]) for r in st["loco"]["per_cluster"]]
    assert lo["best_ranked_always_PSOBV"] and lo["posthoc_only_PSOC_GA_always"], ns_sets
    s1 = copy.deepcopy(stored)
    for k in KEYS:
        s1[k] = B9[k]
    s1["equivalence_levels"]["reproduction"] = rep
    J["ga"] = {"case_level": B9["case_level"]["PSOBV-GA"], "cluster": B9["cluster"]["PSOBV-GA"],
               "loco": B9["loco"]["pairs"]["PSOBV-GA"], "hodges_lehmann": B9["hodges_lehmann"]["PSOBV-GA"],
               "equivalence_levels": B9["equivalence_levels"]["pairs"]["PSOBV-GA"],
               "loco_per_cluster": [dict(dropped=r["dropped"], **r["tests"]["PSOBV-GA"]) for r in B9["loco"]["per_cluster"]]}
    with Pairs9(pool=False):             # MAIN8 stays the 8-method list here: X-modelshift (not output) uses it
        fn1 = os.path.join(SCR, "supp1.tex"); X.write_tables(s1, fn1)
        fm1 = os.path.join(SCR, "macros1.tex"); X.write_macros(s1, fm1)
    t_new = tables_of(read(fn1))
    new = {lab: ls_edit(lab, t_new[lab]) for lab in LABELS}
    # captions: pool-dependent wording
    new["tab:X-loco"] = rep1(new["tab:X-loco"],
                             "The best-ranked method of the eight-method Friedman ranking is PSO-VNS in every LOCO run.",
                             "The best-ranked method of the nine-method Friedman ranking (with GA) is PSO-VNS in every LOCO "
                             "run, and PSO and GA are the only methods not significantly different from it in rank.")
    new["tab:X-equiv-levels"] = rep1(new["tab:X-equiv-levels"], "of the pairs of Table~\\ref{tab:equivalence} at three levels",
                                     "of the pairs of Table~\\ref{tab:equivalence} and of PSO-VNS vs.\\ GA at three levels")
    new["tab:X-bayes"] = rep1(new["tab:X-bayes"], "as Table~\\ref{tab:equivalence}).",
                              "as Table~\\ref{tab:equivalence}; with PSO-VNS vs.\\ GA added).")
    # macros: changed values and new macros
    def mac(text):
        return {m.group(1): m.group(2) for m in re.finditer(r"^\\newcommand\{\\(NX\w+)\}\{(.*)\}(?:   %.*)?$", text, re.M)}
    M0, M1 = mac(read(fm)), mac(read(fm1))
    assert set(M0) <= set(M1), sorted(set(M0) - set(M1))
    chg = {k: (M0[k], M1[k]) for k in M0 if M0[k] != M1[k]}
    newm = {k: M1[k] for k in M1 if k not in M0}
    J["macro_changes"] = {k: dict(old=v[0], new=v[1]) for k, v in chg.items()}
    J["macro_new"] = newm
    # macros of rev5_ga_inference_extra_macros.tex (threshold / multiplicity with GA) must not be touched here
    prev = set(re.findall(r"\\renewcommand\{\\(NX\w+)\}", read(os.path.join(HERE, "rev5_ga_inference_extra_macros.tex"))))
    assert not (set(chg) & prev), sorted(set(chg) & prev)
    J["counts"] = {k: dict(old=old_v, new=new_v) for k, old_v, new_v in
                   [("pairs", len(ALL8), len(ALL9)), ("eq_pairs", len(EQ8), len(EQ9))]}
    return old, new, chg, newm


# ====================================================================== tab:cost
def cost_rows(G, MAINP):
    """the statements of mpce_results.main, section 'computational cost' (copied): table lines and the summary."""
    rt = (G.Seconds / G.Calls * 1000).groupby(G.Algorithm).mean().reindex(MAINP)
    rn = G.groupby(["Algorithm", "Turbines"]).Seconds.mean().unstack().reindex(MAINP)
    me = (G.Seconds / G.Calls * 1000).groupby(G.Algorithm).median().reindex(MAINP)
    lines = [f"{MR.LAB[a]} & {rt[a]:.3f} & {rn.loc[a].min():.1f}--{rn.loc[a].max():.1f} \\\\" for a in MAINP]
    return lines, rt, rn, me


COST_CAP = ("Computational cost at 6,030 objective calls: mean wall-clock time per objective call (including the optimizer "
            "overhead) and range over $N$ of the mean time per run (single core). Times come from the machine that produced "
            "each result file.")


def cost_table(J):
    MR.set_focus("PSOBV")
    ALL = MR.load(HERE, False)[0]
    R6 = ALL[(ALL.Budget == 6030) & (ALL.Init == "random")]
    G = R6[R6.Dataset.isin(["1", "2"]) & R6.Algorithm.isin(MR.MAIN)].reset_index(drop=True)
    MAINP = [a for a in MR.MAIN if a in set(G.Algorithm)]
    assert MAINP == MAIN8, MAINP
    lines, rt, rn, me = cost_rows(G, MAINP)
    t0 = MR.evals(MR.table("table", COST_CAP, "tab:cost", "lcc", "Method & ms per call & s per run", lines))
    stored = tables_of(read(os.path.join(HERE, "mpce_supplementary.tex")))["tab:cost"]
    must_equal(tables_of(t0)["tab:cost"], stored, "S-tab:cost of analysis/mpce_supplementary.tex")
    must_equal(tables_of(t0)["tab:cost"], tables_of(read(os.path.join(LS, "analysis", "mpce_supplementary.tex")))["tab:cost"],
               "S-tab:cost of latex_source/analysis/mpce_supplementary.tex")
    sj = json.load(open(os.path.join(HERE, "mpce_summary.json")))
    assert {a: round(float(v), 4) for a, v in rt.items()} == sj["cost"]["ms_per_call"]
    assert {a: float(me[a]) for a in MAINP} == sj["cost_per_eval"]["median_ms"]
    log("  reproduced exactly: mpce_summary.json cost.ms_per_call and cost_per_eval.median_ms")
    # GA: stored rev2_ga runs (68 cases, 6,030, random), as the Type I supplement pool of mpce_results.main
    GAX = MR.load_ga_secondary(HERE, {})
    GA68 = GAX[(GAX.Source == "rev2_ga") & GAX.Dataset.isin(["1", "2"]) & (GAX.Budget == 6030) & (GAX.Init == "random")]
    assert len(GA68) == 2040 and GA68.Calls.eq(6030).all() and GA68.Seconds.notna().all()
    SUP = MR.with_ga(MAINP)
    GG = pd.concat([G, GA68], ignore_index=True)
    lines9, rt9, rn9, me9 = cost_rows(GG, SUP)
    for a in MAINP:
        assert lines9[SUP.index(a)] == lines[MAINP.index(a)]
    i = SUP.index("GA")
    lines9[i] = lines9[i].replace("GA & ", "GA$^{a}$ & ", 1)
    t1 = MR.evals(MR.table("table", COST_CAP + " $^{a}$GA: runs of a separate batch (rev2\\_ga), so its times are "
                           "indicative only.", "tab:cost", "lcc", "Method & ms per call & s per run", lines9))
    meta9 = [a for a in SUP if a != "SLSQP"]
    J["cost"] = dict(ga_ms_per_eval_mean=float(rt9["GA"]), ga_ms_per_eval_median=float(me9["GA"]),
                     ga_s_per_run_range=[float(rn9.loc["GA"].min()), float(rn9.loc["GA"].max())],
                     ga_s_per_run_by_N={int(k): float(v) for k, v in rn9.loc["GA"].items()}, ga_runs=int(len(GA68)),
                     meta_median_range_old=[float(me[[a for a in MAINP if a != "SLSQP"]].min()),
                                            float(me[[a for a in MAINP if a != "SLSQP"]].max())],
                     meta_median_range_new=[float(me9[meta9].min()), float(me9[meta9].max())],
                     meta_mean_range_new=[float(rt9[meta9].min()), float(rt9[meta9].max())],
                     ga_line=lines9[i])
    return tables_of(t0)["tab:cost"], tables_of(t1)["tab:cost"]


# ====================================================================== main
def main():
    t00 = time.time()
    J = dict(generated=time.strftime("%Y-%m-%d %H:%M:%S"), script="rev5_ga_infer2.py", pool9=MAIN9,
             all_pairs9=[X.pk(a, b) for a, b in ALL9], eq_pairs9=EQ9)
    outs = ("rev5_ga_inference_extra2_tables.tex", "rev5_ga_inference_extra2_macros.tex", "rev5_ga_cost_tables.tex")
    for f in outs:
        p = os.path.join(HERE, f)
        assert not os.path.exists(p) or open(p).readline().startswith("% Generated by analysis/rev5_ga_infer2.py"), \
            f"{f} exists and was not written by this script"
    log("==== data")
    G8, G9, S8, S9, clusters, ga = load_all()
    log(f"  {len(G8)} runs of the 8 main methods + controls, {len(ga)} GA runs (rev2_ga); clusters {clusters}")
    log("==== S-tab:X-cluster / X-loco / X-equiv-levels / X-bayes / X-hl")
    old, new, mchg, mnew = x_tables(J, G8, G9, S8, S9, clusters)
    log("==== S-tab:cost")
    c_old, c_new = cost_table(J)
    hd = "% Generated by analysis/rev5_ga_infer2.py -- do not edit by hand.\n"
    w = lambda fn, txt: open(os.path.join(HERE, fn), "w").write(txt)
    w("rev5_ga_inference_extra2_tables.tex", hd + "% GA versions of S-tab:X-cluster, S-tab:X-loco, S-tab:X-equiv-levels, "
      "S-tab:X-bayes and S-tab:X-hl (pair PSO-VNS vs GA added; complete replacements of those tables of "
      "mpce_supp_inference.tex, same labels, latex_source wording).\n\n" + "\n\n".join(new[lab] for lab in LABELS) + "\n")
    w("rev5_ga_inference_extra2_macros.tex", hd + "% \\NX... macros of mpce_numbers_extra.tex whose value changes when the pair "
      "PSO-VNS vs GA enters the cluster / LOCO / equivalence / Hodges-Lehmann analyses, and the new PSO-VNS-vs-GA macros; "
      "input after mpce_numbers_extra.tex and rev5_ga_inference_extra_macros.tex.\n"
      + "".join(f"\\renewcommand{{\\{k}}}{{{v[1]}}}   % was: {v[0]}\n" for k, v in mchg.items())
      + "".join(f"\\newcommand{{\\{k}}}{{{v}}}\n" for k, v in mnew.items()))
    w("rev5_ga_cost_tables.tex", hd + "% GA version of S-tab:cost (replaces that table of mpce_supplementary.tex).\n\n" + c_new + "\n")
    J["changes"] = {f"S-{lab}": dict(cells=cell_changes(old[lab], new[lab]), text=text_changes(old[lab], new[lab])) for lab in LABELS}
    J["changes"]["S-tab:cost"] = dict(cells=cell_changes(c_old, c_new), text=text_changes(c_old, c_new))
    for k, v in J["changes"].items():
        log(f"\n---- {k}: {len(v['cells']['changed'])} changed cells, {len(v['cells']['added_rows'])} new rows")
        for r in v["cells"]["added_rows"]:
            log(f"  + {r}")
        for c in v["cells"]["changed"]:
            log(f"  {c['row']} | {c['col']}: {c['old']} -> {c['new']}")
        for c in v["text"]:
            log(f"  text: '{c['old']}' -> '{c['new']}'")
    log("\n---- changed \\NX macros")
    for k, v in mchg.items():
        log(f"  \\{k}: {v[0]} -> {v[1]}")
    log(f"\n---- new \\NX macros ({len(mnew)})")
    for k, v in mnew.items():
        log(f"  \\{k}: {v}")
    log(f"\n  LOCO ranking with GA: {J['loco_ranking_with_ga']}")
    log(f"  cost: {J['cost']}")
    J["seconds"] = round(time.time() - t00, 1)
    json.dump(X.clean(J), open(os.path.join(HERE, "rev5_ga_infer2.json"), "w"), indent=1)
    open(os.path.join(HERE, "rev5_ga_infer2.log"), "w").write("\n".join(LOG) + "\n")
    log(f"done in {J['seconds']} s")


if __name__ == "__main__":
    main()
