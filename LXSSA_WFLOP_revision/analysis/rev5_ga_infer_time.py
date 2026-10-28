"""GA in the supplement tables that still use the eight-method pool: inference sensitivity (S-tab:X-threshold,
S-tab:X-multiplicity), imputation sensitivity (S-tab:S-sa-imputation) and elapsed time (S-tab:S-sa-time,
S-tab:S-sa-time-sites).

Usage (from any directory):  python3 analysis/rev5_ga_infer_time.py
Outputs (analysis/; no existing file is modified, nothing under SWEVO_rev2/ is written, no optimizer is run):
  rev5_ga_inference_extra_tables.tex   S-tab:X-threshold and S-tab:X-multiplicity with GA (complete replacements,
                                       same labels; replace those tables of mpce_supp_inference.tex)
  rev5_ga_inference_extra_macros.tex   \\renewcommand of every \\NX... macro of mpce_numbers_extra.tex whose value changes
                                       with GA in the pool (only those; to be input after mpce_numbers_extra.tex)
  rev5_ga_imputation_tables.tex        S-tab:S-sa-imputation with GA (Holm over 8), label / wording of the latex_source copy
  rev5_ga_time_tables.tex              S-tab:S-sa-time and S-tab:S-sa-time-sites with the GA cells, labels / wording of the
                                       latex_source copies; GA timings of the rev5_ga batch marked as indicative
  rev5_ga_infer_time.json, rev5_ga_infer_time.log

Procedure (the generating functions are imported and called):
  X-threshold     mpce_inference_extra.threshold_block on the per-run data of mpce_inference_extra.load + the GA records
                  rev2_ga_s*of10.csv (68 cases, seeds 1-30, 6,030 evaluations, random; std_cols as all other records),
                  module constant MAIN8 set to the nine-method pool (GA after DE, as tab:friedman68-full) for the call;
                  table text by mpce_inference_extra.write_tables. The pool-dependent conclusion 'only PSO not significantly
                  different in rank' becomes 'only PSO and GA ...' (the main-text conclusion with GA), and 'PSO-VNS vs GA
                  not significant (unadjusted and Holm over 8)' is added; columns Rank GA and PSO-VNS/GA are added.
  X-multiplicity  mpce_inference_extra.multiplicity_block (same patch): PSO-VNS vs GA enters the all-family Holm / BH.
  S-sa-imputation rev3_inference.imputation_block on the case statistics with GA, main family PSO-VNS vs 8 methods;
                  table text by rev3_inference.tables.
  S-sa-time       rev3_sites.block_time / time_table; the GA row (not in the local method lists of block_time /
                  time_table) is computed with the statements of block_time copied into bm_entry / bud_entry / time_row,
                  which are first checked to reproduce every existing cell exactly.
  S-sa-time-sites rev3_sites.block_time with rev3_sites.load_hr extended by rev5_ga_hr_s0of1.csv (GA at Horns Rev 6,030
                  feasible / 120,030 random); table by rev3_sites.time_table.
Step 0 (reproduction): every existing table and JSON block is regenerated with the same functions and asserted to equal
the stored outputs (mpce_summary_extra.json, mpce_supp_inference.tex, mpce_numbers_extra.tex, rev3_inference.json,
rev3_inference_tables.tex, rev3_sites.json, rev3_sites_tables.tex) and, after the label / wording edits of the copies,
the latex_source copies.
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
import rev3_inference as RI                 # noqa: E402
import rev3_sites as RS                     # noqa: E402
import mpce_results as MR                   # noqa: E402

LS = os.path.join(ROOT, "SWEVO_rev2", "latex_source")
SCR = tempfile.mkdtemp(prefix="rev5_ga_infer_time_")
LOG = []
MAIN8 = list(X.MAIN8)
MAIN9 = MAIN8[:6] + ["GA"] + MAIN8[6:]          # PSOBV PSOC SSABV SSA LXSSA DE GA BVNS SLSQP (tab:friedman68-full)
assert MAIN8[5] == "DE"


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


def cell_changes(old, new, key_cols=1, header=None):
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
                ch.append(dict(row=" / ".join(k), col=(header[i] if header and i < len(header) else i), old=u, new=v))
    removed = [" & ".join(r) for k, r in O.items() if k not in N]
    return dict(changed=ch, added_rows=added, removed_rows=removed)


def text_changes(old, new):
    strip = lambda t: [w for line in t.splitlines() if not (line.strip().endswith("\\\\") and "&" in line) for w in line.split()]
    a, b = strip(old), strip(new)
    sm = difflib.SequenceMatcher(a=a, b=b, autojunk=False)
    return [dict(old=" ".join(a[i1:i2]), new=" ".join(b[j1:j2])) for op, i1, i2, j1, j2 in sm.get_opcodes() if op != "equal"]


class Pool9:
    """MAIN8 / MAIN_PAIRS of mpce_inference_extra set to the nine-method pool for the duration of a call."""
    def __enter__(self):
        self.old = (X.MAIN8, X.MAIN_PAIRS, dict(X.LAB), dict(X.MAC))
        X.LAB["GA"] = "GA"; X.MAC["GA"] = "GA"
        X.MAIN8 = list(MAIN9)
        X.MAIN_PAIRS = [(X.FOCUS, b) for b in MAIN9 if b != X.FOCUS]
        return self

    def __exit__(self, *a):
        X.MAIN8, X.MAIN_PAIRS = self.old[:2]
        X.LAB.clear(); X.LAB.update(self.old[2]); X.MAC.clear(); X.MAC.update(self.old[3])


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
    Sp = A[A.Algorithm.isin(X.SPM)]
    Sp = Sp[[(d, r, n) in X.SPLITCASES for d, r, n in zip(Sp.Dataset, Sp.Radius, Sp.Turbines)]]
    ga = load_ga()
    G9 = pd.concat([G8, ga[[c for c in G8.columns if c in ga.columns]]], ignore_index=True)
    assert set(G8.columns) <= set(ga.columns) | {"LossPct"}
    SP = X.case_stats(Sp, X.SPM, X.BASE_THR)
    S8 = X.case_stats(G8, sorted(set(MAIN8) | set(X.ABL)), X.BASE_THR)
    S9 = X.case_stats(G9, sorted(set(MAIN9) | set(X.ABL)), X.BASE_THR)
    I = S9.pivot_table(index=X.CASE, columns="Algorithm", values="Ideal")
    assert float(((I["GA"] - I["PSOBV"]).abs() / I["PSOBV"]).max()) < 1e-9
    return A, G8, G9, Sp, SP, S8, S9, ga


# ====================================================================== X-threshold / X-multiplicity
def threshold_summary(TH):
    """the 'threshold_summary' dict of mpce_inference_extra.main (same statements)."""
    base = TH[X.BASE_THR]
    concl_names = list(base["conclusions"])
    T = X.THRESHOLDS
    return dict(
        conclusions=concl_names,
        unchanged={c: bool(all(TH[t]["conclusions"][c] for t in T)) for c in concl_names},
        all_unchanged=bool(all(TH[t]["all_unchanged"] for t in T)),
        best_ranked={str(t): TH[t]["friedman"]["best_ranked"] for t in T},
        rank_PSOBV={str(t): TH[t]["friedman"]["avg_rank"]["PSOBV"] for t in T},
        rank_PSOC={str(t): TH[t]["friedman"]["avg_rank"]["PSOC"] for t in T},
        p_PSOBV_PSOC={str(t): TH[t]["main_case_mean"]["PSOC"]["p"] for t in T},
        order_identical=bool(all(TH[t]["friedman"]["order"] == base["friedman"]["order"] for t in T)),
        orders={str(t): TH[t]["friedman"]["order"] for t in T})


def ga_conclusions(TH):
    """Conclusions of the paper with GA in the pool (main text: PSO-VNS best-ranked, not significantly different from PSO
    or GA in rank and in the case-mean Wilcoxon test, significantly better than the other six). The eight-method
    conclusion 'posthoc_only_PSOC_nonsig' is kept in the JSON (as evaluated on the nine-method ranking) under
    conclusions_8method_definition."""
    for t, x in TH.items():
        c0 = dict(x["conclusions"])
        fr, mc = x["friedman"], x["main_case_mean"]
        new = {}
        for k, v in c0.items():
            if k == "posthoc_only_PSOC_nonsig":
                new["posthoc_only_PSOC_GA_nonsig"] = fr["posthoc_nonsig"] == ["GA", "PSOC"]
            else:
                new[k] = v
            if k == "PSOBV-PSOC_ns":
                new["PSOBV-GA_ns"] = bool(mc["GA"]["p_holm"] >= 0.05 and mc["GA"]["p"] >= 0.05)
        x["conclusions_8method_definition"] = c0
        x["conclusions"] = new
        x["all_unchanged"] = bool(all(new.values()))
    return TH


CAP_EDITS = [
    ("Unq.: case--method cells of the eight-method comparison that do not qualify",
     "Unq.: case--method cells of the nine-method comparison (with GA) that do not qualify"),
    ("Best-ranked method and average ranks of the Friedman ranking of Table~\\ref{M-tab:friedman68};",
     "Best-ranked method and average ranks of the nine-method Friedman ranking of Table~\\ref{M-tab:friedman68};"),
    ("(with the same imputation rule at each threshold) for PSO-VNS vs.\\ PSO and the component-analysis contrasts",
     "(with the same imputation rule at each threshold) for PSO-VNS vs.\\ PSO, PSO-VNS vs.\\ GA and the component-analysis contrasts"),
    ("only PSO not significantly different in rank; PSO-VNS vs.\\ PSO not significant;",
     "only PSO and GA not significantly different in rank; PSO-VNS vs.\\ PSO and PSO-VNS vs.\\ GA not significant;"),
]


def add_ga_columns(tab, TH):
    """Columns Rank GA (after Rank PSO) and PSO-VNS/GA (after PSO-VNS/PSO), formatted as write_tables formats the
    neighbouring columns (f'{:.2f}', mpce_inference_extra.tp)."""
    # two more columns: the table is scaled to the text width as mpce_inference_extra.table(..., resize=True) does
    tab = rep1(tab, "\\begin{tabular}{" + "c" * 12 + "}", "\\resizebox{\\textwidth}{!}{%\n\\begin{tabular}{" + "c" * 14 + "}")
    tab = rep1(tab, "\\end{tabular}\n\\end{table*}", "\\end{tabular}}\n\\end{table*}")
    rk = " & " + X.stack("Rank", "PSO")
    tab = rep1(tab, rk, rk + " & " + X.stack("Rank", "GA"))
    pp = " & " + X.stack("PSO-VNS/", "PSO")
    tab = rep1(tab, pp, pp + " & " + X.stack("PSO-VNS/", "GA"))
    out = []
    for line in tab.splitlines():
        m = re.match(r"^(\d+)( \(paper\))? & ", line)
        if m:
            t = int(m.group(1)); x = TH[str(t)]
            c = line[:-3].split(" & ")
            assert len(c) == 12 and c[4] == f"{x['friedman']['avg_rank']['PSOC']:.2f}" and c[6] == X.tp(x["main_case_mean"]["PSOC"]["p"])
            c = c[:5] + [f"{x['friedman']['avg_rank']['GA']:.2f}"] + c[5:7] + [X.tp(x["main_case_mean"]["GA"]["p"])] + c[7:]
            line = " & ".join(c) + " \\\\"
        out.append(line)
    return "\n".join(out)


def x_tables(J, G8, G9, Sp, SP, S8, S9):
    stored = json.load(open(os.path.join(HERE, "mpce_summary_extra.json")))
    # ---- step 0: recompute the eight-method blocks from the data and compare with the stored JSON
    TH0 = X.threshold_block(G8, Sp)
    same(jround({str(k): v for k, v in TH0.items()}), stored["threshold"], "threshold")
    log("  reproduced exactly: mpce_summary_extra.json threshold (6 thresholds, recomputed from the run records)")
    same(jround(threshold_summary(TH0)), stored["threshold_summary"], "threshold_summary")
    log("  reproduced exactly: mpce_summary_extra.json threshold_summary")
    MU0 = X.multiplicity_block(S8, G8, SP)
    same(jround(MU0), stored["multiplicity"], "multiplicity")
    log("  reproduced exactly: mpce_summary_extra.json multiplicity (28 tests)")
    fn = os.path.join(SCR, "supp0.tex"); X.write_tables(stored, fn)
    must_equal(read(fn), read(os.path.join(HERE, "mpce_supp_inference.tex")), "analysis/mpce_supp_inference.tex (all 15 tables)")
    t_old = tables_of(read(fn))
    lsx = tables_of(read(os.path.join(LS, "analysis", "mpce_supp_inference.tex")))
    for lab in ("tab:X-threshold", "tab:X-multiplicity"):
        must_equal(t_old[lab], lsx[lab], f"S-{lab} of latex_source/analysis/mpce_supp_inference.tex")
    fm = os.path.join(SCR, "macros0.tex"); X.write_macros(stored, fm)
    must_equal(read(fm), read(os.path.join(HERE, "mpce_numbers_extra.tex")), "analysis/mpce_numbers_extra.tex")
    # ---- GA in the pool
    with Pool9():
        TH1 = X.threshold_block(G9, Sp)
        MU1 = X.multiplicity_block(S9, G9, SP)
        TH1raw = copy.deepcopy(TH1)
        TH1 = ga_conclusions(TH1)
        # leave-one-cluster-out ranking statement of the caption of tab:X-loco (report only; table not regenerated)
        cases = S9[X.CASE].drop_duplicates()
        clusters = sorted({(d, int(r)) for d, r in zip(cases.Dataset, cases.Radius)})
        loco = []
        for c in clusters:
            Sx = S9[~((S9.Dataset == c[0]) & (S9.Radius == c[1]))]
            FR = X.friedman(Sx, X.MAIN8)
            loco.append(dict(dropped=f"{c[0]}-{c[1]}", best_ranked=FR["best_ranked"], posthoc_nonsig=FR["posthoc_nonsig"],
                             avg_rank={a: FR["avg_rank"][a] for a in ("PSOBV", "PSOC", "GA")}))
        J["loco_ranking_with_ga"] = jround(loco)
    s1 = copy.deepcopy(stored)
    s1["threshold"] = jround({str(k): v for k, v in TH1.items()})
    s1["threshold_summary"] = jround(threshold_summary(TH1))
    s1["multiplicity"] = jround(MU1)
    J["threshold_old"] = stored["threshold"]; J["threshold_new"] = s1["threshold"]
    J["threshold_new_raw_conclusions"] = {str(t): jround(TH1raw[t]["conclusions"]) for t in TH1raw}
    J["threshold_summary_old"] = stored["threshold_summary"]; J["threshold_summary_new"] = s1["threshold_summary"]
    J["multiplicity_old"] = stored["multiplicity"]; J["multiplicity_new"] = s1["multiplicity"]
    lab0 = dict(X.LAB); mac0 = dict(X.MAC)
    X.LAB["GA"] = "GA"; X.MAC["GA"] = "GA"
    try:
        fn1 = os.path.join(SCR, "supp1.tex"); X.write_tables(s1, fn1)
        fm1 = os.path.join(SCR, "macros1.tex"); X.write_macros(s1, fm1)
    finally:
        X.LAB.clear(); X.LAB.update(lab0); X.MAC.clear(); X.MAC.update(mac0)
    t_new = tables_of(read(fn1))
    th = t_new["tab:X-threshold"]
    for a, b in CAP_EDITS:
        th = rep1(th, a, b)
    th = add_ga_columns(th, s1["threshold"])
    # the one conclusion that changes with GA (asserted), stated in the caption with the values of the table / JSON
    fails = [(t, c) for t in X.THRESHOLDS for c, v in s1["threshold"][str(t)]["conclusions"].items() if not v]
    assert fails == [(30, "PSOBV-GA_ns")], fails
    g30 = s1["threshold"]["30"]["main_case_mean"]["GA"]
    J["threshold_failing_conclusions"] = fails
    th = rep1(th, "and so do those of the controls in Table~\\ref{tab:X-threshold-controls}.}",
              "and so do those of the controls in Table~\\ref{tab:X-threshold-controls}. With GA in the pool, the only change "
              f"is at threshold 30: PSO-VNS vs.\\ GA becomes significant ($p={X.tp(g30['p']).strip('$')}$; Holm over 8: "
              f"$p={X.tp(g30['p_holm']).strip('$')}$), because GA has fewer than 30 feasible runs in {g30['imputed']} cases, which "
              "enter the test with an imputed maximal difference in favour of PSO-VNS.}")
    mu = t_new["tab:X-multiplicity"]
    # macros: changed values only
    def mac(text):
        return {m.group(1): m.group(2) for m in re.finditer(r"^\\newcommand\{\\(NX\w+)\}\{(.*)\}(?:   %.*)?$", text, re.M)}
    M0, M1 = mac(read(fm)), mac(read(fm1))
    assert set(M0) <= set(M1), sorted(set(M0) - set(M1))
    chg = {k: (M0[k], M1[k]) for k in M0 if M0[k] != M1[k]}
    newm = {k: M1[k] for k in M1 if k not in M0}
    J["macro_changes"] = {k: dict(old=v[0], new=v[1]) for k, v in chg.items()}
    J["macro_new"] = newm
    return t_old["tab:X-threshold"], th, t_old["tab:X-multiplicity"], mu, chg, newm


# ====================================================================== S-sa-imputation
def imputation(J, S8, S9):
    R = json.load(open(os.path.join(HERE, "rev3_inference.json")))
    t0 = RI.tables(R)
    must_equal(t0, read(os.path.join(HERE, "rev3_inference_tables.tex")), "analysis/rev3_inference_tables.tex")
    b5m = RI.imputation_block(S8, "tab:friedman68 (PSO-VNS vs 7 methods)", X.MAIN_PAIRS)
    b5a = RI.imputation_block(S8, "tab:ablation (15 contrasts)", X.ABL_CONTR)
    same(jround(b5m), R["B5"]["main"], "B5.main"); same(jround(b5a), R["B5"]["ablation"], "B5.ablation")
    log("  reproduced exactly: rev3_inference.json B5 (main and ablation families, recomputed from the run records)")
    ls_edit = lambda t: rep1(t, "\\label{tab:S-r3-imputation}", "\\label{tab:S-sa-imputation}")
    old = ls_edit(tables_of(t0)["tab:S-r3-imputation"])
    must_equal(old, tables_of(read(os.path.join(LS, "analysis", "rev3_inference_tables.tex")))["tab:S-sa-imputation"],
               "S-tab:S-sa-imputation of latex_source/analysis/rev3_inference_tables.tex")
    pairs9 = [(X.FOCUS, b) for b in MAIN9 if b != X.FOCUS]
    b5m9 = RI.imputation_block(S9, "tab:friedman68 (PSO-VNS vs 8 methods, incl. GA)", pairs9)
    b5a9 = RI.imputation_block(S9, "tab:ablation (15 contrasts)", X.ABL_CONTR)
    same(jround(b5a9), R["B5"]["ablation"], "B5.ablation with GA in the case statistics")
    R1 = copy.deepcopy(R)
    R1["B5"]["main"] = jround(b5m9)
    lab0 = dict(X.LAB); X.LAB["GA"] = "GA"
    try:
        t1 = RI.tables(R1)
    finally:
        X.LAB.clear(); X.LAB.update(lab0)
    new = tables_of(t1)["tab:S-r3-imputation"]
    new = rep1(new, "PSO-VNS vs.\\ each method, Holm over 7)", "PSO-VNS vs.\\ each method, Holm over 8)")
    new = ls_edit(new)
    J["imputation_old"] = R["B5"]["main"]; J["imputation_new"] = jround(b5m9)
    return old, new, R["B5"]["main"], jround(b5m9)


# ====================================================================== S-sa-time / S-sa-time-sites
def rc(med, ratio):
    """rev3_sites.time_table.rc (local function there), verbatim."""
    if med is None or not np.isfinite(med):
        return "--"
    m = f"{med:.1f}" if med >= 10 else f"{med:.2f}"
    return m if ratio is None or not np.isfinite(ratio) else f"{m} ({ratio:.2f})"


def bm_entry(x, src):
    """the per-method dict of rev3_sites.block_time (benchmark), statements copied."""
    byN = x.groupby("Turbines").Seconds.median()
    msN = x.groupby("Turbines").ms.median()
    return dict(source=src, runs=int(len(x)), median_s_per_run=float(x.Seconds.median()),
                median_s_by_N={int(k): float(v) for k, v in byN.items()},
                median_ms_per_eval=float(x.ms.median()), mean_ms_per_eval=float(x.ms.mean()),
                pooled_ms_per_eval=float(1000 * x.Seconds.sum() / x.Calls.sum()),
                median_ms_by_N={int(k): float(v) for k, v in msN.items()},
                total_cpu_h=float(x.Seconds.sum() / 3600))


def add_ratios(v, pv):
    v["ratio_median_s_to_psovns"] = v["median_s_per_run"] / pv["median_s_per_run"]
    v["ratio_by_N_to_psovns"] = {n: v["median_s_by_N"][n] / pv["median_s_by_N"][n] for n in v["median_s_by_N"]}
    return v


def bud_entry(g, a):
    """median_s / ratio_to_psovns / median_ms_per_eval of one method in one (budget, init) group of block_time."""
    med = g.groupby("Algorithm").Seconds.median()
    return (float(med[a]), float(med[a] / med["PSOBV"]),
            float((1000 * g[g.Algorithm == a].Seconds / g[g.Algorithm == a].Calls).median()))


def time_row(a, bm, bf, label):
    """one row of tab:S-r3-time (rev3_sites.time_table, statements copied)."""
    v = bm[a]
    c = [f"{v['median_ms_per_eval']:.2f}", f"{v['mean_ms_per_eval']:.2f}",
         rc(v["median_s_by_N"].get(5), None), rc(v["median_s_by_N"].get(10), None),
         rc(v["median_s_by_N"].get(15), v["ratio_by_N_to_psovns"].get(15))]
    for k in ("6030_feasible", "30030_random", "120030_random"):
        s_ = bf.get(k, {})
        c.append(rc(s_.get("median_s", {}).get(a), s_.get("ratio_to_psovns", {}).get(a)) if a in s_.get("median_s", {}) else "n/r")
    return f"{label} & " + " & ".join(c) + " \\\\"


GA_REV5 = {"30030_random": "rev5_ga_b30k_s0of1.csv", "120030_random": "rev5_ga_b120k_s0of1.csv",
           "6030_feasible": "rev5_ga_feas_s0of1.csv"}
LARGE = [(ds, r, n) for ds in ("1", "2") for r, n in ((500, 10), (750, 12), (1000, 15))]


def rd_ga(fn):
    d = MR.std_cols(pd.read_csv(os.path.join(HERE, fn), usecols=lambda c: c not in ("Coordinates", "Curve")))
    return d.drop_duplicates(MR.KEY, keep="last").reset_index(drop=True)


def time_tables(J, ga):
    stored = json.load(open(os.path.join(HERE, "rev3_sites.json")))
    # ---- step 0: block_time recomputed and compared with the stored JSON; tables from the stored JSON
    s0 = {}
    RS.block_time(s0, [])
    same(jround(s0["time"]), stored["time"], "time")
    log("  reproduced exactly: rev3_sites.json time (block_time recomputed from the run records)")
    summ0 = dict(stored, time=s0["time"])          # in-memory blocks (integer N keys, as in rev3_sites.main)
    tex0 = []
    RS.time_table(summ0, tex0)
    gen = read(os.path.join(HERE, "rev3_sites_tables.tex"))
    for t in tex0:
        assert t in gen, "time table not in rev3_sites_tables.tex"
    log("  reproduced exactly: tab:S-r3-time and tab:S-r3-time-sites of analysis/rev3_sites_tables.tex")
    T0 = {k: v for t in tex0 for k, v in tables_of(t).items()}
    lsx = tables_of(read(os.path.join(LS, "analysis", "rev3_sites_tables.tex")))

    def ls_time(t):
        t = rep1(t, "\\label{tab:S-r3-time}", "\\label{tab:S-sa-time}")
        t = rep1(t, "Elapsed time per evaluation and per run in the original benchmark records",
                 "Elapsed time per evaluation and per run in the stored benchmark records")
        t = rep1(t, "Elapsed (wall-clock) seconds of a run as stored in the original records",
                 "Elapsed (wall-clock) seconds of a run as stored in the run records")
        m = re.search(r"Total serial compute \(sum of Seconds\): ([\d.]+) CPU-hours for the ([\d{},]+) original records and "
                      r"([\d.]+) for the ([\d{},]+) additional-experiment records; the runs were executed in four parallel worker "
                      r"processes, so the elapsed time of a study was about a quarter of its total\.", t)
        assert m
        return t.replace(m.group(0), f"Aggregate worker-hours (sum of the per-run wall-clock Seconds; process CPU time was not recorded "
                                     f"for the optimization runs): {m.group(1)} for the {m.group(2)} records of the main study and "
                                     f"{m.group(3)} for the {m.group(4)} records of the additional experiments; the runs were executed "
                                     f"in four parallel worker processes, so the elapsed time of a study was about a quarter of its "
                                     f"total. The initialization times in brackets are measured process CPU times.")

    def ls_sites(t):
        t = rep1(t, "\\label{tab:S-r3-time-sites}", "\\label{tab:S-sa-time-sites}")
        return rep1(t, "in the original and additional-experiment records (median seconds; ratio to PSO-VNS)",
                    "in the records of the main study and the additional experiments (median seconds; ratio to PSO-VNS)")
    old_t = ls_time(T0["tab:S-r3-time"]); old_s = ls_sites(T0["tab:S-r3-time-sites"])
    must_equal(old_t, lsx["tab:S-sa-time"], "S-tab:S-sa-time of latex_source/analysis/rev3_sites_tables.tex")
    must_equal(old_s, lsx["tab:S-sa-time-sites"], "S-tab:S-sa-time-sites of latex_source/analysis/rev3_sites_tables.tex")
    # ---- check of the copied statements on every existing method
    T = s0["time"]; bm = T["benchmark"]["methods"]; bf = T["budget_feasinit"]["settings"]
    A = X.load(HERE).assign(ms=lambda d: 1000 * d.Seconds / d.Calls)
    for a, v in bm.items():
        e = add_ratios(bm_entry(A[A.Algorithm == a], v["source"]), bm["PSOBV"])
        same(jround(e), jround(v), f"bm_entry {a}")
    for a in ["PSOBV", "PSOC", "SSABV", "SSA", "LXSSA", "DE", "BVNS", "SLSQP", "LXBV", "RSVNS", "RSDVNS"]:
        assert time_row(a, bm, bf, RS.LAB[a]) in T0["tab:S-r3-time"], a
    isL = lambda d: pd.Series([(x, y, z) in LARGE for x, y, z in zip(d.Dataset, d.Radius, d.Turbines)], index=d.index)
    bud = pd.concat([d for d in (RS.read_orig(e) for e in ("b30k", "b30kp", "b120k", "b120kp", "feas", "feasp")) if d is not None],
                    ignore_index=True)
    bud = bud[bud.Dataset.isin(["1", "2"])]
    BL = pd.concat([A[isL(A)].assign(Source="main 6,030"), bud], ignore_index=True)
    BL = BL[isL(BL)]
    for (b, ini), g in BL.groupby(["Budget", "Init"]):
        k = f"{b}_{ini}"
        for a in bf[k]["median_s"]:
            m, r, ms = bud_entry(g, a)
            assert (m, r, ms) == (bf[k]["median_s"][a], bf[k]["ratio_to_psovns"][a], bf[k]["median_ms_per_eval"][a]), (k, a)
    log("  copied statements (bm_entry, bud_entry, time_row, rc) reproduce every stored benchmark / budget value and row")
    # ---- GA: benchmark (rev2_ga, 68 cases) and six largest cases (rev2_ga at 6,030 random; rev5_ga feas / b30k / b120k)
    gaA = ga.assign(ms=lambda d: 1000 * d.Seconds / d.Calls)
    bm1 = copy.deepcopy(bm)
    bm1["GA"] = add_ratios(bm_entry(gaA, "rev2_ga"), bm["PSOBV"])
    bf1 = copy.deepcopy(bf)
    gl = {"6030_random": gaA[isL(gaA)]}
    for k, fn in GA_REV5.items():
        d = rd_ga(fn); d = d[d.Dataset.isin(["1", "2"])]
        gl[k] = d[isL(d)]
    for k, d in gl.items():
        b, ini = k.split("_")
        assert len(d) == 180 and (d.Budget == int(b)).all() and (d.Init == ini).all() and d.Calls.eq(int(b)).all(), k
        g = pd.concat([BL[(BL.Budget == int(b)) & (BL.Init == ini)], d], ignore_index=True)
        m, r, ms = bud_entry(g, "GA")
        bf1[k]["median_s"]["GA"] = m; bf1[k]["ratio_to_psovns"]["GA"] = r; bf1[k]["median_ms_per_eval"]["GA"] = ms
    # ---- Horns Rev: load_hr + rev5_ga_hr (GA at 6,030 feasible and 120,030 random)
    hr5 = rd_ga("rev5_ga_hr_s0of1.csv")
    hr5 = hr5[(hr5.Algorithm == "GA") & (hr5.Dataset == "HR") & (hr5.Turbines == 16)]
    assert sorted(map(tuple, hr5[["Budget", "Init"]].drop_duplicates().values.tolist())) == [(6030, "feasible"), (120030, "random")]
    orig_load_hr = RS.load_hr
    RS.load_hr = lambda: pd.concat([orig_load_hr(), hr5], ignore_index=True)
    try:
        s1 = {}
        RS.block_time(s1, [])
    finally:
        RS.load_hr = orig_load_hr
    hs0, hs1 = stored["time"]["hr"]["settings"], s1["time"]["hr"]["settings"]
    for k in hs0:                                                       # only GA cells added
        for f in ("median_s", "ratio_to_psovns", "median_ms_per_eval"):
            same({a: v for a, v in hs1[k][f].items() if not (a == "GA" and a not in hs0[k][f])}, hs0[k][f], f"hr {k} {f}")
    same({k: v for k, v in jround(s1["time"]).items() if k not in ("hr", "compute")},
         {k: v for k, v in stored["time"].items() if k not in ("hr", "compute")}, "time without hr")
    # ---- totals: the 1,360 later GA runs (rev4_galg / rev5_ga; listed separately in the archive section of the
    #      supplement) are stated next to the 11,640 additional-experiment records, not merged into them
    extra = ["rev4_galg_s0of1.csv", "rev4_galgb_s0of1.csv", "rev5_ga_b30k_s0of1.csv", "rev5_ga_b120k_s0of1.csv",
             "rev5_ga_feas_s0of1.csv", "rev5_ga_hr_s0of1.csv", "rev5_ga_spacing_s0of1.csv"]
    dx = pd.concat([pd.read_csv(os.path.join(HERE, f), usecols=["Seconds"]) for f in extra], ignore_index=True)
    ga_extra = dict(files=extra, records=int(len(dx)), worker_h=float(dx.Seconds.sum() / 3600))
    assert ga_extra["records"] == 1360
    cp1 = T["compute"]
    summ1 = dict(stored, time=copy.deepcopy(s1["time"]))
    summ1["time"]["benchmark"]["methods"] = bm1
    summ1["time"]["budget_feasinit"]["settings"] = bf1
    summ1["time"]["hr"]["settings"] = hs1
    summ1["time"]["compute"] = cp1
    tex1 = []
    RS.time_table(summ1, tex1)
    T1 = {k: v for t in tex1 for k, v in tables_of(t).items()}
    # GA row in tab:S-r3-time, after DE (order of tab:friedman68-full); GA timings of the rev5_ga batch marked
    t = T1["tab:S-r3-time"]
    de = time_row("DE", bm1, bf1, "DE")
    garow = time_row("GA", bm1, bf1, "GA$^{b}$")
    t = rep1(t, de + "\n", de + "\n" + garow + "\n")
    t = rep1(t, "$^{a}$Same method, runs of the preliminary comparison (fresh\\_grid). n/r: not run.",
             "$^{a}$Same method, runs of the preliminary comparison (fresh\\_grid). $^{b}$GA: 68-case runs of a separate run "
             "batch (rev2\\_ga); its six-largest-case runs with feasible initialization, 30{,}030 and 120{,}030 evaluations "
             "(rev5\\_ga) were made later on a machine shared with another job, so these three times are inflated and "
             "indicative only. n/r: not run.")
    new_t = ls_time(t)
    o = (f"{cp1['additional_cpu_h']:.1f} for the {cp1['additional_records']:,} records of the additional experiments;").replace(",", "{,}")
    n_x = f"{ga_extra['records']:,}".replace(",", "{,}")
    new_t = rep1(new_t, o, o[:-1] + f" and {ga_extra['worker_h']:.1f} for the {n_x} later GA runs (rev4\\_galg, rev5\\_ga);")
    s = T1["tab:S-r3-time-sites"]
    lines = s.split("\n")
    for i, line in enumerate(lines):
        if line.startswith("GA & "):
            c = line[:-3].split(" & ")
            assert c[2] != "n/r" and c[4] != "n/r"
            c[2] += "$^{a}$"; c[4] += "$^{a}$"
            lines[i] = " & ".join(c) + " \\\\"
    s = "\n".join(lines)
    s = rep1(s, "Batches: Horns Rev~1 mpce\\_hrfix (ten methods, one batch) and rev2\\_gahr (GA);",
             "Batches: Horns Rev~1 mpce\\_hrfix (ten methods, one batch), rev2\\_gahr (GA, random initialization, 6{,}030 and "
             "30{,}030) and rev5\\_ga\\_hr (GA, 6{,}030 F and 120{,}030; $^{a}$run on a machine shared with another job, so "
             "these two times are inflated and indicative only);")
    new_s = ls_sites(s)
    J["time"] = dict(ga_benchmark=jround(bm1["GA"]), ga_budget_feasinit={k: dict(median_s=bf1[k]["median_s"]["GA"],
                     ratio_to_psovns=bf1[k]["ratio_to_psovns"]["GA"], median_ms_per_eval=bf1[k]["median_ms_per_eval"]["GA"]) for k in bf1},
                     ga_hr={k: dict(median_s=hs1[k]["median_s"].get("GA"), ratio_to_psovns=hs1[k]["ratio_to_psovns"].get("GA"),
                                    median_ms_per_eval=hs1[k]["median_ms_per_eval"].get("GA")) for k in hs1},
                     compute=T["compute"], later_ga_runs=ga_extra,
                     meta_ms_per_eval_median_range_old=[min(bm[a]["median_ms_per_eval"] for a in MAIN8 if a != "SLSQP"),
                                                        max(bm[a]["median_ms_per_eval"] for a in MAIN8 if a != "SLSQP")],
                     meta_ms_per_eval_median_range_new=[min(bm1[a]["median_ms_per_eval"] for a in MAIN9 if a != "SLSQP"),
                                                        max(bm1[a]["median_ms_per_eval"] for a in MAIN9 if a != "SLSQP")],
                     meta_ms_per_eval_mean_range_new=[min(bm1[a]["mean_ms_per_eval"] for a in MAIN9 if a != "SLSQP"),
                                                      max(bm1[a]["mean_ms_per_eval"] for a in MAIN9 if a != "SLSQP")])
    return old_t, new_t, old_s, new_s


# ====================================================================== main
HDR = {"thr": ["Thresh.", "Unq.", "Best", "Rank PSO-VNS", "Rank PSO", "pz ns", "PSO-VNS/PSO", "SSA-VNS/RS-VNS",
               "LX-SSA-VNS/RS-VNS", "SSA-VNS/LX-SSA-VNS", "PSO-VNS/RS-VNS", "Same"],
       "mult": ["Test", "n", "W/L", "dL", "p", "p_Holm", "q_BH", "Holm sig."],
       "imp": ["Pair", "Only 1st/2nd", "Zeros", "Imputed Wilcoxon", "Imputed Pratt", "n", "Dropped Wilcoxon", "Dropped Pratt", "Verdicts"],
       "time": ["Method", "ms/eval median", "ms/eval mean", "N=5", "N=10", "N=15", "6,030 feas. init.", "30,030", "120,030"],
       "sites": ["Method", "HR 6,030 R", "HR 6,030 F", "HR 30,030", "HR 120,030", "IEA 16/6,030", "IEA 16/30,030",
                 "IEA 36/6,030", "IEA 36/30,030"]}


def main():
    t00 = time.time()
    J = dict(generated=time.strftime("%Y-%m-%d %H:%M:%S"), script="rev5_ga_infer_time.py", pool9=MAIN9)
    outs = ("rev5_ga_inference_extra_tables.tex", "rev5_ga_inference_extra_macros.tex", "rev5_ga_imputation_tables.tex",
            "rev5_ga_time_tables.tex")
    for f in outs:
        p = os.path.join(HERE, f)
        assert not os.path.exists(p) or open(p).readline().startswith("% Generated by analysis/rev5_ga_infer_time.py"), \
            f"{f} exists and was not written by this script"
    log("==== data")
    A, G8, G9, Sp, SP, S8, S9, ga = load_all()
    log(f"  {len(G8)} runs of the 8 main methods + controls, {len(ga)} GA runs (rev2_ga), {len(Sp)} split runs")
    log("==== S-tab:X-threshold / S-tab:X-multiplicity")
    th_old, th_new, mu_old, mu_new, mchg, mnew = x_tables(J, G8, G9, Sp, SP, S8, S9)
    log("==== S-tab:S-sa-imputation")
    im_old, im_new, b5_old, b5_new = imputation(J, S8, S9)
    log("==== S-tab:S-sa-time / S-tab:S-sa-time-sites")
    tt_old, tt_new, ts_old, ts_new = time_tables(J, ga)
    hd = "% Generated by analysis/rev5_ga_infer_time.py -- do not edit by hand.\n"
    w = lambda fn, txt: open(os.path.join(HERE, fn), "w").write(txt)
    w("rev5_ga_inference_extra_tables.tex", hd + "% GA versions of S-tab:X-threshold and S-tab:X-multiplicity (nine-method pool; "
      "complete replacements of those tables of mpce_supp_inference.tex, same labels).\n\n" + th_new + "\n\n" + mu_new + "\n")
    w("rev5_ga_inference_extra_macros.tex", hd + "% \\NX... macros of mpce_numbers_extra.tex whose value changes with GA in the main "
      "pool (S-tab:X-threshold, S-tab:X-multiplicity); input after mpce_numbers_extra.tex.\n"
      + "".join(f"\\renewcommand{{\\{k}}}{{{v[1]}}}   % was: {v[0]}\n" for k, v in mchg.items())
      + "".join(f"\\newcommand{{\\{k}}}{{{v}}}\n" for k, v in mnew.items()))
    w("rev5_ga_imputation_tables.tex", hd + "% GA version of S-tab:S-sa-imputation (main family PSO-VNS vs 8 methods, Holm over 8; "
      "replaces that table of rev3_inference_tables.tex).\n\n" + im_new + "\n")
    w("rev5_ga_time_tables.tex", hd + "% GA versions of S-tab:S-sa-time and S-tab:S-sa-time-sites (replace those tables of "
      "rev3_sites_tables.tex); labels and wording of the latex_source versions.\n\n" + tt_new + "\n\n" + ts_new + "\n")
    J["changes"] = {
        "S-tab:X-threshold": dict(cells=cell_changes(th_old, th_new, 1, None), text=text_changes(th_old, th_new)),
        "S-tab:X-multiplicity": dict(cells=cell_changes(mu_old, mu_new, 1, HDR["mult"]), text=text_changes(mu_old, mu_new)),
        "S-tab:S-sa-imputation": dict(cells=cell_changes(im_old, im_new, 1, HDR["imp"]), text=text_changes(im_old, im_new)),
        "S-tab:S-sa-time": dict(cells=cell_changes(tt_old, tt_new, 1, HDR["time"]), text=text_changes(tt_old, tt_new)),
        "S-tab:S-sa-time-sites": dict(cells=cell_changes(ts_old, ts_new, 1, HDR["sites"]), text=text_changes(ts_old, ts_new))}
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
    for k, v in mnew.items():
        log(f"  new \\{k}: {v}")
    log("\n---- verdicts")
    for t in X.THRESHOLDS:
        x = J["threshold_new"][str(t)]
        log(f"  thr {t}: best {x['friedman']['best_ranked']}, order {x['friedman']['order']}, ns {x['friedman']['posthoc_nonsig']}, "
            f"conclusions {x['conclusions']}")
    log(f"  imputation verdict changes (main, GA): {b5_new['verdict_changes']}; old: {b5_old['verdict_changes']}")
    log(f"  LOCO ranking with GA: {J['loco_ranking_with_ga']}")
    J["seconds"] = round(time.time() - t00, 1)
    json.dump(X.clean(J), open(os.path.join(HERE, "rev5_ga_infer_time.json"), "w"), indent=1)
    open(os.path.join(HERE, "rev5_ga_infer_time.log"), "w").write("\n".join(LOG) + "\n")
    log(f"done in {J['seconds']} s")


if __name__ == "__main__":
    main()
