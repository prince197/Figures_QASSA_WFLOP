"""GA in the supplement tables that still exclude it: Lillgrund site table and 1-deg / PyWake re-evaluation, and the
all-run endpoint tables (benchmark pooled and seed level, Lillgrund paired outcome and reliability).

Usage (from any directory):  python3 analysis/rev5_ga_supp_tables.py [--recompute-check]
Outputs (analysis/; no existing file is modified, nothing under SWEVO_rev2/ is written, no optimizer is run):
  rev5_ga_sites_reeval.csv               the 60 GA Lillgrund layouts (rev4_galg / rev4_galgb) re-evaluated exactly as
                                         rev3_sites.reeval_layouts re-evaluates the other Lillgrund layouts (5-deg replay,
                                         1-deg bins centred at 0.5 deg, PyWake 2.6.20 NOJ k = 0.04 at the 1-deg bins);
                                         same columns and float format as rev3_sites_reeval.csv
  rev5_ga_sites_tables.tex               S-tab:S-xp-site and S-tab:S-sa-lg-eval with GA (labels and wording of the
                                         versions in SWEVO_rev2/latex_source)
  rev5_ga_inference_tables.tex           S-tab:S-sa-allrun with GA (Holm over 8)
  rev5_ga_dependence_tables.tex          S-tab:S-sa-dep-bench with GA (Holm over 8)
  rev5_ga_archive_reliability_tables.tex S-tab:archive-reliability and S-tab:archive-paired with GA (Holm over 18);
                                         complete replacement of latex_source/analysis/archive_reliability_tables.tex
  rev5_ga_supp_tables.json, rev5_ga_supp_tables.log

Procedure (the generating functions are imported and called, not re-implemented):
  S-xp-site    rev2_analysis.block_site (tab:R2-site), pool rev2_analysis.SP9 (+ GA after DE); the GA records are added
               through rev2_analysis.read_rev2(..., prefix="rev4"). The printed supplement table is tab:R2-site with a
               hand-edited caption / note; its body rows and Friedman p values are substituted.
  S-sa-lg-eval rev3_sites.block_lg, pool rev3_sites.SP9 (+ GA); the per-layout re-evaluations of the other methods are the
               cached rev3_sites_reeval.csv, those of GA are computed with the code of rev3_sites.reeval_layouts (LG part,
               copied verbatim into lg_rows below); MinSpacing of the GA runs from their records (as add_minspacing).
  S-sa-allrun  rev3_inference.tables on rev3_inference.json, the GA pair and Holm-8 values from rev5_ga_allrun.json.
  S-sa-dep-bench  rev3_dependence.tables on rev3_dependence.json, the GA pair and Holm-8 values from rev5_ga_allrun.json.
  archive      audit_archive.reliability / latex_reliability on rev2_lg16* (+ rev4_galg*), the exact Clopper-Pearson
               column with rev3_dependence.cp (the latex_source version of the file adds this column by hand).
Step 0 (reproduction): every existing table is regenerated with the same functions and asserted to equal, character
for character, the generated file in analysis/ (rev3_sites_tables.tex, rev3_inference_tables.tex,
rev3_dependence_tables.tex, SWEVO_rev2/experiments_docs/rev2_tables.tex, SWEVO_rev2/validation/reliability_tables.tex)
and, after the label / wording edits found between those files and their latex_source copies, the latex_source copy
(and the inline S-xp-site table of SWEVO_supplement.tex). --recompute-check additionally recomputes the 540 stored
Lillgrund re-evaluations of rev3_sites_reeval.csv with lg_rows and asserts equality (slow part, ~1 min).
One process; OMP / BLAS threads = 1.
"""
import os, sys, re, glob, json, time, copy, argparse, difflib
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import mpce_results as MR                   # noqa: E402
import rev2_analysis as R2                  # noqa: E402
import rev3_sites as RS                     # noqa: E402
import rev3_inference as RI                 # noqa: E402
import rev3_dependence as RD                # noqa: E402
import audit_archive as AA                  # noqa: E402

LS = os.path.join(ROOT, "SWEVO_rev2", "latex_source")
DOCS = os.path.join(ROOT, "SWEVO_rev2", "experiments_docs")
VALID = os.path.join(ROOT, "SWEVO_rev2", "validation")
REEVAL_CACHE = os.path.join(HERE, "rev3_sites_reeval.csv")
GA_REEVAL = os.path.join(HERE, "rev5_ga_sites_reeval.csv")
LOG = []


def log(msg=""):
    print(msg, flush=True)
    LOG.append(str(msg))


def read(p):
    return open(p, encoding="utf-8").read()


def tables_of(text):
    """label -> complete table environment (\\begin{table...} ... \\end{table...})."""
    out = {}
    for m in re.finditer(r"\\begin\{(table\*?)\}.*?\\end\{\1\}", text, re.S):
        lab = re.search(r"\\label\{([^}]*)\}", m.group(0))
        out[lab.group(1) if lab else f"nolabel{len(out)}"] = m.group(0)
    return out


def must_equal(a, b, what):
    if a != b:
        d = "\n".join(list(difflib.unified_diff(b.splitlines(), a.splitlines(), "expected", "regenerated", lineterm=""))[:40])
        raise AssertionError(f"{what}: regenerated text differs\n{d}")
    log(f"  reproduced exactly: {what}")


def jround(o):
    return json.loads(json.dumps(RS.clean(o)))


def rows_of(tab):
    """table body rows 'first & ... \\\\' -> list of (key, cells)."""
    out = []
    for line in tab.splitlines():
        s = line.strip()
        if s.endswith("\\\\") and "&" in s and not s.startswith(("\\multicolumn", "&", "Pair", "Method", "Layout", "Budget")):
            cells = [c.strip() for c in s[:-2].split("&")]
            out.append(cells)
    return out


def cell_changes(old, new, key_cols=1, header=None):
    """Cell-level old -> new changes between two versions of a table (rows keyed by their first key_cols cells)."""
    def keyed(rows):                            # a key repeated in a second panel gets a suffix
        out, seen = {}, {}
        for r in rows:
            k = tuple(r[:key_cols]); seen[k] = seen.get(k, 0) + 1
            out[k if seen[k] == 1 else k + (f"(panel {seen[k]})",)] = r
        return out
    O, N = keyed(rows_of(old)), keyed(rows_of(new))
    ch, added = [], []
    for k, r in N.items():
        if k not in O:
            added.append(" & ".join(r)); continue
        o = O[k]
        for i, (u, v) in enumerate(zip(o, r)):
            if u != v:
                ch.append(dict(row=" / ".join(k), col=(header[i] if header and i < len(header) else i), old=u, new=v))
    removed = [" & ".join(r) for k, r in O.items() if k not in N]
    return dict(changed=ch, added_rows=added, removed_rows=removed)


def text_changes(old, new):
    """word-level differences outside the table rows (captions / notes / headers)."""
    strip = lambda t: [w for line in t.splitlines() if not (line.strip().endswith("\\\\") and "&" in line) for w in line.split()]
    a, b = strip(old), strip(new)
    sm = difflib.SequenceMatcher(a=a, b=b, autojunk=False)
    return [dict(old=" ".join(a[i1:i2]), new=" ".join(b[j1:j2])) for op, i1, i2, j1, j2 in sm.get_opcodes() if op != "equal"]


# ====================================================================== Lillgrund re-evaluation (rev3_sites code)
COLS = ["Algorithm", "Dataset", "Turbines", "Seed", "Budget", "Init", "Feasible", "Objective", "MinSpacing", "Seconds", "Coordinates"]


def lg_rows(Y, with_installed):
    """The Lillgrund part of rev3_sites.reeval_layouts, verbatim (same bins, model, centring, columns)."""
    import rev2_site_model as lg
    import py_wake
    from py_wake import NOJ
    from py_wake.examples.data import lillgrund as L
    rows = []
    wd1, f1, p1 = lg.bins(1.0, 0.5)
    lm = NOJ(L.LillgrundSite(), L.SWT23(), k=0.04)
    gx, gy = lg.WT_X[lg.I16].mean(), lg.WT_Y[lg.I16].mean()
    if with_installed:
        xyl, _ = lg.site(16)
        li = dict(**{"5deg_2.5": lg.aep_gwh(xyl), "1deg_0.5": lg.aep_gwh(xyl, True, wd1, f1, p1),
                     "PyWakeNOJ_1deg": float(lm(lg.WT_X[lg.I16], lg.WT_Y[lg.I16], wd=wd1, ws=lg.WS).aep().sum())})
        rows.append(dict(Site="LG", Algorithm="INSTALLED", Seed=0, Budget=0, Init="--", Feasible=True, Objective=li["5deg_2.5"],
                         Source="rev2_site_model.site(16)", **li))
    for r in Y.itertuples(index=False):
        xy = RS.parse_xy(r.Coordinates)
        e = {"5deg_2.5": lg.aep_gwh(xy), "1deg_0.5": lg.aep_gwh(xy, True, wd1, f1, p1),
             "PyWakeNOJ_1deg": float(lm(xy[:, 0] + gx, xy[:, 1] + gy, wd=wd1, ws=lg.WS).aep().sum())}
        rows.append(dict(Site="LG", Algorithm=r.Algorithm, Seed=r.Seed, Budget=r.Budget, Init=r.Init, Feasible=r.Feasible,
                         Objective=r.Objective, Source=r.Source, **e))
    out = pd.DataFrame(rows)
    out["Feasible"] = out.Feasible.astype(str).str.lower().isin(["true", "1", "1.0"])
    out.insert(0, "PyWakeVersion", py_wake.__version__)
    return out


def roundtrip(df, path):
    """write with the float format of rev3_sites_reeval.csv and read back (same rounding as the cached rows)."""
    df.to_csv(path, index=False, float_format="%.10g")
    return pd.read_csv(path)


def ga_layouts():
    Y = pd.concat([RS.raw_shards("rev4_galg_s*of*.csv", COLS).assign(Source="rev4_galg"),
                   RS.raw_shards("rev4_galgb_s*of*.csv", COLS).assign(Source="rev4_galgb")], ignore_index=True)
    Y = Y[Y.Dataset.astype(str) == "LG"].drop_duplicates(["Algorithm", "Seed", "Budget", "Init"], keep="last")
    assert len(Y) == 60 and set(Y.Algorithm) == {"GA"} and sorted(Y.Budget.unique()) == [6030, 30030]
    return Y


def reeval_block(J, recompute_check, scratch):
    E = pd.read_csv(REEVAL_CACHE)
    if recompute_check:
        t0 = time.time()
        cols = COLS
        Y = pd.concat([RS.raw_shards("rev2_lg16_s*of*.csv", cols).assign(Source="rev2_lg16"),
                       RS.raw_shards("rev2_lg16b_s*of*.csv", cols).assign(Source="rev2_lg16b")], ignore_index=True)
        Y = Y[Y.Dataset.astype(str) == "LG"].drop_duplicates(["Algorithm", "Seed", "Budget", "Init"], keep="last")
        chk = roundtrip(lg_rows(Y, True), os.path.join(scratch, "lg_recompute_check.csv"))
        ref = E[E.Site == "LG"].reset_index(drop=True)[list(chk.columns)]
        assert len(chk) == len(ref) == 541
        pd.testing.assert_frame_equal(chk.reset_index(drop=True), ref, check_exact=True, check_dtype=False)
        log(f"  reproduced exactly: the 541 Lillgrund rows (installed + 540 layouts) of rev3_sites_reeval.csv "
            f"recomputed with lg_rows ({time.time() - t0:.0f} s)")
        J["reeval_recompute_check"] = "541 LG rows of rev3_sites_reeval.csv recomputed: identical after the %.10g round trip"
    t0 = time.time()
    G = roundtrip(lg_rows(ga_layouts(), False), GA_REEVAL)
    F = G[G.Feasible]
    d = (F["5deg_2.5"] - F.Objective).abs()
    J["ga_reeval"] = dict(file=os.path.basename(GA_REEVAL), layouts=int(len(G)), feasible=int(len(F)),
                          replay_5deg_max_abs_dev_gwh=float(d.max()) if len(F) else None,
                          seconds=round(time.time() - t0, 1))
    log(f"  GA re-evaluated: {len(G)} layouts ({len(F)} feasible), 5-deg replay max |diff| {d.max():.2e} GWh/yr "
        f"({time.time() - t0:.0f} s) -> {os.path.basename(GA_REEVAL)}")
    return E, G


# ====================================================================== S-sa-lg-eval (rev3_sites.block_lg)
SP10 = ["PSOBV", "PSOC", "SSABV", "SSA", "LXSSA", "DE", "GA", "BVNS", "SLSQP", "RSDVNS"]


def lg_eval_latex_edits(t):
    """edits between analysis/rev3_sites_tables.tex and its latex_source copy (this table)."""
    return t.replace("\\label{tab:S-r3-lg-eval}", "\\label{tab:S-sa-lg-eval}").replace(
        "Table~\\ref{tab:S-rev-site}", "Table~\\ref{tab:S-xp-site}")


def lg_eval_ga_wording(t):
    """the hard-coded counts of block_lg's caption / note for the pool with GA."""
    rep = [("re-evaluated with 1$^\\circ$ direction bins and PyWake NOJ (nine methods,",
            "re-evaluated with 1$^\\circ$ direction bins and PyWake NOJ (ten methods,"),
           ("replayed on the stored coordinates (rounded to 1~cm; maximum absolute difference",
            "replayed on the stored coordinates (rounded to 1~cm, GA: 17 significant digits; maximum absolute difference"),
           ("Holm over its 8 comparisons", "Holm over its 9 comparisons")]
    for a, b in rep:
        assert t.count(a) == 1, a
        t = t.replace(a, b)
    return t


def lg_eval_block(J, E, G):
    gen = read(os.path.join(HERE, "rev3_sites_tables.tex"))
    lsx = tables_of(read(os.path.join(LS, "analysis", "rev3_sites_tables.tex")))
    ref = json.load(open(os.path.join(HERE, "rev3_sites.json")))["lg_eval"]
    # ---- step 0: existing table
    s0, t0 = {}, []
    RS.block_lg(s0, t0, False)
    assert t0[0] in gen, "lg table not in rev3_sites_tables.tex"
    log("  reproduced exactly: tab:S-r3-lg-eval of analysis/rev3_sites_tables.tex")
    assert jround(s0["lg_eval"]) == ref
    log("  reproduced exactly: rev3_sites.json lg_eval")
    old = lg_eval_latex_edits(t0[0].rstrip("\n"))
    must_equal(old, lsx["tab:S-sa-lg-eval"], "S-tab:S-sa-lg-eval of latex_source/analysis/rev3_sites_tables.tex")
    # ---- GA: cached rows + GA rows, GA MinSpacing from its records, pool SP9 + GA
    EG = pd.concat([E, G], ignore_index=True)
    orig = (RS.reeval_layouts, RS.add_minspacing, RS.SP9)

    def add_ms(Ex, site):
        assert site == "LG"
        parts = [R2.read_rev2(e, HERE)[0] for e in ("lg16", "lg16b")] + \
                [R2.read_rev2(e, HERE, prefix="rev4")[0] for e in ("galg", "galgb")]
        Hm = pd.concat(parts, ignore_index=True)[["Algorithm", "Seed", "Budget", "Init", "MinSpacing"]]
        return Ex.merge(Hm, on=["Algorithm", "Seed", "Budget", "Init"], how="left")
    try:
        RS.reeval_layouts = lambda recompute, cache: EG
        RS.add_minspacing = add_ms
        RS.SP9 = SP10
        s1, t1 = {}, []
        RS.block_lg(s1, t1, False)
    finally:
        RS.reeval_layouts, RS.add_minspacing, RS.SP9 = orig
    new = lg_eval_latex_edits(lg_eval_ga_wording(t1[0].rstrip("\n")))
    J["lg_eval_old"] = jround(s0["lg_eval"])
    J["lg_eval_new"] = jround(s1["lg_eval"])
    return old, new, s0["lg_eval"], s1["lg_eval"]


# ====================================================================== S-xp-site (rev2_analysis.block_site)
def supp_xp_site():
    s = read(os.path.join(LS, "SWEVO_supplement.tex"))
    i = s.index("\\label{tab:S-xp-site}")
    a = s.rindex("\\begin{table}", 0, i)
    b = s.index("\\end{table}", i) + len("\\end{table}")
    return s[a:b]


def body_rows(t):
    a = t.index("Installed block")
    b = t.index("\\bottomrule")
    return t[a:b]


def site_block(J):
    orig_read, orig_sp9 = R2.read_rev2, R2.SP9

    def run(with_ga):
        summ, tex = dict(data={}), []
        if with_ga:
            def rd(exp, data_dir, prefix="rev2"):
                d, st = orig_read(exp, data_dir, prefix)
                if prefix == "rev2" and exp in ("lg16", "lg16b"):
                    g, sg = orig_read({"lg16": "galg", "lg16b": "galgb"}[exp], data_dir, "rev4")
                    assert g is not None and len(g) == 30
                    d = pd.concat([d, g], ignore_index=True)
                    st = dict(st, rows=int(len(d)), ga=sg)
                return d, st
            R2.read_rev2, R2.SP9 = rd, SP10
        try:
            R2.block_site(HERE, summ, tex)
        finally:
            R2.read_rev2, R2.SP9 = orig_read, orig_sp9
        return summ, tex[0]
    s0, t0 = run(False)
    docs = read(os.path.join(DOCS, "rev2_tables.tex"))
    assert t0 in docs, "tab:R2-site differs from experiments_docs/rev2_tables.tex"
    log("  reproduced exactly: tab:R2-site of SWEVO_rev2/experiments_docs/rev2_tables.tex")
    ref = json.load(open(os.path.join(DOCS, "rev2_summary.json")))["lillgrund"]
    assert jround(s0["lillgrund"]) == ref
    log("  reproduced exactly: rev2_summary.json lillgrund")
    supp = supp_xp_site()
    assert body_rows(t0) == body_rows(supp), "S-xp-site body differs from tab:R2-site"
    fr0 = {b: R2.fp(s0["lillgrund"]["budgets"][b]["run_level_friedman_p"]).strip("$") for b in (6030, 30030)}
    fsup = f"$p={fr0[6030]}$ (6{{,}}030) and ${fr0[30030]}$ (30{{,}}030)"
    assert supp.count(fsup) == 1, fsup
    log("  reproduced exactly: body rows and Friedman p of the inline S-tab:S-xp-site (SWEVO_supplement.tex)")
    s1, t1 = run(True)
    fr1 = {b: R2.fp(s1["lillgrund"]["budgets"][b]["run_level_friedman_p"]).strip("$") for b in (6030, 30030)}
    new = supp.replace(body_rows(supp), body_rows(t1)).replace(fsup, f"$p={fr1[6030]}$ (6{{,}}030) and ${fr1[30030]}$ (30{{,}}030)")
    J["xp_site_old"] = jround(s0["lillgrund"])
    J["xp_site_new"] = jround(s1["lillgrund"])
    return supp, new, s0["lillgrund"], s1["lillgrund"]


# ====================================================================== archive tables (audit_archive)
ARCH_NOTE_OLD = ("CI: 20,000 paired-seed percentile bootstrap samples, seed $20261004+B$. Two-sided exact sign tests exclude "
                 "ties; Holm correction covers all 16 comparisons in this table. This additional endpoint does not replace the "
                 "original signed-rank tests and gives no ordering between failed layouts.")
ARCH_NOTE_LS = ("95\\% CI: percentile interval of the superiority from 20,000 paired-seed bootstrap samples, seed $20261004+B$; "
                "it is degenerate when no loss or tie occurs (30/0/0: every resample consists of wins only, hence [1.000, 1.000]), "
                "which does not mean that the win probability is known exactly. Exact 95\\% CI: Clopper--Pearson interval of the "
                "conditional win probability $W/(W+L)$ over the discordant seed pairs (independent seeds), defined also at the "
                "boundary (30/0/0: [0.884, 1.000]); in the rows without ties (30,030 evaluations) it also serves as an interval "
                "for the superiority; with ties, $W/(W+L)$ and $(W+T/2)/30$ are different quantities. Two-sided exact sign tests "
                "exclude ties; Holm correction covers all 16 comparisons in this table. This additional endpoint does not replace "
                "the signed-rank tests and gives no ordering between failed layouts.")


def rel_line(r6, r30, label):
    """one row of tab:archive-reliability, the cell code of audit_archive.latex_reliability."""
    cells = [label]
    for r in (r6, r30):
        cells.extend([f"{r['feasible']}/{r['runs']}", f"[{r['wilson95_lower']:.2f}, {r['wilson95_upper']:.2f}]",
                      f"{r['conditional_mean_GWh']:.3f}" if r['conditional_mean_GWh'] is not None else '--'])
    return ' & '.join(cells) + r' \\'


def archive_ls(tex, pairs, nfam):
    """latex_source form: Clopper-Pearson column (rev3_dependence.cp), wording of the copy."""
    t = tex.replace("Section~\\ref{sec:S-r3-precision}", "Section~\\ref{sec:S-sa-precision}")
    t = t.replace(r"\begin{tabular}{rlcccc}", r"\begin{tabular}{rlccccc}").replace(
        r"Budget & Compared with & W/T/L & Superiority & 95\% CI & Holm $p$ \\",
        r"Budget & Compared with & W/T/L & Superiority & 95\% CI & Exact 95\% CI, $W/(W+L)$ & Holm $p$ \\")
    for p in pairs:
        a = (f"{p['budget']:,} & {p['label']} & {p['wins']}/{p['ties']}/{p['losses']} & {p['superiority']:.3f} & "
             f"[{p['bootstrap95_lower']:.3f}, {p['bootstrap95_upper']:.3f}] & ")
        c = RD.cp(p["wins"], p["wins"] + p["losses"])
        assert t.count(a) == 1, a
        t = t.replace(a, a + f"[{c[0]:.3f}, {c[1]:.3f}] & ")
    assert t.count(ARCH_NOTE_OLD.replace("16", str(nfam))) == 1
    return t.replace(ARCH_NOTE_OLD.replace("16", str(nfam)), ARCH_NOTE_LS.replace("16", str(nfam)))


def archive_block(J):
    fr = {}
    for fam, pat in (("lg16", "rev2_lg16_s*of4.csv"), ("lg16b", "rev2_lg16b_s*of30.csv")):
        fr[fam] = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(os.path.join(HERE, pat)))], ignore_index=True)
        assert len(fr[fam]) == AA.FAMILIES[fam]
    rows0, pairs0 = AA.reliability(fr)
    t0 = AA.latex_reliability(rows0, pairs0)
    must_equal(t0, read(os.path.join(VALID, "reliability_tables.tex")), "SWEVO_rev2/validation/reliability_tables.tex")
    by0 = {(r["algorithm"], r["budget"]): r for r in rows0}
    for a in ['PSOBV', 'SLSQP', 'PSOC', 'BVNS', 'SSABV', 'RSDVNS', 'SSA', 'LXSSA', 'DE']:
        assert rel_line(by0[a, 6030], by0[a, 30030], AA.LABELS[a]) in t0.splitlines()
    old = archive_ls(t0, pairs0, 16)
    must_equal(old, read(os.path.join(LS, "analysis", "archive_reliability_tables.tex")),
               "latex_source/analysis/archive_reliability_tables.tex (tab:archive-reliability, tab:archive-paired)")
    # GA
    ga = {"lg16": pd.read_csv(os.path.join(HERE, "rev4_galg_s0of1.csv")), "lg16b": pd.read_csv(os.path.join(HERE, "rev4_galgb_s0of1.csv"))}
    fr2 = {f: pd.concat([fr[f], ga[f][fr[f].columns]], ignore_index=True) for f in fr}
    labels0 = dict(AA.LABELS)
    AA.LABELS["GA"] = "GA"
    try:
        rows1, pairs1 = AA.reliability(fr2)          # Holm over all 18 rows (stored under the key 'holm16_p')
    finally:
        AA.LABELS.clear(); AA.LABELS.update(labels0)
    assert len(pairs1) == 18
    t1 = AA.latex_reliability(rows1, pairs1)
    by1 = {(r["algorithm"], r["budget"]): r for r in rows1}
    de = rel_line(by1["DE", 6030], by1["DE", 30030], "DE")
    assert t1.count(de + "\n") == 1
    t1 = t1.replace(de + "\n", de + "\n" + rel_line(by1["GA", 6030], by1["GA", 30030], "GA") + "\n")
    t1 = t1.replace("Holm correction covers all 16 comparisons in this table", "Holm correction covers all 18 comparisons in this table")
    new = archive_ls(t1, pairs1, 18)
    # cross-check with rev5_ga_allrun.json (Holm-18 values and the GA rows)
    R5 = json.load(open(os.path.join(HERE, "rev5_ga_allrun.json")))["lillgrund"]["all_pairs_holm18"]
    for p, q in zip(pairs1, R5):
        assert (p["budget"], p["second"]) == (q["budget"], q["second"]) and p["holm16_p"] == q["holm18_p"]
        assert (p["wins"], p["ties"], p["losses"], p["superiority"]) == (q["wins"], q["ties"], q["losses"], q["superiority"])
    log("  Holm-18 values and GA rows identical to rev5_ga_allrun.json")
    J["archive_old"] = dict(rows=jround(rows0), pairs=jround(pairs0))
    J["archive_new"] = dict(rows=jround(rows1), pairs=[dict(jround(p), holm18_p=p["holm16_p"]) for p in pairs1])
    return old, new


# ====================================================================== S-sa-allrun / S-sa-dep-bench
def with_main(main):
    old = RI.MAIN7
    RI.MAIN7 = main
    RD.MAIN7 = main
    return old


def allrun_block(J):
    R = json.load(open(os.path.join(HERE, "rev3_inference.json")))
    t0 = RI.tables(R)
    must_equal(t0, read(os.path.join(HERE, "rev3_inference_tables.tex")), "analysis/rev3_inference_tables.tex")
    old = tables_of(t0)["tab:S-r3-allrun"]
    ins = ("pairs treated as independent), Holm-adjusted",
           "pairs treated as independent, which the seeds synchronized across cases do not guarantee; seed-level inference: "
           "Section~\\ref{sec:S-sa-dep}), Holm-adjusted")
    ls_edit = lambda t: t.replace("\\label{tab:S-r3-allrun}", "\\label{tab:S-sa-allrun}").replace(*ins)
    old = ls_edit(old)
    must_equal(old, tables_of(read(os.path.join(LS, "analysis", "rev3_inference_tables.tex")))["tab:S-sa-allrun"],
               "S-tab:S-sa-allrun of latex_source/analysis/rev3_inference_tables.tex")
    G5 = json.load(open(os.path.join(HERE, "rev5_ga_allrun.json")))["benchmark"]
    R1 = copy.deepcopy(R)
    P = R1["B1"]["pairs"]
    P["PSOBV-GA"] = copy.deepcopy(G5["ga_pair"]["pooled"])
    for k, adj in G5["adjusted_main8"].items():
        P[k]["main"]["p_sign_holm_main7"] = adj["pooled_sign_holm"]          # the column printed for the PSO-VNS block
        P[k]["violation"]["p_sign_holm_main7"] = adj["pooled_sign_holm_violation"]
    main8 = RI.MAIN7[:5] + ["GA"] + RI.MAIN7[5:]                              # PSOC SSABV SSA LXSSA DE GA BVNS SLSQP
    assert RI.MAIN7[4] == "DE"
    lab0 = dict(RI.LAB)
    prev = with_main(main8)
    RI.LAB["GA"] = "GA"
    try:
        t1 = RI.tables(R1)
    finally:
        with_main(prev)
        RI.LAB.clear(); RI.LAB.update(lab0)
    new = tables_of(t1)["tab:S-r3-allrun"]
    for a, b in (("Holm-adjusted over the 7 comparisons of PSO-VNS", "Holm-adjusted over the 8 comparisons of PSO-VNS"),
                 ("PSO-VNS vs.\\ each method (Holm over 7)", "PSO-VNS vs.\\ each method (Holm over 8)")):
        assert new.count(a) == 1, a
        new = new.replace(a, b)
    new = ls_edit(new)
    ga_row = [l for l in new.splitlines() if l.startswith("PSO-VNS vs.\\ GA &")][0]
    assert ga_row == json.load(open(os.path.join(HERE, "rev5_ga_allrun.json")))["latex_rows"]["S_sa_allrun"], ga_row
    log("  GA row identical to rev5_ga_allrun.json latex_rows.S_sa_allrun")
    return old, new


def dep_block(J):
    R = json.load(open(os.path.join(HERE, "rev3_dependence.json")))
    t0 = RD.tables(R)
    must_equal(t0, read(os.path.join(HERE, "rev3_dependence_tables.tex")), "analysis/rev3_dependence_tables.tex")
    old = tables_of(t0)["tab:S-sa-dep-bench"]
    must_equal(old, tables_of(read(os.path.join(LS, "analysis", "rev3_dependence_tables.tex")))["tab:S-sa-dep-bench"],
               "S-tab:S-sa-dep-bench of latex_source/analysis/rev3_dependence_tables.tex")
    G5 = json.load(open(os.path.join(HERE, "rev5_ga_allrun.json")))["benchmark"]
    R1 = copy.deepcopy(R)
    P = R1["benchmark"]["pairs"]
    P["PSOBV-GA"] = copy.deepcopy(G5["ga_pair"]["seed"])
    for k, adj in G5["adjusted_main8"].items():
        P[k]["seed"]["p_signflip_holm"] = adj["signflip_holm"]
        P[k]["seed"]["p_wilcoxon_holm"] = adj["wilcoxon_holm"]
        P[k]["group"]["p_signflip_groups_holm"] = adj["groups_holm"]
    main8 = RI.MAIN7[:5] + ["GA"] + RI.MAIN7[5:]
    prev = with_main(main8)
    try:
        t1 = RD.tables(R1)
    finally:
        with_main(prev)
    new = tables_of(t1)["tab:S-sa-dep-bench"]
    for a, b in (("Holm-adjusted over the 7 comparisons of PSO-VNS", "Holm-adjusted over the 8 comparisons of PSO-VNS"),
                 ("PSO-VNS vs.\\ each method (Holm over 7)", "PSO-VNS vs.\\ each method (Holm over 8)")):
        assert new.count(a) == 1, a
        new = new.replace(a, b)
    ga_row = [l for l in new.splitlines() if l.startswith("PSO-VNS vs.\\ GA &")][0]
    assert ga_row == json.load(open(os.path.join(HERE, "rev5_ga_allrun.json")))["latex_rows"]["S_sa_dep_bench"], ga_row
    log("  GA row identical to rev5_ga_allrun.json latex_rows.S_sa_dep_bench")
    return old, new


# ====================================================================== main-text numbers (tab:farm-III-tests)
POOL9 = ["PSOBV", "PSOC", "SSABV", "SSA", "LXSSA", "DE", "GA", "BVNS", "SLSQP"]       # rev4_farm_tables.POOL


def farm3_tests(J, E, G):
    """tab:farm-III-tests: 1-deg / PyWake means and the 1-deg run-level test of PSO-VNS under three Holm families:
    old (SP9, 8 comparisons; as printed), SP9 + GA (9 comparisons; the regenerated S-sa-lg-eval) and the nine-method
    pool of tab:farm-III (no RSD-VNS, 8 comparisons)."""
    EG = pd.concat([E, G], ignore_index=True)
    L = EG[(EG.Site == "LG") & (EG.Algorithm != "INSTALLED")].copy()
    parts = [R2.read_rev2(e, HERE)[0] for e in ("lg16", "lg16b")] + [R2.read_rev2(e, HERE, prefix="rev4")[0] for e in ("galg", "galgb")]
    Hm = pd.concat(parts, ignore_index=True)[["Algorithm", "Seed", "Budget", "Init", "MinSpacing"]]
    X = L.merge(Hm, on=["Algorithm", "Seed", "Budget", "Init"], how="left")
    X["Seed"] = X.Seed.astype(int); X["Budget"] = X.Budget.astype(int)
    assert X.MinSpacing.notna().all()
    out = {}
    for b in RS.BUDS:
        s = X[X.Budget == b]
        means = {}
        for a in SP10:
            f = s[(s.Algorithm == a) & s.Feasible]
            means[a] = dict(feasible=int(len(f)), **{c: (float(f[c].mean()) if len(f) else None)
                                                     for c in ("5deg_2.5", "1deg_0.5", "PyWakeNOJ_1deg")})
        fams = {}
        for name, pool in (("SP9_holm8_printed", RS.SP9), ("SP9+GA_holm9", SP10), ("farmIII_pool9_holm8", POOL9)):
            fams[name] = {c: RS.paper_tests(s[s.Algorithm.isin(pool)], "PSOBV", pool, c) for c in ("1deg_0.5", "PyWakeNOJ_1deg", "5deg_2.5")}
        # GA vs PSO-VNS on jointly feasible seeds (unadjusted Wilcoxon), as tab:S-sa-hr-eval for Horns Rev
        pairs = {c: RS.jf_pair(s, "PSOBV", "GA", c) for c in ("5deg_2.5", "1deg_0.5", "PyWakeNOJ_1deg")}
        # feasibility, PSO-VNS vs GA: exact McNemar (two-sided binomial test on the discordant seed pairs), as for PSO
        from scipy.stats import binomtest
        fa = s[s.Algorithm == "PSOBV"].set_index("Seed").Feasible
        fb = s[s.Algorithm == "GA"].set_index("Seed").Feasible.reindex(fa.index)
        n10, n01 = int((fa & ~fb).sum()), int((~fa & fb).sum())
        mc = dict(psovns_only=n10, ga_only=n01, p_exact=float(binomtest(n10, n10 + n01).pvalue) if n10 + n01 else 1.0)
        out[str(b)] = dict(means=means, tests=fams, psovns_minus_ga_joint_feasible=pairs, mcnemar_psovns_vs_ga=mc)
    J["farm_III_tests"] = jround(out)
    return out


# ====================================================================== output
HDR = {"lg": ["Layout / method", "Feas. 6,030", "5deg 6,030", "1deg 6,030", "PyWake 6,030", "", "Feas. 30,030", "5deg 30,030",
              "1deg 30,030", "PyWake 30,030", "PSO-VNS vs., 1deg (30,030)"],
       "xp": ["Layout / method", "AEP 6,030", "Feas. 6,030", "Loss 6,030", "Dinst 6,030", "Rank 6,030", "AEP 30,030",
              "Feas. 30,030", "Loss 30,030", "Dinst 30,030", "Rank 30,030"],
       "allrun": ["Pair", "Feas. (%)", "dL", "W/T/L", "Score", "Case 95% CI", "Two-stage 95% CI", "p_sign (Holm)", "Score by V"],
       "dep": ["Pair", "W/T/L", "Score", "95% CI", "Seeds +/-", "p_flip (Holm)", "p_W (Holm)", "Groups +/-", "p_grp (Holm)"],
       "rel": ["Method", "Feas. 6,030", "95% CI 6,030", "Mean AEP 6,030", "Feas. 30,030", "95% CI 30,030", "Mean AEP 30,030"],
       "paired": ["Budget", "Compared with", "W/T/L", "Superiority", "95% CI", "Exact 95% CI", "Holm p"]}


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--recompute-check", action="store_true", help="recompute the 541 cached Lillgrund re-evaluations")
    ap.add_argument("--scratch", default=os.environ.get("TMPDIR", "/tmp"))
    a = ap.parse_args(argv)
    t00 = time.time()
    J = dict(generated=time.strftime("%Y-%m-%d %H:%M:%S"), script="rev5_ga_supp_tables.py")
    for f in ("rev5_ga_sites_tables.tex", "rev5_ga_inference_tables.tex", "rev5_ga_dependence_tables.tex",
              "rev5_ga_archive_reliability_tables.tex"):
        assert not os.path.exists(os.path.join(HERE, f)) or open(os.path.join(HERE, f)).readline().startswith(
            "% Generated by analysis/rev5_ga_supp_tables.py"), f"{f} exists and was not written by this script"
    log("==== Lillgrund re-evaluation of the GA layouts")
    E, G = reeval_block(J, a.recompute_check, a.scratch)
    log("==== S-tab:S-sa-lg-eval")
    lg_old, lg_new, lgs0, lgs1 = lg_eval_block(J, E, G)
    log("==== S-tab:S-xp-site")
    xp_old, xp_new, xps0, xps1 = site_block(J)
    log("==== S-tab:archive-reliability / S-tab:archive-paired")
    ar_old, ar_new = archive_block(J)
    log("==== S-tab:S-sa-allrun")
    al_old, al_new = allrun_block(J)
    log("==== S-tab:S-sa-dep-bench")
    de_old, de_new = dep_block(J)
    log("==== main text tab:farm-III-tests")
    f3 = farm3_tests(J, E, G)
    hd = "% Generated by analysis/rev5_ga_supp_tables.py -- do not edit by hand.\n"
    w = lambda fn, txt: open(os.path.join(HERE, fn), "w").write(txt)
    w("rev5_ga_sites_tables.tex", hd + "% GA versions of S-tab:S-xp-site (replaces the inline table of SWEVO_supplement.tex) and "
      "S-tab:S-sa-lg-eval (replaces that table of rev3_sites_tables.tex); labels and wording of the latex_source versions.\n\n"
      + xp_new + "\n\n" + lg_new + "\n")
    w("rev5_ga_inference_tables.tex", hd + "% GA version of S-tab:S-sa-allrun (replaces that table of rev3_inference_tables.tex).\n\n" + al_new + "\n")
    w("rev5_ga_dependence_tables.tex", hd + "% GA version of S-tab:S-sa-dep-bench (replaces that table of rev3_dependence_tables.tex).\n\n" + de_new + "\n")
    w("rev5_ga_archive_reliability_tables.tex", hd + "% GA version of archive_reliability_tables.tex (both tables; complete replacement).\n"
      + ar_new)
    # ---- changes
    old_r = tables_of(ar_old); new_r = tables_of(ar_new)
    J["changes"] = {
        "S-tab:S-xp-site": dict(cells=cell_changes(xp_old, xp_new, 1, HDR["xp"]), text=text_changes(xp_old, xp_new)),
        "S-tab:S-sa-lg-eval": dict(cells=cell_changes(lg_old, lg_new, 1, HDR["lg"]), text=text_changes(lg_old, lg_new)),
        "S-tab:S-sa-allrun": dict(cells=cell_changes(al_old, al_new, 1, HDR["allrun"]), text=text_changes(al_old, al_new)),
        "S-tab:S-sa-dep-bench": dict(cells=cell_changes(de_old, de_new, 1, HDR["dep"]), text=text_changes(de_old, de_new)),
        "S-tab:archive-reliability": dict(cells=cell_changes(old_r["tab:archive-reliability"], new_r["tab:archive-reliability"], 1, HDR["rel"]),
                                          text=text_changes(old_r["tab:archive-reliability"], new_r["tab:archive-reliability"])),
        "S-tab:archive-paired": dict(cells=cell_changes(old_r["tab:archive-paired"], new_r["tab:archive-paired"], 2, HDR["paired"]),
                                     text=text_changes(old_r["tab:archive-paired"], new_r["tab:archive-paired"]))}
    for k, v in J["changes"].items():
        log(f"\n---- {k}: {len(v['cells']['changed'])} changed cells, {len(v['cells']['added_rows'])} new rows")
        for r in v["cells"]["added_rows"]:
            log(f"  + {r}")
        for c in v["cells"]["changed"]:
            log(f"  {c['row']} | {c['col']}: {c['old']} -> {c['new']}")
        for c in v["text"]:
            log(f"  text: '{c['old']}' -> '{c['new']}'")
    # ---- summary numbers quoted in the text
    for tag, s in (("old", lgs0), ("new", lgs1)):
        log(f"lg_eval {tag}: replay {s['replay_5deg']}; above installed {s['above_installed_total']}; mean drop 5->1 "
            f"{s['mean_drop_5to1_pct_all_feasible']:.4f} %")
    for b in ("6030", "30030"):
        r = f3[b]
        log(f"farm-III-tests {b}: GA feasible {r['means']['GA']['feasible']}, 1deg {r['means']['GA']['1deg_0.5']}, PyWake "
            f"{r['means']['GA']['PyWakeNOJ_1deg']}; PSO-VNS 1deg {r['means']['PSOBV']['1deg_0.5']:.4f} PyWake "
            f"{r['means']['PSOBV']['PyWakeNOJ_1deg']:.4f}")
        for fam, t in r["tests"].items():
            log(f"   {fam} 1deg: " + ", ".join(f"{k} {v['p_holm']:.3g} ({v['outcome']})" for k, v in t["1deg_0.5"].items()))
        log(f"   PSO-VNS - GA joint feasible: " + json.dumps(r["psovns_minus_ga_joint_feasible"]))
        log(f"   McNemar PSO-VNS vs GA: {r['mcnemar_psovns_vs_ga']}")
    J["seconds"] = round(time.time() - t00, 1)
    json.dump(RS.clean(J), open(os.path.join(HERE, "rev5_ga_supp_tables.json"), "w"), indent=1)
    open(os.path.join(HERE, "rev5_ga_supp_tables.log"), "w").write("\n".join(LOG) + "\n")
    log(f"done in {J['seconds']} s")


if __name__ == "__main__":
    main()
