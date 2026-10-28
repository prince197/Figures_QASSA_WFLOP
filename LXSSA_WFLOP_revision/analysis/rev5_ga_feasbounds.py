"""GA in the benchmark pool of the two coordinate-precision robustness statements of the manuscript:
(1) worst case of the open labels (rev3_feasbounds.block_bench_worst: every open MS-SLSQP label infeasible, 68-case
    feasibility-aware Friedman ranks) and (2) the replay from the rounded coordinates (feasible layouts re-evaluated
    with the optimization model, 15-deg bins = mpce_direction J1; ranks recomputed; Horns Rev 1 / Lillgrund site means
    from the precision-audit replay).

Usage (from any directory):  python3 analysis/rev5_ga_feasbounds.py [--out-dir DIR]
Outputs (analysis/ by default): rev5_ga_feasbounds.json, rev5_ga_feasbounds_tables.tex (rev3_feasbounds_tables.tex
regenerated with the nine-method benchmark worst case), rev5_ga_feasbounds.log (stdout, if redirected).
No optimizer is run; no existing file is modified.

Step 0 (reproduction, asserted):
  * rev3_feasbounds.block_bench_worst rerun -> equals rev3_feasbounds.json["benchmark_worst_case"] exactly;
    rev3_feasbounds.write_tex(rev3_feasbounds.json) -> byte-identical to rev3_feasbounds_tables.tex.
  * eight-method replay: mpce_direction.with_values(G, layouts, "J1", "I1") + mpce_results.case_stats / rank_matrix /
    friedman_block -> avg ranks of mpce_summary_dir.json["benchmark"]["res"]["rec"] and ["J1"] exactly; the count of
    layouts with |J1 - Objective| > 0.05 pp of the ideal value (233 of 21,716) from mpce_direction_layouts.csv.gz.
Step 1 (GA): GA = rev2_ga_s*of10.csv (mpce_results.std_cols, 68 cases, duplicates on KEY: last kept, as
  mpce_results.load_ga_secondary). GA has no open label (its 3 class-(iii) records were decided by full-precision
  reruns), so the worst case flips the same 448 MS-SLSQP labels. GA's J1 (15-deg replay of its rounded coordinates) is
  recomputed here with mpce_direction.bench_objective(xy, ds, 1) for all its feasible layouts and compared with
  rev5_ga_direction_layouts.csv.gz.
Site means (Horns Rev 1, Lillgrund): per budget / initialization, the mean over the labelled-feasible runs of the stored
  objective and of the objective replayed from the rounded coordinates (rev3_precision_audit_records.csv: ObjectiveStored
  + ReplayDiff; Lillgrund GA rev4_galg*: rev3_precision_audit.audit_file on the 17-digit coordinates).
"""
import os, sys, glob, json, time, argparse, io, contextlib
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import mpce_results as MR                   # noqa: E402
import mpce_direction as D                  # noqa: E402
import rev3_feasbounds as FB                # noqa: E402
import rev3_precision_audit as PA           # noqa: E402

KEY = ["Algorithm", "Dataset", "Radius", "Turbines", "Seed", "Budget", "Init"]


def same(x, y, path=""):
    if isinstance(x, dict):
        assert set(x) == set(y), (path, set(x) ^ set(y))
        for k in x:
            same(x[k], y[k], path + "/" + str(k))
        return
    if isinstance(x, (list, tuple)):
        assert len(x) == len(y), path
        for i, (u, v) in enumerate(zip(x, y)):
            same(u, v, f"{path}[{i}]")
        return
    if isinstance(x, float) and isinstance(y, float) and np.isnan(x) and np.isnan(y):
        return
    assert x == y, (path, x, y)


def load_ga():
    files = sorted(glob.glob(os.path.join(HERE, "rev2_ga_s*of10.csv")))
    assert len(files) == 10
    d = MR.std_cols(pd.concat([pd.read_csv(f) for f in files], ignore_index=True))
    d = d[(d.Algorithm == "GA") & d.Dataset.isin(["1", "2"])].drop_duplicates(KEY, keep="last").assign(Source="rev2_ga")
    assert len(d) == 2040
    return d


# ------------------------------------------------------------------ (1) worst case
def worst_case(o, ALL, ga, MP):
    """rev3_feasbounds.block_bench_worst (benchmark part) with the pool MP; GA rows appended when GA in MP."""
    R6 = ALL[(ALL.Budget == 6030) & (ALL.Init == "random")]
    G = R6[R6.Dataset.isin(["1", "2"]) & R6.Algorithm.isin(MR.MAIN)].reset_index(drop=True)
    if "GA" in MP:
        G = pd.concat([G, ga[[c for c in G.columns if c in ga.columns]]], ignore_index=True)
    cols = ["Algorithm", "Dataset", "Radius", "Turbines", "Seed"]
    os_ = FB.open_keys(o, "mpce_slsqp", "SLSQP").copy()
    os_["Dataset"] = os_.Ds
    os_["Radius"] = os_.Radius.astype(float).astype(int)
    ks = [tuple(x) for x in os_[cols].itertuples(index=False)]
    G1, n1 = FB.flip(G, ks, cols)
    res = {}
    for name, g in (("stored", G), ("msslsqp_open_infeasible", G1)):
        S = MR.case_stats(g, MP)
        R = MR.rank_matrix(S, MP)
        FR = MR.friedman_block(R, MP)
        sl = S[S.Algorithm == "SLSQP"]
        res[name] = dict(avg_rank=FR["avg_rank"], order=sorted(FR["avg_rank"], key=FR["avg_rank"].get),
                         p_holm_vs_focus=FR.get("p_holm_vs_focus") or FR.get("pz_focus"), chi2=FR["chi2"], p=FR["p"],
                         best_cases={a: v for a, v in FR.get("best", {}).items()} if isinstance(FR.get("best"), dict) else None,
                         msslsqp_feasible_pct=float(g[g.Algorithm == "SLSQP"].Feasible.mean() * 100),
                         msslsqp_cases_unqualified=int((~sl.Qualified).sum()),
                         msslsqp_mean_loss_qualified=float(sl[sl.Qualified].Loss.mean()),
                         allrun_psovns_vs_msslsqp=FB.allrun(g, "PSOBV", "SLSQP"))
        if "GA" in MP:
            res[name]["allrun_psovns_vs_ga"] = FB.allrun(g, "PSOBV", "GA")
            res[name]["ga_feasible_pct"] = float(g[g.Algorithm == "GA"].Feasible.mean() * 100)
    st, wc = res["stored"], res["msslsqp_open_infeasible"]
    S0 = MR.case_stats(G, MP); S1 = MR.case_stats(G1, MP)
    rk0 = S0.pivot_table(index=MR.CASE, columns="Algorithm", values="Rank")
    rk1 = S1.pivot_table(index=MR.CASE, columns="Algorithm", values="Rank")
    ch = (rk0 != rk1).any(axis=1)
    n_open_case = os_.groupby(["Dataset", "Radius", "Turbines"]).size()
    res["flipped"] = n1
    res["cases_with_open_msslsqp"] = int(len(n_open_case))
    res["max_open_per_case"] = int(n_open_case.max())
    res["cases_rank_changed"] = int(ch.sum())
    res["order_unchanged"] = st["order"] == wc["order"]
    res["avg_rank_change"] = {a: wc["avg_rank"][a] - st["avg_rank"][a] for a in MP}
    return res


# ------------------------------------------------------------------ (2) replay
def replay_ranks(G, F, MP):
    Gm = G[G.Algorithm.isin(MP)]
    out = {}
    for name, cols in (("rec", None), ("J1", ("J1", "I1"))):
        X = Gm if cols is None else D.with_values(Gm, F, *cols)
        S = MR.case_stats(X[X.Algorithm.isin(MP)], MP)
        FR = MR.friedman_block(MR.rank_matrix(S, MP), MP, focus="PSOBV")
        d = D.case_diff(S, "PSOBV", "PSOC")
        out[name] = dict(avg_rank=FR["avg_rank"], order=sorted(MP, key=lambda a: FR["avg_rank"][a]),
                         p_holm_vs_focus=FR["p_holm_vs_focus"], psovns_minus_pso_pp=float(d.mean()),
                         psovns_minus_pso_cases=int(len(d)))
    out["order_unchanged"] = out["rec"]["order"] == out["J1"]["order"]
    return out


def layouts(fn):
    F = pd.read_csv(os.path.join(HERE, fn))
    F["Dataset"] = F.Dataset.astype(str)
    return F


def pp_dev(F):
    return 100 * (F.J1 - F.Objective).abs() / F.Ideal


def site_means(rec, extra):
    """Per site / budget / init: method means of the labelled-feasible runs, stored and replayed; orders."""
    R = rec[(rec.Role != "duplicate") & rec.Study.isin(["mpce_hrfix", "rev2_gahr", "rev2_lg16", "rev2_lg16b"])
            & (rec.Turbines == 16)].copy()
    R = pd.concat([R, extra], ignore_index=True)
    R = R[R.Feasible]
    R["Site"] = np.where(R.Kind == "lg", "LG", "HR")
    R["Replay"] = R.ObjectiveStored + R.ReplayDiff
    out = {}
    for (site, b, ini), g in R.groupby(["Site", "Budget", "Init"]):
        m = g.groupby("Algorithm")[["ObjectiveStored", "Replay"]].mean()
        o_all = (list(m.ObjectiveStored.sort_values(ascending=False).index), list(m.Replay.sort_values(ascending=False).index))
        mg = m.drop(index="GA", errors="ignore")
        o_wo = (list(mg.ObjectiveStored.sort_values(ascending=False).index), list(mg.Replay.sort_values(ascending=False).index))
        out[f"{site}_{b}_{ini}"] = dict(
            n_feasible={a: int(v) for a, v in g.groupby("Algorithm").size().items()},
            mean_stored={a: float(v) for a, v in m.ObjectiveStored.items()},
            mean_replay={a: float(v) for a, v in m.Replay.items()},
            order_stored=o_all[0], order_replay=o_all[1], order_unchanged=o_all[0] == o_all[1],
            ga_present="GA" in m.index, order_unchanged_without_ga=o_wo[0] == o_wo[1],
            min_gap_stored=float(np.diff(np.sort(m.ObjectiveStored.values)).min()) if len(m) > 1 else None)
    return out


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", default=HERE)
    a = ap.parse_args(argv)
    t0 = time.time()
    OUT = dict(generated=time.strftime("%Y-%m-%d %H:%M:%S"), script="rev5_ga_feasbounds.py")
    REF = json.load(open(os.path.join(HERE, "rev3_feasbounds.json")))
    # ---------------- step 0a: tables reproduction from the stored JSON
    tmp = os.path.join(a.out_dir, "rev5_ga_feasbounds_tables.tex")
    FB.write_tex(REF, tmp)
    assert open(tmp).read() == open(os.path.join(HERE, "rev3_feasbounds_tables.tex")).read(), "write_tex reproduction"
    print("reproduced rev3_feasbounds_tables.tex byte for byte from rev3_feasbounds.json", flush=True)
    # ---------------- step 0b: block_bench_worst reproduction
    MR.set_focus("PSOBV")
    r, f = FB.load_audit()
    o = r[r.Open]
    assert len(o) == 920
    tmpout = {}
    with contextlib.redirect_stdout(io.StringIO()):
        ALL = FB.block_bench_worst(o, tmpout)
    same(json.loads(json.dumps(FB.clean(tmpout["benchmark_worst_case"]))), REF["benchmark_worst_case"], "bench_worst")
    print("reproduced rev3_feasbounds.json benchmark_worst_case exactly", flush=True)
    ga = load_ga()
    ga_open = int(r[(r.Study == "rev2_ga") & r.Open].shape[0])
    assert ga_open == 0
    m8 = worst_case(o, ALL, ga, list(MR.MAIN))
    same(json.loads(json.dumps(FB.clean(m8))), REF["benchmark_worst_case"]["main8_68cases"], "main8 via worst_case")
    MP9 = MR.with_ga(list(MR.MAIN))
    m9 = worst_case(o, ALL, ga, MP9)
    OUT["worst_case"] = dict(pool9=MP9, ga_open_labels=ga_open, main8=m8, main9=m9)
    for nm, res in (("8", m8), ("9", m9)):
        st, wc = res["stored"], res["msslsqp_open_infeasible"]
        print(f"[1] {nm} methods: flipped {res['flipped']}; order unchanged {res['order_unchanged']}; cases with rank change "
              f"{res['cases_rank_changed']}", flush=True)
        print("    stored: " + ", ".join(f"{x} {st['avg_rank'][x]:.3f}" for x in st["order"]), flush=True)
        print("    worst : " + ", ".join(f"{x} {wc['avg_rank'][x]:.3f}" for x in wc["order"]), flush=True)
    st9, wc9 = m9["stored"], m9["msslsqp_open_infeasible"]
    print(f"    all-run PSO-VNS vs GA stored {st9['allrun_psovns_vs_ga']} -> worst {wc9['allrun_psovns_vs_ga']}", flush=True)
    print(f"    all-run PSO-VNS vs MS-SLSQP stored {st9['allrun_psovns_vs_msslsqp']} -> worst {wc9['allrun_psovns_vs_msslsqp']}", flush=True)
    # regenerated table (nine-method benchmark worst case; all other blocks from the stored JSON)
    out9 = json.loads(json.dumps(REF))
    out9["benchmark_worst_case"]["main8_68cases"] = json.loads(json.dumps(FB.clean(m9)))
    FB.write_tex(out9, tmp)
    txt = open(tmp).read()
    assert txt.count("Benchmark (68 cases, eight methods)") == 1
    txt = txt.replace("Benchmark (68 cases, eight methods)", "Benchmark (68 cases, nine methods incl.\\ GA)")
    txt = txt.replace("%% generated by rev3_feasbounds.py -- do not edit by hand",
                      "%% generated by rev5_ga_feasbounds.py (rev3_feasbounds.write_tex; benchmark worst case with GA, "
                      "nine methods; all other blocks from rev3_feasbounds.json) -- do not edit by hand", 1)
    open(tmp, "w").write(txt)
    old = open(os.path.join(HERE, "rev3_feasbounds_tables.tex")).read().splitlines()
    new = txt.splitlines()
    OUT["tables_changed_lines"] = [dict(line=i + 1, old=u, new=v) for i, (u, v) in enumerate(zip(old, new)) if u != v]
    assert len(old) == len(new)
    # ---------------- (2) replay
    G = ALL[(ALL.Budget == 6030) & (ALL.Init == "random") & ALL.Dataset.isin(["1", "2"])
            & ALL.Algorithm.isin(D.ALLM)].reset_index(drop=True)
    G9 = pd.concat([G, ga[[c for c in G.columns if c in ga.columns]]], ignore_index=True)
    F8 = layouts("mpce_direction_layouts.csv.gz")
    F9 = layouts("rev5_ga_direction_layouts.csv.gz")
    SD = json.load(open(os.path.join(HERE, "mpce_summary_dir.json")))["benchmark"]["res"]
    R8 = replay_ranks(G, F8, list(D.MAIN8))
    for nm in ("rec", "J1"):
        same(R8[nm]["avg_rank"], SD[nm]["friedman"]["avg_rank"], "replay8 " + nm)
        assert R8[nm]["order"] == SD[nm]["rank_order"]
    dev8 = pp_dev(F8)
    n8 = int((dev8 > 0.05).sum())
    assert (n8, len(F8)) == (233, 21716), (n8, len(F8))
    print(f"reproduced 8-method replay ranks (mpce_summary_dir.json rec / J1) and 233 of 21,716 layouts > 0.05 pp", flush=True)
    # same layouts for the eleven methods in both files; GA's J1 recomputed from its rounded coordinates
    k = ["Algorithm", "Dataset", "Radius", "Turbines", "Seed"]
    m = F8.merge(F9[F9.Algorithm != "GA"], on=k, suffixes=("", "_9"))
    assert len(m) == len(F8) and all(np.array_equal(m[c].values, m[c + "_9"].values) for c in ("Objective", "J1", "I1"))
    gf = ga[ga.Feasible].copy()
    gf["J1re"], gf["I1re"] = zip(*[D.bench_objective(D.W.parse_coords(c), ds, 1) for c, ds in zip(gf.Coordinates, gf.Dataset)])
    gm = gf.merge(F9[F9.Algorithm == "GA"], on=k, suffixes=("", "_f"))
    assert len(gm) == len(gf) == 2012
    ga_j1_maxdiff = float((gm.J1re - gm.J1).abs().max())
    ga_i1_maxdiff = float((gm.I1re - gm.I1).abs().max())
    print(f"GA J1 recomputed from rounded coordinates: max |diff| vs rev5_ga_direction_layouts {ga_j1_maxdiff:.3g} (ideal {ga_i1_maxdiff:.3g})", flush=True)
    MP9d = list(D.MAIN8[:6]) + ["GA"] + list(D.MAIN8[6:])
    R9 = replay_ranks(G9, F9, MP9d)
    SD9 = json.load(open(os.path.join(HERE, "rev5_ga_summary_dir.json")))["benchmark"]["res"]
    for nm in ("rec", "J1"):
        same(R9[nm]["avg_rank"], SD9[nm]["friedman"]["avg_rank"], "replay9 vs rev5_ga_summary_dir " + nm)
    # the layout files store J1 / I1 rounded to 0.01 objective units: independent full-precision check, J1 / I1 of
    # EVERY feasible layout (eleven methods + GA) recomputed from its rounded run-record coordinates
    Ff = G9[G9.Feasible][k + ["Coordinates", "Objective", "Ideal"]].copy()
    Ff["J1"], Ff["I1"] = zip(*[D.bench_objective(D.W.parse_coords(c), ds, 1) for c, ds in zip(Ff.Coordinates, Ff.Dataset)])
    assert len(Ff) == len(F9)
    R9x = replay_ranks(G9, Ff, MP9d)
    R8x = replay_ranks(G, Ff[Ff.Algorithm != "GA"], list(D.MAIN8))
    r9x_maxdiff = max(abs(R9x["J1"]["avg_rank"][x] - R9["J1"]["avg_rank"][x]) for x in MP9d)
    r8x_maxdiff = max(abs(R8x["J1"]["avg_rank"][x] - SD["J1"]["friedman"]["avg_rank"][x]) for x in D.MAIN8)
    print(f"full-precision J1 of all {len(Ff)} layouts: 9-method order {R9x['J1']['order']} (max |avg rank diff| vs file "
          f"{r9x_maxdiff:.3g}); 8-method max |diff| vs mpce_summary_dir {r8x_maxdiff:.3g}", flush=True)
    dev9 = pp_dev(F9)
    devga = pp_dev(F9[F9.Algorithm == "GA"])
    OUT["replay"] = dict(
        main8=R8, main9=R9, main9_ga_full_precision_j1=R9x, main9_fullprec_max_avg_rank_diff=r9x_maxdiff, main8_full_precision_j1=R8x, main8_fullprec_max_avg_rank_diff_vs_summary=r8x_maxdiff,
        layout_files_rounding="J1 / I1 / Objective stored with 2 decimals in the *_direction_layouts.csv.gz files", ga_j1_recomputed_max_abs_diff=ga_j1_maxdiff, ga_i1_recomputed_max_abs_diff=ga_i1_maxdiff,
        layouts_gt_005pp=dict(eight_plus_components=dict(n=n8, of=int(len(F8))),
                              with_ga=dict(n=int((dev9 > 0.05).sum()), of=int(len(F9))),
                              ga=dict(n=int((devga > 0.05).sum()), of=int(len(devga)), max_pp=float(devga.max()))))
    for nm, R in (("8", R8), ("9", R9)):
        print(f"[2] {nm} methods: order unchanged {R['order_unchanged']}; PSO-VNS - PSO {R['rec']['psovns_minus_pso_pp']:+.4f} -> "
              f"{R['J1']['psovns_minus_pso_pp']:+.4f} pp", flush=True)
        print("    stored  : " + ", ".join(f"{x} {R['rec']['avg_rank'][x]:.3f}" for x in R["rec"]["order"]), flush=True)
        print("    replayed: " + ", ".join(f"{x} {R['J1']['avg_rank'][x]:.3f}" for x in R["J1"]["order"]), flush=True)
    print(f"    layouts > 0.05 pp: {OUT['replay']['layouts_gt_005pp']}", flush=True)
    # ---------------- site means
    rec = pd.read_csv(os.path.join(HERE, "rev3_precision_audit_records.csv"), low_memory=False, dtype={"Dataset": str, "Init": str})
    ex = []
    for st_id, fn in (("rev4_galg", "rev4_galg_s0of1.csv"), ("rev4_galgb", "rev4_galgb_s0of1.csv")):
        ex += PA.audit_file(os.path.join(HERE, fn), st_id, "used")
    ex = pd.DataFrame(ex)
    sm = site_means(rec, ex[[c for c in rec.columns if c in ex.columns]])
    OUT["site_means"] = sm
    for kk, v in sm.items():
        print(f"[2] site {kk}: GA {'in' if v['ga_present'] else 'not in'} pool; order unchanged {v['order_unchanged']} "
              f"(without GA {v['order_unchanged_without_ga']}); order {v['order_replay']}", flush=True)
    OUT["seconds"] = round(time.time() - t0, 1)
    json.dump(FB.clean(OUT), open(os.path.join(a.out_dir, "rev5_ga_feasbounds.json"), "w"), indent=1)
    for c in OUT["tables_changed_lines"]:
        print("changed table line", c["line"], flush=True)
    print(f"done in {OUT['seconds']} s", flush=True)


if __name__ == "__main__":
    main()
