"""Bounds and sensitivity of the strict feasibility labels that the stored (rounded) coordinates cannot decide.

Usage (from analysis/):  python3 rev3_feasbounds.py [--skip-transfer] [--out-dir DIR]
Outputs: rev3_feasbounds.json, rev3_feasbounds_tables.tex (labels tab:S-sa-feasbounds, tab:S-sa-feasworst,
         tab:S-sa-roundtransfer), rev3_feasbounds.log (stdout copy).

Reads only stored files (no optimizer is run): rev3_precision_audit_records.csv (per-record slacks of the rounded
coordinates, written by rev3_precision_audit.py), the merged full-precision rerun files rev3_fullprec_<study>.csv
(rev3_precision_rerun.py collect) and the stored run records. Re-evaluation of stored layouts with the existing
evaluators only.

Open labels: audited records of class (ii)/(iii) whose strict 1e-6 m label the full-precision reruns did not decide
(LabelDecided != True); 920 records, all SLSQP-based and all labelled feasible.

1  Residual bounds. For every open label, from the rounded coordinates (d decimals, eps = 0.5 10^-d m):
     boundary excess e = max over turbines of the signed distance outside the site (circle of the case radius:
     benchmark r = 500 / 750 / 1000 m, IEA37 r = 1300 m (16 turbines) / 2000 m (36); Horns Rev 1 parallelogram of
     hornsrev_model.site(N)), and spacing deficit s = S_min - min pair distance (S_min = 4D = 308 m benchmark,
     2D = 260 m IEA37 (D = 130 m), 4D = 320 m Horns Rev 1);
     guaranteed bounds on the TRUE (unrounded) residuals: e_true <= e + sqrt(2) eps, s_true <= s + 2 sqrt(2) eps
     (+ 1e-9 m float pad, as the audit). The stored MinSpacing column (written with round-trip precision from the
     unrounded layout) gives the true spacing deficit S_min - MinSpacing directly (reported separately).
2  Tolerance sensitivity at tau = 1e-6, 1e-4, 1e-3, 1e-2 m: decided feasible if (raw slack - bound) >= -tau for both
     constraints, decided infeasible if (raw slack + bound) < -tau for either, otherwise undecidable; raw slacks
     S_min - s and -e (m). Rule R: both constraints from the rounded coordinates; rule M: spacing from the stored
     full-precision MinSpacing (bound 1e-9 m), boundary from the rounded coordinates. Also the smallest tolerance
     tau* that decides every label feasible.
3  Worst case for the conclusions: the IEA37 13-method pool (rev3_sites.load_iea / pool_analysis, the paper's
     feasibility-aware rank rule and the best-method seed-paired Wilcoxon / Holm family) recomputed with
     (W1) every open exact-gradient label infeasible, (W2) W1 plus every open MS-SLSQP (finite-difference) label
     infeasible, and (A) an adversarial case that keeps every exact-gradient method qualified: in each setting the
     open runs with the highest AEP of each exact-gradient method are made infeasible until 15 of its 30 runs remain
     feasible (the lowest mean it can have while qualified). Benchmark: the 68-case eight-method comparison
     (mpce_results.case_stats / rank_matrix / friedman_block) and the all-run PSO-VNS vs MS-SLSQP outcome with every
     open MS-SLSQP label infeasible; Horns Rev 1 pool (rev3_sites.load_hr) with the open MS-SLSQP labels infeasible;
     budget / initialization groups: MS-SLSQP qualification before and after.
4  Transfer sensitivity of re-evaluations to coordinate rounding: layouts with a full-precision rerun reproduced bit
     for bit or to floating-point noise (Reproduction in {bit, noise}) that belong to a re-evaluation set
     (benchmark: the feasible 6,030-evaluation random-start runs of the 11 benchmark methods of mpce_direction.py,
     1-deg bins mpce_direction.bench_objective(xy, ds, 15), plus the Gaussian wake at 1 deg; Horns Rev 1: mpce_hrfix +
     rev2_gahr, 1-deg bins centred at 0.5 deg, mpce_direction.hr_aep, and PyWake NOJ 1 deg as rev3_sites.py;
     Lillgrund: rev2_lg16 / rev2_lg16b, 1-deg bins and PyWake NOJ 1 deg as rev3_sites.py): re-evaluated objective
     of the full-precision layout minus that of the stored rounded layout, in pp of the evaluator's ideal value.
     All eligible layouts are evaluated (no subsampling).
"""
import os, sys, re, glob, json, math, time, argparse, warnings
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
warnings.filterwarnings("ignore", category=RuntimeWarning)

TOL = 1e-6
FLOAT_PAD = 1e-9
TAUS = (1e-6, 1e-4, 1e-3, 1e-2)
GRADM = ("SLSQPX", "PSOSLSQPX")
LOG = []


def log(m=""):
    print(m, flush=True)
    LOG.append(str(m))


# ====================================================================== data
def load_audit():
    r = pd.read_csv(os.path.join(HERE, "rev3_precision_audit_records.csv"), low_memory=False,
                    dtype={"Dataset": str, "Init": str})
    r = r[r.Role != "duplicate"].copy()
    fp = [p for p in sorted(glob.glob(os.path.join(HERE, "rev3_fullprec_*.csv")))
          if not re.search(r"_(s\d+of\d+|local)\.csv$", os.path.basename(p))]
    f = pd.concat([pd.read_csv(p, low_memory=False, dtype={"Dataset": str, "Init": str, "Spacing": str}) for p in fp],
                  ignore_index=True)
    assert not f.duplicated(["File", "Row"]).any()
    r = r.merge(f[["File", "Row", "LabelDecided", "Reproduction", "StrictLabel"]], on=["File", "Row"], how="left")
    r["Open"] = (r.Class != "i") & (r.LabelDecided != True)
    r["Eps"] = 0.5 * 10.0 ** (-r.Decimals.astype(float))
    r["Bsp"] = 2 * math.sqrt(2) * r.Eps + FLOAT_PAD
    r["Bbd"] = math.sqrt(2) * r.Eps + FLOAT_PAD
    r["Exc"] = TOL - r.SlackBoundary                     # boundary excess of the rounded layout (m, signed)
    r["Def"] = TOL - r.SlackSpacing                      # spacing deficit of the rounded layout (m, signed)
    r["DefM"] = r.SMin - r.MinSpacingStored              # true spacing deficit from the stored MinSpacing column
    r["Ds"] = r.Dataset.astype(str).str.replace(r"\.0$", "", regex=True)
    return r, f


def decide(g, tau, rule):
    raw_bd, b_bd = -g.Exc, g.Bbd
    if rule == "R":
        raw_sp, b_sp = -g.Def, g.Bsp
    else:
        raw_sp, b_sp = -g.DefM, FLOAT_PAD
    feas = ((raw_sp - b_sp) >= -tau) & ((raw_bd - b_bd) >= -tau)
    inf = ((raw_sp + b_sp) < -tau) | ((raw_bd + b_bd) < -tau)
    return feas, inf & ~feas, ~feas & ~inf


def tau_star(g, rule):
    sp = (g.Def + g.Bsp) if rule == "R" else (g.DefM + FLOAT_PAD)
    return float(np.maximum(sp, g.Exc + g.Bbd).max())


GEOM = {"grid": "circle r = 500/750/1000 m; S_min = 4D = 308 m (D = 77 m)",
        "iea": "circle r = 1300 m (16 turbines) / 2000 m (36); S_min = 2D = 260 m (D = 130 m)",
        "hr": "Horns Rev 1 parallelogram (hornsrev_model.site); S_min = 4D = 320 m (D = 80 m)"}
STUDY_LAB = {"mpce_slsqp": "Benchmark, 68 cases", "mpce_feas": "Feasible init., 6{,}030", "mpce_b30k": "30{,}030 eval.",
             "mpce_b120k": "120{,}030 eval.", "mpce_hrfix": "HR1, all budgets", "mpce_iea": "IEA37",
             "rev2_grad": "IEA37", "fresh_hr": "HR1, earlier binning"}
MULTI_KIND = ("mpce_feas", "mpce_b30k", "mpce_b120k")
ALG_LAB = {"SLSQP": "MS-SLSQP", "SLSQPX": "MS-SLSQP, exact grad.", "PSOSLSQPX": "PSO-SLSQP, exact grad."}
KIND_LAB = {"grid": "bench.", "iea": "IEA37", "hr": "HR1"}


# ====================================================================== 1 + 2
def block_bounds(r, out):
    o = r[r.Open].copy()
    log(f"[1] open labels: {len(o)} (labelled feasible: {int(o.Feasible.sum())}); by study: "
        + ", ".join(f"{k} {v}" for k, v in o.groupby("Study").size().items()))
    rows = []
    for (st, alg, kind), g in o.groupby(["Study", "Algorithm", "Kind"], sort=False):
        rad = g.Radius.astype(float)
        d = dict(study=st, algorithm=alg, kind=kind, role=g.Role.iloc[0], n=int(len(g)),
                 labelled_feasible=int(g.Feasible.sum()), decimals=sorted(int(x) for x in g.Decimals.unique()),
                 eps_m=float(g.Eps.max()), radii_m=sorted(float(x) for x in rad.unique()),
                 smin_m=sorted(float(x) for x in g.SMin.unique()), geometry=GEOM[kind],
                 max_boundary_excess_rounded_m=float(g.Exc.max()),
                 max_spacing_deficit_rounded_m=float(g.Def.max()),
                 bound_true_boundary_excess_m=float((g.Exc + g.Bbd).max()),
                 bound_true_spacing_deficit_m=float((g.Def + g.Bsp).max()),
                 true_spacing_deficit_from_minspacing_m=float(g.DefM.max()),
                 bound_true_boundary_rel_radius=(float(((g.Exc + g.Bbd) / rad).max()) if kind != "hr" else None),
                 bound_true_boundary_rel_smin=float(((g.Exc + g.Bbd) / g.SMin).max()),
                 bound_true_spacing_rel_smin=float(((g.Def + g.Bsp) / g.SMin).max()),
                 n_boundary_active=int((g.Exc.abs() <= g.Bbd).sum()),
                 n_spacing_active_rounded=int((g.Def.abs() <= g.Bsp).sum()),
                 n_spacing_satisfied_minspacing=int((g.DefM <= TOL).sum()),
                 n_rounded_outside=int((g.Exc > TOL).sum()), n_rounded_spacing_short=int((g.Def > TOL).sum()),
                 reproduction=g.Reproduction.value_counts().to_dict())
        rows.append(d)
    out["open_by_study_method"] = rows
    tot = dict(n=int(len(o)), labelled_feasible=int(o.Feasible.sum()),
               max_boundary_excess_rounded_m={k: float(g.Exc.max()) for k, g in o.groupby("Kind")},
               max_spacing_deficit_rounded_m={k: float(g.Def.max()) for k, g in o.groupby("Kind")},
               bound_true_boundary_excess_m={k: float((g.Exc + g.Bbd).max()) for k, g in o.groupby("Kind")},
               bound_true_spacing_deficit_m={k: float((g.Def + g.Bsp).max()) for k, g in o.groupby("Kind")},
               true_spacing_deficit_from_minspacing_m={k: float(g.DefM.max()) for k, g in o.groupby("Kind")},
               n_spacing_satisfied_minspacing=int((o.DefM <= TOL).sum()),
               n_boundary_active=int((o.Exc.abs() <= o.Bbd).sum()))
    out["open_total"] = tot
    for d in rows:
        log(f"  {d['study']:11s} {d['algorithm']:10s} {d['kind']:4s} n={d['n']:4d} d={d['decimals']} exc_r={d['max_boundary_excess_rounded_m']:.6f} "
            f"def_r={d['max_spacing_deficit_rounded_m']:.6f} | bound exc={d['bound_true_boundary_excess_m']:.6f} "
            f"def={d['bound_true_spacing_deficit_m']:.6f} | defM={d['true_spacing_deficit_from_minspacing_m']:.2e} "
            f"bd_active={d['n_boundary_active']} sp_active={d['n_spacing_active_rounded']} spM_ok={d['n_spacing_satisfied_minspacing']}")
    # ---- 2: tolerance sensitivity
    eg = r[r.Study == "rev2_grad"].copy()
    assert len(eg) == 240, len(eg)
    groups = {"open920": o, "open_exact_gradient": o[o.Study == "rev2_grad"], "exact_gradient240": eg,
              "open_benchmark_circle": o[o.Kind == "grid"], "open_iea": o[o.Kind == "iea"], "open_hr": o[o.Kind == "hr"]}
    sens = {}
    for name, g in groups.items():
        sens[name] = dict(n=int(len(g)), tau_star_R=tau_star(g, "R"), tau_star_M=tau_star(g, "M"), tol={})
        for tau in TAUS:
            row = {}
            for rule in ("R", "M"):
                fe, inf, und = decide(g, tau, rule)
                row[rule] = dict(feasible=int(fe.sum()), infeasible=int(inf.sum()), undecidable=int(und.sum()))
                if name == "exact_gradient240":
                    full = g.LabelDecided == True
                    row[rule + "_with_reruns"] = dict(feasible=int((fe | (full & (g.StrictLabel == True))).sum()),
                                                      infeasible=int((inf & ~full).sum()),
                                                      undecidable=int((und & ~full).sum()))
            sens[name]["tol"][f"{tau:g}"] = row
        log(f"[2] {name} (n={len(g)}): tau*_R={sens[name]['tau_star_R']:.6f} m, tau*_M={sens[name]['tau_star_M']:.6f} m; "
            + "; ".join(f"{t}: R {v['R']['feasible']}/{v['R']['infeasible']}/{v['R']['undecidable']} "
                        f"M {v['M']['feasible']}/{v['M']['infeasible']}/{v['M']['undecidable']}"
                        for t, v in sens[name]["tol"].items()))
    out["tolerance_sensitivity"] = sens
    return o


# ====================================================================== 3
def flip(y, keys, cols):
    """Set Feasible = False for the rows of y whose key tuple (cols) is in keys; returns (copy, number flipped)."""
    y = y.copy()
    k = pd.MultiIndex.from_frame(y[cols].astype(str))
    kk = pd.MultiIndex.from_tuples([tuple(map(str, t)) for t in keys], names=cols)
    m = k.isin(kk)
    y.loc[m, "Feasible"] = False
    return y, int(m.sum())


def open_keys(o, study, alg=None):
    g = o[o.Study == study]
    if alg is not None:
        g = g[g.Algorithm.isin(alg if isinstance(alg, (list, tuple)) else [alg])]
    return g


def summarize_pool(res, methods):
    S = {a: dict(feasible=v["feasible"], runs=v["runs"], mean=v["mean"], rank=v["rank"], qualified=v["qualified"])
         for a, v in res["methods"].items() if a in methods}
    fam = res["families"]["best_vs_all"]
    ns = [b for b, t in fam["tests"].items() if t["outcome"] == "T"]
    worse = [b for b, t in fam["tests"].items() if t["outcome"] not in ("T",) and False]
    return dict(best=res["best"], best_mean=res["methods"][res["best"]]["mean"], methods=S,
                best_family=dict(focus=fam["focus"], not_significant=ns,
                                 outcomes={b: t["outcome"] for b, t in fam["tests"].items()},
                                 p_holm={b: t["p_holm"] for b, t in fam["tests"].items()}))


def block_iea_worst(o, out):
    import rev3_sites as S
    IE = S.load_iea().copy()
    IE["Ds"] = IE.Dataset.astype(str)
    cols = ["Algorithm", "Turbines", "Seed", "Budget", "Init"]
    assert not IE.duplicated(cols).any()
    og = open_keys(o, "rev2_grad", list(GRADM))
    os_ = open_keys(o, "mpce_iea", "SLSQP")
    kg = [tuple(x) for x in og[cols].itertuples(index=False)]
    ks = [tuple(x) for x in os_[cols].itertuples(index=False)]
    res = {}
    show = list(GRADM) + ["SLSQP", "PSOBV", "PSOC", "GA", "SSABV", "BVNS"]
    for n in (16, 36):
        for b in S.BUDS:
            y0 = IE[(IE.Turbines == n) & (IE.Budget == b) & (IE.Init == "random")].copy()
            fam = {"best_vs_all": (lambda best: best, None)}
            sc = {}
            y1, n1 = flip(y0, kg, cols)
            y2, n2 = flip(y1, ks, cols)
            # adversarial but qualified: per exact-gradient method, flip its highest-AEP open runs while >= 15 stay feasible
            ya = y0.copy()
            na = {}
            for a in GRADM:
                ko = set(map(lambda t: tuple(map(str, t)), [k for k in kg if k[0] == a and k[1] == n and k[3] == b]))
                ga = ya[(ya.Algorithm == a) & ya.Feasible]
                ga = ga[[tuple(map(str, t)) in ko for t in ga[cols].itertuples(index=False)]].sort_values("Objective", ascending=False)
                nf = int(((ya.Algorithm == a) & ya.Feasible).sum())
                kflip = max(0, min(len(ga), nf - 15))
                ya.loc[ga.index[:kflip], "Feasible"] = False
                na[a] = dict(open=int(len(ga)), feasible_stored=nf, flipped=kflip)
            for name, y, nfl in (("stored", y0, 0), ("W1_exactgrad_open_infeasible", y1, n1),
                                 ("W2_plus_msslsqp_open_infeasible", y2, n1 + n2),
                                 ("A_adversarial_qualified", ya, sum(v["flipped"] for v in na.values()))):
                r = S.pool_analysis(y, S.IEA13, fam)
                free = [a for a in S.FREE if a in r["methods"]]
                bf = sorted(free, key=lambda a: (r["methods"][a]["rank"], -(r["methods"][a]["mean"] or -1e18)))[0]
                d = summarize_pool(r, show)
                d.update(flipped=nfl, best_gradient_free=bf, best_free_mean=r["methods"][bf]["mean"],
                         exact_gradient_first=r["best"] in GRADM,
                         ranks_all={a: v["rank"] for a, v in r["methods"].items()},
                         feasible_all={a: v["feasible"] for a, v in r["methods"].items()},
                         mean_all={a: v["mean"] for a, v in r["methods"].items()})
                if name == "A_adversarial_qualified":
                    d["adversarial_detail"] = na
                sc[name] = d
            res[f"{n}T_{b}"] = sc
            for name, d in sc.items():
                log(f"[3] IEA {n}T {b} {name:32s} flipped={d['flipped']:3d} best={d['best']} ({d['best_mean']:.0f}) "
                    f"bestfree={d['best_gradient_free']} ({d['best_free_mean']:.0f}) "
                    + " ".join(f"{a}:{d['methods'][a]['feasible']}/{d['methods'][a]['rank']:g}/{(d['methods'][a]['mean'] or 0):.0f}"
                               for a in list(GRADM) + ["SLSQP", "PSOBV"])
                    + f" | n.s. vs best: {d['best_family']['not_significant']}")
    # consistency with the published pool (rev3_sites.json)
    try:
        J = json.load(open(os.path.join(HERE, "rev3_sites.json")))["iea_pool"]["settings"]
        chk = all(J[k]["best"] == v["stored"]["best"] and all(
            J[k]["methods"][a]["feasible"] == v["stored"]["feasible_all"][a] and J[k]["methods"][a]["rank"] == v["stored"]["ranks_all"][a]
            for a in v["stored"]["ranks_all"]) for k, v in res.items())
    except Exception as e:                      # noqa
        chk = f"not checked: {e!r}"
    log(f"[3] stored scenario reproduces rev3_sites.json iea_pool: {chk}")
    out["iea_worst_case"] = dict(settings=res, reproduces_rev3_sites=chk,
                                 open_exact_gradient=int(len(og)), open_msslsqp_iea=int(len(os_)),
                                 open_by_setting={f"{n}T_{b}": {a: int(((og.Algorithm == a) & (og.Turbines == n) & (og.Budget == b)).sum())
                                                                 for a in GRADM} | {"SLSQP": int(((os_.Turbines == n) & (os_.Budget == b)).sum())}
                                                  for n in (16, 36) for b in (6030, 30030)})


def block_hr_worst(o, out):
    import rev3_sites as S
    H = S.load_hr().copy()
    cols = ["Algorithm", "Seed", "Budget", "Init"]
    oh = open_keys(o, "mpce_hrfix", "SLSQP")
    kh = [tuple(x) for x in oh[cols].itertuples(index=False)]
    res = {}
    for b in S.BUDS:
        y0 = H[(H.Budget == b) & (H.Init == "random")].copy()
        y1, n1 = flip(y0, kh, cols)
        sc = {}
        for name, y in (("stored", y0), ("msslsqp_open_infeasible", y1)):
            r = S.pool_analysis(y, S.HR11, {"best_vs_all": (lambda best: best, None)})
            d = summarize_pool(r, ["SLSQP", "PSOBV", "PSOC", "SSABV"])
            d["ranks_all"] = {a: v["rank"] for a, v in r["methods"].items()}
            d["order_without_msslsqp"] = [a for a in sorted(d["ranks_all"], key=d["ranks_all"].get) if a != "SLSQP"]
            sc[name] = d
        sc["order_of_other_methods_unchanged"] = sc["stored"]["order_without_msslsqp"] == sc["msslsqp_open_infeasible"]["order_without_msslsqp"]
        sc["flipped"] = n1
        res[str(b)] = sc
        log(f"[3] HR1 {b}: flipped {n1}; MS-SLSQP feasible/rank/mean stored "
            f"{sc['stored']['methods']['SLSQP']['feasible']}/{sc['stored']['methods']['SLSQP']['rank']:g}/{sc['stored']['methods']['SLSQP']['mean']:.3f} "
            f"-> {sc['msslsqp_open_infeasible']['methods']['SLSQP']['feasible']}/{sc['msslsqp_open_infeasible']['methods']['SLSQP']['rank']:g}/"
            f"{(sc['msslsqp_open_infeasible']['methods']['SLSQP']['mean'] or 0):.3f}; best {sc['stored']['best']} -> {sc['msslsqp_open_infeasible']['best']}; "
            f"order of the other methods unchanged: {sc['order_of_other_methods_unchanged']}")
    out["hr_worst_case"] = res


def allrun(G, a, b):
    """All-run paired outcome a vs b (rev3_inference B1 main rule: feasible beats infeasible, two feasible runs by the
    stored objective with tie |d| <= 1e-9, two infeasible runs tie)."""
    k = ["Dataset", "Radius", "Turbines", "Seed"]
    A = G[G.Algorithm == a].set_index(k)[["Feasible", "Objective"]]
    B = G[G.Algorithm == b].set_index(k)[["Feasible", "Objective"]]
    j = A.index.intersection(B.index)
    A, B = A.loc[j], B.loc[j]
    fa, fb = A.Feasible.values, B.Feasible.values
    d = A.Objective.values - B.Objective.values
    W = (fa & ~fb) | (fa & fb & (d > 1e-9))
    L = (~fa & fb) | (fa & fb & (d < -1e-9))
    Wn, Ln = int(W.sum()), int(L.sum())
    T = len(j) - Wn - Ln
    from scipy.stats import binomtest
    return dict(n=int(len(j)), W=Wn, T=int(T), L=Ln, score=(Wn + T / 2) / len(j),
                p_sign=float(binomtest(Wn, Wn + Ln).pvalue) if Wn + Ln else 1.0)


def block_bench_worst(o, out):
    import mpce_results as MR
    MR.set_focus("PSOBV")
    ALL, _, _, _ = MR.load(HERE, False)
    out_b = {}
    R6 = ALL[(ALL.Budget == 6030) & (ALL.Init == "random")]
    G = R6[R6.Dataset.isin(["1", "2"]) & R6.Algorithm.isin(MR.MAIN)].reset_index(drop=True)
    cols = ["Algorithm", "Dataset", "Radius", "Turbines", "Seed"]
    os_ = open_keys(o, "mpce_slsqp", "SLSQP").copy()
    os_["Dataset"] = os_.Ds
    os_["Radius"] = os_.Radius.astype(float).astype(int)
    ks = [tuple(x) for x in os_[cols].itertuples(index=False)]
    G1, n1 = flip(G, ks, cols)
    MP = list(MR.MAIN)
    res = {}
    for name, g in (("stored", G), ("msslsqp_open_infeasible", G1)):
        S = MR.case_stats(g, MP)
        R = MR.rank_matrix(S, MP)
        FR = MR.friedman_block(R, MP)
        sl = S[S.Algorithm == "SLSQP"]
        res[name] = dict(avg_rank=FR["avg_rank"], order=sorted(FR["avg_rank"], key=FR["avg_rank"].get),
                         p_holm_vs_focus=FR.get("p_holm_vs_focus") or FR.get("pz_focus"),
                         chi2=FR["chi2"], p=FR["p"], best_cases={a: v for a, v in FR.get("best", {}).items()} if isinstance(FR.get("best"), dict) else None,
                         msslsqp_feasible_pct=float(g[g.Algorithm == "SLSQP"].Feasible.mean() * 100),
                         msslsqp_cases_unqualified=int((~sl.Qualified).sum()),
                         msslsqp_mean_loss_qualified=float(sl[sl.Qualified].Loss.mean()),
                         allrun_psovns_vs_msslsqp=allrun(g, "PSOBV", "SLSQP"))
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
    out_b["main8_68cases"] = res
    log(f"[3] benchmark 68 cases: flipped {n1} MS-SLSQP labels in {len(n_open_case)} cases (max {res['max_open_per_case']} per case); "
        f"MS-SLSQP feasible {st['msslsqp_feasible_pct']:.1f}% -> {wc['msslsqp_feasible_pct']:.1f}%, unqualified cases "
        f"{st['msslsqp_cases_unqualified']} -> {wc['msslsqp_cases_unqualified']}; cases with any rank change {int(ch.sum())}")
    log("    avg ranks stored: " + ", ".join(f"{a} {st['avg_rank'][a]:.2f}" for a in st["order"]))
    log("    avg ranks worst : " + ", ".join(f"{a} {wc['avg_rank'][a]:.2f}" for a in wc["order"]))
    log(f"    all-run PSO-VNS vs MS-SLSQP stored {st['allrun_psovns_vs_msslsqp']} -> worst {wc['allrun_psovns_vs_msslsqp']}")
    # budget / initialization groups with open MS-SLSQP labels (6 largest cases; Horns Rev 1 via hrfix)
    grp = []
    for study in ("mpce_feas", "mpce_b30k", "mpce_b120k", "mpce_hrfix"):
        g = open_keys(o, study, "SLSQP").copy()
        if study != "mpce_hrfix":
            g = g[g.Kind != "hr"]                     # HR rows of feas/b30k/b120k: earlier binning, replaced by mpce_hrfix
        if not len(g):
            continue
        g["Dataset"] = g.Ds
        g["Radius"] = g.Radius.astype(float).astype(int)
        for (ds, rad, n, b, ini), gg in g.groupby(["Dataset", "Radius", "Turbines", "Budget", "Init"]):
            A = ALL[(ALL.Algorithm == "SLSQP") & (ALL.Dataset == ds) & (ALL.Radius == rad) & (ALL.Turbines == n)
                    & (ALL.Budget == b) & (ALL.Init == ini)]
            nf = int(A.Feasible.sum())
            grp.append(dict(study=study, dataset=ds, radius=int(rad), turbines=int(n), budget=int(b), init=ini,
                            runs=int(len(A)), feasible=nf, open=int(len(gg)), feasible_worst=nf - int(len(gg)),
                            qualified_stored=bool(nf >= math.ceil(len(A) / 2)),
                            qualified_worst=bool(nf - len(gg) >= math.ceil(len(A) / 2))))
    out_b["budget_init_groups"] = grp
    nchg = sum(1 for d in grp if d["qualified_stored"] != d["qualified_worst"])
    out_b["budget_init_groups_qualification_lost"] = nchg
    out_b["budget_init_groups_n"] = len(grp)
    log(f"[3] budget/init groups with open MS-SLSQP labels: {len(grp)}; qualification lost under the worst case: {nchg}")
    for d in grp:
        log(f"    {d['study']:10s} {d['dataset']}/{d['radius']}/{d['turbines']} {d['budget']} {d['init']}: feasible {d['feasible']}/{d['runs']}, open {d['open']}")
    out["benchmark_worst_case"] = out_b
    return ALL


# ====================================================================== 4
def parse_xy(c):
    return np.array([[float(v) for v in p.split()] for p in c.split(";")])


_ROWS = {}


def stored_coords(fname, row):
    if fname not in _ROWS:
        _ROWS[fname] = pd.read_csv(os.path.join(HERE, fname), usecols=["Coordinates"]).Coordinates.values
    return _ROWS[fname][int(row)]


def stats_pp(d):
    d = np.abs(np.asarray(d, float))
    if not len(d):
        return dict(n=0)
    return dict(n=int(len(d)), max_pp=float(d.max()), p95_pp=float(np.percentile(d, 95)), median_pp=float(np.median(d)),
                n_gt_0_05pp=int((d > 0.05).sum()), n_gt_0_01pp=int((d > 0.01).sum()), n_exact_zero=int((d == 0).sum()))


def block_transfer(f, ALL, out):
    import mpce_direction as D
    from record_io import decode_coordinates
    t0 = time.time()
    F = f[f.Reproduction.isin(["bit", "noise"])].copy()
    res = {}
    # ---------------- benchmark
    G = ALL[(ALL.Budget == 6030) & (ALL.Init == "random") & ALL.Dataset.isin(["1", "2"]) & ALL.Algorithm.isin(D.ALLM)
            & ALL.Feasible].copy()
    key = ["Algorithm", "Dataset", "Radius", "Turbines", "Seed"]
    G = G.set_index(key)
    Fb = F[F.Study.str.startswith(("fresh_", "mpce_")) & (F.Budget == 6030) & (F.Init.fillna("random") == "random")].copy()
    Fb["Dataset"] = Fb.Dataset.astype(str).str.replace(r"\.0$", "", regex=True)
    Fb = Fb[Fb.Dataset.isin(["1", "2"])]
    rows = []
    for r in Fb.itertuples(index=False):
        k = (r.Algorithm, r.Dataset, int(float(r.Radius)), int(r.Turbines), int(r.Seed))
        if k not in G.index:
            continue
        sc = stored_coords(r.File, r.Row)
        gc = G.loc[k, "Coordinates"]
        if isinstance(gc, pd.Series):
            gc = gc.iloc[0]
        if gc != sc:                                  # this record is not the one the pipeline uses
            continue
        xf = np.array(decode_coordinates(r.CoordinatesFull)); xr = parse_xy(sc)
        jf, i15 = D.bench_objective(xf, r.Dataset, 15); jr, _ = D.bench_objective(xr, r.Dataset, 15)
        gf, _ = D.bench_objective(xf, r.Dataset, 15, wake="gauss"); gr, _ = D.bench_objective(xr, r.Dataset, 15, wake="gauss")
        j1f, i1 = D.bench_objective(xf, r.Dataset, 1); j1r, _ = D.bench_objective(xr, r.Dataset, 1)
        rows.append(dict(File=r.File, Row=r.Row, Algorithm=r.Algorithm, Reproduction=r.Reproduction,
                         J15=100 * (jf - jr) / i15, G15=100 * (gf - gr) / i15, J1=100 * (j1f - j1r) / i1))
    B = pd.DataFrame(rows)
    res["benchmark"] = dict(set="feasible 6,030-evaluation random-start runs of the 11 benchmark methods (mpce_direction.ALLM)",
                            n_layouts=int(len(B)), by_algorithm=B.groupby("Algorithm").size().to_dict() if len(B) else {},
                            reproduction=B.Reproduction.value_counts().to_dict() if len(B) else {},
                            jensen_1deg=stats_pp(B.J15), gauss_1deg=stats_pp(B.G15), paper_15deg=stats_pp(B.J1))
    log(f"[4] benchmark: {len(B)} layouts; 1 deg Jensen {res['benchmark']['jensen_1deg']}; Gauss 1 deg {res['benchmark']['gauss_1deg']}; "
        f"15 deg {res['benchmark']['paper_15deg']} ({time.time() - t0:.0f} s)")
    # ---------------- Horns Rev 1
    import hornsrev_model as hrm
    from py_wake import NOJ
    from py_wake.examples.data.hornsrev1 import Hornsrev1Site, V80, wt16_x, wt16_y
    import py_wake
    b1 = D.hr_bins(1.0, 0.5)
    b5 = D.hr_bins(*D.HR_BINS["5deg_2.5"])
    mdl = NOJ(Hornsrev1Site(), V80(), k=0.04)
    cx, cy = np.mean(wt16_x), np.mean(wt16_y)
    Fh = F[F.Study.isin(["mpce_hrfix", "rev2_gahr"]) & (F.Dataset.astype(str) == "HR") & (F.Turbines == 16)]
    rows = []
    for r in Fh.itertuples(index=False):
        sc = stored_coords(r.File, r.Row)
        xf = np.array(decode_coordinates(r.CoordinatesFull)); xr = parse_xy(sc)
        i1 = D.hr_aep(xf, b1, with_wake=False)
        a1 = D.hr_aep(xf, b1) - D.hr_aep(xr, b1)
        a5 = D.hr_aep(xf, b5) - D.hr_aep(xr, b5)
        i5 = D.hr_aep(xf, b5, with_wake=False)
        pf = float(mdl(xf[:, 0] + cx, xf[:, 1] + cy, wd=b1["wd"], ws=hrm.WS).aep().sum())
        pr = float(mdl(xr[:, 0] + cx, xr[:, 1] + cy, wd=b1["wd"], ws=hrm.WS).aep().sum())
        pi = float(mdl(xf[:, 0] + cx, xf[:, 1] + cy, wd=b1["wd"], ws=hrm.WS).aep(with_wake_loss=False).sum())
        rows.append(dict(File=r.File, Row=r.Row, Algorithm=r.Algorithm, Feasible=r.StoredLabel, Reproduction=r.Reproduction,
                         H1=100 * a1 / i1, H5=100 * a5 / i5, H1gwh=a1, PW=100 * (pf - pr) / pi, PWgwh=pf - pr))
    Hh = pd.DataFrame(rows)
    res["horns_rev"] = dict(set="mpce_hrfix + rev2_gahr, 16 turbines, all stored layouts (as rev3_sites.reeval_layouts)",
                            n_layouts=int(len(Hh)), n_feasible=int(Hh.Feasible.sum()) if len(Hh) else 0,
                            reproduction=Hh.Reproduction.value_counts().to_dict() if len(Hh) else {},
                            jensen_1deg=stats_pp(Hh.H1), paper_5deg=stats_pp(Hh.H5), pywake_noj_1deg=stats_pp(Hh.PW),
                            max_abs_gwh_1deg=float(Hh.H1gwh.abs().max()) if len(Hh) else None,
                            max_abs_gwh_pywake=float(Hh.PWgwh.abs().max()) if len(Hh) else None,
                            pywake_version=py_wake.__version__)
    log(f"[4] Horns Rev 1: {len(Hh)} layouts; 1 deg {res['horns_rev']['jensen_1deg']}; PyWake {res['horns_rev']['pywake_noj_1deg']} "
        f"({time.time() - t0:.0f} s)")
    # ---------------- Lillgrund
    import rev2_site_model as lg
    from py_wake.examples.data import lillgrund as L
    wd1, f1, p1 = lg.bins(1.0, 0.5)
    lm = NOJ(L.LillgrundSite(), L.SWT23(), k=0.04)
    gx, gy = lg.WT_X[lg.I16].mean(), lg.WT_Y[lg.I16].mean()
    Fl = F[F.Study.isin(["rev2_lg16", "rev2_lg16b"])]
    rows = []
    for r in Fl.itertuples(index=False):
        sc = stored_coords(r.File, r.Row)
        xf = np.array(decode_coordinates(r.CoordinatesFull)); xr = parse_xy(sc)
        i1 = lg.aep_gwh(xf, False, wd1, f1, p1)
        a1 = lg.aep_gwh(xf, True, wd1, f1, p1) - lg.aep_gwh(xr, True, wd1, f1, p1)
        a5 = lg.aep_gwh(xf) - lg.aep_gwh(xr)
        pf = float(lm(xf[:, 0] + gx, xf[:, 1] + gy, wd=wd1, ws=lg.WS).aep().sum())
        pr = float(lm(xr[:, 0] + gx, xr[:, 1] + gy, wd=wd1, ws=lg.WS).aep().sum())
        pi = float(lm(xf[:, 0] + gx, xf[:, 1] + gy, wd=wd1, ws=lg.WS).aep(with_wake_loss=False).sum())
        rows.append(dict(File=r.File, Row=r.Row, Algorithm=r.Algorithm, Feasible=r.StoredLabel, Reproduction=r.Reproduction,
                         L1=100 * a1 / i1, L5=100 * a5 / lg.aep_gwh(xf, False), PW=100 * (pf - pr) / pi))
    Ll = pd.DataFrame(rows)
    res["lillgrund"] = dict(set="rev2_lg16 + rev2_lg16b (as rev3_sites.reeval_layouts)", n_layouts=int(len(Ll)),
                            reproduction=Ll.Reproduction.value_counts().to_dict() if len(Ll) else {},
                            jensen_1deg=stats_pp(Ll.L1), paper_5deg=stats_pp(Ll.L5), pywake_noj_1deg=stats_pp(Ll.PW))
    log(f"[4] Lillgrund: {len(Ll)} layouts; 1 deg {res['lillgrund']['jensen_1deg']}; PyWake {res['lillgrund']['pywake_noj_1deg']} "
        f"({time.time() - t0:.0f} s)")
    res["seconds"] = time.time() - t0
    out["transfer"] = res


# ====================================================================== LaTeX
def num(v, d=3):
    return "--" if v is None or (isinstance(v, float) and not np.isfinite(v)) else f"{v:.{d}f}"


def sci(v):
    if v is None or not np.isfinite(v):
        return "--"
    if v == 0:
        return "0"
    m, e = f"{v:.1e}".split("e")
    return f"${m}\\times10^{{{int(e)}}}$"


def mm(v):
    return "--" if v is None else f"{1000 * v:.2f}"


def mmf(v):
    """mm with two decimals, or scientific notation below 0.01 mm."""
    if v is None:
        return "--"
    return f"{1000 * v:.2f}" if (v == 0 or abs(1000 * v) >= 0.01) else sci(1000 * v)


def ppf(v):
    return f"{v:.4f}" if (v == 0 or abs(v) >= 1e-3) else sci(v)


def write_tex(out, path):
    T = []
    # ---- table 1: residual bounds + tolerance sensitivity
    lines = []
    order = ["rev2_grad", "mpce_iea", "mpce_slsqp", "mpce_feas", "mpce_b30k", "mpce_b120k", "mpce_hrfix", "fresh_hr"]
    rows = sorted(out["open_by_study_method"], key=lambda d: (order.index(d["study"]) if d["study"] in order else 99,
                                                               d["algorithm"], d["kind"]))
    for d in rows:
        rel_r = sci(d["bound_true_boundary_rel_radius"]) if d["bound_true_boundary_rel_radius"] is not None else "--"
        sl = STUDY_LAB.get(d['study'], d['study']) + (f", {KIND_LAB[d['kind']]}" if d['study'] in MULTI_KIND else "")
        lines.append(f"{sl} & {ALG_LAB.get(d['algorithm'], d['algorithm'])} & "
                     f"{d['n']:,} & {','.join(map(str, d['decimals']))} & {mm(d['max_boundary_excess_rounded_m'])} & "
                     f"{mm(d['bound_true_boundary_excess_m'])} & {rel_r} & {mm(d['max_spacing_deficit_rounded_m'])} & "
                     f"{mm(d['bound_true_spacing_deficit_m'])} & {sci(d['bound_true_spacing_rel_smin'])} & "
                     f"{mmf(d['true_spacing_deficit_from_minspacing_m'])} \\\\")
    tot = out["open_total"]
    lines.append("\\midrule")
    lines.append(f"Total & & {tot['n']:,} & & & & & & & & \\\\")
    S = out["tolerance_sensitivity"]
    lines2 = []
    for name, lab in (("open920", "All open labels"), ("open_exact_gradient", "Open, exact gradient"),
                      ("exact_gradient240", "All exact-gradient runs$^b$")):
        g = S[name]
        c = []
        for t in TAUS:
            v = g["tol"][f"{t:g}"]
            vr = v["R_with_reruns"] if name == "exact_gradient240" else v["R"]
            vm = v["M_with_reruns"] if name == "exact_gradient240" else v["M"]
            c.append(f"{vr['feasible']}/{vr['infeasible']}/{vr['undecidable']}")
            c.append(f"{vm['feasible']}/{vm['infeasible']}/{vm['undecidable']}")
        lines2.append(f"{lab} ({g['n']}) & " + " & ".join(c) + f" & {mm(g['tau_star_R'])} & {mm(g['tau_star_M'])} \\\\")
    cap = ("Strict feasibility labels ($10^{-6}$~m) that neither the stored rounded coordinates nor the full-precision "
           "reruns decide (all SLSQP-based runs, all labelled feasible): residuals implied by the stored coordinates and "
           "guaranteed bounds on the residuals of the unrounded layouts (top), and the labels decided at larger tolerances (bottom).")
    note = ("Top: Exc.: largest boundary excess of the rounded layout (signed distance outside the site, negative = inside); "
            "Spc.: largest spacing deficit of the rounded layout ($S_{\\min}$ minus the smallest pair distance; negative = "
            "spacing larger than $S_{\\min}$); bound: guaranteed upper bound on the same quantity for the unrounded layout, "
            "rounded value plus $\\sqrt2\\,\\varepsilon$ (boundary) or $2\\sqrt2\\,\\varepsilon$ (spacing), "
            "$\\varepsilon=0.5\\cdot10^{-d}$~m for $d$ stored decimals; rel.: bound divided by the farm radius $r$ (boundary) "
            "or by $S_{\\min}$ (spacing); MinSp.: true spacing deficit from the stored minimum-spacing value, which the "
            "drivers wrote at full precision from the unrounded layout. All lengths in mm. Geometry: benchmark circle "
            "$r=500$, 750, 1000~m, $S_{\\min}=4D=308$~m; IEA37 Case Study~1 circle $r=1300$~m (16 turbines) and 2000~m (36), "
            "$S_{\\min}=2D=260$~m; Horns Rev~1 (HR1) parallelogram of the installed block, $S_{\\min}=4D=320$~m (relative "
            "boundary bound not defined). HR1 rows of the 30,030 and 120,030-evaluation and feasible-initialization "
            "studies and the last row use the earlier direction binning and do not enter the Horns Rev~1 results. "
            "Bottom: labels decided feasible / decided infeasible / undecidable at tolerance $\\tau$: decided feasible if both "
            "slacks minus their rounding bounds are at least $-\\tau$, decided infeasible if a slack plus its bound is below "
            "$-\\tau$. R: both constraints from the rounded coordinates; M: spacing from the stored full-precision minimum "
            "spacing, boundary from the rounded coordinates. $\\tau^\\ast$ (mm): smallest tolerance that decides every label "
            "of the row feasible. $^b$All 240 exact-gradient runs, the 20 labels decided by bit-identical or "
            "floating-point-noise reruns counted as decided at every $\\tau$.")
    head = ("Study & Method & $n$ & $d$ & Exc. & bound & rel. & Spc. & bound & rel. & MinSp.")
    tab1 = (f"\\begin{{table*}}[!htb]\n\\centering\n\\caption{{{cap}}}\n\\label{{tab:S-sa-feasbounds}}\n"
            "\\scriptsize\\setlength{\\tabcolsep}{2.5pt}\n"
            "\\begin{tabular}{llrrrrrrrrr}\n\\toprule\n"
            "& & & & \\multicolumn{3}{c}{Boundary (mm)} & \\multicolumn{3}{c}{Spacing, rounded coord. (mm)} & \\\\\n"
            "\\cmidrule(lr){5-7}\\cmidrule(lr){8-10}\n" + head + " \\\\\n\\midrule\n" + "\n".join(lines) +
            "\n\\bottomrule\n\\end{tabular}\n\n\\vspace{4pt}\n"
            "\\begin{tabular}{lrrrrrrrrrr}\n\\toprule\n"
            "& \\multicolumn{2}{c}{$\\tau=10^{-6}$~m} & \\multicolumn{2}{c}{$10^{-4}$~m} & \\multicolumn{2}{c}{$10^{-3}$~m} & "
            "\\multicolumn{2}{c}{$10^{-2}$~m} & \\multicolumn{2}{c}{$\\tau^\\ast$ (mm)} \\\\\n"
            "\\cmidrule(lr){2-3}\\cmidrule(lr){4-5}\\cmidrule(lr){6-7}\\cmidrule(lr){8-9}\\cmidrule(lr){10-11}\n"
            "Labels ($n$) & R & M & R & M & R & M & R & M & R & M \\\\\n\\midrule\n" + "\n".join(lines2) +
            "\n\\bottomrule\n\\end{tabular}\n"
            f"\\par\\vspace{{2pt}}\\parbox{{\\textwidth}}{{\\scriptsize {note}}}\n\\end{{table*}}\n")
    T.append(tab1)
    # ---- table 2: worst case
    W = out["iea_worst_case"]["settings"]
    lines = []
    scen = [("stored", "Stored labels"), ("A_adversarial_qualified", "Adversarial, qualified$^a$"),
            ("W1_exactgrad_open_infeasible", "Open exact-grad.\\ infeasible"),
            ("W2_plus_msslsqp_open_infeasible", "+ open MS-SLSQP infeasible")]
    for n in (16, 36):
        for b in (6030, 30030):
            k = f"{n}T_{b}"
            for i, (s, lab) in enumerate(scen):
                d = W[k][s]
                m = d["methods"]
                best = LAB_SHORT.get(d["best"], d["best"])
                ns = [LAB_SHORT.get(x, x) for x in d["best_family"]["not_significant"]]
                c = [f"{m[a]['feasible']} / {m[a]['rank']:g}" for a in ("PSOSLSQPX", "SLSQPX", "SLSQP")]
                pre = (f"\\multirow{{4}}{{*}}{{{n}, {b:,}}}".replace(",", "{,}").replace("{,} ", ", ", 1) if i == 0 else "")
                lines.append(f"{pre} & {lab} & " + " & ".join(c) + f" & {best} & {big0(d['best_mean'])} & "
                             f"{', '.join(ns) if ns else 'none'} \\\\")
            if not (n == 36 and b == 30030):
                lines.append("\\midrule")
    BW = out["benchmark_worst_case"]["main8_68cases"]
    st, wc = BW["stored"], BW["msslsqp_open_infeasible"]
    HW = out["hr_worst_case"]
    hr_txt = "; ".join(
        f"{int(b):,}".replace(",", "{,}") + f" evaluations: {HW[b]['flipped']} labels, MS-SLSQP {HW[b]['stored']['methods']['SLSQP']['feasible']}"
        f"$\\to${HW[b]['msslsqp_open_infeasible']['methods']['SLSQP']['feasible']} feasible runs, rank "
        f"{HW[b]['stored']['methods']['SLSQP']['rank']:g}$\\to${HW[b]['msslsqp_open_infeasible']['methods']['SLSQP']['rank']:g}, best method "
        f"{LAB_SHORT.get(HW[b]['msslsqp_open_infeasible']['best'])}" for b in HW)
    gq = out["benchmark_worst_case"]
    note = ("IEA37 Case Study~1, 13-method pool (scenario: turbines, evaluations; random starts, 30 seed-paired runs). "
            "Feas./Rank: feasible runs and feasibility-aware rank of the mean AEP (fewer than 15 feasible runs = ranked below "
            "the qualified methods). Best: rank-1 method and its mean AEP (MWh) over its feasible runs; n.s.: methods not "
            "significantly different from the best method (seed-paired Wilcoxon signed-rank, infeasible runs below every "
            "feasible run, Holm over its 12 comparisons, $\\alpha=0.05$). Scenarios: stored labels; $^a$in each scenario the "
            "open runs with the highest AEP of each exact-gradient method are taken as infeasible until only 15 of its 30 runs "
            "remain feasible (the lowest mean it can have while ranked by its mean); every open exact-gradient label taken "
            "as infeasible; in addition every open label of MS-SLSQP with finite differences taken as infeasible. "
            f"Benchmark (68 cases, eight methods): with all {BW['flipped']} open MS-SLSQP labels infeasible "
            f"({BW['cases_with_open_msslsqp']} cases, at most {BW['max_open_per_case']} per case), MS-SLSQP has "
            f"{num(wc['msslsqp_feasible_pct'], 1)}\\% feasible runs and fewer than 15 feasible runs in "
            f"{wc['msslsqp_cases_unqualified']} case{'s' if wc['msslsqp_cases_unqualified'] != 1 else ''}; average ranks "
            + ", ".join(f"{LAB_SHORT.get(a, a)} {st['avg_rank'][a]:.2f}$\\to${wc['avg_rank'][a]:.2f}" for a in st["order"])
            + f" (order {'unchanged' if BW['order_unchanged'] else 'changed'}); all-run outcome PSO-VNS vs.\\ MS-SLSQP "
            f"{st['allrun_psovns_vs_msslsqp']['W']}/{st['allrun_psovns_vs_msslsqp']['T']}/{st['allrun_psovns_vs_msslsqp']['L']}"
            f"$\\to${wc['allrun_psovns_vs_msslsqp']['W']}/{wc['allrun_psovns_vs_msslsqp']['T']}/{wc['allrun_psovns_vs_msslsqp']['L']} "
            f"(wins/ties/losses). Horns Rev~1 pool (11 methods), open MS-SLSQP labels infeasible: {hr_txt}. "
            f"Budget and initialization studies: {gq['budget_init_groups_n']} case groups contain open MS-SLSQP labels; "
            f"MS-SLSQP keeps at least half of its runs feasible in "
            f"{gq['budget_init_groups_n'] - gq['budget_init_groups_qualification_lost']} of them and loses that qualification in "
            f"{gq['budget_init_groups_qualification_lost']} (Horns Rev~1).")
    cap = ("Worst-case effect of the undecided strict labels on the IEA37 comparison: feasible runs and feasibility-aware "
           "ranks of the SLSQP-based methods, best method and the methods not significantly worse than it, when undecided "
           "labels are taken as infeasible.")
    tab2 = (f"\\begin{{table*}}[!htb]\n\\centering\n\\caption{{{cap}}}\n\\label{{tab:S-sa-feasworst}}\n"
            "\\scriptsize\\setlength{\\tabcolsep}{2.5pt}\n\\begin{tabular}{llccclrp{5.2cm}}\n\\toprule\n"
            "& & PSO-SLSQP, exact & MS-SLSQP, exact & MS-SLSQP & & & \\\\\n"
            "Scenario & Labels & Feas.\\ / Rank & Feas.\\ / Rank & Feas.\\ / Rank & Best & Mean AEP & n.s.\\ vs.\\ best \\\\\n"
            "\\midrule\n" + "\n".join(lines) + "\n\\bottomrule\n\\end{tabular}\n"
            f"\\par\\vspace{{2pt}}\\parbox{{\\textwidth}}{{\\scriptsize {note}}}\n\\end{{table*}}\n")
    T.append(tab2)
    # ---- table 3: transfer
    if "transfer" in out:
        X = out["transfer"]
        lines = []

        def row(site, ev, s):
            if not s.get("n"):
                return f"{site} & {ev} & 0 & -- & -- & -- \\\\"
            return (f"{site} & {ev} & {big0(s['n'])} & {ppf(s['p95_pp'])} & {ppf(s['max_pp'])} & {s['n_gt_0_05pp']} \\\\")
        bm, hrr, lgg = X["benchmark"], X["horns_rev"], X["lillgrund"]
        lines += [row("Benchmark", "15$^\\circ$ (optimization model)", bm["paper_15deg"]),
                  row("", "1$^\\circ$, Jensen", bm["jensen_1deg"]), row("", "1$^\\circ$, Gaussian", bm["gauss_1deg"]),
                  row("Horns Rev~1", "5$^\\circ$ (optimization model)", hrr["paper_5deg"]),
                  row("", "1$^\\circ$, Jensen", hrr["jensen_1deg"]), row("", "PyWake NOJ, 1$^\\circ$", hrr["pywake_noj_1deg"]),
                  row("Lillgrund", "5$^\\circ$ (optimization model)", lgg["paper_5deg"]),
                  row("", "1$^\\circ$, Jensen", lgg["jensen_1deg"]), row("", "PyWake NOJ, 1$^\\circ$", lgg["pywake_noj_1deg"])]
        note = ("Layouts whose optimization run was rerun with a 17-digit writer and reproduced the stored record bit for bit "
                "or to floating-point noise (so that the unrounded final layout is known) and that belong to a re-evaluation "
                "set: benchmark, feasible runs of the 11 benchmark methods (68 cases, 6{,}030 evaluations, random starts); "
                "Horns Rev~1, the 16-turbine runs of the ten methods and GA; Lillgrund, the 16-turbine runs at both budgets. "
                "$|\\Delta|$: absolute difference between the re-evaluated objective of the unrounded and of the rounded "
                "(stored) layout, in percentage points of the evaluator's wake-free value; P95: 95th percentile; "
                "$>0.05$: layouts with $|\\Delta|>0.05$~pp. All eligible layouts evaluated (no subsampling). The "
                "recovered layouts are those of runs whose strict label the rounded coordinates could not decide (a constraint "
                "active at the optimum) and a small determinism sample; they are not a random sample of all layouts.")
        cap = ("Sensitivity of the objective and of the re-evaluations to the rounding of the stored coordinates, from "
               "layouts recovered at full precision.")
        tab3 = (f"\\begin{{table}}[!htb]\n\\centering\n\\caption{{{cap}}}\n\\label{{tab:S-sa-roundtransfer}}\n"
                "\\scriptsize\\setlength{\\tabcolsep}{3pt}\n\\begin{tabular}{llrrrr}\n\\toprule\n"
                "Site & Evaluator & $n$ & P95 $|\\Delta|$ & Max $|\\Delta|$ & $>0.05$ \\\\\n\\midrule\n" + "\n".join(lines) +
                "\n\\bottomrule\n\\end{tabular}\n"
                f"\\par\\vspace{{2pt}}\\parbox{{\\columnwidth}}{{\\scriptsize {note}}}\n\\end{{table}}\n")
        T.append(tab3)
    open(path, "w").write("%% generated by rev3_feasbounds.py -- do not edit by hand\n" + "\n".join(T))


LAB_SHORT = {"PSOBV": "PSO-VNS", "PSOC": "PSO", "SSABV": "SSA-VNS", "SSA": "SSA", "LXSSA": "LX-SSA", "DE": "DE",
             "BVNS": "VNS", "SLSQP": "MS-SLSQP", "LXBV": "LX-SSA-VNS", "RSVNS": "RS-VNS", "RSDVNS": "RSD-VNS", "GA": "GA",
             "SLSQPX": "MS-SLSQP (exact)", "PSOSLSQPX": "PSO-SLSQP (exact)"}


def big0(v):
    return "--" if v is None else f"{v:,.0f}".replace(",", "{,}")


def clean(o):
    if isinstance(o, dict):
        return {str(k): clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [clean(v) for v in o]
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return None if not np.isfinite(o) else float(o)
    if isinstance(o, (np.bool_,)):
        return bool(o)
    if isinstance(o, float) and not np.isfinite(o):
        return None
    return o


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--skip-transfer", action="store_true")
    ap.add_argument("--out-dir", default=HERE)
    a = ap.parse_args(argv)
    t0 = time.time()
    r, f = load_audit()
    out = dict(script="rev3_feasbounds.py", generated=time.strftime("%Y-%m-%d %H:%M:%S"), tolerance_m=TOL,
               float_pad_m=FLOAT_PAD, taus_m=list(TAUS),
               bounds="boundary: e + sqrt(2) eps; spacing: s + 2 sqrt(2) eps; eps = 0.5 10^-d (d stored decimals)",
               geometry=GEOM, audited_records=int(len(r)))
    o = block_bounds(r, out)
    block_iea_worst(o, out)
    block_hr_worst(o, out)
    ALL = block_bench_worst(o, out)
    if not a.skip_transfer:
        block_transfer(f, ALL, out)
    elif os.path.exists(os.path.join(a.out_dir, "rev3_feasbounds.json")):
        prev = json.load(open(os.path.join(a.out_dir, "rev3_feasbounds.json")))
        if "transfer" in prev:
            out["transfer"] = prev["transfer"]
    out["seconds"] = time.time() - t0
    json.dump(clean(out), open(os.path.join(a.out_dir, "rev3_feasbounds.json"), "w"), indent=1)
    write_tex(clean(out), os.path.join(a.out_dir, "rev3_feasbounds_tables.tex"))
    open(os.path.join(a.out_dir, "rev3_feasbounds.log"), "w").write("\n".join(LOG) + "\n")
    log(f"done in {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
