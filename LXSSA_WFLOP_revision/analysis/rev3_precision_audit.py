"""Precision audit: audit of the stored (rounded) coordinates of ALL run records.

Usage:  python3 rev3_precision_audit.py [--out-dir DIR]
Output: rev3_precision_audit.json            (summary per study and per file)
        rev3_precision_audit_records.csv     (one row per record: keys, decimals, slacks, class, replay difference)

For every record with coordinates (original: fresh_*.csv, mpce_*_s*of*.csv; additional experiments: rev2_*_s*of*.csv) the script
  1. reads the number of decimals d actually stored in that record's Coordinates (max over its tokens) and sets the
     per-coordinate rounding bound eps = 0.5 * 10^-d (+ 1e-9 m for binary64 parsing / arithmetic);
  2. recomputes, from the rounded coordinates, the two constraint slacks of the strict label used by every driver
       spacing slack  s_sp = MinSpacing - (s_min - 1e-6)          (label needs s_sp >= 0)
       boundary slack s_bd = 1e-6 - max outside distance          (label needs s_bd >= 0; signed, < 0 inside)
     with the correct site geometry (circle of the case radius for the benchmark and IEA37; the Horns Rev 1
     parallelogram of hornsrev_model.site(N) and the Lillgrund parallelogram of rev2_site_model.site(16)) and the
     minimum spacing of the study (4D = 8R benchmark; SMin column of the spacing study (5D / 6D); 2D = 260 m IEA37;
     4D = 320 m Horns Rev 1; 3D = 279 m Lillgrund);
  3. bounds the effect of rounding: a turbine moves by at most sqrt(2) eps, so a pair distance changes by at most
     b_sp = 2 sqrt(2) eps and a distance to the circle or to an edge line by at most b_bd = sqrt(2) eps;
  4. classifies the stored label
       (i)   confirmed      : label feasible and s_sp > b_sp and s_bd > b_bd, or label infeasible and
                              (s_sp < -b_sp or s_bd < -b_bd)
       (ii)  contradicted   : the opposite certificate (label feasible but s_sp < -b_sp or s_bd < -b_bd, or label
                              infeasible but both slacks exceed their bounds)
       (iii) undecidable    : otherwise (the true slack sign cannot be decided within the rounding bound);
     a second classification ("hybrid") uses the stored MinSpacing column, which pandas wrote with round-trip
     (17-digit) precision from the unrounded layout, for the spacing constraint (exact) and the rounded coordinates
     only for the boundary; it is reported separately and is NOT used to select reruns;
  5. replays the objective from the rounded coordinates with the evaluator of the record (benchmark: wflop_model;
     IEA37: iea37_model; Lillgrund: rev2_site_model; Horns Rev 1: hornsrev_model for mpce_hrfix and rev2_gahr, and the
     alternative direction binning (bins centred at 0, 5, ..., 355 deg, nearest sector by np.round, as those
     runs used) for all other Horns Rev records) and reports replay - stored (objective units and pp of
     the ideal value).
Nothing is relabelled; existing files are only read.
"""
import os, re, sys, glob, json, math, argparse
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import wflop_model as wm
import hornsrev_model as hr
import iea37_model as iea
import rev2_site_model as lg

TOL = 1e-6
FLOAT_PAD = 1e-9

# ------------------------------------------------------------------ alternative Horns Rev 1 direction binning (replay only)
_WD_OLD = np.arange(0.0, 360.0, 5.0)
_SEC_OLD = (np.round(_WD_OLD / 30.0).astype(int)) % 12
_F_OLD = hr.SEC_F[_SEC_OLD] / 6.0
_A_OLD, _K_OLD = hr.SEC_A[_SEC_OLD], hr.SEC_K[_SEC_OLD]
_PWS_OLD = (np.exp(-(hr._lo[None] / _A_OLD[:, None]) ** _K_OLD[:, None])
            - np.exp(-(hr._hi[None] / _A_OLD[:, None]) ** _K_OLD[:, None]))
_TOW_OLD = np.deg2rad(270.0 - _WD_OLD)


def hr_aep(xy, model="new", with_wake=True):
    """hornsrev_model.aep_gwh with the current ("new") or the alternative ("old") direction bins."""
    if model == "new":
        return hr.aep_gwh(xy, with_wake)
    xy = np.asarray(xy, float)
    n = len(xy)
    if with_wake:
        dx = xy[:, 0][:, None] - xy[:, 0][None, :]
        dy = xy[:, 1][:, None] - xy[:, 1][None, :]
        c = np.cos(_TOW_OLD)[:, None, None]; s = np.sin(_TOW_OLD)[:, None, None]
        x = dx[None] * c + dy[None] * s
        lat = np.abs(-dx[None] * s + dy[None] * c)
        inw = (x > 0) & (lat < hr.RR + hr.KW * x)
        g4 = np.where(inw, (hr.RR / (hr.RR + hr.KW * np.maximum(x, 0))) ** 4, 0.0)
        G = np.sqrt(g4.sum(2))
        u = hr.WS[None, :, None] * (1 - hr.A_CT[None, :, None] * G[:, None, :])
    else:
        u = np.broadcast_to(hr.WS[None, :, None], (len(_WD_OLD), len(hr.WS), n))
    p = np.interp(u, hr.WS_TAB, hr.P_TAB)
    return np.einsum("dsn,ds,d->", p, _PWS_OLD, _F_OLD) * hr.HOURS / 1e6


# ------------------------------------------------------------------ studies (file families)
# (study id, glob, description, role)   role: "used" = feeds a published table, "superseded", "duplicate"
STUDIES = [
    ("fresh_grid", "fresh_grid.csv", "Benchmark 68 cases, 6 methods (LX-SSA, SSA, PSO, DE, mod. VNS, MS-SLSQP)", "used"),
    ("fresh_bgrid", "fresh_bgrid.csv", "Benchmark 68 cases, LX-SSA-VNS and SSA-VNS", "used"),
    ("fresh_vgrid", "fresh_vgrid.csv", "Benchmark 68 cases, VNS", "used"),
    ("fresh_bsplit", "fresh_bsplit.csv", "Split 25/75, LX-SSA-VNS (12 cases)", "used"),
    ("fresh_hgrid", "fresh_hgrid.csv", "Benchmark 68 cases, LX-VNS (not used in the results)", "superseded"),
    ("fresh_hsplit", "fresh_hsplit.csv", "Split, LX-VNS (not used in the results)", "superseded"),
    ("fresh_hr", "fresh_*hr*.csv", "Horns Rev 1, 16/80 turbines, preliminary studies (alternative binning)", "superseded"),
    ("mpce_psobv", "mpce_psobv_s*of2.csv", "PSO-VNS, 68 cases (+ HR16 alternative binning)", "used"),
    ("mpce_psoc", "mpce_psoc_s*of1.csv", "PSO (constriction), 68 cases", "used"),
    ("mpce_slsqp", "mpce_slsqp_s*of1.csv", "MS-SLSQP, 68 cases", "used"),
    ("mpce_rsvns", "mpce_rsvns_s*of1.csv", "RS-VNS, 68 cases", "used"),
    ("mpce_rsdisc", "mpce_rsdisc_s*of1.csv", "RSD-VNS, 68 cases", "used"),
    ("mpce_psosplit", "mpce_psosplit_s*of1.csv", "PSO-VNS split 25/75 (12 cases)", "used"),
    ("mpce_omega90", "mpce_omega90_s*of1.csv", "PSO-VNS split 0.9 (12 cases)", "used"),
    ("mpce_csweep", "mpce_csweep_s*of1.csv", "PSO coefficient sweep (12 cases)", "used"),
    ("mpce_feas", "mpce_feasx_s*of8.csv", "Feasible initialization, 8 methods (6 cases + HR16 alternative binning)", "used"),
    ("mpce_feasp", "mpce_feasp_s*of1.csv", "Feasible initialization, PSO-VNS", "used"),
    ("mpce_b30k", "mpce_b30k_s*of3.csv", "30,030 evaluations, 9 methods (6 cases + HR16 alternative binning)", "used"),
    ("mpce_b30kp", "mpce_b30kp_s*of1.csv", "30,030 evaluations, PSO-VNS", "used"),
    ("mpce_b120k", "mpce_b120k_s*of8.csv", "120,030 evaluations, 9 methods (6 cases + HR16 alternative binning)", "used"),
    ("mpce_b120kp", "mpce_b120kp_s*of4.csv", "120,030 evaluations, PSO-VNS", "used"),
    ("mpce_hr16new", "mpce_hr16new_s*of1.csv", "HR16 PSO and RS-VNS (alternative binning)", "superseded"),
    ("mpce_hrfix", "mpce_hrfix_s*of24.csv", "Horns Rev 1 16 turbines, 10 methods, all budgets", "used"),
    ("mpce_iea", "mpce_iea[13]6*_s*of1.csv", "IEA37 CS1 16/36 turbines, 10 methods", "used"),
    ("rev2_ga", "rev2_ga_s*of10.csv", "GA, 68 cases", "used"),
    ("rev2_gahr", "rev2_gahr_s*of2.csv", "GA, Horns Rev 1", "used"),
    ("rev2_gaiea", "rev2_gaiea_s*of2.csv", "GA, IEA37", "used"),
    ("rev2_grad", "rev2_grad_s*of40.csv", "Exact-gradient SLSQP / PSO-SLSQP, IEA37", "used"),
    ("rev2_laplace", "rev2_laplace_s*of12.csv", "Laplace ablation (12 cases)", "used"),
    ("rev2_spacing", "rev2_spacing_s*of40.csv", "Spacing 5D/6D (12 cases)", "used"),
    ("rev2_lg16", "rev2_lg16_s*of4.csv", "Lillgrund, 6,030 evaluations", "used"),
    ("rev2_lg16b", "rev2_lg16b_s*of30.csv", "Lillgrund, 30,030 evaluations", "used"),
    ("mpce_feas_merged", "mpce_feas_s0of1.csv", "merge of the mpce_feasx shards (duplicate, not counted)", "duplicate"),
]
HR_NEW_FILES = ("mpce_hrfix", "rev2_gahr")
KEYS = ["File", "Row", "Algorithm", "Dataset", "Radius", "Turbines", "Seed", "Budget", "Init", "Spacing"]


def decimals(s):
    return max((len(t.split(".")[1]) if "." in t else 0) for t in re.split(r"[ ;]", s.strip()) if t)


def parse(s):
    return np.array([[float(v) for v in p.split()] for p in s.split(";")])


def pair_min(xy):
    d = np.sqrt(((xy[:, None, :] - xy[None, :, :]) ** 2).sum(-1))
    return d[np.triu_indices(len(xy), 1)].min() if len(xy) > 1 else np.inf


_POLY = {}


def geometry(kind, n):
    if kind == "hr":
        if n not in _POLY:
            _POLY[n] = hr.site(n)[1]
        return _POLY[n]
    if kind == "lg":
        if "lg" not in _POLY:
            _POLY["lg"] = lg.site(16)[1]
        return _POLY["lg"]
    return None


def record_kind(row, fam):
    ds = str(row.get("Dataset", "HR" if fam.startswith("fresh_") and "hr" in fam else ""))
    if ds == "HR":
        return "hr"
    if ds == "LG":
        return "lg"
    if ds.startswith("IEA37"):
        return "iea"
    return "grid"


def smin_of(row, kind):
    if kind == "hr":
        return 4 * hr.D
    if kind == "lg":
        return float(row["Spacing"]) if "Spacing" in row and pd.notna(row["Spacing"]) else lg.SMIN
    if kind == "iea":
        return iea.SMIN
    if "SMin" in row and pd.notna(row.get("SMin")):
        return float(row["SMin"])
    return 8 * wm.R


def evaluate(xy, row, kind, hrmodel):
    """(stored objective, replayed objective, ideal) in the record's objective units."""
    if kind == "grid":
        obj, ideal = wm.farm_objective(xy, int(row["Dataset"]))
        return float(row["Objective"]), obj, float(row["Ideal"])
    if kind == "iea":
        return float(row["Objective"]), iea.aep(xy), float(row["Ideal"])
    if kind == "lg":
        return float(row["Objective"]), lg.aep_gwh(xy), float(row["Ideal"])
    stored = float(row["AEP"]) if "AEP" in row else float(row["Objective"])
    ideal = float(row["IdealAEP"]) if "IdealAEP" in row else float(row["Ideal"])
    return stored, hr_aep(xy, hrmodel), ideal


def slacks(xy, kind, row, smin):
    """(spacing slack, boundary slack) of the strict 1e-6 label, signed (>= 0 means satisfied)."""
    n = len(xy)
    s_sp = (pair_min(xy) - (smin - TOL)) if n > 1 else np.inf
    if kind in ("grid", "iea"):
        rad = float(row["Radius"])
        viol = np.sqrt((xy ** 2).sum(1)).max() - rad
    else:
        poly = geometry(kind, n)
        a, b = poly, np.roll(poly, -1, axis=0)
        e = b - a
        nrm = np.c_[e[:, 1], -e[:, 0]] / np.linalg.norm(e, axis=1)[:, None]
        viol = np.einsum("ijk,jk->ij", xy[:, None, :] - a[None], nrm).max()   # signed (max over turbines of max edge)
    return s_sp, TOL - viol


def classify(label, s_sp, s_bd, b_sp, b_bd):
    cert_feas = (s_sp > b_sp) and (s_bd > b_bd)
    cert_inf = (s_sp < -b_sp) or (s_bd < -b_bd)
    if label:
        return "i" if cert_feas else ("ii" if cert_inf else "iii")
    return "i" if cert_inf else ("ii" if cert_feas else "iii")


def audit_file(path, study, role):
    fam = os.path.basename(path)
    df = pd.read_csv(path)
    out = []
    for r, row in enumerate(df.to_dict("records")):
        kind = record_kind(row, fam)
        hrmodel = "new" if (kind == "hr" and fam.startswith(HR_NEW_FILES)) else "old"
        s = row["Coordinates"]
        d = decimals(s)
        eps = 0.5 * 10.0 ** (-d)
        b_sp, b_bd = 2 * math.sqrt(2) * eps + FLOAT_PAD, math.sqrt(2) * eps + FLOAT_PAD
        xy = parse(s)
        n = int(row["Turbines"])
        assert len(xy) == n, (path, r)
        smin = smin_of(row, kind)
        s_sp, s_bd = slacks(xy, kind, row, smin)
        label = bool(row["Feasible"])
        cls = classify(label, s_sp, s_bd, b_sp, b_bd)
        ms = float(row["MinSpacing"])
        s_sp_col = ms - (smin - TOL)
        # hybrid: spacing exact from the stored full-precision MinSpacing, boundary from the rounded coordinates
        clsh = classify(label, s_sp_col, s_bd, FLOAT_PAD, b_bd)
        stored, rep, ideal = evaluate(xy, row, kind, hrmodel)
        out.append(dict(File=os.path.basename(path), Row=r, Study=study, Role=role, Kind=kind,
                        Algorithm=row["Algorithm"], Dataset=row.get("Dataset", "HR"), Radius=row.get("Radius", 0),
                        Turbines=n, Seed=int(row["Seed"]), Budget=int(row.get("Budget", 6030)),
                        Init=row.get("Init", "random"), Spacing=row.get("Spacing", ""), SMin=smin,
                        Feasible=label, Decimals=d, Eps=eps, SlackSpacing=s_sp, SlackBoundary=s_bd,
                        BoundSpacing=b_sp, BoundBoundary=b_bd, Class=cls,
                        MinSpacingStored=ms, MinSpacingRounded=s_sp + smin - TOL, SlackSpacingStored=s_sp_col,
                        ClassHybrid=clsh, HRModel=hrmodel if kind == "hr" else "",
                        ObjectiveStored=stored, ObjectiveReplayRounded=rep, Ideal=ideal,
                        ReplayDiff=rep - stored, ReplayDiffPP=100.0 * (stored - rep) / ideal,
                        Seconds=float(row["Seconds"])))
    return out


def summarize(rec):
    def block(g):
        lab = g[g.Feasible]
        uns = g[~g.Feasible]
        d = dict(records=int(len(g)), labelled_feasible=int(len(lab)), labelled_infeasible=int(len(uns)),
                 decimals=sorted(int(x) for x in g.Decimals.unique()),
                 class_i=int((g.Class == "i").sum()), class_ii=int((g.Class == "ii").sum()),
                 class_iii=int((g.Class == "iii").sum()),
                 feasible_i=int((lab.Class == "i").sum()), feasible_ii=int((lab.Class == "ii").sum()),
                 feasible_iii=int((lab.Class == "iii").sum()),
                 infeasible_i=int((uns.Class == "i").sum()), infeasible_ii=int((uns.Class == "ii").sum()),
                 infeasible_iii=int((uns.Class == "iii").sum()),
                 hybrid_i=int((g.ClassHybrid == "i").sum()), hybrid_ii=int((g.ClassHybrid == "ii").sum()),
                 hybrid_iii=int((g.ClassHybrid == "iii").sum()),
                 rounded_strict_fail_of_feasible=int(((lab.SlackSpacing < 0) | (lab.SlackBoundary < 0)).sum()),
                 rounded_strict_fail_spacing=int((lab.SlackSpacing < 0).sum()),
                 rounded_strict_fail_boundary=int((lab.SlackBoundary < 0).sum()),
                 stored_minspacing_contradicts_label=int(((g.SlackSpacingStored < 0) & g.Feasible).sum()),
                 max_abs_minspacing_stored_minus_rounded=float((g.MinSpacingStored - g.MinSpacingRounded).abs().max()),
                 minspacing_within_bound=bool(((g.MinSpacingStored - g.MinSpacingRounded).abs() <= g.BoundSpacing).all()),
                 max_abs_replay_diff=float(g.ReplayDiff.abs().max()),
                 max_abs_replay_diff_pp=float(g.ReplayDiffPP.abs().max()),
                 max_abs_replay_diff_feasible=float(lab.ReplayDiff.abs().max()) if len(lab) else None,
                 max_abs_replay_diff_pp_feasible=float(lab.ReplayDiffPP.abs().max()) if len(lab) else None,
                 replay_exact=int((g.ReplayDiff == 0).sum()),
                 cpu_seconds_class_ii_iii=float(g[g.Class != "i"].Seconds.sum()),
                 cpu_seconds_all=float(g.Seconds.sum()))
        if len(g):
            i = g.ReplayDiff.abs().idxmax()
            d["largest_replay"] = dict(File=g.File[i], Algorithm=g.Algorithm[i], Seed=int(g.Seed[i]),
                                       Turbines=int(g.Turbines[i]), Budget=int(g.Budget[i]),
                                       diff=float(g.ReplayDiff[i]), diff_pp=float(g.ReplayDiffPP[i]))
        return d
    res = {}
    for st, g in rec.groupby("Study", sort=False):
        res[st] = block(g.reset_index(drop=True))
        res[st]["role"] = g.Role.iloc[0]
        res[st]["files"] = int(g.File.nunique())
        res[st]["by_kind"] = {k: block(gg.reset_index(drop=True)) for k, gg in g.groupby("Kind")} \
            if g.Kind.nunique() > 1 else None
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", default=HERE)
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    rows = []
    for st, pat, desc, role in STUDIES:
        files = sorted(glob.glob(os.path.join(HERE, pat)))
        if st == "fresh_hr":
            files = [f for f in files if re.match(r"fresh_[bvh]?hr(16|80)\.csv$", os.path.basename(f))]
        assert files, pat
        for f in files:
            rows += audit_file(f, st, role)
        print(st, len(files), "files", flush=True)
    rec = pd.DataFrame(rows)
    drop = ["Eps", "BoundSpacing", "BoundBoundary", "MinSpacingRounded", "ObjectiveReplayRounded"]   # derivable
    rec.drop(columns=drop).to_csv(os.path.join(a.out_dir, "rev3_precision_audit_records.csv"), index=False,
                                  float_format="%.10g")
    summ = summarize(rec)
    cnt = rec[rec.Role != "duplicate"]
    tot = dict(records=int(len(cnt)), original=int(cnt.File.str.match(r"(fresh|mpce)_").sum()),
               additional=int(cnt.File.str.startswith("rev2_").sum()),
               class_i=int((cnt.Class == "i").sum()), class_ii=int((cnt.Class == "ii").sum()),
               class_iii=int((cnt.Class == "iii").sum()),
               hybrid_i=int((cnt.ClassHybrid == "i").sum()), hybrid_ii=int((cnt.ClassHybrid == "ii").sum()),
               hybrid_iii=int((cnt.ClassHybrid == "iii").sum()),
               cpu_hours_class_ii_iii=float(cnt[cnt.Class != "i"].Seconds.sum() / 3600))
    desc = {s[0]: s[2] for s in STUDIES}
    out = dict(description=__doc__.split("\n")[0], tolerance_m=TOL, float_pad_m=FLOAT_PAD,
               bounds="b_sp = 2 sqrt(2) eps, b_bd = sqrt(2) eps, eps = 0.5 10^-d per record",
               total=tot, studies={k: dict(description=desc[k], **v) for k, v in summ.items()})
    with open(os.path.join(a.out_dir, "rev3_precision_audit.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    pd.set_option("display.width", 250)
    t = pd.DataFrame({k: {c: v[c] for c in ("records", "labelled_feasible", "class_i", "class_ii", "class_iii",
                                              "hybrid_ii", "hybrid_iii", "rounded_strict_fail_of_feasible",
                                              "max_abs_replay_diff_pp", "cpu_seconds_class_ii_iii")}
                      for k, v in summ.items()}).T
    print(t.to_string())
    print(json.dumps(tot))


if __name__ == "__main__":
    main()
