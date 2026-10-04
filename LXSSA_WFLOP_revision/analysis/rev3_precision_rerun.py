"""Revision 3, item C4: full-precision reruns of stored records through the UNCHANGED experiment drivers.

Usage
  python3 rev3_precision_rerun.py sample  [--procs=2]                 determinism sample (one record per study x
                                                                       site kind x method; class (ii)/(iii) preferred)
  python3 rev3_precision_rerun.py run STUDY [SHARD NSHARDS] [--procs=2] [--all-classes]
                                                                       rerun the class (ii)/(iii) records of STUDY
                                                                       (rev3_precision_audit_records.csv)
  python3 rev3_precision_rerun.py keys FILE.json [--procs=2]           rerun an explicit list [[File, Row], ...]
  python3 rev3_precision_rerun.py plan                                 write rev3_precision_plan.json
  python3 rev3_precision_rerun.py collect                              write rev3_fullprec_<study>.csv (+ summary)
Results are appended, one JSON line per record, to $REV3_PRECISION_RESDIR/<tag>.jsonl (default
analysis/rev3_precision_results/) (resumable: records already
present are skipped); `collect` merges every results file it finds there.

How a record is rerun. The record (File, Row) of rev3_precision_audit_records.csv is mapped to the task tuple and
runner of the driver that produced it (full_grid_experiments.run_grid / run_hr, mpce_experiments.run_grid / run_hr /
run_grid_x / run_grid_cs, iea37_experiments.run_iea, rev2_ga._call, rev2_gradient.run_iea_g,
rev2_laplace.run_grid_lap, rev2_spacing.run_grid_sp, rev2_site.run_lg). Nothing in the drivers is edited:
  * the final layout is captured by wrapping the module-level name `min_spacing` that each runner calls last, as
    MinSpacing=min_spacing(xy), on the returned layout (the wrapper returns the original value unchanged); the capture
    is verified: the driver's own formatting of the captured layout (.3f / .2f for the legacy drivers, 17 significant
    digits for the rev2 drivers) must equal the Coordinates string of the rerun row, and min_spacing of it must equal
    the rerun MinSpacing bit for bit;
  * the full-precision layout is written with record_io.encode_coordinates (17 significant digits);
  * init_hook.GEN is reset to None before every task (full_grid_experiments.run_grid never sets it);
  * Horns Rev 1 records produced before the binning fix of 2026-09-28 (every Horns Rev record except mpce_hrfix and
    rev2_gahr) are rerun with the pre-fix direction bins (centres 0, 5, ..., 355 deg, nearest sector by np.round),
    set in the worker by overwriting the module constants of hornsrev_model (restored afterwards);
  * a rerun is "bit-identical" if Objective, Calls, Feasible, MinSpacing and the rounded Coordinates and Curve strings
    equal the stored ones exactly (float equality of the parsed CSV values, which pandas wrote with repr precision).
The strict 1e-6 m label is then decided from the full-precision layout with exactly the drivers' formula
(min spacing >= s_min - 1e-6 and max distance outside the circle / parallelogram <= 1e-6), and the objective is
replayed from it with the record's evaluator.
"""
import os, sys, json, time, math, glob, argparse
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
os.chdir(HERE)
RESDIR = os.environ.get("REV3_PRECISION_RESDIR", os.path.join(HERE, "rev3_precision_results"))
TOL = 1e-6

import init_hook
import wflop_model as wm
import hornsrev_model as hr
import iea37_model as iea
import rev2_site_model as lgm
from record_io import encode_coordinates
import full_grid_experiments as fg
import mpce_experiments as mx
import iea37_experiments as ix
import rev2_ga as rga
import rev2_gradient as rgr
import rev2_laplace as rlap
import rev2_spacing as rsp
import rev2_site as rsite
import rev3_precision_audit as aud

# ------------------------------------------------------------------ capture of the final layout
_CAP = {}


def _recorder(orig, reshape=False):
    def ms(xy):
        a = np.array(xy, dtype=float, copy=True)
        _CAP["xy"] = a.reshape(-1, 2) if reshape else a
        return orig(xy)
    ms._orig = orig
    return ms


for _mod in (fg, mx, rsp, rsite):
    _mod.min_spacing = _recorder(_mod.min_spacing)
iea.min_spacing = _recorder(iea.min_spacing, reshape=True)
_ORIG_MS = wm.min_spacing

# ------------------------------------------------------------------ Horns Rev 1 direction binning (old / new)
_HR_NEW = {k: getattr(hr, k) for k in ("WD", "_SEC", "F_WD", "A_WD", "K_WD", "TOWARD", "P_WS")}
_HR_OLD = dict(WD=aud._WD_OLD, _SEC=aud._SEC_OLD, F_WD=aud._F_OLD, A_WD=aud._A_OLD, K_WD=aud._K_OLD,
               TOWARD=aud._TOW_OLD, P_WS=aud._PWS_OLD)


def set_hr_model(model):
    for k, v in (_HR_OLD if model == "old" else _HR_NEW).items():
        setattr(hr, k, v)


# ------------------------------------------------------------------ record -> (runner, task, formatter)
def fmt(xy, d):
    return ";".join(f"{a:.{d}f} {b:.{d}f}" for a, b in xy)


def task_of(rec):
    """(callable, task tuple, coordinate formatter, HR model) for one stored record (a row of the CSV file)."""
    f = rec["File"]
    alg, seed = rec["Algorithm"], int(rec["Seed"])
    n = int(rec["Turbines"])
    ds = str(rec.get("Dataset"))
    budget = int(rec.get("Budget", 6030)) if pd.notna(rec.get("Budget", np.nan)) else 6030
    init = rec.get("Init", "random") if isinstance(rec.get("Init"), str) else "random"
    hrm = "new" if f.startswith(aud.HR_NEW_FILES) else "old"
    full = encode_coordinates
    if f.startswith("fresh_"):
        if ds == "HR" or "Dataset" not in rec or pd.isna(rec.get("Dataset")):
            return fg.run_hr, (alg, n, seed), lambda xy: fmt(xy, 2), hrm
        return fg.run_grid, (alg, int(ds), int(rec["Radius"]), n, seed), lambda xy: fmt(xy, 3), None
    if f.startswith("mpce_iea"):
        return ix.run_iea, (alg, n, seed, budget, init), lambda xy: fmt(xy, 3), None
    if f.startswith("mpce_"):
        if ds == "HR":
            return mx.run_hr, (alg, n, seed, budget, init), lambda xy: fmt(xy, 2), hrm
        t = (alg, int(ds), int(rec["Radius"]), n, seed, budget, init)
        if f.startswith("mpce_csweep"):
            return mx.run_grid_cs, t, lambda xy: fmt(xy, 3), None
        if f.startswith(("mpce_omega90", "mpce_rsdisc")):
            return mx.run_grid_x, t, lambda xy: fmt(xy, 3), None
        return mx.run_grid, t, lambda xy: fmt(xy, 3), None
    if f.startswith(("rev2_ga_", "rev2_gahr", "rev2_gaiea")):
        if ds == "HR":
            fn, t = mx.run_hr, (alg, n, seed, budget, init)
        elif ds.startswith("IEA37"):
            fn, t = ix.run_iea, (alg, n, seed, budget, init)
        else:
            fn, t = mx.run_grid, (alg, int(ds), int(rec["Radius"]), n, seed, budget, init)
        return (lambda task, fn=fn: rga._call((fn, task))), t, full, hrm if ds == "HR" else None
    if f.startswith("rev2_grad"):
        return rgr.run_iea_g, (alg, n, seed, budget, init), full, None
    if f.startswith("rev2_laplace"):
        return rlap.run_grid_lap, (alg, int(ds), int(rec["Radius"]), n, seed, budget, init), full, None
    if f.startswith("rev2_spacing"):
        return rsp.run_grid_sp, (alg, int(ds), int(rec["Radius"]), n, seed, budget, init, rec["Spacing"]), full, None
    if f.startswith("rev2_lg16"):
        return rsite.run_lg, (alg, n, seed, budget, init), full, None
    raise ValueError(f)


# ------------------------------------------------------------------ strict label and replay from a layout
def strict(xy, kind, smin, radius):
    n = len(xy)
    ms = _ORIG_MS(xy) if n > 1 else np.inf
    if kind in ("grid", "iea"):
        viol = np.sqrt((xy ** 2).sum(1)).max() - radius
    else:
        poly = hr.site(n)[1] if kind == "hr" else lgm.site(16)[1]
        viol = (hr.outside_distance(xy, poly) if kind == "hr" else lgm.outside_distance(xy, poly)).max()
        a, b = poly, np.roll(poly, -1, axis=0)
        e = b - a
        nrm = np.c_[e[:, 1], -e[:, 0]] / np.linalg.norm(e, axis=1)[:, None]
        viol_signed = np.einsum("ijk,jk->ij", xy[:, None, :] - a[None], nrm).max()
        return ms, ms - (smin - TOL), TOL - viol_signed, bool((ms >= smin - TOL) and (viol <= TOL))
    return ms, ms - (smin - TOL), TOL - viol, bool((ms >= smin - TOL) and (viol <= TOL))


def replay(xy, kind, ds, hrm):
    if kind == "grid":
        return wm.farm_objective(xy, int(ds))[0]
    if kind == "iea":
        return iea.aep(xy)
    if kind == "lg":
        return lgm.aep_gwh(xy)
    return aud.hr_aep(xy, hrm)


# ------------------------------------------------------------------ one record
_STORED = {}


def stored_row(fname, row):
    if fname not in _STORED:
        _STORED[fname] = pd.read_csv(os.path.join(HERE, fname), dtype={"Dataset": str})
    return _STORED[fname].iloc[int(row)].to_dict()


def rerun(key):
    fname, row = key
    s = stored_row(fname, row)
    s["File"] = fname
    fn, task, form, hrm = task_of(s)
    init_hook.GEN = None
    _CAP.clear()
    if hrm is not None:
        set_hr_model(hrm)
    t0 = time.perf_counter()
    try:
        out = fn(task)
    finally:
        set_hr_model("new")
    sec = time.perf_counter() - t0
    xy = _CAP["xy"]
    cap_ok = (form(xy) == out["Coordinates"]) and (_ORIG_MS(xy) == out["MinSpacing"])
    kind = aud.record_kind(s, fname)
    smin = aud.smin_of(s, kind)
    radius = float(s.get("Radius", 0) or 0)
    ms, ssp, sbd, lab = strict(xy, kind, smin, radius)
    st_obj = float(s["AEP"]) if "AEP" in s else float(s["Objective"])
    re_obj = float(out["AEP"]) if "AEP" in out else float(out["Objective"])
    rep = float(replay(xy, kind, s.get("Dataset"), hrm or "new"))
    d = 2 if (kind in ("hr", "lg")) else 3
    rounded = fmt(xy, d)
    res = dict(File=fname, Row=int(row), Algorithm=s["Algorithm"], Dataset=str(s.get("Dataset", "HR")),
               Radius=s.get("Radius", 0), Turbines=int(s["Turbines"]), Seed=int(s["Seed"]),
               Budget=int(s.get("Budget", 6030)), Init=s.get("Init", "random"), Spacing=s.get("Spacing", ""),
               SMin=smin, Kind=kind, HRModel=hrm or "",
               CoordinatesFull=encode_coordinates(xy), CaptureVerified=bool(cap_ok),
               MinSpacingFull=ms, SpacingSlackFull=ssp, BoundarySlackFull=sbd, StrictLabel=lab,
               StoredLabel=bool(s["Feasible"]), RerunLabel=bool(out["Feasible"]),
               ObjectiveStored=st_obj, ObjectiveRerun=re_obj, ObjectiveReplayFull=rep,
               BitIdentical=bool(re_obj == st_obj), ReplayBitIdentical=bool(rep == st_obj),
               CallsMatch=bool(int(out["Calls"]) == int(s["Calls"])),
               MinSpacingMatch=bool(float(out["MinSpacing"]) == float(s["MinSpacing"])),
               RoundedCoordsMatch=bool(rounded == s["Coordinates"]),
               CurveMatch=bool(out["Curve"] == s["Curve"]) if "Curve" in s else None,
               SecondsStored=float(s["Seconds"]), SecondsRerun=sec)
    res["AllMatch"] = bool(res["BitIdentical"] and res["CallsMatch"] and res["MinSpacingMatch"] and
                           res["RoundedCoordsMatch"] and res["CurveMatch"] is not False and
                           res["RerunLabel"] == res["StoredLabel"] and cap_ok)
    return res


def _safe(key):
    try:
        return rerun(key)
    except Exception as e:                                   # report, never hide
        import traceback
        return dict(File=key[0], Row=int(key[1]), Error=repr(e), Trace=traceback.format_exc()[-2000:])


def run_keys(keys, tag, procs):
    os.makedirs(RESDIR, exist_ok=True)
    path = os.path.join(RESDIR, f"{tag}.jsonl")
    done = set()
    for p in glob.glob(os.path.join(RESDIR, "*.jsonl")):
        for line in open(p):
            r = json.loads(line)
            if "Error" not in r:
                done.add((r["File"], r["Row"]))
    todo = [k for k in keys if (k[0], int(k[1])) not in done]
    print(f"{tag}: {len(keys)} keys, {len(todo)} to run, procs={procs}", flush=True)
    t0 = time.time()
    from multiprocessing import Pool
    with Pool(procs, maxtasksperchild=50) as pool, open(path, "a") as fh:
        for i, r in enumerate(pool.imap_unordered(_safe, todo, chunksize=1), 1):
            fh.write(json.dumps(r, default=lambda o: o.item() if hasattr(o, "item") else str(o)) + "\n")
            fh.flush()
            if i % 10 == 0 or i == len(todo) or "Error" in r:
                print(f"  {i}/{len(todo)} {time.time() - t0:.0f}s last {r['File']}:{r['Row']} "
                      f"{'ERROR ' + r['Error'] if 'Error' in r else ('match' if r['AllMatch'] else 'MISMATCH')}",
                      flush=True)
    print(f"{tag}: done in {time.time() - t0:.0f}s", flush=True)


def audit_records():
    return pd.read_csv(os.path.join(HERE, "rev3_precision_audit_records.csv"), low_memory=False)


def sample_keys():
    r = audit_records()
    r = r[r.Role != "duplicate"].copy()
    r["pref"] = (r.Class == "i").astype(int)
    keys = []
    for (st, kind, alg), g in r.groupby(["Study", "Kind", "Algorithm"]):
        g = g.sort_values(["pref", "Seconds"])
        x = g.iloc[0]
        keys.append((x.File, int(x.Row)))
    return keys


def study_keys(study, all_classes=False):
    r = audit_records()
    r = r[(r.Study == study)]
    if not all_classes:
        r = r[r.Class != "i"]
    return [(f, int(i)) for f, i in zip(r.File, r.Row)]


# ------------------------------------------------------------------ plan
ORDER_PRIORITY = ["rev2_grad", "rev2_lg16", "rev2_lg16b", "rev2_laplace", "rev2_spacing"]


def plan():
    r = audit_records()
    r = r[(r.Role != "duplicate") & (r.Class != "i")]
    done = set()
    for p in glob.glob(os.path.join(RESDIR, "*.jsonl")):
        for line in open(p):
            q = json.loads(line)
            if "Error" not in q:
                done.add((q["File"], q["Row"]))
    studies = []
    for st, g in r.groupby("Study"):
        keys = [(f, int(i)) for f, i in zip(g.File, g.Row)]
        rem = [k for k in keys if k not in done]
        gr = g[[(f, int(i)) not in done for f, i in zip(g.File, g.Row)]]
        studies.append(dict(study=st, role=g.Role.iloc[0], records_class_ii_iii=len(keys), done_locally=len(keys) - len(rem),
                            remaining=len(rem), cpu_hours_stored_seconds=round(g.Seconds.sum() / 3600, 3),
                            cpu_hours_remaining=round(gr.Seconds.sum() / 3600, 3),
                            command=f"python3 rev3_precision_rerun.py run {st} SHARD NSHARDS --procs=P",
                            remaining_keys=rem))
    studies.sort(key=lambda d: d["cpu_hours_remaining"])
    tot = sum(d["cpu_hours_remaining"] for d in studies)
    out = dict(note=("Rerun plan for the class (ii)/(iii) records of rev3_precision_audit_records.csv. Each record is "
                     "identified by (File, Row): the 0-based data row of the stored CSV in analysis/. A worker runs "
                     "'python3 rev3_precision_rerun.py run STUDY SHARD NSHARDS --procs=P' from analysis/ (shard = "
                     "records[SHARD::NSHARDS] of the study's remaining list in this file order) and returns "
                     "analysis/rev3_precision_results/<study>_s<SHARD>of<NSHARDS>.jsonl; the lead copies the jsonl files "
                     "back and runs 'python3 rev3_precision_rerun.py collect'. CPU estimates use the stored Seconds "
                     "column (measured on 4 busy cores); local reruns ran at roughly the speed given in "
                     "rev3_precision_summary.json."),
               environment=dict(python=sys.version.split()[0], numpy=np.__version__,
                                scipy=__import__("scipy").__version__, pandas=pd.__version__),
               remaining_cpu_hours=round(tot, 2), studies=studies)
    with open(os.path.join(HERE, "rev3_precision_plan.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    for d in studies:
        print(f"{d['study']:16s} {d['records_class_ii_iii']:5d} done {d['done_locally']:5d} remaining {d['remaining']:5d} "
              f"{d['cpu_hours_remaining']:7.2f} h")
    print("remaining CPU-hours", round(tot, 2))


# ------------------------------------------------------------------ collect
def collect():
    rows = []
    for p in sorted(glob.glob(os.path.join(RESDIR, "*.jsonl"))):
        for line in open(p):
            rows.append(json.loads(line))
    res = pd.DataFrame(rows)
    err = res[res.get("Error").notna()] if "Error" in res else res.iloc[:0]
    res = res[res.get("Error").isna()] if "Error" in res else res
    res = res.drop_duplicates(["File", "Row"], keep="last")
    r = audit_records().set_index(["File", "Row"])
    res["Study"] = [r.loc[(f, i), "Study"] for f, i in zip(res.File, res.Row)]
    res["Class"] = [r.loc[(f, i), "Class"] for f, i in zip(res.File, res.Row)]
    res["ReplayDiffRounded"] = [r.loc[(f, i), "ReplayDiff"] for f, i in zip(res.File, res.Row)]
    res["LabelChanged"] = res.StrictLabel != res.StoredLabel
    res["ReplayDiffFull"] = res.ObjectiveReplayFull - res.ObjectiveStored
    cols = ["File", "Row", "Study", "Algorithm", "Dataset", "Radius", "Turbines", "Seed", "Budget", "Init", "Spacing",
            "SMin", "Class", "CoordinatesFull", "MinSpacingFull", "SpacingSlackFull", "BoundarySlackFull",
            "StrictLabel", "StoredLabel", "RerunLabel", "LabelChanged", "ObjectiveReplayFull", "ObjectiveStored",
            "ObjectiveRerun", "BitIdentical", "ReplayBitIdentical", "ReplayDiffFull", "ReplayDiffRounded",
            "CallsMatch", "MinSpacingMatch", "RoundedCoordsMatch", "CurveMatch", "CaptureVerified", "AllMatch",
            "HRModel", "SecondsStored", "SecondsRerun"]
    summ = {}
    for st, g in res.groupby("Study"):
        g = g.sort_values(["File", "Row"])
        g[cols].to_csv(os.path.join(HERE, f"rev3_fullprec_{st}.csv"), index=False, float_format="%.17g")
        summ[st] = dict(reruns=int(len(g)), class_iii_reruns=int((g.Class != "i").sum()),
                        bit_identical=int(g.BitIdentical.sum()), all_match=int(g.AllMatch.sum()),
                        capture_verified=int(g.CaptureVerified.sum()),
                        strict_confirmed=int((~g.LabelChanged).sum()), strict_changed=int(g.LabelChanged.sum()),
                        changed=[dict(File=a, Row=int(b), Algorithm=c, Seed=int(d), Stored=bool(e), Strict=bool(h),
                                      SpacingSlack=float(i), BoundarySlack=float(j))
                                 for a, b, c, d, e, h, i, j in g[g.LabelChanged][["File", "Row", "Algorithm", "Seed",
                                     "StoredLabel", "StrictLabel", "SpacingSlackFull", "BoundarySlackFull"]].itertuples(index=False)],
                        replay_full_bit_identical=int(g.ReplayBitIdentical.sum()),
                        max_abs_replay_diff_full=float(g.ReplayDiffFull.abs().max()),
                        min_boundary_slack_full=float(g.BoundarySlackFull.min()),
                        min_spacing_slack_full=float(g.SpacingSlackFull.min()),
                        seconds_rerun=float(g.SecondsRerun.sum()), seconds_stored=float(g.SecondsStored.sum()),
                        mismatches=[dict(File=a, Row=int(b)) for a, b in g[~g.AllMatch][["File", "Row"]].itertuples(index=False)])
    out = dict(studies=summ, errors=[dict(File=a, Row=int(b), Error=c) for a, b, c in
                                     err[["File", "Row", "Error"]].itertuples(index=False)] if len(err) else [])
    with open(os.path.join(HERE, "rev3_precision_rerun_summary.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    for st, d in summ.items():
        print(f"{st:16s} reruns {d['reruns']:5d} bit-identical {d['bit_identical']:5d} all-match {d['all_match']:5d} "
              f"strict changed {d['strict_changed']:3d} rerun s {d['seconds_rerun']:8.0f} stored s {d['seconds_stored']:8.0f}")
    print("errors:", len(out["errors"]))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd")
    ap.add_argument("arg", nargs="*")
    ap.add_argument("--procs", type=int, default=2)
    ap.add_argument("--all-classes", action="store_true")
    a = ap.parse_args()
    if a.cmd == "sample":
        run_keys(sample_keys(), "sample", a.procs)
    elif a.cmd == "run":
        st = a.arg[0]
        keys = study_keys(st, a.all_classes)
        tag = st
        if len(a.arg) >= 3:
            sh, nsh = int(a.arg[1]), int(a.arg[2])
            keys = keys[sh::nsh]; tag = f"{st}_s{sh}of{nsh}"
        run_keys(keys, tag, a.procs)
    elif a.cmd == "keys":
        run_keys([tuple(k) for k in json.load(open(a.arg[0]))], os.path.splitext(os.path.basename(a.arg[0]))[0], a.procs)
    elif a.cmd == "plan":
        plan()
    elif a.cmd == "collect":
        collect()
    else:
        raise SystemExit(__doc__)
