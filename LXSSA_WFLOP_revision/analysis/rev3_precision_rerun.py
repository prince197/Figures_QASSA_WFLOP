"""Revision 3, item C4: full-precision reruns of stored records through the UNCHANGED experiment drivers.

Usage
  python3 rev3_precision_rerun.py sample  [--procs=2]                 determinism sample (one record per study x
                                                                       site kind x method; class (ii)/(iii) preferred)
  python3 rev3_precision_rerun.py run STUDY [--procs=2] [--all-classes]
                                                                       rerun the class (ii)/(iii) records of STUDY
                                                                       (rev3_precision_audit_records.csv)
  python3 rev3_precision_rerun.py run STUDY SHARD NSHARDS [--procs=1]  rerun shard SHARD of the plan (explicit key
                                                                       list in rev3_precision_plan.json); resumable;
                                                                       writes rev3_fullprec_STUDY_sSHARDofNSHARDS.csv
  python3 rev3_precision_rerun.py keys FILE.json [--procs=2]           rerun an explicit list [[File, Row], ...]
  python3 rev3_precision_rerun.py plan [STUDY,STUDY]                   write rev3_precision_plan.json (studies listed
                                                                       are finished locally and get no shards)
  python3 rev3_precision_rerun.py export-local                         write $REV3_PRECISION_RESDIR/*.jsonl as
                                                                       rev3_fullprec_<study>_local.csv
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
for _v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):      # one thread per worker process
    os.environ.setdefault(_v, "1")
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
DRIVER = {
    "fresh_grid": "full_grid_experiments.py grid (run_grid)", "fresh_bgrid": "full_grid_experiments.py bgrid",
    "fresh_vgrid": "full_grid_experiments.py vgrid", "fresh_bsplit": "full_grid_experiments.py bsplit",
    "fresh_hgrid": "full_grid_experiments.py hgrid", "fresh_hsplit": "full_grid_experiments.py hsplit",
    "fresh_hr": "full_grid_experiments.py hr16/hr80/bhr16/bhr80/vhr16/vhr80/hhr16/hhr80 (run_hr, pre-fix HR bins)",
    "mpce_psobv": "mpce_experiments.py psobv", "mpce_psoc": "mpce_experiments.py psoc",
    "mpce_slsqp": "mpce_experiments.py slsqp", "mpce_rsvns": "mpce_experiments.py rsvns",
    "mpce_rsdisc": "mpce_experiments.py rsdisc (run_grid_x)", "mpce_psosplit": "mpce_experiments.py psosplit",
    "mpce_omega90": "mpce_experiments.py omega90 (run_grid_x)", "mpce_csweep": "mpce_experiments.py csweep (run_grid_cs)",
    "mpce_feas": "mpce_experiments.py feasx 0..7 8", "mpce_feasp": "mpce_experiments.py feasp",
    "mpce_b30k": "mpce_experiments.py b30k 0..2 3", "mpce_b30kp": "mpce_experiments.py b30kp",
    "mpce_b120k": "mpce_experiments.py b120k (8 shards)", "mpce_b120kp": "mpce_experiments.py b120kp (4 shards)",
    "mpce_hr16new": "mpce_experiments.py hr16new", "mpce_hrfix": "mpce_experiments.py hrfix 0..23 24",
    "mpce_iea": "iea37_experiments.py iea16/iea36/iea16p/iea36p", "rev2_ga": "rev2_ga.py ga",
    "rev2_gahr": "rev2_ga.py gahr", "rev2_gaiea": "rev2_ga.py gaiea", "rev2_grad": "rev2_gradient.py grad",
    "rev2_laplace": "rev2_laplace.py laplace", "rev2_spacing": "rev2_spacing.py spacing",
    "rev2_lg16": "rev2_site.py lg16", "rev2_lg16b": "rev2_site.py lg16b"}
ORDER_PRIORITY = ["rev2_grad", "rev2_lg16", "rev2_lg16b", "rev2_laplace", "rev2_spacing"]


def shard_files():
    """Raw rerun files in analysis/: cloud shards rev3_fullprec_<study>_s<i>of<k>.csv and the local reruns
    rev3_fullprec_<study>_local.csv (not the per-study summaries rev3_fullprec_<study>.csv written by collect)."""
    import re
    return sorted(p for p in glob.glob(os.path.join(HERE, "rev3_fullprec_*.csv"))
                  if re.search(r"_(s\d+of\d+|local)\.csv$", os.path.basename(p)))


def export_local():
    """Write the results of RESDIR (*.jsonl) as raw files rev3_fullprec_<study>_local.csv in analysis/."""
    rows = []
    for p in sorted(glob.glob(os.path.join(RESDIR, "*.jsonl"))):
        for line in open(p):
            q = json.loads(line)
            if "Error" not in q:
                rows.append(q)
    df = pd.DataFrame(rows).drop_duplicates(["File", "Row"], keep="last")
    r = audit_records().set_index(["File", "Row"])
    df["_study"] = r.loc[list(zip(df.File, df.Row)), "Study"].values
    for st, g in df.groupby("_study"):
        fn = os.path.join(HERE, f"rev3_fullprec_{st}_local.csv")
        g.drop(columns="_study").sort_values(["File", "Row"]).to_csv(fn, index=False)
        print(f"wrote {os.path.basename(fn)}: {len(g)}")


def _done_keys():
    done = set()
    for p in glob.glob(os.path.join(RESDIR, "*.jsonl")):
        for line in open(p):
            q = json.loads(line)
            if "Error" not in q:
                done.add((q["File"], int(q["Row"])))
    for p in shard_files():
        q = pd.read_csv(p, usecols=["File", "Row"])
        done.update((f, int(i)) for f, i in zip(q.File, q.Row))
    return done


def plan(exclude=(), target_hours=1.0):
    """Cost-balanced explicit shards (LPT packing on the stored Seconds, about target_hours CPU-hours per shard)
    of the class (ii)/(iii) records not yet rerun; studies in `exclude` are being finished locally."""
    r = audit_records()
    r = r[(r.Role != "duplicate") & (r.Class != "i")]
    done = _done_keys()
    studies = []
    for st, g in r.groupby("Study"):
        keys = [(f, int(i)) for f, i in zip(g.File, g.Row)]
        gr = g[[(f, int(i)) not in done for f, i in zip(g.File, g.Row)]].sort_values("Seconds", ascending=False)
        hrs = gr.Seconds.sum() / 3600
        d = dict(study=st, role=g.Role.iloc[0], driver=DRIVER.get(st, ""), records_class_ii_iii=len(keys),
                 done_locally=len(keys) - len(gr), remaining=len(gr), cpu_hours_remaining=round(hrs, 3),
                 local=st in exclude, shards=[], shard_hours=[], shard_commands=[], shard_outputs=[])
        if len(gr) and st not in exclude:
            nsh = max(1, int(math.ceil(hrs / target_hours)))
            bins = [[] for _ in range(nsh)]; load = [0.0] * nsh
            for f, i, sec in zip(gr.File, gr.Row, gr.Seconds):              # longest first into the lightest shard
                j = int(np.argmin(load)); bins[j].append([f, int(i)]); load[j] += sec
            d["shards"] = bins
            d["shard_hours"] = [round(x / 3600, 3) for x in load]
            d["shard_commands"] = [f"python3 rev3_precision_rerun.py run {st} {j} {nsh} --procs=1" for j in range(nsh)]
            d["shard_outputs"] = [f"analysis/rev3_fullprec_{st}_s{j}of{nsh}.csv" for j in range(nsh)]
        studies.append(d)
    studies.sort(key=lambda d: -d["cpu_hours_remaining"])
    cloud = [d for d in studies if d["shards"]]
    out = dict(note=("Rerun plan for the class (ii)/(iii) records of rev3_precision_audit_records.csv not yet rerun. "
                     "A record is (File, Row) = 0-based data row of the stored CSV in analysis/. Each shard is an explicit "
                     "cost-balanced key list (LPT packing on the stored Seconds). A worker, in a checkout of analysis/ "
                     "(needs the stored run CSVs, rev3_precision_audit_records.csv, this plan and the model/optimizer "
                     "modules), runs the shard command from analysis/; it is resumable (finished records, kept one JSON "
                     "line each in analysis/rev3_precision_results/<study>_s<i>of<k>.jsonl, are skipped on restart) and "
                     "ends by writing analysis/rev3_fullprec_<study>_s<i>of<k>.csv. Copy the shard CSVs back into analysis/ "
                     "and run 'python3 rev3_precision_rerun.py collect' and 'python3 rev3_precision_rerun.py table'. "
                     "Hours = sum of the stored Seconds (originally measured with 4 processes on 4 cores); on an idle core "
                     "a rerun took 1.0-1.4x the stored time, on this shared machine 2-4x. SLSQP-based methods "
                     "(rev2_grad, MS-SLSQP rows) are platform dependent: expect 'near'/'diverged' reruns."),
               environment=dict(python=sys.version.split()[0], numpy=np.__version__,
                                scipy=__import__("scipy").__version__, pandas=pd.__version__,
                                env="OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 recommended"),
               remaining_cpu_hours_cloud=round(sum(d["cpu_hours_remaining"] for d in cloud), 2),
               n_shards_cloud=sum(len(d["shards"]) for d in cloud),
               all_shard_commands=[c for d in cloud for c in d["shard_commands"]], studies=studies)
    with open(os.path.join(HERE, "rev3_precision_plan.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    for d in studies:
        print(f"{d['study']:16s} {d['records_class_ii_iii']:5d} done {d['done_locally']:5d} remaining {d['remaining']:5d} "
              f"{d['cpu_hours_remaining']:7.2f} h shards {len(d['shards'])} {'LOCAL' if d['local'] else ''}")
    print("cloud CPU-hours", out["remaining_cpu_hours_cloud"], "shards", out["n_shards_cloud"])


def write_shard_csv(st, sh, nsh, keys):
    rows = []
    for p in glob.glob(os.path.join(RESDIR, "*.jsonl")):
        for line in open(p):
            q = json.loads(line)
            if "Error" not in q:
                rows.append(q)
    want = {(f, int(i)) for f, i in keys}
    df = pd.DataFrame([q for q in rows if (q["File"], int(q["Row"])) in want]).drop_duplicates(["File", "Row"], keep="last")
    fn = os.path.join(HERE, f"rev3_fullprec_{st}_s{sh}of{nsh}.csv")
    df.sort_values(["File", "Row"]).to_csv(fn, index=False)
    print(f"wrote {fn}: {len(df)} of {len(want)} records", flush=True)


# ------------------------------------------------------------------ collect
DELTA = 1e-8          # margin (m) required of a strict-label slack when the rerun is reproduced only to float noise


def _category(g):
    same = g.RoundedCoordsMatch & g.CallsMatch & (g.CurveMatch != False) & (g.RerunLabel == g.StoredLabel) & g.CaptureVerified
    bit = same & g.BitIdentical & g.MinSpacingMatch
    rel = (g.ObjectiveRerun - g.ObjectiveStored).abs() / g.ObjectiveStored.abs()
    noise = same & ~bit & (rel <= 1e-9) & ((g.MinSpacingFull - g.MinSpacingStoredCol).abs() <= 1e-6)
    near = (g.RoundedCoordsMatch & (g.RerunLabel == g.StoredLabel) & g.CaptureVerified & ~bit & ~noise & (rel <= 1e-9)
            & ((g.MinSpacingFull - g.MinSpacingStoredCol).abs() <= 1e-6))
    return np.where(bit, "bit", np.where(noise, "noise", np.where(near, "near", "diverged")))


def load_results():
    rows = []
    for p in sorted(glob.glob(os.path.join(RESDIR, "*.jsonl"))):
        for line in open(p):
            rows.append(json.loads(line))
    res = pd.DataFrame(rows)
    shard = [pd.read_csv(p, dtype={"Dataset": str, "Spacing": str, "Init": str}) for p in shard_files()]
    if shard:
        res = pd.concat([res] + shard, ignore_index=True)
    for c in ("CurveMatch",):
        res[c] = res[c].map(lambda v: None if (v is None or (isinstance(v, float) and np.isnan(v))) else
                            (v if isinstance(v, (bool, np.bool_)) else str(v) == "True"))
    if "Error" not in res:
        res["Error"] = np.nan
    err = res[res.Error.notna()]
    res = res[res.Error.isna()].drop_duplicates(["File", "Row"], keep="last").copy()
    r = audit_records().set_index(["File", "Row"])
    idx = list(zip(res.File, res.Row))
    for c, cc in (("Study", "Study"), ("Class", "Class"), ("ReplayDiffRounded", "ReplayDiff"),
                  ("MinSpacingStoredCol", "MinSpacingStored"), ("Ideal", "Ideal")):
        res[c] = r.loc[idx, cc].values
    res["LabelChanged"] = res.StrictLabel != res.StoredLabel
    res["ReplayDiffFull"] = res.ObjectiveReplayFull - res.ObjectiveStored
    res["Reproduction"] = _category(res)
    m_ok = np.minimum(res.SpacingSlackFull, res.BoundarySlackFull)            # > 0: strictly feasible by m_ok
    m_bad = -np.minimum(res.SpacingSlackFull, res.BoundarySlackFull)         # > 0: infeasible by m_bad
    res["Margin"] = np.where(res.StrictLabel, m_ok, m_bad)
    res["LabelDecided"] = (res.Reproduction == "bit") | ((res.Reproduction == "noise") & (res.Margin > DELTA))
    from record_io import decode_coordinates
    res["WakeFlips"] = [wake_flips(np.array(decode_coordinates(c)), k, h) for c, k, h in
                        zip(res.CoordinatesFull, res.Kind, res.HRModel)]
    return res, err


def wake_flips(full, kind, hrm):
    """Number of (direction, turbine pair) in-wake indicators that differ between the full-precision layout and its
    rounding to the stored decimals (top-hat models: benchmark Jensen, Horns Rev 1, Lillgrund; NaN for the IEA37
    Gaussian wake, which has no in-wake switch)."""
    if kind == "iea":
        return np.nan
    d = 3 if kind == "grid" else 2
    rnd = np.array([[float(f"{v:.{d}f}") for v in p] for p in full])

    def ind(xy):
        dx = xy[:, 0][:, None] - xy[:, 0][None]; dy = xy[:, 1][:, None] - xy[:, 1][None]
        if kind == "grid":
            th = wm.THETA; c, s = np.cos(th)[:, None, None], np.sin(th)[:, None, None]; A = wm.R / wm.K
            beta = np.arccos(np.clip((dx * c + dy * s + A) / np.sqrt((dx + A * c) ** 2 + (dy + A * s) ** 2), -1, 1))
            return beta < np.arctan(wm.K)
        if kind == "lg":
            tow, rr, kw = np.deg2rad(270.0 - lgm.WD), lgm.RR, lgm.KW
        else:
            tow, rr, kw = (aud._TOW_OLD if hrm == "old" else hr.TOWARD), hr.RR, hr.KW
        c = np.cos(tow)[:, None, None]; s = np.sin(tow)[:, None, None]
        x = dx * c + dy * s
        return (x > 0) & (np.abs(-dx * s + dy * c) < rr + kw * x)
    return int((ind(full) != ind(rnd)).sum())


COLS = ["File", "Row", "Study", "Algorithm", "Dataset", "Radius", "Turbines", "Seed", "Budget", "Init", "Spacing",
        "SMin", "Class", "CoordinatesFull", "MinSpacingFull", "SpacingSlackFull", "BoundarySlackFull", "Margin",
        "StrictLabel", "StoredLabel", "RerunLabel", "LabelChanged", "LabelDecided", "ObjectiveReplayFull",
        "ObjectiveStored", "ObjectiveRerun", "BitIdentical", "Reproduction", "ReplayBitIdentical", "ReplayDiffFull",
        "ReplayDiffRounded", "WakeFlips", "CallsMatch", "MinSpacingMatch", "RoundedCoordsMatch", "CurveMatch", "CaptureVerified",
        "AllMatch", "HRModel", "SecondsStored", "SecondsRerun"]


def collect():
    res, err = load_results()
    summ = {}
    for st, g in res.groupby("Study"):
        g = g.sort_values(["File", "Row"])
        g[COLS].to_csv(os.path.join(HERE, f"rev3_fullprec_{st}.csv"), index=False, float_format="%.17g")
        ch = g[g.LabelChanged]
        summ[st] = dict(reruns=int(len(g)), class_ii_iii_reruns=int((g.Class != "i").sum()),
                        bit_identical=int((g.Reproduction == "bit").sum()),
                        float_noise=int((g.Reproduction == "noise").sum()),
                        near=int((g.Reproduction == "near").sum()),
                        diverged=int((g.Reproduction == "diverged").sum()),
                        objective_bit_identical=int(g.BitIdentical.sum()),
                        label_decided=int(g.LabelDecided.sum()),
                        strict_confirmed=int((g.LabelDecided & ~g.LabelChanged).sum()),
                        strict_changed=int((g.LabelDecided & g.LabelChanged).sum()),
                        strict_changed_undecided=int((~g.LabelDecided & g.LabelChanged).sum()),
                        changed=[dict(File=a, Row=int(b), Algorithm=c, Seed=int(d), Stored=bool(e), Strict=bool(f),
                                      SpacingSlack=float(h), BoundarySlack=float(i), Reproduction=j)
                                 for a, b, c, d, e, f, h, i, j in ch[["File", "Row", "Algorithm", "Seed", "StoredLabel",
                                     "StrictLabel", "SpacingSlackFull", "BoundarySlackFull", "Reproduction"]].itertuples(index=False)],
                        min_margin=float(g.Margin.min()),
                        margins_below_1e_8=int((g.Margin <= DELTA).sum()),
                        replay_full_bit_identical=int(g.ReplayBitIdentical.sum()),
                        replay_rounded_nonzero=int((g.ReplayDiffRounded != 0).sum()),
                        replay_rounded_gt_1e_4_pp=int((100 * g.ReplayDiffRounded.abs() / g.Ideal > 1e-4).sum()),
                        replay_rounded_gt_1e_4_pp_with_wake_flip=int(((100 * g.ReplayDiffRounded.abs() / g.Ideal > 1e-4)
                                                                      & (g.WakeFlips > 0)).sum()),
                        max_abs_replay_diff_rounded_pp=float((100 * g.ReplayDiffRounded / g.Ideal).abs().max()),
                        max_abs_replay_diff_full=float(g.ReplayDiffFull.abs().max()),
                        max_abs_replay_diff_full_pp=float((100 * g.ReplayDiffFull / g.Ideal).abs().max()),
                        seconds_rerun=float(g.SecondsRerun.sum()), seconds_stored=float(g.SecondsStored.sum()),
                        not_bit_identical=[dict(File=a, Row=int(b), Algorithm=c, Seed=int(d), Reproduction=e,
                                                RelObjDiff=float((f - h) / abs(h)))
                                           for a, b, c, d, e, f, h in g[g.Reproduction != "bit"][["File", "Row", "Algorithm",
                                               "Seed", "Reproduction", "ObjectiveRerun", "ObjectiveStored"]].itertuples(index=False)])
    out = dict(delta_m=DELTA, studies=summ,
               errors=[dict(File=a, Row=int(b), Error=c) for a, b, c in err[["File", "Row", "Error"]].itertuples(index=False)])
    with open(os.path.join(HERE, "rev3_precision_rerun_summary.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    for st, d in summ.items():
        print(f"{st:16s} reruns {d['reruns']:5d} bit {d['bit_identical']:5d} noise {d['float_noise']:4d} diverged "
              f"{d['diverged']:4d} decided {d['label_decided']:5d} changed {d['strict_changed']:3d} "
              f"rerun s {d['seconds_rerun']:8.0f} stored s {d['seconds_stored']:8.0f}")
    print("errors:", len(out["errors"]))
    return res


# ------------------------------------------------------------------ LaTeX table
SHORT = {
    "fresh_grid": "Bench., 6 methods$^a$", "fresh_bgrid": "Bench., (LX-)SSA-VNS",
    "fresh_vgrid": "Bench., VNS", "fresh_bsplit": "Split, LX-SSA-VNS", "fresh_hgrid": "Bench., LX-VNS$^s$",
    "fresh_hsplit": "Split, LX-VNS$^s$", "fresh_hr": "HR1, earlier$^{s,h}$",
    "mpce_psobv": "Bench., PSO-VNS$^h$", "mpce_psoc": "Bench., PSO", "mpce_slsqp": "Bench., MS-SLSQP",
    "mpce_rsvns": "Bench., RS-VNS", "mpce_rsdisc": "Bench., RSD-VNS", "mpce_psosplit": "Split, PSO-VNS",
    "mpce_omega90": "Split, PSO-VNS, $\\omega=0.9$", "mpce_csweep": "PSO coeff. sweep",
    "mpce_feas": "Feas. init., 8 methods$^h$", "mpce_feasp": "Feas. init., PSO-VNS$^h$",
    "mpce_b30k": "30,030 eval., 9 methods$^h$", "mpce_b30kp": "30,030 eval., PSO-VNS$^h$",
    "mpce_b120k": "120,030 eval., 9 methods$^h$", "mpce_b120kp": "120,030 eval., PSO-VNS$^h$",
    "mpce_hr16new": "HR1, PSO, RS-VNS$^{s,h}$", "mpce_hrfix": "HR1, 10 methods",
    "mpce_iea": "IEA37, 10 methods", "rev2_ga": "GA, bench.", "rev2_gahr": "GA, HR1",
    "rev2_gaiea": "GA, IEA37", "rev2_grad": "Exact gradients, IEA37", "rev2_laplace": "Laplace ablation",
    "rev2_spacing": "Spacing $5D$/$6D$", "rev2_lg16": "Lillgrund, 6,030", "rev2_lg16b": "Lillgrund, 30,030",
}


def table():
    A = json.load(open(os.path.join(HERE, "rev3_precision_audit.json")))["studies"]
    res, _ = load_results()
    lines = []
    tot = np.zeros(12, dtype=float)
    mx_all = 0.0
    for st, d in A.items():
        if d["role"] == "duplicate":
            continue
        g = res[res.Study == st]
        g3 = g[g.Class != "i"]
        nb = int((g3.Reproduction == "bit").sum()); nn = int((g3.Reproduction == "noise").sum())
        nr = int((g3.Reproduction == "near").sum())
        conf = int((g3.LabelDecided & ~g3.LabelChanged).sum()); chg = int((g3.LabelDecided & g3.LabelChanged).sum())
        open_ = d["class_ii"] + d["class_iii"] - conf - chg
        v = [d["records"], d["labelled_feasible"], d["class_i"], d["class_ii"], d["class_iii"], len(g3), nb, nn, nr,
             conf, chg, open_]
        tot += np.array(v, dtype=float)
        mx = d["max_abs_replay_diff_pp"]; mx_all = max(mx_all, mx)
        lines.append(SHORT[st] + " & " + " & ".join(f"{int(x):,}" for x in v) + f" & {mx:.3f} \\\\")
    body = "\n".join(lines)
    tex = r"""\begin{table}[!htbp]
\centering
\caption{Verification of the strict feasibility labels ($10^{-6}$~m) of all stored run records with coordinates
(original and revision studies; all cases, methods, seeds 1--30 or 1--10 and budgets as stored).
Class (i): label confirmed from the rounded coordinates with slack beyond the rounding bound
($2\sqrt2\,\varepsilon$ for spacing, $\sqrt2\,\varepsilon$ for the boundary, $\varepsilon=0.5\cdot10^{-d}$~m for $d$
stored decimals); (ii): contradicted beyond the bound; (iii): undecidable within the bound. Rerun: class (ii)/(iii)
records rerun through the unchanged drivers with a 17-digit writer; bit: objective, minimum spacing, evaluations,
rounded coordinates and convergence curve identical to the stored record; noise: identical except for float noise
($\le10^{-9}$ relative in the objective); near: rounded coordinates, label and objective ($\le10^{-9}$ relative)
reproduced but convergence curve or minimum spacing not (SLSQP-based methods; platform dependent), not used to
decide labels. Strict: label decided from the full-precision layout (bit-identical
rerun, or float-noise rerun with slack $>10^{-8}$~m), confirmed or changed. Open: class (ii)/(iii) labels not yet
decided (reruns pending, see \texttt{rev3\_precision\_plan.json}). Replay: maximum $|$objective recomputed from
the rounded coordinates minus stored objective$|$, in percentage points of the ideal (wake-free) value.
$^a$LX-SSA, SSA, PSO, DE, modified VNS, MS-SLSQP. $^s$Superseded, not used in the manuscript.
Replay in pp. HR1: Horns Rev~1. $^h$Contains Horns Rev~1 runs with the pre-2026-09-28 direction binning (replayed and rerun with that binning).}
\label{tab:S-r3-precision}
\scriptsize\setlength{\tabcolsep}{1.5pt}
\begin{tabular}{lrrrrrrrrrrrrr}
\toprule
& & & \multicolumn{3}{c}{Rounded coordinates} & \multicolumn{4}{c}{Reruns} & \multicolumn{2}{c}{Strict label} & & \\
\cmidrule(lr){4-6}\cmidrule(lr){7-10}\cmidrule(lr){11-12}
Study & Rec. & Feas. & (i) & (ii) & (iii) & Rerun & Bit & Noise & Near & Conf. & Chg. & Open & Replay \\
\midrule
""" + body + r"""
\midrule
Total & """ + " & ".join(f"{int(x):,}" for x in tot) + f" & {mx_all:.3f} \\\\" + r"""
\bottomrule
\end{tabular}
\end{table}
"""
    with open(os.path.join(HERE, "rev3_precision_tables.tex"), "w") as fh:
        fh.write(tex)
    print(tex)


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
        if len(a.arg) >= 3:
            sh, nsh = int(a.arg[1]), int(a.arg[2])
            pl = {d["study"]: d for d in json.load(open(os.path.join(HERE, "rev3_precision_plan.json")))["studies"]}
            shards = pl[st]["shards"]
            if len(shards) != nsh:
                raise SystemExit(f"plan has {len(shards)} shards for {st}, not {nsh}")
            keys = [tuple(k) for k in shards[sh]]
            run_keys(keys, f"{st}_s{sh}of{nsh}", a.procs)
            write_shard_csv(st, sh, nsh, keys)
        else:
            run_keys(keys, st, a.procs)
    elif a.cmd == "keys":
        run_keys([tuple(k) for k in json.load(open(a.arg[0]))], os.path.splitext(os.path.basename(a.arg[0]))[0], a.procs)
    elif a.cmd == "plan":
        plan(exclude=tuple(a.arg[0].split(",")) if a.arg else ())
    elif a.cmd == "collect":
        collect()
    elif a.cmd == "export-local":
        export_local()
    elif a.cmd == "table":
        table()
    else:
        raise SystemExit(__doc__)
