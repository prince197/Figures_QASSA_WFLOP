"""Horns Rev 1 wake model (hornsrev_model.py) against DTU PyWake at identical direction and speed bins (R4 issue 2).

OPTIONAL dependency (not in requirements.txt):  python3 -m pip install py_wake==2.6.20
(it pulls xarray, netCDF4, h5netcdf, autograd, numpy_financial, tqdm, pooch, joblib). Run:
    python3 analysis/pywake_check.py [--layouts]
The results are stored in analysis/pywake_check.csv (and, with --layouts, analysis/pywake_check_hr16runs.csv/.json);
mpce_direction.py reads them, so the rest of the pipeline does not need PyWake.

What is compared (installed 80-turbine farm and the 16-turbine north-west block, i.e. PyWake's wt16_x/wt16_y = our
hornsrev_model.I16):
  * bins "ours_5deg_2.5": the paper's 72 direction bins of 5 deg centred at 2.5, 7.5, ... and speed bins 3..25 m/s;
    "1deg_0.5": 360 bins of 1 deg centred at 0.5, 1.5, ...; "pywake_default_1deg_0": PyWake's default (0, 1, ..., 359).
    PyWake's Hornsrev1Site (UniformWeibullSite, 'nearest' sector interpolation) is evaluated at exactly these wd/ws
    values; the column PMaxAbsDiff gives max |P_PyWake(wd, ws) - P_ours(wd, ws)| (identical bins and probabilities).
  * models: PyWake NOJ(site, V80(), k=0.04) with its defaults (rotor-area overlap AreaOverlapAvgModel, Madsen a(C_T),
    C_T at the effective speed of the wake-generating turbine, SquaredSum superposition, PropagateDownwind);
    NOJ with RotorCenter() (hub-centre test as ours); NOJ with RotorCenter() and the 1-D momentum a(C_T)
    (ct2a_mom1d; our 1 - sqrt(1 - C_T) = 2a); our model (mpce_direction.hr_aep = hornsrev_model.aep_gwh at the paper's
    bins) and a DIAGNOSTIC variant of our model with C_T at the local (waked) speed (mpce_direction.hr_aep(local_ct=True)).
  * --layouts: every run of experiment hrfix (16 turbines) and the installed block with PyWake NOJ at the 1-deg bins
    "1deg_0.5" (coordinates shifted by the centre of the installed block, as hornsrev_model.site()).
"""
import os, sys, glob, json, time, argparse
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import hornsrev_model as H
import mpce_direction as D

BINS = {"ours_5deg_2.5": (5.0, 2.5), "1deg_0.5": (1.0, 0.5), "pywake_default_1deg_0": (1.0, 0.0)}


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--layouts", action="store_true", help="also re-evaluate the 970 hrfix runs with PyWake NOJ (1 deg)")
    ap.add_argument("--out-dir", default=HERE)
    args = ap.parse_args(argv)
    try:
        import py_wake
        from py_wake.examples.data.hornsrev1 import Hornsrev1Site, V80, wt_x, wt_y, wt16_x, wt16_y
        from py_wake import NOJ
        from py_wake.rotor_avg_models import RotorCenter
        from py_wake.deficit_models.utils import ct2a_mom1d
    except ImportError as e:
        print(f"PyWake not available ({e}); install with  python3 -m pip install py_wake==2.6.20  -- no output written.")
        sys.exit(2)
    ver = py_wake.__version__
    site, wt = Hornsrev1Site(), V80()
    models = {"NOJ_k0.04": NOJ(site, wt, k=0.04),
              "NOJ_k0.04_RotorCenter": NOJ(site, wt, k=0.04, rotorAvgModel=RotorCenter()),
              "NOJ_k0.04_RotorCenter_momentum": NOJ(site, wt, k=0.04, rotorAvgModel=RotorCenter(), ct2a=ct2a_mom1d)}
    settings = {"NOJ_k0.04": "NOJ(Hornsrev1Site(), V80(), k=0.04): AreaOverlapAvgModel, ct2a_madsen, SquaredSum, PropagateDownwind",
                "NOJ_k0.04_RotorCenter": "as NOJ_k0.04 with rotorAvgModel=RotorCenter()",
                "NOJ_k0.04_RotorCenter_momentum": "as NOJ_k0.04 with rotorAvgModel=RotorCenter(), ct2a=ct2a_mom1d",
                "paper_model": "hornsrev_model: Jensen top-hat k=0.04, hub-centre test, RSS, C_T at free-stream speed",
                "paper_model_localCT": "hornsrev_model equations with C_T at the local (waked) speed of the upstream turbine (diagnostic)"}
    farms = {"HR80": (np.asarray(wt_x, float), np.asarray(wt_y, float)), "HR16": (np.asarray(wt16_x, float), np.asarray(wt16_y, float))}
    # consistency of the layouts with hornsrev_model
    assert np.allclose(farms["HR80"][0], H.WT_X) and np.allclose(farms["HR80"][1], H.WT_Y)
    assert np.allclose(farms["HR16"][0], H.WT_X[H.I16]) and np.allclose(farms["HR16"][1], H.WT_Y[H.I16])
    rows = []
    for bname, (step, off) in BINS.items():
        b = D.hr_bins(step, off)
        P_ours = b["f"][:, None] * b["pws"]
        for farm, (x, y) in farms.items():
            xy = np.c_[x, y] - np.c_[x, y].mean(0)
            ideal = D.hr_aep(xy, b, with_wake=False)
            for mname, lc in (("paper_model", False), ("paper_model_localCT", True)):
                a = D.hr_aep(xy, b, local_ct=lc)
                rows.append(dict(Farm=farm, Bins=bname, Model=mname, Source="ours", Version="--", Settings=settings[mname],
                                 WdStep=step, WdFirst=off, NWd=len(b["wd"]), AEP_GWh=a, IdealAEP_GWh=ideal,
                                 WakeLossPct=100 * (1 - a / ideal), PMaxAbsDiff=np.nan, Seconds=np.nan))
            for mname, mdl in models.items():
                t = time.time()
                sim = mdl(x, y, wd=b["wd"], ws=H.WS)
                a = float(sim.aep().sum())
                a0 = float(sim.aep(with_wake_loss=False).sum())
                P = np.asarray(sim.P.transpose("wd", "ws").values if hasattr(sim.P, "transpose") else sim.P)
                if P.ndim == 3:             # (i, wd, ws) for some site types
                    P = P[0]
                rows.append(dict(Farm=farm, Bins=bname, Model=mname, Source="PyWake", Version=ver, Settings=settings[mname],
                                 WdStep=step, WdFirst=off, NWd=len(b["wd"]), AEP_GWh=a, IdealAEP_GWh=a0,
                                 WakeLossPct=100 * (1 - a / a0), PMaxAbsDiff=float(np.abs(P - P_ours).max()),
                                 Seconds=time.time() - t))
                print(f"{farm} {bname:22s} {mname:32s} AEP {a:9.3f}  ideal {a0:9.3f}  loss {100 * (1 - a / a0):6.3f} %  "
                      f"|dP| {rows[-1]['PMaxAbsDiff']:.1e}", flush=True)
            for r in rows[-5:-3]:
                print(f"{farm} {bname:22s} {r['Model']:32s} AEP {r['AEP_GWh']:9.3f}  ideal {r['IdealAEP_GWh']:9.3f}  loss {r['WakeLossPct']:6.3f} %")
    out = pd.DataFrame(rows)
    out.insert(0, "Generated", time.strftime("%Y-%m-%d %H:%M:%S"))
    out.to_csv(os.path.join(args.out_dir, "pywake_check.csv"), index=False, float_format="%.10g")
    print(f"wrote pywake_check.csv ({len(out)} rows), PyWake {ver}")

    if args.layouts:
        files = sorted(glob.glob(os.path.join(HERE, "mpce_hrfix_s*of24.csv")))
        df = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)
        df = df[df.Dataset.astype(str) == "HR"].drop_duplicates(["Algorithm", "Seed", "Budget", "Init"], keep="last")
        x16, y16 = farms["HR16"]
        cx, cy = x16.mean(), y16.mean()
        b = D.hr_bins(1.0, 0.5)
        m = models["NOJ_k0.04"]
        inst = float(m(x16, y16, wd=b["wd"], ws=H.WS).aep().sum())
        t = time.time()
        vals = []
        for i, c in enumerate(df.Coordinates):
            xy = np.array([[float(v) for v in p.split()] for p in c.split(";")])
            vals.append(float(m(xy[:, 0] + cx, xy[:, 1] + cy, wd=b["wd"], ws=H.WS).aep().sum()))
            if i % 100 == 0:
                print(f"  layouts {i}/{len(df)} {time.time() - t:.0f} s", flush=True)
        df["PyWakeNOJ_1deg"] = vals
        df[["Algorithm", "Seed", "Budget", "Init", "Feasible", "Objective", "PyWakeNOJ_1deg"]].to_csv(
            os.path.join(args.out_dir, "pywake_check_hr16runs.csv"), index=False, float_format="%.6f")
        json.dump(dict(pywake_version=ver, model=settings["NOJ_k0.04"], bins="1 deg centred at 0.5, 1.5, ...; ws 3..25 m/s",
                       installed_aep_gwh=inst, n_runs=int(len(df)), seconds=time.time() - t,
                       coordinates="run coordinates (local) + centre of the installed 16-turbine block"),
                  open(os.path.join(args.out_dir, "pywake_check_hr16runs.json"), "w"), indent=1)
        print(f"wrote pywake_check_hr16runs.csv ({len(df)} runs), installed {inst:.3f} GWh/yr")


if __name__ == "__main__":
    main()
