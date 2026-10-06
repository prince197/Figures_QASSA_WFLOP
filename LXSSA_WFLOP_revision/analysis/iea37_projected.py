"""IEA37 Case Study 1: published layouts that violate the boundary slightly, projected onto it.

Usage:  python3 analysis/iea37_projected.py      (also run from mpce_direction.py, which writes the \\NFIEA... macros)

The published comparison of the paper counts a published layout as feasible if no turbine lies more than 1 mm
outside the boundary circle and no pair is closer than 2D = 260 m (column Feasible_tol1e-3m of
iea37_published_results.csv; "strict" convention). Several submissions place turbines a few millimetres (or, for
participant 12 with 16 turbines, 3.5 m) outside the circle. Here every turbine outside the circle is moved
radially onto it (x <- x R / |x|), the spacing is re-checked (>= 260 m - 1e-6 m) and the AEP of the projected layout
is recomputed with iea37_model.aep (which reproduces the official calculator). A projected layout is a feasible layout
in its own right, so the "projected" convention uses it as a published reference whenever it is feasible (layouts
that violate the spacing cannot be repaired this way and stay infeasible).

The gaps of our layouts (PSO-VNS: best and mean of 30 runs; best run of the eight main methods) to the best feasible
published layout are then recomputed under both conventions, in % of AEP and in percentage points of wake loss
(wake loss = 1 - AEP / wake-free AEP; the wake-free AEP is every turbine at rated power in every direction).

Outputs (analysis/): iea37_projected.csv (one row per published layout), iea37_projected.json.
"""
import os, sys, json
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import iea37_model as IM

PARTICIPANTS = [f"par{i}" for i in range(1, 13)]
STRICT_TOL = 1e-3          # m, published-layout tolerance of the paper (Feasible_tol1e-3m)
OUR_TOL = 1e-6             # m, tolerance of our own layouts / of the projected layouts
MAIN8 = ["PSOBV", "PSOC", "SSABV", "SSA", "LXSSA", "DE", "BVNS", "SLSQP"]
FOCUS = "PSOBV"


def project(xy, radius):
    r = np.sqrt((xy ** 2).sum(1))
    s = np.where(r > radius, radius / r, 1.0)
    return xy * s[:, None]


def published(log=print):
    rows = []
    for n in (16, 36):
        R = IM.RADIUS[n]
        ideal = IM.ideal_aep(n)
        for p in PARTICIPANTS:
            fn = os.path.join(IM.DATA, "iea37-cs1-results", f"iea37-par{p[3:]}-opt{n}.yaml")
            xy, rep = IM.load_submission(p[3:], n)
            a = IM.aep(xy)
            exc = float(np.sqrt((xy ** 2).sum(1)).max() - R)
            sp = float(IM.min_spacing(xy))
            strict = bool(exc <= STRICT_TOL and sp >= IM.SMIN - STRICT_TOL)
            xyp = project(xy, R)
            spp = float(IM.min_spacing(xyp))
            ap = float(IM.aep(xyp))
            featp = bool(spp >= IM.SMIN - OUR_TOL and np.sqrt((xyp ** 2).sum(1)).max() <= R + OUR_TOL)
            rows.append(dict(turbines=n, participant=p, file=os.path.relpath(fn, HERE), aep_published=float(rep), aep_ours=float(a),
                             max_excess_m=exc, n_outside_strict=int((np.sqrt((xy ** 2).sum(1)) > R + STRICT_TOL).sum()),
                             min_spacing_m=sp, feasible_strict=strict, projected_ok=True, aep_projected=ap,
                             min_spacing_projected_m=spp, feasible_projected=featp,
                             loss_pct=100 * (1 - a / ideal), loss_projected_pct=100 * (1 - ap / ideal),
                             aep_change_projection=ap - a))
    P = pd.DataFrame(rows)
    return P


def ref(P, n, conv):
    """Best feasible published layout of scenario n under a convention ('strict' or 'projected')."""
    Q = P[P.turbines == n]
    if conv == "strict":
        Q = Q[Q.feasible_strict]
        a = Q.aep_ours
    else:
        Q = Q[Q.feasible_strict | Q.feasible_projected]
        a = np.where(Q.feasible_strict, Q.aep_ours, Q.aep_projected)
    Q = Q.assign(A=a)
    b = Q.sort_values("A").iloc[-1]
    return dict(best_by=b.participant, best_aep=float(b.A), best_loss_pct=100 * (1 - float(b.A) / IM.ideal_aep(n)),
                n_feasible=int(len(Q)), aeps=sorted(map(float, Q.A)),
                best_max_excess_m=float(b.max_excess_m), best_min_spacing_m=float(b.min_spacing_projected_m if conv == "projected" else b.min_spacing_m))


def run(M=None, ALL=None, log=print, out_dir=HERE):
    P = published(log)
    # consistency with the published-results file used by mpce_results.py
    pub = pd.read_csv(os.path.join(HERE, "iea37_published_results.csv"))
    pub = pub[pub.Participant.str.startswith("par")]
    mm = P.merge(pub, left_on=["turbines", "participant"], right_on=["Turbines", "Participant"])
    chk = dict(n_matched=int(len(mm)),
               aep_matches_file=bool(np.allclose(mm.aep_ours, mm.AEP_MWh_iea37_model, rtol=0, atol=1e-3)),
               projected_aep_matches_file=bool(np.allclose(mm.aep_projected, mm.AEP_MWh_boundary_projected, rtol=0, atol=1e-2)),
               strict_matches_file=bool((mm.feasible_strict == mm["Feasible_tol1e-3m"]).all()))
    changed = P[(~P.feasible_strict) & P.feasible_projected]
    still = P[(~P.feasible_strict) & (~P.feasible_projected)]
    out = dict(rule=dict(strict_tol_m=STRICT_TOL, projected_tol_m=OUR_TOL, projection="radial, only turbines outside the circle",
                         spacing_m=IM.SMIN),
               consistency=chk,
               participants=P.to_dict(orient="records"),
               changed_status=[dict(participant=r.participant, turbines=int(r.turbines), max_excess_m=float(r.max_excess_m),
                                    aep_published=float(r.aep_ours), aep_projected=float(r.aep_projected),
                                    min_spacing_projected_m=float(r.min_spacing_projected_m)) for r in changed.itertuples()],
               still_infeasible=[dict(participant=r.participant, turbines=int(r.turbines), min_spacing_m=float(r.min_spacing_m),
                                      max_excess_m=float(r.max_excess_m)) for r in still.itertuples()],
               scenarios={})
    if ALL is None:
        import mpce_results as M
        ALL, _, _, _ = M.load(HERE, False)
    IE = ALL[ALL.Dataset.astype(str).str.startswith("IEA37") & (ALL.Init == "random")]
    for n in (16, 36):
        ideal = IM.ideal_aep(n)
        sc = dict(strict=ref(P, n, "strict"), projected=ref(P, n, "projected"), ideal_aep=ideal, ours={})
        for b in (6030, 30030):
            Y = IE[(IE.Turbines == n) & (IE.Budget == b) & IE.Feasible]
            f = Y[Y.Algorithm == FOCUS].Objective
            if not len(f):
                continue
            ob = Y[Y.Algorithm.isin(MAIN8)].sort_values("Objective").iloc[-1]
            best, mean = float(f.max()), float(f.mean())
            g = dict(psovns_best=best, psovns_mean=mean, psovns_feasible=int(len(f)),
                     psovns_best_loss_pct=100 * (1 - best / ideal), psovns_mean_loss_pct=100 * (1 - mean / ideal),
                     our_best_method=str(ob.Algorithm), our_best_aep=float(ob.Objective))
            for conv in ("strict", "projected"):
                r = sc[conv]
                g[conv] = dict(ref_by=r["best_by"], ref_aep=r["best_aep"], ref_loss_pct=r["best_loss_pct"],
                               n_feasible_published=r["n_feasible"],
                               best_gap_pct=100 * (best / r["best_aep"] - 1), mean_gap_pct=100 * (mean / r["best_aep"] - 1),
                               best_ratio_pct=100 * best / r["best_aep"], mean_ratio_pct=100 * mean / r["best_aep"],
                               best_gap_loss_pp=g["psovns_best_loss_pct"] - r["best_loss_pct"],
                               mean_gap_loss_pp=g["psovns_mean_loss_pct"] - r["best_loss_pct"],
                               best_rank=int(1 + sum(p > best for p in r["aeps"])), mean_rank=int(1 + sum(p > mean for p in r["aeps"])),
                               our_best_gap_pct=100 * (g["our_best_aep"] / r["best_aep"] - 1),
                               our_best_gap_loss_pp=100 * (1 - g["our_best_aep"] / ideal) - r["best_loss_pct"],
                               our_best_rank=int(1 + sum(p > g["our_best_aep"] for p in r["aeps"])))
            sc["ours"][str(b)] = g
            log(f"  IEA37 {n}T {b}: PSO-VNS best {best:,.1f} (loss {g['psovns_best_loss_pct']:.2f} %); strict ref "
                f"{g['strict']['ref_by']} {g['strict']['ref_aep']:,.1f} gap {g['strict']['best_gap_pct']:+.2f} % / "
                f"{g['strict']['best_gap_loss_pp']:+.2f} pp rank {g['strict']['best_rank']}; projected ref {g['projected']['ref_by']} "
                f"{g['projected']['ref_aep']:,.1f} gap {g['projected']['best_gap_pct']:+.2f} % / {g['projected']['best_gap_loss_pp']:+.2f} pp "
                f"rank {g['projected']['best_rank']}; best of ours {g['our_best_method']} {g['projected']['our_best_gap_pct']:+.2f} %")
        out["scenarios"][str(n)] = sc
    log("  changed status (feasible only after projection): " +
        ", ".join(f"{c['participant']}-{c['turbines']} ({c['max_excess_m']:.4f} m)" for c in out["changed_status"]))
    log("  still infeasible: " + ", ".join(f"{c['participant']}-{c['turbines']} (spacing {c['min_spacing_m']:.1f} m)" for c in out["still_infeasible"]))
    log(f"  consistency with iea37_published_results.csv: {chk}")
    P.to_csv(os.path.join(out_dir, "iea37_projected.csv"), index=False, float_format="%.6f")

    def cl(o):
        if isinstance(o, dict):
            return {k: cl(v) for k, v in o.items()}
        if isinstance(o, list):
            return [cl(v) for v in o]
        if isinstance(o, (np.integer,)):
            return int(o)
        if isinstance(o, (np.floating,)):
            return float(o)
        if isinstance(o, np.bool_):
            return bool(o)
        return o
    out = cl(out)
    json.dump(out, open(os.path.join(out_dir, "iea37_projected.json"), "w"), indent=1)
    return out


if __name__ == "__main__":
    run()
