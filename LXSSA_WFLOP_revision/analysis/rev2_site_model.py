"""Second measured-wind site (R1.9, AE.8): a 16-turbine block of the Lillgrund offshore wind farm (Sweden).

Data (DTU PyWake 2.6.20, py_wake/examples/data/lillgrund.py, MIT licence), copied below as constants so that PyWake
is not needed at run time:
  * installed layout of the 48 Siemens SWT-2.3-93 turbines (UTM easting/northing, m), wt_x / wt_y;
  * SWT-2.3-93 power (kW) and thrust-coefficient curves at 3, 4, ..., 25 m/s (D = 93 m, hub height 65 m);
  * the 12-sector Weibull climate of LillgrundSite (sector frequencies normalised to one). Provenance, as stated in
    the PyWake source: "The Weibull parameters are based on 7-months of measurements between 06/2012 - 01/2013,
    Lillgrund 61m met-mast", with reference doi:10.1016/j.renene.2016.07.038 (T. Goecmen, G. Giebel, Estimation of
    turbulence intensity using rotor effective wind speed in Lillgrund and Horns Rev-I offshore wind farms,
    Renewable Energy 99 (2016) 524-532). It is therefore a sector-Weibull fit to met-mast measurements at the site,
    like the Horns Rev 1 climate -- but from only seven months (June-January, so seasonally unbalanced) and at 61 m
    rather than the 65-m hub height (no shear correction is applied, as in PyWake).
Model: identical to hornsrev_model (Jensen top-hat wake with k = 0.04, hub-centre in-wake test, root-sum-square
superposition, thrust coefficient at the free-stream speed, 5-degree direction bins centred at 2.5, 7.5, ...,
357.5 deg with exactly six bins per 30-degree sector, 1 m/s speed bins from 3 to 25 m/s, AEP in GWh/yr).

Block choice. The farm is a skewed lattice with basis vectors U = (-266, -299) m (|U| = 400 m = 4.3 D, the
within-row spacing) and W = (-266, 152) m (|W| = 306 m = 3.3 D, the between-row spacing) and has a gap of two
missing positions in its centre. Only two complete 4 x 4 blocks exist; we take the one at the eastern corner of the
farm (lattice indices 0..3 along U and W; PyWake turbine indices I16 below), the analogue of the corner block of
Horns Rev 1. The boundary is the parallelogram through its four corner turbines (0, 3, 26, 23), enlarged by 0.2 %
about the block centre (0.1 % is not enough: turbine 1 lies 0.69 m outside the exact corner parallelogram because the
as-built positions deviate from the lattice by a few decimetres); all 16 installed turbines lie >= 0.2 m inside.

Minimum spacing. The Horns Rev 1 block uses l_min = 4 D. That is IMPOSSIBLE here: the installed block has
3.3 D x 4.3 D spacing, and by Oler's inequality (N <= A / (sqrt(3)/2 d^2) + P / (2 d) + 1 for points at mutual
distance >= d in a convex region of area A and perimeter P) at most 15.7 turbines fit at d = 4 D = 372 m
(A = 1.081 km^2, P = 4.24 km); multistart max-min-distance optimisation finds at most 3.64 D. We therefore use
l_min = 3 D = 279 m (SMIN below), just under the installed minimum spacing (3.29 D), so that the installed block is
feasible and serves as the reference, as the installed Horns Rev 1 block does (Oler bound at 3 D: 24.6 turbines).
make_objective(n, smin=...) accepts another value.

Validation against PyWake (run  python3 rev2_site_model.py --pywake ; needs py_wake==2.6.20; writes
rev2_site_pywake_check.csv).
"""
import sys
import numpy as np

D = 93.0
RR = D / 2
KW = 0.04
HOURS = 8760.0
SMIN = 3 * D                      # see the module docstring (4 D is infeasible for 16 turbines in this block)
ENLARGE = 1.002

WS_TAB = np.arange(3.0, 26.0)
P_TAB = np.array([0.0, 65.0, 180.0, 352.0, 590.0, 906.0, 1308.0, 1767.0, 2085.0, 2234.0, 2283.0, 2296.0,
                  2299.0] + [2300.0] * 10)                                                    # kW
CT_TAB = np.array([0.0, 0.81, 0.84, 0.83, 0.85, 0.86, 0.87, 0.79, 0.67, 0.45, 0.34, 0.26, 0.21, 0.17, 0.14,
                   0.12, 0.10, 0.09, 0.07, 0.07, 0.06, 0.05, 0.05])

SEC_F = np.array([3.8, 4.5, 0.4, 2.8, 8.3, 7.5, 9.9, 14.8, 14.3, 17.0, 12.6, 4.1])
SEC_F = SEC_F / SEC_F.sum()
SEC_A = np.array([4.5, 4.7, 3.0, 7.2, 8.8, 8.2, 8.4, 9.5, 9.2, 9.9, 10.3, 6.7])
SEC_K = np.array([1.69, 1.78, 1.82, 1.70, 1.97, 2.49, 2.72, 2.70, 2.88, 3.34, 2.84, 2.23])

WT_X = np.array([361469, 361203, 360936, 360670, 360404, 360137, 359871, 361203, 360936, 360670, 360404, 360137,
                 359871, 359604, 359338, 360936, 360670, 360404, 360137, 359871, 359604, 359338, 359071, 360670,
                 360404, 360137, 359871, 359338, 359071, 358805, 360390, 360137, 359871, 359604, 359071, 358805,
                 359871, 359604, 359338, 359071, 358805, 359604, 359338, 359071, 358805, 359338, 359071, 358805], float)
WT_Y = np.array([6154543, 6154244, 6153946, 6153648, 6153349, 6153051, 6152753, 6154695, 6154396, 6154098, 6153800,
                 6153501, 6153203, 6152905, 6152606, 6154847, 6154548, 6154250, 6153952, 6153653, 6153355, 6153057,
                 6152758, 6154999, 6154701, 6154402, 6154104, 6153507, 6153209, 6152910, 6155136, 6154853, 6154554,
                 6154256, 6153659, 6153361, 6155005, 6154706, 6154408, 6154110, 6153811, 6155157, 6154858, 6154560,
                 6154262, 6155309, 6155010, 6154712], float)
# 4 x 4 block: lattice (i along U, j along W), ordered i-major: (0,0),(0,1),(0,2),(0,3),(1,0),...,(3,3)
I16 = [0, 7, 15, 23, 1, 8, 16, 24, 2, 9, 17, 25, 3, 10, 18, 26]
CORNERS = [0, 3, 15, 12]          # positions in I16 of the lattice corners (0,0), (0,3), (3,3), (3,0)

WD = np.arange(2.5, 360.0, 5.0)                        # meteorological "from" direction (bin centres)
_SEC = (np.floor((WD + 15.0) / 30.0).astype(int)) % 12  # exactly 6 bins per sector
F_WD = SEC_F[_SEC] / 6.0
A_WD, K_WD = SEC_A[_SEC], SEC_K[_SEC]
TOWARD = np.deg2rad(270.0 - WD)
WS = np.arange(3.0, 26.0)
_lo, _hi = WS - 0.5, WS + 0.5
P_WS = (np.exp(-(_lo[None] / A_WD[:, None]) ** K_WD[:, None])
        - np.exp(-(_hi[None] / A_WD[:, None]) ** K_WD[:, None]))       # (n_wd, n_ws)
A_CT = 1 - np.sqrt(1 - np.interp(WS, WS_TAB, CT_TAB))                  # (n_ws,)


def site(n_turbines=16):
    """Installed 16-turbine block (local coordinates, centred on the block mean) and its outline (counter-clockwise)."""
    if n_turbines != 16:
        raise ValueError("only the 16-turbine block is defined")
    xy = np.c_[WT_X[I16], WT_Y[I16]]
    xy = xy - xy.mean(0)
    poly = xy[CORNERS]
    area = 0.5 * np.sum(poly[:, 0] * np.roll(poly[:, 1], -1) - np.roll(poly[:, 0], -1) * poly[:, 1])
    if area < 0:
        poly = poly[::-1]
    return xy, poly * ENLARGE


def outside_distance(xy, poly):
    """Distance (m) of each point outside the convex polygon (0 inside)."""
    a, b = poly, np.roll(poly, -1, axis=0)
    e = b - a
    nrm = np.c_[e[:, 1], -e[:, 0]] / np.linalg.norm(e, axis=1)[:, None]
    s = np.einsum("ijk,jk->ij", xy[:, None, :] - a[None], nrm)
    return np.maximum(s.max(1), 0.0)


def aep_gwh(xy, with_wake=True, wd=None, f_wd=None, p_ws=None, local_ct=False):
    """Annual energy production (GWh/yr) of layout xy (N x 2, m). wd/f_wd/p_ws override the paper's 5-degree bins
    (used only by the PyWake check); local_ct=True takes C_T at the waked speed of the upstream turbine (diagnostic,
    processing turbines from upstream to downstream)."""
    xy = np.asarray(xy, float)
    n = len(xy)
    wd = WD if wd is None else wd
    f_wd = F_WD if f_wd is None else f_wd
    p_ws = P_WS if p_ws is None else p_ws
    if with_wake:
        tow = np.deg2rad(270.0 - wd)
        dx = xy[:, 0][:, None] - xy[:, 0][None, :]
        dy = xy[:, 1][:, None] - xy[:, 1][None, :]
        c = np.cos(tow)[:, None, None]; s = np.sin(tow)[:, None, None]
        x = dx[None] * c + dy[None] * s               # downstream distance of i behind j
        lat = np.abs(-dx[None] * s + dy[None] * c)
        inw = (x > 0) & (lat < RR + KW * x)
        g = np.where(inw, (RR / (RR + KW * np.maximum(x, 0))) ** 2, 0.0)   # (n_wd, i, j)
        if not local_ct:
            G = np.sqrt((g ** 2).sum(2))                                   # (n_wd, n)
            u = WS[None, :, None] * (1 - A_CT[None, :, None] * G[:, None, :])
        else:
            # as mpce_direction.hr_aep(local_ct=True), which reproduces PyWake NOJ + RotorCenter + ct2a_mom1d
            u = np.empty((len(wd), len(WS), n))
            for d in range(len(wd)):
                proj = xy @ np.array([np.cos(tow[d]), np.sin(tow[d])])
                a_up = np.zeros((len(WS), n))
                for i in np.argsort(proj, kind="stable"):
                    ui = WS * (1 - np.sqrt(((a_up * g[d, i][None, :]) ** 2).sum(1)))
                    u[d, :, i] = ui
                    a_up[:, i] = 1 - np.sqrt(1 - np.interp(ui, WS_TAB, CT_TAB))
    else:
        u = np.broadcast_to(WS[None, :, None], (len(wd), len(WS), n))
    p = np.interp(u, WS_TAB, P_TAB)
    e = np.einsum("dsn,ds,d->", p, p_ws, f_wd)
    return e * HOURS / 1e6


def make_objective(n_turbines=16, smin=SMIN, penalty=1e10):
    """Minimisation objective in the authors' form: wake loss (GWh) + quadratic penalties (as hornsrev_model)."""
    _, poly = site(n_turbines)
    ideal = aep_gwh(np.zeros((n_turbines, 2)), with_wake=False)

    def f(x):
        xy = x.reshape(-1, 2)
        loss = ideal - aep_gwh(xy)
        g = outside_distance(xy, poly)
        iu = np.triu_indices(len(xy), 1)
        d = np.sqrt(((xy[:, None] - xy[None]) ** 2).sum(-1))[iu]
        gs = smin - d
        return loss + ((1 + penalty * g[g > 0]) ** 2).sum() + ((1 + penalty * gs[gs > 0]) ** 2).sum()
    return f, ideal, poly


def bins(step, first):
    """Direction bins of width step centred at first, first + step, ... (sector by nearest centre, frequency shared
    equally by the bins of a sector, as PyWake's 'nearest' sector lookup times the bin width)."""
    wd = np.arange(first, 360.0, step)
    sec = (np.floor((wd + 15.0) / 30.0).astype(int)) % 12
    cnt = np.bincount(sec, minlength=12)
    f = SEC_F[sec] / cnt[sec]
    pws = (np.exp(-(_lo[None] / SEC_A[sec][:, None]) ** SEC_K[sec][:, None])
           - np.exp(-(_hi[None] / SEC_A[sec][:, None]) ** SEC_K[sec][:, None]))
    return wd, f, pws


def pywake_check(out="rev2_site_pywake_check.csv"):
    import time, pandas as pd
    import py_wake
    from py_wake.examples.data import lillgrund as L
    from py_wake import NOJ
    from py_wake.rotor_avg_models import RotorCenter
    from py_wake.deficit_models.utils import ct2a_mom1d
    assert np.allclose(L.wt_x, WT_X) and np.allclose(L.wt_y, WT_Y)
    assert np.allclose(L.power_curve[:, 1], P_TAB) and np.allclose(L.ct_curve[:, 1], CT_TAB)
    s0, wt = L.LillgrundSite(), L.SWT23()
    assert abs(wt.diameter() - D) < 1e-9
    models = {"NOJ_k0.04": NOJ(s0, wt, k=0.04),
              "NOJ_k0.04_RotorCenter": NOJ(s0, wt, k=0.04, rotorAvgModel=RotorCenter()),
              "NOJ_k0.04_RotorCenter_momentum": NOJ(s0, wt, k=0.04, rotorAvgModel=RotorCenter(), ct2a=ct2a_mom1d)}
    farms = {"LG16": (WT_X[I16], WT_Y[I16]), "LG48": (WT_X, WT_Y)}
    rows = []
    for bname, (step, first) in {"ours_5deg_2.5": (5.0, 2.5), "1deg_0.5": (1.0, 0.5)}.items():
        wd, f, pws = bins(step, first)
        for farm, (x, y) in farms.items():
            xy = np.c_[x, y] - np.c_[x, y].mean(0)
            ideal = aep_gwh(xy, False, wd, f, pws)
            for mname, lc in (("paper_model", False), ("paper_model_localCT", True)):
                a = aep_gwh(xy, True, wd, f, pws, local_ct=lc)
                rows.append(dict(Farm=farm, Bins=bname, Model=mname, Source="ours", AEP_GWh=a, IdealAEP_GWh=ideal,
                                 WakeLossPct=100 * (1 - a / ideal), PMaxAbsDiff=np.nan))
            for mname, m in models.items():
                t = time.time()
                sim = m(x, y, wd=wd, ws=WS)
                a, a0 = float(sim.aep().sum()), float(sim.aep(with_wake_loss=False).sum())
                P = np.asarray(sim.P.transpose("wd", "ws").values)
                rows.append(dict(Farm=farm, Bins=bname, Model=mname, Source=f"PyWake {py_wake.__version__}",
                                 AEP_GWh=a, IdealAEP_GWh=a0, WakeLossPct=100 * (1 - a / a0),
                                 PMaxAbsDiff=float(np.abs(P - f[:, None] * pws).max())))
    df = pd.DataFrame(rows)
    ref = df[df.Model == "paper_model"].set_index(["Farm", "Bins"]).AEP_GWh
    df["OursMinusThisPct"] = [100 * (ref[(r.Farm, r.Bins)] / r.AEP_GWh - 1) for r in df.itertuples()]
    df.to_csv(out, index=False, float_format="%.8g")
    with __import__("pandas").option_context("display.width", 200, "display.max_columns", 20):
        print(df)
    return df


if __name__ == "__main__":
    xy, poly = site(16)
    f, ideal, _ = make_objective(16)
    a = aep_gwh(xy)
    d = np.sqrt(((xy[:, None] - xy[None]) ** 2).sum(-1))[np.triu_indices(16, 1)]
    print(f"installed block: AEP {a:.4f} GWh/yr, wake-free {ideal:.4f}, wake loss {ideal - a:.4f} GWh/yr "
          f"({100 * (1 - a / ideal):.3f} %), min spacing {d.min():.1f} m = {d.min() / D:.3f} D, "
          f"max outside {outside_distance(xy, poly).max():.3g} m, F_p {f(xy.ravel()):.4f}, "
          f"half-width {np.abs(poly).max():.1f} m, sum F_WD {F_WD.sum():.12f}")
    if "--pywake" in sys.argv:
        pywake_check()
