"""IEA Wind Task 37 WFLO Case Study 1 ("optimization only") model, 16- and 36-turbine scenarios.

Source: https://github.com/IEAWindTask37/iea37-wflo-casestudies (mirror of
https://github.com/byuflowlab/iea37-wflo-casestudies), folder cs1-2, commit 267f6e5 (2020-11-13).
The official files used here are copied, unmodified, into iea37_data/.

Definitions (iea37-wflocs-announcement.pdf, iea37-wakemodel.pdf, iea37-335mw.yaml, iea37-windrose.yaml):
  * turbine: IEA37 3.35 MW onshore reference turbine, D = 130 m, hub height 110 m, cut-in 4 m/s,
    rated 9.8 m/s, cut-out 25 m/s, P = Prated ((V - Vci)/(Vr - Vci))^3 for Vci <= V < Vr,
    Prated for Vr <= V < Vco, 0 otherwise;
  * wind: 16 direction bins (0, 22.5, ..., 337.5 deg; meteorological, 0 = North = +y, clockwise),
    constant free-stream speed 9.8 m/s, TI 0.075;
  * wake: simplified Bastankhah Gaussian, CT = 8/9, ky = 0.0324555, sigma = ky dx + D/sqrt(8),
    root-sum-square superposition;
  * boundary: circle centred at (0, 0), radius 1300 m (16 turbines), 2000 m (36), 3000 m (64);
  * spacing: no turbine closer than two rotor diameters (260 m) to any other.
Note: in the official naming, "Case Study 1" is the optimization-only study with the three farm
sizes 16/36/64; "Case Study 2" is the combined (own wake model) 9-turbine study.

AEP is in MWh and is computed exactly as in iea37-aepcalc.py (vectorised over directions and pairs).
"""
import os
import math
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "iea37_data")

D = 130.0                       # rotor diameter (m)
HUB_HEIGHT = 110.0              # m (not used by the wake model)
CUT_IN, RATED_WS, CUT_OUT = 4.0, 9.8, 25.0
RATED_PWR = 3350000.0           # W
WIND_SPEED = 9.8                # m/s, constant
WIND_DIR = np.arange(16) * 22.5
WIND_FREQ = np.array([.025, .024, .029, .036, .063, .065, .100, .122,
                      .063, .038, .039, .083, .213, .046, .032, .022])
CT = 4.0 * 1. / 3. * (1.0 - 1. / 3.)
K = 0.0324555
HRS = 365. * 24.
RADIUS = {16: 1300.0, 36: 2000.0, 64: 3000.0}
SMIN = 2 * D
BASELINE_AEP = {16: 366941.57116, 36: 737883.09851}   # 'default' AEP in iea37-ex16.yaml / iea37-ex36.yaml

# Rotation constants computed exactly as in WindFrame() of iea37-aepcalc.py (math.radians, cos/sin(-rad))
_rad = np.array([math.radians(270. - d) for d in WIND_DIR])
_COS = np.cos(-_rad)
_SIN = np.sin(-_rad)
_PAIRS = {}


def dir_power(xy):
    """Farm power (W) for each of the 16 directions; xy is (n, 2) in metres.

    Same equations as GaussianWake()/DirPower() of iea37-aepcalc.py, evaluated once per unordered
    pair: for a pair (i, j), the wake acts on whichever turbine is downwind (x_i - x_j > 0 means i
    is downwind of j; |x_i - x_j| is bit-identical to the official x), none if they are abeam."""
    xy = np.asarray(xy, float).reshape(-1, 2)
    n = len(xy)
    if n not in _PAIRS:
        iu, ju = np.triu_indices(n, 1)
        _PAIRS[n] = (iu, ju, np.arange(len(WIND_DIR))[:, None] * n)
    iu, ju, off = _PAIRS[n]
    x, y = xy[:, 0], xy[:, 1]
    fx = x[None, :] * _COS[:, None] - y[None, :] * _SIN[:, None]      # (16, n) downwind coordinate
    fy = x[None, :] * _SIN[:, None] + y[None, :] * _COS[:, None]      # (16, n) crosswind coordinate
    dx = fx[:, iu] - fx[:, ju]                                         # (16, pairs)
    dy = fy[:, iu] - fy[:, ju]
    adx = np.abs(dx)
    sigma = K * adx + D / np.sqrt(8.)
    d2 = ((1. - np.sqrt(1. - CT / (8. * sigma ** 2 / D ** 2))) * np.exp(-0.5 * (dy / sigma) ** 2)) ** 2
    d2[adx == 0.] = 0.0                                                # abeam: no wake either way
    rec = np.where(dx > 0., iu, ju) + off                              # index of the waked turbine
    loss = np.sqrt(np.bincount(rec.ravel(), d2.ravel(), minlength=len(WIND_DIR) * n)).reshape(-1, n)
    v = WIND_SPEED * (1. - loss)                                       # root-sum-square superposition
    p = np.where((CUT_IN <= v) & (v < RATED_WS), RATED_PWR * ((v - CUT_IN) / (RATED_WS - CUT_IN)) ** 3, 0.0)
    p = np.where((RATED_WS <= v) & (v < CUT_OUT), RATED_PWR, p)
    return p.sum(-1)


def aep_binned(xy):
    """AEP (MWh) per direction bin, as printed by iea37-aepcalc.py."""
    return HRS * (WIND_FREQ * dir_power(np.asarray(xy, float).reshape(-1, 2))) / 1.E6


def aep(xy):
    """Farm AEP (MWh)."""
    return aep_binned(xy).sum()


def ideal_aep(n):
    """Wake-free AEP (MWh): every turbine at rated power in every direction (9.8 m/s = rated)."""
    return (HRS * (WIND_FREQ * n * RATED_PWR) / 1.E6).sum()


def min_spacing(xy):
    xy = np.asarray(xy, float).reshape(-1, 2)
    d = np.sqrt(((xy[:, None] - xy[None]) ** 2).sum(-1))
    return d[np.triu_indices(len(xy), 1)].min()


def is_feasible(xy, n=None, tol=1e-6):
    xy = np.asarray(xy, float).reshape(-1, 2)
    r = RADIUS[len(xy) if n is None else n]
    return min_spacing(xy) >= SMIN - tol and np.sqrt((xy ** 2).sum(1)).max() <= r + tol


def make_objective(case, penalty=1e10):
    """case = 16 or 36 (or 64). Returns (f, ideal_aep, radius, smin).

    f(x) = (ideal - AEP) [MWh] + sum over turbines outside the circle of (1 + 1e10 g)^2 with
    g = x^2 + y^2 - r^2 (m^2) + sum over pairs closer than smin of (1 + 1e10 (smin - d))^2
    (the authors' penalty convention, authors_objective.py)."""
    n = int(case)
    radius, smin, ideal = RADIUS[n], SMIN, ideal_aep(n)
    iu = np.triu_indices(n, 1)

    def f(x):
        xy = np.asarray(x, float).reshape(-1, 2)
        loss = ideal - aep(xy)
        gb = (xy ** 2).sum(1) - radius ** 2
        gs = smin - np.sqrt(((xy[:, None] - xy[None]) ** 2).sum(-1))[iu]
        return loss + ((1 + penalty * gb[gb > 0]) ** 2).sum() + ((1 + penalty * gs[gs > 0]) ** 2).sum()
    return f, ideal, radius, smin


def _yaml(fn):
    import yaml
    with open(fn) as fh:
        return yaml.safe_load(fh)["definitions"]


def load_layout(fn):
    """(xy (n, 2), reported AEP MWh, reported binned AEP) from an official-format .yaml file."""
    d = _yaml(fn)
    xy = np.c_[d["position"]["items"]["xc"], d["position"]["items"]["yc"]].astype(float)
    e = d["plant_energy"]["properties"]["annual_energy_production"]
    return xy, float(e["default"]), np.asarray(e.get("binned", []), float)


def load_baseline(n):
    """Official example (baseline) layout iea37-ex<n>.yaml: (xy, official AEP MWh)."""
    xy, a, _ = load_layout(os.path.join(DATA, f"iea37-ex{n}.yaml"))
    return xy, a


def load_submission(par, n):
    """Participant submission iea37-par<par>-opt<n>.yaml: (xy, reported AEP MWh)."""
    xy, a, _ = load_layout(os.path.join(DATA, "iea37-cs1-results", f"iea37-par{par}-opt{n}.yaml"))
    return xy, a


if __name__ == "__main__":
    for n in (16, 36, 64):
        xy, ref = load_baseline(n)
        _, _, binned = load_layout(os.path.join(DATA, f"iea37-ex{n}.yaml"))
        a = aep(xy)
        print(f"ex{n}: AEP {a:.6f} official {ref:.5f} rel.err {(a - ref) / ref:.2e}; "
              f"max bin err {np.abs(aep_binned(xy) - binned).max():.2e}; ideal {ideal_aep(n):.4f}; "
              f"feasible {is_feasible(xy)}")
