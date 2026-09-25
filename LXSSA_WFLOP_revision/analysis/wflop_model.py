"""Re-implementation of the Kusiak--Song Jensen/Weibull WFLOP benchmark objective.

Used only to re-evaluate stored layouts and to run the sensitivity checks reported in the
revised manuscript. Parameters follow Section V-A of the manuscript.
"""
import numpy as np

R = 38.5            # rotor radius (m)
CT = 0.8            # thrust coefficient
K = 0.075           # Jensen wake-spreading constant
S_CI, S_R, S_CO = 3.5, 14.0, 25.0
P_RATED = 1500.0    # kW
LAM, ETA = 140.86, -500.0
BIN_WIDTH = 15.0    # legacy directional-bin-width factor (degrees)
THETA = np.deg2rad(np.arange(7.5, 360.0, 15.0))   # 24 bin mid-points

DS1_PSI = np.full(24, 13.0)
DS1_W = np.array([0.0] + [0.01] * 4 + [0.2, 0.6] + [0.01] * 16 + [0.0])
DS2_PSI = np.array([7, 5, 5, 5, 5, 4, 5, 6, 7, 7, 8, 9.5, 10, 8.5, 8.5, 6.5, 4.6, 2.6,
                    8, 5, 6.4, 5.2, 4.5, 3.9], dtype=float)
DS2_W = np.array([0.0002, 0.008, 0.0227, 0.0242, 0.0225, 0.0339, 0.0423, 0.029, 0.0617,
                  0.0813, 0.0994, 0.1394, 0.1839, 0.1115, 0.0765, 0.008, 0.0051, 0.0019,
                  0.0012, 0.001, 0.0017, 0.0031, 0.0097, 0.0317])
DATASETS = {1: (np.full(24, 2.0), DS1_PSI, DS1_W), 2: (np.full(24, 2.0), DS2_PSI, DS2_W)}


def _wcdf_c(s, k, psi):
    """Weibull survival function exp(-(s/psi)^k)."""
    return np.exp(-(s / psi) ** k)


def expected_power_linear(k, psi):
    """Kusiak--Song discretised expected power (kW) of one turbine, elementwise in psi."""
    s = np.arange(S_CI, S_R + 1e-9, 0.5)
    lo, hi = s[:-1], s[1:]
    mid = 0.5 * (lo + hi)
    k = np.asarray(k)[..., None]
    psi = np.asarray(psi)[..., None]
    lin = LAM * np.sum(mid * (_wcdf_c(lo, k, psi) - _wcdf_c(hi, k, psi)), axis=-1)
    lin += ETA * (_wcdf_c(S_CI, k, psi) - _wcdf_c(S_R, k, psi))[..., 0]
    return lin + P_RATED * _wcdf_c(S_R, k, psi)[..., 0]


def _power_curve_cubic(s, cutout):
    p = np.where(s < S_CI, 0.0,
                 np.where(s <= S_R, P_RATED * (s ** 3 - S_CI ** 3) / (S_R ** 3 - S_CI ** 3), P_RATED))
    if cutout:
        p = np.where(s > S_CO, 0.0, p)
    return p


def expected_power_cubic(k, psi, cutout=True, ds=0.01):
    """Expected power (kW) for a cubic (nonlinear) power curve, fine-grid quadrature."""
    s = np.arange(0.0, 40.0 + ds, ds)
    lo, hi = s[:-1], s[1:]
    p = _power_curve_cubic(0.5 * (lo + hi), cutout)
    k = np.asarray(k)[..., None]
    psi = np.asarray(psi)[..., None]
    return np.sum(p * (_wcdf_c(lo, k, psi) - _wcdf_c(hi, k, psi)), axis=-1)


def jensen_deficits(xy, theta=THETA):
    """Combined Jensen deficit delta_i(theta), shape (n_theta, N) (Lemmas 1-2, eq. 5)."""
    xy = np.asarray(xy, dtype=float)
    dx = xy[:, 0][:, None] - xy[:, 0][None, :]   # rho_i - rho_j
    dy = xy[:, 1][:, None] - xy[:, 1][None, :]
    c, s = np.cos(theta)[:, None, None], np.sin(theta)[:, None, None]
    a = R / K
    num = dx * c + dy * s + a
    den = np.sqrt((dx + a * c) ** 2 + (dy + a * s) ** 2)
    beta = np.arccos(np.clip(num / den, -1.0, 1.0))
    d = np.abs(dx * c + dy * s)
    n = len(xy)
    inwake = (beta < np.arctan(K)) & ~np.eye(n, dtype=bool)[None]
    dij = (1 - np.sqrt(1 - CT)) / (1 + K * d / R) ** 2
    return np.sqrt(np.sum(np.where(inwake, dij, 0.0) ** 2, axis=2))


def gaussian_deficits(xy, theta=THETA, kstar=0.04):
    """Bastankhah--Porte-Agel (2014) Gaussian deficit with the same RSS superposition."""
    xy = np.asarray(xy, dtype=float)
    D = 2 * R
    dx = xy[:, 0][:, None] - xy[:, 0][None, :]
    dy = xy[:, 1][:, None] - xy[:, 1][None, :]
    c, s = np.cos(theta)[:, None, None], np.sin(theta)[:, None, None]
    x = dx * c + dy * s                      # downstream distance of i behind j
    rr = np.abs(-dx * s + dy * c)            # radial offset from wake centreline
    beta = 0.5 * (1 + np.sqrt(1 - CT)) / np.sqrt(1 - CT)
    eps = 0.2 * np.sqrt(beta)
    sig = kstar * np.maximum(x, 0) / D + eps
    arg = np.clip(1 - CT / (8 * sig ** 2), 0.0, 1.0)
    dij = (1 - np.sqrt(arg)) * np.exp(-(rr / D) ** 2 / (2 * sig ** 2))
    n = len(xy)
    dij = np.where((x > 0) & ~np.eye(n, dtype=bool)[None], dij, 0.0)
    return np.sqrt(np.sum(dij ** 2, axis=2))


def farm_objective(xy, dataset, wake="jensen", curve="linear", **kw):
    """Benchmark objective (legacy units: 15 x expected farm power in kW) and ideal value."""
    k, psi, w = DATASETS[dataset]
    delta = jensen_deficits(xy) if wake == "jensen" else gaussian_deficits(xy, **kw)
    psi_i = psi[:, None] * (1 - delta)
    kk = np.broadcast_to(k[:, None], psi_i.shape)
    if curve == "linear":
        ep, ep0 = expected_power_linear(kk, psi_i), expected_power_linear(k, psi)
    else:
        cut = curve == "cubic_cutout"
        ep, ep0 = expected_power_cubic(kk, psi_i, cut), expected_power_cubic(k, psi, cut)
    n = len(xy)
    obj = BIN_WIDTH * np.sum(w[:, None] * ep)
    ideal = BIN_WIDTH * n * np.sum(w * ep0)
    return obj, ideal


def parse_coords(sval):
    return np.array([[float(v) for v in p.split()] for p in sval.split(";")])


def min_spacing(xy):
    d = np.sqrt(((xy[:, None, :] - xy[None, :, :]) ** 2).sum(-1))
    return d[np.triu_indices(len(xy), 1)].min()
