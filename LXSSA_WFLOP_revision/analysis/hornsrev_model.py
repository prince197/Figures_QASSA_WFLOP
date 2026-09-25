"""Horns Rev 1 case with the real turbine, wind climate and farm outline (R1-9).

Data (DTU PyWake 2.6.20, py_wake/examples/data/hornsrev1.py, MIT licence): installed layout of
80 Vestas V80 turbines, the V80 power and thrust-coefficient curves, and the 12-sector Weibull
wind climate. Model: Jensen top-hat wake with k = 0.04, hub-centre in-wake test, root-sum-square
superposition, thrust coefficient at the free-stream speed, 5-degree direction bins (sector
parameters of the nearest 30-degree sector) and 1 m/s speed bins from 3 to 25 m/s.
"""
import numpy as np

D = 80.0
RR = D / 2
KW = 0.04
HOURS = 8760.0

WS_TAB = np.arange(3.0, 26.0)
P_TAB = np.array([0.0, 66.6, 154.0, 282.0, 460.0, 696.0, 996.0, 1341.0, 1661.0, 1866.0, 1958.0,
                  1988.0, 1997.0, 1999.0] + [2000.0] * 9)                     # kW
CT_TAB = np.array([0.0, 0.818, 0.806, 0.804, 0.805, 0.806, 0.807, 0.793, 0.739, 0.709, 0.409,
                   0.314, 0.249, 0.202, 0.167, 0.14, 0.119, 0.102, 0.088, 0.077, 0.067, 0.06, 0.053])

SEC_F = np.array([3.597152, 3.948682, 5.167395, 7.000154, 8.364547, 6.43485,
                  8.643194, 11.77051, 15.15757, 14.73792, 10.01205, 5.165975])
SEC_F = SEC_F / SEC_F.sum()
SEC_A = np.array([9.176929, 9.782334, 9.531809, 9.909545, 10.04269, 9.593921,
                  9.584007, 10.51499, 11.39895, 11.68746, 11.63732, 10.08803])
SEC_K = np.array([2.392578, 2.447266, 2.412109, 2.591797, 2.755859, 2.595703,
                  2.583984, 2.548828, 2.470703, 2.607422, 2.626953, 2.326172])

WD = np.arange(0.0, 360.0, 5.0)                        # meteorological "from" direction
_SEC = (np.round(WD / 30.0).astype(int)) % 12          # nearest sector
F_WD = SEC_F[_SEC] / 6.0
A_WD, K_WD = SEC_A[_SEC], SEC_K[_SEC]
TOWARD = np.deg2rad(270.0 - WD)                         # direction the wind blows towards
WS = np.arange(3.0, 26.0)
_lo, _hi = WS - 0.5, WS + 0.5
P_WS = (np.exp(-(_lo[None] / A_WD[:, None]) ** K_WD[:, None])
        - np.exp(-(_hi[None] / A_WD[:, None]) ** K_WD[:, None]))       # (n_wd, n_ws)
A_CT = 1 - np.sqrt(1 - np.interp(WS, WS_TAB, CT_TAB))                  # (n_ws,)

WT_X = np.array([423974, 424042, 424111, 424179, 424247, 424315, 424384, 424452, 424534,
                 424602, 424671, 424739, 424807, 424875, 424944, 425012, 425094, 425162,
                 425231, 425299, 425367, 425435, 425504, 425572, 425654, 425722, 425791,
                 425859, 425927, 425995, 426064, 426132, 426214, 426282, 426351, 426419,
                 426487, 426555, 426624, 426692, 426774, 426842, 426911, 426979, 427047,
                 427115, 427184, 427252, 427334, 427402, 427471, 427539, 427607, 427675,
                 427744, 427812, 427894, 427962, 428031, 428099, 428167, 428235, 428304,
                 428372, 428454, 428522, 428591, 428659, 428727, 428795, 428864, 428932,
                 429014, 429082, 429151, 429219, 429287, 429355, 429424, 429492], float)
WT_Y = np.tile([6151447, 6150891, 6150335, 6149779, 6149224, 6148668, 6148112, 6147556], 10).astype(float)
I16 = [0, 1, 2, 3, 8, 9, 10, 11, 16, 17, 18, 19, 24, 25, 26, 27]


def site(n_turbines):
    """Installed layout (local coordinates, centred) and its outline polygon (counter-clockwise)."""
    idx = I16 if n_turbines == 16 else list(range(80))
    xy = np.c_[WT_X[idx], WT_Y[idx]]
    c = xy.mean(0)
    xy = xy - c
    rows = 4 if n_turbines == 16 else 10               # columns of 8 (or 4) turbines
    col = 4 if n_turbines == 16 else 8
    grid = xy.reshape(rows, col, 2)
    poly = np.array([grid[0, 0], grid[0, -1], grid[-1, -1], grid[-1, 0]])
    area = 0.5 * np.sum(poly[:, 0] * np.roll(poly[:, 1], -1) - np.roll(poly[:, 0], -1) * poly[:, 1])
    if area < 0:
        poly = poly[::-1]
    poly = poly * 1.001        # enlarge by 0.1% so the as-built (integer) coordinates lie inside
    return xy, poly


def outside_distance(xy, poly):
    """Distance (m) of each point outside the convex polygon (0 inside)."""
    a, b = poly, np.roll(poly, -1, axis=0)
    e = b - a
    nrm = np.c_[e[:, 1], -e[:, 0]] / np.linalg.norm(e, axis=1)[:, None]   # outward normals
    s = np.einsum("ijk,jk->ij", xy[:, None, :] - a[None], nrm)
    return np.maximum(s.max(1), 0.0)


def aep_gwh(xy, with_wake=True):
    """Annual energy production (GWh) of layout xy (N x 2, metres)."""
    xy = np.asarray(xy, float)
    n = len(xy)
    if with_wake:
        dx = xy[:, 0][:, None] - xy[:, 0][None, :]
        dy = xy[:, 1][:, None] - xy[:, 1][None, :]
        c = np.cos(TOWARD)[:, None, None]; s = np.sin(TOWARD)[:, None, None]
        x = dx[None] * c + dy[None] * s               # downstream distance of i behind j
        lat = np.abs(-dx[None] * s + dy[None] * c)
        inw = (x > 0) & (lat < RR + KW * x)
        g4 = np.where(inw, (RR / (RR + KW * np.maximum(x, 0))) ** 4, 0.0)
        G = np.sqrt(g4.sum(2))                         # (n_wd, n)
        u = WS[None, :, None] * (1 - A_CT[None, :, None] * G[:, None, :])   # (n_wd, n_ws, n)
    else:
        u = np.broadcast_to(WS[None, :, None], (len(WD), len(WS), n))
    p = np.interp(u, WS_TAB, P_TAB)
    e = np.einsum("dsn,ds,d->", p, P_WS, F_WD)       # expected farm power, kW
    return e * HOURS / 1e6


def make_objective(n_turbines, smin=4 * D, penalty=1e10):
    """Minimisation objective in the authors' form: wake loss (GWh) + quadratic penalties."""
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
