"""Verifier for wind-layout-yield.

Ground truth (tests/truth.npz) is the clean hourly hub-height wind speed and direction that
the mast record was synthesised from (authoring/provenance/generate.py); it is never produced
by a solver. The benchmark power in tests/benchmark.json is the best true mean power found by
a truth-informed optimisation (authoring/evidence/bench.py) or by any calibration solver.
"""
import csv
import json
import math
import os
import re

import numpy as np
import pytest

OUT = "/app/output"
HERE = os.path.dirname(os.path.abspath(__file__))
SITE = json.load(open(os.path.join(HERE, "site.json")))
BENCH = json.load(open(os.path.join(HERE, "benchmark.json")))
METRICS = "/logs/verifier/metrics.json"

GEOM_TOL_M = 0.1            # stated in the instruction
MAX_SHORTFALL = 0.0075      # layout quality: 1 - P_true / P_best
MAX_YIELD_ERR = 0.020       # yield accuracy: |P_reported / P_true - 1|

D = SITE["rotor_diameter_m"]
R0 = D / 2
CT = SITE["thrust_coefficient"]
K = SITE["wake_decay_constant"]
A0 = 1.0 - math.sqrt(1.0 - CT)

_metrics = {}


def _dump():
    os.makedirs(os.path.dirname(METRICS), exist_ok=True)
    with open(METRICS, "w") as fh:
        json.dump(_metrics, fh, indent=2)


def _norm(s):
    return re.sub(r"[^a-z0-9]", "", s.strip().lower())


def _num(v):
    if isinstance(v, bool):
        raise ValueError("boolean is not a number")
    x = float(str(v).strip())
    if not math.isfinite(x):
        raise ValueError("non-finite")
    return x


def load_layout():
    path = os.path.join(OUT, "layout.csv")
    with open(path, newline="") as fh:
        rows = [r for r in csv.reader(fh) if any(c.strip() for c in r)]
    head = [_norm(c) for c in rows[0]]
    xi = next(i for i, h in enumerate(head) if h in ("xm", "x", "easting", "eastingm"))
    yi = next(i for i, h in enumerate(head) if h in ("ym", "y", "northing", "northingm"))
    return np.array([[_num(r[xi]), _num(r[yi])] for r in rows[1:]], float)


def load_yield():
    obj = json.load(open(os.path.join(OUT, "yield.json")))
    keys = {_norm(k): v for k, v in obj.items()}
    return _num(keys["netmeanpowerkw"]), (_num(keys["grossmeanpowerkw"]) if "grossmeanpowerkw" in keys else None)


def load_curve():
    d = np.loadtxt(os.path.join(HERE, "power_curve.csv"), delimiter=",", skiprows=1)
    return d[:, 0], d[:, 1]


def pcurve(u, curve):
    v, p = curve
    out = np.interp(u, v, p)
    out[u > v[-1]] = 0.0
    return out


def true_mean_power(xy, chunk=3000):
    t = np.load(os.path.join(HERE, "truth.npz"), allow_pickle=False)
    U, theta = t["U"], t["theta"]
    curve = load_curve()
    dx = xy[None, :, 0] - xy[:, None, 0]     # [j, i] = x_i - x_j
    dy = xy[None, :, 1] - xy[:, None, 1]
    tot = 0.0
    for a in range(0, len(U), chunk):
        th = np.deg2rad(theta[a:a + chunk])[:, None, None]
        ex, ey = -np.sin(th), -np.cos(th)      # unit vector of travel for wind FROM theta
        down = dx * ex + dy * ey
        lat = np.abs(dy * ex - dx * ey)
        rw = R0 + K * down
        inw = (down > 0) & (lat <= rw)
        d = np.where(inw, A0 * (R0 / np.where(inw, rw, 1.0)) ** 2, 0.0)
        u = U[a:a + chunk, None] * (1.0 - np.sqrt((d ** 2).sum(axis=1)))
        tot += pcurve(u, curve).sum()
    gross = len(xy) * pcurve(U, curve).mean()
    return tot / len(U), gross


@pytest.fixture(scope="module")
def layout():
    return load_layout()


@pytest.fixture(scope="module")
def truth(layout):
    net, gross = true_mean_power(layout)
    _metrics.update(true_net_kw=net, true_gross_kw=gross, wake_efficiency=net / gross,
                    benchmark_kw=BENCH["best_true_net_mean_power_kw"])
    _dump()
    return net, gross


def test_outputs_parse(layout):
    net, gross = load_yield()
    assert layout.shape == (SITE["n_turbines"], 2), f"expected {SITE['n_turbines']} turbines, got {layout.shape[0]}"
    assert net > 0


def test_layout_feasible(layout):
    r = np.sqrt((layout ** 2).sum(axis=1))
    dist = np.sqrt(((layout[:, None] - layout[None]) ** 2).sum(-1))
    dist[np.diag_indices(len(layout))] = np.inf
    _metrics.update(max_radius_m=float(r.max()), min_spacing_m=float(dist.min()))
    _dump()
    assert r.max() <= SITE["boundary_radius_m"] + GEOM_TOL_M, f"turbine outside lease: r={r.max():.2f} m"
    assert dist.min() >= SITE["min_spacing_m"] - GEOM_TOL_M, f"spacing violated: {dist.min():.2f} m"


def test_layout_quality(layout, truth):
    shortfall = 1.0 - truth[0] / BENCH["best_true_net_mean_power_kw"]
    _metrics.update(layout_shortfall=shortfall)
    _dump()
    assert shortfall <= MAX_SHORTFALL, f"layout delivers {100 * shortfall:.2f}% less than the best known layout"


def test_yield_accuracy(layout, truth):
    net, gross = load_yield()
    err = net / truth[0] - 1.0
    _metrics.update(reported_net_kw=net, yield_error=err, reported_gross_kw=gross,
                    gross_error=(gross / truth[1] - 1.0) if gross else None)
    _dump()
    assert abs(err) <= MAX_YIELD_ERR, f"reported net mean power is off by {100 * err:+.2f}%"
