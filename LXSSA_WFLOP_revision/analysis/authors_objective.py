"""Minimisation objective in the form used by the authors' optimizers: wake loss (ideal farm
objective minus benchmark objective) plus a large constraint penalty. The authors' own
objective.py was not supplied; the penalty form below is a reconstruction (see calibration)."""
import numpy as np
from wflop_model import farm_objective, R

PENALTY = 1e10


def make_objective(dataset, radius, smin=8 * R, penalty="count", curve="linear"):
    def f(x):
        xy = x.reshape(-1, 2)
        obj, ideal = farm_objective(xy, dataset, curve=curve)
        n = len(xy)
        dd = np.sqrt(((xy[:, None] - xy[None]) ** 2).sum(-1))[np.triu_indices(n, 1)]
        rr = np.sqrt((xy ** 2).sum(1))
        if penalty == "count":
            pen = PENALTY * ((dd < smin - 1e-9).sum() + (rr > radius + 1e-9).sum())
        else:
            pen = PENALTY * (np.maximum(0, smin - dd).sum() + np.maximum(0, rr - radius).sum())
        return ideal - obj + pen
    return f
