"""Minimisation objective used with the authors' optimizers: wake loss (ideal farm objective minus
benchmark objective) plus a constraint penalty.

penalty="authors" (default) is the exact penalty of the authors' objective.py
(objective_original.py): for every turbine outside the circle, (1 + 1e10 g)^2 with
g = x^2 + y^2 - r^2 (m^2), and for every pair closer than smin, (1 + 1e10 (smin - d))^2.
"linear" and "count" are the reconstructions used before objective.py was available."""
import numpy as np
from wflop_model import farm_objective, R

PENALTY = 1e10


def make_objective(dataset, radius, smin=8 * R, penalty="authors", curve="linear"):
    def f(x):
        xy = x.reshape(-1, 2)
        obj, ideal = farm_objective(xy, dataset, curve=curve)
        n = len(xy)
        dd = np.sqrt(((xy[:, None] - xy[None]) ** 2).sum(-1))[np.triu_indices(n, 1)]
        rr = np.sqrt((xy ** 2).sum(1))
        if penalty == "authors":
            gb = (xy ** 2).sum(1) - radius ** 2
            gs = smin - dd
            pen = (((1 + PENALTY * gb[gb > 0]) ** 2).sum() + ((1 + PENALTY * gs[gs > 0]) ** 2).sum())
        elif penalty == "count":
            pen = PENALTY * ((dd < smin - 1e-9).sum() + (rr > radius + 1e-9).sum())
        else:
            pen = PENALTY * (np.maximum(0, smin - dd).sum() + np.maximum(0, rr - radius).sum())
        return ideal - obj + pen
    return f
