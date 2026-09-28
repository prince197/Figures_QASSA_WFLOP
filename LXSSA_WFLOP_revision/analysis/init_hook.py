"""Initial-population hook shared by all optimizers.

GEN is None (default): uniform random coordinates in the bounds, exactly the original
np.random.uniform(lb, ub, (n, dim)) call, so earlier results are reproduced bit for bit.
GEN set to a callable gen(n, dim) -> (n, dim) array: feasibility-preserving initialization
(feasible_init.py). The generator draws from the global NumPy stream, so runs stay seed-paired.
"""
import numpy as np

GEN = None


def init_pop(n, dim, lb, ub):
    if GEN is None:
        return np.random.uniform(lb, ub, (n, dim))
    return np.asarray(GEN(n, dim), dtype=float)
