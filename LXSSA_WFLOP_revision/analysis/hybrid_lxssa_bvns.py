"""LX-SSA-VNS (final version): two-phase hybrid of the Laplacian Salp Swarm Algorithm and the
original basic variable neighbourhood search (BVNS, original_vns.py).

Phase 1 (global exploration): the original LX-SSA (authors' implementation, phi = 0, chi = 1,
population 30) runs for the first `split` fraction of the objective-call budget
(T1 = round((split*B - Np) / (2 Np)) iterations; one LX-SSA iteration costs 2 Np calls).
Phase 2 (intensification): BVNS starts from the LX-SSA food source H - a complete
best-improvement local search from H, then shaking in nested l_inf neighbourhoods of the whole
layout, each followed by the complete local search - until the budget is exhausted.

`phase1="SSA"` gives the ablation variant SSA-VNS (standard SSA, T1 = round((split*B - Np)/Np)).
The random stream is seeded once and phase 1 draws the initial population first, so for a
given seed all methods start from the same initial population (paired comparison).
"""
import numpy as np
from authors_optimizers import LXSSA, SSA
from original_vns import BVNS, BudgetExhausted


class HybridBVNS:
    def __init__(self, pop_size=30, budget=6030, split=0.5, phi=0.0, chi=1.0, radius=500.0,
                 phase1="LXSSA", seed=None):
        self.pop_size = pop_size; self.budget = budget; self.split = split
        self.phi = phi; self.chi = chi; self.radius = radius; self.phase1 = phase1
        per_iter = 2 * pop_size if phase1 == "LXSSA" else pop_size
        self.iters1 = max(1, int(round((split * budget - pop_size) / per_iter)))
        if seed is not None:
            np.random.seed(seed)

    def optimize(self, obj_fun, dim, lb, ub):
        calls = [0]

        def f(x):
            if calls[0] >= self.budget:
                raise BudgetExhausted
            calls[0] += 1
            return obj_fun(x)

        # ---- Phase 1: LX-SSA or SSA (no re-seeding: continues the seeded random stream) ----
        if self.phase1 == "LXSSA":
            opt = LXSSA(self.pop_size, self.iters1, self.phi, self.chi, seed=None)
        else:
            opt = SSA(self.pop_size, self.iters1, seed=None)
        h, fh, _ = opt.optimize(f, dim, lb, ub)
        self.phase1_calls = calls[0]
        # ---- Phase 2: original basic VNS from the food source ----
        return BVNS(self.pop_size, self.budget, self.radius).search(
            f, np.array(h, dtype=float).copy(), fh, dim, lb, ub)
