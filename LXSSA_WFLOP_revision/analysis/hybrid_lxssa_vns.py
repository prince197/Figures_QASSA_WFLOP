"""LX-SSA-VNS: a two-phase hybrid of the Laplacian Salp Swarm Algorithm and variable
neighbourhood search for continuous WFLOP.

Phase 1 (global exploration): the original LX-SSA (authors' implementation, phi = 0, chi = 1,
population 30) runs for the first `split` fraction of the objective-call budget.
Phase 2 (intensification): VNS starts from the LX-SSA food source H and spends the remaining
budget. Shaking moves k in {1, 2, 3} randomly chosen turbines by Gaussian steps with standard
deviation 0.05r, 0.15r, 0.40r; every shake is followed by a ten-call first-improvement local
search of single-turbine Gaussian moves with an adaptive step (x1.5 on success, x0.8 on
failure). An improvement is accepted and resets k = 1, otherwise k increases cyclically.

Rationale: in the controlled study LX-SSA stagnates after about half of the budget, whereas VNS
keeps improving from a good incumbent, and LX-SSA explores many regions of the layout space
while VNS refines one. The random-number stream is seeded once, so for a given seed the
initial population equals that of all other methods (paired comparison).
"""
import numpy as np
from authors_optimizers import LXSSA


class LXSSAVNS:
    def __init__(self, pop_size=30, budget=6030, split=0.5, phi=0.0, chi=1.0, radius=500.0,
                 k_max=3, shake_frac=(0.05, 0.15, 0.40), ls_steps=10, seed=None):
        self.pop_size = pop_size; self.budget = budget; self.split = split
        self.phi = phi; self.chi = chi; self.radius = radius
        self.k_max = k_max; self.shake_frac = shake_frac; self.ls_steps = ls_steps
        # one LX-SSA iteration costs 2 * pop_size calls; the initial population pop_size calls
        self.iters1 = max(1, int(round((split * budget - pop_size) / (2 * pop_size))))
        if seed is not None:
            np.random.seed(seed)

    def optimize(self, obj_fun, dim, lb, ub):
        calls = [0]

        def f(x):
            calls[0] += 1
            return obj_fun(x)

        # ---- Phase 1: LX-SSA (no re-seeding: continues the seeded random stream) ----
        x, fx, _ = LXSSA(self.pop_size, self.iters1, self.phi, self.chi, seed=None).optimize(f, dim, lb, ub)
        x = np.array(x, dtype=float).copy()
        self.phase1_calls = calls[0]

        # ---- Phase 2: VNS intensification from the LX-SSA food source ----
        n = dim // 2
        step = 0.1 * self.radius
        k = 1
        while calls[0] < self.budget:
            y = x.copy()
            for t in np.random.choice(n, min(k, n), replace=False):
                y[2 * t:2 * t + 2] += np.random.normal(0, self.shake_frac[k - 1] * self.radius, 2)
            y = np.clip(y, lb, ub)
            fy = f(y)
            for _ in range(self.ls_steps):
                if calls[0] >= self.budget:
                    break
                z = y.copy()
                t = np.random.randint(n)
                z[2 * t:2 * t + 2] += np.random.normal(0, step, 2)
                z = np.clip(z, lb, ub)
                fz = f(z)
                if fz < fy:
                    y, fy = z, fz
                    step = min(step * 1.5, 0.5 * self.radius)
                else:
                    step = max(step * 0.8, 1e-3 * self.radius)
            if fy < fx:
                x, fx, k = y, fy, 1
            else:
                k = k + 1 if k < self.k_max else 1
        return x, fx, None
