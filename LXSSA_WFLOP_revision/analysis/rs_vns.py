"""RS-VNS: random-multistart control for the two-phase hybrids.

Phase 1 replaces the swarm by pure random sampling with the same number of calls
(the seeded initial population of 30 followed by further random layouts, split*B calls in
total, drawn with init_hook.init_pop so that feasible initialization applies as well); the best
sampled layout is then refined by the same original basic VNS (original_vns.BVNS) for the rest
of the budget. Comparing SSA-VNS with RS-VNS isolates the value of the swarm phase.
"""
import numpy as np
from init_hook import init_pop
from original_vns import BVNS, BudgetExhausted


class RSVNS:
    def __init__(self, pop_size=30, budget=6030, split=0.5, radius=500.0, seed=None):
        self.pop_size = pop_size; self.budget = budget; self.split = split; self.radius = radius
        self.n1 = int(round(split * budget))
        if seed is not None:
            np.random.seed(seed)

    def optimize(self, obj_fun, dim, lb, ub):
        calls = [0]

        def f(x):
            if calls[0] >= self.budget:
                raise BudgetExhausted
            calls[0] += 1
            return obj_fun(x)

        pop = init_pop(self.pop_size, dim, lb, ub)
        fit = np.array([f(p) for p in pop])
        best, fbest = pop[np.argmin(fit)].copy(), fit.min()
        while calls[0] < self.n1:
            x = init_pop(1, dim, lb, ub)[0]
            fx = f(x)
            if fx < fbest:
                best, fbest = x.copy(), fx
        self.phase1_calls = calls[0]
        return BVNS(self.pop_size, self.budget, self.radius).search(f, best, fbest, dim, lb, ub)
