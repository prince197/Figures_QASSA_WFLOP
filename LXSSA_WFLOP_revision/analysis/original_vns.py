"""Original basic variable neighbourhood search (BVNS).

Framework of Mladenovic and Hansen (1997) / Hansen and Mladenovic (2001), in the continuous
form of Mladenovic, Drazic, Kovacevic-Vujcic and Cangalovic (EJOR 191, 2008), with constraints
handled by the exterior penalty F_p (as in that paper):

    x <- initial solution; k <- 1
    repeat until the budget is exhausted
        x'  <- Shake(x, k)               random point, uniform in the k-th l_inf shell
                                         rho_{k-1} < ||x' - x||_inf <= rho_k around x
        x'' <- BestImprovement(x')       local search run to a local minimum
        if F(x'') < F(x): x <- x'', k <- 1    (NeighbourhoodChange)
        else:             k <- k + 1, and k <- 1 after k_max

Differences from the modified VNS used elsewhere in this study (extra_baselines.VNS):
  * shaking perturbs ALL 2N coordinates (nested l_inf neighbourhoods of the whole layout),
    not k selected turbines;
  * the local search is a complete best-improvement descent to a local minimum, not ten
    trial moves: a derivative-free compass (coordinate) search that evaluates all 2*dim moves
    +-h e_i, moves to the best improving one, and halves h when none improves, until
    h < h_min. (Gradient descent, used for smooth problems in Mladenovic et al. 2008, is not
    applicable because the Jensen top-hat wake makes F_p discontinuous.)

Defaults (not tuned): k_max = 5, rho_k = 0.1 k r, initial step h0 = 0.05 r, h_min = 1e-3 r.
The initial solution is the best member of the same seeded initial population of 30 that all
other methods draw, so runs are seed-paired; every objective call counts against the budget.
"""
import numpy as np


class BudgetExhausted(Exception):
    pass


class BVNS:
    def __init__(self, pop_size=30, budget=6030, radius=500.0, k_max=5, rho_frac=0.1,
                 h0_frac=0.05, hmin_frac=1e-3, seed=None):
        self.pop_size = pop_size; self.budget = budget; self.radius = radius
        self.k_max = k_max; self.rho = [rho_frac * k * radius for k in range(k_max + 1)]
        self.h0 = h0_frac * radius; self.hmin = hmin_frac * radius
        if seed is not None:
            np.random.seed(seed)

    def optimize(self, obj_fun, dim, lb, ub):
        calls = [0]

        def f(x):
            if calls[0] >= self.budget:
                raise BudgetExhausted
            calls[0] += 1
            return obj_fun(x)

        pop = np.random.uniform(lb, ub, (self.pop_size, dim))
        try:
            fit = np.array([f(p) for p in pop])
        except BudgetExhausted:
            raise RuntimeError("budget smaller than the initial population")
        return self.search(f, pop[np.argmin(fit)].copy(), fit.min(), dim, lb, ub)

    def search(self, f, x, fx, dim, lb, ub):
        """BVNS from the incumbent (x, fx). `f` must raise BudgetExhausted when the budget is
        used up. Returns the best point evaluated, its value and None."""
        best = [x.copy(), fx]                               # best point evaluated so far

        def shake(x, k):
            lo, hi = self.rho[k - 1], self.rho[k]
            while True:                       # uniform in the l_inf shell (rejection)
                d = np.random.uniform(-hi, hi, dim)
                if np.abs(d).max() > lo:
                    return np.clip(x + d, lb, ub)

        def best_improvement(y, fy):
            h = self.h0
            while h >= self.hmin:
                cand_x, cand_f = None, fy
                for i in range(dim):
                    for s in (h, -h):
                        z = y.copy(); z[i] = min(max(z[i] + s, lb), ub)
                        fz = f(z)
                        if fz < best[1]:
                            best[0], best[1] = z.copy(), fz
                        if fz < cand_f:
                            cand_x, cand_f = z, fz
                if cand_x is None:
                    h *= 0.5                  # no improving move: refine the step
                else:
                    y, fy = cand_x, cand_f    # best-improvement move
            return y, fy

        try:
            x, fx = best_improvement(x, fx)   # descend from the initial solution
            k = 1
            while True:
                y = shake(x, k)
                fy = f(y)
                if fy < best[1]:
                    best[0], best[1] = y.copy(), fy
                y, fy = best_improvement(y, fy)
                if fy < fx:
                    x, fx, k = y, fy, 1
                else:
                    k = k + 1 if k < self.k_max else 1
        except BudgetExhausted:
            pass
        # return the best point evaluated (equals x unless the budget ran out mid-descent)
        return best[0], best[1], None
