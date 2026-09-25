"""Additional baselines requested by Reviewer 3 (R3-3), with the same interface, seeding and
initial population as the authors' optimizers.

VNS   : continuous variable neighbourhood search (shake k turbines with a k-dependent radius,
        followed by a short first-improvement local search), in the spirit of basic VNS
        (Mladenovic and Hansen) adapted to continuous turbine coordinates.
MSSLSQP: multistart SLSQP, a gradient-based constrained solver. Wake loss is minimised with
        explicit spacing and boundary constraints; gradients are forward differences, and every
        objective evaluation (including those for gradients) counts against the budget.

Both draw the same initial population as the authors' optimizers for a given seed
(np.random.seed(seed); np.random.uniform(lb, ub, (pop, dim))), so runs are seed-paired.
"""
import numpy as np
from scipy.optimize import minimize


class BudgetExhausted(Exception):
    pass


class VNS:
    def __init__(self, pop_size=30, budget=3030, radius=500.0, k_max=3,
                 shake_frac=(0.05, 0.15, 0.40), ls_steps=10, seed=None):
        self.pop_size = pop_size; self.budget = budget; self.radius = radius
        self.k_max = k_max; self.shake_frac = shake_frac; self.ls_steps = ls_steps
        if seed is not None:
            np.random.seed(seed)

    def optimize(self, obj_fun, dim, lb, ub):
        pop = np.random.uniform(lb, ub, (self.pop_size, dim))
        fit = np.array([obj_fun(p) for p in pop])
        calls = self.pop_size
        x = pop[np.argmin(fit)].copy(); fx = fit.min()
        n = dim // 2
        step = 0.1 * self.radius
        k = 1
        while calls < self.budget:
            # shaking: move k randomly chosen turbines within radius shake_frac[k-1] * r
            y = x.copy()
            for t in np.random.choice(n, min(k, n), replace=False):
                y[2 * t:2 * t + 2] += np.random.normal(0, self.shake_frac[k - 1] * self.radius, 2)
            y = np.clip(y, lb, ub)
            fy = obj_fun(y); calls += 1
            # first-improvement local search: single-turbine moves with an adaptive step
            for _ in range(self.ls_steps):
                if calls >= self.budget:
                    break
                z = y.copy()
                t = np.random.randint(n)
                z[2 * t:2 * t + 2] += np.random.normal(0, step, 2)
                z = np.clip(z, lb, ub)
                fz = obj_fun(z); calls += 1
                if fz < fy:
                    y, fy = z, fz
                    step = min(step * 1.5, 0.5 * self.radius)
                else:
                    step = max(step * 0.8, 1e-3 * self.radius)
            if fy < fx:
                x, fx, k = y, fy, 1          # move and restart neighbourhoods
            else:
                k = k + 1 if k < self.k_max else 1
        return x, fx, None


class MSSLSQP:
    """Multistart SLSQP on the constrained problem. `wake_fun` returns the wake loss (no
    penalty); starts are the initial population sorted by the penalised objective, then new
    uniform points once the population is exhausted."""

    def __init__(self, wake_fun, radius, smin, pop_size=30, budget=3030, eps=1.0, seed=None,
                 boundary_cons=None):
        self.wake_fun = wake_fun; self.radius = radius; self.smin = smin
        self.boundary_cons = boundary_cons      # optional g(xy) >= 0 for a non-circular site
        self.pop_size = pop_size; self.budget = budget; self.eps = eps
        if seed is not None:
            np.random.seed(seed)

    def _cons(self, x):
        p = x.reshape(-1, 2)
        n = len(p)
        iu = np.triu_indices(n, 1)
        d2 = ((p[:, None] - p[None]) ** 2).sum(-1)[iu]
        if self.boundary_cons is not None:
            bnd = self.boundary_cons(p)
        else:
            bnd = (self.radius ** 2 - (p ** 2).sum(1)) / self.radius ** 2
        return np.concatenate([(d2 - self.smin ** 2) / self.smin ** 2, bnd])

    def optimize(self, obj_fun, dim, lb, ub):
        pop = np.random.uniform(lb, ub, (self.pop_size, dim))
        fit = np.array([obj_fun(p) for p in pop])
        self.calls = self.pop_size
        order = list(np.argsort(fit))
        best = [pop[order[0]].copy(), fit[order[0]]]

        def feasible(x):
            return self._cons(x).min() >= -1e-9

        def wrapped(x):
            if self.calls >= self.budget:
                raise BudgetExhausted
            self.calls += 1
            w = self.wake_fun(x)
            if feasible(x):
                if best[1] > 1e9 or w < best[1]:
                    best[0], best[1] = x.copy(), w
            return w

        starts = [pop[i] for i in order]
        while self.calls < self.budget:
            x0 = starts.pop(0) if starts else np.random.uniform(lb, ub, dim)
            try:
                minimize(wrapped, x0, method="SLSQP", bounds=[(lb, ub)] * dim,
                         constraints=[{"type": "ineq", "fun": self._cons}],
                         options=dict(maxiter=200, ftol=1e-9, eps=self.eps))
            except BudgetExhausted:
                break
        return best[0], best[1], None
