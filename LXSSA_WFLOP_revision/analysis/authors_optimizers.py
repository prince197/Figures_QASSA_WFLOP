"""Authors' optimizer implementations (GA, PSO, DE, SSA, LX-SSA), copied from the code supplied
by the authors in September 2026. Only the classes needed for the reviewer experiments are kept;
the update logic and the order of random-number calls are unchanged. All minimise obj_fun."""
import numpy as np


class GA:
    def __init__(self, pop_size=50, max_iter=100, pc=0.8, pm=0.1, seed=None):
        self.pop_size = pop_size; self.max_iter = max_iter; self.pc = pc; self.pm = pm
        if seed is not None:
            np.random.seed(seed)

    def optimize(self, obj_fun, dim, lb, ub):
        pop = np.random.uniform(lb, ub, (self.pop_size, dim))
        fitness = np.array([obj_fun(ind) for ind in pop])
        best_idx = np.argmin(fitness)
        best_pos = pop[best_idx].copy(); best_score = fitness[best_idx]
        curve = [best_score]
        for _ in range(self.max_iter):
            new_pop = []
            while len(new_pop) < self.pop_size:
                i1, i2 = np.random.choice(self.pop_size, 2, replace=False)
                p1 = pop[i1] if fitness[i1] < fitness[i2] else pop[i2]
                i1, i2 = np.random.choice(self.pop_size, 2, replace=False)
                p2 = pop[i1] if fitness[i1] < fitness[i2] else pop[i2]
                c1 = p1.copy(); c2 = p2.copy()
                if np.random.rand() < self.pc:
                    point = np.random.randint(1, dim)
                    c1[:point] = p1[:point]; c1[point:] = p2[point:]
                    c2[:point] = p2[:point]; c2[point:] = p1[point:]
                for child in [c1, c2]:
                    mask = np.random.rand(dim) < self.pm
                    child[mask] = np.random.uniform(lb, ub, np.sum(mask))
                    child = np.clip(child, lb, ub)
                    new_pop.append(child)
                    if len(new_pop) >= self.pop_size:
                        break
            pop = np.array(new_pop)
            fitness = np.array([obj_fun(ind) for ind in pop])
            idx = np.argmin(fitness)
            if fitness[idx] < best_score:
                best_score = fitness[idx]; best_pos = pop[idx].copy()
            curve.append(best_score)
        return best_pos, best_score, np.array(curve)


class PSO:
    def __init__(self, pop_size=50, max_iter=100, w=0.7, c1=2.0, c2=2.0, seed=None):
        self.pop_size = pop_size; self.max_iter = max_iter
        self.w = w; self.c1 = c1; self.c2 = c2
        if seed is not None:
            np.random.seed(seed)

    def optimize(self, obj_fun, dim, lb, ub):
        X = np.random.uniform(lb, ub, (self.pop_size, dim))
        V = np.zeros((self.pop_size, dim))
        pbest = X.copy()
        pbest_score = np.array([obj_fun(x) for x in X])
        gbest_idx = np.argmin(pbest_score)
        gbest = pbest[gbest_idx].copy(); gbest_score = pbest_score[gbest_idx]
        curve = [gbest_score]
        for _ in range(self.max_iter):
            for i in range(self.pop_size):
                r1 = np.random.rand(dim); r2 = np.random.rand(dim)
                V[i] = (self.w * V[i] + self.c1 * r1 * (pbest[i] - X[i])
                        + self.c2 * r2 * (gbest - X[i]))
                X[i] = X[i] + V[i]
                X[i] = np.clip(X[i], lb, ub)
                score = obj_fun(X[i])
                if score < pbest_score[i]:
                    pbest_score[i] = score; pbest[i] = X[i].copy()
                if score < gbest_score:
                    gbest_score = score; gbest = X[i].copy()
            curve.append(gbest_score)
        return gbest, gbest_score, np.array(curve)


class DE:
    def __init__(self, pop_size=50, max_iter=100, F=0.5, CR=0.9, seed=None):
        self.pop_size = pop_size; self.max_iter = max_iter; self.F = F; self.CR = CR
        if seed is not None:
            np.random.seed(seed)

    def optimize(self, obj_fun, dim, lb, ub):
        pop = np.random.uniform(lb, ub, (self.pop_size, dim))
        fitness = np.array([obj_fun(ind) for ind in pop])
        best_idx = np.argmin(fitness)
        best_pos = pop[best_idx].copy(); best_score = fitness[best_idx]
        curve = [best_score]
        for _ in range(self.max_iter):
            for i in range(self.pop_size):
                idxs = list(range(self.pop_size)); idxs.remove(i)
                r1, r2, r3 = np.random.choice(idxs, 3, replace=False)
                mutant = pop[r1] + self.F * (pop[r2] - pop[r3])
                mutant = np.clip(mutant, lb, ub)
                trial = pop[i].copy()
                j_rand = np.random.randint(dim)
                for j in range(dim):
                    if (np.random.rand() < self.CR) or (j == j_rand):
                        trial[j] = mutant[j]
                trial = np.clip(trial, lb, ub)
                trial_fit = obj_fun(trial)
                if trial_fit < fitness[i]:
                    pop[i] = trial; fitness[i] = trial_fit
                    if trial_fit < best_score:
                        best_score = trial_fit; best_pos = trial.copy()
            curve.append(best_score)
        return best_pos, best_score, np.array(curve)


class SSA:
    def __init__(self, pop_size=50, max_iter=100, seed=None):
        self.pop_size = pop_size; self.max_iter = max_iter
        if seed is not None:
            np.random.seed(seed)

    def optimize(self, obj_fun, dim, lb, ub):
        pop = np.random.uniform(lb, ub, (self.pop_size, dim))
        fitness = np.array([obj_fun(ind) for ind in pop])
        idx = np.argsort(fitness); pop = pop[idx]; fitness = fitness[idx]
        food_position = pop[0].copy(); food_fitness = fitness[0]
        curve = [food_fitness]
        for t in range(1, self.max_iter + 1):
            c1 = 2 * np.exp(-(4 * t / self.max_iter) ** 2)
            new_pop = pop.copy()
            for i in range(self.pop_size // 2):                       # leaders
                for j in range(dim):
                    c2 = np.random.rand(); c3 = np.random.rand()
                    step = c1 * ((ub - lb) * c2 + lb)
                    if c3 < 0.5:
                        new_pop[i, j] = food_position[j] + step
                    else:
                        new_pop[i, j] = food_position[j] - step
            for i in range(self.pop_size // 2, self.pop_size):         # followers
                new_pop[i] = (new_pop[i - 1] + pop[i]) / 2
            new_pop = np.clip(new_pop, lb, ub)
            pop = new_pop
            fitness = np.array([obj_fun(ind) for ind in pop])
            idx = np.argsort(fitness); pop = pop[idx]; fitness = fitness[idx]
            if fitness[0] < food_fitness:
                food_fitness = fitness[0]; food_position = pop[0].copy()
            curve.append(food_fitness)
        return food_position, food_fitness, np.array(curve)


class LXSSA:
    def __init__(self, pop_size=50, max_iter=100, phi=0.0, chi=1.0, seed=None):
        self.pop_size = pop_size; self.max_iter = max_iter; self.phi = phi; self.chi = chi
        if seed is not None:
            np.random.seed(seed)

    def laplace_random(self):
        z = np.random.rand()
        if z <= 0.5:
            gamma = self.phi + self.chi * np.log(2 * z)
        else:
            gamma = self.phi - self.chi * np.log(2 * (1 - z))
        return gamma

    def optimize(self, obj_fun, dim, lb, ub):
        pop = np.random.uniform(lb, ub, (self.pop_size, dim))
        fitness = np.array([obj_fun(ind) for ind in pop])
        idx = np.argsort(fitness); pop = pop[idx]; fitness = fitness[idx]
        food_position = pop[0].copy(); food_fitness = fitness[0]
        curve = [food_fitness]
        for l in range(1, self.max_iter + 1):
            r1 = 2 * np.exp(-(4 * l / self.max_iter) ** 2)
            new_pop = pop.copy()
            half = self.pop_size // 2
            for i in range(half):                                      # leaders
                for j in range(dim):
                    r2 = np.random.rand(); r3 = np.random.rand()
                    step = r1 * ((ub - lb) * r2 + lb)
                    if r3 >= 0.5:
                        new_pop[i, j] = food_position[j] + step
                    else:
                        new_pop[i, j] = food_position[j] - step
            for i in range(half, self.pop_size):                       # followers
                follower = (pop[i] + new_pop[i - 1]) / 2.0
                gamma = self.laplace_random()
                current_salp = pop[i]
                candidate = current_salp + gamma * (food_position - current_salp)
                follower = np.clip(follower, lb, ub)
                candidate = np.clip(candidate, lb, ub)
                if obj_fun(candidate) < obj_fun(follower):
                    new_pop[i] = candidate
                else:
                    new_pop[i] = follower
            new_pop = np.clip(new_pop, lb, ub)
            pop = new_pop
            fitness = np.array([obj_fun(ind) for ind in pop])
            idx = np.argsort(fitness); pop = pop[idx]; fitness = fitness[idx]
            if fitness[0] < food_fitness:
                food_fitness = fitness[0]; food_position = pop[0].copy()
            curve.append(food_fitness)
        return food_position, food_fitness, np.array(curve)
