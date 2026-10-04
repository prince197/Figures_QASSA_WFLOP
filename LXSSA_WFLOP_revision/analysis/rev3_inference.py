"""Revision 3, topic `inference` (items B1, B3, B4, B5, B7 of the remaining-work report).

Usage (from the repository root):  python3 LXSSA_WFLOP_revision/analysis/rev3_inference.py [--quick]
                                    [--null-reps 2000] [--null-boot-case 10000] [--null-boot-seed 2000]

Data: the 68 benchmark cases (data sets I/II x r = 500/750/1000 m x N = 2..10/12/15), 6,030 evaluations, random
starts, seeds 1-30, loaded with mpce_inference_extra.load (the same per-run records as the paper's pipeline);
coordinates (rounded to 3 decimals in the legacy CSVs) are read from the same files for the violation ordering.
No optimizer is run. Helpers of mpce_inference_extra.py (case_stats, case_diffs, cm_test-like imputation, wil,
holm, seed_level, cr2, wild_level, webb_all, case_boot_means, boot_summary) are imported, not modified.

B1  All-run paired outcome (R4). Per case and seed: feasible beats infeasible (stored labels); two feasible runs are
    compared by the stored benchmark objective (tie if |difference| <= 1e-9 objective units -- the tie rule of the
    paper's run-level tests, mpce_results.wil; objective = 15 x expected farm power in kW, ~10^4); two infeasible
    runs tie (main variant) or, as a sensitivity, the run with the smaller normalized total violation
    V = sum_i max(0, x_i^2 + y_i^2 - r^2) / r^2 + sum_{i<j} max(0, l_min - l_ij) / l_min  (l_min = 308 m) wins
    (tie if |dV| <= 1e-5; V recomputed from the stored coordinates, which are rounded to 1e-3 m).
    Score (W + T/2) / n over the n = 68 x 30 = 2,040 seed pairs. 95% percentile intervals: (a) cluster bootstrap of
    the cases (all 30 seeds of a resampled case together), (b) two-stage (cases, then seeds within resampled cases),
    (c, JSON only) seed level (cases fixed, seeds resampled within cases); 20,000 resamples, fixed seed.
    Exact two-sided sign test over the discordant pairs (W vs. L; pairs treated as independent), Holm over the 7
    comparisons of PSO-VNS with the other main methods, and separately over the 5 sampling-control contrasts.
B3  Bootstrap TOST audit (R7): p_TOST (max of the add-one tail shares beyond -/+m) vs. 90% percentile-interval
    inclusion for all 18 pairs at the seed (independent and joint), case and cluster (CR2, wild) levels; null
    simulation at Delta = -m and +m for the two primary pairs at the case and the seed level.
B4  Synchronized seed resampling (R7): seed-level bootstrap with one seed-index vector for all cases vs. independent
    per case (all 18 pairs), per-seed t(29) on the benchmark-average difference, and a test of cross-case seed
    dependence (variance ratio of the per-seed benchmark average; correlation of data set I/II twins).
B5  Imputation sensitivity (R4): the case-mean Wilcoxon (Holm family of tab:friedman68, and of tab:ablation) with
    imputed maximal differences vs. with those cases dropped, each with zero differences dropped (paper) or with
    Pratt's zero handling.
B7  Equal-cluster target: Table S-eqclus recomputed from the run records, plus the restricted wild-cluster
    bootstrap-t (Webb six-point, all 6^6 draws enumerated) for the equal-cluster mean.

Outputs: analysis/rev3_inference.json, analysis/rev3_inference_tables.tex, analysis/rev3_inference.log (stdout).
"""
import os, re, sys, glob, json, math, time, argparse
import numpy as np, pandas as pd
from scipy.stats import wilcoxon, binomtest, rankdata, t as t_dist

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import mpce_inference_extra as X  # noqa: E402

CASE, KEY = X.CASE, X.KEY
M = X.EQ_MARGIN                                  # 0.05 pp (read from mpce_summary.json by mpce_inference_extra)
LMIN = 308.0                                     # 4D = 8R, R = 38.5 m
OBJ_TIE = 1e-9                                   # objective units (benchmark objective = 15 x expected power, kW)
V_TIE = 1e-5                                     # normalized violation
B1_BOOT, B1_SEED = 20000, 20261004
NULL_SEED = 20261005
SYNC_DIAG_PERM, SYNC_DIAG_SEED = 20000, 20261006
MAIN7 = [b for b in X.MAIN8 if b != "PSOBV"]     # PSOC, SSABV, SSA, LXSSA, DE, BVNS, SLSQP
CTRL5 = [("SSABV", "RSDVNS"), ("SSABV", "RSVNS"), ("RSDVNS", "RSVNS"), ("LXBV", "RSDVNS"), ("PSOBV", "RSDVNS")]
B1_PAIRS = [("PSOBV", b) for b in MAIN7] + CTRL5
PRIMARY = [("PSOBV", "PSOC"), ("SSABV", "RSDVNS")]
EQCLUS_PAIRS = [("PSOBV", "PSOC"), ("SSABV", "RSDVNS"), ("SSABV", "RSVNS"), ("RSDVNS", "RSVNS"), ("LXBV", "RSDVNS"),
                ("PSOBV", "RSDVNS")]
# Table S-eqclus as printed in SWEVO_supplement.tex (revision 2; reconstructed from rounded table means)
EQCLUS_PRINTED = {"PSOBV-PSOC": (-0.018, -0.067, 0.031, -0.014, -0.063, 0.035, [0.064, 0.043, -0.008, -0.063, -0.031, -0.090]),
                  "SSABV-RSDVNS": (0.024, -0.051, 0.099, 0.012, -0.067, 0.092, [-0.153, 0.016, 0.025, -0.013, 0.064, 0.137]),
                  "SSABV-RSVNS": (-0.068, -0.126, -0.011, -0.079, -0.142, -0.017, [-0.212, -0.055, -0.018, -0.119, -0.064, -0.008]),
                  "RSDVNS-RSVNS": (-0.093, -0.131, -0.054, -0.092, -0.125, -0.058, [-0.059, -0.071, -0.042, -0.106, -0.128, -0.145]),
                  "LXBV-RSDVNS": (0.087, 0.004, 0.170, 0.077, -0.005, 0.158, [-0.067, 0.041, 0.059, 0.053, 0.152, 0.222]),
                  "PSOBV-RSDVNS": (-0.306, -0.442, -0.171, -0.332, -0.478, -0.185, [-0.621, -0.233, -0.131, -0.454, -0.315, -0.235])}
LAB = X.LAB


def pk(a, b):
    return f"{a}-{b}"


def plab(a, b):
    return f"{LAB[a]} vs.\\ {LAB[b]}"


# ------------------------------------------------------------------ data
def load_coords(data_dir):
    """KEY -> Coordinates string, from exactly the files and with the duplicate rules of mpce_inference_extra.load."""
    parts = []
    use = lambda c: c != "Curve"
    for fn, drop in (("fresh_grid.csv", ["VNS", "PSO", "SLSQP"]), ("fresh_vgrid.csv", []), ("fresh_bgrid.csv", [])):
        d = X.std_cols(pd.read_csv(os.path.join(HERE, fn), usecols=use))
        parts.append(d[~d.Algorithm.isin(drop)])
    for exp in ("rsvns", "psoc", "psobv", "slsqp", "psosplit", "omega90", "rsdisc"):
        files = glob.glob(os.path.join(data_dir, f"mpce_{exp}_s*of*.csv"))
        by_k = {}
        for fn in files:
            m = re.search(rf"mpce_{exp}_s(\d+)of(\d+)\.csv$", os.path.basename(fn))
            if m:
                by_k.setdefault(int(m.group(2)), {})[int(m.group(1))] = fn
        complete = {k: v for k, v in by_k.items() if set(v) == set(range(k))}
        k = max(complete, key=lambda kk: sum(os.path.getsize(f) for f in complete[kk].values()))
        df = pd.concat([pd.read_csv(f, usecols=use) for f in sorted(complete[k].values())], ignore_index=True)
        parts.append(X.std_cols(df).drop_duplicates(KEY, keep="last"))
    A = pd.concat(parts, ignore_index=True).drop_duplicates(KEY, keep="first")
    A = A[(A.Budget == 6030) & (A.Init == "random") & A.Dataset.isin(["1", "2"])].reset_index(drop=True)
    return A[KEY + ["Objective", "Coordinates"]]


def violation(coord, r):
    xy = np.array([[float(v) for v in p.split()] for p in coord.split(";")])
    gb = np.maximum(0.0, (xy ** 2).sum(1) - r * r).sum() / (r * r)
    D = np.sqrt(((xy[:, None, :] - xy[None, :, :]) ** 2).sum(-1))
    iu = np.triu_indices(len(xy), 1)
    gs = np.maximum(0.0, LMIN - D[iu]).sum() / LMIN
    return float(gb + gs), float(gb), float(gs)


def run_matrix(G, alg, cases):
    """(n_cases, 30) arrays of objective, loss %, feasible, violation for one method (cases in the given order)."""
    R = G[G.Algorithm == alg].set_index(CASE + ["Seed"]).sort_index()
    obj = np.empty((len(cases), 30)); loss = np.empty_like(obj); fe = np.zeros_like(obj, bool); V = np.full_like(obj, np.nan)
    for i, c in enumerate(cases):
        g = R.loc[c]
        assert list(g.index) == list(range(1, 31)), (alg, c)
        obj[i] = g.Objective.values; loss[i] = g.LossPct.values; fe[i] = g.Feasible.values; V[i] = g.V.values
    return dict(obj=obj, loss=loss, feas=fe, V=V)


# ------------------------------------------------------------------ B1
def outcome(Ra, Rb, viol=False):
    """per case x seed: 1 (first wins), 0.5 (tie), 0 (first loses)."""
    fa, fb = Ra["feas"], Rb["feas"]
    O = np.full(fa.shape, 0.5)
    O[fa & ~fb] = 1.0; O[~fa & fb] = 0.0
    both = fa & fb
    d = Ra["obj"] - Rb["obj"]                                       # objective: larger is better
    O[both & (d > OBJ_TIE)] = 1.0; O[both & (d < -OBJ_TIE)] = 0.0
    if viol:
        nb = ~fa & ~fb
        dv = Ra["V"] - Rb["V"]                                      # violation: smaller is better
        O[nb & (dv < -V_TIE)] = 1.0; O[nb & (dv > V_TIE)] = 0.0
    return O


def boot_scores(O, B, seed, chunk=1000):
    """95% percentile intervals of the mean score: (a) cases resampled, (b) cases then seeds within, (c) seeds within
    fixed cases."""
    rng = np.random.default_rng(seed); nc, ns = O.shape
    cm = O.mean(1)
    a = np.empty(B); b = np.empty(B); c = np.empty(B)
    for s in range(0, B, chunk):
        k = min(chunk, B - s)
        ci = rng.integers(0, nc, (k, nc))
        a[s:s + k] = cm[ci].mean(1)
        si = rng.integers(0, ns, (k, nc, ns))
        b[s:s + k] = O[ci[:, :, None], si].mean((1, 2))
        si2 = rng.integers(0, ns, (k, nc, ns))
        c[s:s + k] = O[np.arange(nc)[None, :, None], si2].mean((1, 2))
    q = lambda v: [float(x) for x in np.quantile(v, [0.025, 0.975])]
    return dict(case=q(a), two_stage=q(b), seed=q(c))


def b1_block(G, S, cases, B=B1_BOOT):
    RM = {a: run_matrix(G, a, cases) for a in sorted({x for p in B1_PAIRS for x in p})}
    res = {}
    for j, (a, b) in enumerate(B1_PAIRS):
        Ra, Rb = RM[a], RM[b]
        row = dict(a=a, b=b, n_pairs=int(Ra["feas"].size),
                   feas_a_pct=float(100 * Ra["feas"].mean()), feas_b_pct=float(100 * Rb["feas"].mean()))
        d = X.case_diffs(S, a, b)
        row["cond_mean_dloss_pp"] = float(d.mean()); row["cond_n_cases"] = int(len(d))
        row["both_infeasible_pairs"] = int((~Ra["feas"] & ~Rb["feas"]).sum())
        for var, viol in (("main", False), ("violation", True)):
            O = outcome(Ra, Rb, viol)
            W, T, L = int((O == 1).sum()), int((O == 0.5).sum()), int((O == 0).sum())
            r = dict(W=W, T=T, L=L, score=float(O.mean()),
                     p_sign=float(binomtest(W, W + L, 0.5).pvalue) if W + L > 0 else 1.0,
                     case_score_gt_half=int((O.mean(1) > 0.5 + 1e-12).sum()),
                     case_score_lt_half=int((O.mean(1) < 0.5 - 1e-12).sum()))
            r["ci95"] = boot_scores(O, B, B1_SEED + 10 * j + (1 if viol else 0))
            row[var] = r
        res[pk(a, b)] = row
        m = row["main"]
        print(f"B1 {pk(a, b):14s} feas {row['feas_a_pct']:.1f}/{row['feas_b_pct']:.1f} W/T/L {m['W']}/{m['T']}/{m['L']} "
              f"score {m['score']:.3f} case {np.round(m['ci95']['case'], 3)} 2st {np.round(m['ci95']['two_stage'], 3)} "
              f"p {m['p_sign']:.3g} | viol score {row['violation']['score']:.3f}", flush=True)
    for var in ("main", "violation"):
        for fam, pairs in (("main7", [pk("PSOBV", b) for b in MAIN7]), ("controls5", [pk(a, b) for a, b in CTRL5])):
            for k, h in zip(pairs, X.holm([res[k][var]["p_sign"] for k in pairs])):
                res[k][var]["p_sign_holm_" + fam] = float(h)
    # comparison with conditional-quality ranks (tab:friedman68) and case-mean tests
    fr = X.friedman(S, X.MAIN8)
    order_rank = [a for a in fr["order"] if a != "PSOBV"]
    order_score = sorted(MAIN7, key=lambda b: res[pk("PSOBV", b)]["main"]["score"])     # lowest score = strongest rival
    order_score_v = sorted(MAIN7, key=lambda b: res[pk("PSOBV", b)]["violation"]["score"])
    cmp = dict(avg_rank=fr["avg_rank"], rank_order_rivals=order_rank, score_order_rivals_main=order_score,
               score_order_rivals_violation=order_score_v,
               same_order_main=order_rank == order_score, same_order_violation=order_rank == order_score_v)
    return res, cmp, RM


# ------------------------------------------------------------------ B3
def tost_counts(bm, m=M):
    B = len(bm); lo, hi = np.quantile(bm, [0.05, 0.95])
    k_lo, k_hi = int((bm <= -m).sum()), int((bm >= m).sum())
    p = max(k_lo + 1, k_hi + 1) / (B + 1)
    return dict(B=B, k_lo=k_lo, k_hi=k_hi, p=float(p), ci90=[float(lo), float(hi)], ci_in=bool(-m < lo and hi < m),
                p_rule=bool(p <= 0.05))


def audit_block(G, S, EL):
    """p <= 0.05 vs CI inclusion for every pair and level (recomputed bootstrap means; CR2 and wild from the stored
    summary of mpce_inference_extra, recomputed where cheap)."""
    rows = {}
    for key in X.EQ_PAIRS:
        a, b = key.split("-")
        d = X.case_diffs(S, a, b); dv = d.values
        bm_c = X.case_boot_means(dv)
        case = tost_counts(bm_c)
        stored = EL["pairs"][key]
        assert abs(case["p"] - stored["case"]["p_tost"]) < 1e-12 and max(abs(u - v) for u, v in zip(case["ci90"], stored["case"]["ci90"])) < 1e-12
        sl, bm_s = X.seed_level(G, a, b, d)
        seed = tost_counts(bm_s)
        assert abs(seed["p"] - sl["p_tost"]) < 1e-12
        jt = sl["joint"]; Bj = sl["resamples"]
        kj = round(jt["p_tost"] * (Bj + 1)) - 1
        joint = dict(B=Bj, k_max=int(kj), p=jt["p_tost"], ci90=jt["ci90"], ci_in=jt["equivalent"], p_rule=bool(jt["p_tost"] <= 0.05))
        cl = stored["cluster"]
        cr = dict(p=cl["cr2"]["p_tost"], ci90=cl["cr2"]["ci90"], ci_in=cl["cr2"]["equivalent"], p_rule=bool(cl["cr2"]["p_tost"] <= 0.05))
        wd = dict(p=cl["wild"]["p_tost"], ci90=cl["wild"]["ci90"], ci_in=bool(-M < cl["wild"]["ci90"][0] and cl["wild"]["ci90"][1] < M),
                  p_rule=bool(cl["wild"]["p_tost"] <= 0.05))
        rows[key] = dict(seed=seed, seed_joint=joint, case=case, cluster_cr2=cr, cluster_wild=wd, seed_level_full=sl)
    summ = {}
    for lev in ("seed", "seed_joint", "case", "cluster_cr2", "cluster_wild"):
        v = [r[lev] for r in rows.values()]
        summ[lev] = dict(n=len(v), equivalent_ci=sum(x["ci_in"] for x in v), p_le_05=sum(x["p_rule"] for x in v),
                         discordant=[k for k, r in rows.items() if r[lev]["ci_in"] != r[lev]["p_rule"]])
        if lev in ("seed", "case"):
            summ[lev]["k_equal_500"] = [k for k, r in rows.items() if 500 in (r[lev]["k_lo"], r[lev]["k_hi"])]
            # distance of the binding tail count from the critical count 500 (B = 10,000)
            summ[lev]["min_abs_k_minus_500"] = int(min(min(abs(r[lev]["k_lo"] - 500), abs(r[lev]["k_hi"] - 500)) for r in rows.values()))
    return rows, summ


def null_case(d, m_true, R, Bi, rng):
    """case-level null simulation: population = observed case differences shifted to mean m_true; a replicate draws n
    cases iid from it and applies the case-level percentile bootstrap TOST (Bi resamples)."""
    d0 = d - d.mean() + m_true; n = len(d0)
    out = dict(eq_ci=0, eq_p=0, rej_low=0, rej_high=0, rej_low_p=0, rej_high_p=0, t_eq=0)
    for _ in range(R):
        x = d0[rng.integers(0, n, n)]
        bm = x[rng.integers(0, n, (Bi, n))].mean(1)
        lo, hi = np.quantile(bm, [0.05, 0.95])
        k_lo, k_hi = (bm <= -M).sum(), (bm >= M).sum()
        pl, ph = (k_lo + 1) / (Bi + 1), (k_hi + 1) / (Bi + 1)
        out["eq_ci"] += bool(-M < lo and hi < M); out["eq_p"] += bool(max(pl, ph) <= 0.05)
        out["rej_low"] += bool(lo > -M); out["rej_high"] += bool(hi < M)
        out["rej_low_p"] += bool(pl <= 0.05); out["rej_high_p"] += bool(ph <= 0.05)
        se = x.std(ddof=1) / math.sqrt(n)
        pt = max(t_dist.sf((x.mean() + M) / se, n - 1), t_dist.cdf((x.mean() - M) / se, n - 1))
        out["t_eq"] += bool(pt <= 0.05)
    return out


def null_seed(Ra_loss, Ra_f, Rb_loss, Rb_f, m_true, R, Bi, rng, chunk=500):
    """seed-level null simulation: population per case = its 30 observed seed pairs, with the first method's losses
    shifted by a constant so that the benchmark estimand (mean over cases of the feasible-run case-mean difference)
    equals m_true; a replicate draws 30 seed pairs iid per case and applies the seed-level percentile bootstrap TOST
    (Bi resamples of the seed pairs within every case, independently per case)."""
    nc, ns = Ra_loss.shape
    fm = lambda L, F: np.array([L[i][F[i]].mean() for i in range(nc)])
    est = (fm(Ra_loss, Ra_f) - fm(Rb_loss, Rb_f)).mean()
    La = Ra_loss + (m_true - est)
    assert abs((fm(La, Ra_f) - fm(Rb_loss, Rb_f)).mean() - m_true) < 1e-12
    FLa, FLb = np.where(Ra_f, La, 0.0), np.where(Rb_f, Rb_loss, 0.0)
    Fa, Fb = Ra_f.astype(float), Rb_f.astype(float)
    off = (np.arange(nc) * ns)[:, None]
    out = dict(eq_ci=0, eq_p=0, rej_low=0, rej_high=0, rej_low_p=0, rej_high_p=0, fallback=0)
    t0 = time.time()
    for rep in range(R):
        if rep and rep % 250 == 0:
            print(f"   null seed level: {rep}/{R} replications, {time.time() - t0:.0f} s", flush=True)
        J = rng.integers(0, ns, (nc, ns)) + off                    # replicate sample (flat indices)
        xa, xfa, xb, xfb = FLa.ravel()[J], Fa.ravel()[J], FLb.ravel()[J], Fb.ravel()[J]
        full_a = xa.sum(1) / np.maximum(xfa.sum(1), 1); full_b = xb.sum(1) / np.maximum(xfb.sum(1), 1)
        bm = np.empty(Bi)
        for s in range(0, Bi, chunk):
            k = min(chunk, Bi - s)
            I = rng.integers(0, ns, (k, nc, ns)) + off[None]
            sa, na = xa.ravel()[I].sum(2), xfa.ravel()[I].sum(2)
            sb, nb = xb.ravel()[I].sum(2), xfb.ravel()[I].sum(2)
            out["fallback"] += int((na == 0).sum() + (nb == 0).sum())
            ma = np.where(na > 0, sa / np.maximum(na, 1), full_a[None]); mb = np.where(nb > 0, sb / np.maximum(nb, 1), full_b[None])
            bm[s:s + k] = (ma - mb).mean(1)
        lo, hi = np.quantile(bm, [0.05, 0.95])
        k_lo, k_hi = (bm <= -M).sum(), (bm >= M).sum()
        pl, ph = (k_lo + 1) / (Bi + 1), (k_hi + 1) / (Bi + 1)
        out["eq_ci"] += bool(-M < lo and hi < M); out["eq_p"] += bool(max(pl, ph) <= 0.05)
        out["rej_low"] += bool(lo > -M); out["rej_high"] += bool(hi < M)
        out["rej_low_p"] += bool(pl <= 0.05); out["rej_high_p"] += bool(ph <= 0.05)
    return out


def cp(k, n):
    ci = binomtest(int(k), int(n)).proportion_ci(0.95, method="exact")
    return [float(ci.low), float(ci.high)]


def null_block(G, S, RM_all, R, Bc, Bs, cases, Rc=None):
    res = {}
    for j, (a, b) in enumerate(PRIMARY):
        d = X.case_diffs(S, a, b)
        assert len(d) == len(cases)
        Ra, Rb = RM_all[a], RM_all[b]
        for side, mt in (("minus", -M), ("plus", M)):
            for lev in ("case", "seed"):
                rng = np.random.default_rng(NULL_SEED + 100 * j + (0 if side == "minus" else 10) + (0 if lev == "case" else 1))
                t0 = time.time()
                if lev == "case":
                    reps = Rc or R
                    o = null_case(d.values, mt, reps, Bc, rng); Bi = Bc
                else:
                    reps = R
                    o = null_seed(Ra["loss"], Ra["feas"], Rb["loss"], Rb["feas"], mt, R, Bs, rng); Bi = Bs
                binding = "rej_low" if side == "minus" else "rej_high"
                r = dict(reps=reps, inner_resamples=Bi, true_mean=mt, seconds=round(time.time() - t0, 1), **o,
                         rate_eq_ci=o["eq_ci"] / reps, rate_eq_p=o["eq_p"] / reps, ci95_eq_ci=cp(o["eq_ci"], reps),
                         ci95_eq_p=cp(o["eq_p"], reps), rate_binding=o[binding] / reps, ci95_binding=cp(o[binding], reps))
                if lev == "case":
                    r["rate_eq_t"] = o["t_eq"] / reps; r["ci95_eq_t"] = cp(o["t_eq"], reps)
                res[f"{pk(a, b)}|{lev}|{side}"] = r
                print(f"NULL {pk(a, b)} {lev:4s} {side:5s}: TOST rejection (CI rule) {r['rate_eq_ci']:.4f} {np.round(r['ci95_eq_ci'], 4)}"
                      f" p-rule {r['rate_eq_p']:.4f} binding one-sided {r['rate_binding']:.4f} [{r['seconds']} s]", flush=True)
    return res


# ------------------------------------------------------------------ B4
def sync_diag(RM_all, cases, rng):
    """cross-case dependence of the seed effects of the paired run differences (both feasible; PSO-VNS - PSO all
    feasible): E[c, s] = (L_a - L_b)[c, s] - case mean. Statistic 1: variance ratio of the per-seed benchmark average,
    Var_s(mean_c E) / ((1/n^2) sum_c Var_s E_c) (1 under independence across cases). Statistic 2: mean Pearson
    correlation over seeds between the twin cases (same r and N in data sets I and II, identical initial populations).
    Permutation p values: seeds permuted independently within every case (SYNC_DIAG_PERM permutations)."""
    out = {}
    twins = [(i, cases.index(("2", c[1], c[2]))) for i, c in enumerate(cases) if c[0] == "1" and ("2", c[1], c[2]) in cases]
    for a, b in PRIMARY:
        Ra, Rb = RM_all[a], RM_all[b]
        D = Ra["loss"] - Rb["loss"]; ok = Ra["feas"] & Rb["feas"]
        keep = ok.all(1)                                          # cases in which all 30 pairs are feasible
        E = D[keep] - D[keep].mean(1, keepdims=True)
        n = E.shape[0]
        def stats(E):
            vr = E.mean(0).var(ddof=1) / (E.var(1, ddof=1).sum() / n ** 2)
            Ek = {k: e for k, e in zip(np.flatnonzero(keep), E)}
            cs = [np.corrcoef(Ek[i], Ek[j])[0, 1] for i, j in twins if i in Ek and j in Ek and Ek[i].std() > 0 and Ek[j].std() > 0]
            return vr, float(np.mean(cs)), len(cs)
        vr, rc, ntw = stats(E)
        P = np.empty((SYNC_DIAG_PERM, 2))
        for k in range(SYNC_DIAG_PERM):
            Ep = np.take_along_axis(E, rng.permuted(np.tile(np.arange(30), (n, 1)), axis=1), axis=1)
            P[k] = stats(Ep)[:2]
        out[pk(a, b)] = dict(cases_all_feasible=int(n), variance_ratio=float(vr), p_variance_ratio=float((1 + (P[:, 0] >= vr).sum()) / (1 + SYNC_DIAG_PERM)),
                             perm_q95_variance_ratio=float(np.quantile(P[:, 0], 0.95)),
                             twin_mean_corr=rc, twin_pairs=ntw, p_twin_corr=float((1 + (P[:, 1] >= rc).sum()) / (1 + SYNC_DIAG_PERM)))
        print(f"SYNC {pk(a, b)}: variance ratio {vr:.3f} (perm p {out[pk(a, b)]['p_variance_ratio']:.3g}), twin corr {rc:+.3f} "
              f"over {ntw} twins (p {out[pk(a, b)]['p_twin_corr']:.3g})", flush=True)
    return out


def replay_initial_populations():
    """the seeded initial population (init_hook.init_pop after np.random.seed(seed), as in every driver): identical
    for the two data sets at the same (r, N) and an exact scaling by r for the same N."""
    import init_hook
    init_hook.GEN = None
    rep = dict(identical_dataset_twins=True, scaled_across_radii=True, checked=0)
    for s in (1, 7, 30):
        for n in (6, 10):
            pops = {}
            for r in (500, 750, 1000):
                np.random.seed(s); pops[r] = init_hook.init_pop(30, 2 * n, -r, r)
            for r in (750, 1000):
                rep["scaled_across_radii"] &= bool(np.allclose(pops[r] / r, pops[500] / 500, atol=1e-12))
            rep["checked"] += 1
    rep["note"] = ("the wind data set does not enter init_pop; the drivers (mpce_experiments.run_grid, full_grid_experiments) "
                   "seed the global MT19937 stream with the run's seed (np.random.seed) in the optimizer constructor and draw "
                   "the initial population first, so seed k gives the same uniform draws in every case")
    return rep


# ------------------------------------------------------------------ B5
def imputation_block(S, family, pairs):
    P = S.pivot_table(index=CASE, columns="Algorithm", values="Loss")
    Q = S.pivot_table(index=CASE, columns="Algorithm", values="Qualified").astype(float)
    rows = {}
    for a, b in pairs:
        qa, qb = Q[a] > 0.5, Q[b] > 0.5
        both = qa & qb
        d = (P[b] - P[a])[both].values                     # > 0: a better (as cm_test / case_mean_wilcoxon)
        big = 10 * (np.nanmax(np.abs(d)) if len(d) else 1.0) + 1.0
        dimp = np.r_[d, np.full(int((qa & ~qb).sum()), big), np.full(int((~qa & qb).sum()), -big)]
        r = dict(n_cases_imputed_version=int(len(dimp)), n_both=int(both.sum()), only_a=int((qa & ~qb).sum()),
                 only_b=int((~qa & qb).sum()), neither=int((~qa & ~qb).sum()), zeros=int((np.abs(d) <= 1e-9).sum()),
                 mean_dloss_pp=float(np.mean(-d)))
        for nm, v in (("imputed", dimp), ("dropped", d)):
            p, rb, nn = X.wil(v)
            vz = np.where(np.abs(v) <= 1e-9, 0.0, v)
            pp = float(wilcoxon(vz, zero_method="pratt").pvalue) if (np.abs(vz) > 0).any() else 1.0
            r[nm] = dict(p_wilcox=p, rb=rb, n_nonzero=nn, p_pratt=pp, wins=int((v > 1e-9).sum()), losses=int((v < -1e-9).sum()))
        rows[pk(a, b)] = r
    keys = list(rows)
    for nm in ("imputed", "dropped"):
        for z in ("wilcox", "pratt"):
            for k, h in zip(keys, X.holm([rows[k][nm]["p_" + z] for k in keys])):
                rows[k][nm]["p_" + z + "_holm"] = float(h)
                rows[k][nm]["verdict_" + z] = ("ns" if h >= 0.05 else ("A" if rows[k][nm]["wins"] > rows[k][nm]["losses"] else "B"))
    changed = [k for k in keys if len({rows[k][nm]["verdict_" + z] for nm in ("imputed", "dropped") for z in ("wilcox", "pratt")}) > 1]
    return dict(family=family, pairs=rows, verdict_changes=changed)


# ------------------------------------------------------------------ B7
def eqclus_block(S, clusters, EL):
    W = X.webb_all(len(clusters))
    res = {}
    for a, b in EQCLUS_PAIRS:
        d = X.case_diffs(S, a, b); dv = d.values
        cid = np.array([clusters.index((c[0], c[1])) for c in d.index])
        cm = np.array([dv[cid == g].mean() for g in range(len(clusters))])
        G = len(cm); est = cm.mean(); se = cm.std(ddof=1) / math.sqrt(G); q = t_dist.ppf(0.95, G - 1)
        t5 = [float(est - q * se), float(est + q * se)]
        pt = float(max(t_dist.sf((est + M) / se, G - 1), t_dist.cdf((est - M) / se, G - 1)))
        wild = X.wild_level(cm, np.arange(G), G, W)                  # one 'case' per cluster: CR2 = t5 variance (see notes)
        cw = X.cr2(dv, cid, G)
        cww = EL["pairs"][pk(a, b)]["cluster"]["wild"]
        pr = EQCLUS_PRINTED[pk(a, b)]
        r = dict(n_cases=int(len(dv)), case_weighted=dict(mean=float(dv.mean()), cr2_ci90=cw["cr2"]["ci90"], cr2_p_tost=cw["cr2"]["p_tost"],
                                                       wild_ci90=cww["ci90"], wild_p_tost=cww["p_tost"]),
                 equal_cluster=dict(mean=float(est), cluster_means=[float(v) for v in cm], t5_ci90=t5, t5_p_tost=pt, se=float(se),
                                    wild_ci90=wild["ci90"], wild_p_tost=wild["p_tost"], wild_min_margin=wild["min_margin_pp"],
                                    equivalent_t5=bool(-M < t5[0] and t5[1] < M), equivalent_wild=bool(wild["p_tost"] <= 0.05),
                                    nonsuperior_t5=bool(t5[0] > -M), nonsuperior_wild=bool(wild["ci90"][0] > -M)))
        got = [r["case_weighted"]["mean"], *cw["cr2"]["ci90"], est, *t5]
        r["printed_rev2"] = dict(values=list(pr[:6]), cluster_means=pr[6],
                                 max_abs_diff=float(max(abs(u - v) for u, v in zip(got, pr[:6]))),
                                 max_abs_diff_cluster_means=float(max(abs(u - v) for u, v in zip(cm, pr[6]))),
                                 matches_to_3_decimals=bool(all(f"{u:.3f}" == f"{v:.3f}" or abs(u - v) < 0.0011 for u, v in zip(got, pr[:6]))))
        res[pk(a, b)] = r
        print(f"EQCLUS {pk(a, b):14s} cw {dv.mean():+.4f} CR2 {np.round(cw['cr2']['ci90'], 4)} | eq {est:+.4f} t5 {np.round(t5, 4)} "
              f"wild {np.round(wild['ci90'], 4)} p {wild['p_tost']:.3f} | printed diff {r['printed_rev2']['max_abs_diff']:.4f}", flush=True)
    return res


# ------------------------------------------------------------------ LaTeX
def f3(v, sign=False):
    s = f"{v:+.3f}" if sign else f"{v:.3f}"
    return s.replace("-", "$-$")


def ci(v, d=3):
    return f"[{f3(v[0])}, {f3(v[1])}]"


def ptex(p):
    if p < 1e-3:
        m_, e = f"{p:.1e}".split("e")
        return f"${m_}\\times10^{{{int(e)}}}$"
    return f"{p:.3f}"


def tables(R):
    out = ["% Generated by analysis/rev3_inference.py (revision 3, topic inference). Do not edit by hand.", ""]
    # ---- B1
    b1 = R["B1"]["pairs"]
    L = []
    for k in [pk("PSOBV", b) for b in MAIN7] + ["mid"] + [pk(a, b) for a, b in CTRL5]:
        if k == "mid":
            L.append("\\midrule\n\\multicolumn{10}{@{}l}{\\emph{Primary pair SSA-VNS vs.\\ RSD-VNS and the other sampling-control contrasts of Table~\\ref{M-tab:equiv-main}}} \\\\")
            continue
        r = b1[k]; m = r["main"]; v = r["violation"]
        fam = "p_sign_holm_main7" if k.startswith("PSOBV-") and k.split("-")[1] in MAIN7 else "p_sign_holm_controls5"
        L.append(f"{plab(r['a'], r['b'])} & {r['feas_a_pct']:.1f}/{r['feas_b_pct']:.1f} & {f3(r['cond_mean_dloss_pp'], True)} & "
                 f"{m['W']}/{m['T']}/{m['L']} & {m['score']:.3f} & {ci(m['ci95']['case'])} & {ci(m['ci95']['two_stage'])} & "
                 f"{ptex(m[fam])} & {v['score']:.3f} & {ci(v['ci95']['case'])} \\\\")
    out += ["\\begin{table}[!htbp]", "\\centering",
            "\\caption{All-run paired outcome on the 68 benchmark cases (6,030 evaluations, random starts, seeds 1--30; "
            "$n=2{,}040$ seed pairs per comparison). Per case and seed, a feasible run beats an infeasible one, two feasible runs "
            "are compared by the benchmark objective (tie if the difference is at most $10^{-9}$ objective units), two infeasible "
            "runs tie. Score $(W+T/2)/n$ of the first method (0.5 = no difference; above 0.5 = first method better) with 95\\% "
            "percentile intervals from 20{,}000 resamples (fixed seed) of the cases (all seeds of a case together) and two-stage "
            "(cases, then seeds within cases). $p_{\\rm sign}$: exact two-sided sign test over the discordant pairs ($W$ vs.\\ $L$, "
            "pairs treated as independent), Holm-adjusted over the 7 comparisons of PSO-VNS (top) or over the 5 contrasts below. "
            "Separate outcomes: Feas., feasible runs (\\%) of the first/second method; $\\overline{\\Delta L}$, conditional mean "
            "wake-loss difference (pp, first minus second; negative = first better) over the cases in which both are qualified "
            "(Table~\\ref{M-tab:friedman68}). Last two columns: two infeasible runs ordered by the normalized violation "
            "$V=\\sum_i\\max(0,x_i^2+y_i^2-r^2)/r^2+\\sum_{i<j}\\max(0,\\ell_{\\min}-\\ell_{ij})/\\ell_{\\min}$ recomputed from "
            "the stored coordinates (rounded to $10^{-3}$~m; tie if $|\\Delta V|\\le10^{-5}$), with the case-bootstrap interval.}",
            "\\label{tab:S-r3-allrun}", "\\scriptsize\\setlength{\\tabcolsep}{2pt}",
            "\\begin{tabular}{@{}lccccccccc@{}}", "\\toprule",
            "& & & \\multicolumn{5}{c}{Infeasible pairs tie} & \\multicolumn{2}{c}{Ordered by $V$} \\\\",
            "\\cmidrule(lr){4-8}\\cmidrule(l){9-10}",
            "Pair (first vs.\\ second) & Feas.\\ (\\%) & $\\overline{\\Delta L}$ & $W/T/L$ & Score & Case 95\\% CI & Two-stage 95\\% CI & $p_{\\rm sign}$ & Score & Case 95\\% CI \\\\",
            "\\midrule", "\\multicolumn{10}{@{}l}{\\emph{Main comparison: PSO-VNS vs.\\ each method (Holm over 7)}} \\\\"] + L + [
            "\\bottomrule", "\\end{tabular}", "\\end{table}", ""]
    # ---- B3
    A = R["B3"]["audit_summary"]; NS = R["B3"]["null"]
    names = dict(seed="Seed, independent per case", seed_joint="Seed, one joint seed vector", case="Case (bootstrap over cases)",
                 cluster_cr2="Cluster, CR2 $t_5$", cluster_wild="Cluster, wild bootstrap-$t$")
    L = []
    for lev, nm in names.items():
        s = A[lev]
        L.append(f"{nm} & {s['n']} & {s['equivalent_ci']} & {s['p_le_05']} & {len(s['discordant'])} \\\\")
    L2 = []
    for (a, b) in PRIMARY:
        for lev in ("case", "seed"):
            cells = []
            for side in ("minus", "plus"):
                r = NS[f"{pk(a, b)}|{lev}|{side}"]
                cells.append(f"{100 * r['rate_eq_ci']:.2f} [{100 * r['ci95_eq_ci'][0]:.2f}, {100 * r['ci95_eq_ci'][1]:.2f}]")
            r0 = NS[f"{pk(a, b)}|{lev}|minus"]
            th = lambda v: f"{v:,}".replace(",", "{,}")
            L2.append(f"{plab(a, b)} & {lev} & {th(r0['reps'])} & {th(r0['inner_resamples'])} & {cells[0]} & {cells[1]} \\\\")
    out += ["\\begin{table}[!htbp]", "\\centering",
            "\\caption{Audit of the bootstrap TOST (margin $m=0.05$~pp; 68 cases, 6,030 evaluations, random starts, seeds 1--30; "
            "$\\Delta L$ first minus second method). Top: the 18 pairs of Tables~\\ref{tab:equivalence} and~\\ref{tab:X-equiv-levels} "
            "at each level; Eq.\\ (CI): 90\\% interval strictly inside $(-m,m)$; $p\\le0.05$: reported $p_{\\rm TOST}=\\max(k_-+1,k_++1)/(B+1)$, "
            "$k_-$ ($k_+$) = bootstrap means $\\le-m$ ($\\ge m$), $B=10{,}000$ (wild: enumerated tail shares, CR2: $t_5$). Bottom: "
            "empirical size of the equivalence decision (90\\% percentile interval inside $(-m,m)$) under a true mean "
            "difference of exactly $-m$ or $+m$, in \\% with exact 95\\% binomial intervals (nominal: at most 5\\%). Case level: "
            "the observed case differences shifted to mean $\\pm m$, 68 cases drawn with replacement per replication; seed "
            "level: the first method's run losses shifted by a constant so that the benchmark mean of the case-mean "
            "differences is $\\pm m$, 30 seed pairs drawn with replacement within every case; each replication then applies "
            "the procedure of the paper with the stated number of inner resamples.}",
            "\\label{tab:S-r3-tostaudit}", "\\scriptsize\\setlength{\\tabcolsep}{3pt}",
            "\\begin{tabular}{@{}lcccc@{}}", "\\toprule",
            "Level & Pairs & Eq.\\ (CI) & $p_{\\rm TOST}\\le0.05$ & Discordant \\\\", "\\midrule"] + L + [
            "\\bottomrule", "\\end{tabular}", "\\par\\medskip",
            "\\begin{tabular}{@{}llcccc@{}}", "\\toprule",
            "Pair & Level & Replications & Inner resamples & Size at $\\Delta=-m$ (\\%) & Size at $\\Delta=+m$ (\\%) \\\\", "\\midrule"] + L2 + [
            "\\bottomrule", "\\end{tabular}", "\\end{table}", ""]
    # ---- B4
    S4 = R["B4"]["pairs"]; D4 = R["B4"]["dependence"]
    L = []
    for key in X.EQ_PAIRS:
        a, b = key.split("-"); r = S4[key]
        pst = r.get("per_seed_t")
        L.append(f"{plab(a, b)} & {f3(r['mean'], True)} & {ci(r['indep']['ci90'])} & {r['indep']['min_margin']:.3f} & "
                 f"{'yes' if r['indep']['eq'] else 'no'} & {ci(r['joint']['ci90'])} & {r['joint']['min_margin']:.3f} & "
                 f"{'yes' if r['joint']['eq'] else 'no'} & " + (f"{ci(pst['ci90'])}" if pst else "--") + " \\\\")
        if key == "RSDVNS-RSVNS":
            L.append("\\midrule")
    dd = "; ".join(f"{plab(*k.split('-'))}: variance ratio {v['variance_ratio']:.2f} (permutation $p={v['p_variance_ratio']:.3f}$), "
                   f"mean twin correlation {v['twin_mean_corr']:+.3f} ($p={v['p_twin_corr']:.3f}$, {v['twin_pairs']} twins)".replace("-", "$-$", 0)
                   for k, v in D4.items())
    out += ["\\begin{table}[!htbp]", "\\centering",
            "\\caption{Seed-level (fixed-benchmark) equivalence with independent and with synchronized seed resampling (68 cases, "
            "6,030 evaluations, random starts, seeds 1--30; $\\overline{\\Delta L}$ first minus second method, pp; negative = first "
            "better; $m=0.05$~pp; 10{,}000 resamples, fixed seeds of Table~\\ref{tab:X-equiv-levels}). Independent: the 30 seed "
            "pairs resampled independently in every case (Table~\\ref{tab:X-equiv-levels}); joint: one seed-index vector drawn "
            "for all cases, so every case uses the same resampled seeds (seed $k$ seeds the same random stream, hence the same "
            "initial population draws, in every case). Per-seed $t$: 90\\% $t_{29}$ interval of the 30 per-seed benchmark "
            "averages (pairs with all runs feasible). $m_{\\min}$: smallest margin for equivalence; Eq.: 90\\% interval inside "
            "$\\pm m$. Cross-case dependence of the seed effects (permutation of seeds within cases, 20{,}000 permutations): "
            + dd + ".}",
            "\\label{tab:S-r3-syncseed}", "\\scriptsize\\setlength{\\tabcolsep}{2pt}",
            "\\begin{tabular}{@{}lcccccccc@{}}", "\\toprule",
            "& & \\multicolumn{3}{c}{Independent per case} & \\multicolumn{3}{c}{Joint seed vector} & \\\\",
            "\\cmidrule(lr){3-5}\\cmidrule(lr){6-8}",
            "Pair (first vs.\\ second) & $\\overline{\\Delta L}$ & 90\\% CI & $m_{\\min}$ & Eq. & 90\\% CI & $m_{\\min}$ & Eq. & Per-seed $t$ 90\\% CI \\\\",
            "\\midrule"] + L + ["\\bottomrule", "\\end{tabular}", "\\end{table}", ""]
    # ---- B5
    L = []
    for fam_key, title in (("main", "Main comparison (Table~\\ref{M-tab:friedman68}; PSO-VNS vs.\\ each method, Holm over 7)"),
                           ("ablation", "Component analysis (Table~\\ref{M-tab:ablation}; Holm over 15), contrasts with imputed cases")):
        F = R["B5"][fam_key]
        L.append(f"\\multicolumn{{9}}{{@{{}}l}}{{\\emph{{{title}}}}} \\\\")
        for k, r in F["pairs"].items():
            if fam_key == "ablation" and r["only_a"] + r["only_b"] == 0:
                continue
            a, b = k.split("-"); im, dr = r["imputed"], r["dropped"]
            L.append(f"{plab(a, b)} & {r['only_a']}/{r['only_b']} & {r['zeros']} & {ptex(im['p_wilcox_holm'])} & {ptex(im['p_pratt_holm'])} & "
                     f"{r['n_both']} & {ptex(dr['p_wilcox_holm'])} & {ptex(dr['p_pratt_holm'])} & "
                     f"{im['verdict_wilcox']}/{im['verdict_pratt']}/{dr['verdict_wilcox']}/{dr['verdict_pratt']} \\\\")
        if fam_key == "main":
            L.append("\\midrule")
    out += ["\\begin{table}[!htbp]", "\\centering",
            "\\caption{Sensitivity of the case-mean Wilcoxon signed-rank tests to the imputation of maximal differences and to "
            "the handling of zero differences (68 cases, 6,030 evaluations, random starts, seeds 1--30; case means of the "
            "feasible runs, qualification at 15 of 30 feasible runs). Imputed: the paper's test (a case in which only one "
            "method qualifies counts as a maximal difference in its favour; $n=68$). Dropped: only the cases in which both "
            "qualify. Wilcoxon: zero differences ($|d|\\le10^{-9}$~pp) dropped (paper); Pratt: zeros kept in the ranking. All "
            "$p$ values Holm-adjusted within the stated family and version. Only 1st/2nd: cases in which only the first / "
            "second method qualifies; Zeros: cases with $|d|\\le10^{-9}$ among those in which both qualify. Verdicts "
            "(imputed-Wilcoxon / imputed-Pratt / dropped-Wilcoxon / dropped-Pratt): A, first method significantly better; "
            "B, second; ns, not significant ($\\alpha=0.05$).}",
            "\\label{tab:S-r3-imputation}", "\\scriptsize\\setlength{\\tabcolsep}{2.5pt}",
            "\\begin{tabular}{@{}lcccccccc@{}}", "\\toprule",
            "& & & \\multicolumn{2}{c}{Imputed} & \\multicolumn{3}{c}{Dropped} & \\\\",
            "\\cmidrule(lr){4-5}\\cmidrule(lr){6-8}",
            "Pair (first vs.\\ second) & Only 1st/2nd & Zeros & Wilcoxon & Pratt & $n$ & Wilcoxon & Pratt & Verdicts \\\\",
            "\\midrule"] + L + ["\\bottomrule", "\\end{tabular}", "\\end{table}", ""]
    # ---- B7
    E = R["B7"]
    L = []
    for a, b in EQCLUS_PAIRS:
        r = E[pk(a, b)]; cw, eq = r["case_weighted"], r["equal_cluster"]
        L.append(f"{plab(a, b)} & {f3(cw['mean'], True)} & {ci(cw['cr2_ci90'])} & {ci(cw['wild_ci90'])} & {f3(eq['mean'], True)} & "
                 f"{ci(eq['t5_ci90'])} & {ci(eq['wild_ci90'])} & {eq['wild_p_tost']:.3f} \\\\")
        L.append(f"\\multicolumn{{8}}{{@{{}}r@{{}}}}{{\\emph{{cluster means:}} " + ", ".join(f3(v) for v in eq["cluster_means"]) + "} \\\\")
    out += ["\\begin{table}[!htbp]", "\\centering",
            "\\caption{Case-weighted and equal-cluster estimates at the cluster level, recomputed from the run records (68 "
            "cases in six (data set, radius) clusters, 6,030 evaluations, random starts, seeds 1--30; $\\overline{\\Delta L}$ first "
            "minus second method, pp; negative = first better; margin $m=0.05$~pp). Case-weighted: mean over the cases with "
            "the CR2 $t_5$ interval and the restricted wild-cluster bootstrap-$t$ interval of Table~\\ref{tab:X-equiv-levels}. "
            "Equal-cluster: mean of the six cluster means with the $t_5$ interval over these means (this equals the CR2 interval "
            "for the equal-cluster weighting) and the restricted wild-cluster bootstrap-$t$ interval (CR2-studentized, all "
            "$6^6=46{,}656$ Webb six-point draws enumerated; interval by test inversion); $p_{\\rm TOST}$ of the wild bootstrap. "
            "Cluster means in the order Data Set~I, $r=500$, 750, 1000~m; Data Set~II, $r=500$, 750, 1000~m.}",
            "\\label{tab:S-r3-eqclus}", "\\scriptsize\\setlength{\\tabcolsep}{2.5pt}",
            "\\begin{tabular}{@{}lccccccc@{}}", "\\toprule",
            "& \\multicolumn{3}{c}{Case-weighted} & \\multicolumn{4}{c}{Equal-cluster} \\\\",
            "\\cmidrule(lr){2-4}\\cmidrule(l){5-8}",
            "Pair (first vs.\\ second) & $\\overline{\\Delta L}$ & CR2 90\\% CI & Wild 90\\% CI & $\\overline{\\Delta L}$ & $t_5$ 90\\% CI & Wild 90\\% CI & Wild $p_{\\rm TOST}$ \\\\",
            "\\midrule"] + L + ["\\bottomrule", "\\end{tabular}", "\\end{table}", ""]
    return "\n".join(out)


def clean(o):
    if isinstance(o, dict):
        return {str(k): clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [clean(v) for v in o]
    if isinstance(o, (np.bool_,)):
        return bool(o)
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, np.floating):
        return float(o)
    return o


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default=HERE)
    ap.add_argument("--out-dir", default=HERE)
    ap.add_argument("--null-reps", type=int, default=2000)
    ap.add_argument("--null-reps-case", type=int, default=5000)
    ap.add_argument("--null-boot-case", type=int, default=10000)
    ap.add_argument("--null-boot-seed", type=int, default=10000)
    ap.add_argument("--b1-boot", type=int, default=B1_BOOT)
    ap.add_argument("--skip-null", action="store_true")
    args = ap.parse_args(argv)
    t0 = time.time()
    A = X.load(args.data_dir)
    G = A[A.Algorithm.isin(set(X.MAIN8) | set(X.ABL))].copy()
    S = X.case_stats(G, sorted(set(X.MAIN8) | set(X.ABL)), X.BASE_THR)
    cases = [tuple(c) for c in S[CASE].drop_duplicates().sort_values(CASE).itertuples(index=False)]
    clusters = sorted({(d, int(r)) for d, r, _ in cases})
    assert len(cases) == 68 and len(clusters) == 6
    # coordinates -> normalized violation of the infeasible runs (feasible runs: V not needed, set 0)
    C = load_coords(args.data_dir)
    G = G.merge(C[KEY + ["Coordinates", "Objective"]].rename(columns={"Objective": "ObjC"}), on=KEY, how="left", validate="1:1")
    assert G.Coordinates.notna().all() and np.allclose(G.ObjC, G.Objective, rtol=0, atol=0)
    G["V"] = 0.0
    inf = ~G.Feasible
    vv = [violation(c, r) for c, r in zip(G.Coordinates[inf], G.Radius[inf])]
    G.loc[inf, "V"] = [v[0] for v in vv]
    vinfo = dict(infeasible_runs=int(inf.sum()), recomputed_V_zero=int(sum(v[0] == 0 for v in vv)),
                 boundary_only=int(sum(v[1] > 0 and v[2] == 0 for v in vv)), spacing_only=int(sum(v[1] == 0 and v[2] > 0 for v in vv)),
                 both=int(sum(v[1] > 0 and v[2] > 0 for v in vv)), median_V=float(np.median([v[0] for v in vv])))
    # feasible runs recomputed from rounded coordinates (precision caveat): count labelled-feasible with V > 1e-6 relative
    print(f"loaded {len(G)} runs, {len(cases)} cases; infeasible runs {vinfo}", flush=True)
    EL = json.load(open(os.path.join(HERE, "mpce_summary_extra.json")))["equivalence_levels"]
    R = dict(generated=time.strftime("%Y-%m-%d %H:%M:%S"), script="rev3_inference.py", margin_pp=M,
             data=dict(n_runs=int(len(G)), n_cases=68, seeds="1-30", budget=6030, init="random", methods=sorted(set(G.Algorithm)),
                       loader="mpce_inference_extra.load", coordinates="legacy CSVs, rounded to 1e-3 m"),
             conventions=dict(direction="Delta L = first minus second (pp of the ideal objective); score = share of seed pairs won by the first method, ties 1/2",
                              objective_tie=OBJ_TIE, objective_tie_units="benchmark objective (15 x expected farm power in kW)",
                              violation=dict(definition="V = sum_i max(0, x_i^2+y_i^2-r^2)/r^2 + sum_{i<j} max(0, l_min-l_ij)/l_min, l_min = 308 m",
                                             tie=V_TIE, **vinfo)))
    # ---- B1
    t1 = time.time()
    b1, cmp, RM = b1_block(G, S, cases, args.b1_boot)
    R["B1"] = dict(pairs=b1, comparison_with_ranks=cmp, resamples=args.b1_boot, seed=B1_SEED,
                   ci_levels=dict(case="cases resampled (cluster bootstrap; all seeds of a case together)",
                                  two_stage="cases resampled, then the 30 seeds within each resampled case",
                                  seed="cases fixed, seeds resampled within cases (JSON only)"),
                   seconds=round(time.time() - t1, 1))
    print(f"B1 rank order of rivals {cmp['rank_order_rivals']} vs score order {cmp['score_order_rivals_main']}", flush=True)
    # ---- B5
    R["B5"] = dict(main=imputation_block(S, "tab:friedman68 (PSO-VNS vs 7 methods)", X.MAIN_PAIRS),
                   ablation=imputation_block(S, "tab:ablation (15 contrasts)", X.ABL_CONTR))
    MS = json.load(open(os.path.join(HERE, "mpce_summary.json")))
    dp = max(abs(math.log10(R["B5"]["main"]["pairs"][pk("PSOBV", b)]["imputed"]["p_wilcox"]) -
                 math.log10(MS["main"]["case_mean_wilcoxon"][b]["p"])) for b in MAIN7)
    dpa = max(abs(math.log10(R["B5"]["ablation"]["pairs"][k]["imputed"]["p_wilcox_holm"]) - math.log10(MS["ablation"]["case_mean"][k]["p_holm"]))
              for k in R["B5"]["ablation"]["pairs"])
    R["B5"]["reproduction"] = dict(max_abs_diff_log10p_main=float(dp), max_abs_diff_log10p_holm_ablation=float(dpa), ok=bool(max(dp, dpa) < 1e-6))
    print(f"B5 reproduction {R['B5']['reproduction']}; verdict changes main {R['B5']['main']['verdict_changes']}, "
          f"ablation {R['B5']['ablation']['verdict_changes']}", flush=True)
    # ---- B3 audit + B4 (seed level independent vs joint; all 18 pairs)
    t1 = time.time()
    rows, summ = audit_block(G, S, EL)
    R["B3"] = dict(audit_summary=summ, audit=dict((k, {l: v for l, v in r.items() if l != "seed_level_full"}) for k, r in rows.items()),
                   convention=("H0-: Delta <= -m is rejected iff the share k-/B of bootstrap means <= -m satisfies (k- + 1)/(B + 1) <= 0.05; "
                               "H0+: Delta >= m likewise with means >= +m; ties at the margin count against equivalence; p_TOST = max; "
                               "the CI rule requires -m < q_0.05 and q_0.95 < m (numpy linear interpolation). With B = 10,000, "
                               "p <= 0.05 iff k <= 499, which implies the CI rule; the two can differ only if exactly k = 500 "
                               "(then p = 0.0501 but q_0.05 interpolates between the 500th and 501st order statistics)."))
    R["B4"] = dict(pairs={}, resamples=X.BOOT_N, seed_independent=X.SEED_BOOT_SEED, seed_joint=X.SEED_BOOT_SEED + 1)
    for k, r in rows.items():
        sl = r["seed_level_full"]
        R["B4"]["pairs"][k] = dict(mean=sl["mean_dloss_pp"],
                                   indep=dict(ci90=sl["ci90"], min_margin=sl["min_margin_pp"], eq=sl["equivalent"], p_tost=sl["p_tost"],
                                              width=sl["ci90"][1] - sl["ci90"][0]),
                                   joint=dict(ci90=sl["joint"]["ci90"], min_margin=sl["joint"]["min_margin_pp"], eq=sl["joint"]["equivalent"],
                                              p_tost=sl["joint"]["p_tost"], width=sl["joint"]["ci90"][1] - sl["joint"]["ci90"][0]),
                                   per_seed_t=sl.get("per_seed_t"), fallback_resamples=sl["fallback_resamples"])
        q = R["B4"]["pairs"][k]
        q["width_ratio_joint_to_indep"] = q["joint"]["width"] / q["indep"]["width"]
        q["verdict_changed"] = q["joint"]["eq"] != q["indep"]["eq"]
        q["nonsup_changed"] = (q["joint"]["ci90"][0] > -M) != (q["indep"]["ci90"][0] > -M)
    R["B4"]["any_verdict_changed"] = [k for k, q in R["B4"]["pairs"].items() if q["verdict_changed"] or q["nonsup_changed"]]
    R["B4"]["width_ratio_range"] = [min(q["width_ratio_joint_to_indep"] for q in R["B4"]["pairs"].values()),
                                    max(q["width_ratio_joint_to_indep"] for q in R["B4"]["pairs"].values())]
    R["B4"]["dependence"] = sync_diag(RM, cases, np.random.default_rng(SYNC_DIAG_SEED))
    R["B4"]["initial_population_replay"] = replay_initial_populations()
    R["B4"]["seconds"] = round(time.time() - t1, 1)
    print(f"B3 audit {json.dumps(summ)}", flush=True)
    print(f"B4 verdict changes {R['B4']['any_verdict_changed']}; width ratio {np.round(R['B4']['width_ratio_range'], 3)}", flush=True)
    # ---- B7
    R["B7"] = eqclus_block(S, clusters, EL)
    # ---- B3 null simulation
    if not args.skip_null:
        t1 = time.time()
        R["B3"]["null"] = null_block(G, S, RM, args.null_reps, args.null_boot_case, args.null_boot_seed, cases,
                                    Rc=args.null_reps_case)
        R["B3"]["null_design"] = dict(reps_seed=args.null_reps, reps_case=args.null_reps_case, inner_case=args.null_boot_case, inner_seed=args.null_boot_seed, seed=NULL_SEED,
                                      seconds=round(time.time() - t1, 1))
    R["seconds"] = round(time.time() - t0, 1)
    json.dump(clean(R), open(os.path.join(args.out_dir, "rev3_inference.json"), "w"), indent=1)
    if not args.skip_null:
        open(os.path.join(args.out_dir, "rev3_inference_tables.tex"), "w").write(tables(R))
    print(f"done in {R['seconds']} s", flush=True)


if __name__ == "__main__":
    main()
