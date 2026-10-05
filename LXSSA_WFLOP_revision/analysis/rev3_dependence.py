"""Dependence-aware inference for the all-run paired endpoint (synchronized seeds) and boundary-aware intervals.

Usage (from any directory):  python3 analysis/rev3_dependence.py [--calib-reps 5000] [--calib-inner 10000]
                                                                 [--skip-calib] [--tables-only] [--out-dir DIR]
Outputs (analysis/ by default): rev3_dependence.json, rev3_dependence_tables.tex (tab:S-sa-dep-bench,
tab:S-sa-dep-studies), rev3_dependence.log (stdout, if redirected). No optimizer is run; only stored run records are
read, through the loaders of the scripts that generate the paper's tables (imported, not modified):
  benchmark    mpce_inference_extra.load + rev3_inference.run_matrix / outcome      (tab:S-sa-allrun)
  constraint   rev3_constraint_analysis.load / cube / outcome                        (tab:S-sa-constraint)
  1-deg study  rev3_fine_analysis.read_shards / in_cases / arm (tie rule 1e-6)       (tab:S-sa-fine)
  Lillgrund    analysis/rev2_lg16_s*of*.csv, rev2_lg16b_s*of*.csv (rule of audit_archive.reliability; checked
               against SWEVO_rev2/validation/paired_outcomes.csv)                    (tab:archive-paired)

All-run paired outcome (unchanged): per case and seed, a feasible run beats an infeasible one, two feasible runs are
compared by the objective (tie within the tie rule of the generating script), two infeasible runs tie; score
(W + T/2)/n. The W/T/L counts are asserted to equal those of the generating scripts' JSON outputs.

Seed-level (fixed-benchmark) inference. Seed k seeds the same random stream in every case of a study (benchmark: seeds
1-30 over 68 cases; constraint study: 31-60 over 6 cases; 1-deg study: 31-60 over 8 cases), so the outcomes of
different cases with the same seed are not independent; the independent replicates are the seed vectors. For seed k
the seed score is s_k = mean over the cases of the outcome; the score is mean_k s_k (= (W + T/2)/n).
  * 95% percentile interval: the complete seed vectors resampled jointly (B = 10,000, fixed RNG seed).
  * Exact sign-flip randomization test of mean_k (s_k - 1/2) = 0, all 2^K sign patterns enumerated (meet in the middle;
    integer arithmetic). Assumption (exchangeability): under H0, swapping the two methods' labels within a seed leaves
    the joint distribution of that seed's run pairs unchanged, so s_k - 1/2 is symmetric about 0 and independent over
    seeds; this is the null of no difference between the methods, not merely a zero mean score.
  * Wilcoxon signed-rank test of the K seed scores against 1/2 (check; zero differences dropped).
  * Holm adjustment within the family of the generating table: benchmark, the 7 comparisons of PSO-VNS and the 5
    sampling-control contrasts; constraint study, the 28 tests of each variant (4 methods x 6 cases and pooled), in
    which the pooled sign test is replaced by the seed-level test and the per-case exact sign tests (30 independent
    seeds each) are kept; 1-deg study, the 4 comparisons of PSO-VNS per arm.
Designed-group sensitivity (benchmark): the six (data set, radius) groups of rev3_inference / mpce_inference_extra;
group score = outcome mean over the cases and seeds of a group; exact sign-flip test of the case-weighted score over
the 2^6 = 64 group sign patterns (smallest attainable two-sided p = 1/32), t_5 interval of the equal-group mean.
Matched 1-deg comparison: direct 1-deg arm minus the 15-deg control arm (same seeds 31-60 and cases), wake loss under
the 1-deg objective over the jointly feasible seed pairs (case means, then mean over the 8 cases), joint-seed 95%
interval and sign-flip test of the per-seed case-average difference; matched all-run score (direct vs control run of
the same method and seed).
Calibration (optional): the joint-seed percentile-bootstrap TOST of tab:S-sa-syncseed for PSO-VNS vs PSO and SSA-VNS
vs RSD-VNS under a shifted null (first method's losses shifted so the benchmark estimand is -m or +m), replicate
samples drawn as 30 complete seed vectors with replacement (the dependence structure kept), inner joint bootstrap as
in mpce_inference_extra.seed_level (fallback to the replicate's case mean if a resample has no feasible run).
Lillgrund: exact (Clopper-Pearson) 95% interval of the conditional win probability W/(W+L) for every row; with ties
this differs from the score (W + T/2)/n.
"""
import os, sys, json, math, time, glob, argparse, itertools
import numpy as np, pandas as pd
from scipy.stats import wilcoxon, binomtest, t as t_dist

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import mpce_inference_extra as X            # noqa: E402
import rev3_inference as RI                 # noqa: E402
import rev3_constraint_analysis as RCA      # noqa: E402
import rev3_fine_analysis as RF             # noqa: E402

BOOT_B, BOOT_SEED = 10000, 20261007
CALIB_SEED = 20261008
M = X.EQ_MARGIN
ALPHA = 0.05
LAB = dict(X.LAB, GA="GA")
VALID = os.path.join(ROOT, "SWEVO_rev2", "validation")


# ------------------------------------------------------------------ core procedures
def flip_sums(v):
    s = np.zeros(1, dtype=v.dtype)
    for x in v:
        s = np.concatenate([s + x, s - x])
    return s


def signflip_exact(a, tol=0.0):
    """Exact two-sided sign-flip p of sum(a) (all 2^K patterns, meet in the middle). a: integers (exact) or floats
    (tol absorbs rounding). Returns (p_two_sided, p_greater, p_less)."""
    a = np.asarray(a)
    K = len(a); h = K // 2
    s1 = flip_sums(a[:h]); s2 = np.sort(flip_sums(a[h:]))
    t = a.sum(); N = float(2 ** K)
    ge = lambda thr: (len(s2) - np.searchsorted(s2, thr - s1 - tol, side="left")).sum()   # #(x + y >= thr)
    le = lambda thr: np.searchsorted(s2, thr - s1 + tol, side="right").sum()            # #(x + y <= thr)
    p_gt = ge(t) / N; p_lt = le(t) / N
    at = abs(t)
    p2 = 1.0 if at <= tol else (ge(at) + le(-at)) / N
    return float(min(1.0, p2)), float(p_gt), float(p_lt)


def wil_p(x):
    x = np.asarray(x, float); x = x[np.abs(x) > 1e-12]
    if len(x) == 0:
        return 1.0
    return float(wilcoxon(x).pvalue)


def seed_inference(O, ncase_units, boot_seed, B=BOOT_B):
    """O: outcomes [case, seed] in {0, 1/2, 1}. Seed scores, joint-seed percentile CI, exact sign flip, Wilcoxon."""
    nc, ns = O.shape
    s = O.mean(0)
    a = np.rint(2 * O.sum(0)).astype(np.int64) - nc                     # 2 nc (s_k - 1/2), integer
    assert np.allclose(a, 2 * nc * (s - 0.5))
    p2, pg, pl = signflip_exact(a)
    bm = s[np.random.default_rng(boot_seed).integers(0, ns, (B, ns))].mean(1)
    se = s.std(ddof=1) / math.sqrt(ns); q = t_dist.ppf(0.975, ns - 1)
    return dict(n_seeds=int(ns), n_cases=int(nc), score=float(O.mean()), seed_scores=[float(v) for v in s],
                seed_score_sd=float(s.std(ddof=1)), seeds_above=int((s > 0.5 + 1e-12).sum()),
                seeds_below=int((s < 0.5 - 1e-12).sum()), seeds_equal=int((np.abs(s - 0.5) <= 1e-12).sum()),
                ci95_joint_seed=[float(v) for v in np.quantile(bm, [0.025, 0.975])], boot_B=B, boot_seed=boot_seed,
                t29_ci95=[float(s.mean() - q * se), float(s.mean() + q * se)],
                p_signflip=p2, p_signflip_greater=pg, p_signflip_less=pl, p_signflip_method=f"exact, 2^{ns} patterns",
                p_wilcoxon=wil_p(s - 0.5), seeds_all_won=int((s >= 1 - 1e-12).sum()), seeds_all_lost=int((s <= 1e-12).sum()),
                ci_degenerate=bool(np.ptp(bm) == 0),
                # boundary-aware bound: E[s_k] >= P(s_k = 1), exact 95% interval of P(s_k = 1) over the K seed vectors
                cp95_prob_seed_all_won=cp(int((s >= 1 - 1e-12).sum()), ns))


def wtl(O):
    return int((O == 1).sum()), int((O == 0.5).sum()), int((O == 0).sum())


def cp(k, n):
    if n == 0:
        return [None, None]
    c = binomtest(int(k), int(n)).proportion_ci(0.95, method="exact")
    return [float(c.low), float(c.high)]


# ------------------------------------------------------------------ 1. benchmark
MAIN7, CTRL5 = RI.MAIN7, RI.CTRL5
BENCH_PAIRS = RI.B1_PAIRS


def group_inference(O, cases):
    groups = sorted({(c[0], int(c[1])) for c in cases})
    gid = np.array([groups.index((c[0], int(c[1]))) for c in cases])
    G = len(groups)
    contrib = np.array([np.rint(2 * (O[gid == g] - 0.5).sum()).astype(np.int64) for g in range(G)])   # 2 sum(O - 1/2)
    gs = np.array([O[gid == g].mean() for g in range(G)])
    signs = np.array(list(itertools.product([1, -1], repeat=G)))
    obs = abs(contrib.sum())
    p = 1.0 if obs == 0 else float(np.mean(np.abs(signs @ contrib) >= obs))
    est = gs.mean(); se = gs.std(ddof=1) / math.sqrt(G); q = t_dist.ppf(0.975, G - 1)
    return dict(groups=[f"{d}-{r}" for d, r in groups], n_cases=[int((gid == g).sum()) for g in range(G)],
                group_scores=[float(v) for v in gs], groups_above=int((gs > 0.5 + 1e-12).sum()),
                groups_below=int((gs < 0.5 - 1e-12).sum()), case_weighted_score=float(O.mean()),
                p_signflip_groups=p, p_min_attainable=float(2 / 2 ** G), equal_group_score=float(est),
                equal_group_t5_ci95=[float(est - q * se), float(est + q * se)])


def bench_block():
    A = X.load(HERE)
    G = A[A.Algorithm.isin(set(X.MAIN8) | set(X.ABL))].copy()
    S = X.case_stats(G, sorted(set(X.MAIN8) | set(X.ABL)), X.BASE_THR)
    cases = [tuple(c) for c in S[X.CASE].drop_duplicates().sort_values(X.CASE).itertuples(index=False)]
    assert len(cases) == 68
    G["V"] = 0.0                                                    # violation ordering not used (main variant)
    RM = {a: RI.run_matrix(G, a, cases) for a in sorted({x for p in BENCH_PAIRS for x in p})}
    REF = json.load(open(os.path.join(HERE, "rev3_inference.json")))["B1"]["pairs"]
    res = {}
    for j, (a, b) in enumerate(BENCH_PAIRS):
        O = RI.outcome(RM[a], RM[b])
        W, T, L = wtl(O)
        ref = REF[RI.pk(a, b)]["main"]
        assert (W, T, L) == (ref["W"], ref["T"], ref["L"]), (a, b, (W, T, L), ref)
        r = dict(a=a, b=b, W=W, T=T, L=L, n=int(O.size), pooled_sign_p=ref["p_sign"],
                 pooled_sign_p_holm=ref.get("p_sign_holm_main7", ref.get("p_sign_holm_controls5")))
        r["seed"] = seed_inference(O, len(cases), BOOT_SEED + j)
        r["group"] = group_inference(O, cases)
        res[RI.pk(a, b)] = r
    for fam, pairs in (("main7", [RI.pk("PSOBV", b) for b in MAIN7]), ("controls5", [RI.pk(a, b) for a, b in CTRL5])):
        for key, sub, out in (("p_signflip", "seed", "p_signflip_holm"), ("p_wilcoxon", "seed", "p_wilcoxon_holm"),
                              ("p_signflip_groups", "group", "p_signflip_groups_holm")):
            for k, h in zip(pairs, X.holm([res[k][sub][key] for k in pairs])):
                res[k][sub][out] = float(h)
        for k in pairs:
            res[k]["family"] = fam
            res[k]["pooled_sign_p_holm"] = REF[k]["main"]["p_sign_holm_" + fam]
    for k, r in res.items():
        sd = r["seed"]
        r["pooled_significant"] = bool(r["pooled_sign_p_holm"] < ALPHA)
        r["seed_significant"] = bool(sd["p_signflip_holm"] < ALPHA)
        r["seed_ci_excludes_half"] = bool(sd["ci95_joint_seed"][0] > 0.5 or sd["ci95_joint_seed"][1] < 0.5)
        r["group_significant"] = bool(r["group"]["p_signflip_groups_holm"] < ALPHA)
        print(f"BENCH {k:14s} {r['W']}/{r['T']}/{r['L']} score {sd['score']:.3f} seedCI {np.round(sd['ci95_joint_seed'], 3)} "
              f"seeds +{sd['seeds_above']}/-{sd['seeds_below']} p_flip {sd['p_signflip']:.3g} (Holm {sd['p_signflip_holm']:.3g}) "
              f"p_W {sd['p_wilcoxon']:.3g} | pooled Holm {r['pooled_sign_p_holm']:.3g} | groups +{r['group']['groups_above']}/"
              f"-{r['group']['groups_below']} p_grp {r['group']['p_signflip_groups']:.3g}", flush=True)
    return res, G, S, cases, RM


# ------------------------------------------------------------------ 5. calibration of the joint-seed TOST
def joint_tost_counts(La, Fa, Lb, Fb, Bi, rng, chunk=2500):
    """Joint-seed percentile bootstrap of the benchmark estimand mean_c [mean feasible L_a - mean feasible L_b]:
    one seed-index vector for all cases (as mpce_inference_extra.seed_level, joint); counts-matrix implementation."""
    nc, ns = La.shape
    FLa, FLb = np.where(Fa, La, 0.0), np.where(Fb, Lb, 0.0)
    fa, fb = Fa.astype(float), Fb.astype(float)
    full_a = FLa.sum(1) / np.maximum(fa.sum(1), 1); full_b = FLb.sum(1) / np.maximum(fb.sum(1), 1)
    bm = np.empty(Bi); fbk = 0
    for s in range(0, Bi, chunk):
        k = min(chunk, Bi - s)
        I = rng.integers(0, ns, (k, ns))
        C = np.bincount((I + ns * np.arange(k)[:, None]).ravel(), minlength=k * ns).reshape(k, ns).astype(float)
        sa, na, sb, nb = C @ FLa.T, C @ fa.T, C @ FLb.T, C @ fb.T
        fbk += int((na == 0).sum() + (nb == 0).sum())
        ma = np.where(na > 0, sa / np.maximum(na, 1), full_a[None]); mb = np.where(nb > 0, sb / np.maximum(nb, 1), full_b[None])
        bm[s:s + k] = (ma - mb).mean(1)
    return bm, fbk


def calib_block(G, S, RM, cases, R, Bi):
    out = {}
    for j, (a, b) in enumerate(RI.PRIMARY):
        d = X.case_diffs(S, a, b)
        assert len(d) == len(cases)
        sl = X.seed_level(G, a, b, d)[0]
        obs = sl["joint"]["ci90"]
        Ra, Rb = RM[a], RM[b]
        La0, Fa, Lb, Fb = Ra["loss"], Ra["feas"], Rb["loss"], Rb["feas"]
        # implementation check on the observed data (different RNG stream from the paper's)
        bm, _ = joint_tost_counts(La0, Fa, Lb, Fb, X.BOOT_N, np.random.default_rng(CALIB_SEED - 1))
        chk = [float(v) for v in np.quantile(bm, [0.05, 0.95])]
        fm = lambda L, F: np.array([L[i][F[i]].mean() for i in range(len(L))])
        est = float((fm(La0, Fa) - fm(Lb, Fb)).mean())
        assert abs(est - sl["joint"]["mean_dloss_pp"]) < 1e-12
        for side, mt in (("minus", -M), ("plus", M)):
            rng = np.random.default_rng(CALIB_SEED + 10 * j + (0 if side == "minus" else 1))
            La = La0 + (mt - est)
            assert abs((fm(La, Fa) - fm(Lb, Fb)).mean() - mt) < 1e-12
            t0 = time.time(); LIM = np.empty((R, 2)); fb = 0
            eq = rl = rh = 0
            for rep in range(R):
                J = rng.integers(0, La.shape[1], La.shape[1])          # 30 complete seed vectors, with replacement
                bm, f_ = joint_tost_counts(La[:, J], Fa[:, J], Lb[:, J], Fb[:, J], Bi, rng)
                fb += f_
                lo, hi = np.quantile(bm, [0.05, 0.95]); LIM[rep] = lo, hi
                eq += bool(-M < lo and hi < M); rl += bool(lo > -M); rh += bool(hi < M)
            k = int((LIM[:, 0] >= obs[0]).sum()) if side == "minus" else int((LIM[:, 1] <= obs[1]).sum())
            bind = rl if side == "minus" else rh
            r = dict(reps=R, inner_resamples=Bi, true_mean=mt, rate_eq_ci=eq / R, ci95_eq_ci=cp(eq, R),
                     rate_binding=bind / R, ci95_binding=cp(bind, R), observed_joint_ci90=obs,
                     p_calibrated_binding=(k + 1) / (R + 1), fallback_resamples=fb, seconds=round(time.time() - t0, 1))
            out[f"{RI.pk(a, b)}|{side}"] = r
            print(f"CALIB {RI.pk(a, b)} {side}: size (equivalence decision) {100 * r['rate_eq_ci']:.2f}% "
                  f"{np.round(100 * np.array(r['ci95_eq_ci']), 2)}; one-sided {100 * r['rate_binding']:.2f}%; calibrated p "
                  f"{r['p_calibrated_binding']:.4f} [{r['seconds']} s]", flush=True)
        out[f"{RI.pk(a, b)}|implementation_check"] = dict(paper_joint_ci90=obs, counts_impl_ci90=chk,
                                                          estimate=est)
    return out


# ------------------------------------------------------------------ 2. constraint study
def constraint_block():
    REF = json.load(open(os.path.join(HERE, "rev3_constraint.json")))["variants"]
    seeds = list(RCA.RC.SEEDS)
    fi = RCA.METH.index(RCA.FOCUS)
    out = {}
    for vi, v in enumerate(RCA.VARS):
        d = RCA.load(v)
        F = RCA.cube(d, "Feasible", seeds).astype(bool)
        Ob = RCA.cube(d, "Objective", seeds).astype(float)
        L = RCA.cube(d, "LossPct", seeds).astype(float)
        res = {}; tests = []
        for mi, m in enumerate(RCA.OTHERS):
            j = RCA.METH.index(m)
            O = RCA.outcome(F[fi], F[j], Ob[fi], Ob[j], L[fi], L[j])
            W, T, Lo = wtl(O)
            ref = REF[v]["comparisons"][m]
            assert (W, T, Lo) == (ref["pooled"]["W"], ref["pooled"]["T"], ref["pooled"]["L"]), (v, m)
            r = dict(W=W, T=T, L=Lo, n=int(O.size), pooled_sign_p=ref["pooled"]["sign_p"],
                     pooled_sign_p_holm=ref["pooled"]["sign_p_holm"], paper_category=ref["pooled"]["category"])
            r["seed"] = seed_inference(O, O.shape[0], BOOT_SEED + 100 + 10 * vi + mi)
            r["cases"] = {}
            for k, c in enumerate(RCA.CASES):
                cl = RCA.CLAB[c]; w, t, l = wtl(O[k])
                assert (w, t, l) == (ref["cases"][cl]["W"], ref["cases"][cl]["T"], ref["cases"][cl]["L"])
                r["cases"][cl] = dict(W=w, T=t, L=l, score=float(O[k].mean()), sign_p=ref["cases"][cl]["sign_p"],
                                      sign_p_holm_paper=ref["cases"][cl]["sign_p_holm"])
                tests.append((m, cl, ref["cases"][cl]["sign_p"]))
            tests.append((m, "pooled", r["seed"]["p_signflip"]))
            res[m] = r
        adj = X.holm([t[2] for t in tests])
        assert len(tests) == 28
        for (m, cl, _), h in zip(tests, adj):
            if cl == "pooled":
                res[m]["seed"]["p_signflip_holm28"] = float(h)
            else:
                res[m]["cases"][cl]["sign_p_holm28_seedpooled"] = float(h)
        pooled_only = X.holm([res[m]["seed"]["p_signflip"] for m in RCA.OTHERS])
        for m, h in zip(RCA.OTHERS, pooled_only):
            r = res[m]; sd = r["seed"]
            sd["p_signflip_holm4"] = float(h)
            sd["p_wilcoxon_holm4"] = float(X.holm([res[x]["seed"]["p_wilcoxon"] for x in RCA.OTHERS])[RCA.OTHERS.index(m)])
            r["pooled_significant"] = bool(r["pooled_sign_p_holm"] < ALPHA)
            r["seed_significant"] = bool(sd["p_signflip_holm28"] < ALPHA)
            r["seed_category"] = ("PSO-VNS better" if sd["score"] > 0.5 else "PSO-VNS worse") if r["seed_significant"] else "ns"
            r["case_flags_changed"] = [cl for cl in r["cases"] if (r["cases"][cl]["sign_p_holm28_seedpooled"] < ALPHA)
                                       != (r["cases"][cl]["sign_p_holm_paper"] < ALPHA)]
            print(f"CONSTR {v:4s} {m:6s} {r['W']}/{r['T']}/{r['L']} score {sd['score']:.3f} seedCI {np.round(sd['ci95_joint_seed'], 3)} "
                  f"seeds +{sd['seeds_above']}/-{sd['seeds_below']} p_flip {sd['p_signflip']:.3g} Holm28 {sd['p_signflip_holm28']:.3g} "
                  f"| paper pooled Holm {r['pooled_sign_p_holm']:.3g} ({r['paper_category']}) -> {r['seed_category']}", flush=True)
        out[v] = res
    return out


# ------------------------------------------------------------------ 3. direct 1-deg study
def fine_mats(Gm, a):
    R = Gm[Gm.Algorithm == a].set_index(RF.CASE + ["Seed"]).sort_index()
    seeds = list(range(31, 61))
    obj = np.empty((len(RF.CASES), 30)); loss = np.empty_like(obj); fe = np.zeros_like(obj, bool)
    for i, c in enumerate(RF.CASES):
        g = R.loc[c]
        assert list(g.index) == seeds, (a, c)
        obj[i] = g.Objective.values; loss[i] = g.LossPct.values; fe[i] = g.Feasible.values.astype(bool)
    return dict(obj=obj, loss=loss, feas=fe)


def fine_outcome(A, B):
    O = np.full(A["feas"].shape, 0.5)
    O[A["feas"] & ~B["feas"]] = 1.0; O[~A["feas"] & B["feas"]] = 0.0
    both = A["feas"] & B["feas"]; d = A["obj"] - B["obj"]
    O[both & (d > RF.TIE_OBJ)] = 1.0; O[both & (d < -RF.TIE_OBJ)] = 0.0
    return O


def matched_loss(D, C, boot_seed, B=BOOT_B):
    """direct minus control, wake loss (%) under the 1-deg objective, jointly feasible seed pairs; per case mean,
    then mean over cases; joint-seed bootstrap and sign flip of the per-seed case-average difference."""
    both = D["feas"] & C["feas"]
    nc, ns = both.shape
    dd = np.where(both, D["loss"] - C["loss"], 0.0); bf = both.astype(float)
    cmD = np.array([D["loss"][i][both[i]].mean() for i in range(nc)]); cmC = np.array([C["loss"][i][both[i]].mean() for i in range(nc)])
    est = float((cmD - cmC).mean())
    full = dd.sum(1) / np.maximum(bf.sum(1), 1)
    rng = np.random.default_rng(boot_seed)
    I = rng.integers(0, ns, (B, ns))
    Cn = np.bincount((I + ns * np.arange(B)[:, None]).ravel(), minlength=B * ns).reshape(B, ns).astype(float)
    sd_, nd_ = Cn @ dd.T, Cn @ bf.T
    bm = np.where(nd_ > 0, sd_ / np.maximum(nd_, 1), full[None]).mean(1)
    # per seed: mean over the cases in which the pair is jointly feasible
    ps = np.array([dd[:, k][both[:, k]].mean() for k in range(ns)])
    p2, pg, pl = signflip_exact(ps, tol=1e-12 * np.abs(ps).sum())
    fdm = lambda X_: float(np.mean([X_["loss"][i][X_["feas"][i]].mean() for i in range(nc)]))
    return dict(jointly_feasible_pairs=int(both.sum()), n_cases=int(nc), loss_direct=float(cmD.mean()),
                loss_control=float(cmC.mean()), diff=est, ci95_joint_seed=[float(v) for v in np.quantile(bm, [0.025, 0.975])],
                ci90_joint_seed=[float(v) for v in np.quantile(bm, [0.05, 0.95])], per_seed_mean=float(ps.mean()),
                seeds_lower=int((ps < 0).sum()), seeds_higher=int((ps > 0).sum()), p_signflip=p2, p_wilcoxon=wil_p(ps),
                loss_direct_all_feasible=fdm(D), loss_control_all_feasible=fdm(C),
                feasible_direct=int(D["feas"].sum()), feasible_control=int(C["feas"].sum()), runs=int(D["feas"].size))


def fine_block():
    REF = json.load(open(os.path.join(HERE, "rev3_fine.json")))
    Dn, _ = RF.read_shards("rev3_fine", HERE)
    Dn = RF.in_cases(Dn[Dn.Rose == "1deg"])
    Dc, _ = RF.read_shards("rev3_fine15", HERE)
    Dc = RF.in_cases(Dc[Dc.Rose == "15deg"])
    arms = dict(direct=RF.arm(Dn, "Objective", "Ideal"), ctrl1=RF.arm(Dc, "ObjFine", "IdealFine"),
                ctrl15=RF.arm(Dc, "Objective", "Ideal"))
    MT = {n: {a: fine_mats(Gm, a) for a in RF.METHODS} for n, Gm in arms.items()}
    out = dict(allrun={}, matched={})
    others = [b for b in RF.METHODS if b != RF.FOCUS]
    for ai, (name, mats) in enumerate(MT.items()):
        res = {}
        for bi, b in enumerate(others):
            O = fine_outcome(mats[RF.FOCUS], mats[b])
            W, T, L = wtl(O)
            ref = REF["arms"][name]["allrun_vs_focus"][b]
            assert (W, T, L) == (ref["W"], ref["T"], ref["L"]), (name, b)
            res[b] = dict(W=W, T=T, L=L, n=int(O.size), pooled_sign_p=ref["sign_test_p"], pooled_sign_p_holm=ref["sign_test_p_holm4"],
                          case_boot_ci95=ref["ci95"], seed=seed_inference(O, O.shape[0], BOOT_SEED + 200 + 10 * ai + bi))
        for b, h in zip(others, X.holm([res[b]["seed"]["p_signflip"] for b in others])):
            res[b]["seed"]["p_signflip_holm4"] = float(h)
        for b, h in zip(others, X.holm([res[b]["seed"]["p_wilcoxon"] for b in others])):
            res[b]["seed"]["p_wilcoxon_holm4"] = float(h)
        for b in others:
            r = res[b]; sd = r["seed"]
            r["pooled_significant"] = bool(r["pooled_sign_p_holm"] < ALPHA); r["seed_significant"] = bool(sd["p_signflip_holm4"] < ALPHA)
            print(f"FINE {name:6s} PSO-VNS vs {b:6s} {r['W']}/{r['T']}/{r['L']} score {sd['score']:.3f} seedCI "
                  f"{np.round(sd['ci95_joint_seed'], 3)} seeds +{sd['seeds_above']}/-{sd['seeds_below']} p_flip {sd['p_signflip']:.3g} "
                  f"Holm4 {sd['p_signflip_holm4']:.3g} | pooled Holm {r['pooled_sign_p_holm']:.3g}", flush=True)
        out["allrun"][name] = res
    ref_m = REF["compare"]["direct_minus_ctrl1_seed_paired"]
    for ai, a in enumerate(RF.METHODS):
        D, C = MT["direct"][a], MT["ctrl1"][a]
        r = matched_loss(D, C, BOOT_SEED + 300 + ai)
        assert abs(r["diff"] - ref_m[a]["mean_pp"]) < 1e-9 and r["jointly_feasible_pairs"] == ref_m[a]["n_pairs"], a
        O = fine_outcome(D, C)
        W, T, L = wtl(O)
        r["allrun_direct_vs_control"] = dict(W=W, T=T, L=L, seed=seed_inference(O, O.shape[0], BOOT_SEED + 400 + ai))
        out["matched"][a] = r
        sd = r["allrun_direct_vs_control"]["seed"]
        print(f"MATCHED {a:6s} feas {r['feasible_direct']}/{r['feasible_control']} L direct {r['loss_direct']:.3f} control "
              f"{r['loss_control']:.3f} diff {r['diff']:+.3f} CI95 {np.round(r['ci95_joint_seed'], 3)} seeds lower {r['seeds_lower']}/30 "
              f"p_flip {r['p_signflip']:.3g} | all-run {W}/{T}/{L} score {sd['score']:.3f} {np.round(sd['ci95_joint_seed'], 3)}", flush=True)
    ps = X.holm([out["matched"][a]["p_signflip"] for a in RF.METHODS])
    for a, h in zip(RF.METHODS, ps):
        out["matched"][a]["p_signflip_holm5"] = float(h)
    # PSO-VNS - PSO conditional contrast, joint-seed version (estimand of the lower panel of tab:S-sa-fine)
    out["psovns_pso_joint"] = {}
    import mpce_results as MR
    for name in ("direct", "ctrl1"):
        Gm = arms[name]
        S = MR.case_stats(Gm, RF.METHODS)
        d = X.case_diffs(S, "PSOBV", "PSOC")
        sl = X.seed_level(Gm, "PSOBV", "PSOC", d)[0]
        out["psovns_pso_joint"][name] = dict(n_cases=sl["n_cases"], mean=sl["mean_dloss_pp"], indep_ci90=sl["ci90"],
                                             joint_ci90=sl["joint"]["ci90"], joint_ci95=sl["joint"]["ci95"],
                                             joint_equivalent=sl["joint"]["equivalent"], joint_min_margin=sl["joint"]["min_margin_pp"])
        print(f"FINE PSO-VNS - PSO {name}: {sl['mean_dloss_pp']:+.4f} indep90 {np.round(sl['ci90'], 4)} joint90 "
              f"{np.round(sl['joint']['ci90'], 4)}", flush=True)
    return out


# ------------------------------------------------------------------ 4. Lillgrund boundary-aware intervals
def lillgrund_block():
    rows = []
    for fam in ("lg16", "lg16b"):
        files = sorted(glob.glob(os.path.join(HERE, f"rev2_{fam}_s*of*.csv")))
        df = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)
        df["Feasible"] = df.Feasible.astype(str).str.lower().isin(["true", "1"])
        budget = int(df.Budget.iloc[0])
        focus = df[df.Algorithm == "PSOBV"].set_index("Seed")
        for alg in sorted(set(df.Algorithm) - {"PSOBV"}):
            other = df[df.Algorithm == alg].set_index("Seed")
            sc = []
            for s in sorted(focus.index):
                a, b = focus.loc[s], other.loc[s]
                if a.Feasible != b.Feasible:
                    sc.append(float(a.Feasible))
                elif not a.Feasible:
                    sc.append(0.5)
                else:
                    sc.append(0.5 if abs(a.Objective - b.Objective) <= 1e-9 else float(a.Objective > b.Objective))
            sc = np.array(sc); W, T, L = int((sc == 1).sum()), int((sc == 0.5).sum()), int((sc == 0).sum())
            rows.append(dict(budget=budget, second=alg, label=LAB.get(alg, alg), W=W, T=T, L=L, score=float(sc.mean()),
                             win_prob_conditional=(W / (W + L)) if W + L else None, cp95_conditional=cp(W, W + L),
                             score_equals_conditional=bool(T == 0)))
    ref = pd.read_csv(os.path.join(VALID, "paired_outcomes.csv"))
    for r in rows:
        q = ref[(ref.budget == r["budget"]) & (ref.second == r["second"])].iloc[0]
        assert (r["W"], r["T"], r["L"]) == (q.wins, q.ties, q.losses), r
        r["bootstrap95_paper"] = [float(q.bootstrap95_lower), float(q.bootstrap95_upper)]
        r["holm16_p_paper"] = float(q.holm16_p)
        r["bootstrap_degenerate"] = bool(r["L"] == 0 and r["T"] == 0)
    for r in rows:
        print(f"LG {r['budget']:5d} {r['label']:8s} {r['W']}/{r['T']}/{r['L']} score {r['score']:.3f} W/(W+L) "
              f"{r['win_prob_conditional']:.3f} CP95 [{r['cp95_conditional'][0]:.4f}, {r['cp95_conditional'][1]:.4f}]", flush=True)
    return rows


# ------------------------------------------------------------------ LaTeX
def f3(v):
    return f"{v:.3f}".replace("-", "$-$")


def fci(c, d=3):
    return f"[{c[0]:.{d}f}, {c[1]:.{d}f}]".replace("-", "$-$")


def ptex(p):
    if p >= 0.001:
        return f"{p:.3f}"
    m_, e = f"{p:.1e}".split("e")
    return f"${m_}\\times10^{{{int(e)}}}$"


def plab(a, b):
    return f"{LAB[a]} vs.\\ {LAB[b]}"


def tables(R):
    out = ["% Generated by analysis/rev3_dependence.py. Do not edit by hand.", ""]
    BN = R["benchmark"]["pairs"]
    L = []
    for k in [RI.pk("PSOBV", b) for b in MAIN7] + ["mid"] + [RI.pk(a, b) for a, b in CTRL5]:
        if k == "mid":
            L.append("\\midrule\n\\multicolumn{9}{@{}l}{\\emph{Sampling-control contrasts (Holm over 5)}} \\\\")
            continue
        r = BN[k]; s = r["seed"]; g = r["group"]
        L.append(f"{plab(r['a'], r['b'])} & {r['W']}/{r['T']}/{r['L']} & {s['score']:.3f} & {fci(s['ci95_joint_seed'])} & "
                 f"{s['seeds_above']}/{s['seeds_below']} & {ptex(s['p_signflip_holm'])} & {ptex(s['p_wilcoxon_holm'])} & "
                 f"{g['groups_above']}/{g['groups_below']} & {ptex(g['p_signflip_groups_holm'])} \\\\")
    cal = R.get("calibration")
    L2 = []
    if cal:
        for a, b in RI.PRIMARY:
            c0, c1 = cal[f"{RI.pk(a, b)}|minus"], cal[f"{RI.pk(a, b)}|plus"]
            cell = lambda c: f"{100 * c['rate_eq_ci']:.2f} [{100 * c['ci95_eq_ci'][0]:.2f}, {100 * c['ci95_eq_ci'][1]:.2f}]"
            pc = lambda p: "$<$0.001" if p < 0.001 else f"{p:.3f}"
            th = lambda v: f"{v:,}".replace(",", "{,}")
            L2.append(f"{plab(a, b)} & {fci(c0['observed_joint_ci90'])} & {th(c0['reps'])} & {th(c0['inner_resamples'])} & "
                      f"{cell(c0)} & {cell(c1)} & {pc(c0['p_calibrated_binding'])} / {pc(c1['p_calibrated_binding'])} \\\\")
    cap = ("Seed-level (fixed-benchmark) inference for the all-run paired outcome of Table~\\ref{tab:S-sa-allrun} (68 cases, "
           "6,030 evaluations, random starts, seeds 1--30; feasible beats infeasible, two feasible runs by the objective, two "
           "infeasible runs tie). Seed $k$ seeds the same random stream in all 68 cases, so the 30 seed vectors, not the "
           "2,040 seed pairs, are the independent replicates. Seed score $s_k$: outcome mean of seed $k$ over the 68 cases; "
           "Score $=(W+T/2)/n=\\overline{s}$. 95\\% CI: percentile interval from 10{,}000 joint resamples of the 30 complete "
           "seed vectors (fixed seed). Seeds $+/-$: seeds with $s_k$ above/below 0.5. $p_{\\rm flip}$: exact sign-flip "
           "randomization test of $\\overline{s}-0.5$ over all $2^{30}$ sign patterns, valid if, under $H_0$, exchanging "
           "the two methods within a seed leaves the distribution of that seed's run pairs unchanged ($s_k-0.5$ symmetric "
           "about 0, independent over seeds); $p_{\\rm W}$: Wilcoxon signed-rank test of the 30 seed scores against 0.5. "
           "Groups: the six designed (data set, radius) groups; $+/-$ groups with score above/below 0.5, $p_{\\rm grp}$ exact "
           "sign-flip test of the case-weighted score over the $2^6=64$ group sign patterns (smallest attainable two-sided "
           "$p=1/32$; a sensitivity with six units only). All $p$ Holm-adjusted over the 7 comparisons of PSO-VNS or over the 5 "
           "contrasts.")
    if cal:
        cap += (" Bottom: empirical size (\\%, exact 95\\% binomial interval; nominal at most 5\\%) of the equivalence "
                "decision of the joint-seed percentile-bootstrap TOST of Table~\\ref{tab:S-sa-syncseed} (90\\% interval "
                "inside $(-m,m)$, $m=0.05$~pp) when the benchmark mean difference is exactly $-m$ or $+m$ (first method's run "
                "losses shifted by a constant); each replication draws 30 complete seed vectors with replacement and applies "
                "the joint bootstrap with the stated number of inner resamples. Calibrated $p_-$ ($p_+$): share of the "
                "replications at $-m$ ($+m$) whose lower (upper) 90\\% limit is at least as far inside the margin as the "
                "observed one (add-one).")
    out += ["\\begin{table}[!htbp]", "\\centering", "\\caption{" + cap + "}", "\\label{tab:S-sa-dep-bench}",
            "\\scriptsize\\setlength{\\tabcolsep}{2pt}", "\\begin{tabular}{@{}lcccccccc@{}}", "\\toprule",
            "& & & \\multicolumn{4}{c}{Seed level (30 seed vectors)} & \\multicolumn{2}{c}{Six groups} \\\\",
            "\\cmidrule(lr){4-7}\\cmidrule(l){8-9}",
            "Pair (first vs.\\ second) & $W/T/L$ & Score & 95\\% CI & Seeds $+/-$ & $p_{\\rm flip}$ & $p_{\\rm W}$ & $+/-$ & $p_{\\rm grp}$ \\\\",
            "\\midrule", "\\multicolumn{9}{@{}l}{\\emph{PSO-VNS vs.\\ each method (Holm over 7)}} \\\\"] + L + ["\\bottomrule", "\\end{tabular}"]
    if cal:
        out += ["\\par\\medskip", "\\begin{tabular}{@{}lcccccc@{}}", "\\toprule",
                "Pair & Observed 90\\% CI & Repl. & Inner $B$ & Size at $-m$ (\\%) & Size at $+m$ (\\%) & Calibrated $p_-$ / $p_+$ \\\\",
                "\\midrule"] + L2 + ["\\bottomrule", "\\end{tabular}"]
    out += ["\\end{table}", ""]
    # ---- studies
    C = R["constraint"]; FN = R["fine"]
    L = []
    for v in RCA.VARS:
        for i, m in enumerate(RCA.OTHERS):
            r = C[v][m]; s = r["seed"]
            L.append(f"{RCA.VSHORT[v] if i == 0 else ''} & {LAB[m]} & {r['W']}/{r['T']}/{r['L']} & {s['score']:.3f} & "
                     f"{fci(s['ci95_joint_seed'])} & {s['seeds_above']}/{s['seeds_below']} & {ptex(s['p_signflip_holm28'])} & "
                     f"{ptex(r['pooled_sign_p_holm'])} \\\\")
        L.append("\\midrule")
    L.append("\\multicolumn{8}{@{}l}{\\emph{Direct 1$^\\circ$ optimization, 8 cases, seeds 31--60 (Holm over 4)}} \\\\")
    for i, b in enumerate([x for x in RF.METHODS if x != RF.FOCUS]):
        r = FN["allrun"]["direct"][b]; s = r["seed"]
        L.append(f"1$^\\circ$ & {LAB[b]} & {r['W']}/{r['T']}/{r['L']} & {s['score']:.3f} & {fci(s['ci95_joint_seed'])} & "
                 f"{s['seeds_above']}/{s['seeds_below']} & {ptex(s['p_signflip_holm4'])} & {ptex(r['pooled_sign_p_holm'])} \\\\")
    L3 = []
    for a in RF.METHODS:
        r = FN["matched"][a]; ar = r["allrun_direct_vs_control"]; s = ar["seed"]
        L3.append(f"{LAB[a]} & {r['feasible_direct']}/{r['feasible_control']} & {r['jointly_feasible_pairs']} & {f3(r['loss_direct'])} & "
                  f"{f3(r['loss_control'])} & {f3(r['diff'])} & {fci(r['ci95_joint_seed'])} & {r['seeds_lower']}/{30 - r['seeds_lower']} & "
                  f"{ptex(r['p_signflip_holm5'])} & {ar['W']}/{ar['T']}/{ar['L']} & {s['score']:.3f} " +
                  (f"[{s['cp95_prob_seed_all_won'][0]:.2f}, 1]$^{{\\dagger}}$" if s["ci_degenerate"] else fci(s['ci95_joint_seed'], 2)) + " \\\\")
    cap = ("Seed-level inference for the pooled all-run tests of the constraint-handling study (Table~\\ref{tab:S-sa-constraint}; "
           "variants (i) penalty with box clipping, (ii) Deb's rules with box clipping, (iii) Deb's rules with radial projection; "
           "six cases, seeds 31--60, 180 seed pairs per comparison) and of the direct 1$^\\circ$ study "
           "(Table~\\ref{tab:S-sa-fine}; eight cases, seeds 31--60, 240 seed pairs), all-run score of PSO-VNS against the "
           "method. Seed $k$ seeds every case of a study, so the 30 seed vectors are the independent replicates; definitions "
           "as in Table~\\ref{tab:S-sa-dep-bench} (95\\% CI from 10{,}000 joint resamples of the seed vectors; $p_{\\rm flip}$ "
           "exact over all $2^{30}$ sign patterns). Constraint study: $p_{\\rm flip}$ Holm-adjusted in the family of 28 tests "
           "per variant of Table~\\ref{tab:S-sa-constraint} (4 methods $\\times$ (6 per-case exact sign tests and the pooled "
           "test, here the seed-level test)); 1$^\\circ$ study: Holm over the 4 comparisons. $p_{\\rm sign}$: the exact sign "
           "test over all discordant seed pairs (pairs treated as independent) with the Holm adjustment of the source table. "
           "Bottom: matched comparison of the direct 1$^\\circ$ arm with the 15$^\\circ$ control arm (same methods, cases and "
           "seeds 31--60), wake loss $L$ (\\%) under the 1$^\\circ$ objective. Feas.: feasible runs of 240 (direct/control); "
           "$n$: seed pairs with both runs feasible; $\\bar L$: mean over the 8 cases of the case means over these pairs; "
           "$\\Delta\\bar L$ direct minus control (pp; negative = direct better) with 95\\% joint-seed interval; Seeds $-/+$: "
           "seeds whose case-average difference is negative/positive; $p_{\\rm flip}$ exact sign-flip test of the 30 per-seed "
           "differences, Holm over the 5 methods. All-run: direct run vs.\\ control run of the same method and seed "
           "(feasible beats infeasible, two feasible runs by the 1$^\\circ$ objective), score with 95\\% joint-seed interval; "
           "$^{\\dagger}$all pairs won, so the percentile interval is degenerate ([1, 1]); given instead is the exact "
           "(Clopper--Pearson) 95\\% interval for the probability that a seed vector is won in all 8 cases (30 of 30 "
           "seeds), whose lower limit bounds the score from below.")
    out += ["\\begin{table}[!htbp]", "\\centering", "\\caption{" + cap + "}", "\\label{tab:S-sa-dep-studies}",
            "\\scriptsize\\setlength{\\tabcolsep}{2.5pt}", "\\begin{tabular}{@{}llcccccc@{}}", "\\toprule",
            "Study & Method & $W/T/L$ & Score & 95\\% CI & Seeds $+/-$ & $p_{\\rm flip}$ & $p_{\\rm sign}$ \\\\", "\\midrule",
            "\\multicolumn{8}{@{}l}{\\emph{Constraint handling, variants (i)--(iii) (Holm over 28 per variant)}} \\\\"]
    out += L + ["\\bottomrule", "\\end{tabular}", "\\par\\medskip", "\\setlength{\\tabcolsep}{2pt}",
                "\\begin{tabular}{@{}lcccccccccc@{}}", "\\toprule",
                "& & & \\multicolumn{6}{c}{Wake loss, jointly feasible seed pairs} & \\multicolumn{2}{c}{All-run, direct vs.\\ control} \\\\",
                "\\cmidrule(lr){4-9}\\cmidrule(l){10-11}",
                "Method & Feas. & $n$ & $\\bar L$ direct & $\\bar L$ control & $\\Delta\\bar L$ & 95\\% CI & Seeds $-/+$ & $p_{\\rm flip}$ & $W/T/L$ & Score [95\\% CI] \\\\",
                "\\midrule"] + L3 + ["\\bottomrule", "\\end{tabular}", "\\end{table}", ""]
    return "\n".join(out)


def clean(o):
    if isinstance(o, dict):
        return {str(k): clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [clean(v) for v in o]
    if isinstance(o, np.bool_):
        return bool(o)
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, np.floating):
        return float(o)
    return o


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", default=HERE)
    ap.add_argument("--calib-reps", type=int, default=5000)
    ap.add_argument("--calib-inner", type=int, default=10000)
    ap.add_argument("--skip-calib", action="store_true")
    ap.add_argument("--tables-only", action="store_true", help="only rewrite the .tex from rev3_dependence.json")
    args = ap.parse_args(argv)
    if args.tables_only:                                # rewrite the LaTeX tables from an existing JSON
        R = json.load(open(os.path.join(args.out_dir, "rev3_dependence.json")))
        open(os.path.join(args.out_dir, "rev3_dependence_tables.tex"), "w").write(tables(R))
        return
    t0 = time.time()
    # self-check of the exact sign-flip routine against full enumeration (K = 12)
    rng = np.random.default_rng(1)
    v = rng.integers(-5, 6, 12)
    tot = np.array([np.dot(s, v) for s in itertools.product([1, -1], repeat=12)])
    assert abs(signflip_exact(v)[0] - (np.abs(tot) >= abs(v.sum())).mean()) < 1e-15
    assert abs(signflip_exact(v)[1] - (tot >= v.sum()).mean()) < 1e-15
    R = dict(generated=time.strftime("%Y-%m-%d %H:%M:%S"), script="rev3_dependence.py", margin_pp=M,
             conventions=dict(score="(W + T/2)/n of the first method; seed score s_k = mean over the cases of seed k",
                              ci="95% percentile interval, complete seed vectors resampled jointly",
                              boot=dict(B=BOOT_B, seed=BOOT_SEED, note="seed + offset per comparison"),
                              test=("exact two-sided sign-flip randomization test of mean_k(s_k - 1/2), all 2^K patterns; "
                                    "H0: within a seed the two methods are exchangeable (label swap leaves the joint "
                                    "distribution of the seed's run pairs unchanged), seeds independent"),
                              check="Wilcoxon signed-rank of the K seed scores against 1/2 (zeros dropped)",
                              alpha=ALPHA))
    t1 = time.time()
    bench, G, S, cases, RM = bench_block()
    R["benchmark"] = dict(pairs=bench, n_cases=68, seeds="1-30", families=dict(main7=[RI.pk("PSOBV", b) for b in MAIN7],
                          controls5=[RI.pk(a, b) for a, b in CTRL5]), seconds=round(time.time() - t1, 1))
    if not args.skip_calib:
        t1 = time.time()
        R["calibration"] = calib_block(G, S, RM, cases, args.calib_reps, args.calib_inner)
        R["calibration_design"] = dict(reps=args.calib_reps, inner=args.calib_inner, seed=CALIB_SEED,
                                       seconds=round(time.time() - t1, 1))
    R["constraint"] = constraint_block()
    R["fine"] = fine_block()
    R["lillgrund"] = lillgrund_block()
    # summary of verdicts
    summ = []
    for k, r in bench.items():
        summ.append(dict(study="benchmark", comparison=k, family=r["family"], pooled_holm=r["pooled_sign_p_holm"],
                         seed_holm=r["seed"]["p_signflip_holm"], pooled_sig=r["pooled_significant"], seed_sig=r["seed_significant"],
                         group_holm=r["group"]["p_signflip_groups_holm"]))
    for v, res in R["constraint"].items():
        for m, r in res.items():
            summ.append(dict(study=f"constraint-{v}", comparison=f"PSOBV-{m}", pooled_holm=r["pooled_sign_p_holm"],
                             seed_holm=r["seed"]["p_signflip_holm28"], pooled_sig=r["pooled_significant"], seed_sig=r["seed_significant"],
                             case_flags_changed=r["case_flags_changed"]))
    for b, r in R["fine"]["allrun"]["direct"].items():
        summ.append(dict(study="fine-direct", comparison=f"PSOBV-{b}", pooled_holm=r["pooled_sign_p_holm"],
                         seed_holm=r["seed"]["p_signflip_holm4"], pooled_sig=r["pooled_significant"], seed_sig=r["seed_significant"]))
    R["verdicts"] = summ
    R["verdict_changes"] = [s for s in summ if s["pooled_sig"] != s["seed_sig"]]
    R["seconds"] = round(time.time() - t0, 1)
    json.dump(clean(R), open(os.path.join(args.out_dir, "rev3_dependence.json"), "w"), indent=1)
    open(os.path.join(args.out_dir, "rev3_dependence_tables.tex"), "w").write(tables(R))
    print("verdict changes (pooled sign test vs seed level, Holm 0.05):", json.dumps(clean(R["verdict_changes"])), flush=True)
    print(f"done in {R['seconds']} s", flush=True)


if __name__ == "__main__":
    main()
