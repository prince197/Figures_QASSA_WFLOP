"""[rev3 audit copy of make_theory_figures.py; see the comment at FIG below]
Theory numbers, figure and checks of Section S-theory (supplement) and optA/drafts/theory_main.tex.

Usage (from analysis/):
    python make_theory_figures.py            compute everything, write ../figures_mpce/theory_stability.pdf, print
    python make_theory_figures.py --check    additionally verify that every number typed in optA/supp_theory.tex
                                             and optA/drafts/theory_main.tex equals the value computed here
                                             (checks T01...; exit code 1 on any FAIL)
    python make_theory_figures.py --check --swevo   the same checks against the SWEVO methods section (optA/sw/04_methods.tex)
    python make_theory_figures.py --fast     smaller Monte Carlo (for a quick look; --check needs the full run)

Everything here is mathematics or geometry (no objective evaluations, no run data):
  T01-T06  PSO stability (Proposition S-pso-stability): order-1 / order-2 bounds, spectral radii, margins, the
           stationary variance factor, and a numerical confirmation that the analytic order-2 region equals
           {spectral radius of the second-moment matrix < 1} (random points, general c1 != c2, and a c1 = c2 grid).
  T07-T09  box sampling (Proposition S-box): (pi/4)^N, the rigorous upper bound with spacing, Monte Carlo estimates
           of the probability that a uniform layout in the bounding square is feasible (fixed seed).
  T10      the lens area is nonincreasing in the distance of the centres (used in the proof of the bound).
  T11-T12  DE crossover count and SSA leader radius; compass-search final step (Proposition S-ls).
Single process; runtime about 2-3 minutes (Monte Carlo 10^7 layouts per case).
"""
import os, re, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
# rev3 numbers audit copy (analysis/rev3_numbers_theory_check.py): identical computations and checks as
# make_theory_figures.py, but (1) the figure goes to a scratch directory (env REV3_FIG_DIR, default /tmp), never to
# figures_mpce/, and (2) the checked text is the revision-2 supplement theory and methods sections
# (SWEVO_rev2/latex_source/optA/sw/supp_theory.tex and optA/sw/04_methods.tex); with --repo-text the repository's
# SWEVO copies (optA/sw/supp_theory.tex, optA/sw/04_methods.tex) are checked instead, for comparison.
FIG = os.environ.get("REV3_FIG_DIR", "/tmp/rev3_numbers_theory_fig")
_R2 = os.path.join(ROOT, "SWEVO_rev2", "latex_source", "optA", "sw")
TEX = [os.path.join(_R2, "supp_theory.tex"), os.path.join(_R2, "04_methods.tex")]
if "--repo-text" in sys.argv:
    TEX = [os.path.join(ROOT, "optA", "sw", "supp_theory.tex"), os.path.join(ROOT, "optA", "sw", "04_methods.tex")]
print("rev3 theory check; text files:", TEX, "; figure dir:", FIG)
FAST = "--fast" in sys.argv
SEED = 20260928

INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e4e3df"
BLUE, BLUE_XL = "#2a78d6", "#e6f0fc"
ORANGE, GREEN = "#eb6834", "#008300"
plt.rcParams.update({"font.family": "STIXGeneral", "mathtext.fontset": "stix", "font.size": 8.5,
                     "text.color": INK, "axes.labelcolor": INK, "xtick.color": MUTED, "ytick.color": MUTED,
                     "pdf.fonttype": 42})

SETTINGS = {"old": (0.7, 2.0, 2.0), "constriction": (0.7298, 1.49618, 1.49618)}


# ----------------------------------------------------------------------------------------------------------
# PSO stability under stagnation, no bounds: x_{t+1} = (1+w-phi_t) x_t - w x_{t-1} + phi1 p + phi2 g,
# phi_k = c_k U(0,1) independent for every coordinate and move (authors_optimizers.PSO).
def order1_bound(w):
    return 4 * (1 + w)                                     # c1 + c2 < 4(1+w)  (E phi = (c1+c2)/2 < 2(1+w))


def order2_bound(w):
    return 24 * (1 - w * w) / (7 - 5 * w)                  # c1 = c2: c1 + c2 < 24(1-w^2)/(7-5w)  [Poli 2009]


def M2(w, c1, c2):
    """Second-moment matrix acting on (E y_t^2, E y_t y_{t-1}, E y_{t-1}^2)."""
    mu = (c1 + c2) / 2; s2 = (c1 * c1 + c2 * c2) / 12; m = 1 + w - mu; a = m * m + s2
    return np.array([[a, -2 * w * m, w * w], [m, -w, 0.0], [1.0, 0.0, 0.0]])


def rho1(w, c1, c2):
    return max(abs(np.roots([1.0, -(1 + w - (c1 + c2) / 2), w])))


def rho2(w, c1, c2):
    return max(abs(np.linalg.eigvals(M2(w, c1, c2))))


def chi1(w, c1, c2):
    """chi(1) = det(I - M): order-2 stable iff |w| < 1 and chi(1) > 0 (Proposition S-pso-stability)."""
    mu = (c1 + c2) / 2; s2 = (c1 * c1 + c2 * c2) / 12
    return 2 * mu * (1 - w * w) - (1 - w) * mu * mu - (1 + w) * s2


def stat_var_factor(w, c):
    """c1 = c2 = c: lim Var x_t / (p - g)^2 = c(1+w) / (4[12(1-w^2) - c(7-5w)])  [Poli 2009]."""
    return c * (1 + w) / (4 * (12 * (1 - w * w) - c * (7 - 5 * w)))


def stability_numbers():
    out = {}
    for k, (w, c1, c2) in SETTINGS.items():
        s = c1 + c2
        out[k] = dict(w=w, sum=s, b1=order1_bound(w), b2=order2_bound(w), rho1=rho1(w, c1, c2),
                      rho2=rho2(w, c1, c2), margin2=order2_bound(w) - s, rel2=s / order2_bound(w),
                      var=stat_var_factor(w, c1))
    # analytic region == numerical region (points closer than 1e-3 to the boundary are skipped)
    rng = np.random.default_rng(SEED)
    bad = n = 0
    for _ in range(20000 if FAST else 200000):
        w = rng.uniform(-0.999, 0.999); c1 = rng.uniform(0, 6); c2 = rng.uniform(0, 6)
        g = chi1(w, c1, c2)
        if abs(g) < 1e-3:
            continue
        n += 1; bad += (rho2(w, c1, c2) < 1) != (g > 0)
    out["region_random"] = (bad, n)
    bad = n = 0
    for w in np.linspace(-0.999, 0.999, 201 if FAST else 601):
        for s in np.linspace(0.005, 8, 201 if FAST else 601):
            if abs(s - order2_bound(w)) < 1e-3:
                continue
            n += 1; bad += (rho2(w, s / 2, s / 2) < 1) != (s < order2_bound(w))
    out["region_grid"] = (bad, n)
    return out


# ----------------------------------------------------------------------------------------------------------
# Box sampling: N points uniform in [-r, r]^2; all inside the circle and pairwise >= s apart.
def lens(r, rho, d=None):
    """Area of B(c, rho) intersected with B(0, r), |c| = d (default d = r); rho <= 2r."""
    d = r if d is None else d
    if d <= abs(r - rho):
        return np.pi * min(r, rho) ** 2
    if d >= r + rho:
        return 0.0
    k = (-d + r + rho) * (d + r - rho) * (d - r + rho) * (d + r + rho)
    return (r * r * np.arccos((d * d + r * r - rho * rho) / (2 * d * r))
            + rho * rho * np.arccos((d * d + rho * rho - r * r) / (2 * d * rho)) - 0.5 * np.sqrt(k))


def spacing_bound(N, r, s):
    """Upper bound on P(pairwise spacing >= s | N iid uniform points in the disc of radius r)."""
    A = np.pi * r * r; a1 = lens(r, s); a2 = lens(r, s / 2); P = 1.0
    for k in range(1, N):
        P *= max(0.0, 1 - max(a1, k * a2) / A)
    return P


def mc_spacing(N, r, s, n, seed=SEED, chunk=500000):
    """Monte Carlo: number of the n iid uniform-in-disc layouts that satisfy the spacing (sequential rejection)."""
    rng = np.random.default_rng(seed); hits = 0
    for start in range(0, n, chunk):
        m = min(chunk, n - start)
        ang = rng.uniform(0, 2 * np.pi, (m, N)); rad = r * np.sqrt(rng.random((m, N)))
        x = rad * np.cos(ang); y = rad * np.sin(ang)
        alive = np.arange(m)
        for k in range(1, N):
            ok = np.ones(len(alive), bool)
            for j in range(k):
                ok &= (x[alive, k] - x[alive, j]) ** 2 + (y[alive, k] - y[alive, j]) ** 2 >= s * s
            alive = alive[ok]
            if len(alive) == 0:
                break
        hits += len(alive)
    return hits


def clopper_pearson(k, n, level=0.95):
    from scipy.stats import beta
    a = 1 - level
    lo = 0.0 if k == 0 else beta.ppf(a / 2, k, n - k + 1)
    hi = beta.ppf(1 - a / 2, k + 1, n - k) if k < n else 1.0
    return lo, hi


# rows of Table S-th-box: (N, r, s, label)
BOX_ROWS = [(2, 500, 308.0, "benchmark, smallest case"), (10, 500, 308.0, "benchmark"), (12, 750, 308.0, "benchmark"),
            (15, 1000, 308.0, "benchmark"), (16, 1300, 260.0, "IEA37"), (20, 1000, 308.0, "hypothetical"),
            (36, 2000, 260.0, "IEA37")]
RS_DRAWS = 3015               # uniform layouts of the RS-VNS Phase 1 at B = 6,030 (incl. the initial population)


def box_numbers():
    rows = []
    n = 10 ** 6 if FAST else 10 ** 7
    for N, r, s, lab in BOX_ROWS:
        q = (np.pi / 4) ** N
        ub = q * spacing_bound(N, r, s)
        k = mc_spacing(N, r, s, n)
        lo, hi = clopper_pearson(k, n)
        rows.append(dict(N=N, r=r, s=s, lab=lab, inside=q, bound=ub, hits=k, n=n, est=q * k / n,
                         lo=q * lo, hi=q * hi, exp_rs=RS_DRAWS * q * k / n, exp_rs_hi=RS_DRAWS * q * hi))
    # N = 2 exact value (distance density of two uniform points in a disc) as a check of the Monte Carlo
    r, s = 500.0, 308.0
    t = np.linspace(0, s, 200001)
    f = 4 * t / (np.pi * r * r) * (np.arccos(t / (2 * r)) - t / (2 * r) * np.sqrt(1 - t * t / (4 * r * r)))
    exact2 = (np.pi / 4) ** 2 * (1 - np.trapezoid(f, t) if hasattr(np, "trapezoid") else 1 - np.trapz(f, t))
    # lens area nonincreasing in the distance of the centres (T10)
    mono = all(lens(r0, rho, d1) >= lens(r0, rho, d2) - 1e-9
               for r0 in (500., 1000.) for rho in (154., 308.)
               for d1, d2 in zip(np.linspace(1, r0, 400)[:-1], np.linspace(1, r0, 400)[1:]))
    return rows, exact2, mono


def fmt_sci(x, digits=2):
    """LaTeX scientific notation with `digits` significant digits (0 -> '0')."""
    if x == 0:
        return "0"
    e = int(np.floor(np.log10(abs(x))))
    m = x / 10 ** e
    if round(m, digits - 1) >= 10:
        m /= 10; e += 1
    if -2 <= e <= 0:
        return f"{x:.{digits - 1 - e}f}"
    return f"{m:.{digits - 1}f}\\cdot10^{{{e}}}"


def box_rows_tex(rows):
    out = []
    for r in rows:
        est = (f"${fmt_sci(r['est'])}$ $[{fmt_sci(r['lo'])},\\,{fmt_sci(r['hi'])}]$" if r["hits"] else
               f"$0$ $[0,\\,{fmt_sci(r['hi'])}]$")
        ex = f"{r['exp_rs']:,.0f}".replace(",", "{,}") if r["exp_rs"] >= 1 else (
            f"${fmt_sci(r['exp_rs'])}$" if r["hits"] else f"$<{fmt_sci(r['exp_rs_hi'])}$")
        rr = f"{r['r']:,}".replace(",", "{,}")
        out.append(f"{r['N']} & {rr} & {r['s']:.0f} & ${fmt_sci(r['inside'])}$ & ${fmt_sci(r['bound'])}$ & {est} & {ex} \\\\")
    return out


# ----------------------------------------------------------------------------------------------------------
def figure(st):
    fig, ax = plt.subplots(figsize=(3.45, 2.45))
    w = np.linspace(-1, 1, 801)
    ax.fill_between(w, 0, order1_bound(w), color=BLUE_XL, lw=0, zorder=1)
    ax.fill_between(w, 0, order2_bound(w), color=BLUE, alpha=0.30, lw=0, zorder=2)
    ax.plot(w, order1_bound(w), color=BLUE, lw=1.0, ls=(0, (4, 2)), zorder=3)
    ax.plot(w, order2_bound(w), color=BLUE, lw=1.4, zorder=3)
    ax.text(0.68, 5.4, "order-1 only", ha="center", va="center", fontsize=7.5, color=INK)
    ax.text(-0.84, 1.22, r"$c_1+c_2=4(1+w)$", ha="left", va="bottom", fontsize=7.0, color=MUTED,
            rotation=np.degrees(np.arctan(4.0)), rotation_mode="anchor", transform_rotates_text=True)
    ax.text(-0.22, 1.05, "order-1 and order-2\n" r"$c_1+c_2<\frac{24(1-w^2)}{7-5w}$", ha="center", va="center",
            fontsize=7.5, color=INK)
    ax.text(-0.93, 7.55, "not order-1 stable", ha="left", va="center", fontsize=7.2, color=MUTED)
    o, c = st["old"], st["constriction"]
    ax.scatter([o["w"]], [o["sum"]], s=34, marker="s", color=ORANGE, edgecolor="white", lw=0.8, zorder=6)
    ax.scatter([c["w"]], [c["sum"]], s=40, marker="o", color=GREEN, edgecolor="white", lw=0.8, zorder=6)
    ax.annotate(f"old setting\n$w=0.7$, $c_1+c_2=4$\n(order-2 bound {o['b2']:.2f})", (o["w"], o["sum"]),
                xytext=(-0.33, 6.45), fontsize=7.2, ha="center", va="center", color=INK,
                arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.6, shrinkA=1, shrinkB=4))
    ax.annotate(f"constriction\n$w=0.7298$, $c_1+c_2=2.99$\n(order-2 bound {c['b2']:.2f})", (c["w"], c["sum"]),
                xytext=(0.55, 1.0), fontsize=7.2, ha="center", va="center", color=INK,
                arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.6, shrinkA=1, shrinkB=4))
    ax.set_xlim(-1, 1); ax.set_ylim(0, 8)
    ax.set_xlabel(r"inertia weight $w$"); ax.set_ylabel(r"$c_1+c_2$")
    ax.set_xticks([-1, -0.5, 0, 0.5, 1]); ax.set_yticks([0, 2, 4, 6, 8])
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    for sp in ("left", "bottom"):
        ax.spines[sp].set_color(MUTED); ax.spines[sp].set_linewidth(0.6)
    ax.tick_params(length=2.5, width=0.6)
    os.makedirs(FIG, exist_ok=True)
    fig.savefig(os.path.join(FIG, "theory_stability.pdf"), bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)


# ----------------------------------------------------------------------------------------------------------
def main():
    st = stability_numbers()
    rows, exact2, mono = box_numbers()
    figure(st)
    o, c = st["old"], st["constriction"]
    de = {n: 1 + 0.9 * (n - 1) for n in (20, 24, 30, 32)}
    xi_half = 2 * np.exp(-4.0)                     # SSA xi_1 at t = T/2
    hbar = 0.05 / 2 ** 5                           # compass search: last step h0/2^5 (in units of r)
    assert 0.05 / 2 ** 5 >= 1e-3 > 0.05 / 2 ** 6

    print("PSO stability (stagnation, no bounds):")
    for k in ("old", "constriction"):
        d = st[k]
        print(f"  {k:12s} w={d['w']} c1+c2={d['sum']:.5f}  order-1 bound {d['b1']:.4f}  order-2 bound {d['b2']:.4f}"
              f"  margin {d['margin2']:+.4f} (ratio {d['rel2']:.4f})  rho1 {d['rho1']:.4f}  rho2 {d['rho2']:.4f}"
              f"  var factor {d['var']:.4f}")
    print("  region check (random general c1,c2): mismatches %d of %d" % st["region_random"])
    print("  region check (c1 = c2 grid):         mismatches %d of %d" % st["region_grid"])
    print("Box sampling (uniform in the square):")
    for r in rows:
        print(f"  N={r['N']:2d} r={r['r']:5d} s={r['s']:.0f}  (pi/4)^N={r['inside']:.3e}  bound={r['bound']:.3e}"
              f"  MC hits {r['hits']}/{r['n']}  P_N={r['est']:.3e} [{r['lo']:.2e}, {r['hi']:.2e}]"
              f"  E[#feasible of {RS_DRAWS}]={r['exp_rs']:.2e} (<= {r['exp_rs_hi']:.2e})  ({r['lab']})")
    for line in box_rows_tex(rows):     # LaTeX rows of Table S-th-box (typed in the supplement; T13)
        print("    " + line)
    print(f"  N=2 exact P_2 (r=500, s=308) = {exact2:.4f};  lens monotone: {mono}")
    print("DE: expected number of mutant coordinates 1 + CR(n-1):", {k: round(v, 1) for k, v in de.items()})
    print(f"SSA: xi_1(T/2) = {xi_half:.4f};  compass search final step h = {hbar} r")

    if "--check" not in sys.argv:
        return 0
    if FAST:
        print("--check needs the full Monte Carlo (run without --fast)"); return 1
    txt = "\n".join(open(f, encoding="utf-8").read() for f in TEX if os.path.exists(f))
    txt_main = open(TEX[1], encoding="utf-8").read() if os.path.exists(TEX[1]) else ""
    checks = []

    def has(cid, desc, s, where=txt):
        checks.append((cid, desc + f"  [{s}]", s in where))
    has("T01", "old: order-1 bound 6.8 and order-2 bound 3.50", f"{o['b2']:.2f}")
    has("T01", "old: order-1 bound", f"{o['b1']:.1f}")
    has("T02", "constriction: c1+c2 and order-2 bound", f"{c['sum']:.2f}<{c['b2']:.2f}".replace("<", "<"))
    has("T02", "constriction: order-1 bound", f"{c['b1']:.2f}")
    has("T03", "spectral radii rho2 (old, constriction)", f"{o['rho2']:.2f}")
    has("T03", "spectral radius rho2 constriction", f"{c['rho2']:.2f}")
    has("T03", "order-1 radii sqrt(w)", f"{o['rho1']:.2f}")
    has("T03", "order-1 radius constriction", f"{c['rho1']:.2f}")
    has("T04", "margins", f"{c['margin2']:.2f}")
    has("T04", "old excess", f"{-o['margin2']:.2f}")
    has("T04", "constriction bound 4 digits", f"{c['b2']:.4f}")
    has("T04", "old bound 4 digits", f"{o['b2']:.4f}")
    has("T04", "constriction relative margin", f"{100 * (1 - c['rel2']):.1f}\\% below")
    has("T04", "old relative excess", f"{100 * (o['rel2'] - 1):.1f}\\% above")
    has("T05", "stationary variance factor and SD (constriction)", f"{c['var']:.2f}")
    has("T05", "stationary SD factor", f"{np.sqrt(c['var']):.2f}")
    checks.append(("T06", "analytic order-2 region equals numerical region (no mismatch)",
                   st["region_random"][0] == 0 and st["region_grid"][0] == 0))
    for r in rows:
        has("T07", f"(pi/4)^N N={r['N']}", fmt_sci(r["inside"]))
        has("T08", f"upper bound N={r['N']} r={r['r']}", fmt_sci(r["bound"]))
        if r["hits"] > 0:
            has("T09", f"MC estimate N={r['N']} r={r['r']}", fmt_sci(r["est"]))
        has("T09", f"MC upper CI N={r['N']} r={r['r']}", fmt_sci(r["hi"]))
    for line in box_rows_tex(rows):
        checks.append(("T13", f"table row verbatim in supp_theory.tex  [{line[:40]}...]", line in txt))
    has("T09", "N=2 exact value", f"{exact2:.3f}")
    checks.append(("T10", "lens area nonincreasing in the distance of the centres", mono))
    has("T11", "DE mutant coordinates n=30", f"{de[30]:.1f}")
    has("T11", "SSA xi_1(T/2)", f"{xi_half:.3f}")
    has("T12", "compass search final step", f"{hbar}")
    has("T12", "main text: order-2 bound old", f"{o['b2']:.2f}", txt_main)
    has("T12", "main text: order-2 bound constriction", f"{c['b2']:.2f}", txt_main)
    fails = 0
    for cid, desc, ok in checks:
        print(f"{cid} {'PASS' if ok else 'FAIL'}  {desc}")
        fails += not ok
    print(f"{len(checks) - fails} PASS, {fails} FAIL")
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
