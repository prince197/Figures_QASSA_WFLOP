"""Graphical abstract for the SWEVO submission (Elsevier: single panel, >= 531 x 1328 px, 5 x 13 cm).

Usage (from anywhere):  python3 analysis/make_graphical_abstract.py
Writes figures_final/graphical_abstract.{pdf,png,tif} and copies them to
SWEVO_rev2/latex_source/figures_final/.

Every number shown is read from the analysis JSON files and checked against the value printed in the paper:
  * PSO setting (Table 4, Section 7.1): mpce_summary.json /baseline/{old,constriction}/avg_rank
  * Sampling control (Tables 8/9, Section 8.3): mpce_summary.json /ablation/contrasts/SSABV-{RSVNS,RSDVNS}/dloss_pp
  * Constraint handling (Section 10.4, Table S79): rev3_constraint.json /variants/{pen,deb}/ranks/*/case_rule
  * Direction bins (Section 10.5, Table S81): rev3_fine.json /arms/{ctrl1,direct}/contrast_psovns_pso/mean_dloss_pp
    and /arms/{ctrl1,direct}/leader
The canvas is 130 x 52 mm (aspect 2.5 = 1328/531); all text is >= 9 pt at that size.
"""
import json
import os
import shutil

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(ROOT, "figures_final")
COPY_TO = os.path.join(ROOT, "SWEVO_rev2", "latex_source", "figures_final")


def load(name):
    with open(os.path.join(HERE, name)) as f:
        return json.load(f)


S = load("mpce_summary.json")
C = load("rev3_constraint.json")
F = load("rev3_fine.json")

# ---- numbers (with the value printed in the paper as a check) ----
old = S["baseline"]["old"]["avg_rank"]
con = S["baseline"]["constriction"]["avg_rank"]
pso_old, pso_con = old["PSO"], con["PSO"]            # 5.04 -> 1.53 (Table 4)
ssa_old, ssa_con = old["SSABV"], con["SSABV"]       # 1.99 -> 2.77 (Table 4)
pos_old = S["baseline"]["old"]["pso_position"]      # 5
pos_con = S["baseline"]["constriction"]["pso_position"]  # 1
ab = S["ablation"]["contrasts"]
gain_sq = -ab["SSABV-RSVNS"]["dloss_pp"]            # SSA-VNS - RS-VNS = -0.068 -> gain +0.068
gain_disc = -ab["SSABV-RSDVNS"]["dloss_pp"]         # SSA-VNS - RSD-VNS = +0.024 -> gain -0.024
rk = {v: C["variants"][v]["ranks"] for v in ("pen", "deb")}
pv_pen, ga_pen = rk["pen"]["PSOBV"]["case_rule"], rk["pen"]["GA"]["case_rule"]   # 1.67, 2.50
pv_deb, ga_deb = rk["deb"]["PSOBV"]["case_rule"], rk["deb"]["GA"]["case_rule"]   # 2.17, 1.83
adv15 = -F["arms"]["ctrl1"]["contrast_psovns_pso"]["mean_dloss_pp"]   # 0.086 (opt. 15 deg, eval. 1 deg)
adv1 = -F["arms"]["direct"]["contrast_psovns_pso"]["mean_dloss_pp"]   # 0.194 (opt. 1 deg)
lead15, lead1 = F["arms"]["ctrl1"]["leader"], F["arms"]["direct"]["leader"]
N_CASES = S["baseline"]["old"]["n_cases"]           # 68
MARGIN = 0.05                                        # post hoc equivalence margin (pp), Section 6.2

checks = [(pso_old, 5.04, 2), (pso_con, 1.53, 2), (ssa_old, 1.99, 2), (ssa_con, 2.77, 2),
          (gain_sq, 0.068, 3), (gain_disc, -0.024, 3), (pv_pen, 1.67, 2), (ga_pen, 2.50, 2),
          (pv_deb, 2.17, 2), (ga_deb, 1.83, 2), (adv15, 0.086, 3), (adv1, 0.194, 3)]
for val, paper, nd in checks:
    assert round(val, nd) == paper, (val, paper)
assert (pos_old, pos_con, N_CASES) == (5, 1, 68) and lead15 == lead1 == ["PSOBV"]

# ---- style ----
# method colours as in the paper's figures (mpce_results.COL): PSO-VNS blue, SSA-VNS sky blue, PSO green
OI = dict(blue="#0072B2", orange="#E69F00", green="#009E73", verm="#D55E00", sky="#56B4E9",
          purple="#CC79A7", yellow="#F0E442", black="#000000")
GREY, LIGHT, INK = "#6E6E6E", "#EFEFEF", "#1A1A1A"
plt.rcParams.update({
    "font.family": "Liberation Sans", "font.size": 9, "mathtext.fontset": "custom", "mathtext.rm": "Liberation Sans",
    "mathtext.bf": "Liberation Sans:bold", "mathtext.it": "Liberation Sans:italic",
    "mathtext.sf": "Liberation Sans", "mathtext.cal": "Liberation Sans", "mathtext.tt": "Liberation Mono", "pdf.fonttype": 42, "ps.fonttype": 42,
    "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "xtick.major.size": 2, "ytick.major.size": 2, "xtick.major.pad": 1.5, "ytick.major.pad": 1.5,
    "text.color": INK, "axes.labelcolor": INK, "xtick.color": INK, "ytick.color": INK,
})
W, H = 130.0, 52.0                     # mm
MM = 1 / 25.4
fig = plt.figure(figsize=(W * MM, H * MM))
fig.patch.set_facecolor("white")


def ax_mm(x, y, w, h, **kw):
    return fig.add_axes([x / W, y / H, w / W, h / H], **kw)


def text_mm(x, y, s, **kw):
    return fig.text(x / W, y / H, s, **kw)


def box_mm(x, y, w, h, fc, ec="none", lw=0.0, r=1.5):
    fig.patches.append(FancyBboxPatch((x / W, y / H), w / W, h / H, transform=fig.transFigure,
                                      boxstyle=f"round,pad=0,rounding_size={r / H}",
                                      fc=fc, ec=ec, lw=lw, zorder=-10, mutation_aspect=W / H))


def arrow_mm(x0, x1, y):
    fig.patches.append(FancyArrowPatch((x0 / W, y / H), (x1 / W, y / H), transform=fig.transFigure,
                                       arrowstyle="-|>,head_length=3,head_width=1.8", color=GREY,
                                       lw=1.4, zorder=3))


# ---- layout (mm) ----
LX, LW = 0.8, 24.2                # test bed
MX, MW = 27.8, 73.4               # four protocol choices (2 x 2)
RX, RW = 104.0, 25.2              # take-away
TOP = 51.0

# ---- left: test bed ----
box_mm(LX, 1.0, LW, 50.0, LIGHT)
text_mm(LX + LW / 2, TOP - 1.0, "Test bed", ha="center", va="top", fontsize=10, fontweight="bold")
ia = ax_mm(LX + 1.0, 21.0, LW - 2.0, 24.5)
ia.set_xlim(-1.35, 1.15)
ia.set_ylim(-1.15, 1.15)
ia.set_aspect("equal")
ia.axis("off")
ia.add_patch(Circle((0, 0), 1.0, fc="white", ec=GREY, lw=0.9, ls=(0, (3, 1.5))))
turbines = [(-0.62, 0.48), (0.05, 0.72), (0.62, 0.42), (-0.70, -0.18), (-0.12, 0.10), (0.55, -0.25),
            (-0.35, -0.70), (0.25, -0.72), (0.80, 0.05)]
for (x, y) in turbines:          # top view: rotor bar perpendicular to the wind (from the west)
    ia.plot([x, x], [y - 0.15, y + 0.15], color=OI["blue"], lw=1.6, solid_capstyle="round")
    ia.plot([x, x + 0.13], [y, y], color=OI["blue"], lw=1.0)
for yy in (-0.45, 0.0, 0.45):
    ia.annotate("", xy=(-1.04, yy), xytext=(-1.38, yy),
                arrowprops=dict(arrowstyle="-|>,head_length=0.35,head_width=0.2", color=GREY, lw=0.9))
lines = [("68", "layout cases"), ("30", "paired seeds"), ("equal", "budgets")]
for i, (big, small) in enumerate(lines):
    y = 17.0 - i * 5.6
    text_mm(LX + LW / 2, y, f"$\\bf{{{big}}}$ {small}" if big != "equal" else "equal budgets",
            ha="center", va="center", fontsize=9)

# ---- middle: four choices ----
text_mm(MX + MW / 2, TOP - 0.6, "Protocol choices change the conclusions",
        ha="center", va="top", fontsize=10, fontweight="bold")
PW, PH = 36.3, 21.8
pos = [(MX, 23.6), (MX + MW - PW, 23.6), (MX, 1.0), (MX + MW - PW, 1.0)]


def panel(k, title, verdict, vcolor):
    x, y = pos[k]
    box_mm(x, y, PW, PH, "white", ec="#BDBDBD", lw=0.6)
    text_mm(x + 1.4, y + PH - 1.0, title, ha="left", va="top", fontsize=9, fontweight="bold")
    text_mm(x + PW / 2, y + 1.0, verdict, ha="center", va="bottom", fontsize=9, fontweight="bold",
            color=vcolor)
    return x, y


def slope(ax, a0, a1, b0, b1, ca, cb, la, lb, xt, ylim, yticks, right_axis=True):
    """Two-point slope chart; method labels left of the first point, rank ticks on the right."""
    for v0, v1, c in ((a0, a1, ca), (b0, b1, cb)):
        ax.plot([0, 1], [v0, v1], color=c, lw=2.0, marker="o", ms=3.6, solid_capstyle="round", zorder=3)
    ax.text(-0.13, a0, la, color=ca, va="center", ha="right", fontsize=9, fontweight="bold")
    ax.text(-0.13, b0, lb, color=cb, va="center", ha="right", fontsize=9, fontweight="bold")
    ax.set_xlim(-0.08, 1.12)
    ax.set_ylim(*ylim)
    ax.set_xticks([0, 1], xt)
    ax.yaxis.tick_right()
    ax.set_yticks(yticks if right_axis else [])
    for sp in ("top", "left") if right_axis else ("top", "left", "right"):
        ax.spines[sp].set_visible(False)
    ax.tick_params(axis="x", length=0)


# A: PSO baseline setting
x, y = panel(0, "PSO baseline setting", f"PSO rank {pso_old:.1f} → {pso_con:.1f}", OI["green"])
a = ax_mm(x + 15.0, y + 8.6, 16.0, 7.6)
slope(a, pso_old, pso_con, ssa_old, ssa_con, OI["green"], OI["sky"], "PSO", "SSA-VNS",
      ["old", "constr."], (5.6, 0.8), [1, 3, 5])

# B: random-sampling control
x, y = panel(1, "Sampling control (4D)", "no gain vs disc", OI["orange"])
b = ax_mm(x + 11.5, y + 8.0, 21.5, 9.6)
vals = [gain_disc, gain_sq]
b.barh([0, 1], vals, height=0.62, color=[GREY, OI["orange"]], zorder=2)
b.axvline(0, color=INK, lw=0.6, zorder=3)
b.set_yticks([0, 1], ["disc", "square"])
b.set_ylim(-0.55, 1.55)
b.set_xlim(-0.06, 0.13)
b.set_xticks([])
for v, yy in zip(vals, (0, 1)):
    b.text(max(v, 0) + 0.006, yy, f"{v:+.2f}".replace("-", "−"), va="center", ha="left",
           fontsize=9)
for s in ("top", "right", "bottom", "left"):
    b.spines[s].set_visible(False)
b.tick_params(axis="y", length=0)
text_mm(x + 11.5 + 21.5 * 0.06 / 0.19, y + 6.4, "SSA-phase gain (pp)", ha="center", va="center",
        fontsize=9, color=GREY)

# C: constraint handling
x, y = panel(2, "Constraint handling", "leader changes", OI["verm"])
c = ax_mm(x + 16.3, y + 8.6, 13.4, 7.6)
slope(c, pv_pen, pv_deb, ga_pen, ga_deb, OI["blue"], OI["verm"], "PSO-VNS", "GA",
      ["penalty", "Deb+clip"], (2.75, 1.45), [1.5, 2.5], right_axis=False)
c.tick_params(axis="x", labelsize=8.5)

# D: wind-direction bins
x, y = panel(3, "Wind bins (eval. 1°)", "leader stays, gap ≈×2", OI["blue"])
d = ax_mm(x + 12.5, y + 8.0, 20.5, 9.6)
vals = [adv1, adv15]
d.barh([0, 1], vals, height=0.62, color=[OI["blue"], "#7FB3D9"], zorder=2)
d.axvline(0, color=INK, lw=0.6, zorder=3)
d.set_yticks([0, 1], ["opt. 1°", "opt. 15°"])
d.set_ylim(-0.55, 1.55)
d.set_xlim(0, 0.25)
d.set_xticks([])
for v, yy in zip(vals, (0, 1)):
    d.text(v + 0.008, yy, f"{v:.2f}", va="center", ha="left", fontsize=9)
for s in ("top", "right", "bottom", "left"):
    d.spines[s].set_visible(False)
d.tick_params(axis="y", length=0)
text_mm(x + PW / 2 + 1.5, y + 6.4, "gain over PSO (pp)", ha="center", va="center", fontsize=9,
        color=GREY)

# ---- right: take-away ----
box_mm(RX, 1.0, RW, 50.0, "#E3F0F8")
text_mm(RX + RW / 2, TOP - 1.0, "Report", ha="center", va="top", fontsize=10, fontweight="bold")
items = [("baseline", "settings"), ("in-farm", "controls"), ("all-run", "reliability"),
         ("constraint &", "bin checks")]
for i, (l1, l2) in enumerate(items):
    yy = 41.5 - i * 10.4
    text_mm(RX + 1.6, yy, "✓", ha="left", va="center", fontsize=10, fontweight="bold", fontfamily="DejaVu Sans",
            color=OI["blue"])
    text_mm(RX + 5.0, yy + 1.9, l1, ha="left", va="center", fontsize=9)
    text_mm(RX + 5.0, yy - 1.9, l2, ha="left", va="center", fontsize=9)

arrow_mm(LX + LW + 0.4, MX - 0.6, 26.0)
arrow_mm(MX + MW + 0.4, RX - 0.6, 26.0)

os.makedirs(OUT, exist_ok=True)
pdf = os.path.join(OUT, "graphical_abstract.pdf")
png = os.path.join(OUT, "graphical_abstract.png")
tif = os.path.join(OUT, "graphical_abstract.tif")
meta = {"Title": "Graphical abstract", "Creator": "analysis/make_graphical_abstract.py", "CreationDate": None}
fig.savefig(pdf, metadata=meta)
fig.savefig(png, dpi=600)          # 3071 x 1228 px
try:
    from PIL import Image
    fig.savefig(tif + ".tmp.png", dpi=300)
    im = Image.open(tif + ".tmp.png").convert("RGB")
    im.save(tif, compression="tiff_lzw", dpi=(300, 300))
    os.remove(tif + ".tmp.png")
except Exception as e:  # PIL without TIFF support
    print("TIFF not written:", e)
    tif = None
outs = [p for p in (pdf, png, tif) if p]
if os.path.isdir(COPY_TO):
    for p in outs:
        shutil.copy2(p, COPY_TO)
for p in outs:
    print(p)
