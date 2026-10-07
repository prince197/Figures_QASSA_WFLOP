"""Graphical abstract for the SWEVO submission (Elsevier: single panel, >= 531 x 1328 px, 5 x 13 cm).

Usage (from anywhere):  python3 analysis/make_graphical_abstract.py
Writes figures_final/graphical_abstract.{pdf,png,tif} and copies them to
SWEVO_rev2/latex_source/figures_final/.

Left: the two phases of PSO-VNS; middle: the evaluation; right: the findings, including where PSO-VNS does not lead.
Every number shown is read from the analysis JSON files and checked against the value printed in the paper:
  * Average ranks of the main comparison: mpce_summary.json /main/friedman/avg_rank
  * Direction bins: rev3_fine.json /arms/{ctrl1,direct}/contrast_psovns_pso/mean_dloss_pp and /arms/{ctrl1,direct}/leader
  * Horns Rev 1 feasibility at 6,030 evaluations: mpce_summary.json /hr16/methods/{PSOBV,PSOC}/feasible
The canvas is 130 x 52 mm (aspect 2.5 = 1328/531); all text is >= 9 pt at that size.
"""
import json
import os
import shutil

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(ROOT, "figures_final")
COPY_TO = os.path.join(ROOT, "SWEVO_rev2", "latex_source", "figures_final")


def load(name):
    with open(os.path.join(HERE, name)) as f:
        return json.load(f)


S = load("mpce_summary.json")
F = load("rev3_fine.json")

# ---- numbers (with the value printed in the paper as a check) ----
rank = S["main"]["friedman"]["avg_rank"]
top = sorted(rank, key=rank.get)[:4]                                   # PSO-VNS, PSO, SSA-VNS, VNS
adv15 = -F["arms"]["ctrl1"]["contrast_psovns_pso"]["mean_dloss_pp"]   # 0.086 (opt. 15 deg, eval. 1 deg)
adv1 = -F["arms"]["direct"]["contrast_psovns_pso"]["mean_dloss_pp"]   # 0.194 (opt. 1 deg)
lead15, lead1 = F["arms"]["ctrl1"]["leader"], F["arms"]["direct"]["leader"]
hr_pv, hr_pso = S["hr16"]["methods"]["PSOBV"]["feasible"], S["hr16"]["methods"]["PSOC"]["feasible"]   # 30, 19
checks = [(rank["PSOBV"], 1.82, 2), (rank["PSOC"], 1.97, 2), (rank["SSABV"], 3.45, 2), (rank["BVNS"], 4.40, 2),
          (adv15, 0.086, 3), (adv1, 0.194, 3)]
for val, paper, nd in checks:
    assert round(val, nd) == paper, (val, paper)
assert top == ["PSOBV", "PSOC", "SSABV", "BVNS"] and lead15 == lead1 == ["PSOBV"] and (hr_pv, hr_pso) == (30, 19)
NAME = {"PSOBV": "PSO-VNS", "PSOC": "PSO", "SSABV": "SSA-VNS", "BVNS": "VNS"}

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
LX, LW = 0.8, 34.0                # proposed method
MX, MW = 38.0, 59.6               # evaluation
RX, RW = 100.4, 28.8              # findings
TOP = 51.0

# ---- left: PSO-VNS ----
box_mm(LX, 1.0, LW, 50.0, LIGHT)
text_mm(LX + LW / 2, TOP - 1.0, "PSO-VNS", ha="center", va="top", fontsize=10, fontweight="bold",
        color=OI["blue"])


def phase(y, h, title, l1, l2, l3, fc):
    box_mm(LX + 1.6, y, LW - 3.2, h, fc, ec="#BDBDBD", lw=0.6)
    cx = LX + LW / 2
    text_mm(cx, y + h - 1.4, title, ha="center", va="top", fontsize=9, fontweight="bold")
    for i, t in enumerate((l1, l2, l3)):
        text_mm(cx, y + h - 6.2 - i * 3.9, t, ha="center", va="top", fontsize=9,
                color=GREY if i == 2 else INK)


phase(25.4, 18.6, "Phase 1: PSO", "swarm reaches the", "feasible region", "ωB evaluations", "white")
phase(2.6, 18.6, "Phase 2: VNS", "compass search and", "shaking refine it", "(1−ω)B evaluations", "white")
fig.patches.append(FancyArrowPatch(((LX + LW / 2) / W, 25.2 / H), ((LX + LW / 2) / W, 21.4 / H),
                                   transform=fig.transFigure, arrowstyle="-|>,head_length=2.5,head_width=1.5",
                                   color=OI["blue"], lw=1.4, zorder=3))
text_mm(LX + LW / 2 + 2.0, 23.3, "best layout", ha="left", va="center", fontsize=8, color=GREY)

# ---- middle: evaluation ----
text_mm(MX + MW / 2, TOP - 0.6, "Controlled evaluation", ha="center", va="top", fontsize=10, fontweight="bold")
text_mm(MX + MW / 2, 44.6, "2 wind data sets × 3 farm sizes", ha="center", va="top", fontsize=9)
text_mm(MX + MW / 2, 40.6, "Horns Rev 1 · Lillgrund · IEA37", ha="center", va="top", fontsize=9)
PW, PH = 29.0, 34.0
for k, (x, title, verdict) in enumerate(((MX, "Average rank", "1st of 8"),
                                         (MX + MW - PW, "Gain vs PSO (pp)", "≈×2 with 1° bins"))):
    box_mm(x, 1.0, PW, PH, "white", ec="#BDBDBD", lw=0.6)
    text_mm(x + PW / 2, 1.0 + PH - 1.2, title, ha="center", va="top", fontsize=9, fontweight="bold")
    text_mm(x + PW / 2, 2.2, verdict, ha="center", va="bottom", fontsize=9, fontweight="bold", color=OI["blue"])

a = ax_mm(MX + 12.6, 8.0, 12.0, 20.0)
vals = [rank[m] for m in top]
cols = [OI["blue"], OI["green"], OI["sky"], GREY]
a.barh(range(4), vals, height=0.62, color=cols, zorder=2)
a.set_yticks(range(4), [NAME[m] for m in top])
for t, c in zip(a.get_yticklabels(), cols):
    t.set_color(c)
    t.set_fontweight("bold")
a.set_ylim(3.6, -0.6)
a.set_xlim(0, 6.6)
a.set_xticks([])
for v, yy in zip(vals, range(4)):
    a.text(v + 0.25, yy, f"{v:.2f}", va="center", ha="left", fontsize=8.5)
for sp in ("top", "right", "bottom"):
    a.spines[sp].set_visible(False)
a.tick_params(axis="y", length=0, labelsize=8.5)

d = ax_mm(MX + MW - PW + 11.0, 11.0, 12.0, 14.0)
vals = [adv15, adv1]
d.barh([0, 1], vals, height=0.62, color=["#7FB3D9", OI["blue"]], zorder=2)
d.axvline(0, color=INK, lw=0.6, zorder=3)
d.set_yticks([0, 1], ["opt. 15°", "opt. 1°"])
d.set_ylim(1.6, -0.6)
d.set_xlim(0, 0.27)
d.set_xticks([])
for v, yy in zip(vals, (0, 1)):
    d.text(v + 0.01, yy, f"{v:.2f}", va="center", ha="left", fontsize=8.5)
for sp in ("top", "right", "bottom", "left"):
    d.spines[sp].set_visible(False)
d.tick_params(axis="y", length=0, labelsize=8.5)
text_mm(MX + MW - PW / 2, 8.6, "evaluated at 1°", ha="center", va="center", fontsize=8, color=GREY)

# ---- right: findings ----
box_mm(RX, 1.0, RW, 50.0, "#E3F0F8")
text_mm(RX + RW / 2, TOP - 1.0, "Findings", ha="center", va="top", fontsize=10, fontweight="bold")
items = [("✓", "ranks first", "of 8 methods"), ("✓", "+0.19 pp over", "PSO with 1° bins"),
         ("✓", f"feasible {hr_pv}/30", f"vs PSO {hr_pso}/30*"), ("–", "GA, SLSQP lead", "in 2 settings")]
for i, (mark, l1, l2) in enumerate(items):
    yy = 41.5 - i * 10.4
    text_mm(RX + 1.4, yy, mark, ha="left", va="center", fontsize=10, fontweight="bold", fontfamily="DejaVu Sans",
            color=OI["blue"] if mark == "✓" else GREY)
    text_mm(RX + 4.6, yy + 1.9, l1, ha="left", va="center", fontsize=9)
    text_mm(RX + 4.6, yy - 1.9, l2, ha="left", va="center", fontsize=9)
text_mm(RX + RW / 2, 2.0, "*Horns Rev 1", ha="center", va="bottom", fontsize=8, color=GREY)

arrow_mm(LX + LW + 0.3, MX - 0.5, 26.0)
arrow_mm(MX + MW + 0.3, RX - 0.5, 26.0)

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
