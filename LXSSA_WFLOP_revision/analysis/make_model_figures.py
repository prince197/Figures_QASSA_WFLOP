"""Model figures of Section III (vector PDF + PNG in ../figures_final):
fig_wind_farm   circular farm, direction convention and an example layout with 4D exclusion zones
fig_wake_model  Jensen top-hat wake behind a rotor (side view)
fig_half_cone   wake half-cone geometry: virtual vertex A, angle alpha, test angle beta_ij
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Polygon, FancyArrowPatch, Arc, Wedge, Rectangle, Ellipse

OUT = "../figures_final"
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e4e3df"
BLUE, BLUE_L, BLUE_XL = "#2a78d6", "#86b6ef", "#e6f0fc"
ORANGE, GREEN = "#eb6834", "#008300"
plt.rcParams.update({"font.family": "STIXGeneral", "mathtext.fontset": "stix", "font.size": 9,
                     "text.color": INK, "pdf.fonttype": 42, "axes.linewidth": 0})


def save(fig, name):
    os.makedirs(OUT, exist_ok=True)
    fig.savefig(f"{OUT}/{name}.pdf", bbox_inches="tight", pad_inches=0.02)
    fig.savefig(f"{OUT}/{name}.png", dpi=220, bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)


def arrow(ax, p, q, color=INK, lw=0.9, style="-|>", ms=7, both=False, z=4):
    ax.add_patch(FancyArrowPatch(p, q, arrowstyle=("<|-|>" if both else style), mutation_scale=ms,
                                 color=color, lw=lw, shrinkA=0, shrinkB=0, zorder=z))


def dim(ax, p, q, text, off, color=MUTED, rot=0, tpos=0.5, ha="center", va="center", fs=9):
    arrow(ax, p, q, color=color, lw=0.7, both=True, ms=6)
    mx, my = p[0] + tpos * (q[0] - p[0]), p[1] + tpos * (q[1] - p[1])
    ax.text(mx + off[0], my + off[1], text, color=INK, rotation=rot, ha=ha, va=va, fontsize=fs)


# ------------------------------------------------------------------ Fig. 1: circular wind farm
def wind_farm():
    r = 750.0
    B = pd.read_csv("final_best_layouts_maxN.csv")
    row = B[(B.Dataset == 1) & (B.Radius == 750) & (B.Algorithm == "LXBV")].sort_values("Objective").iloc[-1]
    xy = np.array([[float(v) for v in p.split()] for p in row.Coordinates.split(";")])
    fig, ax = plt.subplots(figsize=(3.45, 3.55))
    ax.set_aspect("equal"); ax.axis("off")
    ax.add_patch(Circle((0, 0), r, facecolor="#f7f9fc", edgecolor=INK, lw=1.1, zorder=1))
    for rr in (0.25 * r, 0.5 * r, 0.75 * r):
        ax.add_patch(Circle((0, 0), rr, fill=False, edgecolor=GRID, lw=0.6, zorder=1))
    for a in range(0, 360, 30):
        t = np.radians(a)
        ax.plot([0, r * np.cos(t)], [0, r * np.sin(t)], color=GRID, lw=0.6, zorder=1)
    for a, lab in ((0, "0°\nE"), (90, "90° N"), (180, "180°\nW"), (270, "270° S")):
        t = np.radians(a)
        ha = {0: "left", 180: "right"}.get(a, "center"); va = {90: "bottom", 270: "top"}.get(a, "center")
        ax.text(1.05 * r * np.cos(t), 1.05 * r * np.sin(t), lab, ha=ha, va=va, fontsize=8.5, color=INK,
                linespacing=1.0)
    # layout with 4D exclusion zones (radius 2D = half the minimum spacing)
    for x, y in xy:
        ax.add_patch(Circle((x, y), 154, facecolor=BLUE, alpha=0.07, edgecolor="none", zorder=2))
        ax.add_patch(Circle((x, y), 154, fill=False, edgecolor=BLUE_L, lw=0.6, ls=(0, (2, 1.5)), zorder=2))
    ax.scatter(xy[:, 0], xy[:, 1], s=18, color=BLUE, edgecolor="white", lw=0.6, zorder=4)
    # direction angle theta and wind direction
    th = np.radians(135)
    ax.add_patch(Arc((0, 0), 0.36 * r, 0.36 * r, theta1=0, theta2=135, color=ORANGE, lw=1.1, zorder=5))
    ax.text(0.25 * r * np.cos(np.radians(62)), 0.25 * r * np.sin(np.radians(62)), r"$\theta$", color=ORANGE,
            fontsize=11, ha="center", va="center", zorder=6)
    arrow(ax, (0, 0), (0.62 * r * np.cos(th), 0.62 * r * np.sin(th)), color=ORANGE, lw=1.5, ms=10, z=6)
    ax.plot([0, r], [0, 0], color=ORANGE, lw=0.7, ls=(0, (2, 1.5)), zorder=5)
    # radius
    ang = np.radians(-62)
    arrow(ax, (0, 0), (r * np.cos(ang), r * np.sin(ang)), color=INK, lw=0.8, ms=6, z=5)
    ax.text(0.5 * r * np.cos(ang) + 30, 0.5 * r * np.sin(ang), r"$r$", fontsize=11, ha="left", va="center")
    ax.scatter([0], [0], s=8, color=INK, zorder=7)
    # key
    ky = -1.3 * r
    ax.scatter([-0.95 * r], [ky], s=18, color=BLUE, edgecolor="white", lw=0.6, clip_on=False)
    ax.text(-0.87 * r, ky, "turbine", va="center", fontsize=8)
    ax.add_patch(Circle((-0.22 * r, ky), 55, facecolor=BLUE, alpha=0.07, edgecolor="none", clip_on=False))
    ax.add_patch(Circle((-0.22 * r, ky), 55, fill=False, edgecolor=BLUE_L, lw=0.6, ls=(0, (2, 1.5)), clip_on=False))
    ax.text(-0.12 * r, ky, r"radius $2D$ (half the $4D$ spacing)", va="center", fontsize=8)
    arrow(ax, (-0.95 * r - 30, ky - 0.13 * r), (-0.72 * r, ky - 0.13 * r), color=ORANGE, lw=1.2, ms=8)
    ax.text(-0.66 * r, ky - 0.13 * r, r"wind direction $\theta$ (from east, counter-clockwise)", va="center",
            fontsize=8)
    ax.set_xlim(-1.3 * r, 1.3 * r); ax.set_ylim(-1.5 * r, 1.17 * r)
    save(fig, "fig_wind_farm")


# ------------------------------------------------------------------ Fig. 2: Jensen wake (side view)
def turbine(ax, x, R, hub=0.0, width=0.22):
    ax.add_patch(Ellipse((x, hub), width * R, 2 * R, facecolor="#c9d6e6", edgecolor=INK, lw=0.9, zorder=5))
    ax.add_patch(Rectangle((x + 0.06 * R, hub - 0.1 * R), 0.45 * R, 0.2 * R, facecolor="#8a8f96",
                           edgecolor=INK, lw=0.6, zorder=6))
    ax.add_patch(Circle((x, hub), 0.07 * R, facecolor=INK, zorder=7))


def wake_model():
    R, K, d = 1.0, 0.2, 8.0             # K exaggerated for readability (benchmark K = 0.075)
    fig, ax = plt.subplots(figsize=(3.45, 2.05))
    ax.set_aspect("equal"); ax.axis("off")
    x0, xe = 0.0, d
    Rw = R + K * d
    ax.add_patch(Polygon([(x0, R), (xe, Rw), (xe, -Rw), (x0, -R)], closed=True, facecolor=BLUE_XL,
                         edgecolor="none", zorder=1))
    ax.plot([x0, xe], [R, Rw], color=BLUE, lw=0.9, ls=(0, (4, 2)), zorder=2)
    ax.plot([x0, xe], [-R, -Rw], color=BLUE, lw=0.9, ls=(0, (4, 2)), zorder=2)
    ax.plot([xe, xe], [-Rw, Rw], color=BLUE, lw=0.8, zorder=2)
    turbine(ax, x0, R, width=0.3)
    for y in np.linspace(-1.6, 1.6, 7):                      # free stream
        arrow(ax, (-4.6, y), (-2.4, y), color=MUTED, lw=0.8, ms=6)
    ax.text(-3.5, 1.85, r"$s_{\mathrm{up}}$", ha="center", va="bottom", fontsize=10)
    for y in np.linspace(-1.2, 1.2, 5):                      # reduced speed in the wake
        arrow(ax, (xe - 3.0, y), (xe - 1.6, y), color=BLUE, lw=0.9, ms=6)
    ax.text(xe - 2.3, 1.45, r"$s_{\mathrm{down}}$", ha="center", va="bottom", fontsize=10, color=BLUE)
    for y in (Rw + 0.45, -Rw - 0.45):                        # bypass flow
        arrow(ax, (xe - 3.0, y), (xe - 0.2, y), color=MUTED, lw=0.8, ms=6)
    ax.text(xe - 3.2, Rw + 0.45, r"$s_{\mathrm{up}}$", ha="right", va="center", fontsize=9)
    dim(ax, (x0 - 0.75, -R), (x0 - 0.75, R), r"$2R$", (-0.2, -0.45), ha="right", fs=9)
    dim(ax, (xe + 0.5, -Rw), (xe + 0.5, Rw), r"$2(R+Kd)$", (0.2, 0), ha="left")
    dim(ax, (x0, -Rw - 1.0), (xe, -Rw - 1.0), r"$d$", (0, -0.25), va="top")
    ax.plot([x0, x0], [-R - 0.15, -Rw - 1.15], color=GRID, lw=0.6, zorder=0)
    ax.plot([xe, xe], [-Rw - 0.1, -Rw - 1.15], color=GRID, lw=0.6, zorder=0)
    ax.text(x0 + 2.4, 0.0, "wake", color=BLUE, fontsize=9, ha="center", va="center", style="italic")
    ax.set_xlim(-5.0, xe + 3.3); ax.set_ylim(-Rw - 1.65, Rw + 0.85)
    save(fig, "fig_wake_model")


# ------------------------------------------------------------------ Fig. 3: half-cone geometry
def half_cone():
    R, K = 1.0, 0.25                    # K exaggerated for readability (benchmark K = 0.075)
    xA, d = -R / K, 6.0
    alpha = np.degrees(np.arctan(K))
    fig, ax = plt.subplots(figsize=(3.45, 2.15))
    ax.set_aspect("equal"); ax.axis("off")
    xe, Rw = d, R + K * d
    ax.add_patch(Polygon([(xA, 0), (xe, Rw), (xe, -Rw)], closed=True, facecolor=BLUE_XL, edgecolor="none", zorder=1))
    ax.plot([xA, xe], [0, Rw], color=BLUE, lw=0.9, ls=(0, (4, 2)), zorder=2)
    ax.plot([xA, xe], [0, -Rw], color=BLUE, lw=0.9, ls=(0, (4, 2)), zorder=2)
    ax.plot([xA, xe + 0.4], [0, 0], color=MUTED, lw=0.7, ls=(0, (1, 1.5)), zorder=2)
    turbine(ax, 0.0, R, width=0.3)
    ax.text(0.0, R + 0.2, r"turbine $j$", ha="center", va="bottom", fontsize=8.5)
    ax.scatter([xA], [0], s=16, color=INK, zorder=6)
    ax.text(xA - 0.2, 0, r"$A$", ha="right", va="center", fontsize=10)
    ax.add_patch(Arc((xA, 0), 3.0, 3.0, theta1=0, theta2=alpha, color=BLUE, lw=1.0, zorder=3))
    ax.text(xA + 1.65, 0.28, r"$\alpha$", color=BLUE, fontsize=10, va="center")
    ti = (4.2, -1.05)                                        # downstream turbine i (inside the wake)
    ax.plot([xA, ti[0]], [0, ti[1]], color=ORANGE, lw=0.9, zorder=3)
    ax.scatter([ti[0]], [ti[1]], s=24, color=ORANGE, edgecolor="white", lw=0.6, zorder=6)
    ax.text(ti[0] + 0.25, ti[1] - 0.1, r"turbine $i$", fontsize=8.5, va="top")
    beta = np.degrees(np.arctan2(ti[1], ti[0] - xA))
    ax.add_patch(Arc((xA, 0), 5.6, 5.6, theta1=beta, theta2=0, color=ORANGE, lw=1.0, zorder=3))
    ax.text(xA + 2.95, -0.36, r"$\beta_{ij}$", color=ORANGE, fontsize=10, va="center")
    dim(ax, (0.55, 0.14), (0.55, R), r"$R$", (0.2, 0.05), ha="left", fs=9)
    dim(ax, (xA, -Rw - 0.45), (0, -Rw - 0.45), r"$R/K$", (0, -0.25), va="top")
    dim(ax, (0, -Rw - 0.45), (xe, -Rw - 0.45), r"$d$", (0, -0.25), va="top")
    dim(ax, (xe + 0.55, -Rw), (xe + 0.55, Rw), r"$2(R+Kd)$", (0.25, 0), ha="left")
    for x, y0 in ((xA, -0.1), (0, -R - 0.1), (xe, -Rw)):
        ax.plot([x, x], [y0, -Rw - 0.6], color=GRID, lw=0.6, zorder=0)
    ax.text(xA + 0.2, Rw + 0.1, r"turbine $i$ is waked by $j$ if $\beta_{ij}<\alpha=\arctan K$", fontsize=8,
            color=MUTED, ha="left", va="bottom")
    ax.set_xlim(xA - 0.8, xe + 3.1); ax.set_ylim(-Rw - 1.1, Rw + 0.55)
    save(fig, "fig_half_cone")


if __name__ == "__main__":
    wind_farm(); wake_model(); half_cone()
    print("model figures written")
