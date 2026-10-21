"""Figures for Cone Matching and Integral Constraints in Relative Langlands Duality."""
from __future__ import annotations

from fractions import Fraction as F
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

import cones
import examples as ex

INK, GRID, BLUE, RED, GRAY, AMBER = "#1a1a1a", "#d9d9d9", "#2563eb", "#b3202c", "#78716c", "#b45309"


def _style(ax):
    ax.grid(True, color=GRID, lw=0.5)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.tick_params(labelsize=10)


def plot_cones(path, z=6):
    """A2 slice at fixed z: the automorphic cone (A), equal to (S_1), and the spectral cone (S_3) of GL3."""
    fig, ax = plt.subplots(figsize=(5.6, 6.0))
    # (A) slice: max(x1, x2) <= z <= x1 + x2 ; (S_d) slice: max(x1, x2) <= z/d <= x1 + x2
    ax.add_patch(Polygon([(z, 0), (z, z), (0, z)], closed=True, facecolor=BLUE, alpha=0.18, edgecolor=BLUE, lw=1.2,
                         label="cone (A) = cone (S) with $d=1$"))
    w = z / 3
    ax.add_patch(Polygon([(w, 0), (w, w), (0, w)], closed=True, facecolor=RED, alpha=0.18, edgecolor=RED, lw=1.2,
                         label="cone (S) with $d=3$"))
    pts = {"A": [], "S": [], "both": []}
    for x1 in range(z + 1):
        for x2 in range(z + 1):
            if (z - x1 - 2 * x2) % 3:  # GL3: z = a + b + c = 3c + x1 + 2 x2
                continue
            x = [F(x1), F(x2)]
            a, s = cones.in_A(ex.A2, x, z), cones.in_S(ex.A2, x, z, 3)
            if a and s:
                pts["both"].append((x1, x2))
            elif a:
                pts["A"].append((x1, x2))
            elif s:
                pts["S"].append((x1, x2))
    for key, col, mk, lab in (("A", BLUE, "o", "GL$_3$ label in (A) only"), ("S", RED, "s", "GL$_3$ label in (S), $d=3$, only")):
        if pts[key]:
            ax.plot(*zip(*pts[key]), mk, color=col, ms=6, label=lab)
    ax.set_xlim(-0.4, z + 0.6)
    ax.set_ylim(-0.4, z + 0.6)
    ax.set_aspect("equal")
    ax.set_xlabel(r"$x_1=\langle\chi,\alpha_1\rangle$", fontsize=11)
    ax.set_ylabel(r"$x_2=\langle\chi,\alpha_2\rangle$", fontsize=11)
    ax.set_title(rf"type $A_2$, slice $z=\langle\chi,\theta_G\rangle={z}$", fontsize=11, loc="left")
    ax.legend(fontsize=9, frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.13), ncol=2)
    _style(ax)
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def plot_integrality(r, path, z=6, p1=2):
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 5.2))
    ax = axes[0]
    good, bad = [], []
    for p2 in range(z + 1):
        for p3 in range(z + 1):
            p = [F(p1), F(p2), F(p3)]
            if cones.in_A(ex.A1_3, p, z):
                ((good if (p1 + p2 + p3 - z) % 2 == 0 else bad)).append((p2, p3))
    ax.plot(*zip(*good), "o", color=BLUE, ms=7, label="in both supports")
    ax.plot(*zip(*bad), "x", color=RED, ms=8, mew=2, label=r"in the cone, $(p_1+p_2+p_3-z)/2\notin\mathbb{Z}$")
    ax.set_xlim(-0.5, z + 0.5)
    ax.set_ylim(-0.5, z + 0.5)
    ax.set_aspect("equal")
    ax.set_xlabel(r"$p_2$", fontsize=11)
    ax.set_ylabel(r"$p_3$", fontsize=11)
    ax.set_title(rf"a. labels with $p_1={p1}$, $z={z}$", fontsize=11, loc="left")
    ax.legend(fontsize=9, frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.14), ncol=1)
    _style(ax)
    ax = axes[1]
    rows = r["examples"]["gl2_pgl2_pgl2"]
    zs = [row["z"] for row in rows]
    ax.plot(zs, [row["labels_A"] for row in rows], "o-", color=BLUE, ms=4, lw=1.4, label="automorphic support")
    ax.plot(zs, [row["labels_S"] for row in rows], "s-", color=RED, ms=4, lw=1.4, label="spectral support with the integrality\ncondition of Proposition 5")
    ax.set_xlabel(r"central degree $z$", fontsize=11)
    ax.set_ylabel("number of labels", fontsize=11)
    ax.set_title("b. labels by degree", fontsize=11, loc="left")
    ax.legend(fontsize=9, frameon=False, loc="upper left")
    _style(ax)
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def make_all(r, out):
    out = Path(out)
    plot_cones(out / "cones.png")
    plot_integrality(r, out / "integrality.png")
