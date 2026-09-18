"""CONTINUATION -- one figure: the canonical magnitude/comparator split of g(zeta).

Palette: same validated reference instance as the parent audit
(categorical slots 1-3, light surface, recessive grid; light mode only).

Run:  ../.venv/bin/python src/comparator_figure.py
"""

from __future__ import annotations

import pathlib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

_ROOT = next(p for p in pathlib.Path(__file__).resolve().parents
             if (p / "src").is_dir())
HERE = _ROOT
FIG = HERE / "figures"
FIG.mkdir(exist_ok=True)

SURFACE, INK, INK2, MUTED, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#8a8985", "#e6e5e1"
C1, C2, C3 = "#2a78d6", "#eb6834", "#1baf7a"
plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "font.size": 9.5, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "text.color": INK,
    "xtick.color": INK2, "ytick.color": INK2, "axes.linewidth": 0.8,
    "grid.color": GRID, "grid.linewidth": 0.7, "lines.linewidth": 2.0,
    "legend.frameon": False, "figure.dpi": 150, "axes.titlesize": 10,
    "axes.titleweight": "bold", "axes.spines.top": False, "axes.spines.right": False,
})

z = np.linspace(-0.97, 0.97, 800)
h = 0.5 * np.sqrt(1 - z ** 2)
xi = np.arctanh(z)
g = (1 + z) * np.sqrt((1 - z) / 2)
g_sym = 2 * h ** 1.5
g_cmp = np.exp(xi / 2)

fig, ax = plt.subplots(1, 2, figsize=(10.0, 3.9))

ax[0].plot(z, g, color=C1, label=r"$g(\zeta)$  (allocation factor)")
ax[0].plot(z, g_sym, color=C2, ls="--", label=r"$g_{\rm sym}=2h^{3/2}$  (even)")
ax[0].plot(z, g_cmp, color=C3, ls="-.", label=r"$g_{\rm cmp}=e^{\xi/2}$  (comparator)")
ax[0].plot([1 / 3], [(4 / 3) * np.sqrt(1 / 3)], "o", ms=8, color=C1, mec=SURFACE,
           mew=1.5, zorder=5)
ax[0].plot([0.0], [2 * 0.5 ** 1.5], "s", ms=7, color=C2, mec=SURFACE, mew=1.5, zorder=5)
ax[0].annotate(r"$\zeta=\frac{1}{3}$", (1 / 3, (4 / 3) * np.sqrt(1 / 3)),
               textcoords="offset points", xytext=(6, 8), fontsize=9, color=C1)
ax[0].annotate(r"apex $\zeta=0$", (0.0, 2 * 0.5 ** 1.5), textcoords="offset points",
               xytext=(-38, 8), fontsize=8.6, color=C2)
ax[0].set_xlabel(r"allocation $\zeta$")
ax[0].set_ylabel("factor")
ax[0].set_ylim(0, 2.1)
ax[0].set_title(r"$g(\zeta)=2h^{3/2}\,e^{\xi/2}$: unique even $\times$ odd split")
ax[0].legend(fontsize=8.4, loc="upper left")
ax[0].grid(True, alpha=0.6)

d_sym = -3 * z / (2 * (1 - z ** 2))
d_cmp = 1 / (2 * (1 - z ** 2))
ax[1].plot(z, d_sym, color=C2, ls="--", label=r"$d\log g_{\rm sym}/d\zeta$")
ax[1].plot(z, d_cmp, color=C3, ls="-.", label=r"$d\log g_{\rm cmp}/d\zeta$")
ax[1].plot(z, d_sym + d_cmp, color=C1, label=r"sum $=(1-3\zeta)/[2(1-\zeta^2)]$")
ax[1].axhline(0, color=MUTED, lw=0.9)
ax[1].axvline(1 / 3, color=C1, lw=1.2, ls=":")
ax[1].text(-0.12, 2.75, "magnitude pulls to the apex,\ncomparator pushes to rotation:\n"
           r"they cancel at $\zeta=1/3$", ha="center", fontsize=8.2, color=INK2)
ax[1].set_xlabel(r"allocation $\zeta$")
ax[1].set_ylabel("logarithmic derivative")
ax[1].set_ylim(-4, 4)
ax[1].set_title("The landmark is the balance of the two layers")
ax[1].legend(fontsize=8.4, loc="lower left")
ax[1].grid(True, alpha=0.6)

fig.text(0.005, -0.045, "EXACT (continuation, src/comparator_exact.py). The Thales "
         "altitude family sees only the dashed curve; the comparator supplies the "
         "dash-dotted one.", fontsize=6.6, color=MUTED, ha="left", va="top")
fig.savefig(FIG / "fig_comparator_split.png", bbox_inches="tight")
print(f"wrote {FIG / 'fig_comparator_split.png'}")
