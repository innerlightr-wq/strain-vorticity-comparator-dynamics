"""XI EVOLUTION AUDIT -- one figure: the size and locality of the three terms.

Panel (a): distributions of the local, pressure and viscous contributions to
           D xi/Dt in the DNS snapshot.
Panel (b): median |pressure| / |local| as a function of the allocation zeta,
           i.e. where in the strain-rotation balance the nonlocal obstruction
           matters most.

Palette: the project's validated reference instance (slots 1-3, light surface).

Run:  ../../../.venv/bin/python src/xi_figure.py
"""

from __future__ import annotations

import json
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

d = np.load(HERE / "results" / "xi_dns_terms.npz")
meta = json.loads((HERE / "results" / "xi_dns_pressure.json").read_text())
L, Pr, V, zeta = d["local"], d["pressure"], d["viscous"], d["zeta"]

fig, ax = plt.subplots(1, 2, figsize=(10.0, 3.9))

bins = np.linspace(-1.5, 1.5, 161)
for arr, col, lab in ((L, C1, r"local  $\|S\|[(1+\frac{1}{2}e^{2\xi})A-s/\sqrt{6}]$"),
                      (Pr, C2, r"pressure  $S{:}H_{\rm dev}/\|S\|^2$"),
                      (V, C3, r"viscous")):
    x = arr / np.std(L)
    inside = float(np.mean((x >= bins[0]) & (x <= bins[-1])))
    ax[0].hist(x, bins=bins, histtype="step", color=col, lw=2.0, density=True,
               label=lab + f"   [{inside:.0%} shown]")
ax[0].set_yscale("log")
ax[0].set_xlabel(r"contribution to $D\xi/Dt$, in units of $\mathrm{std(local)}$")
ax[0].set_ylabel("density")
ax[0].set_title("The three terms of the comparator equation")
ax[0].legend(fontsize=8.0, loc="upper right")
ax[0].grid(True, alpha=0.6)

edges = np.quantile(zeta, np.linspace(0, 1, 21))
mid, med, frac = [], [], []
for i in range(len(edges) - 1):
    m = (zeta >= edges[i]) & (zeta < edges[i + 1])
    if m.sum() > 100:
        mid.append(0.5 * (edges[i] + edges[i + 1]))
        med.append(np.median(np.abs(Pr[m]) / (np.abs(L[m]) + 1e-300)))
        frac.append(np.mean(np.abs(Pr[m]) > np.abs(L[m])))
ax[1].plot(mid, med, color=C2, marker="o", ms=5,
           label=r"median $|{\rm pressure}|/|{\rm local}|$")
ax[1].plot(mid, frac, color=C1, marker="s", ms=5, ls="--",
           label=r"fraction with $|{\rm pressure}|>|{\rm local}|$")
ax[1].axhline(1.0, color=MUTED, lw=0.9, ls=":")
ax[1].axvline(0.0, color=MUTED, lw=0.9)
ax[1].set_xlabel(r"allocation $\zeta$")
ax[1].set_ylabel("ratio / fraction")
ax[1].set_title("The nonlocal obstruction across the balance axis")
ax[1].set_ylim(0.4, 1.95)
ax[1].legend(fontsize=8.2, loc="upper right")
ax[1].grid(True, alpha=0.6)

fig.text(0.005, -0.045,
         f"EMPIRICAL (unforced Taylor-Green DNS, {meta['parameters']['n']}^3, "
         f"Re = {meta['parameters']['re']:.0f}, t = {meta['parameters']['t']}). "
         "The pressure term is the part restricted Euler discards.",
         fontsize=6.6, color=MUTED, ha="left", va="top")
fig.tight_layout()
fig.savefig(FIG / "fig_xi_terms.png", bbox_inches="tight")
print(f"wrote {FIG / 'fig_xi_terms.png'}")
