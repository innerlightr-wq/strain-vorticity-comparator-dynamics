"""Figures for the audit.  Reads results/*.json and results/*.npz only.

Palette: the validated reference instance of the dataviz skill.
  categorical slots 1-3 (all-pairs safe): blue #2a78d6, orange #eb6834, aqua #1baf7a
  diverging: blue <-> red with a neutral gray midpoint (signed production)
  ordinal (ordered zeta / ordered time): single blue hue, steps 250-700
  light chart surface #fcfcfb, ink #0b0b0b / #52514e
One y-axis per panel throughout; panels that share a quantity share its scale.
Static report figures are rendered light-surface only.

Run:  .venv/bin/python src/figures.py
"""

from __future__ import annotations

import json
import pathlib
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

# repo-root finder: works from src/ and from any experiments/ subdirectory
_ROOT = next(p for p in pathlib.Path(__file__).resolve().parents
             if (p / "src").is_dir())
sys.path.insert(0, str(_ROOT / "src"))
import svcore as sv  # noqa: E402

ROOT = _ROOT
RES = ROOT / "results"
FIG = ROOT / "figures"
FIG.mkdir(exist_ok=True)

SURFACE = "#fcfcfb"
INK, INK2, MUTED, GRID = "#0b0b0b", "#52514e", "#8a8985", "#e6e5e1"
C1, C2, C3 = "#2a78d6", "#eb6834", "#1baf7a"
DIVERGING = LinearSegmentedColormap.from_list(
    "bl_gy_rd", ["#104281", "#2a78d6", "#f0efec", "#e34948", "#8f1f1f"])
BLUE_ORDINAL = ["#86b6ef", "#6da7ec", "#3987e5", "#2a78d6",
                "#256abf", "#184f95", "#0d366b"]

plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
    "savefig.facecolor": SURFACE, "font.size": 9.5,
    "axes.edgecolor": MUTED, "axes.labelcolor": INK, "text.color": INK,
    "xtick.color": INK2, "ytick.color": INK2, "axes.linewidth": 0.8,
    "grid.color": GRID, "grid.linewidth": 0.7, "lines.linewidth": 2.0,
    "legend.frameon": False, "figure.dpi": 150, "axes.titlesize": 10,
    "axes.titleweight": "bold", "axes.spines.top": False,
    "axes.spines.right": False,
})


def finish(fig, path, note=None):
    if note:
        fig.text(0.005, -0.045, note, fontsize=6.6, color=MUTED, ha="left", va="top")
    fig.savefig(FIG / path, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote figures/{path}")


npz_ens = np.load(RES / "ensemble_samples.npz")
ens = lambda tag, key: npz_ens[f"{tag}__{key}"]
print("figures:")

# ---------------------------------------------------------------- figure 1
fig, ax = plt.subplots(1, 2, figsize=(9.6, 3.9))
th = np.linspace(0, np.pi, 400)
ax[0].plot(0.5 * np.cos(th), 0.5 * np.sin(th), color=INK2, lw=1.6)
ax[0].plot([-0.5, 0.5], [0, 0], color=MUTED, lw=0.9)
for z0, lab, col, dy in [(+1.0, "solid-body rotation", C1, 13),
                         (-1.0, "pure strain", C2, 13),
                         (0.0, "simple shear = Thales apex", C3, 11),
                         (1 / 3, r"$\zeta=1/3$", "#4a3aa7", -18)]:
    L, h = sv.L_of_zeta(z0), sv.h_of_zeta(z0)
    ax[0].plot([L], [h], "o", ms=8, color=col, mec=SURFACE, mew=1.5, zorder=5)
    ax[0].annotate(lab, (L, h), textcoords="offset points", xytext=(0, dy),
                   ha="center", fontsize=8, color=col)
ax[0].set_xlabel(r"$L=\zeta/2$   (lateral)")
ax[0].set_ylabel(r"$h=\sqrt{ab}$   (altitude)")
ax[0].set_title("Thales semicircle: the image of $\\zeta$")
ax[0].set_xlim(-0.74, 0.74)
ax[0].set_ylim(-0.08, 0.66)
ax[0].set_aspect("equal")
ax[0].grid(True, alpha=0.6)
ax[0].text(0.0, 0.17, "pure rotation and pure strain\nland on the SAME point",
           ha="center", fontsize=7.6, color=INK2)

zz = np.linspace(-1, 1, 600)
ax[1].plot(zz, sv.h_of_zeta(zz), color=C1, label=r"$h=\frac{1}{2}\sqrt{1-\zeta^2}$")
ax[1].plot(zz, sv.L_of_zeta(zz), color=C2, label=r"$L=\zeta/2$")
ax[1].plot(zz, sv.D_of_zeta(zz), color=C3,
           label=r"$D=\frac{1}{2}(1-\sqrt{1-\zeta^2})$")
ax[1].axvline(0, color=MUTED, lw=0.8, ls=":")
ax[1].set_ylim(-0.60, 0.80)
ax[1].set_xlabel(r"$\zeta=(\|\Omega\|_F^2-\|S\|_F^2)/\|\nabla u\|_F^2$")
ax[1].set_ylabel("Thales coordinate")
ax[1].set_title("Exact reparameterizations of $\\zeta$")
ax[1].legend(loc="upper center", ncol=3, fontsize=8, columnspacing=1.0,
             handlelength=1.4)
ax[1].grid(True, alpha=0.6)
ax[1].text(0.32, -0.34, "$h$ and $D$ are EVEN: 2-to-1 in $\\zeta$\n$L$ is a bijection",
           ha="center", fontsize=7.6, color=INK2)
finish(fig, "fig1_thales_map.png",
       "EXACT. Every Thales coordinate is a deterministic function of zeta alone; "
       "h and D additionally discard sign(zeta).")

# ---------------------------------------------------------------- figure 2
z_s, p_s = ens("E3_zeta_designed", "zeta"), ens("E3_zeta_designed", "p_norm")
fig, ax = plt.subplots(figsize=(6.8, 4.4))
sub = np.random.default_rng(0).choice(len(z_s), 25_000, replace=False)
ax.scatter(z_s[sub], p_s[sub], s=1.2, c=MUTED, alpha=0.18, lw=0)
zz = np.linspace(-1, 1, 800)
env = sv.production_envelope(zz)
ax.plot(zz, env, color=C1)
ax.plot(zz, -env, color=C1)
ax.axhline(0, color=MUTED, lw=0.8)
ax.axvline(1 / 3, color=C2, lw=1.3, ls="--")
ax.plot([1 / 3], [sv.SHARP_P_CONST], "o", ms=9, color=C2, mec=SURFACE, mew=1.5,
        zorder=6)
ax.annotate(r"$\zeta=\frac{1}{3}$:" "\n"
            r"$P/\|\nabla u\|_F^3 \leq \frac{4\sqrt{2}}{9}=0.6285$",
            (1 / 3, sv.SHARP_P_CONST), textcoords="offset points", xytext=(11, -4),
            ha="left", va="top", fontsize=8.4, color=C2)
ax.plot([0.0], [0.0], "s", ms=8, color=C3, mec=SURFACE, mew=1.5, zorder=6)
ax.annotate("Thales apex $\\zeta=0$\n(simple shear: $P=0$)", (0, 0),
            textcoords="offset points", xytext=(-8, -36), ha="center",
            fontsize=8.2, color=C3)
ax.set_xlim(-1.05, 1.42)
ax.set_xlabel(r"allocation $\zeta$")
ax.set_ylabel(r"normalized production $P/\|\nabla u\|_F^3$")
ax.set_title("At fixed allocation, production fills its whole envelope")
ax.legend(handles=[Line2D([], [], color=C1, lw=2,
                          label=r"exact envelope $\pm\sqrt{2/3}\,g(\zeta)$"),
                   Line2D([], [], color=MUTED, marker="o", ls="", ms=6,
                          label="random states (isotropic orientation)")],
          loc="lower left", fontsize=8.4)
ax.grid(True, alpha=0.6)
finish(fig, "fig2_production_envelope.png",
       "EXACT envelope (curve) + COMPUTATIONAL sample (points). zeta bounds |P| and "
       "determines nothing inside the bound.")

# ---------------------------------------------------------------- figure 3
npz_dns = np.load(RES / "dns_samples.npz")
dns_tag = sorted({f.split("_")[0] for f in npz_dns.files},
                 key=lambda t: float(t[1:]))[-1]
panels = [("gauss", "Gaussian solenoidal field (no alignment preference)",
           ens("E1_gaussian_solenoidal", "zeta"), ens("E1_gaussian_solenoidal", "A"),
           ens("E1_gaussian_solenoidal", "p_norm")),
          ("dns", f"Navier-Stokes DNS, Taylor-Green ({dns_tag})",
           npz_dns[f"{dns_tag}_zeta"], npz_dns[f"{dns_tag}_A"],
           npz_dns[f"{dns_tag}_P"] / npz_dns[f"{dns_tag}_Q"] ** 1.5)]
vmax = max(float(np.nanpercentile(np.abs(p), 99)) for _, _, _, _, p in panels)
fig, axs = plt.subplots(1, 2, figsize=(10.6, 4.1))
for k, (_, title, zt, At, pt) in enumerate(panels):
    a = axs[k]
    idx = np.random.default_rng(1).choice(len(zt), min(40_000, len(zt)), replace=False)
    sc = a.scatter(zt[idx], At[idx], c=pt[idx], s=1.6, lw=0, cmap=DIVERGING,
                   norm=TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax))
    for sgn, va in ((+1, "bottom"), (-1, "top")):
        a.axhline(sgn * sv.SQRT_2_3, color=INK2, lw=1.2, ls="--")
        a.text(-0.97, sgn * sv.SQRT_2_3 + sgn * 0.02,
               f"$A={'+' if sgn > 0 else '-'}\\sqrt{{2/3}}$ (sharp)",
               fontsize=7.8, color=INK2, va=va)
    a.set_xlim(-1.02, 1.02)
    a.set_ylim(-0.95, 0.95)
    a.set_xlabel(r"allocation $\zeta$")
    a.set_ylabel(r"alignment $A$" if k == 0 else "")
    a.set_title(title, fontsize=9.2)
    a.grid(True, alpha=0.5)
cb = fig.colorbar(sc, ax=axs, pad=0.015, fraction=0.035)
cb.set_label(r"$P/\|\nabla u\|_F^3$  (shared scale)", fontsize=8.5)
cb.outline.set_visible(False)
finish(fig, "fig3_joint_state_space.png",
       "The admissible set is the full rectangle (EXACT); which part of it a flow "
       "occupies is dynamics, not geometry. Both panels share one color scale.")

# ---------------------------------------------------------------- figure 4
cf = json.loads((RES / "canonical_flows.json").read_text())
B = cf["burgers"]
a_vals = np.array(B["sweep_fixed_r"]["1.0"]["a"])
fig, ax = plt.subplots(2, 2, figsize=(9.8, 6.4))

ax[0, 0].plot(a_vals, B["sweep_fixed_r"]["1.0"]["zeta"], color=C1, label=r"$\zeta$")
ax[0, 0].plot(a_vals, B["sweep_fixed_r"]["1.0"]["A"], color=C2, label=r"$A$")
ax[0, 0].axhline(sv.SQRT_2_3, color=MUTED, lw=0.9, ls=":")
ax[0, 0].text(1.2e-2, sv.SQRT_2_3 + 0.05, r"$\sqrt{2/3}$", fontsize=8, color=INK2)
ax[0, 0].set_xscale("log")
ax[0, 0].set_ylim(-1.3, 1.1)
ax[0, 0].set_xlabel("axial strain $a$")
ax[0, 0].set_ylabel("coordinate")
ax[0, 0].set_title("fixed physical radius $r=1$", fontsize=9.4)
ax[0, 0].legend(fontsize=8.4, loc="lower left", ncol=2)
ax[0, 0].grid(True, alpha=0.6)

ax[0, 1].plot(a_vals, B["sweep_fixed_eta"]["1.0"]["zeta"], color=C1, label=r"$\zeta$")
ax[0, 1].plot(a_vals, B["sweep_fixed_eta"]["1.0"]["A"], color=C2, label=r"$A$")
ax[0, 1].set_xscale("log")
ax[0, 1].set_ylim(-1.3, 1.1)
ax[0, 1].set_xlabel("axial strain $a$")
ax[0, 1].set_title(r"fixed similarity radius $\eta=1$", fontsize=9.4)
ax[0, 1].legend(fontsize=8.4, loc="lower left", ncol=2)
ax[0, 1].grid(True, alpha=0.6)
ax[0, 1].text(1.0, -0.66, "exactly constant in $a$", ha="center", fontsize=8.4,
              color=INK2)

Pr = np.array(B["sweep_fixed_r"]["1.0"]["P"])
Pe = np.array(B["sweep_fixed_eta"]["1.0"]["P"])
ax[1, 0].plot(a_vals, Pe, color=C1, label=r"fixed $\eta=1$   ($P\propto a^3$)")
ax[1, 0].plot(a_vals, np.maximum(Pr, 1e-16), color=C2, label=r"fixed $r=1$")
ax[1, 0].set_xscale("log")
ax[1, 0].set_yscale("log")
ax[1, 0].set_ylim(1e-8, 1e8)
ax[1, 0].set_xlabel("axial strain $a$")
ax[1, 0].set_ylabel(r"production $P$")
ax[1, 0].set_title("production over the same two sweeps", fontsize=9.4)
ax[1, 0].legend(fontsize=8.4, loc="upper left")
ax[1, 0].grid(True, alpha=0.6, which="major")

prof = B["profile_at_a_1"]
ax[1, 1].plot(prof["eta"], prof["zeta"], color=C1, label=r"$\zeta(\eta)$")
ax[1, 1].plot(prof["eta"], prof["A"], color=C2, label=r"$A(\eta)$")
ax[1, 1].axhline(sv.SQRT_2_3, color=MUTED, lw=0.9, ls=":")
ax[1, 1].text(3.4, sv.SQRT_2_3 + 0.06, r"$\sqrt{2/3}$ (bound)", fontsize=8, color=INK2)
ax[1, 1].set_ylim(-1.3, 1.1)
ax[1, 1].set_xlabel(r"similarity radius $\eta=r/\sqrt{4\nu/a}$")
ax[1, 1].set_title("radial profile: the core saturates the bound", fontsize=9.4)
ax[1, 1].legend(fontsize=8.4, loc="lower left", ncol=2)
ax[1, 1].grid(True, alpha=0.6)
fig.suptitle("Burgers vortex: allocation, alignment and production "
             r"($Re_\Gamma=100$)", fontsize=10.5, fontweight="bold")
fig.tight_layout(rect=(0, 0, 1, 0.96))
finish(fig, "fig4_burgers.png",
       "EXACT symbolic sweeps. On the axis A = sqrt(2/3) exactly for every a, Gamma, "
       "nu; at fixed eta both coordinates are a-independent while P spans 12 decades.")

# ---------------------------------------------------------------- figure 5
re = json.loads((RES / "restricted_euler.json").read_text())
npz_re = np.load(RES / "restricted_euler_samples.npz")
z0, A0, amp = npz_re["zeta0"], npz_re["A0"], npz_re["amplification"]
fig, ax = plt.subplots(1, 2, figsize=(10.0, 3.9))
for i, zv in enumerate(np.unique(z0)):
    m = z0 == zv
    Au = np.unique(A0[m])
    ax[0].plot(Au, [amp[m][A0[m] == av].mean() for av in Au],
               color=BLUE_ORDINAL[i % len(BLUE_ORDINAL)], marker="o", ms=4.5,
               label=f"$\\zeta={zv:+.2f}$")
ax[0].axhline(0, color=MUTED, lw=0.8)
ax[0].set_xlabel(r"initial alignment $A$")
ax[0].set_ylabel(r"$\ln\,[e_\omega(T)/e_\omega(0)]$")
ax[0].set_title("Amplification at matched $\\zeta$ and scale", fontsize=9.6)
ax[0].legend(fontsize=7.6, ncol=2, loc="upper left")
ax[0].grid(True, alpha=0.6)

sens = re["A_sensitivity_vs_exact_law"]
zs = np.array([s["zeta"] for s in sens])
sl = np.array([s["slope_dAmp_dA"] for s in sens])
pred = np.array([s["predicted_shape"] for s in sens])
ax[1].plot(zs, sl, color=C1, marker="o", ms=6,
           label=r"measured $\partial$amp$/\partial A$")
ax[1].plot(zs, pred * np.mean(sl / pred), color=C2, ls="--",
           label=r"exact law $\propto\sqrt{2(1-\zeta)}$")
ax[1].set_xlabel(r"allocation $\zeta$")
ax[1].set_ylabel("sensitivity to alignment")
ax[1].set_title("Sensitivity follows the exact allocation law", fontsize=9.6)
ax[1].legend(fontsize=8.4)
ax[1].grid(True, alpha=0.6)
finish(fig, "fig5_restricted_euler.png",
       "COMPUTATIONAL (restricted Euler: no pressure Hessian, no viscosity). "
       f"Scatter of slope/law across zeta: {re['A_sensitivity_ratio_cv']:.1%}.")

# ---------------------------------------------------------------- figure 6
rj = json.loads((RES / "dns_reanalysis.json").read_text())
dj = json.loads((RES / "dns_taylor_green.json").read_text())
times = sorted([t for t in rj if t.startswith("t") and t[1:].replace(".", "").isdigit()],
               key=lambda t: float(t[1:]))
fig, ax = plt.subplots(1, 3, figsize=(13.4, 3.9))
last = times[-1]
rows = rj[last]["by_A_quartile"]
xs = np.arange(4)
w = 0.26
ax[0].bar(xs - w, [r["mean_G_P"] for r in rows], w, color=C1, label=r"production $P$")
ax[0].bar(xs, [r["mean_G_visc"] for r in rows], w, color=C2,
          label=r"viscous $-\epsilon_\omega+T_\nu$")
ax[0].bar(xs + w, [r["mean_G_tot"] for r in rows], w, color=C3,
          label=r"net $D_t e_\omega$")
ax[0].axhline(0, color=MUTED, lw=0.9)
ax[0].set_xticks(xs)
ax[0].set_xticklabels([f"Q{q+1}\n$A$={r['mean_A']:+.2f}" for q, r in enumerate(rows)],
                      fontsize=8)
ax[0].set_ylabel(r"term $/\,\|\nabla u\|_F^3$")
ax[0].set_title(f"DNS {last}: balance by $A$-quartile\nwithin matched "
                r"$(\zeta,$ scale$)$ cells", fontsize=9.4)
ax[0].legend(fontsize=8.2, loc="lower left")
ax[0].grid(True, alpha=0.5, axis="y")

for i, t in enumerate(times):
    fr = [q["frac_Dt_e_positive"] for q in rj[t]["by_A_quartile"]]
    ax[1].plot(xs, fr, color=BLUE_ORDINAL[2 + i], marker="o", ms=6,
               label=f"$t={t[1:]}$")
ax[1].axhline(0.5, color=MUTED, lw=0.9, ls=":")
ax[1].set_xticks(xs)
ax[1].set_xticklabels([f"Q{q+1}" for q in range(4)])
ax[1].set_xlabel(r"$A$-quartile within matched $(\zeta,$ scale$)$ cells")
ax[1].set_ylabel(r"fraction with $D_t e_\omega>0$")
ax[1].set_title("Probability of net enstrophy growth", fontsize=9.4)
ax[1].legend(fontsize=8.2, ncol=2)
ax[1].grid(True, alpha=0.5)

lag = dj.get("lagrangian") or {}
if lag:
    lr = lag["by_A_quartile"]
    ax[2].bar(xs, [r["mean_amplification"] for r in lr], 0.55, color=C1)
    for q, r in enumerate(lr):
        v = r["mean_amplification"]
        ax[2].annotate(f"{v:+.3f}", (q, v), ha="center",
                       va="bottom" if v >= 0 else "top", fontsize=8, color=INK2,
                       textcoords="offset points", xytext=(0, 4 if v >= 0 else -4))
    ax[2].axhline(0, color=MUTED, lw=0.9)
    ax[2].set_xticks(xs)
    ax[2].set_xticklabels([f"Q{q+1}\n$A$={r['mean_A']:+.2f}"
                           for q, r in enumerate(lr)], fontsize=8)
    ax[2].set_ylabel(r"$\ln\,[e_\omega(t_1)/e_\omega(t_0)]$")
    ax[2].set_title(f"Lagrangian amplification, $t={lag['t_seed']:.1f}"
                    f"\\rightarrow{lag['t_final']:.1f}$", fontsize=9.4)
    ax[2].grid(True, alpha=0.5, axis="y")
fig.tight_layout()
finish(fig, "fig6_dns_conditional.png",
       "EMPIRICAL (DNS, unforced Taylor-Green, 96^3, Re=400). Production and net rate "
       "are reported separately: P > 0 is not net growth.")

# ---------------------------------------------------------------- figure 7
en = json.loads((RES / "ensembles.json").read_text())
scan = en["phase6_leakage_resolution_scan"]
fig, ax = plt.subplots(figsize=(6.6, 4.2))
nb = np.array([s["n_zeta_bins"] for s in scan], float)
ax.plot(nb, np.abs([s["partial_corr_h"] for s in scan]), color=C1, marker="o", ms=6,
        label=r"$|\mathrm{partial\ corr}(h,\,|p|\,\mid\,\zeta)|$")
ax.plot(nb, np.abs([s["partial_corr_A"] for s in scan]), color=C2, marker="s", ms=6,
        label=r"$|\mathrm{partial\ corr}(|A|,\,|p|\,\mid\,\zeta)|$")
ax.plot(nb, 1.2 / nb, color=MUTED, ls=":", lw=1.4, label=r"$\propto 1/n_{\rm bins}$")
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel(r"number of $\zeta$ strata used for conditioning")
ax.set_ylabel(r"partial correlation with $|P|/\|\nabla u\|_F^3$")
ax.set_title("The Thales altitude's apparent signal is discretization leakage")
ax.legend(fontsize=8.4, loc="center left")
ax.grid(True, alpha=0.6, which="both")
finish(fig, "fig7_leakage_scan.png",
       "CONTROL. h's partial correlation decays as the conditioning is refined "
       "(it is an exact function of zeta); A's does not.")

print("done")
