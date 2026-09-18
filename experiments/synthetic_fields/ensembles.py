"""PHASE 5, 6, 9 -- synthetic ensembles, the Thales test, adversarial controls.

Ensembles
  E1  Gaussian isotropic solenoidal field: gradients read off a divergence-free
      Gaussian random field (the natural measure on "random incompressible local
      velocity-gradient states").
  E2  Naive independent construction (as suggested in the brief): random
      traceless symmetric S, independent isotropic omega with independent
      magnitude -- sweeps zeta much more widely than E1.
  E3  zeta-designed: zeta drawn uniformly on (-1,1), orientation isotropic.
      Used to answer "at fixed zeta, how much variation remains?" without the
      confound of an ensemble-specific zeta distribution.
  E4  Preferentially aligned: omega tilted towards a chosen strain eigenvector.

Every correlation reported here is checked against the exact algebra of
DERIVATIONS.md; quantities that are deterministic functions of one another are
labelled as such and are NOT reported as empirical findings.

Run:  .venv/bin/python src/ensembles.py
"""

from __future__ import annotations

import json
import pathlib
import sys

import numpy as np

# repo-root finder: works from src/ and from any experiments/ subdirectory
_ROOT = next(p for p in pathlib.Path(__file__).resolve().parents
             if (p / "src").is_dir())
sys.path.insert(0, str(_ROOT / "src"))
import statlib as st  # noqa: E402
import svcore as sv  # noqa: E402

OUT = _ROOT / "results"
OUT.mkdir(exist_ok=True)
FAILURES: list[str] = []
report: dict[str, object] = {}


def require(name, ok, detail=""):
    print(f"[{'PASS' if ok else 'FAIL':4s}] {name}" + (f"   {detail}" if detail else ""))
    if not ok:
        FAILURES.append(name)


# generators live in src/generators.py so that other scripts can import them
# without executing this analysis
from generators import (gaussian_solenoidal_gradients, ensemble_naive,  # noqa: E402
                        ensemble_zeta_designed, ensemble_aligned, isotropic_dirs,
                        random_traceless_symmetric)

# ---------------------------------------------------------------- build data
print("=" * 72)
print("PHASE 5: synthetic ensembles")
print("=" * 72)

G1 = gaussian_solenoidal_gradients(n_grid=64, n_real=4, seed=11)
c1 = sv.coordinates_from_gradient(G1)
S1 = sv.sym(G1)
w1 = sv.vorticity_from_skew(sv.skew(G1))
print(f"E1 Gaussian solenoidal: n = {len(G1)} gradient samples")
require("E1 is incompressible to round-off",
        np.max(np.abs(np.trace(S1, axis1=-2, axis2=-1))) < 1e-9,
        f"max |div u| = {np.max(np.abs(np.trace(S1, axis1=-2, axis2=-1))):.2e}")

S2, w2 = ensemble_naive(400_000, seed=12)
c2 = sv.coordinates(S2, w2)
S3, w3, zeta_design = ensemble_zeta_designed(400_000, seed=13)
c3 = sv.coordinates(S3, w3)
require("E3 realizes the designed zeta exactly",
        np.max(np.abs(c3["zeta"] - zeta_design)) < 1e-12,
        f"max error {np.max(np.abs(c3['zeta'] - zeta_design)):.2e}")
S4, w4 = ensemble_aligned(400_000, seed=14, which=2, kappa=4.0)
c4 = sv.coordinates(S4, w4)
S4b, w4b = ensemble_aligned(400_000, seed=15, which=1, kappa=4.0)
c4b = sv.coordinates(S4b, w4b)

ens = {"E1_gaussian_solenoidal": (S1, w1, c1),
       "E2_naive_independent": (S2, w2, c2),
       "E3_zeta_designed": (S3, w3, c3),
       "E4_aligned_max": (S4, w4, c4),
       "E4_aligned_mid": (S4b, w4b, c4b)}

# -------------------------------------------------- bounds hold in every sample
print("\n--- exact bounds, checked on every sample ---")
bounds = {}
for tag, (S, w, c) in ens.items():
    good = np.isfinite(c["A"])
    maxA = float(np.nanmax(np.abs(c["A"][good])))
    maxp = float(np.nanmax(np.abs(c["p_norm"][good])))
    env_ok = bool(np.all(np.abs(c["p_norm"][good]) <=
                         sv.production_envelope(c["zeta"][good]) + 1e-10))
    bounds[tag] = {"n": int(good.sum()), "max_abs_A": maxA, "max_abs_p_norm": maxp,
                   "envelope_respected": env_ok,
                   "frac_at_bound_1pct": float(np.mean(np.abs(c["A"][good]) >
                                                       0.99 * sv.SQRT_2_3))}
    require(f"{tag}: |A| <= sqrt(2/3) and |p| <= envelope(zeta)",
            maxA <= sv.SQRT_2_3 + 1e-10 and env_ok and maxp <= sv.SHARP_P_CONST + 1e-10,
            f"max|A| = {maxA:.9f}, max|p| = {maxp:.9f} (sharp: {sv.SHARP_P_CONST:.9f})")

# ------------------------------------------- Phase 5 core question: at fixed zeta
print("\n--- at fixed zeta, how much variation remains in A and P? ---")


def conditional_spread(c, nz=20):
    z = c["zeta"]
    A = c["A"]
    p = c["p_norm"]
    m = np.isfinite(A) & np.isfinite(p) & np.isfinite(z)
    z, A, p = z[m], A[m], p[m]
    idx, edges = st.quantile_bins(z, nz)
    rows = []
    for b in range(idx.max() + 1):
        sel = idx == b
        if sel.sum() < 100:
            continue
        rows.append(dict(zeta_lo=float(edges[b]), zeta_hi=float(edges[b + 1]),
                         n=int(sel.sum()), A_mean=float(A[sel].mean()),
                         A_std=float(A[sel].std()),
                         A_min=float(A[sel].min()), A_max=float(A[sel].max()),
                         p_std=float(p[sel].std()),
                         frac_P_pos=float((p[sel] > 0).mean())))
    total = dict(A_std_total=float(A.std()), p_std_total=float(p.std()),
                 frac_P_pos_total=float((p > 0).mean()))
    within_A = float(np.sqrt(np.average([r["A_std"] ** 2 for r in rows],
                                        weights=[r["n"] for r in rows])))
    return dict(bins=rows, totals=total, within_bin_A_std=within_A,
                variance_of_A_explained_by_zeta=float(
                    1.0 - within_A ** 2 / A.std() ** 2))


cond = {}
for tag, (S, w, c) in ens.items():
    cond[tag] = conditional_spread(c)
    t = cond[tag]
    print(f"  {tag}: std(A) = {t['totals']['A_std_total']:.4f}, "
          f"within-zeta-bin std(A) = {t['within_bin_A_std']:.4f}  "
          f"=> zeta explains {t['variance_of_A_explained_by_zeta']:.2%} of var(A)")

require("E3: conditioning on zeta removes essentially none of the variance of A",
        abs(cond["E3_zeta_designed"]["variance_of_A_explained_by_zeta"]) < 0.01,
        f"{cond['E3_zeta_designed']['variance_of_A_explained_by_zeta']:.4%} "
        "(designed ensemble: zeta and orientation are independent by construction)")
require("E1: conditioning on zeta removes little of the variance of A in a "
        "physically generated Gaussian field",
        abs(cond["E1_gaussian_solenoidal"]["variance_of_A_explained_by_zeta"]) < 0.10,
        f"{cond['E1_gaussian_solenoidal']['variance_of_A_explained_by_zeta']:.4%}")

# sign of P at fixed zeta
print("\n--- sign of P: what can zeta predict? ---")
sign_tests = {}
for tag, (S, w, c) in ens.items():
    m = np.isfinite(c["A"])
    z, p = c["zeta"][m], c["p_norm"][m]
    idx, edges = st.quantile_bins(z, 20)
    base = max((p > 0).mean(), (p <= 0).mean())
    acc = 0.0
    for b in range(idx.max() + 1):
        sel = idx == b
        if sel.sum() == 0:
            continue
        frac = (p[sel] > 0).mean()
        acc += sel.mean() * max(frac, 1.0 - frac)
    sign_tests[tag] = {"majority_baseline": float(base),
                       "best_accuracy_from_zeta_bins": float(acc),
                       "gain_over_baseline": float(acc - base),
                       "accuracy_from_sign_A": 1.0}
    print(f"  {tag}: majority baseline {base:.4f}, best from 20 zeta bins {acc:.4f} "
          f"(gain {acc - base:+.4f}); sign(A) gives 1.0000 by identity")

require("zeta gives at most a marginal gain over the majority baseline for sign(P)",
        max(v["gain_over_baseline"] for v in sign_tests.values()) < 0.05,
        "sign(P) = sign(A) is an identity, not an empirical result")

# ------------------------------------------- Phase 6: the Thales geometry test
print("\n" + "=" * 72)
print("PHASE 6: Thales geometry test")
print("=" * 72)
thales = {}
for tag, (S, w, c) in ens.items():
    m = np.isfinite(c["A"])
    z, h, L, D = c["zeta"][m], c["h"][m], c["L"][m], c["D"][m]
    thales[tag] = {
        "max_err_L_minus_zeta_over_2": float(np.max(np.abs(L - z / 2))),
        "max_err_h_minus_formula": float(np.max(np.abs(h - sv.h_of_zeta(z)))),
        "max_err_L2_minus_D_1_minus_D": float(np.max(np.abs(L ** 2 - D * (1 - D)))),
        "max_err_L2_plus_h2_minus_quarter": float(np.max(np.abs(L ** 2 + h ** 2 - 0.25))),
    }
require("Thales identities hold to machine precision in every ensemble",
        max(max(v.values()) for v in thales.values()) < 1e-12,
        f"worst residual {max(max(v.values()) for v in thales.values()):.2e}")

# h is a deterministic 2-to-1 function of zeta: R^2 of h on a flexible f(zeta) = 1
det = {}
for tag, (S, w, c) in ens.items():
    m = np.isfinite(c["A"])
    z, h = c["zeta"][m], c["h"][m]
    idxz, _ = st.quantile_bins(z, 40)
    r2_bins = 1.0 - float(np.sum(st.group_mean_residual(h, idxz) ** 2)) / \
        float(np.sum((h - h.mean()) ** 2))
    det[tag] = {"r2_h_on_40_zeta_bins": r2_bins,
                "r2_h_on_exact_formula": st.r2_of(h, sv.h_of_zeta(z))}
    print(f"  {tag}: R^2(h | 40 zeta bins) = {r2_bins:.8f}, "
          f"R^2(h | (1/2)sqrt(1-zeta^2)) = {det[tag]['r2_h_on_exact_formula']:.10f}")
require("h is a deterministic function of zeta (R^2 = 1 against the exact formula)",
        all(v["r2_h_on_exact_formula"] > 1 - 1e-12 for v in det.values()))

# does h/L/D predict A or P beyond zeta?
print("\n--- conditional information: does h, L or D add anything beyond zeta? ---")
cmi = {}
for tag in ("E1_gaussian_solenoidal", "E3_zeta_designed", "E4_aligned_mid"):
    S, w, c = ens[tag]
    m = np.isfinite(c["A"])
    sub = slice(None, 200_000)
    z = c["zeta"][m][sub]
    A = c["A"][m][sub]
    p = c["p_norm"][m][sub]
    h = c["h"][m][sub]
    D = c["D"][m][sub]
    cmi[tag] = {
        "I(p ; A | zeta)": st.cmi_with_null(A, p, z, seed=1),
        "I(p ; h | zeta)": st.cmi_with_null(h, p, z, seed=2),
        "I(p ; D | zeta)": st.cmi_with_null(D, p, z, seed=3),
        "I(A ; h | zeta)": st.cmi_with_null(h, A, z, seed=4),
    }
    for k, v in cmi[tag].items():
        print(f"  {tag}: {k} = {v['cmi']:.5f} nats, shuffled null "
              f"{v['null_mean']:.5f} +- {v['null_std']:.5f}, excess {v['excess']:+.5f} "
              f"(z = {v['z_score']:+.1f})")

for tag, v in cmi.items():
    require(f"{tag}: I(p ; h | zeta) is two orders of magnitude below I(p ; A | zeta)",
            abs(v["I(p ; h | zeta)"]["excess"]) < 0.05 * v["I(p ; A | zeta)"]["excess"],
            f"excess {v['I(p ; h | zeta)']['excess']:+.5f} nats vs "
            f"I(p ; A | zeta) excess {v['I(p ; A | zeta)']['excess']:+.5f} nats")

# The residual I(p;h|zeta) above is NOT information: h = (1/2)sqrt(1-zeta^2) is an
# exact function of zeta, so the true conditional information is exactly zero and any
# nonzero estimate is leakage from conditioning on FINITE zeta bins (h resolves zeta
# inside a bin).  The signature of leakage is that it vanishes as the conditioning is
# refined, while genuine information does not.  Scan it.
print("\n--- resolution scan: leakage vs information under refined conditioning ---")
S, w, c = ens["E3_zeta_designed"]
m = np.isfinite(c["A"])
zz, hh, pp, AA = (c[k][m][:300_000] for k in ("zeta", "h", "p_norm", "A"))
scan = []
for nz in (8, 16, 32, 64, 128, 256, 512):
    # bias-corrected: subtract the within-stratum-shuffled null, which is the
    # plug-in estimator's own floor and grows as the strata get small
    ex_h = st.cmi_with_null(hh, pp, zz, nb=8, nz=nz, n_shuffle=4, seed=101)["excess"]
    ex_A = st.cmi_with_null(AA, pp, zz, nb=8, nz=nz, n_shuffle=4, seed=102)["excess"]
    bh = st.binned_partial(hh, np.abs(pp), zz, nz)
    bA = st.binned_partial(np.abs(AA), np.abs(pp), zz, nz)
    ph = st.binned_partial(hh, pp, zz, nz)
    pA = st.binned_partial(AA, pp, zz, nz)
    pc_h, pc_A = bh["partial_corr"], bA["partial_corr"]
    dr2_h, dr2_A = ph["delta_r2"], pA["delta_r2"]
    scan.append(dict(n_zeta_bins=nz, cmi_h=ex_h, cmi_A=ex_A,
                     partial_corr_h=pc_h, partial_corr_A=pc_A,
                     delta_r2_h=dr2_h, delta_r2_A=dr2_A))
    print(f"  nz = {nz:4d}:  I(p;h|z) = {ex_h:.5f}  I(p;A|z) = {ex_A:.5f} | "
          f"partial corr h = {pc_h:+.2e}  A = {pc_A:+.4f} | "
          f"dR2 h = {dr2_h:.2e}  A = {dr2_A:.4f}")
require("h's apparent conditional signal is discretization leakage: it decays with "
        "conditioning resolution while A's does not",
        scan[-1]["cmi_h"] < 0.35 * scan[0]["cmi_h"] and
        abs(scan[-1]["partial_corr_h"]) < 0.3 * abs(scan[0]["partial_corr_h"]) and
        scan[-1]["cmi_A"] > 0.8 * scan[0]["cmi_A"] and
        abs(scan[-1]["partial_corr_A"]) > 0.8 * abs(scan[0]["partial_corr_A"]),
        f"I(p;h|z): {scan[0]['cmi_h']:.4f} -> {scan[-1]['cmi_h']:.4f}; "
        f"I(p;A|z): {scan[0]['cmi_A']:.4f} -> {scan[-1]['cmi_A']:.4f}")

# ------------------------------------- the Appendix-B residual-altitude test trap
print("\n--- control: the residual-altitude test applied to a target with"
      " provably zero independent h information ---")
S, w, c = ens["E3_zeta_designed"]
m = np.isfinite(c["A"])
a_frac = c["a"][m]
z = c["zeta"][m]
h = c["h"][m]
# target built to be an exact function of the allocation alone, plus noise that is
# independent of everything: by construction h carries NO information beyond a.
rng = np.random.default_rng(99)
target = sv.allocation_factor(z) + 0.05 * rng.normal(size=len(z))
h_lin = st.ols(st.design(a_frac), h)
h_resid = h - st.design(a_frac) @ h_lin["beta"]
print(f"  linear fit of h on a: R^2 = {h_lin['r2']:.4f} "
      f"(so h_resid is the nonlinear part of h)")
trap_linear = st.nested_r2(target, [a_frac], [h_resid])
trap_flexible = st.binned_partial(h_resid, target, a_frac, 40)
print(f"  Appendix-B procedure (linear control in a): delta_R2 = "
      f"{trap_linear['delta_r2']:.4f}, t = {trap_linear['t_extra'][0]:.1f}, "
      f"p = {trap_linear['p_extra'][0]:.3e}  -> would be declared 'independent info'")
print(f"  same test with a FLEXIBLE control in a (40 bins): delta_R2 = "
      f"{trap_flexible['delta_r2']:.3e}, t = {trap_flexible['t']:.2f}, "
      f"p = {trap_flexible['p']:.3f}")
trap_scan = []
for nb in (2, 5, 10, 20, 40, 80, 160, 320, 640):
    tr = st.binned_partial(h_resid, target, a_frac, nb)
    trap_scan.append({"n_bins": nb, "delta_r2": tr["delta_r2"],
                      "t": tr["t"], "p": tr["p"]})
    print(f"    control = {nb:3d} bins in a: delta_R2 = {tr['delta_r2']:.3e}, "
          f"t = {tr['t']:8.2f}")
require("FALSIFIER: the residual-altitude test passes on a target with exactly "
        "zero independent h content when the control in a is linear",
        trap_linear["delta_r2"] > 0.01 and trap_linear["p_extra"][0] < 0.01,
        "the test detects nonlinearity in the target, not independent information in h")
require("the same delta_R2 collapses towards zero as the control in a is made "
        "flexible, confirming there was no information to find",
        trap_scan[-1]["delta_r2"] < 0.01 * trap_linear["delta_r2"],
        f"delta_R2: {trap_linear['delta_r2']:.4f} (linear control) -> "
        f"{trap_scan[-1]['delta_r2']:.2e} (320-bin control)")

# ------------------------------------------------------- Phase 9 controls
print("\n" + "=" * 72)
print("PHASE 9: adversarial controls")
print("=" * 72)
rng = np.random.default_rng(7)
controls: dict[str, object] = {}

# 9.1 / 9.2 rotate omega in a fixed strain eigenframe (all magnitudes fixed)
S0 = np.diag([2.0, -1.0, -1.0]) / np.sqrt(6.0)
w0mag = 2.0
dirs = isotropic_dirs(200_000, rng)
cR = sv.coordinates(np.broadcast_to(S0, (len(dirs), 3, 3)), dirs * w0mag)
controls["9.1_rotate_omega"] = {
    "zeta_spread": float(np.ptp(cR["zeta"])), "h_spread": float(np.ptp(cR["h"])),
    "A_mean": float(cR["A"].mean()), "A_std": float(cR["A"].std()),
    "A_min": float(cR["A"].min()), "A_max": float(cR["A"].max()),
    "P_min": float(cR["P"].min()), "P_max": float(cR["P"].max()),
    "frac_P_negative": float((cR["P"] < 0).mean()),
}
print(f"  9.1 isotropic rotation of omega at fixed magnitudes: zeta spread "
      f"{np.ptp(cR['zeta']):.2e}, A in [{cR['A'].min():+.4f}, {cR['A'].max():+.4f}], "
      f"mean A = {cR['A'].mean():+.6f}, P < 0 for {(cR['P'] < 0).mean():.1%} of orientations")
require("9.1 rotating omega leaves every magnitude coordinate exactly fixed",
        np.ptp(cR["zeta"]) < 1e-12 and np.ptp(cR["h"]) < 1e-12)
require("9.1 EXACT: isotropic orientation gives mean A = 0 for any traceless S",
        abs(cR["A"].mean()) < 5e-3,
        f"sample mean {cR['A'].mean():+.2e}; exactly 0 because E[cos^2] = 1/3 and "
        "tr S = 0 -- so nonzero mean production is purely an alignment effect")

# 9.3 preserve A while changing magnitude
scales = np.exp(rng.normal(0, 1.5, size=50_000))
Sc = S0[None] * scales[:, None, None]
wc = (np.array([1.0, 0.3, 0.2]) / np.linalg.norm([1.0, 0.3, 0.2]))[None] * \
     (w0mag * scales[:, None])
cS = sv.coordinates(Sc, wc)
controls["9.3_preserve_A_change_scale"] = {
    "A_spread": float(np.ptp(cS["A"])), "zeta_spread": float(np.ptp(cS["zeta"])),
    "P_ratio_max_over_min": float(cS["P"].max() / cS["P"].min()),
    "p_norm_spread": float(np.ptp(cS["p_norm"])),
}
require("9.3 scaling S and omega together fixes A, zeta and p_norm while P moves "
        "over orders of magnitude",
        np.ptp(cS["A"]) < 1e-12 and np.ptp(cS["zeta"]) < 1e-12 and
        np.ptp(cS["p_norm"]) < 1e-12,
        f"P spans a factor {cS['P'].max() / cS['P'].min():.3e}")

# 9.4 randomize eigenvectors, keep eigenvalues
lam_fixed = np.array([-0.4082, -0.4082, 0.8165])
Qs = np.linalg.qr(rng.normal(size=(50_000, 3, 3)))[0]
Se = Qs @ (lam_fixed[None, :, None] * np.swapaxes(Qs, -1, -2))
w_fixed = np.array([0.0, 0.0, 2.0])
cE = sv.coordinates(Se, np.broadcast_to(w_fixed, (len(Se), 3)))
controls["9.4_randomize_eigenvectors"] = {
    "zeta_spread": float(np.ptp(cE["zeta"])), "A_min": float(cE["A"].min()),
    "A_max": float(cE["A"].max()), "A_mean": float(cE["A"].mean()),
    "frac_P_negative": float((cE["P"] < 0).mean())}
require("9.4 randomizing strain eigenvectors at fixed eigenvalues leaves zeta fixed "
        "and moves A over its full range",
        np.ptp(cE["zeta"]) < 1e-10 and np.ptp(cE["A"]) > 1.0,
        f"A in [{cE['A'].min():+.4f}, {cE['A'].max():+.4f}], "
        f"P < 0 for {(cE['P'] < 0).mean():.1%}")

# 9.5 isotropic vs preferentially aligned
iso_A = cond["E3_zeta_designed"]["totals"]
controls["9.5_isotropic_vs_aligned"] = {
    "E3_isotropic_mean_A": float(np.nanmean(ens["E3_zeta_designed"][2]["A"])),
    "E4_aligned_max_mean_A": float(np.nanmean(ens["E4_aligned_max"][2]["A"])),
    "E4_aligned_mid_mean_A": float(np.nanmean(ens["E4_aligned_mid"][2]["A"])),
    "E1_gaussian_mean_A": float(np.nanmean(ens["E1_gaussian_solenoidal"][2]["A"])),
    "E1_gaussian_mean_P": float(np.nanmean(ens["E1_gaussian_solenoidal"][2]["P"])),
    "E1_gaussian_mean_p_norm": float(np.nanmean(ens["E1_gaussian_solenoidal"][2]["p_norm"])),
}
cv = controls["9.5_isotropic_vs_aligned"]
print(f"  9.5 mean A: isotropic {cv['E3_isotropic_mean_A']:+.5f}, "
      f"aligned-to-max {cv['E4_aligned_max_mean_A']:+.5f}, "
      f"aligned-to-mid {cv['E4_aligned_mid_mean_A']:+.5f}, "
      f"Gaussian field {cv['E1_gaussian_mean_A']:+.5f}")
require("9.5 mean production is zero for isotropic orientation and positive only "
        "under preferential alignment",
        abs(cv["E3_isotropic_mean_A"]) < 5e-3 and cv["E4_aligned_max_mean_A"] > 0.3)
require("9.5 CONTROL: a Gaussian solenoidal field has zero mean production",
        abs(cv["E1_gaussian_mean_p_norm"]) < 0.01,
        f"mean p_norm = {cv['E1_gaussian_mean_p_norm']:+.2e} -- Gaussian fields have "
        "no alignment preference, so <P> > 0 in turbulence is dynamical, not kinematic")

# 9.6 does any Thales correlation survive conditioning on zeta?
S, w, c = ens["E1_gaussian_solenoidal"]
m = np.isfinite(c["A"])
z, h, p, A = c["zeta"][m], c["h"][m], c["p_norm"][m], c["A"][m]
absp = np.abs(p)
raw_h = float(np.corrcoef(h, absp)[0, 1])
part_scan = []
for nz in (10, 40, 160, 640, 2560):
    part_scan.append({"n_zeta_bins": nz,
                      "partial_corr_h": st.binned_partial(h, absp, z, nz)["partial_corr"],
                      "partial_corr_absA": st.binned_partial(np.abs(A), absp, z,
                                                             nz)["partial_corr"]})
    print(f"  9.6 {nz:4d} zeta bins: partial corr(h, |p|) = "
          f"{part_scan[-1]['partial_corr_h']:+.2e}, "
          f"partial corr(|A|, |p|) = {part_scan[-1]['partial_corr_absA']:+.4f}")
controls["9.6_partial_correlations"] = {
    "corr(h, |p|)_unconditional": raw_h, "scan": part_scan}
print(f"  9.6 unconditional corr(h, |p_norm|) = {raw_h:+.4f}")
require("9.6 the apparent h correlation is mediated entirely by zeta and decays to "
        "zero as the conditioning is refined",
        abs(part_scan[-1]["partial_corr_h"]) < 0.5 * abs(part_scan[0]["partial_corr_h"])
        and abs(part_scan[-1]["partial_corr_absA"]) > 0.3,
        f"partial corr(h,|p|): {part_scan[0]['partial_corr_h']:+.4f} -> "
        f"{part_scan[-1]['partial_corr_h']:+.2e} while |A| holds at "
        f"{part_scan[-1]['partial_corr_absA']:+.4f}")

# ------------------------------------------------------------------ persist
report = {
    "n_samples": {k: int(len(v[2]["zeta"])) for k, v in ens.items()},
    "bounds": bounds,
    "phase5_conditional_spread": cond,
    "phase5_sign_prediction": sign_tests,
    "phase6_thales_identities": thales,
    "phase6_determinism": det,
    "phase6_conditional_mutual_information": cmi,
    "phase6_leakage_resolution_scan": scan,
    "phase6_residual_altitude_trap": {
        "linear_control": trap_linear, "flexible_control": trap_flexible,
        "control_flexibility_scan": trap_scan,
        "interpretation": "the Appendix-B residual test is a nonlinearity detector, "
                          "not an independent-information test"},
    "phase9_controls": controls,
    "failures": FAILURES,
}
(OUT / "ensembles.json").write_text(json.dumps(report, indent=2, default=float))

# subsample for figures
def sub(c, n=60_000, seed=0):
    m = np.isfinite(c["A"]) & np.isfinite(c["p_norm"])
    idx = np.flatnonzero(m)
    r = np.random.default_rng(seed)
    if len(idx) > n:
        idx = r.choice(idx, n, replace=False)
    return {k: np.asarray(c[k])[idx] for k in ("zeta", "A", "P", "p_norm", "h", "L", "D",
                                               "Qtot", "nS2", "w2")}


np.savez_compressed(OUT / "ensemble_samples.npz",
                    **{f"{tag}__{k}": v for tag, (S, w, c) in ens.items()
                       for k, v in sub(c).items()})
print(f"\n{'ALL CHECKS PASSED' if not FAILURES else 'FAILURES: ' + str(FAILURES)}")
print(f"written: {OUT / 'ensembles.json'}, {OUT / 'ensemble_samples.npz'}")
sys.exit(0 if not FAILURES else 1)
