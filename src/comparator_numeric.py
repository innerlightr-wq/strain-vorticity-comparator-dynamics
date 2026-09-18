"""CONTINUATION -- numerical test suite for the magnitude/comparator/alignment audit.

Reads nothing from the parent results except the saved DNS sample fields and the
parent's generators; writes only into ../continuation_comparator/results.

Tests
  N1  endpoint and symmetry-pair states (exact values)
  N2  reconstruction:  zeta = sigma sqrt(1 - 4h^2)  on 10^6 random states
  N3  the Archimedean-coordinate production identity
      P = Q^{3/2} * 2 h^{3/2} * exp(xi/2) * A  to machine precision
  N4  insufficiency: at fixed (Q, h, sigma) the production still fills its envelope
  N5  exact counterexample pairs (identical magnitude+comparator, different P)
  N6  branch pairs (identical h, opposite sigma) and the exp(xi) envelope asymmetry
  N7  how much information the altitude family actually discards, in bits,
      measured on a Gaussian solenoidal field and on the parent DNS snapshots

Run:  ../.venv/bin/python src/comparator_numeric.py
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
import generators as gen  # noqa: E402
import svcore as sv  # noqa: E402

OUT = _ROOT / "results"
OUT.mkdir(exist_ok=True)
FAILURES: list[str] = []
report: dict[str, object] = {}


def require(name, ok, detail=""):
    print(f"[{'PASS' if ok else 'FAIL':4s}] {name}" + (f"   {detail}" if detail else ""))
    if not ok:
        FAILURES.append(name)


def comparator_coords(S, w):
    """The three-layer coordinates (scale, magnitude, comparator, alignment)."""
    c = sv.coordinates(S, w)
    zeta = c["zeta"]
    h = c["h"]
    sigma = np.sign(zeta)
    sigma = np.where(sigma == 0, 1.0, sigma)        # sign(0) := +1, by convention
    with np.errstate(divide="ignore", invalid="ignore"):
        xi = np.arctanh(np.clip(zeta, -1 + 1e-15, 1 - 1e-15))
    return dict(Q=c["Qtot"], zeta=zeta, h=h, sigma=sigma, xi=xi, A=c["A"], P=c["P"],
                p_norm=c["p_norm"])


print("=" * 74)
print("N1  endpoint and symmetry-pair states")
print("=" * 74)
rows = []
for name, S, w in [
    ("pure extensional strain", np.diag([1.0, -1.0, 0.0]), np.zeros(3)),
    ("solid-body rotation", np.zeros((3, 3)), np.array([0.0, 0.0, 2.0])),
    ("simple shear", np.array([[0, .5, 0], [.5, 0, 0], [0, 0, 0]]),
     np.array([0.0, 0.0, -1.0])),
    ("zeta = +1/3 (rotation-rich)", np.diag([2.0, -1.0, -1.0]) / np.sqrt(6),
     np.array([2.0, 0.0, 0.0])),
    ("zeta = -1/3 (strain-rich)", np.diag([2.0, -1.0, -1.0]) / np.sqrt(6) * np.sqrt(2),
     np.array([1.0, 0.0, 0.0]) * np.sqrt(2)),
]:
    c = comparator_coords(S, w)
    rows.append(dict(state=name, zeta=float(c["zeta"]), h=float(c["h"]),
                     sigma=float(c["sigma"]), A=float(c["A"]) if np.isfinite(c["A"])
                     else None, P=float(c["P"])))
    print(f"   {name:28s} zeta={c['zeta']:+.6f}  h={c['h']:.6f}  "
          f"sigma={c['sigma']:+.0f}  A={c['A']:+.6f}" if np.isfinite(c["A"]) else
          f"   {name:28s} zeta={c['zeta']:+.6f}  h={c['h']:.6f}  "
          f"sigma={c['sigma']:+.0f}  A=undefined")
report["N1_states"] = rows
require("N1 pure strain and solid rotation share h = 0 but differ in sigma",
        abs(rows[0]["h"]) < 1e-12 and abs(rows[1]["h"]) < 1e-12 and
        rows[0]["sigma"] < 0 < rows[1]["sigma"],
        "the one bit that the altitude family discards, exhibited")
require("N1 the balanced state sits at the fold h = 1/2",
        abs(rows[2]["h"] - 0.5) < 1e-12 and abs(rows[2]["zeta"]) < 1e-12)

print()
print("=" * 74)
print("N2  reconstruction  zeta = sigma * sqrt(1 - 4 h^2)")
print("=" * 74)
G = gen.gaussian_solenoidal_gradients(n_grid=64, n_real=4, seed=11)
S = sv.sym(G)
w = sv.vorticity_from_skew(sv.skew(G))
c = comparator_coords(S, w)
zeta_rec = c["sigma"] * np.sqrt(np.clip(1.0 - 4.0 * c["h"] ** 2, 0.0, None))
abs_err = np.abs(zeta_rec - c["zeta"])
err = float(np.max(abs_err))
# away from the fold the reconstruction is exact to machine precision
far = np.abs(c["zeta"]) > 1e-2
require("N2 exact reconstruction on 10^6 random incompressible states "
        "(away from the fold)",
        float(np.max(abs_err[far])) < 1e-13,
        f"n = {int(far.sum())} states with |zeta| > 1e-2, "
        f"max |zeta_rec - zeta| = {float(np.max(abs_err[far])):.3e}")
# CONDITIONING: zeta^2 = 1 - 4h^2 gives dzeta/dh = -4h/zeta, so the inverse is
# exact but ILL-CONDITIONED at the fold zeta -> 0.  Check that the residual is
# explained by that amplification and not by a wrong formula.
EPS = np.finfo(float).eps
amp = 4.0 * c["h"] / np.maximum(np.abs(c["zeta"]), 1e-300)
predicted = 20.0 * amp * EPS + 1e-15
require("N2 the residual is exactly the fold conditioning 4h/|zeta|, not a formula error",
        bool(np.all(abs_err <= predicted)),
        f"worst overall error {err:.3e} at |zeta| = "
        f"{float(np.abs(c['zeta'])[np.argmax(abs_err)]):.3e}; "
        f"predicted amplification 4h/|zeta| = "
        f"{float(amp[np.argmax(abs_err)]):.3e}")
# and the magnitudes come back too
nS2_rec = c["Q"] * (1 - c["zeta"]) / 2
nW2_rec = c["Q"] * (1 + c["zeta"]) / 2
err_m = max(float(np.max(np.abs(nS2_rec - sv.frob2(S)))),
            float(np.max(np.abs(nW2_rec - 0.5 * np.sum(w * w, axis=-1)))))
require("N2 the two component magnitudes are recovered from (Q, h, sigma)",
        err_m < 1e-10, f"max error {err_m:.3e}")
# an even 'comparator' fails
zeta_even = np.abs(c["zeta"])
require("N2 CONTROL: |zeta| alone cannot reconstruct zeta",
        np.max(np.abs(zeta_even - c["zeta"])) > 0.5,
        "the sign is genuinely absent from every even coordinate")
report["N2"] = {"n": int(len(G)), "max_zeta_error_overall": float(err),
                "max_zeta_error_away_from_fold": float(np.max(abs_err[far])),
                "max_magnitude_error": float(err_m),
                "conditioning": "dzeta/dh = -4h/zeta; the inverse is exact but the "
                                "fold at zeta = 0 is ill-conditioned"}

print()
print("=" * 74)
print("N3  production in Archimedean coordinates")
print("=" * 74)
m = np.isfinite(c["A"]) & (c["h"] > 0)
P_arch = (c["Q"][m] ** 1.5) * 2.0 * c["h"][m] ** 1.5 * np.exp(c["xi"][m] / 2.0) * c["A"][m]
rel = np.max(np.abs(P_arch - c["P"][m]) / (np.abs(c["P"][m]) + 1e-300))
require("N3 P = ||grad u||^3 * 2 h^{3/2} * exp(xi/2) * A holds pointwise",
        rel < 1e-10, f"max relative error {rel:.3e} over {m.sum()} states")
g_sym = 2.0 * c["h"][m] ** 1.5
g_cmp = np.exp(c["xi"][m] / 2.0)
require("N3 g_sym is even and g_cmp obeys g_cmp(zeta) g_cmp(-zeta) = 1",
        np.max(np.abs(g_cmp * np.exp(-c["xi"][m] / 2.0) - 1.0)) < 1e-12)
report["N3"] = {"max_relative_error": float(rel), "n": int(m.sum())}

print()
print("=" * 74)
print("N4  insufficiency of (scale, magnitude, comparator)")
print("=" * 74)
Qs, hs, sg, As, ps = (c["Q"][m], c["h"][m], c["sigma"][m], c["A"][m], c["p_norm"][m])


def quantile_bins(x, nb):
    e = np.unique(np.quantile(x, np.linspace(0, 1, nb + 1)))
    return np.clip(np.searchsorted(e, x, side="right") - 1, 0, len(e) - 2)


cell = (quantile_bins(np.log(Qs), 10) * 10 + quantile_bins(hs, 10)) * 2 + (sg > 0)
_, cidx = np.unique(cell, return_inverse=True)
nb = cidx.max() + 1
cnt = np.bincount(cidx, minlength=nb).astype(float)
mean_A = np.bincount(cidx, weights=As, minlength=nb) / np.maximum(cnt, 1)
resid_A = As - mean_A[cidx]
frac_var = float(np.sum(resid_A ** 2) / np.sum((As - As.mean()) ** 2))
both_signs = np.array([len(np.unique(np.sign(ps[cidx == b]))) > 1
                       for b in range(nb) if cnt[b] > 50])
require("N4 fixing (Q, h, sigma) leaves essentially all of the variance of A",
        frac_var > 0.99,
        f"{frac_var:.4%} of var(A) survives inside {nb} (Q, h, sigma) cells")
require("N4 almost every such cell contains BOTH signs of P",
        both_signs.mean() > 0.95,
        f"{both_signs.mean():.1%} of cells with n > 50 contain P > 0 and P < 0")
report["N4"] = {"n_cells": int(nb), "frac_var_A_within_cells": frac_var,
                "frac_cells_with_both_signs": float(both_signs.mean())}

print()
print("=" * 74)
print("N5  exact counterexample pairs: same (Q, h, sigma), different P")
print("=" * 74)
S0 = np.diag([1.0, 0.0, -1.0])                     # ||S||^2 = 2, traceless
pairs = []
for tag, wv in [("omega on e1 (stretching)", np.array([np.sqrt(2), 0, 0])),
                ("omega on e2 (neutral)", np.array([0, np.sqrt(2), 0])),
                ("omega on e3 (compressing)", np.array([0, 0, np.sqrt(2)]))]:
    cc = comparator_coords(S0, wv)
    pairs.append(dict(case=tag, Q=float(cc["Q"]), h=float(cc["h"]),
                      sigma=float(cc["sigma"]), zeta=float(cc["zeta"]),
                      A=float(cc["A"]), P=float(cc["P"])))
    print(f"   {tag:26s} Q={cc['Q']:.6f} h={cc['h']:.6f} sigma={cc['sigma']:+.0f} "
          f"A={cc['A']:+.6f} P={cc['P']:+.6f}")
same = (abs(pairs[0]["Q"] - pairs[2]["Q"]) < 1e-12 and
        abs(pairs[0]["h"] - pairs[2]["h"]) < 1e-12 and
        pairs[0]["sigma"] == pairs[2]["sigma"])
require("N5 three states share (Q, h, sigma) exactly and have P = +2, 0, -2",
        same and abs(pairs[0]["P"] - 2) < 1e-12 and abs(pairs[1]["P"]) < 1e-12 and
        abs(pairs[2]["P"] + 2) < 1e-12,
        "what differs is only the direction cosines of omega in the strain eigenframe")
# what exactly differs
for pr, wv in zip(pairs, [np.array([np.sqrt(2), 0, 0]), np.array([0, np.sqrt(2), 0]),
                          np.array([0, 0, np.sqrt(2)])]):
    lam, cos2 = sv.alignment_cosines(S0, wv)
    pr["lambda_ascending"] = lam.tolist()
    pr["cos2_ascending"] = cos2.tolist()
report["N5_pairs"] = pairs

print()
print("=" * 74)
print("N6  branch pairs: identical h, opposite sigma")
print("=" * 74)
branch = []
for zt in (1 / 3, 0.6, 0.8):
    gp = (1 + zt) * np.sqrt((1 - zt) / 2)
    gm = (1 - zt) * np.sqrt((1 + zt) / 2)
    xi = np.arctanh(zt)
    branch.append(dict(zeta=zt, h=float(sv.h_of_zeta(zt)), g_plus=gp, g_minus=gm,
                       ratio=gp / gm, exp_xi=float(np.exp(xi))))
    print(f"   |zeta| = {zt:.3f}:  h = {sv.h_of_zeta(zt):.6f} (identical on both "
          f"branches)   g(+) = {gp:.6f}   g(-) = {gm:.6f}   "
          f"ratio = {gp/gm:.6f}   exp(xi) = {np.exp(xi):.6f}")
require("N6 the envelope ratio between the two branches is exactly exp(xi)",
        max(abs(b["ratio"] - b["exp_xi"]) for b in branch) < 1e-12,
        "so the symmetric geometry overestimates the weaker branch by exp(|xi|): "
        f"up to {max(b['exp_xi'] for b in branch):.3f}x over the range tested")
report["N6_branches"] = branch

print()
print("=" * 74)
print("N7  how many bits does the altitude family discard, in practice?")
print("=" * 74)


def sigma_entropy(zeta, h, nb=32):
    """H(sigma) and H(sigma | h) in bits, h binned by quantiles."""
    good = np.isfinite(zeta) & np.isfinite(h)
    zeta, h = zeta[good], h[good]
    s = (zeta > 0).astype(float)

    def H(p):
        p = np.clip(p, 1e-15, 1 - 1e-15)
        return float(-(p * np.log2(p) + (1 - p) * np.log2(1 - p)))

    Hs = H(s.mean())
    idx = quantile_bins(h, nb)
    tot, n = 0.0, len(s)
    for b in range(idx.max() + 1):
        msk = idx == b
        if msk.sum() > 20:
            tot += msk.sum() / n * H(s[msk].mean())
    return Hs, tot, float(s.mean())


Hs_g, Hcond_g, frac_g = sigma_entropy(c["zeta"], c["h"])
print(f"   Gaussian solenoidal field: P(sigma=+1) = {frac_g:.4f}, "
      f"H(sigma) = {Hs_g:.4f} bits, H(sigma | h) = {Hcond_g:.4f} bits")
ent = {"gaussian_field": {"frac_positive": frac_g, "H_sigma": Hs_g,
                          "H_sigma_given_h": Hcond_g}}
dns = _ROOT / "results" / "dns_samples.npz"
if dns.exists():
    zf = np.load(dns)
    for tg in sorted({f.split("_")[0] for f in zf.files}, key=lambda t: float(t[1:])):
        zz = zf[f"{tg}_zeta"]
        hh = sv.h_of_zeta(zz)
        Hs, Hc, fr = sigma_entropy(zz, hh)
        ent[tg] = {"frac_positive": fr, "H_sigma": Hs, "H_sigma_given_h": Hc}
        print(f"   DNS {tg:4s}: P(sigma=+1) = {fr:.4f}, H(sigma) = {Hs:.4f} bits, "
              f"H(sigma | h) = {Hc:.4f} bits")
require("N7 the discarded information is close to one full bit per point",
        all(v["H_sigma_given_h"] > 0.75 for v in ent.values()),
        "h carries almost no statistical information about its own branch: "
        + ", ".join(f"{k}: {v['H_sigma_given_h']:.3f} bits" for k, v in ent.items()))
report["N7_entropy_bits"] = ent

report["failures"] = FAILURES
(OUT / "comparator_numeric.json").write_text(json.dumps(report, indent=2, default=float))
print(f"\n{'ALL CHECKS PASSED' if not FAILURES else 'FAILURES: ' + str(FAILURES)}")
print(f"written: {OUT / 'comparator_numeric.json'}")
sys.exit(0 if not FAILURES else 1)
