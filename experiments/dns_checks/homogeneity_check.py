"""PHASE 8 -- numerical confirmation that the MEAN allocation of a homogeneous
incompressible flow sits exactly at the Thales apex.

Exactly, ||S||_F^2 - ||W||_F^2 = tr((grad u)^2), and for a homogeneous
incompressible field the mean of the right-hand side vanishes (it is a divergence
plus a term proportional to div u).  Hence

    <||S||_F^2> = <||W||_F^2>        =>        <a> = <b> = 1/2 ,

for EVERY homogeneous incompressible flow, whatever its dynamics.  Apex occupancy
is therefore forced by homogeneity and carries no dynamical information -- in the
vocabulary of the Archimedean note itself, occupancy without meaning.

Run:  .venv/bin/python src/homogeneity_check.py
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
report = {}
FAILURES = []


def require(name, ok, detail=""):
    print(f"[{'PASS' if ok else 'FAIL':4s}] {name}" + (f"   {detail}" if detail else ""))
    if not ok:
        FAILURES.append(name)


print("=" * 72)
print("PHASE 8: is apex occupancy dynamical or kinematic?")
print("=" * 72)
G = gen.gaussian_solenoidal_gradients(n_grid=64, n_real=2, seed=11)
S, W = sv.sym(G), sv.skew(G)
nS2, nW2 = sv.frob2(S), sv.frob2(W)
trA2 = np.trace(G @ G, axis1=-2, axis2=-1)
ratio = float(nW2.mean() / nS2.mean())
rel = float(trA2.mean() / (nS2.mean() + nW2.mean()))
print(f"  Gaussian solenoidal field, n = {len(G)}:")
print(f"    <||S||^2> = {nS2.mean():.6g}   <||W||^2> = {nW2.mean():.6g}   "
      f"ratio = {ratio:.6f}")
print(f"    <tr(grad u ^2)> / <||grad u||^2> = {rel:+.3e}")
print(f"    <a> = {float((nS2 / (nS2 + nW2)).mean()):.6f} (pointwise mean of the ratio) "
      f"vs {float(nS2.mean() / (nS2.mean() + nW2.mean())):.6f} (ratio of means)")
require("Gaussian field: <||W||^2> = <||S||^2> to sampling accuracy",
        abs(ratio - 1.0) < 0.02, f"ratio {ratio:.6f}")
report["gaussian_field"] = {"n": int(len(G)), "mean_nS2": float(nS2.mean()),
                            "mean_nW2": float(nW2.mean()), "ratio": ratio,
                            "mean_trA2_normalized": rel,
                            "mean_of_pointwise_a": float((nS2 / (nS2 + nW2)).mean()),
                            "a_from_ratio_of_means":
                                float(nS2.mean() / (nS2.mean() + nW2.mean()))}

dns = OUT / "dns_samples.npz"
if dns.exists():
    z = np.load(dns)
    tags = sorted({f.split("_")[0] for f in z.files})
    report["dns"] = {}
    print("\n  Navier-Stokes DNS snapshots (subsampled points):")
    for tg in tags:
        Q = z[f"{tg}_Q"]
        ew = z[f"{tg}_e_w"]           # = ||W||_F^2 exactly
        nS2d = Q - ew
        r = float(ew.mean() / nS2d.mean())
        print(f"    {tg}: <||W||^2>/<||S||^2> = {r:.4f}, "
              f"ratio-of-means a = {float(nS2d.mean() / Q.mean()):.4f}, "
              f"mean pointwise zeta = {float(z[f'{tg}_zeta'].mean()):+.4f}")
        report["dns"][tg] = {"ratio_nW2_over_nS2": r,
                             "a_ratio_of_means": float(nS2d.mean() / Q.mean()),
                             "mean_pointwise_zeta": float(z[f"{tg}_zeta"].mean())}
    rs = [v["ratio_nW2_over_nS2"] for v in report["dns"].values()]
    require("DNS: the ratio of means sits at 1 (apex) at every snapshot",
            all(abs(r - 1.0) < 0.15 for r in rs),
            f"ratios {[round(r, 4) for r in rs]} (subsampled, so a few percent scatter)")
    require("the POINTWISE mean of zeta is NOT zero, so the apex is a property of "
            "the means only",
            any(abs(v["mean_pointwise_zeta"]) > 0.02 for v in report["dns"].values()),
            "pointwise zeta means: " +
            ", ".join(f"{v['mean_pointwise_zeta']:+.4f}" for v in report["dns"].values()))

(OUT / "homogeneity_check.json").write_text(json.dumps(report, indent=2))
print(f"\n{'ALL CHECKS PASSED' if not FAILURES else 'FAILURES: ' + str(FAILURES)}")
print(f"written: {OUT / 'homogeneity_check.json'}")
sys.exit(0 if not FAILURES else 1)
