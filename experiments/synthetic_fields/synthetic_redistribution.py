"""WEIGHTED COMPARATOR AUDIT -- section N: the redistribution test.

Three genuine periodic incompressible velocity fields with DIFFERENT spatial
distributions of the comparator xi, all rescaled to the same total enstrophy
int E_W dx (which, by the kinematic identity int E_S = int E_W, also matches the
total strain energy exactly).

Question: can any exactly pressure-free functional distinguish them through the
strain-vorticity balance?

The pressure-free class is  Phi = c E_S + psi(E_W)  (proved in weighted_exact.py),
so its ENTIRE dependence on the strain sector is the single number c * int E_S,
which is pinned to c * int E_W.  This script shows that number is identical across
the three fields while the comparator functionals differ by large factors.

Run:  ../../../../.venv/bin/python src/synthetic_redistribution.py
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
sys.path.insert(0, str(_ROOT / "experiments" / "dns_checks"))
import dns_taylor_green as dns  # noqa: E402
import svcore as sv  # noqa: E402

OUT = _ROOT / "results"
OUT.mkdir(exist_ok=True)
FAILURES: list[str] = []


def require(name, ok, detail=""):
    print(f"[{'PASS' if ok else 'FAIL':4s}] {name}" + (f"   {detail}" if detail else ""))
    if not ok:
        FAILURES.append(name)


N = 64
sp = dns.Spectral(N, nu=0.0)
X, Y, Z = np.meshgrid(sp.x, sp.x, sp.x, indexing="ij")
mean = lambda a: float(np.mean(a))


def project(u):
    return sp.project(sp.vfft(u)) * sp.mask


def band_limited(seed, kmin, kmax):
    """A random solenoidal field with energy in a wavenumber shell."""
    rng = np.random.default_rng(seed)
    kk = np.sqrt(sp.K2)
    amp = ((kk >= kmin) & (kk <= kmax)).astype(float)
    uh = np.stack([sp.fft(rng.normal(size=(N,) * 3)) * amp for _ in range(3)])
    return sp.project(uh) * sp.mask


def abc_flow(A=1.0, B=1.0, C=1.0, k=2):
    """ABC flow: strongly rotation-dominated (Beltrami, omega = k u)."""
    u = np.stack([A * np.sin(k * Z) + C * np.cos(k * Y),
                  B * np.sin(k * X) + A * np.cos(k * Z),
                  C * np.sin(k * Y) + B * np.cos(k * X)])
    return project(u)


def strain_cellular(k=2):
    """A cellular field dominated by plane strain over most of the volume."""
    u = np.stack([np.sin(k * X) * np.cos(k * Y) * np.cos(k * Z),
                  -np.cos(k * X) * np.sin(k * Y) * np.cos(k * Z),
                  np.zeros_like(X)])
    return project(u)


def bundle(uh):
    G = sp.velocity_gradient(uh)
    S, W = sv.sym(G), sv.skew(G)
    E_S = np.sum(S * S, axis=(-2, -1))
    E_W = np.sum(W * W, axis=(-2, -1))
    f = E_W - E_S
    fh = sp.fft(f)
    K2 = np.where(sp.K2 == 0, 1.0, sp.K2)
    ph = -fh / K2
    H = np.empty_like(G)
    for i in range(3):
        for j in range(3):
            H[..., i, j] = sp.ifft(-sp.K[i] * sp.K[j] * ph)
    trH = np.trace(H, axis1=-2, axis2=-1)
    H_dev = H - np.eye(3) * (trH / 3.0)[..., None, None]
    return dict(E_S=E_S, E_W=E_W, xi=0.5 * np.log(E_W / E_S),
                SH=np.sum(S * H_dev, axis=(-2, -1)))


# three structurally different fields
u_abc = sp.vifft(abc_flow())          # Beltrami: maximal helicity, symmetric
u_str = sp.vifft(strain_cellular())        # Taylor-Green cellular, k = 2
u_str1 = sp.vifft(strain_cellular(k=1))    # same family at k = 1
raw = {
    "random broadband (k 1-8)": band_limited(1, 1.0, 8.0),
    "random small-scale (k 6-12)": band_limited(7, 6.0, 12.0),
    "ABC (Beltrami, symmetric)": abc_flow(),
    "Taylor-Green cellular (symmetric)": strain_cellular(),
    # mixing DIFFERENT wavenumbers breaks the symmetry that makes the
    # xi-weighted pressure integral cancel in each pure field
    "ABC(k=2) + 2 x cellular(k=1) (asymmetric)": project(u_abc + 2.0 * u_str1),
}

# rescale each to a common total enstrophy
TARGET = 1.0
fields = {}
for tag, uh in raw.items():
    b = bundle(uh)
    fields[tag] = bundle(uh * np.sqrt(TARGET / mean(b["E_W"])))

print("=" * 78)
print("N  REDISTRIBUTION TEST: matched sector totals, different xi distributions")
print("=" * 78)
rows = {}
for tag, F in fields.items():
    xi = F["xi"]
    xSH = xi * F["SH"]
    rows[tag] = {
        "int_E_S": mean(F["E_S"]), "int_E_W": mean(F["E_W"]),
        "xi_mean": mean(xi), "xi_std": float(np.std(xi)),
        "xi_q10": float(np.quantile(xi, 0.10)), "xi_q90": float(np.quantile(xi, 0.90)),
        "frac_rotation_dominated": float(np.mean(xi > 0)),
        "Xi_1": mean(F["E_S"] * xi),
        "xi_bar_E_S": mean(F["E_S"] * xi) / mean(F["E_S"]),
        "int_E_W_sq": mean(F["E_W"] ** 2),
        "pressure_term_alpha1": -2.0 * mean(xSH),
        "pressure_rel": (abs(mean(xSH)) / mean(np.abs(xSH))
                         if mean(np.abs(xSH)) > 0 else 0.0),
        "mean_SH_rel": (abs(mean(F["SH"])) / mean(np.abs(F["SH"]))
                        if mean(np.abs(F["SH"])) > 0 else 0.0),
    }
    r = rows[tag]
    print(f"\n   {tag}")
    print(f"      int E_S = {r['int_E_S']:.6f}   int E_W = {r['int_E_W']:.6f}")
    print(f"      xi:  mean {r['xi_mean']:+.4f}  std {r['xi_std']:.4f}  "
          f"[q10 {r['xi_q10']:+.3f}, q90 {r['xi_q90']:+.3f}]  "
          f"frac(xi>0) {r['frac_rotation_dominated']:.3f}")
    print(f"      Xi_1 = <E_S xi> = {r['Xi_1']:+.6f}    <xi>_E_S = "
          f"{r['xi_bar_E_S']:+.6f}    <E_W^2> = {r['int_E_W_sq']:.4f}")
    print(f"      pressure part of dXi_1/dt = {r['pressure_term_alpha1']:+.6e}  "
          f"(relative {r['pressure_rel']:.2e});  int S:H_dev relative "
          f"{r['mean_SH_rel']:.1e}")

vals_ES = [r["int_E_S"] for r in rows.values()]
stds = [r["xi_std"] for r in rows.values()]
xibars = [r["xi_bar_E_S"] for r in rows.values()]
require("N the five fields have identical sector totals int E_S = int E_W = 1",
        (max(vals_ES) - min(vals_ES)) / np.mean(vals_ES) < 1e-10,
        f"int E_S spread {(max(vals_ES) - min(vals_ES)):.2e}")
require("N the xi distributions are genuinely different",
        max(stds) / min(stds) > 3.0,
        f"std(xi) ranges over [{min(stds):.3f}, {max(stds):.3f}], a factor "
        f"{max(stds)/min(stds):.1f}")
require("N the comparator functionals differ correspondingly",
        max(xibars) - min(xibars) > 0.2,
        f"<xi>_E_S ranges over [{min(xibars):+.4f}, {max(xibars):+.4f}]")
require("N the unweighted identity holds in every field",
        all(r["mean_SH_rel"] < 1e-12 for r in rows.values()),
        f"worst relative {max(r['mean_SH_rel'] for r in rows.values()):.1e}")
require("N BUT the entire strain-sector content of the pressure-free class "
        "(c * int E_S) is IDENTICAL across all five",
        (max(vals_ES) - min(vals_ES)) / np.mean(vals_ES) < 1e-10,
        "the only E_S-dependence exact pressure freedom permits is c*int E_S, "
        "kinematically pinned to c*int E_W: it cannot see the redistribution at all")
generic = ["random broadband (k 1-8)", "random small-scale (k 6-12)",
           "ABC(k=2) + 2 x cellular(k=1) (asymmetric)"]
symmetric = ["ABC (Beltrami, symmetric)", "Taylor-Green cellular (symmetric)"]
require("N the alpha = 1 pressure term is nonzero for generic fields",
        all(rows[t]["pressure_rel"] > 1e-4 for t in generic),
        "relative sizes " + ", ".join(f"{t.split()[0]}: {rows[t]['pressure_rel']:.1e}"
                                      for t in generic))
require("N and can vanish by SYMMETRY for special fields (an accident, not a law)",
        all(rows[t]["pressure_rel"] < 1e-4 for t in symmetric),
        "ABC and Taylor-Green are symmetric enough that int xi S:H_dev cancels; "
        "this is a property of those fields, not of the functional")

out = {"target_enstrophy": TARGET, "fields": rows, "failures": FAILURES}
(OUT / "synthetic_redistribution.json").write_text(json.dumps(out, indent=2,
                                                             default=float))
print(f"\n{'ALL CHECKS PASSED' if not FAILURES else 'FAILURES: ' + str(FAILURES)}")
print(f"written: {OUT / 'synthetic_redistribution.json'}")
sys.exit(0 if not FAILURES else 1)
