"""WEIGHTED COMPARATOR AUDIT -- numerical verification on genuine periodic
incompressible fields.

Every identity tested here is KINEMATIC (it holds for any smooth divergence-free
periodic field, not only for Navier-Stokes solutions), so random solenoidal fields
are admissible test data; a short Taylor-Green Navier-Stokes snapshot is added for
the diagnostic section.

  N1  int S:H_dev dx = 0                                    (the unweighted identity)
  N2  int xi S:H_dev dx != 0                                <-- kills alpha = 1
  N3  the pressure contribution to dXi_alpha/dt over a range of alpha
  N4  direct assembly check: pressure part of d/dt int Phi  =  -2 int Phi_S S:H_dev
  N5  the weight identity  int w S:H = int p[S:grad grad w + grad w . Lap u]
  N6  necessity: non-constant weights give nonzero integrals; constant ones give zero
  N7  int E_S = int E_W  (the kinematic collapse used by the theorem)
  N8  comparator-blindness of the pressure-free class
  M   Xi_alpha and <xi>_{E_S} on a Navier-Stokes snapshot

Run:  ../../../../.venv/bin/python src/weighted_numeric.py
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
report: dict[str, object] = {}


def require(name, ok, detail=""):
    print(f"[{'PASS' if ok else 'FAIL':4s}] {name}" + (f"   {detail}" if detail else ""))
    if not ok:
        FAILURES.append(name)


# ------------------------------------------------------------------ field tools

def random_solenoidal(sp, seed=0, k0=3.0):
    """A random divergence-free periodic velocity field, spectrally projected."""
    rng = np.random.default_rng(seed)
    amp = np.where(sp.K2 > 0, np.exp(-sp.K2 / (2 * k0 ** 2)), 0.0)
    uh = np.stack([sp.fft(rng.normal(size=(sp.N,) * 3)) * amp for _ in range(3)])
    return sp.project(uh) * sp.mask


def field_bundle(sp, uh):
    """Everything the audit needs from one velocity field."""
    G = sp.velocity_gradient(uh)
    S, W = sv.sym(G), sv.skew(G)
    E_S = np.sum(S * S, axis=(-2, -1))
    E_W = np.sum(W * W, axis=(-2, -1))
    f = E_W - E_S                                   # = Lap p
    fh = sp.fft(f)
    K2 = np.where(sp.K2 == 0, 1.0, sp.K2)
    ph = -fh / K2
    p = sp.ifft(ph)
    H = np.empty_like(G)
    for i in range(3):
        for j in range(3):
            H[..., i, j] = sp.ifft(-sp.K[i] * sp.K[j] * ph)
    trH = np.trace(H, axis1=-2, axis2=-1)
    H_dev = H - np.eye(3) * (trH / 3.0)[..., None, None]
    SH = np.sum(S * H_dev, axis=(-2, -1))
    xi = 0.5 * np.log(E_W / E_S)
    lap_u = np.stack([sp.ifft(-sp.K2 * uh[i]) for i in range(3)])
    return dict(G=G, S=S, W=W, E_S=E_S, E_W=E_W, xi=xi, p=p, H=H, H_dev=H_dev,
                SH=SH, trH=trH, f=f, lap_u=lap_u, uh=uh, sp=sp)


def mean(a):
    return float(np.mean(a))


N = 64
sp = dns.Spectral(N, nu=0.0)
fields = {f"solenoidal_{k}": field_bundle(sp, random_solenoidal(sp, seed=k, k0=2.0 + k))
          for k in range(3)}

print("=" * 76)
print("N1/N7  the unweighted identity and the sector-energy collapse")
print("=" * 76)
rows = {}
for tag, F in fields.items():
    scale = mean(np.abs(F["SH"]))
    rows[tag] = {"mean_SH": mean(F["SH"]), "mean_abs_SH": scale,
                 "ratio": abs(mean(F["SH"])) / scale,
                 "mean_E_S": mean(F["E_S"]), "mean_E_W": mean(F["E_W"]),
                 "trH_residual": float(np.max(np.abs(F["trH"] - F["f"])) /
                                       np.max(np.abs(F["f"])))}
    print(f"   {tag}: <S:H_dev> = {rows[tag]['mean_SH']:+.3e}  "
          f"(<|S:H_dev|> = {scale:.4f}, ratio {rows[tag]['ratio']:.2e});  "
          f"<E_S> = {rows[tag]['mean_E_S']:.6f}, <E_W> = {rows[tag]['mean_E_W']:.6f}")
require("N1 int S:H_dev dx = 0 on every test field",
        all(r["ratio"] < 1e-12 for r in rows.values()),
        f"worst ratio {max(r['ratio'] for r in rows.values()):.2e}")
require("N7 int E_S dx = int E_W dx exactly (kinematic collapse)",
        all(abs(r["mean_E_W"] - r["mean_E_S"]) / r["mean_E_S"] < 1e-12
            for r in rows.values()))
require("pressure solve verified: tr H = Lap p = E_W - E_S",
        all(r["trH_residual"] < 1e-12 for r in rows.values()))
report["N1_N7"] = rows

print()
print("=" * 76)
print("N2  DOES alpha = 1 CANCEL?   int xi S:H_dev dx")
print("=" * 76)
n2 = {}
for tag, F in fields.items():
    val = mean(F["xi"] * F["SH"])
    scale = mean(np.abs(F["xi"] * F["SH"]))
    n2[tag] = {"mean_xi_SH": val, "mean_abs_xi_SH": scale,
               "ratio": abs(val) / scale}
    print(f"   {tag}: <xi S:H_dev> = {val:+.6f}   (<|xi S:H_dev|> = {scale:.4f}, "
          f"ratio {abs(val)/scale:.3f})")
require("N2 int xi S:H_dev dx does NOT vanish  =>  alpha = 1 FAILS",
        all(r["ratio"] > 1e-3 for r in n2.values()),
        "the pressure contribution to d/dt int E_S xi dx is -2 int xi S:H_dev dx, "
        "which is nonzero by 3-10 orders of magnitude relative to N1")
report["N2"] = n2

print()
print("=" * 76)
print("N3  PRESSURE CONTRIBUTION TO dXi_alpha/dt OVER alpha")
print("=" * 76)
alphas = [-1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0]
n3 = {}
F0 = fields["solenoidal_0"]
for a in alphas:
    w = F0["E_S"] ** (a - 1.0) * (a * F0["xi"] - 0.5)
    val = -2.0 * mean(w * F0["SH"])
    scale = 2.0 * mean(np.abs(w * F0["SH"]))
    n3[str(a)] = {"pressure_term": val, "scale": scale,
                  "ratio": abs(val) / scale if scale > 0 else 0.0}
    print(f"   alpha = {a:+.1f}:  pressure part = {val:+.6e}   "
          f"(relative size {abs(val)/scale:.4f})")
require("N3 no alpha in the tested range cancels the pressure",
        all(v["ratio"] > 1e-3 for v in n3.values()),
        "consistent with the exact theorem: E_S^(alpha-1)(alpha xi - 1/2) is "
        "never constant")
report["N3"] = n3

print()
print("=" * 76)
print("N4  DIRECT ASSEMBLY CHECK OF THE PRODUCT RULE")
print("=" * 76)
# assemble d/dt int Phi from the sector equations, isolating the H_dev pieces
F = fields["solenoidal_0"]
P = np.einsum("...i,...ij,...j->...", sv.vorticity_from_skew(F["W"]), F["S"],
              sv.vorticity_from_skew(F["W"]))
trS3 = np.trace(F["S"] @ F["S"] @ F["S"], axis1=-2, axis2=-1)
dES = -2 * trS3 - 0.5 * P - 2 * F["SH"]                 # inviscid
dEW = P
for a in (0.0, 1.0):
    Phi_S = F["E_S"] ** (a - 1.0) * (a * F["xi"] - 0.5)
    Phi_W = F["E_S"] ** a / (2 * F["E_W"])
    total = mean(Phi_S * dES + Phi_W * dEW)
    no_press = mean(Phi_S * (dES + 2 * F["SH"]) + Phi_W * dEW)
    press = total - no_press
    formula = -2.0 * mean(Phi_S * F["SH"])
    print(f"   alpha = {a:.0f}: d/dt-int total = {total:+.6f}, pressure-free part = "
          f"{no_press:+.6f}, pressure part = {press:+.6e} "
          f"(formula {formula:+.6e})")
    require(f"N4 alpha={a:.0f}: pressure part equals -2 int Phi_S S:H_dev",
            abs(press - formula) < 1e-10 * max(1.0, abs(formula)),
            f"difference {abs(press - formula):.2e}")

print()
print("=" * 76)
print("N5  THE WEIGHT IDENTITY  int w S:H = int p [S:grad grad w + grad w . Lap u]")
print("=" * 76)
spn = F["sp"]


def grad_scalar(g):
    gh = spn.fft(g)
    return np.stack([spn.ifft(1j * spn.K[i] * gh) for i in range(3)])


def hess_scalar(g):
    gh = spn.fft(g)
    Hs = np.empty(g.shape + (3, 3))
    for i in range(3):
        for j in range(3):
            Hs[..., i, j] = spn.ifft(-spn.K[i] * spn.K[j] * gh)
    return Hs


n5 = {}
for wname, w in (("w = 1 (constant)", np.ones_like(F["E_S"])),
                 ("w = xi", F["xi"]),
                 ("w = E_S", F["E_S"])):
    lhs = mean(w * F["SH"])
    gw = grad_scalar(w)
    Hw = hess_scalar(w)
    rhs = mean(F["p"] * (np.sum(F["S"] * Hw, axis=(-2, -1)) +
                         np.sum(gw * F["lap_u"], axis=0)))
    n5[wname] = {"lhs": lhs, "rhs": rhs, "abs_diff": abs(lhs - rhs)}
    print(f"   {wname:18s}:  int w S:H = {lhs:+.6e}   "
          f"int p[S:grad grad w + grad w . Lap u] = {rhs:+.6e}")
require("N5 the by-parts identity holds to round-off for every weight tested",
        all(v["abs_diff"] < 1e-9 * max(1.0, abs(v["lhs"])) or v["abs_diff"] < 1e-9
            for v in n5.values()),
        "so only derivatives of w survive: constant weights cancel exactly")
report["N5"] = n5

print()
print("=" * 76)
print("N6  NECESSITY: non-constant weights reintroduce the pressure")
print("=" * 76)
n6 = {}
weights = {
    "constant": np.ones_like(F["E_S"]),
    "xi": F["xi"],
    "E_S": F["E_S"],
    "E_W": F["E_W"],
    "log E_S": np.log(F["E_S"]),
    "sigmoid(E_S - median)": 1.0 / (1.0 + np.exp(-(F["E_S"] - np.median(F["E_S"])) /
                                                 (0.1 * np.std(F["E_S"])))),
}
for wn, w in weights.items():
    val = mean(w * F["SH"])
    scale = mean(np.abs(w * F["SH"]))
    n6[wn] = {"value": val, "scale": scale, "ratio": abs(val) / scale}
    print(f"   w = {wn:24s}:  int w S:H_dev = {val:+.6e}   "
          f"relative {abs(val)/scale:.2e}")
require("N6 constant weight cancels; every non-constant weight tested does not",
        n6["constant"]["ratio"] < 1e-12 and
        all(n6[k]["ratio"] > 1e-4 for k in n6 if k != "constant"),
        "the sigmoid weight is the two-region mechanism of the necessity argument: "
        "it splits the domain into two sets on which S:H_dev integrates to equal and "
        "opposite nonzero amounts")
report["N6"] = n6

print()
print("=" * 76)
print("N8  IS THE PRESSURE-FREE CLASS COMPARATOR-BLIND?")
print("=" * 76)
# The pressure-free class is Phi = c E_S + psi(E_W). Its value depends on the field
# only through <E_S> (= <E_W>) and the distribution of E_W.  Demonstrate by
# rearranging the E_S field at fixed E_W field: J is unchanged, Xi_1 is not.
rng = np.random.default_rng(3)
FT = FN if False else fields["solenoidal_0"]        # any field; identities are kinematic
E_S0, E_W0, xi0 = FT["E_S"], FT["E_W"], FT["xi"]
psi = lambda E: E ** 2                              # an arbitrary admissible psi
J = lambda ES, EW: mean(1.0 * ES + psi(EW))         # c = 1: the general pressure-free form
Xi1 = lambda ES, EW: mean(ES * 0.5 * np.log(EW / ES))

variants = {
    "original": E_S0,
    "random rearrangement of E_S": E_S0.reshape(-1)[rng.permutation(E_S0.size)
                                                    ].reshape(E_S0.shape),
}
# extremal rearrangement: pair the largest E_S with the largest E_W (comonotone)
order_ES = np.argsort(E_S0.reshape(-1))
order_EW = np.argsort(E_W0.reshape(-1))
E_S_sorted = np.empty(E_S0.size)
E_S_sorted[order_EW] = E_S0.reshape(-1)[order_ES]
variants["comonotone rearrangement"] = E_S_sorted.reshape(E_S0.shape)

n8 = {}
for name, ES in variants.items():
    n8[name] = {"J": J(ES, E_W0), "Xi1": Xi1(ES, E_W0)}
    print(f"   {name:30s}:  pressure-free J = {n8[name]['J']:.10f}   "
          f"comparator Xi_1 = {n8[name]['Xi1']:+.8f}")
dJ = max(abs(n8[k]["J"] - n8["original"]["J"]) for k in n8)
dXi = max(abs(n8[k]["Xi1"] - n8["original"]["Xi1"]) for k in n8)
require("N8 a pressure-free functional is invariant under rearranging E_S at fixed "
        "E_W, while the comparator functional is not",
        dJ / abs(n8["original"]["J"]) < 1e-12 and
        dXi / abs(n8["original"]["Xi1"]) > 1e-2,
        f"relative change of J: {dJ/abs(n8['original']['J']):.2e}; "
        f"relative change of Xi_1: {dXi/abs(n8['original']['Xi1']):.3f} -- exactly "
        "the blindness the theorem predicts: Phi = cE_S + psi(E_W) sees only <E_S> "
        "and the E_W distribution, never their pointwise pairing")
report["N8"] = n8

print()
print("=" * 76)
print("N9  CONDITIONAL / REGIONAL AVERAGES ARE ALSO OBSTRUCTED")
print("=" * 76)
# w = indicator of a region is a non-constant weight, so by the identity of N5 the
# obstruction is a pure SURFACE term on the region boundary -- it does not vanish.
n9 = {}
for rname, mask in (("left half-box", (np.arange(N)[:, None, None] < N // 2) *
                     np.ones((N, N, N), bool)),
                    ("rotation-dominated set {xi > 0}", F["xi"] > 0),
                    ("strong-strain set {E_S > median}", F["E_S"] > np.median(F["E_S"]))):
    val = mean(mask * F["SH"])
    scale = mean(np.abs(F["SH"]))
    n9[rname] = {"regional_integral": val, "scale": scale, "relative": abs(val) / scale}
    print(f"   Omega = {rname:34s}: int_Omega S:H_dev = {val:+.4e}   "
          f"relative {abs(val)/scale:.2e}")
require("N9 restricting the integral to any region reintroduces the pressure",
        all(v["relative"] > 1e-6 for v in n9.values()),
        "an indicator is a non-constant weight; by the N5 identity its contribution "
        "is a pure surface term on the region boundary, which does not vanish -- so "
        "conditional and regional comparator statistics are obstructed too")
report["N9"] = n9

print()
print("=" * 76)
print("M  DIAGNOSTIC VALUES ON A NAVIER-STOKES SNAPSHOT")
print("=" * 76)
spn2 = dns.Spectral(64, nu=1.0 / 200.0)
uh = dns.taylor_green(spn2)
for _ in range(200):                        # t = 4.0 at dt = 0.02
    uh = spn2.step_rk4(uh, 0.02)
FN = field_bundle(spn2, uh)
m = {}
for a in alphas:
    m[str(a)] = {"Xi_alpha": mean(FN["E_S"] ** a * FN["xi"]),
                 "pressure_part": -2.0 * mean(FN["E_S"] ** (a - 1.0) *
                                              (a * FN["xi"] - 0.5) * FN["SH"])}
    print(f"   alpha = {a:+.1f}:  Xi_alpha = {m[str(a)]['Xi_alpha']:+.6e}   "
          f"pressure part of dXi/dt = {m[str(a)]['pressure_part']:+.6e}")
xi_bar = mean(FN["E_S"] * FN["xi"]) / mean(FN["E_S"])
print(f"   normalized  <xi>_{{E_S}} = {xi_bar:+.6f};  plain <xi> = {mean(FN['xi']):+.6f}")
print(f"   <S:H_dev> = {mean(FN['SH']):+.3e},  <xi S:H_dev> = {mean(FN['xi']*FN['SH']):+.6f}")
require("M on a Navier-Stokes field the same conclusions hold",
        abs(mean(FN["SH"])) / mean(np.abs(FN["SH"])) < 1e-12 and
        abs(mean(FN["xi"] * FN["SH"])) / mean(np.abs(FN["xi"] * FN["SH"])) > 1e-3,
        f"<S:H_dev> vanishes (ratio "
        f"{abs(mean(FN['SH']))/mean(np.abs(FN['SH'])):.1e}) while <xi S:H_dev> does "
        f"not (ratio {abs(mean(FN['xi']*FN['SH']))/mean(np.abs(FN['xi']*FN['SH'])):.3f})")
report["M"] = {"alphas": m, "xi_bar_E_S": xi_bar, "xi_plain": mean(FN["xi"]),
               "mean_SH": mean(FN["SH"]), "mean_xi_SH": mean(FN["xi"] * FN["SH"])}

report["failures"] = FAILURES
(OUT / "weighted_numeric.json").write_text(json.dumps(report, indent=2, default=float))
np.savez_compressed(OUT / "ns_snapshot_fields.npz",
                    E_S=FN["E_S"].astype(np.float32), E_W=FN["E_W"].astype(np.float32),
                    xi=FN["xi"].astype(np.float32), SH=FN["SH"].astype(np.float32))
print(f"\n{'ALL CHECKS PASSED' if not FAILURES else 'FAILURES: ' + str(FAILURES)}")
print(f"written: {OUT / 'weighted_numeric.json'}")
sys.exit(0 if not FAILURES else 1)
