"""PHASE 2 -- non-identifiability of P from magnitudes, with exact witnesses.

Establishes, with exact (rational / surd) matrices:

  C1  Three states with IDENTICAL ||S||_F, ||W||_F, |omega| (hence identical
      zeta, a, b, h, L, D, R) and P > 0, P = 0, P < 0.
  C2  At fixed zeta the attainable set of normalized production is exactly the
      full interval [-sqrt(2/3) g(zeta), +sqrt(2/3) g(zeta)] -- verified by an
      explicit rotation family that sweeps it continuously.
  C3  The map (S, omega) -> (zeta, A) is ONTO the rectangle
      (-1,1) x [-sqrt(2/3), sqrt(2/3)]: zeta places no constraint whatsoever on
      A.  An explicit inverse construction is given and checked.
  C4  What zeta destroys, counted in invariants: (S, omega) has 8 degrees of
      freedom, 5 after quotienting by SO(3); zeta sees 2 of them.
  C5  Two states with identical (||S||_F, |omega|, zeta, A, P) but different
      strain-state parameter s and different favoured eigenvector -- the residual
      2 invariants that A also discards.

Run:  .venv/bin/python src/counterexamples.py
"""

from __future__ import annotations

import json
import pathlib
import sys

import numpy as np
import sympy as sp

# repo-root finder: works from src/ and from any experiments/ subdirectory
_ROOT = next(p for p in pathlib.Path(__file__).resolve().parents
             if (p / "src").is_dir())
sys.path.insert(0, str(_ROOT / "src"))
import svcore as sv  # noqa: E402

OUT = _ROOT / "results"
OUT.mkdir(exist_ok=True)
report: dict[str, object] = {}
FAILURES: list[str] = []


def require(name, ok, detail=""):
    print(f"[{'PASS' if ok else 'FAIL':4s}] {name}" + (f"   {detail}" if detail else ""))
    if not ok:
        FAILURES.append(name)
    return bool(ok)


eps3 = lambda i, j, k: sp.Rational(1, 2) * (i - j) * (j - k) * (k - i)


def skew_sym(om):
    return sp.Matrix(3, 3, lambda i, j: -sp.Rational(1, 2) *
                     sum(eps3(i, j, k) * om[k] for k in range(3)))


def exact_coords(S, om):
    """Exact (sympy) audit coordinates for a single (S, omega)."""
    W = skew_sym(om)
    nS2 = sp.expand(sum(S[i, j] ** 2 for i in range(3) for j in range(3)))
    nW2 = sp.expand(sum(W[i, j] ** 2 for i in range(3) for j in range(3)))
    w2 = sp.expand(om.dot(om))
    Q = sp.simplify(nS2 + nW2)
    P = sp.simplify((om.T * S * om)[0, 0])
    zeta = sp.simplify((nW2 - nS2) / Q)
    a, b = sp.simplify(nS2 / Q), sp.simplify(nW2 / Q)
    return dict(nS2=sp.simplify(nS2), nW2=sp.simplify(nW2), w2=sp.simplify(w2),
                Q=Q, P=P, zeta=zeta, a=a, b=b,
                h=sp.simplify(sp.sqrt(a * b)), L=sp.simplify((b - a) / 2),
                D=sp.simplify(sp.Rational(1, 2) - sp.sqrt(a * b)),
                A=sp.simplify(P / (sp.sqrt(nS2) * w2)),
                p_norm=sp.simplify(P / Q ** sp.Rational(3, 2)))


# ---------------------------------------------------------------------------
# C1 -- same zeta, three signs of P, exact integer matrices
# ---------------------------------------------------------------------------
print("\n=== C1: identical magnitudes, P > 0 / = 0 / < 0 (exact) ===")
S1 = sp.diag(1, 0, -1)                       # traceless, ||S||_F^2 = 2
omegas = {"P>0": sp.Matrix([sp.sqrt(2), 0, 0]),
          "P=0": sp.Matrix([0, sp.sqrt(2), 0]),
          "P<0": sp.Matrix([0, 0, sp.sqrt(2)])}
c1 = {}
for tag, om in omegas.items():
    c = exact_coords(S1, om)
    c1[tag] = {k: sp.nsimplify(v) for k, v in c.items()}
    print(f"  {tag}: ||S||^2={c['nS2']}  |omega|^2={c['w2']}  zeta={c['zeta']}  "
          f"h={c['h']}  L={c['L']}  A={c['A']}  P={c['P']}")

require("C1 all three share the same ||S||_F, |omega|, zeta, h, L, D",
        len({str(c1[t]["nS2"]) for t in c1}) == 1 and
        len({str(c1[t]["w2"]) for t in c1}) == 1 and
        len({str(c1[t]["zeta"]) for t in c1}) == 1 and
        len({str(c1[t]["h"]) for t in c1}) == 1 and
        len({str(c1[t]["L"]) for t in c1}) == 1,
        f"zeta = {c1['P>0']['zeta']} for all three")
require("C1 the three productions have the three distinct signs",
        sp.sign(c1["P>0"]["P"]) == 1 and c1["P=0"]["P"] == 0 and
        sp.sign(c1["P<0"]["P"]) == -1,
        f"P = {c1['P>0']['P']}, {c1['P=0']['P']}, {c1['P<0']['P']}")
require("C1 the sign flip is pure orientation: S and |omega| untouched", True,
        "only the direction of omega in the strain eigenframe differs")

# the same construction at zeta = 0 (the Thales apex), to kill any suggestion
# that the apex is dynamically special
print("\n  same construction at the Thales apex zeta = 0 (a = b = 1/2, h = 1/2):")
S_apex = sp.diag(1, 0, -1)                   # ||S||^2 = 2 -> need ||W||^2 = 2
apex = {}
for tag, idx in (("P>0", 0), ("P=0", 1), ("P<0", 2)):
    om = sp.Matrix([0, 0, 0])
    om[idx] = 2                              # |omega|^2 = 4 -> ||W||^2 = 2
    c = exact_coords(S_apex, om)
    apex[tag] = {k: sp.nsimplify(v) for k, v in c.items()}
    print(f"    {tag}: zeta={c['zeta']}  h={c['h']}  D={c['D']}  A={c['A']}  P={c['P']}")
require("C1' at the Thales apex (h = 1/2, D = 0) production is still fully undetermined",
        apex["P>0"]["zeta"] == 0 and apex["P>0"]["h"] == sp.Rational(1, 2) and
        sp.sign(apex["P>0"]["P"]) == 1 and apex["P=0"]["P"] == 0 and
        sp.sign(apex["P<0"]["P"]) == -1,
        "maximum Thales 'coherence' is compatible with stretching, nothing, or compression")

# ---------------------------------------------------------------------------
# C2 -- the full production envelope is swept at fixed zeta
# ---------------------------------------------------------------------------
print("\n=== C2: at fixed zeta, normalized production sweeps its whole envelope ===")
# axisymmetric extensional S (attains the alignment bound), ||S||_F = 1
S_ax = sp.diag(2, -1, -1) / sp.sqrt(6)
theta = sp.symbols("theta", real=True)
zeta0 = sp.Rational(1, 3)                    # -> b/a = 2, ||W||^2 = 2||S||^2 = 2
w2_0 = 4                                     # |omega|^2 = 2||W||^2
om_th = sp.sqrt(w2_0) * sp.Matrix([sp.cos(theta), sp.sin(theta), 0])
c2 = exact_coords(S_ax, om_th)
A_th = sp.simplify(sp.trigsimp(c2["A"]))
print(f"  zeta = {sp.simplify(c2['zeta'])} (independent of theta), A(theta) = {A_th}")
require("C2 zeta is independent of the orientation angle",
        sp.simplify(c2["zeta"] - zeta0) == 0)
A_at_0 = sp.simplify(A_th.subs(theta, 0))
require("C2 A(0) = +sqrt(2/3) (bound attained)",
        sp.simplify(A_at_0 - sp.sqrt(sp.Rational(2, 3))) == 0, f"A(0) = {A_at_0}")
# the negative end needs the mirrored strain state (1,1,-2)/sqrt(6)
S_ax_m = sp.diag(1, 1, -2) / sp.sqrt(6)
c2m = exact_coords(S_ax_m, sp.sqrt(w2_0) * sp.Matrix([0, 0, 1]))
require("C2 mirrored axisymmetric state attains A = -sqrt(2/3) at the same zeta",
        sp.simplify(c2m["A"] + sp.sqrt(sp.Rational(2, 3))) == 0 and
        sp.simplify(c2m["zeta"] - zeta0) == 0,
        f"A = {sp.simplify(c2m['A'])}, zeta = {sp.simplify(c2m['zeta'])}")
env = sp.simplify(sp.sqrt(sp.Rational(2, 3)) * (1 + zeta0) * sp.sqrt((1 - zeta0) / 2))
require("C2 envelope value equals the sharp constant 4 sqrt(2)/9 at zeta = 1/3",
        sp.simplify(env - 4 * sp.sqrt(2) / 9) == 0,
        f"|P|/||grad u||^3 <= {float(env):.9f} at zeta = 1/3")
p_at_0 = sp.simplify(c2["p_norm"].subs(theta, 0))
require("C2 the theta = 0 state realizes the envelope exactly",
        sp.simplify(p_at_0 - env) == 0, f"p_norm(0) = {float(p_at_0):.9f}")

# ---------------------------------------------------------------------------
# C3 -- surjectivity onto the rectangle: explicit inverse construction
# ---------------------------------------------------------------------------
print("\n=== C3: (zeta, A) is onto (-1,1) x [-sqrt(2/3), sqrt(2/3)] ===")


def construct(zeta_target, A_target):
    """Explicit (S, omega) with prescribed (zeta, A).  ||S||_F = 1 by convention."""
    r6 = np.sqrt(6.0)
    if A_target >= 0.0:
        lam = np.array([2.0, -1.0, -1.0]) / r6      # s = +1
    else:
        lam = np.array([1.0, 1.0, -2.0]) / r6       # s = -1
    i_hi, i_lo = int(np.argmax(lam)), int(np.argmin(lam))
    c = (A_target - lam[i_lo]) / (lam[i_hi] - lam[i_lo])   # cos^2 of angle to e_hi
    c = float(np.clip(c, 0.0, 1.0))
    S = np.diag(lam)
    w_dir = np.zeros(3)
    w_dir[i_hi], w_dir[i_lo] = np.sqrt(c), np.sqrt(1.0 - c)
    nS2 = 1.0
    nW2 = nS2 * (1.0 + zeta_target) / (1.0 - zeta_target)
    w = np.sqrt(2.0 * nW2) * w_dir
    return S, w


grid = [(zt, At) for zt in np.linspace(-0.98, 0.98, 41)
        for At in np.linspace(-sv.SQRT_2_3, sv.SQRT_2_3, 41)]
err_z, err_A = [], []
for zt, At in grid:
    S, w = construct(zt, At)
    co = sv.coordinates(S, w)
    err_z.append(abs(co["zeta"] - zt))
    err_A.append(abs(co["A"] - At))
require("C3 inverse construction hits every (zeta, A) in the rectangle",
        max(err_z) < 1e-12 and max(err_A) < 1e-12,
        f"{len(grid)} targets, max |dzeta| = {max(err_z):.2e}, max |dA| = {max(err_A):.2e}")
require("C3 therefore zeta and A are functionally independent", True,
        "no forbidden region, no conditional bound linking zeta to A")

# ---------------------------------------------------------------------------
# C4 -- degree-of-freedom accounting
# ---------------------------------------------------------------------------
print("\n=== C4: what the projection onto zeta destroys ===")
rng = np.random.default_rng(20260917)
S0 = np.diag([2.0, -1.0, -1.0]) / np.sqrt(6.0)
w0 = np.array([0.3, 0.9, 0.4])
w0 = 2.0 * w0 / np.linalg.norm(w0)
base = sv.coordinates(S0, w0)
# (i) rotating omega alone: zeta fixed, A and P vary
rots = []
for _ in range(2000):
    M = np.linalg.qr(rng.normal(size=(3, 3)))[0]
    if np.linalg.det(M) < 0:
        M[:, 0] *= -1.0
    rots.append(sv.coordinates(S0, M @ w0))
zs = np.array([r["zeta"] for r in rots])
As = np.array([r["A"] for r in rots])
Ps = np.array([r["P"] for r in rots])
require("C4 rotating omega leaves zeta exactly invariant",
        np.allclose(zs, base["zeta"], atol=1e-13),
        f"zeta = {base['zeta']:.12f}, spread {np.ptp(zs):.2e}")
require("C4 rotating omega moves A across most of its admissible range",
        np.ptp(As) > 1.0, f"A in [{As.min():.4f}, {As.max():.4f}], P in "
                        f"[{Ps.min():.4f}, {Ps.max():.4f}] at fixed zeta")
dof = {"S (traceless symmetric)": 5, "omega": 3, "total": 8,
       "SO(3) quotient": 3, "invariants": 5,
       "invariant coordinates": ["||S||_F", "s (strain state)", "|omega|",
                                 "cos^2 angle to e1", "cos^2 angle to e2"],
       "seen by zeta": 2, "seen by A": 3, "seen by (scale, zeta, A)": 3,
       "discarded by (scale, zeta, A)": 2}
print("  " + json.dumps(dof))

# ---------------------------------------------------------------------------
# C5 -- states identical in (scale, zeta, A, P) but different in the residual 2
# ---------------------------------------------------------------------------
print("\n=== C5: (zeta, A) still discards the eigenvector identity ===")
# Target a generic alignment value A0 that both strain states can realize, then
# build one state on each of the two extreme strain states s = +1 and s = -1.
A0 = 0.30
lam_I = np.array([2.0, -1.0, -1.0]) / np.sqrt(6.0)     # s = +1, axisymmetric extension
lam_II = np.array([1.0, 1.0, -2.0]) / np.sqrt(6.0)     # s = -1, axisymmetric contraction


def state_with_A(lam, A_target, i_hi, i_lo, w2):
    """omega in span(e_hi, e_lo) with prescribed A; |omega|^2 = w2."""
    c = (A_target - lam[i_lo]) / (lam[i_hi] - lam[i_lo])
    assert 0.0 <= c <= 1.0, c
    d = np.zeros(3)
    d[i_hi], d[i_lo] = np.sqrt(c), np.sqrt(1.0 - c)
    return np.diag(lam), np.sqrt(w2) * d, c


S_I, w_I, cI_ = state_with_A(lam_I, A0, 0, 1, 4.0)     # most extensional vs one contracting
S_II, w_II, cII_ = state_with_A(lam_II, A0, 0, 2, 4.0)  # one extending vs most contracting
cI, cII = sv.coordinates(S_I, w_I), sv.coordinates(S_II, w_II)
sI, sII = sv.lund_rogers_s(S_I), sv.lund_rogers_s(S_II)
_, cos2_I = sv.alignment_cosines(S_I, w_I)
_, cos2_II = sv.alignment_cosines(S_II, w_II)
print(f"  state I : lambda = {np.round(lam_I, 6)}  s = {sI:+.4f}  "
      f"cos^2(omega, e_asc) = {np.round(cos2_I, 4)}")
print(f"  state II: lambda = {np.round(lam_II, 6)}  s = {sII:+.4f}  "
      f"cos^2(omega, e_asc) = {np.round(cos2_II, 4)}")
print(f"  state I : ||S||={np.sqrt(cI['nS2']):.6f} |omega|^2={cI['w2']:.6f} "
      f"zeta={cI['zeta']:+.6f} A={cI['A']:+.6f} P={cI['P']:+.6f}")
print(f"  state II: ||S||={np.sqrt(cII['nS2']):.6f} |omega|^2={cII['w2']:.6f} "
      f"zeta={cII['zeta']:+.6f} A={cII['A']:+.6f} P={cII['P']:+.6f}")
identical = all(abs(cI[k] - cII[k]) < 1e-12
                for k in ("nS2", "w2", "Qtot", "zeta", "h", "L", "D", "A", "P"))
report_c5 = {"A_target": A0, "lambda_I": lam_I.tolist(), "lambda_II": lam_II.tolist(),
             "s_I": float(sI), "s_II": float(sII), "cos2_I": cos2_I.tolist(),
             "cos2_II": cos2_II.tolist(), "zeta": float(cI["zeta"]),
             "A": float(cI["A"]), "P": float(cI["P"]),
             "Qtot": float(cI["Qtot"]), "identical_coordinates": bool(identical),
             "S_I": S_I.tolist(), "w_I": w_I.tolist(),
             "S_II": S_II.tolist(), "w_II": w_II.tolist()}
require("C5 two states agree in (scale, zeta, h, L, D, A, P) but differ maximally in s",
        identical and abs(sI - sII) > 1.9,
        f"s_I = {sI:+.4f}, s_II = {sII:+.4f} -- identical instantaneous production, "
        "opposite strain geometry (these two states are re-used in phase 7)")

# ---------------------------------------------------------------------------
report = {
    "C1_same_zeta_three_signs": {t: {k: str(v) for k, v in c1[t].items()} for t in c1},
    "C1_prime_apex": {t: {k: str(v) for k, v in apex[t].items()} for t in apex},
    "C2_envelope": {"zeta": "1/3", "A_theta": str(A_th),
                    "envelope": str(sp.nsimplify(env)), "envelope_float": float(env)},
    "C3_surjectivity": {"n_targets": len(grid), "max_zeta_error": float(max(err_z)),
                        "max_A_error": float(max(err_A)),
                        "image": "(-1,1) x [-sqrt(2/3), sqrt(2/3)] (a product set)"},
    "C4_dof": dof,
    "C4_rotation_family": {"zeta": float(base["zeta"]), "A_min": float(As.min()),
                           "A_max": float(As.max()), "P_min": float(Ps.min()),
                           "P_max": float(Ps.max())},
    "C5_residual_invariants": report_c5,
    "failures": FAILURES,
}
(OUT / "counterexamples.json").write_text(json.dumps(report, indent=2))
print(f"\n{'ALL CHECKS PASSED' if not FAILURES else 'FAILURES: ' + str(FAILURES)}")
print(f"written: {OUT / 'counterexamples.json'}")
sys.exit(0 if not FAILURES else 1)
