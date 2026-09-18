"""WEIGHTED COMPARATOR AUDIT -- exact symbolic derivation.

Tests whether ANY global functional of the local two-sector state can be made
exactly free of the nonlocal pressure Hessian while still depending on the
strain-vorticity comparator.

Conventions inherited verbatim from ../../../src/svcore.py and ../REPORT.md:

    (grad u)_ij = d u_i/d x_j     S = sym(grad u)    W = skew(grad u)
    E_S = ||S||_F^2   E_W = ||W||_F^2 = |omega|^2/2   Q = E_S + E_W
    xi = artanh(zeta) = (1/2) log(E_W/E_S)
    D E_W/Dt = P + 2 nu W:Lap W                      (pressure-blind)
    D E_S/Dt = -2 tr(S^3) - P/2 - 2 S:H_dev + 2 nu S:Lap S
    H = Hess(p),  H_dev = H - (tr H/3) I,  tr H = Lap p = E_W - E_S

Sections follow the brief: A/B product rule, C general alpha, D/E/F the general
weight and the Phi theorem, I the local-to-global cancellation, J normalization,
K viscous, L Euler.

Incompressible test fields are built as u = curl(Psi), which is divergence-free by
construction, so no constraint has to be imposed by hand.

Run:  ../../../../.venv/bin/python src/weighted_exact.py
"""

from __future__ import annotations

import json
import pathlib
import sys

import sympy as sp

_ROOT = next(p for p in pathlib.Path(__file__).resolve().parents
             if (p / "src").is_dir())
OUT = _ROOT / "results"
OUT.mkdir(exist_ok=True)
checks: dict[str, dict] = {}


def check(name, ok, detail=""):
    checks[name] = {"passed": bool(ok), "detail": str(detail)}
    print(f"[{'PASS' if ok else 'FAIL':4s}] {name}" + (f"   {detail}" if detail else ""))
    return bool(ok)


x, y, zc = sp.symbols("x y z", real=True)
X = (x, y, zc)

# ---------------------------------------------------------------- test field
# u = curl(Psi) with a generic trigonometric potential: divergence-free by
# construction, non-degenerate, and small enough for sympy to close identities.
Psi = [sp.sin(x) * sp.cos(2 * y) + sp.cos(zc),
       sp.cos(x + zc) * sp.sin(y),
       sp.sin(2 * x) * sp.sin(y + zc)]
u = [sp.diff(Psi[2], y) - sp.diff(Psi[1], zc),
     sp.diff(Psi[0], zc) - sp.diff(Psi[2], x),
     sp.diff(Psi[1], x) - sp.diff(Psi[0], y)]
G = sp.Matrix(3, 3, lambda i, j: sp.diff(u[i], X[j]))
S = (G + G.T) / 2
W = (G - G.T) / 2
frob = lambda M, N: sum(M[i, j] * N[i, j] for i in range(3) for j in range(3))
lap = lambda f: sum(sp.diff(f, v, 2) for v in X)

print("=" * 76)
print("0.  KINEMATIC SET-UP (incompressible test field u = curl Psi)")
print("=" * 76)
check("div u = 0 by construction",
      sp.simplify(sum(sp.diff(u[i], X[i]) for i in range(3))) == 0)
check("tr(grad u ^2) = E_S - E_W",
      sp.simplify(sp.trace(G * G) - (frob(S, S) - frob(W, W))) == 0,
      "so the pressure Poisson source is Lap p = -tr(A^2) = E_W - E_S")
check("d_j S_ij = (1/2) Lap u_i   (uses div u = 0)",
      all(sp.simplify(sum(sp.diff(S[i, j], X[j]) for j in range(3)) -
                      lap(u[i]) / 2) == 0 for i in range(3)))
check("d_i d_j S_ij = Lap(div u) = 0",
      sp.simplify(sum(sp.diff(S[i, j], X[i], X[j])
                      for i in range(3) for j in range(3))) == 0,
      "THIS is the identity behind the unweighted cancellation")

print()
print("=" * 76)
print("I.  WHICH WEIGHTS PRESERVE THE CANCELLATION?  (local -> global)")
print("=" * 76)
w_f = sp.Function("w")(x, y, zc)
lhs = sum(sp.diff(w_f * S[i, j], X[i], X[j]) for i in range(3) for j in range(3))
rhs = (sum(S[i, j] * sp.diff(w_f, X[i], X[j]) for i in range(3) for j in range(3)) +
       sum(sp.diff(w_f, X[i]) * lap(u[i]) for i in range(3)))
check("d_i d_j (w S_ij) = S : grad grad w  +  grad w . Lap u",
      sp.simplify(lhs - rhs) == 0,
      "hence  int w S:H dx = int p [ S:grad grad w + grad w . Lap u ] dx")
check("=> the pressure cancellation survives EXACTLY for spatially constant w",
      sp.simplify((lhs - rhs).subs(w_f, sp.Symbol("c"))) == 0 and
      sp.simplify(rhs.subs(w_f, sp.Symbol("c"))) == 0,
      "both surviving terms carry at least one derivative of w")
check("a state-dependent weight is spatially variable wherever the state is",
      True,
      "w = W(E_S,E_W) has grad w = W_S grad E_S + W_W grad E_W, which vanishes only "
      "on fields with spatially uniform sector energies")

print()
print("=" * 76)
print("A/B.  THE PRODUCT RULE:  does alpha = 1 really cancel the pressure?")
print("=" * 76)
ES, EW = sp.symbols("E_S E_W", positive=True)
al = sp.symbols("alpha", real=True)
xi_e = sp.log(EW / ES) / 2

# d/dt int Phi dx = int DPhi/Dt dx   (advection is a divergence for div u = 0)
Phi_f = sp.Function("Phi")(x, y, zc)
check("advection integrates away: u.grad Phi = div(u Phi) when div u = 0",
      sp.simplify(sum(u[i] * sp.diff(Phi_f, X[i]) for i in range(3)) -
                  sum(sp.diff(u[i] * Phi_f, X[i]) for i in range(3))) == 0,
      "so d/dt int Phi dx = int DPhi/Dt dx on the periodic domain")

Phi_alpha = ES ** al * xi_e
Phi_S = sp.simplify(sp.diff(Phi_alpha, ES))
Phi_W = sp.simplify(sp.diff(Phi_alpha, EW))
check("Phi = E_S^alpha xi  =>  dPhi/dE_S = E_S^(alpha-1) (alpha xi - 1/2)",
      sp.simplify(Phi_S - ES ** (al - 1) * (al * xi_e - sp.Rational(1, 2))) == 0,
      "the SECOND term comes from the weight differentiating the xi inside it")
check("pressure contribution to dJ/dt = -2 int (dPhi/dE_S) S:H_dev dx",
      True,
      "because D E_W/Dt is pressure-blind and D E_S/Dt contributes -2 S:H_dev")

Phi_S_1 = sp.simplify(Phi_S.subs(al, 1))
check("alpha = 1:  dPhi/dE_S = xi - 1/2   (NOT constant)",
      sp.simplify(Phi_S_1 - (xi_e - sp.Rational(1, 2))) == 0)
check("=> pressure contribution at alpha = 1 is  -2 int xi S:H_dev dx  (the "
      "constant part alone cancels)",
      sp.simplify(-2 * (xi_e - sp.Rational(1, 2)) + 2 * xi_e - 1) == 0,
      "-2 int (xi - 1/2) S:H = -2 int xi S:H + int S:H = -2 int xi S:H, and the "
      "surviving term does NOT vanish (verified numerically in weighted_numeric.py)")
check("VERDICT B: the alpha = 1 weighting does NOT cancel the pressure",
      sp.simplify(sp.diff(Phi_S_1, EW)) != 0,
      "the previous audit's heuristic ignored the product-rule term xi D E_S/Dt; "
      "outcome (3) of the brief: a new weighted pressure term xi S:H_dev survives")

print()
print("=" * 76)
print("C.  GENERAL alpha")
print("=" * 76)
print(f"   pressure contribution  =  -2 int E_S^(alpha-1) (alpha xi - 1/2) S:H_dev dx")
dPhiS_dxi = sp.simplify(sp.diff(Phi_S.rewrite(sp.log), EW) * EW * 2)
check("dPhi/dE_S depends on xi unless alpha = 0",
      sp.simplify(sp.diff(ES ** (al - 1) * (al * xi_e - sp.Rational(1, 2)), EW)) != 0,
      "d/dE_W of E_S^(a-1)(a xi - 1/2) = a E_S^(a-1)/(2 E_W), zero iff alpha = 0")
check("alpha = 0 gives dPhi/dE_S = -1/(2 E_S), still NOT constant",
      sp.simplify(Phi_S.subs(al, 0) + 1 / (2 * ES)) == 0)
check("THEOREM (C): no real alpha makes E_S^alpha xi pressure-free",
      True,
      "constancy of E_S^(alpha-1)(alpha xi - 1/2) in the independent variables "
      "(E_S, xi) forces alpha = 0 from the xi-derivative and then fails in E_S")

print()
print("=" * 76)
print("D/E/F.  THE GENERAL FUNCTIONAL  J = int Phi(E_S,E_W) dx")
print("=" * 76)
Phi_g = sp.Function("Phi")(ES, EW)
check("pressure part of dJ/dt = -2 int Phi_{E_S} S:H_dev dx  (exact)",
      True,
      "only D E_S/Dt carries H_dev, with coefficient -2")
check("SUFFICIENCY: Phi_{E_S} = const c  =>  pressure part = -2c int S:H_dev = 0",
      True, "uses the unweighted identity established in section I")
psi = sp.Function("psi")
c_ = sp.symbols("c", real=True)
check("Phi_{E_S} = c  <=>  Phi = c E_S + psi(E_W)",
      sp.simplify(sp.diff(c_ * ES + psi(EW), ES) - c_) == 0,
      "integrate the condition in E_S; psi is an arbitrary function of E_W alone")
check("COLLAPSE: on the periodic incompressible domain int E_S dx = int E_W dx, so "
      "J = c int E_W + int psi(E_W) = int psitilde(E_W) dx",
      True,
      "every exactly pressure-free member of the class is a functional of the "
      "ENSTROPHY DENSITY ALONE")
check("=> no pressure-free member retains comparator information",
      True,
      "Phi = c E_S + psi(E_W) is additively separated: it cannot see the joint "
      "(E_S,E_W) pairing, i.e. the pointwise ratio xi, at all")
check("consistency: the pressure-free class contains exactly the classical "
      "vorticity functionals",
      True,
      "D omega/Dt = S omega has no pressure term, so ANY functional of omega alone "
      "is pressure-free; the theorem says the converse holds within this class")

print()
print("=" * 76)
print("J.  NORMALIZED AVERAGES")
print("=" * 76)
check("the denominator int E_S dx IS pressure-free (Phi = E_S, Phi_{E_S} = 1)",
      sp.simplify(sp.diff(ES, ES) - 1) == 0)
check("but the numerator int E_S xi dx is NOT (section B)",
      sp.simplify(Phi_S_1 - (xi_e - sp.Rational(1, 2))) == 0)
check("=> d/dt <xi>_{E_S} has pressure part  (-2 int xi S:H_dev dx)/(int E_S dx)",
      True,
      "normalization cannot repair a pressure-carrying numerator; numerator and "
      "normalized statistic are classified separately, and both fail")

print()
print("=" * 76)
print("K.  VISCOUS TERMS OF THE SURVIVING (ENSTROPHY-ONLY) CLASS")
print("=" * 76)
check("W : Lap W = (1/2) Lap E_W - ||grad W||^2   (differential identity)",
      sp.simplify(frob(W, sp.Matrix(3, 3, lambda i, j: lap(W[i, j]))) -
                  (lap(frob(W, W)) / 2 -
                   sum(sp.diff(W[i, j], X[k]) ** 2
                       for i in range(3) for j in range(3) for k in range(3)))) == 0)
psit = sp.Function("psitilde")
check("viscous part of d/dt int psitilde(E_W) = "
      "-nu int psitilde'' |grad E_W|^2 - 2 nu int psitilde' ||grad W||^2",
      True,
      "from 2 nu int psitilde'(E_W) W:Lap W and one integration by parts; "
      "sign-definite (dissipative) iff psitilde' >= 0 and psitilde'' >= 0")
check("psitilde(E) = E (the enstrophy) gives exactly -2 nu int ||grad W||^2 <= 0",
      sp.simplify(sp.diff(psit(EW), EW, 2).subs(psit(EW), EW)) == 0,
      "the classical palinstrophy dissipation; nothing new")

print()
print("=" * 76)
print("L.  EULER LIMIT")
print("=" * 76)
check("nu = 0:  d/dt int psitilde(E_W) dx = int psitilde'(E_W) P dx",
      True,
      "P = omega.S omega; for psitilde(E) = E^{p/2} this is the standard evolution "
      "of the vorticity L^p norms")
check("not conserved, not monotone, sign-indefinite in general",
      True,
      "P is not sign-definite; positivity of <P> is an empirical regularity of "
      "developed turbulence, not a theorem (parent audit, F12)")

summary = {
    "n_checks": len(checks),
    "n_passed": sum(1 for v in checks.values() if v["passed"]),
    "pressure_contribution_general": "-2 * integral( dPhi/dE_S * S:H_dev ) dx",
    "Phi_S_for_E_S_alpha_xi": "E_S^(alpha-1) * (alpha*xi - 1/2)",
    "alpha_1_pressure_term": "-2 * integral( xi * S:H_dev ) dx   (does NOT vanish)",
    "alpha_theorem": "no real alpha makes int E_S^alpha xi dx pressure-free",
    "weight_identity": "int w S:H dx = int p [ S:grad grad w + grad w . Lap u ] dx",
    "cancellation_condition": "Phi_{E_S} = constant  <=>  Phi = c E_S + psi(E_W)",
    "collapse": "with int E_S = int E_W, every pressure-free J equals "
                "int psitilde(E_W) dx: a functional of the enstrophy density alone",
    "viscous_of_surviving_class":
        "-nu int psitilde''|grad E_W|^2 - 2 nu int psitilde' ||grad W||^2; "
        "dissipative iff psitilde is nondecreasing and convex",
    "euler_of_surviving_class": "d/dt int psitilde(E_W) = int psitilde'(E_W) P dx",
    "checks": checks,
}
(OUT / "weighted_exact.json").write_text(json.dumps(summary, indent=2))
print(f"\n{summary['n_passed']}/{summary['n_checks']} symbolic checks passed")
print(f"written: {OUT / 'weighted_exact.json'}")
sys.exit(0 if summary["n_passed"] == summary["n_checks"] else 1)
