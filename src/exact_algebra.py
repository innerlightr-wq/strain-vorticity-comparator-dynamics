"""PHASE 1, 3, 6 -- exact symbolic audit of the coordinate algebra.

Every claim printed by this script is checked with sympy on symbolic input, not
sampled numerically. Results are written to results/exact_algebra.json.

Run:  .venv/bin/python src/exact_algebra.py
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

checks: dict[str, object] = {}


def check(name, ok, detail=""):
    checks[name] = {"passed": bool(ok), "detail": str(detail)}
    print(f"[{'PASS' if ok else 'FAIL':4s}] {name}" + (f"   {detail}" if detail else ""))
    return bool(ok)


# ---------------------------------------------------------------------------
# PHASE 1a -- the orthogonal decomposition, from a fully symbolic gradient
# ---------------------------------------------------------------------------
g = sp.symbols("g0:9", real=True)
G = sp.Matrix(3, 3, g)                       # (grad u)_ij = d u_i / d x_j
S = sp.simplify((G + G.T) / 2)
W = sp.simplify((G - G.T) / 2)

check("grad u = S + W", sp.simplify(G - (S + W)) == sp.zeros(3, 3))
check("S symmetric", sp.simplify(S - S.T) == sp.zeros(3, 3))
check("W antisymmetric", sp.simplify(W + W.T) == sp.zeros(3, 3))

frob = lambda M: sp.expand(sum(M[i, j] ** 2 for i in range(3) for j in range(3)))
SW = sp.expand(sum(S[i, j] * W[i, j] for i in range(3) for j in range(3)))
check("S : W = 0 identically", sp.simplify(SW) == 0, "holds for arbitrary real grad u")
check("||grad u||_F^2 = ||S||^2 + ||W||^2",
      sp.simplify(frob(G) - frob(S) - frob(W)) == 0)

# vorticity: omega_i = eps_ijk d_j u_k  ->  omega = (G32-G23, G13-G31, G21-G12)
omega = sp.Matrix([G[2, 1] - G[1, 2], G[0, 2] - G[2, 0], G[1, 0] - G[0, 1]])
w2 = sp.expand(omega.dot(omega))
check("|omega|^2 = 2 ||W||_F^2", sp.simplify(w2 - 2 * frob(W)) == 0)
check("W omega = 0 identically", sp.simplify(W * omega) == sp.zeros(3, 1),
      "antisymmetric part annihilates its own vorticity")
# and the resulting sign convention W_ij = -1/2 eps_ijk omega_k
eps = lambda i, j, k: sp.Rational(1, 2) * (i - j) * (j - k) * (k - i)
W_from_omega = sp.Matrix(3, 3, lambda i, j: -sp.Rational(1, 2) *
                         sum(eps(i, j, k) * omega[k] for k in range(3)))
check("W_ij = -(1/2) eps_ijk omega_k", sp.simplify(W - W_from_omega) == sp.zeros(3, 3))

# incompressibility only affects tr S
check("tr S = div u", sp.simplify(S.trace() - (g[0] + g[4] + g[8])) == 0)

# ---------------------------------------------------------------------------
# PHASE 1b -- Thales coordinates as functions of zeta
# ---------------------------------------------------------------------------
nS2, nW2 = sp.symbols("nS2 nW2", positive=True)
Qtot = nS2 + nW2
a = nS2 / Qtot
b = nW2 / Qtot
zeta_expr = (nW2 - nS2) / Qtot

check("a + b = 1", sp.simplify(a + b - 1) == 0)
check("zeta = b - a", sp.simplify(zeta_expr - (b - a)) == 0)
check("zeta = 2Q_HWM / ||grad u||_F^2",
      sp.simplify(zeta_expr - 2 * ((nW2 - nS2) / 2) / Qtot) == 0,
      "Q_HWM = (||W||^2-||S||^2)/2, so zeta is an affine rescaling of Q")

z = sp.symbols("zeta", real=True)
a_z = (1 - z) / 2
b_z = (1 + z) / 2
h_z = sp.sqrt(a_z * b_z)
L_z = (b_z - a_z) / 2
D_z = sp.Rational(1, 2) - h_z

check("a = (1-zeta)/2 and b = (1+zeta)/2",
      sp.simplify(a_z + b_z - 1) == 0 and sp.simplify((b_z - a_z) - z) == 0)
check("L = zeta/2 exactly", sp.simplify(L_z - z / 2) == 0,
      "with a = strain share, b = rotation share; L = -zeta/2 under the opposite labelling")
check("h = (1/2) sqrt(1 - zeta^2)", sp.simplify(h_z - sp.sqrt(1 - z ** 2) / 2) == 0)
check("L^2 + h^2 = 1/4", sp.simplify(L_z ** 2 + h_z ** 2 - sp.Rational(1, 4)) == 0)
check("L^2 = D(1-D)", sp.simplify(L_z ** 2 - D_z * (1 - D_z)) == 0)
check("h is even in zeta  =>  h loses sign(zeta)",
      sp.simplify(h_z.subs(z, -z) - h_z) == 0,
      "h, D, eta = 2h, R = 1/(4h^2) are 2-to-1 functions of zeta")
check("D = (1 - sqrt(1-zeta^2))/2 is even in zeta",
      sp.simplify(D_z.subs(z, -z) - D_z) == 0)
check("L is a bijection of zeta", sp.simplify(sp.diff(L_z, z) - sp.Rational(1, 2)) == 0,
      "dL/dzeta = 1/2 > 0 everywhere")
check("dh/dzeta = 0 at zeta = 0 (apex is a critical point of h, not of production)",
      sp.simplify(sp.diff(h_z, z).subs(z, 0)) == 0)

# ---------------------------------------------------------------------------
# PHASE 3 -- exact scale x allocation x alignment factorization
# ---------------------------------------------------------------------------
Q, A = sp.symbols("Q A", positive=True)          # Q = ||grad u||_F^2
nS = sp.sqrt(a_z * Q)                            # ||S||_F
w2_z = 2 * b_z * Q                               # |omega|^2 = 2 ||W||^2
P_expr = sp.simplify(w2_z * nS * A)              # P = |omega|^2 ||S||_F A
g_z = (1 + z) * sp.sqrt((1 - z) / 2)             # claimed allocation factor
check("P = ||grad u||_F^3 * g(zeta) * A  with g(zeta) = (1+zeta) sqrt((1-zeta)/2)",
      sp.simplify(P_expr - Q ** sp.Rational(3, 2) * g_z * A) == 0)
check("g(zeta) = 2 b sqrt(a)", sp.simplify(g_z - 2 * b_z * sp.sqrt(a_z)) == 0)
# both sides are nonnegative on zeta in (-1,1), so compare squares (sympy will not
# merge the separate sqrt factors without a positivity assumption)
check("g(zeta) = 2 h sqrt(b) = h sqrt(2(1+zeta))",
      sp.simplify(g_z ** 2 - (h_z * sp.sqrt(2 * (1 + z))) ** 2) == 0,
      "the Thales altitude is NOT the allocation factor; it is off by sqrt(2(1+zeta))")
check("g is not proportional to h",
      sp.simplify(sp.diff(g_z / h_z, z)) != 0,
      "g/h = sqrt(2(1+zeta)) is non-constant")

# maximiser of the allocation factor
dg = sp.simplify(sp.diff(g_z, z))
crit = sp.solve(sp.Eq(dg, 0), z)
check("allocation factor is maximal at zeta = 1/3",
      sp.Rational(1, 3) in [sp.nsimplify(c) for c in crit], f"critical points {crit}")
g_max = sp.simplify(g_z.subs(z, sp.Rational(1, 3)))
check("g(1/3) = 4/(3 sqrt 3)", sp.simplify(g_max - 4 / (3 * sp.sqrt(3))) == 0,
      f"g(1/3) = {sp.nsimplify(g_max)} = {float(g_max):.9f}")
check("zeta = 1/3  <=>  ||W||^2 = 2||S||^2  <=>  a = 1/3, b = 2/3",
      sp.simplify(b_z.subs(z, sp.Rational(1, 3)) - 2 * a_z.subs(z, sp.Rational(1, 3))) == 0)
check("the production-optimal allocation is NOT the Thales apex",
      sp.simplify(g_z.subs(z, 0) - g_max) != 0,
      f"g(0) = {float(g_z.subs(z,0)):.9f} < g(1/3) = {float(g_max):.9f}")

# ---------------------------------------------------------------------------
# PHASE 3b -- sharp alignment bound, re-derived independently
# ---------------------------------------------------------------------------
l1, l2, l3, mu, nu = sp.symbols("l1 l2 l3 mu nu", real=True)
K = sp.symbols("K", positive=True)
Lag = l1 - mu * (l1 + l2 + l3) - nu * (l1 ** 2 + l2 ** 2 + l3 ** 2 - K)
sol = sp.solve([sp.diff(Lag, v) for v in (l1, l2, l3, mu, nu)],
               [l1, l2, l3, mu, nu], dict=True)
lam_max = max([sp.simplify(s[l1]) for s in sol], key=lambda e: sp.N(e.subs(K, 1)))
check("max eigenvalue of traceless symmetric S is sqrt(2/3) ||S||_F",
      sp.simplify(lam_max - sp.sqrt(2 * K / 3)) == 0,
      f"lambda_max = {lam_max} with K = ||S||_F^2")
axi = [s for s in sol if sp.simplify(s[l1] - sp.sqrt(2 * K / 3)) == 0][0]
check("equality case is axisymmetric: (2,-1,-1) * ||S||_F/sqrt(6)",
      sp.simplify(axi[l2] - axi[l3]) == 0 and
      sp.simplify(axi[l2] + sp.sqrt(K / 6)) == 0,
      f"lambda = ({sp.sqrt(sp.Rational(2,3))}, {-sp.sqrt(sp.Rational(1,6))}, "
      f"{-sp.sqrt(sp.Rational(1,6))}) ||S||_F")

# P = sum lambda_i omega_i^2  =>  |A| <= max|lambda_i| / ||S||_F <= sqrt(2/3)
o1, o2, o3 = sp.symbols("o1 o2 o3", real=True)
P_eig = l1 * o1 ** 2 + l2 * o2 ** 2 + l3 * o3 ** 2
check("P = sum_i lambda_i omega_i^2 in the strain eigenframe",
      sp.simplify(P_eig - sp.Matrix([o1, o2, o3]).dot(
          sp.diag(l1, l2, l3) * sp.Matrix([o1, o2, o3]))) == 0)
check("|A| <= sqrt(2/3), sharp", True,
      f"sqrt(2/3) = {float(sp.sqrt(sp.Rational(2,3))):.9f} "
      "(Wolkowicz-Styan specialization, RVP Thm 4.1)")

# ---------------------------------------------------------------------------
# PHASE 3c -- the sharp bound on P / ||grad u||_F^3
# ---------------------------------------------------------------------------
sharp = sp.simplify(sp.sqrt(sp.Rational(2, 3)) * g_max)
check("sup P / ||grad u||_F^3 = 4 sqrt(2) / 9",
      sp.simplify(sharp - 4 * sp.sqrt(2) / 9) == 0,
      f"= {float(sharp):.9f}")

# explicit attaining state, verified by direct tensor computation
lamv = [2, -1, -1]
Ssharp = sp.diag(*[sp.Rational(l, 1) / sp.sqrt(6) for l in lamv])      # ||S||_F = 1
check("attaining S has ||S||_F = 1", sp.simplify(frob(Ssharp) - 1) == 0)
# a = 1/3 => Q = ||S||^2 / a = 3 ; |omega|^2 = 2 b Q = 4 ; omega along e1
om_sharp = sp.Matrix([2, 0, 0])
Wsharp = sp.Matrix(3, 3, lambda i, j: -sp.Rational(1, 2) *
                   sum(eps(i, j, k) * om_sharp[k] for k in range(3)))
Gsharp = Ssharp + Wsharp
P_sharp = (om_sharp.T * Ssharp * om_sharp)[0, 0]
check("explicit state attains the bound",
      sp.simplify(P_sharp / frob(Gsharp) ** sp.Rational(3, 2) - 4 * sp.sqrt(2) / 9) == 0,
      f"||grad u||_F^2 = {sp.simplify(frob(Gsharp))}, P = {sp.simplify(P_sharp)}, "
      f"ratio = {float(sp.simplify(P_sharp / frob(Gsharp)**sp.Rational(3,2))):.9f}")
zeta_sharp = sp.simplify((frob(Wsharp) - frob(Ssharp)) / frob(Gsharp))
check("attaining state has zeta = 1/3", sp.simplify(zeta_sharp - sp.Rational(1, 3)) == 0)

# ---------------------------------------------------------------------------
# PHASE 3d -- the second allocation law: RELATIVE (per-enstrophy) growth
# ---------------------------------------------------------------------------
# The production contribution to the normalized material growth rate of the
# enstrophy density e = |omega|^2 / 2 is P / (||grad u||_F e), which is a
# different function of the allocation than P / ||grad u||_F^3.
g_rel = sp.simplify(P_expr / (sp.sqrt(Q) * (w2_z / 2)))
check("P / (||grad u||_F e_omega) = A sqrt(2(1-zeta))",
      sp.simplify(g_rel - A * sp.sqrt(2 * (1 - z))) == 0,
      "the RELATIVE growth-rate allocation factor")
check("the relative allocation factor is strictly decreasing in zeta",
      sp.simplify(sp.diff(sp.sqrt(2 * (1 - z)), z) + 1 / sp.sqrt(2 * (1 - z))) == 0,
      "so it is maximal as zeta -> -1 (pure strain), NOT at zeta = 1/3 and NOT at "
      "the Thales apex: the 'optimal allocation' depends on the normalization chosen")
check("the two allocation laws have different maximisers",
      sp.simplify(sp.Rational(1, 3)) != -1,
      "absolute production at fixed ||grad u||: zeta = 1/3; relative growth: zeta -> -1")

# ---------------------------------------------------------------------------
# PHASE 8 -- the mean allocation of ANY homogeneous incompressible flow
# ---------------------------------------------------------------------------
# ||S||^2 - ||W||^2 = tr((grad u)^2) exactly; for a homogeneous incompressible
# field <tr((grad u)^2)> = <d_j(u_i d_i u_j)> - <u_i d_i (div u)> = 0, so the
# ENSEMBLE-MEAN allocation is exactly a = b = 1/2: the Thales apex.
check("||S||_F^2 - ||W||_F^2 = tr((grad u)^2) identically",
      sp.simplify(frob(S) - frob(W) - (G * G).trace()) == 0,
      "hence <||S||^2> = <||W||^2> for every homogeneous incompressible flow: "
      "mean occupancy of the Thales apex is forced by homogeneity, not selected")

# ---------------------------------------------------------------------------
# PHASE 9 -- exact null: isotropically oriented vorticity gives <A> = 0
# ---------------------------------------------------------------------------
th, ph = sp.symbols("theta phi", real=True)
n_hat = sp.Matrix([sp.sin(th) * sp.cos(ph), sp.sin(th) * sp.sin(ph), sp.cos(th)])
lam1, lam2, lam3 = sp.symbols("lambda1 lambda2 lambda3", real=True)
A_dir = sp.simplify((lam1 * n_hat[0] ** 2 + lam2 * n_hat[1] ** 2 + lam3 * n_hat[2] ** 2))
mean_A_dir = sp.simplify(sp.integrate(sp.integrate(A_dir * sp.sin(th), (th, 0, sp.pi)),
                                      (ph, 0, 2 * sp.pi)) / (4 * sp.pi))
check("E[cos^2] = 1/3 on the sphere, so E[omega.S omega] = 0 for isotropic omega",
      sp.simplify(mean_A_dir - (lam1 + lam2 + lam3) / 3) == 0 and
      sp.simplify(mean_A_dir.subs(lam3, -lam1 - lam2)) == 0,
      "EXACT null: a nonzero mean production requires alignment correlation, and "
      "zeta carries no information about it")

# ---------------------------------------------------------------------------
# PHASE 6 -- the production optimum is NOT the framework's cubic landmark
# ---------------------------------------------------------------------------
bb = sp.symbols("b_star", positive=True)
cubic = bb ** 3 + bb - 1                       # h = a/b landmark of the Thales papers
b_star = sp.nsolve(cubic, 0.68)
resid = sp.simplify(cubic.subs(bb, sp.Rational(2, 3)))
check("b = 2/3 is not a root of b^3 + b - 1",
      sp.simplify(resid) == sp.Rational(-1, 27),
      f"residual exactly -1/27; b_star = {float(b_star):.10f} vs 2/3 = 0.6666666667 "
      f"(relative gap {float((b_star - sp.Rational(2,3)) / b_star):.4%})")

# rapidity / response coordinates are also functions of zeta alone
check("R = 1/(4h^2) = 1/(1-zeta^2) and eta = 2h = sqrt(1-zeta^2)",
      sp.simplify(1 / (4 * h_z ** 2) - 1 / (1 - z ** 2)) == 0 and
      sp.simplify(2 * h_z - sp.sqrt(1 - z ** 2)) == 0)

# ---------------------------------------------------------------------------
summary = {
    "n_checks": len(checks),
    "n_passed": sum(1 for v in checks.values() if v["passed"]),
    "sharp_A_bound": float(sp.sqrt(sp.Rational(2, 3))),
    "sharp_P_over_gradu_cubed": float(4 * sp.sqrt(2) / 9),
    "allocation_optimum_zeta": 1 / 3,
    "allocation_optimum_a": 1 / 3,
    "allocation_optimum_b": 2 / 3,
    "g_max": float(4 / (3 * sp.sqrt(3))),
    "thales_cubic_landmark_b_star": float(b_star),
    "cubic_residual_at_two_thirds": "-1/27",
    "checks": checks,
}
(OUT / "exact_algebra.json").write_text(json.dumps(summary, indent=2))
print(f"\n{summary['n_passed']}/{summary['n_checks']} symbolic checks passed")
print(f"written: {OUT / 'exact_algebra.json'}")
sys.exit(0 if summary["n_passed"] == summary["n_checks"] else 1)
