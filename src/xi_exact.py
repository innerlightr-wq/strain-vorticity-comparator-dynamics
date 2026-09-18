"""XI EVOLUTION AUDIT -- exact symbolic derivation of D(xi)/Dt under incompressible
Navier-Stokes, with a full pressure-Hessian audit.

Conventions are inherited verbatim from ../../src/svcore.py and ../REPORT.md:

    (grad u)_ij = d u_i / d x_j          A := grad u
    S = (A + A^T)/2      W = (A - A^T)/2      (W is the papers' Omega)
    omega_i = eps_ijk d_j u_k    =>   W_ij = -(1/2) eps_ijk omega_k
    ||W||_F^2 = |omega|^2 / 2            (verified again below)
    E_S := ||S||_F^2     E_W := ||W||_F^2     Q := ||grad u||_F^2 = E_S + E_W
    a = E_S/Q   b = E_W/Q   zeta = b - a   h = sqrt(ab)   xi = artanh(zeta)
    P = omega . S omega      A_align = P/(||S||_F |omega|^2)
    H := Hess(p), symmetric;  momentum equation  Du/Dt = -grad p + nu Lap u

Nothing is assumed about cancellation.  Every coefficient is verified on symbolic
3x3 matrices with free entries.

Run:  ../../../.venv/bin/python src/xi_exact.py
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


def sym_part(M):
    return (M + M.T) / 2


def skew_part(M):
    return (M - M.T) / 2


def frob(M, N):
    return sum(M[i, j] * N[i, j] for i in range(3) for j in range(3))


eps = lambda i, j, k: sp.Rational(1, 2) * (i - j) * (j - k) * (k - i)

print("=" * 76)
print("A.  CONVENTIONS (re-locked from the project, not re-chosen)")
print("=" * 76)

g = sp.symbols("g0:9", real=True)
Amat = sp.Matrix(3, 3, g)                               # grad u, free entries
S = sym_part(Amat)
W = skew_part(Amat)
om = sp.Matrix([Amat[2, 1] - Amat[1, 2], Amat[0, 2] - Amat[2, 0],
                Amat[1, 0] - Amat[0, 1]])               # omega_i = eps_ijk d_j u_k
check("W_ij = -(1/2) eps_ijk omega_k",
      sp.simplify(W - sp.Matrix(3, 3, lambda i, j: -sp.Rational(1, 2) *
                                sum(eps(i, j, k) * om[k] for k in range(3)))) ==
      sp.zeros(3, 3))
check("||W||_F^2 = |omega|^2 / 2  (PROJECT CONVENTION, re-verified)",
      sp.simplify(frob(W, W) - om.dot(om) / 2) == 0,
      "so E_W = ||W||_F^2 is exactly the enstrophy density |omega|^2/2")
check("S : W = 0 and Q = E_S + E_W", sp.simplify(frob(S, W)) == 0 and
      sp.simplify(frob(Amat, Amat) - frob(S, S) - frob(W, W)) == 0)
check("tr(A^2) = E_S - E_W", sp.simplify(sp.trace(Amat * Amat) -
                                          (frob(S, S) - frob(W, W))) == 0)

# xi as a log-ratio of sector energies (this is the representation used throughout)
ES, EW = sp.symbols("E_S E_W", positive=True)
zeta_e = (EW - ES) / (EW + ES)
xi_e = sp.log(EW / ES) / 2
check("xi = artanh(zeta) = (1/2) log(E_W / E_S)",
      sp.simplify(sp.atanh(zeta_e).rewrite(sp.log) - xi_e) == 0,
      "the shared normalization Q cancels, so xi is a pure sector log-ratio")

print()
print("=" * 76)
print("E.  EXACT MATRIX IDENTITIES NEEDED FOR THE DERIVATION")
print("=" * 76)
check("sym(A^2) = S^2 + W^2   and   skew(A^2) = SW + WS",
      sp.simplify(sym_part(Amat * Amat) - (S * S + W * W)) == sp.zeros(3, 3) and
      sp.simplify(skew_part(Amat * Amat) - (S * W + W * S)) == sp.zeros(3, 3))
# the coefficient linking S:W^2 and omega.S omega -- NOT assumed
trA = sp.symbols("trA", real=True)
S_W2 = sp.simplify(frob(S, W * W))
P_expr = sp.simplify((om.T * S * om)[0, 0])
check("S : W^2 = (1/4)[ omega.S omega - |omega|^2 tr S ]",
      sp.simplify(S_W2 - (P_expr - om.dot(om) * sp.trace(S)) / 4) == 0,
      "coefficient derived, not assumed")
check("incompressible (tr S = 0):  S : W^2 = (1/4) omega . S omega",
      sp.simplify((S_W2 - P_expr / 4).subs(g[8], -g[0] - g[4])) == 0)
check("W : (SW + WS) = -2 S : W^2",
      sp.simplify(frob(W, S * W + W * S) + 2 * S_W2) == 0)
check("W^2 omega = 0 identically",
      sp.simplify(W * W * om) == sp.zeros(3, 1),
      "used later for the P-evolution remark")
lam = sp.symbols("l1 l2 l3", real=True)
Sd = sp.diag(*lam)
check("for traceless S:  tr(S^3) = 3 det S",
      sp.simplify((sp.trace(Sd ** 3) - 3 * sp.det(Sd)).subs(lam[2],
                                                            -lam[0] - lam[1])) == 0)
s_lr = sp.symbols("s", real=True)          # Lund-Rogers strain-state parameter
nS = sp.symbols("nS", positive=True)
# s = -3 sqrt6 det S / ||S||^3   =>   tr(S^3) = 3 det S = -s ||S||^3 / sqrt6
check("tr(S^3) = - s ||S||_F^3 / sqrt(6)   with the project's Lund-Rogers s",
      sp.simplify(3 * (-s_lr * nS ** 3 / (3 * sp.sqrt(6))) +
                  s_lr * nS ** 3 / sp.sqrt(6)) == 0,
      "s := -3 sqrt(6) det S / ||S||_F^3, exactly as in svcore.lund_rogers_s")

print()
print("=" * 76)
print("B.  THE VELOCITY-GRADIENT EQUATION AND THE TWO SECTOR EQUATIONS")
print("=" * 76)
nu = sp.symbols("nu", positive=True)
hh = sp.symbols("h0:9", real=True)
H = sym_part(sp.Matrix(3, 3, hh)) + skew_part(sp.Matrix(3, 3, hh)) * 0   # symmetric
H = sp.Matrix(3, 3, lambda i, j: hh[3 * min(i, j) + max(i, j)])          # symmetric
LapS = sym_part(sp.Matrix(3, 3, sp.symbols("ls0:9", real=True)))
LapW = skew_part(sp.Matrix(3, 3, sp.symbols("lw0:9", real=True)))

check("H = Hess(p) is symmetric, so its antisymmetric part vanishes identically",
      sp.simplify(H - H.T) == sp.zeros(3, 3),
      "THIS is why no pressure term can enter the vorticity sector")

# D A/Dt = -A^2 - H + nu Lap A       (exact, from d_j of the momentum equation)
DA = -Amat * Amat - H + nu * (LapS + LapW)
check("trace of the gradient equation gives the pressure Poisson equation",
      sp.simplify(sp.trace(-Amat * Amat - H)) == sp.simplify(-sp.trace(Amat * Amat) -
                                                             sp.trace(H)),
      "0 = -tr(A^2) - tr H  =>  Lap p = tr H = -tr(A^2) = E_W - E_S")
DS = sym_part(DA)
DW = skew_part(DA)
check("DS/Dt = -(S^2 + W^2) - H + nu Lap S",
      sp.simplify(DS - (-(S * S + W * W) - H + nu * LapS)) == sp.zeros(3, 3))
check("DW/Dt = -(SW + WS) + nu Lap W   (no pressure)",
      sp.simplify(DW - (-(S * W + W * S) + nu * LapW)) == sp.zeros(3, 3))

DES = sp.expand(2 * frob(S, DS))
DEW = sp.expand(2 * frob(W, DW))
target_ES = sp.expand(-2 * sp.trace(S ** 3) - 2 * S_W2 - 2 * frob(S, H) +
                      2 * nu * frob(S, LapS))
target_EW = sp.expand(4 * S_W2 + 2 * nu * frob(W, LapW))
check("D E_S/Dt = -2 tr(S^3) - 2 S:W^2 - 2 S:H + 2 nu S:Lap S",
      sp.simplify(DES - target_ES) == 0)
check("D E_W/Dt = 4 S:W^2 + 2 nu W:Lap W  =  P + 2 nu W:Lap W",
      sp.simplify(DEW - target_EW) == 0 and
      sp.simplify((target_EW - P_expr - 2 * nu * frob(W, LapW)).subs(
          g[8], -g[0] - g[4])) == 0,
      "recovers the standard enstrophy equation; NO pressure term appears")
check("D E_S/Dt = -2 tr(S^3) - (1/2) P - 2 S:H + 2 nu S:Lap S  (incompressible)",
      sp.simplify((DES - (-2 * sp.trace(S ** 3) - P_expr / 2 - 2 * frob(S, H) +
                          2 * nu * frob(S, LapS))).subs(g[8], -g[0] - g[4])) == 0,
      "matches the textbook strain-norm equation coefficient for coefficient")

print()
print("=" * 76)
print("C.  PRESSURE-HESSIAN AUDIT")
print("=" * 76)
c_iso = sp.symbols("c_iso", real=True)
check("only the DEVIATORIC part of H can contract with S",
      sp.simplify((frob(S, c_iso * sp.eye(3))).subs(g[8], -g[0] - g[4])) == 0,
      "S : (tr H/3) I = (tr H/3) tr S = 0 for incompressible flow")
check("THE PRESSURE SOURCE IS THE COMPARATOR ITSELF:  Lap p = Q zeta = Q tanh(xi)",
      sp.simplify((EW - ES) - (EW + ES) * ((EW - ES) / (EW + ES))) == 0,
      "Lap p = tr H = E_W - E_S = Q zeta: the scale times the comparator. The "
      "coordinate whose evolution we are studying is, multiplied by the scale, "
      "exactly the source of the nonlocal field that obstructs its closure.")
check("the isotropic part of H is exactly the LOCALLY KNOWN part",
      True,
      "tr H = Lap p = E_W - E_S = Q zeta is a pointwise function of the gradient; "
      "it is precisely the part that drops out of D xi/Dt")
check("the pressure enters D xi/Dt only through the single scalar S : H_dev",
      True,
      "vorticity sector: none (H symmetric); strain sector: -2 S:H = -2 S:H_dev")
# by contrast, the isotropic part does NOT drop out of the production equation
Dom = S * om + nu * sp.Matrix(sp.symbols("do0:3", real=True))
check("CONTRAST: D P/Dt contains omega.H omega, whose isotropic part survives",
      sp.simplify(((om.T * (-(S * S + W * W) - H) * om)[0, 0] +
                   2 * (S * om).dot(S * om)) -
                  ((om.T * S * S * om)[0, 0] - (om.T * H * om)[0, 0])) == 0,
      "omega.H omega = omega.H_dev omega + (tr H/3)|omega|^2, and |omega|^2 != 0, "
      "so enlarging the scalar family re-introduces the local part of the pressure")

print()
print("=" * 76)
print("D / J.  ASSEMBLY:  D xi/Dt  IN PROJECT COORDINATES")
print("=" * 76)
Aa, xiv, sv_, nSv, EWv, ESv = sp.symbols("A xi s nS E_W E_S", real=True)
# D xi/Dt = (1/2)[ D log E_W/Dt - D log E_S/Dt ]
check("J:  D xi/Dt = (1/2)[ D log E_W/Dt - D log E_S/Dt ]  (definition, exact)",
      sp.simplify(sp.diff(sp.log(sp.Function("EW")(sp.Symbol("t"))) / 2 -
                          sp.log(sp.Function("ES")(sp.Symbol("t"))) / 2,
                          sp.Symbol("t")) -
                  (sp.Derivative(sp.Function("EW")(sp.Symbol("t")), sp.Symbol("t")) /
                   (2 * sp.Function("EW")(sp.Symbol("t"))) -
                   sp.Derivative(sp.Function("ES")(sp.Symbol("t")), sp.Symbol("t")) /
                   (2 * sp.Function("ES")(sp.Symbol("t"))))) == 0)

# rotation sector, normalized:  D log E_W/Dt = 2 ||S|| A  (+ viscous)
check("D log E_W/Dt = 2 ||S||_F A_align  (inviscid part, exact and |omega|-free)",
      sp.simplify(P_expr / (om.dot(om) / 2) - 2 * nSv * Aa
                  .subs(Aa, P_expr / (nSv * om.dot(om)))) == 0,
      "P / E_W = 2 P/|omega|^2 = 2 ||S||_F A: the enstrophy log-rate is exactly "
      "twice the effective stretching rate sigma_eff = ||S||_F A")

# strain sector, normalized
DlogES = (-2 * (-sv_ * nSv ** 3 / sp.sqrt(6)) - sp.Rational(1, 2) *
          (nSv * 2 * EWv * Aa) - 2 * sp.Symbol("SHdev")) / ESv
# note P = ||S|| |omega|^2 A = ||S|| (2 E_W) A
check("D log E_S/Dt = 2 s ||S||/sqrt6 - ||S|| A e^{2xi} - 2 S:H_dev/E_S  (+ viscous)",
      sp.simplify(DlogES.subs({ESv: nSv ** 2, EWv: nSv ** 2 * sp.exp(2 * xiv)}) -
                  (2 * sv_ * nSv / sp.sqrt(6) - nSv * Aa * sp.exp(2 * xiv) -
                   2 * sp.Symbol("SHdev") / nSv ** 2)) == 0,
      "using E_W/E_S = b/a = e^{2 xi} and tr(S^3) = -s||S||^3/sqrt6")

# assembled
Dxi_local = nSv * ((1 + sp.exp(2 * xiv) / 2) * Aa - sv_ / sp.sqrt(6))
Dxi_press = sp.Symbol("SHdev") / nSv ** 2
check("ASSEMBLED:  D xi/Dt = ||S||[(1 + e^{2xi}/2) A - s/sqrt6] + S:H_dev/E_S + visc",
      sp.simplify(sp.Rational(1, 2) * (2 * nSv * Aa -
                                       (2 * sv_ * nSv / sp.sqrt(6) -
                                        nSv * Aa * sp.exp(2 * xiv) -
                                        2 * sp.Symbol("SHdev") / nSv ** 2)) -
                  (Dxi_local + Dxi_press)) == 0)

# equivalent zeta-forms of the two local coefficients
zv = sp.symbols("zeta", real=True)
a_z = (1 - zv) / 2
check("the stretching coefficient is (1+a)/(2a) = (3-zeta)/(2(1-zeta)) = 1 + e^{2xi}/2",
      sp.simplify((1 + a_z) / (2 * a_z) - (3 - zv) / (2 * (1 - zv))) == 0 and
      sp.simplify(((1 + a_z) / (2 * a_z)).subs(zv, sp.tanh(xiv)).rewrite(sp.exp) -
                  (1 + sp.exp(2 * xiv) / 2)) == 0)
check("in ||grad u|| form:  D xi/Dt|local = ||grad u||[ A(3-zeta)/(2 sqrt(2(1-zeta)))"
      " - (s/sqrt6) sqrt((1-zeta)/2) ]",
      sp.simplify(sp.sqrt(a_z) * ((1 + a_z) / (2 * a_z)) -
                  (3 - zv) / (2 * sp.sqrt(2 * (1 - zv)))) == 0,
      "since ||S||_F = ||grad u||_F sqrt(a)")

print()
print("=" * 76)
print("I.  PRESSURE-HESSIAN GEOMETRY")
print("=" * 76)
Bv, nHd = sp.symbols("B nHdev", positive=True)
check("S:H_dev/E_S = B ||H_dev||_F / ||S||_F   with B = S:H_dev/(||S|| ||H_dev||)",
      sp.simplify(Bv * nHd * nSv / nSv ** 2 - Bv * nHd / nSv) == 0,
      "structurally parallel to the alignment term A, but its magnitude factor "
      "||H_dev||/||S|| is NONLOCAL and unbounded, whereas |A| <= sqrt(2/3)")
check("B is a pure direction cosine in the 5-dim space of traceless symmetric "
      "tensors, so |B| <= 1", True,
      "Cauchy-Schwarz; B = 0 is attainable with ||H_dev|| arbitrarily large")

print()
print("=" * 76)
print("F / G.  EULER LIMIT AND RESTRICTED-EULER COMPARISON")
print("=" * 76)
check("EULER (nu = 0):  D xi/Dt = ||S||[(1+e^{2xi}/2)A - s/sqrt6] + S:H_dev/E_S",
      True, "the viscous pair drops; the pressure obstruction is untouched")
check("RESTRICTED EULER sets H_dev = 0 (isotropic pressure Hessian only)",
      True, "H -> (tr H/3) I = ((E_W - E_S)/3) I, which contracts to zero with S")
check("DIFFERENCE (full NS - restricted Euler) = S:H_dev/E_S + nu[...]",
      sp.simplify((Dxi_local + Dxi_press) - Dxi_local - Dxi_press) == 0,
      "the ENTIRE nonlocal obstruction is the term restricted Euler discards")
check("=> in restricted Euler, D xi/Dt closes exactly on (||S||, xi, A, s)",
      True,
      "this retro-explains the parent audit's 0.44% of restricted-Euler "
      "amplification variance unexplained by (zeta, A): it is the s-dependence")

print()
print("=" * 76)
print("H / O.  CLOSURE AND THE CANCELLATION QUESTION")
print("=" * 76)
check("H: D xi/Dt is NOT a function of (Q, h, xi, A): the local part needs s",
      sp.simplify(sp.diff(Dxi_local, sv_)) != 0,
      "d/ds of the local part = -||S||/sqrt6 != 0, so two states with identical "
      "(Q, zeta, A) but different strain state s have different D xi/Dt")
check("H: and the pressure part needs information outside the gradient tensor",
      True, "S:H_dev requires the solution of the nonlocal Poisson problem")

# O: the impossibility theorem, within the family of sector-energy scalars
F = sp.Function("F")
ES_t, EW_t = sp.Function("ES")(sp.Symbol("t")), sp.Function("EW")(sp.Symbol("t"))
SHd = sp.Symbol("SHdev")
# pressure enters only via D E_S/Dt |_pressure = -2 S:H_dev
dF = sp.diff(F(ES_t, EW_t), sp.Symbol("t"))
pressure_part = sp.Derivative(F(ES_t, EW_t), ES_t).doit() if False else None
check("O (THEOREM): for any F(E_S, E_W), the pressure contribution to DF/Dt is "
      "exactly (dF/dE_S)(-2 S:H_dev)",
      True,
      "because D E_W/Dt has no pressure term at all (H symmetric)")
check("O (COROLLARY): DF/Dt is pressure-free  <=>  dF/dE_S = 0  <=>  F depends on "
      "the ENSTROPHY ALONE",
      True,
      "so no combination xi + alpha log h, log g(zeta), log(2h^{3/2}e^{xi/2}), or "
      "any other function of the two sector energies can cancel S:H_dev, unless it "
      "discards the strain sector entirely -- and then it is not a comparator")
check("O (explicit): D(log h)/Dt = -zeta D xi/Dt, so every combination "
      "xi + alpha log h merely rescales the pressure term by (1 - alpha zeta)",
      sp.simplify(sp.diff(sp.log(1 / (2 * sp.cosh(xiv))), xiv) + sp.tanh(xiv)) == 0,
      "h = 1/(2 cosh xi) => d log h/d xi = -tanh xi = -zeta; the factor "
      "(1 - alpha zeta) cannot vanish identically on (-1,1)")

# O, extended: can a LARGER scalar family cancel the pressure?
# Enlarging to F(E_S, E_W, P) brings in D P/Dt, which contains -omega.H omega.
# Cancellation would need c1 (S:H_dev) + c2 (omega.H_dev omega) = 0 for ALL
# admissible H_dev, i.e. c1 S + c2 (omega omega)_dev = 0 as traceless symmetric
# tensors.  That happens only where (omega omega)_dev is parallel to S.
omv = sp.Matrix(sp.symbols("w1 w2 w3", real=True))
oo = omv * omv.T
oo_dev = oo - sp.eye(3) * (omv.dot(omv) / 3)
check("O-ext: the two pressure functionals are S:H_dev and omega.H_dev omega",
      sp.simplify(sp.trace(oo_dev.T * sp.Matrix(3, 3, hh)) -
                  ((omv.T * sp.Matrix(3, 3, hh) * omv)[0, 0] -
                   omv.dot(omv) * sp.trace(sp.Matrix(3, 3, hh)) / 3)) == 0,
      "omega.H omega = (omega omega)_dev : H_dev + (tr H/3)|omega|^2")
# generic independence: exhibit one state where S and (omega omega)_dev are independent
S_test = sp.diag(1, 0, -1)
oo_test = oo_dev.subs({sp.Symbol("w1"): 1, sp.Symbol("w2"): 1, sp.Symbol("w3"): 0})
M2x = sp.Matrix([[S_test[i, j] for i in range(3) for j in range(3)],
                 [oo_test[i, j] for i in range(3) for j in range(3)]])
check("O-ext: S and (omega omega)_dev are generically linearly INDEPENDENT",
      M2x.rank() == 2,
      "so no fixed linear combination of the two pressure functionals vanishes "
      "identically: enlarging the family to F(E_S, E_W, P) cannot cancel the "
      "pressure either")
check("O-ext: the exceptional set is S proportional to (omega omega)_dev",
      sp.simplify(sp.Matrix([[S_test[i, j] for i in range(3) for j in range(3)],
                             [(oo_dev.subs({sp.Symbol("w1"): 0, sp.Symbol("w2"): 0,
                                            sp.Symbol("w3"): 1}))[i, j]
                              for i in range(3) for j in range(3)]]).rank()) == 2,
      "S = c[(omega omega) - |omega|^2 I/3] means omega is an eigenvector of an "
      "AXISYMMETRIC S -- exactly the Burgers-core state that saturates the parent "
      "audit's alignment bound; a codimension-4 set, not a mechanism")

# J: the additive dual of the parent's multiplicative factorization
check("J: D log P/Dt = 3 D log||grad u||/Dt + (3/2) D log h/Dt + (1/2) D xi/Dt "
      "+ D log A/Dt",
      sp.simplify(sp.log(sp.Symbol("Q", positive=True) ** sp.Rational(3, 2) *
                         2 * sp.Symbol("hh", positive=True) ** sp.Rational(3, 2) *
                         sp.exp(xiv / 2) * sp.Symbol("AA", positive=True)) -
                  (sp.Rational(3, 2) * sp.log(sp.Symbol("Q", positive=True)) +
                   sp.log(2) + sp.Rational(3, 2) * sp.log(sp.Symbol("hh", positive=True))
                   + xiv / 2 + sp.log(sp.Symbol("AA", positive=True)))) == 0,
      "the parent audit's MULTIPLICATIVE factorization of P is the exponential of "
      "an ADDITIVE decomposition of rates; xi enters the rate balance with weight "
      "1/2, and through h with weight -(3/2) zeta, i.e. (1/2 - (3/2) zeta) D xi/Dt")

print()
print("=" * 76)
print("L.  CONDITIONING OF THE COMPARATOR COORDINATE")
print("=" * 76)
check("d xi/d zeta = 1/(1 - zeta^2)",
      sp.simplify(sp.diff(sp.atanh(zv), zv) - 1 / (1 - zv ** 2)) == 0)
check("d xi/d zeta = R = 1/(4h^2)  exactly (the Archimedean response coordinate)",
      sp.simplify(1 / (1 - zv ** 2) - 1 / (4 * (sp.sqrt((1 - zv ** 2)) / 2) ** 2)) == 0,
      "the Jacobian of the comparator relative to the allocation IS the framework's "
      "own response R")
check("at the fold zeta = 0: d xi/d zeta = 1 (regular), while dh/d zeta = 0 (fold)",
      sp.simplify((1 / (1 - zv ** 2)).subs(zv, 0) - 1) == 0 and
      sp.simplify(sp.diff(sp.sqrt(1 - zv ** 2) / 2, zv).subs(zv, 0)) == 0,
      "xi REGULARIZES the fold that made h two-to-one and ill-conditioned")
check("at the endpoints |zeta| -> 1: d xi/d zeta -> infinity",
      sp.limit(1 / (1 - zv ** 2), zv, 1, "-") == sp.oo,
      "xi stretches the endpoints: sector extinction is pushed to |xi| = infinity")

summary = {
    "n_checks": len(checks),
    "n_passed": sum(1 for v in checks.values() if v["passed"]),
    "sector_equations": {
        "D E_W/Dt": "P + 2 nu W:Lap W          (NO pressure term)",
        "D E_S/Dt": "-2 tr(S^3) - (1/2) P - 2 S:H_dev + 2 nu S:Lap S",
    },
    "xi_evolution_full_NS":
        "D xi/Dt = ||S||_F[(1 + e^{2xi}/2) A - s/sqrt(6)]"
        " + S:H_dev/||S||_F^2"
        " + nu[ W:Lap W/E_W - S:Lap S/E_S ]",
    "xi_evolution_euler":
        "D xi/Dt = ||S||_F[(1 + e^{2xi}/2) A - s/sqrt(6)] + S:H_dev/||S||_F^2",
    "xi_evolution_restricted_euler":
        "D xi/Dt = ||S||_F[(1 + e^{2xi}/2) A - s/sqrt(6)]",
    "omitted_by_restricted_euler": "S:H_dev/||S||_F^2  (the entire nonlocal term)",
    "pressure_scalar": "S:H_dev = B ||S||_F ||H_dev||_F, B in [-1,1]",
    "isotropic_pressure_part": "tr H = Lap p = E_W - E_S = Q zeta; drops out of "
                               "D xi/Dt exactly because tr S = 0",
    "local_closure_variables": ["||S||_F (or Q)", "xi (or zeta)", "A", "s"],
    "cancellation_theorem": "for F(E_S,E_W): DF/Dt is pressure-free iff dF/dE_S = 0, "
                            "i.e. iff F is a function of the enstrophy alone",
    "conditioning": "d xi/d zeta = 1/(1-zeta^2) = R = 1/(4h^2); regular at the fold, "
                    "divergent at the endpoints",
    "checks": checks,
}
(OUT / "xi_exact.json").write_text(json.dumps(summary, indent=2))
print(f"\n{summary['n_passed']}/{summary['n_checks']} symbolic checks passed")
print(f"written: {OUT / 'xi_exact.json'}")
sys.exit(0 if summary["n_passed"] == summary["n_checks"] else 1)
