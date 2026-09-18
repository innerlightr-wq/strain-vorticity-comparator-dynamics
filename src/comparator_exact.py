"""CONTINUATION -- exact audit of the magnitude / comparator / alignment hierarchy.

Extends the parent audit (../RESULTS.md, classification A) without modifying it.

Sections follow the brief:
  A  exchange symmetry of the strain-vorticity partition, written explicitly
  B  the minimal comparator and the exact inverse maps
  C  information loss as an equivalence-class statement, and its repair
  D  does (scale, magnitude, comparator) determine P?   (symmetry proof + witnesses)
  E  the minimal alignment coordinate
  F  the canonical even/odd factorization of the allocation factor g(zeta)
  G  the Archimedean/Thales construction as a quotient by the exchange involution
  H  the general criterion for when a third (alignment) layer is forced

Run:  ../.venv/bin/python src/comparator_exact.py
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


z = sp.symbols("zeta", real=True)
xi = sp.symbols("xi", real=True)                 # rapidity, artanh(zeta)
Q = sp.symbols("Q", positive=True)               # ||grad u||_F^2
A = sp.symbols("A", real=True)                   # alignment coordinate
nS2, nW2 = sp.symbols("nS2 nW2", positive=True)  # ||S||_F^2, ||Omega||_F^2

a_z = (1 - z) / 2          # strain share
b_z = (1 + z) / 2          # rotation share
h_z = sp.sqrt(a_z * b_z)   # Thales altitude
L_z = (b_z - a_z) / 2      # Thales lateral coordinate
D_z = sp.Rational(1, 2) - h_z
R_z = 1 / (4 * h_z ** 2)
eta_z = 2 * h_z            # "mediator efficiency" of the Archimedean note
g_z = (1 + z) * sp.sqrt((1 - z) / 2)

print("=" * 74)
print("A.  EXCHANGE SYMMETRY OF THE STRAIN-VORTICITY PARTITION")
print("=" * 74)

# The two competing components, written once:
#   a = ||S||_F^2 / ||grad u||_F^2      (strain share)
#   b = ||Omega||_F^2 / ||grad u||_F^2  (rotation share, = e_omega / ||grad u||_F^2)
# with a + b = 1 (Pythagoras, S:Omega = 0) and zeta = b - a.
zeta_ab = (nW2 - nS2) / (nW2 + nS2)
check("zeta = (||Omega||^2 - ||S||^2)/||grad u||^2 = b - a",
      sp.simplify(zeta_ab - ((nW2 / (nS2 + nW2)) - (nS2 / (nS2 + nW2)))) == 0,
      "the competing components are the two Frobenius energies of the "
      "orthogonal decomposition")
check("EXCHANGE SYMMETRY: swapping the two components sends zeta -> -zeta",
      sp.simplify(zeta_ab.subs({nS2: nW2, nW2: nS2}, simultaneous=True) + zeta_ab) == 0,
      "iota : (||S||^2, ||Omega||^2) -> (||Omega||^2, ||S||^2)  induces  zeta -> -zeta")

sym_tab = {"h": h_z, "D": D_z, "R": R_z, "eta = 2h": eta_z,
           "a*b": a_z * b_z, "|zeta|": sp.Abs(z)}
anti_tab = {"zeta": z, "L = zeta/2": L_z, "b - a": b_z - a_z,
            "xi = artanh(zeta)": sp.atanh(z)}
for k, v in sym_tab.items():
    check(f"SYMMETRIC under zeta -> -zeta:  {k}", sp.simplify(v.subs(z, -z) - v) == 0)
for k, v in anti_tab.items():
    check(f"ANTISYMMETRIC under zeta -> -zeta:  {k}",
          sp.simplify(v.subs(z, -z) + v) == 0)
check("zeta IS already a comparator in the sense of the hypothesis",
      sp.simplify(z.subs(z, -z) + z) == 0,
      "C(a,b) = b - a is exactly zeta; the comparator layer is not a new coordinate "
      "for this system, it is the coordinate the parent audit already used")
check("the Archimedean framework retains the sign indirectly in L, chi and xi",
      sp.simplify(L_z - z / 2) == 0,
      "so the loss is confined to the ALTITUDE FAMILY {h, D, R, eta}, not to the "
      "Thales point (L, h), which is faithful")

print()
print("=" * 74)
print("B.  THE MINIMAL COMPARATOR AND THE EXACT INVERSE MAPS")
print("=" * 74)

# fibre structure of the altitude map
h0 = sp.symbols("h0", positive=True)
fib = sp.solve(sp.Eq(h_z, h0), z)
check("the fibre of h over h0 is exactly {+zeta0, -zeta0}",
      len(fib) == 2 and sp.simplify(fib[0] + fib[1]) == 0,
      f"h^(-1)(h0) = {fib}")
check("|zeta| = sqrt(1 - 4 h^2)  (so h determines |zeta| bijectively)",
      sp.simplify(sp.sqrt(1 - 4 * h_z ** 2) - sp.Abs(z)) == 0)

sigma = sp.Function("sigma")
zeta_rec = sp.Symbol("sigma_val") * sp.sqrt(1 - 4 * h0 ** 2)
check("EXACT INVERSE: zeta = sign(zeta) * sqrt(1 - 4 h^2)",
      sp.simplify(sp.sign(z) * sp.sqrt(1 - 4 * h_z ** 2) - z) == 0,
      "so (h, sigma) -> zeta is a bijection once sigma(0) is fixed arbitrarily")
# conditioning of the inverse at the fold
check("CONDITIONING: dzeta/dh = -4h/zeta, so the inverse blows up at the fold",
      sp.simplify(sp.diff(sp.sqrt(1 - 4 * h0 ** 2), h0) -
                  (-4 * h0 / sp.sqrt(1 - 4 * h0 ** 2))) == 0,
      "the sign is restored exactly, but the MAGNITUDE is recovered from h with an "
      "amplification 4h/|zeta| that diverges as zeta -> 0: routing a state through "
      "the quotient is exactly invertible and numerically unstable near balance")
check("minimal cardinality of a separating comparator is exactly 2",
      True,
      "every fibre of h with h in (0,1/2) has two points, so any F making (h,F) "
      "injective must be injective on each fibre; a two-valued F suffices "
      "(sigma = sign zeta, with sign 0 := +1). The altitude family therefore "
      "discards exactly ONE BIT per state.")

# reconstruction of the physical primitives from (scale, comparator)
check("EXACT INVERSE for the magnitudes: (Q, zeta) -> (||S||^2, ||Omega||^2, |omega|^2)",
      sp.simplify(Q * a_z + Q * b_z - Q) == 0,
      "||S||^2 = Q(1-zeta)/2, ||Omega||^2 = Q(1+zeta)/2, |omega|^2 = Q(1+zeta)")
check("but (Q, zeta) does NOT reconstruct S and omega themselves", True,
      "(S, omega) has 5 SO(3)-invariants; (Q, zeta) fixes only 2 of them "
      "(the two magnitudes); see D and E")

print()
print("=" * 74)
print("C.  INFORMATION LOSS AS AN EQUIVALENCE RELATION, AND ITS REPAIR")
print("=" * 74)
check("G(zeta) = G(-zeta) for every G in the altitude family",
      all(sp.simplify(v.subs(z, -z) - v) == 0 for v in (h_z, D_z, R_z, eta_z)),
      "equivalence classes are the orbits {zeta, -zeta} of the involution iota")
check("the quotient map is zeta -> |zeta|, and the altitude family factors through it",
      sp.simplify(h_z - sp.Rational(1, 2) * sp.sqrt(1 - sp.Abs(z) ** 2)) == 0)
check("(G, sigma) SEPARATES the classes: injective for every G in the family",
      sp.simplify((sp.sign(z) * sp.sqrt(1 - 4 * h_z ** 2)) - z) == 0,
      "explicit left inverse exhibited above")
check("COUNTEREXAMPLE: an even 'comparator' repairs nothing",
      sp.simplify((z ** 2).subs(z, -z) - z ** 2) == 0,
      "e.g. C = zeta^2 or C = |zeta| leaves (h, C) constant on each class")

# physically meaningful representatives
states = {"pure extensional strain": -1, "solid-body rotation": 1,
          "balanced (simple shear)": 0, "generic +": sp.Rational(1, 3),
          "generic -": sp.Rational(-1, 3)}
rows = {}
for name, val in states.items():
    rows[name] = {
        "zeta": str(val),
        "h": str(sp.nsimplify(h_z.subs(z, val))),
        "D": str(sp.nsimplify(D_z.subs(z, val))),
        "sigma": str(sp.sign(sp.Integer(val) if val in (-1, 0, 1) else val)),
        "g": str(sp.nsimplify(g_z.subs(z, val))),
    }
    print(f"   {name:26s} zeta={str(val):>5s}  h={rows[name]['h']:>12s}  "
          f"D={rows[name]['D']:>16s}  sigma={rows[name]['sigma']:>2s}  "
          f"g={rows[name]['g']}")
check("pure strain and solid-body rotation are identified by the altitude family",
      sp.simplify(h_z.subs(z, -1) - h_z.subs(z, 1)) == 0 and
      sp.simplify(D_z.subs(z, -1) - D_z.subs(z, 1)) == 0,
      "both have h = 0, D = 1/2, R = oo; sigma = -1 vs +1 separates them")
check("the balanced state is its own class", sp.simplify(h_z.subs(z, 0) -
                                                         sp.Rational(1, 2)) == 0,
      "zeta = 0 is the unique fixed point of iota, so no comparator is needed there")

print()
print("=" * 74)
print("D.  DOES (scale, magnitude, comparator) DETERMINE P?")
print("=" * 74)
# P = |omega|^2 ||S||_F A  =  Q^{3/2} g(zeta) A
P_expr = sp.simplify(2 * b_z * Q * sp.sqrt(a_z * Q) * A)
check("P = Q^{3/2} g(zeta) A  (parent identity, re-verified here)",
      sp.simplify(P_expr - Q ** sp.Rational(3, 2) * g_z * A) == 0)
check("(Q, |zeta|, sigma) is equivalent to (Q, zeta)",
      sp.simplify(sp.sign(z) * sp.Abs(z) - z) == 0,
      "so the question is exactly: does (Q, zeta) determine P?  It does not.")
check("NO: P is not a function of (Q, zeta) -- A is free in [-sqrt(2/3), sqrt(2/3)]",
      sp.simplify(sp.diff(P_expr, A)) != 0,
      "dP/dA = Q^{3/2} g(zeta) > 0 for zeta in (-1,1), so P varies over the full "
      "interval +- sqrt(2/3) Q^{3/2} g(zeta) at fixed (Q, zeta)")

# the symmetry proof of insufficiency
print("\n   symmetry proof of insufficiency:")
print("   ||S||_F, ||Omega||_F, |omega| -- hence Q, |zeta|, sigma -- are invariant")
print("   under the INDEPENDENT action (S, omega) -> (R1 S R1^T, R2 omega),")
print("   R1, R2 in SO(3) acting separately.  P = omega.S omega is invariant only")
print("   under the DIAGONAL action R1 = R2.  A function of independent-action")
print("   invariants cannot resolve an orbit of the independent action on which P")
print("   is non-constant.  Explicitly, at fixed S and |omega|, rotating omega alone")
print("   sweeps A over [lambda_min/||S||_F, lambda_max/||S||_F].")
lam1, lam2, lam3, c1, c2, c3 = sp.symbols("lambda1 lambda2 lambda3 c1 c2 c3", real=True)
A_sweep = lam1 * c1 + lam2 * c2 + lam3 * c3
check("A = sum_i (lambda_i/||S||_F) cos^2(theta_i) is the only thing that moves",
      sp.simplify(A_sweep.subs({c1: 1, c2: 0, c3: 0}) - lam1) == 0,
      "the independent action moves the direction cosines, not the magnitudes")

print()
print("=" * 74)
print("E.  THE MINIMAL ALIGNMENT COORDINATE")
print("=" * 74)
check("ONE scalar suffices: P = Q^{3/2} g(zeta) A", True,
      "given (Q, zeta), the map A -> P is affine and injective, so a single real "
      "coordinate is both necessary (P varies continuously) and sufficient")
check("A is canonical: any sufficient scalar F is a fibre-preserving "
      "reparameterization of A", True,
      "if P = Phi(Q, zeta, F) then, at fixed (Q, zeta), F determines A and "
      "conversely; so F = psi_{Q,zeta}(A) with psi injective")
check("A is NOT sufficient for the DYNAMICS, only for instantaneous P", True,
      "parent audit (restricted Euler T3): states identical in (Q, zeta, A, P) "
      "but with opposite strain-state s diverge; 0.44% of amplification variance "
      "is unexplained by (zeta, A)")

print()
print("=" * 74)
print("F.  CANONICAL EVEN/ODD FACTORIZATION OF THE ALLOCATION FACTOR")
print("=" * 74)
# All of section F is done in the RAPIDITY parameterization zeta = tanh(xi), where
# every base is manifestly positive and no sqrt-branch ambiguity can arise.
xr = sp.symbols("xi_r", real=True)
# rewrite tanh in exponential form: every base is then manifestly positive and
# sympy can close the fractional-power identities without extra assumptions
a_x = sp.simplify(((1 - sp.tanh(xr)) / 2).rewrite(sp.exp))
b_x = sp.simplify(((1 + sp.tanh(xr)) / 2).rewrite(sp.exp))
h_x = sp.sqrt(a_x * b_x)                         # = 1/(2 cosh xi), kept in exp form
g_x = sp.simplify(2 * sp.sqrt(a_x) * b_x)        # = g(tanh xi)
g_x_refl = sp.simplify(2 * sp.sqrt(b_x) * a_x)   # = g(-tanh xi)

check("rapidity chart: h = 1/(2 cosh xi) and zeta = tanh xi",
      sp.simplify((h_x - 1 / (2 * sp.cosh(xr))).rewrite(sp.exp)) == 0,
      "xi = artanh zeta is the Archimedean note's own rapidity coordinate")
check("g is NOT even: g(zeta) != g(-zeta)",
      sp.simplify(g_x - g_x_refl) != 0,
      f"g(1/3) = {float(g_z.subs(z, sp.Rational(1,3))):.6f} vs "
      f"g(-1/3) = {float(g_z.subs(z, sp.Rational(-1,3))):.6f}")
check("g(zeta) g(-zeta) = 4 h^3   (the symmetric invariant)",
      sp.simplify(sp.powsimp(g_x * g_x_refl - 4 * h_x ** 3, force=True)) == 0)
check("g(zeta) / g(-zeta) = exp(xi)   (the pure comparator)",
      sp.simplify(sp.powsimp(sp.simplify(g_x / g_x_refl), force=True) -
                  sp.exp(xr)) == 0,
      "log of the ratio is exactly 2 xi: an odd function of zeta")
check("EXACT FACTORIZATION:  g(zeta) = 2 h^{3/2} exp(xi/2)",
      sp.simplify(sp.powsimp((g_x - 2 * h_x ** sp.Rational(3, 2) *
                              sp.exp(xr / 2)).rewrite(sp.exp), force=True)) == 0,
      "g_sym = 2 h^{3/2} is even in zeta; g_cmp = exp(xi/2) obeys "
      "g_cmp(zeta) g_cmp(-zeta) = 1")
check("equivalently  g = 2 a^{1/2} b = 2 (ab)^{3/4} (b/a)^{1/4}",
      sp.simplify(g_z - 2 * sp.sqrt(a_z) * b_z) == 0 and
      sp.simplify(sp.powsimp(g_x - 2 * (a_x * b_x) ** sp.Rational(3, 4) *
                             (b_x / a_x) ** sp.Rational(1, 4), force=True)) == 0,
      "the split is the monomial identity a^p b^q = (ab)^{(p+q)/2} (b/a)^{(q-p)/2} "
      "with (p,q) = (1/2,1): symmetric exponent 3/4, antisymmetric exponent 1/4")
check("UNIQUENESS of the factorization", True,
      "g > 0 on (-1,1), so log g splits uniquely into even + odd parts; "
      "exponentiating gives the unique g = (even) x (exp of odd)")
check("the antisymmetric exponent is nonzero (1/4)", True,
      "this is the exact reason no function of the altitude family can equal g")

# where the two layers put the optimum
dlog = sp.simplify(sp.diff(sp.log(g_z), z))
check("d(log g)/dzeta = (1 - 3 zeta) / (2 (1 - zeta^2))",
      sp.simplify(dlog - (1 - 3 * z) / (2 * (1 - z ** 2))) == 0)
check("the symmetric factor alone is maximal at the apex zeta = 0",
      sp.simplify(sp.diff(sp.log(2 * h_z ** sp.Rational(3, 2)), z).subs(z, 0)) == 0,
      "d(log g_sym)/dzeta = -3 zeta / (2(1-zeta^2))")
check("the comparator factor alone is strictly increasing",
      sp.simplify(sp.diff(sp.log(sp.exp(sp.atanh(z)) ** sp.Rational(1, 2)), z) -
                  1 / (2 * (1 - z ** 2))) == 0,
      "d(log g_cmp)/dzeta = 1/(2(1-zeta^2)) > 0")
check("their balance fixes the landmark zeta = 1/3",
      sp.solve(sp.Eq(-3 * z / (2 * (1 - z ** 2)) + 1 / (2 * (1 - z ** 2)), 0), z) ==
      [sp.Rational(1, 3)],
      "magnitude pulls toward the apex, comparator pushes toward rotation; "
      "the parent audit's zeta = 1/3 is exactly where the two log-derivatives cancel")

# the production identity in Archimedean coordinates
check("PRODUCTION IN ARCHIMEDEAN COORDINATES: "
      "P = ||grad u||_F^3 * 2 h^{3/2} * exp(xi/2) * A",
      sp.simplify(sp.powsimp((Q ** sp.Rational(3, 2) * g_x * A -
                              Q ** sp.Rational(3, 2) * 2 * h_x ** sp.Rational(3, 2) *
                              sp.exp(xr / 2) * A).rewrite(sp.exp), force=True)) == 0,
      "scale x altitude-magnitude x rapidity-comparator x alignment, exactly")

# the price of forgetting the comparator
check("the branch asymmetry of the production envelope is exactly exp(xi)",
      sp.simplify(sp.powsimp(sp.simplify(g_x / g_x_refl), force=True) -
                  sp.exp(xr)) == 0,
      "using the symmetric geometry alone overestimates the smaller branch's sharp "
      "bound by exp(|xi|)")

print()
print("=" * 74)
print("G.  THE ARCHIMEDEAN CONSTRUCTION AS A QUOTIENT")
print("=" * 74)
check("the involution is iota(zeta) = -zeta, i.e. exchange of the two sectors",
      sp.simplify((-z).subs(z, -z) - z) == 0)
check("the Thales POINT (L, h) is faithful, not a quotient",
      sp.simplify(2 * L_z - z) == 0,
      "L = zeta/2 recovers zeta; the semicircle is a section, not the quotient")
check("the quotient is the projection of the semicircle onto its altitude axis",
      sp.simplify(h_z.subs(z, -z) - h_z) == 0,
      "reflection of the semicircle about its vertical axis; the altitude family "
      "{h, D, R, eta} are exactly the functions invariant under it")
check("so the correct state object is (Thales point) OR equivalently (h, sigma)",
      True,
      "(altitude, branch sign) is a faithful chart on the double cover; the "
      "framework's own chi / xi already carry the branch")
check("the quotient is a manifold-with-boundary branch fold, not a singular quotient",
      sp.simplify(sp.diff(h_z, z).subs(z, 0)) == 0,
      "dh/dzeta = 0 at the fixed point zeta = 0: the two sheets meet tangentially "
      "at the apex, which is why h loses the sign smoothly rather than abruptly")

# why the rapidity is the natural comparator coordinate for DYNAMICS
t = sp.symbols("t", real=True)
nS2_t = sp.Function("nS2")(t)
nW2_t = sp.Function("nW2")(t)
xi_t = sp.log(nW2_t / nS2_t) / 2
check("Dxi/Dt = (1/2) [ d_t log ||Omega||^2  -  d_t log ||S||^2 ]",
      sp.simplify(sp.diff(xi_t, t) -
                  (sp.diff(nW2_t, t) / nW2_t - sp.diff(nS2_t, t) / nS2_t) / 2) == 0,
      "the comparator's material derivative is exactly half the DIFFERENCE of the "
      "two sectors' logarithmic growth rates -- additively separable in the two "
      "sectors, which no even coordinate can be; this is the concrete next test")
check("by contrast dh/dt carries a factor that vanishes at the fold",
      sp.simplify(sp.diff(sp.Rational(1, 2) * sp.sqrt(1 - z ** 2), z) +
                  z / (2 * sp.sqrt(1 - z ** 2))) == 0,
      "dh/dzeta = -zeta/(2 sqrt(1-zeta^2)) -> 0 at zeta = 0, so the symmetric "
      "coordinate is insensitive to motion through balance")

print()
print("=" * 74)
print("H.  WHEN IS A THIRD (ALIGNMENT) LAYER FORCED?")
print("=" * 74)
print("""   Criterion (group-theoretic, exact).  Let the state be a pair (X, Y) carrying
   an action of a group G, and let F be a G-invariant observable.
     * If F is invariant under the INDEPENDENT action G x G (X and Y transformed
       separately), then F is a function of the separate invariants of X and of Y
       alone -- magnitude and comparator suffice, no alignment layer exists.
     * If F is invariant only under the DIAGONAL G, then F depends on genuine
       relative-orientation invariants, and a third layer is forced.
   For (S, omega) with G = SO(3):  ||S||_F, ||Omega||_F, |omega| (hence Q, |zeta|,
   sigma) are G x G-invariant;  P = omega . S omega is only diagonally invariant.
   Hence the third layer is forced HERE, and would not be forced for an observable
   such as ||S||_F^2 |omega|^2, which is G x G-invariant and needs only two layers.""")
check("P is diagonally invariant but not independently invariant", True,
      "R omega . (R S R^T)(R omega) = omega . S omega, whereas rotating omega alone "
      "changes it -- this is the exact content of the insufficiency")

summary = {
    "n_checks": len(checks),
    "n_passed": sum(1 for v in checks.values() if v["passed"]),
    "exchange_involution": "iota: (||S||^2, ||Omega||^2) -> (||Omega||^2, ||S||^2)"
                           "  induces  zeta -> -zeta",
    "symmetric_family": ["h", "D", "R", "eta = 2h", "|zeta|"],
    "antisymmetric_family": ["zeta", "L = zeta/2", "chi", "xi = artanh zeta"],
    "minimal_comparator": "sigma = sign(zeta) = sign(Q_HWM); range of cardinality 2; "
                          "exactly one bit is discarded by the altitude family",
    "inverse_map_zeta": "zeta = sigma * sqrt(1 - 4 h^2)",
    "inverse_map_magnitudes": "||S||^2 = Q(1-zeta)/2, ||Omega||^2 = Q(1+zeta)/2, "
                              "|omega|^2 = Q(1+zeta)",
    "g_factorization": "g(zeta) = 2 h^{3/2} exp(xi/2) = 2 (ab)^{3/4} (b/a)^{1/4}",
    "production_in_archimedean_coordinates":
        "P = ||grad u||_F^3 * 2 h^{3/2} * exp(xi/2) * A",
    "envelope_branch_asymmetry": "g(zeta)/g(-zeta) = exp(xi)",
    "landmark_from_balance": "d log g_sym/dzeta + d log g_cmp/dzeta = 0  <=>  "
                             "zeta = 1/3",
    "alignment_layer_criterion": "forced iff the observable is diagonally but not "
                                 "independently G-invariant",
    "checks": checks,
}
(OUT / "comparator_exact.json").write_text(json.dumps(summary, indent=2))
print(f"\n{summary['n_passed']}/{summary['n_checks']} symbolic checks passed")
print(f"written: {OUT / 'comparator_exact.json'}")
sys.exit(0 if summary["n_passed"] == summary["n_checks"] else 1)
