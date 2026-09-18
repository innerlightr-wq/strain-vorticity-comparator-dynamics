"""PHASE 4 -- canonical flows: zeta, A, P and the Thales coordinates.

Every velocity field is written in CARTESIAN components and differentiated
symbolically, so no cylindrical-frame convention can leak into the answer.
Undefined quantities are reported as undefined, not regularized.

Flows: solid-body rotation, pure extensional strain, simple shear,
Lamb-Oseen vortex, Burgers vortex (with a sweep over the axial strain a).

Run:  .venv/bin/python src/canonical_flows.py
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

x, y, z, t = sp.symbols("x y z t", real=True)
Om, gam, aa, Gam, nu = sp.symbols("Omega gamma a Gamma nu", positive=True)
XS = (x, y, z)

FAILURES: list[str] = []


def require(name, ok, detail=""):
    print(f"[{'PASS' if ok else 'FAIL':4s}] {name}" + (f"   {detail}" if detail else ""))
    if not ok:
        FAILURES.append(name)


def analyse(u, subs=None, simplify=True):
    """Exact audit coordinates of a symbolic Cartesian velocity field."""
    G = sp.Matrix(3, 3, lambda i, j: sp.diff(u[i], XS[j]))
    if subs:
        G = G.subs(subs)
    S = (G + G.T) / 2
    W = (G - G.T) / 2
    om = sp.Matrix([G[2, 1] - G[1, 2], G[0, 2] - G[2, 0], G[1, 0] - G[0, 1]])
    f2 = lambda M: sum(M[i, j] ** 2 for i in range(3) for j in range(3))
    nS2, nW2 = f2(S), f2(W)
    w2 = om.dot(om)
    P = (om.T * S * om)[0, 0]
    if simplify:
        S, W, om = sp.simplify(S), sp.simplify(W), sp.simplify(om)
        nS2, nW2, w2, P = (sp.simplify(e) for e in (nS2, nW2, w2, P))
    Q = sp.simplify(nS2 + nW2)
    div = sp.simplify(G.trace())
    out = dict(G=G, S=S, W=W, omega=om, div=div, nS2=nS2, nW2=nW2, w2=w2, Qtot=Q, P=P)
    out["zeta"] = sp.simplify((nW2 - nS2) / Q) if Q != 0 else sp.nan
    out["A"] = sp.simplify(P / (sp.sqrt(nS2) * w2)) if (nS2 != 0 and w2 != 0) else None
    if Q != 0:
        a_ = sp.simplify(nS2 / Q)
        b_ = sp.simplify(nW2 / Q)
        out.update(a=a_, b=b_, h=sp.simplify(sp.sqrt(a_ * b_)),
                   L=sp.simplify((b_ - a_) / 2),
                   D=sp.simplify(sp.Rational(1, 2) - sp.sqrt(a_ * b_)))
    return out


def show(tag, r, notes=""):
    print(f"\n--- {tag} ---")
    print(f"  div u = {r['div']}")
    print(f"  ||S||^2 = {r['nS2']}   ||W||^2 = {r['nW2']}   |omega|^2 = {r['w2']}")
    print(f"  zeta = {r['zeta']}   a = {r.get('a')}   b = {r.get('b')}")
    print(f"  h = {r.get('h')}   L = {r.get('L')}   D = {r.get('D')}")
    print(f"  P = {r['P']}   A = {r['A'] if r['A'] is not None else 'UNDEFINED'}")
    if notes:
        print(f"  note: {notes}")


rows = {}


def record(tag, r, defined_A, sign_P, notes):
    rows[tag] = {
        "zeta": str(r["zeta"]), "a": str(r.get("a")), "b": str(r.get("b")),
        "h": str(r.get("h")), "L": str(r.get("L")), "D": str(r.get("D")),
        "P": str(r["P"]), "A": (str(r["A"]) if r["A"] is not None else "UNDEFINED"),
        "A_defined": defined_A, "sign_P": sign_P, "notes": notes,
    }


# ---------------------------------------------------------------------------
print("=" * 72)
print("PHASE 4: canonical flows")
print("=" * 72)

# 1. solid-body rotation
u1 = [-Om * y, Om * x, sp.Integer(0)]
r1 = analyse(u1)
show("solid-body rotation  u = (-Omega y, Omega x, 0)", r1,
     "S = 0, so A is UNDEFINED (0/0); zeta = +1 saturates the rotation end")
require("solid rotation: S = 0, zeta = +1, P = 0",
        r1["nS2"] == 0 and r1["zeta"] == 1 and r1["P"] == 0)
require("solid rotation: A undefined (||S||_F = 0)", r1["A"] is None)
record("solid-body rotation", r1, False, "0",
       "A undefined: no strain to align with; h = 0 and D = 1/2 (Thales 'fully degenerate')")

# 2. pure extensional strain
u2 = [aa * x, -aa * y, sp.Integer(0)]
r2 = analyse(u2)
show("pure extensional strain  u = (a x, -a y, 0)", r2,
     "omega = 0, so A is UNDEFINED (0/0); zeta = -1 saturates the strain end")
require("pure strain: omega = 0, zeta = -1, P = 0",
        r2["w2"] == 0 and r2["zeta"] == -1 and r2["P"] == 0)
require("pure strain: A undefined (|omega| = 0)", r2["A"] is None)
record("pure extensional strain", r2, False, "0",
       "A undefined: no vorticity to stretch; h = 0 and D = 1/2, IDENTICAL to solid "
       "rotation although the states are physically opposite")

require("FALSIFIER: solid rotation and pure strain share h = 0, D = 1/2, R = infinity",
        r1["h"] == r2["h"] == 0 and r1["D"] == r2["D"] == sp.Rational(1, 2),
        "the Thales altitude/deficit cannot distinguish pure rotation from pure strain; "
        "zeta separates them exactly (+1 vs -1)")

# 3. simple shear
u3 = [gam * y, sp.Integer(0), sp.Integer(0)]
r3 = analyse(u3)
show("simple shear  u = (gamma y, 0, 0)", r3,
     "the Thales apex zeta = 0, h = 1/2, D = 0 -- yet P = 0 and A = 0 exactly")
require("simple shear: zeta = 0 (apex), h = 1/2, D = 0, P = 0, A = 0",
        r3["zeta"] == 0 and r3["h"] == sp.Rational(1, 2) and r3["D"] == 0 and
        r3["P"] == 0 and sp.simplify(r3["A"]) == 0)
record("simple shear", r3, True, "0",
       "sits exactly at the Thales apex (maximal 'coherence' h = 1/2) with zero "
       "enstrophy production: the apex is not a production landmark")

# 4. Lamb-Oseen vortex (unstretched), Cartesian
rr = sp.sqrt(x ** 2 + y ** 2)
u_theta_LO = Gam / (2 * sp.pi * rr) * (1 - sp.exp(-rr ** 2 / (4 * nu * t)))
u4 = [-u_theta_LO * y / rr, u_theta_LO * x / rr, sp.Integer(0)]
# evaluate on the x-axis (WLOG by axisymmetry): y -> 0, x -> r > 0
rsym = sp.symbols("r", positive=True)
r4 = analyse(u4, subs={y: 0, x: rsym}, simplify=True)
show("Lamb-Oseen vortex (unstretched)", r4,
     "P = 0 exactly at every r and t: omega is axial and S acts only in the r-theta plane")
require("Lamb-Oseen: div u = 0", sp.simplify(r4["div"]) == 0)
require("Lamb-Oseen: P = 0 exactly for all r, t", sp.simplify(r4["P"]) == 0)
require("Lamb-Oseen: A = 0 exactly (defined wherever S != 0 and omega != 0)",
        sp.simplify(r4["A"]) == 0)
record("Lamb-Oseen (unstretched)", r4, True, "0",
       "zeta varies with r and t (>0 in the core, <0 in the outer shear region); "
       "P = 0 and A = 0 everywhere, so the vortex only diffuses")

# zeta(r) profile for Lamb-Oseen in similarity variable eta = r / sqrt(4 nu t)
eta = sp.symbols("eta", positive=True)
subs_LO = {rsym: eta * sp.sqrt(4 * nu * t)}
zeta_LO = sp.simplify(r4["zeta"].subs(subs_LO))
zeta_LO = sp.simplify(sp.powsimp(zeta_LO, force=True))
print(f"  zeta_LO(eta) = {zeta_LO}")
f_zeta_LO = sp.lambdify(eta, zeta_LO.subs({nu: 1, t: 1, Gam: 1}), "numpy")
etas = np.linspace(0.02, 6.0, 400)
zLO = np.array([float(f_zeta_LO(e)) for e in etas])
sign_change = np.where(np.diff(np.sign(zLO)))[0]
eta_cross = float(etas[sign_change[0]]) if len(sign_change) else float("nan")
print(f"  zeta_LO changes sign at eta ~ {eta_cross:.4f} "
      f"(zeta > 0 inside, < 0 outside) while P = 0 throughout")
require("Lamb-Oseen: zeta sweeps both signs although P = 0 everywhere",
        zLO.max() > 0 > zLO.min(),
        f"zeta in [{zLO.min():.4f}, {zLO.max():.4f}], sign change at eta = {eta_cross:.4f}")

# 5. Burgers vortex
u_theta_B = Gam / (2 * sp.pi * rr) * (1 - sp.exp(-aa * rr ** 2 / (4 * nu)))
u5 = [-aa * x / 2 - u_theta_B * y / rr,
      -aa * y / 2 + u_theta_B * x / rr,
      aa * z]
r5 = analyse(u5, subs={y: 0, x: rsym}, simplify=True)
print("\n--- Burgers vortex ---")
print(f"  div u = {sp.simplify(r5['div'])}")
print(f"  P = {sp.simplify(r5['P'])}")
require("Burgers: div u = 0", sp.simplify(r5["div"]) == 0)

# exact statements at the axis r -> 0
lim = lambda e: sp.simplify(sp.limit(sp.simplify(e), rsym, 0, "+"))
nS2_0, nW2_0, w2_0, P_0 = (lim(r5[k]) for k in ("nS2", "nW2", "w2", "P"))
zeta_0 = sp.simplify((nW2_0 - nS2_0) / (nW2_0 + nS2_0))
A_0 = sp.simplify(P_0 / (sp.sqrt(nS2_0) * w2_0))
c_sym = Gam / (4 * sp.pi * nu)
print(f"  on the axis: ||S||^2 = {nS2_0}, |omega|^2 = {w2_0}, P = {P_0}")
print(f"  on the axis: zeta = {sp.simplify(zeta_0)}")
print(f"  on the axis: A    = {A_0}")
require("Burgers axis: ||S||_F^2 = (3/2) a^2 (pure background strain, shear vanishes)",
        sp.simplify(nS2_0 - sp.Rational(3, 2) * aa ** 2) == 0)
require("Burgers axis: A = sqrt(2/3) EXACTLY -- the sharp alignment bound is saturated",
        sp.simplify(A_0 - sp.sqrt(sp.Rational(2, 3))) == 0,
        "for every a > 0 and every Gamma, nu")
require("Burgers axis: P = a |omega|^2 > 0",
        sp.simplify(P_0 - aa * w2_0) == 0, f"P = {P_0}")
zeta_0_c = sp.simplify(zeta_0.subs(Gam, 4 * sp.pi * nu * sp.Symbol("c", positive=True)))
require("Burgers axis: zeta is INDEPENDENT of the axial strain a",
        sp.simplify(sp.diff(zeta_0, aa)) == 0,
        f"zeta_axis = {sp.simplify(zeta_0_c)} with c = Gamma/(4 pi nu); "
        "depends only on the circulation Reynolds number")
record("Burgers vortex (axis)", dict(r5, zeta=sp.simplify(zeta_0), P=P_0, A=A_0,
                                     a=sp.simplify(nS2_0 / (nS2_0 + nW2_0)),
                                     b=sp.simplify(nW2_0 / (nS2_0 + nW2_0)),
                                     h=sp.sqrt(sp.simplify(nS2_0 * nW2_0 /
                                                           (nS2_0 + nW2_0) ** 2)),
                                     L=sp.simplify((nW2_0 - nS2_0) / (2 * (nS2_0 + nW2_0))),
                                     D=sp.simplify(sp.Rational(1, 2) -
                                                   sp.sqrt(nS2_0 * nW2_0) /
                                                   (nS2_0 + nW2_0))),
       True, "+", "A = sqrt(2/3) exactly (bound saturated); zeta independent of a; "
                  "P = a|omega|^2 grows as a^3 at fixed similarity radius")

# self-similar coordinate: eta = r / sqrt(4 nu / a)
subs_B = {rsym: eta * sp.sqrt(4 * nu / aa)}
zeta_B = sp.simplify(sp.powsimp(sp.simplify(r5["zeta"].subs(subs_B)), force=True))
A_B = sp.simplify(sp.powsimp(sp.simplify(r5["A"].subs(subs_B)), force=True))
P_B = sp.simplify(sp.powsimp(sp.simplify(r5["P"].subs(subs_B)), force=True))
dz = sp.simplify(sp.diff(zeta_B, aa))
dA = sp.simplify(sp.diff(A_B, aa))
require("Burgers: at fixed similarity radius eta, zeta is independent of a",
        dz == 0, "d zeta / d a = 0 at fixed eta")
require("Burgers: at fixed similarity radius eta, A is independent of a",
        dA == 0, "d A / d a = 0 at fixed eta")
P_scaling = sp.simplify(sp.log(P_B).diff(aa) * aa)
require("Burgers: at fixed eta, P scales exactly as a^3",
        sp.simplify(P_scaling - 3) == 0,
        "so the entire (zeta, A) portrait is a-invariant while P varies as a^3: "
        "the scale coordinate carries all of the a dependence")

# numeric sweeps
Re_G = 100.0        # Gamma / nu
subs_num = {Gam: Re_G, nu: 1.0}
f_zeta_fixed_r = sp.lambdify((aa, rsym), r5["zeta"].subs(subs_num), "numpy")
f_A_fixed_r = sp.lambdify((aa, rsym), r5["A"].subs(subs_num), "numpy")
f_P_fixed_r = sp.lambdify((aa, rsym), r5["P"].subs(subs_num), "numpy")
f_Q_fixed_r = sp.lambdify((aa, rsym), r5["Qtot"].subs(subs_num), "numpy")
f_zeta_eta = sp.lambdify((aa, eta), zeta_B.subs(subs_num), "numpy")
f_A_eta = sp.lambdify((aa, eta), A_B.subs(subs_num), "numpy")
f_P_eta = sp.lambdify((aa, eta), P_B.subs(subs_num), "numpy")

a_vals = np.logspace(-2, 2, 41)
sweep_fixed_r = {}
for r_fix in (0.5, 1.0, 2.0):
    sweep_fixed_r[str(r_fix)] = {
        "a": a_vals.tolist(),
        "zeta": [float(f_zeta_fixed_r(av, r_fix)) for av in a_vals],
        "A": [float(f_A_fixed_r(av, r_fix)) for av in a_vals],
        "P": [float(f_P_fixed_r(av, r_fix)) for av in a_vals],
        "Q": [float(f_Q_fixed_r(av, r_fix)) for av in a_vals],
    }
sweep_eta = {}
for e_fix in (0.01, 0.5, 1.0, 2.0):
    sweep_eta[str(e_fix)] = {
        "a": a_vals.tolist(),
        "zeta": [float(f_zeta_eta(av, e_fix)) for av in a_vals],
        "A": [float(f_A_eta(av, e_fix)) for av in a_vals],
        "P": [float(f_P_eta(av, e_fix)) for av in a_vals],
    }
eta_grid = np.linspace(0.01, 5.0, 300)
profile = {"eta": eta_grid.tolist(),
           "zeta": [float(f_zeta_eta(1.0, e)) for e in eta_grid],
           "A": [float(f_A_eta(1.0, e)) for e in eta_grid],
           "P": [float(f_P_eta(1.0, e)) for e in eta_grid]}

print(f"\n  Burgers sweep at fixed physical radius r = 1.0, Re_Gamma = {Re_G}:")
for i in (0, 10, 20, 30, 40):
    av = a_vals[i]
    print(f"    a = {av:8.3f}   zeta = {sweep_fixed_r['1.0']['zeta'][i]:+.6f}   "
          f"A = {sweep_fixed_r['1.0']['A'][i]:+.6f}   "
          f"P = {sweep_fixed_r['1.0']['P'][i]:.6e}")
print("  Burgers sweep at fixed similarity radius eta = 1.0:")
for i in (0, 10, 20, 30, 40):
    av = a_vals[i]
    print(f"    a = {av:8.3f}   zeta = {sweep_eta['1.0']['zeta'][i]:+.6f}   "
          f"A = {sweep_eta['1.0']['A'][i]:+.6f}   "
          f"P = {sweep_eta['1.0']['P'][i]:.6e}")

zf = np.array(sweep_fixed_r["1.0"]["zeta"])
Pf = np.array(sweep_fixed_r["1.0"]["P"])
require("Burgers at fixed r: zeta decreases while P increases (sign-opposite trends)",
        np.all(np.diff(zf) < 0) and np.all(np.diff(Pf[:len(Pf) // 2]) > 0),
        "reproduces the RVP counterexample; the trend is a core-radius effect")
ze = np.array(sweep_eta["1.0"]["zeta"])
require("Burgers at fixed eta: zeta and A are exactly constant in a while P ~ a^3",
        np.ptp(ze) < 1e-12 and np.ptp(np.array(sweep_eta["1.0"]["A"])) < 1e-12,
        f"zeta spread {np.ptp(ze):.2e}, A spread "
        f"{np.ptp(np.array(sweep_eta['1.0']['A'])):.2e}")

# A(eta) monotone decrease away from the bound, and the Thales coordinates of the core
A_prof = np.array(profile["A"])
require("Burgers: A(eta) <= sqrt(2/3) with equality only on the axis",
        A_prof.max() <= sv.SQRT_2_3 + 1e-12 and
        abs(float(f_A_eta(1.0, 1e-6)) - sv.SQRT_2_3) < 1e-6,
        f"max A on the profile = {A_prof.max():.9f}, sqrt(2/3) = {sv.SQRT_2_3:.9f}")

# ---------------------------------------------------------------------------
summary_table = {
    "flows": rows,
    "lamb_oseen": {"zeta_eta_expr": str(zeta_LO), "eta_sign_change": eta_cross,
                   "zeta_min": float(zLO.min()), "zeta_max": float(zLO.max()),
                   "P": "0 exactly", "A": "0 exactly"},
    "burgers": {"axis_A": "sqrt(2/3) exactly", "axis_zeta_expr": str(sp.simplify(zeta_0_c)),
                "axis_zeta_depends_on_a": False,
                "P_scaling_at_fixed_eta": "a^3",
                "Re_Gamma": Re_G,
                "sweep_fixed_r": sweep_fixed_r, "sweep_fixed_eta": sweep_eta,
                "profile_at_a_1": profile},
    "failures": FAILURES,
}
(OUT / "canonical_flows.json").write_text(json.dumps(summary_table, indent=2))
print(f"\n{'ALL CHECKS PASSED' if not FAILURES else 'FAILURES: ' + str(FAILURES)}")
print(f"written: {OUT / 'canonical_flows.json'}")
sys.exit(0 if not FAILURES else 1)
