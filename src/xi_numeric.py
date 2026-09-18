"""XI EVOLUTION AUDIT -- numerical verification of the exact D(xi)/Dt law.

Verifies, on random admissible states and on the nine special configurations:

  M1  the full identity, assembled from the gradient equation itself
  M2  the vorticity sector carries no pressure term (perturbing H changes nothing)
  M3  only H_dev acts: adding c*I to H leaves D xi/Dt invariant
  M4  restricted-Euler reduction (H_dev = 0)
  M5  CLOSURE COUNTEREXAMPLE 1: identical (Q, h, xi, A), different s
  M6  CLOSURE COUNTEREXAMPLE 2: identical (Q, h, xi, A, s), different H_dev
  M7  finite-difference sanity check (secondary)
  K   the special-states table
  L   conditioning of xi versus h

Run:  ../../../.venv/bin/python src/xi_numeric.py
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
import svcore as sv  # noqa: E402

OUT = _ROOT / "results"
OUT.mkdir(exist_ok=True)
FAILURES: list[str] = []
report: dict[str, object] = {}


def require(name, ok, detail=""):
    print(f"[{'PASS' if ok else 'FAIL':4s}] {name}" + (f"   {detail}" if detail else ""))
    if not ok:
        FAILURES.append(name)


sym = lambda M: 0.5 * (M + np.swapaxes(M, -1, -2))
skew = lambda M: 0.5 * (M - np.swapaxes(M, -1, -2))
ddot = lambda M, N: np.sum(M * N, axis=(-2, -1))


def state_coords(G):
    """(Q, E_S, E_W, zeta, h, xi, A, s, ||S||) for a stack of gradients."""
    S, W = sym(G), skew(G)
    w = sv.vorticity_from_skew(W)
    E_S, E_W = ddot(S, S), ddot(W, W)
    Q = E_S + E_W
    zeta = (E_W - E_S) / Q
    with np.errstate(divide="ignore", invalid="ignore"):
        xi = 0.5 * np.log(E_W / E_S)
        A = np.einsum("...i,...ij,...j->...", w, S, w) / (np.sqrt(E_S) * ddot(W, W) * 2)
    return dict(S=S, W=W, w=w, E_S=E_S, E_W=E_W, Q=Q, zeta=zeta,
                h=np.sqrt((E_S / Q) * (E_W / Q)), xi=xi, A=A,
                s=sv.lund_rogers_s(S), nS=np.sqrt(E_S))


def dxi_formula(c, H, nu=0.0, LapS=None, LapW=None):
    """The derived law:  ||S||[(1+e^{2xi}/2)A - s/sqrt6] + S:H_dev/E_S + viscous."""
    H_dev = H - np.eye(3) * (np.trace(H, axis1=-2, axis2=-1) / 3.0)[..., None, None]
    local = c["nS"] * ((1.0 + 0.5 * np.exp(2.0 * c["xi"])) * c["A"] -
                       c["s"] / np.sqrt(6.0))
    press = ddot(c["S"], H_dev) / c["E_S"]
    visc = 0.0
    if nu and LapS is not None:
        visc = nu * (ddot(c["W"], LapW) / c["E_W"] - ddot(c["S"], LapS) / c["E_S"])
    return local, press, visc


def dxi_direct(G, H, nu=0.0, LapS=None, LapW=None):
    """D xi/Dt assembled straight from  DA/Dt = -A^2 - H + nu Lap A."""
    c = state_coords(G)
    DA = -(G @ G) - H
    if nu and LapS is not None:
        DA = DA + nu * (LapS + LapW)
    dES = 2.0 * ddot(c["S"], sym(DA))
    dEW = 2.0 * ddot(c["W"], skew(DA))
    return 0.5 * (dEW / c["E_W"] - dES / c["E_S"]), dES, dEW


def random_states(n, rng, traceless=True):
    G = rng.normal(size=(n, 3, 3))
    if traceless:
        G -= np.eye(3) * (np.trace(G, axis1=-2, axis2=-1) / 3.0)[:, None, None]
    return G


def poisson_H(G, rng):
    """A random symmetric H obeying the trace constraint tr H = -tr(A^2)."""
    M = rng.normal(size=(len(G), 3, 3))
    H = sym(M)
    tr_target = -np.trace(G @ G, axis1=-2, axis2=-1)
    H = H - np.eye(3) * (np.trace(H, axis1=-2, axis2=-1) / 3.0)[:, None, None]
    return H + np.eye(3) * (tr_target / 3.0)[:, None, None]


print("=" * 76)
print("M1-M4  the exact law on random admissible states")
print("=" * 76)
rng = np.random.default_rng(20260918)
N = 200_000
G = random_states(N, rng)
H = poisson_H(G, rng)
nu = 0.37
LapS = sym(rng.normal(size=(N, 3, 3)))
LapW = skew(rng.normal(size=(N, 3, 3)))

c = state_coords(G)
require("trace constraint tr H = -tr(A^2) = E_W - E_S holds by construction",
        np.max(np.abs(np.trace(H, axis1=-2, axis2=-1) - (c["E_W"] - c["E_S"]))) < 1e-10,
        f"max residual {np.max(np.abs(np.trace(H, axis1=-2, axis2=-1) - (c['E_W'] - c['E_S']))):.2e}")

loc, pre, vis = dxi_formula(c, H, nu, LapS, LapW)
direct, dES, dEW = dxi_direct(G, H, nu, LapS, LapW)
rel = np.max(np.abs(loc + pre + vis - direct) / (np.abs(direct) + 1e-12))
require("M1 D xi/Dt = ||S||[(1+e^{2xi}/2)A - s/sqrt6] + S:H_dev/E_S + nu[...]",
        rel < 1e-9,
        f"max relative error {rel:.3e} over {N} random states (with viscous terms)")

# sector equations separately
P = np.einsum("...i,...ij,...j->...", c["w"], c["S"], c["w"])
dEW_formula = P + 2 * nu * ddot(c["W"], LapW)
trS3 = np.trace(c["S"] @ c["S"] @ c["S"], axis1=-2, axis2=-1)
H_dev = H - np.eye(3) * (np.trace(H, axis1=-2, axis2=-1) / 3.0)[:, None, None]
dES_formula = -2 * trS3 - 0.5 * P - 2 * ddot(c["S"], H) + 2 * nu * ddot(c["S"], LapS)
require("M1a D E_W/Dt = P + 2 nu W:Lap W   (no pressure)",
        np.max(np.abs(dEW - dEW_formula)) / np.max(np.abs(dEW)) < 1e-12,
        f"max relative error {np.max(np.abs(dEW - dEW_formula)) / np.max(np.abs(dEW)):.2e}")
require("M1b D E_S/Dt = -2 tr(S^3) - P/2 - 2 S:H + 2 nu S:Lap S",
        np.max(np.abs(dES - dES_formula)) / np.max(np.abs(dES)) < 1e-12,
        f"max relative error {np.max(np.abs(dES - dES_formula)) / np.max(np.abs(dES)):.2e}")

# M2: perturb H arbitrarily -> the rotation sector must not move
H2 = H + sym(rng.normal(size=(N, 3, 3)))
_, _, dEW2 = dxi_direct(G, H2, nu, LapS, LapW)
require("M2 the vorticity sector is exactly pressure-blind",
        np.max(np.abs(dEW2 - dEW)) < 1e-10,
        "an arbitrary symmetric perturbation of H leaves D E_W/Dt unchanged")

# M3: adding c I to H must not change D xi/Dt
cc = rng.normal(size=N)
H3 = H + np.eye(3) * cc[:, None, None]
d3, _, _ = dxi_direct(G, H3, nu, LapS, LapW)
require("M3 only H_dev acts: H -> H + c I leaves D xi/Dt invariant",
        np.max(np.abs(d3 - direct)) < 1e-9,
        f"max change {np.max(np.abs(d3 - direct)):.2e} for |c| up to "
        f"{np.max(np.abs(cc)):.2f}")

# M4: restricted Euler = drop H_dev (and nu)
H_iso = np.eye(3) * (np.trace(H, axis1=-2, axis2=-1) / 3.0)[:, None, None]
d_re, _, _ = dxi_direct(G, H_iso, 0.0)
loc_only, _, _ = dxi_formula(c, H_iso, 0.0)
require("M4 restricted Euler reduces exactly to the local law",
        np.max(np.abs(d_re - loc_only) / (np.abs(d_re) + 1e-12)) < 1e-9,
        "D xi/Dt|_RE = ||S||[(1+e^{2xi}/2)A - s/sqrt6], closing on (||S||, xi, A, s)")
omitted = np.abs(pre)
require("M4a the term restricted Euler omits is exactly S:H_dev/E_S",
        np.max(np.abs((direct - vis) - (loc + pre))) < 1e-9,
        f"median |omitted| = {np.median(omitted):.4f}, "
        f"median |local| = {np.median(np.abs(loc)):.4f} in this synthetic ensemble")
report["M1_M4"] = {"n": N, "max_rel_error_full_law": float(rel),
                   "median_abs_local": float(np.median(np.abs(loc))),
                   "median_abs_pressure": float(np.median(omitted))}

print()
print("=" * 76)
print("M5/M6  CLOSURE COUNTEREXAMPLES")
print("=" * 76)


def build_state(zeta, A_target, phi, scale=1.0):
    """Exact state with prescribed zeta, A and strain-state angle phi, ||grad u||=1."""
    lam = np.sqrt(2.0 / 3.0) * np.cos(phi + 2.0 * np.pi * np.arange(3) / 3.0)
    a = 0.5 * (1.0 - zeta)
    w2 = 1.0 + zeta
    order = np.argsort(lam)
    lo, mid, hi = order
    # omega in the (e_lo, e_hi) plane with cos^2 weight c on e_hi
    cw = (A_target - lam[lo]) / (lam[hi] - lam[lo])
    if not (-1e-12 <= cw <= 1 + 1e-12):
        return None
    cw = float(np.clip(cw, 0.0, 1.0))
    d = np.zeros(3)
    d[hi], d[lo] = np.sqrt(cw), np.sqrt(1.0 - cw)
    S = np.diag(lam * np.sqrt(a)) * scale
    w = d * np.sqrt(w2) * scale
    return S + sv.skew_from_vorticity(w)


# M5: same (Q, h, xi, A), different s
zeta0, A0 = 0.2, 0.30
GI = build_state(zeta0, A0, 0.0)          # one strain state
GII = build_state(zeta0, A0, np.pi / 3)   # the mirrored strain state
pair = np.stack([GI, GII])
cp = state_coords(pair)
H0 = np.zeros((2, 3, 3))
H0 = H0 + np.eye(3) * ((cp["E_W"] - cp["E_S"]) / 3.0)[:, None, None]   # H_dev = 0
locp, prep, _ = dxi_formula(cp, H0)
print(f"   state I : zeta={cp['zeta'][0]:+.6f} h={cp['h'][0]:.6f} xi={cp['xi'][0]:+.6f} "
      f"A={cp['A'][0]:+.6f} s={cp['s'][0]:+.4f}  ->  D xi/Dt = {locp[0]:+.6f}")
print(f"   state II: zeta={cp['zeta'][1]:+.6f} h={cp['h'][1]:.6f} xi={cp['xi'][1]:+.6f} "
      f"A={cp['A'][1]:+.6f} s={cp['s'][1]:+.4f}  ->  D xi/Dt = {locp[1]:+.6f}")
same_coords = (abs(cp["zeta"][0] - cp["zeta"][1]) < 1e-12 and
               abs(cp["h"][0] - cp["h"][1]) < 1e-12 and
               abs(cp["xi"][0] - cp["xi"][1]) < 1e-12 and
               abs(cp["A"][0] - cp["A"][1]) < 1e-12 and
               abs(cp["Q"][0] - cp["Q"][1]) < 1e-12)
require("M5 identical (Q, h, xi, A) but different s give different D xi/Dt",
        same_coords and abs(locp[0] - locp[1]) > 1e-3,
        f"s = {cp['s'][0]:+.4f} vs {cp['s'][1]:+.4f}; D xi/Dt differs by "
        f"{abs(locp[0] - locp[1]):.6f} with ZERO pressure Hessian -- the local part "
        "alone already fails to close on the parent's coordinate set")

# M6: same (Q, h, xi, A, s), different H_dev
G6 = np.stack([GI, GI])
c6 = state_coords(G6)
rng6 = np.random.default_rng(7)
Hd = sym(rng6.normal(size=(3, 3)))
Hd = Hd - np.eye(3) * np.trace(Hd) / 3.0
H6 = np.stack([np.zeros((3, 3)), Hd])
H6 = H6 + np.eye(3) * ((c6["E_W"] - c6["E_S"]) / 3.0)[:, None, None]
loc6, pre6, _ = dxi_formula(c6, H6)
require("M6 identical (Q, h, xi, A, s) but different H_dev give different D xi/Dt",
        abs(loc6[0] - loc6[1]) < 1e-12 and abs(pre6[0] - pre6[1]) > 1e-3,
        f"local parts identical ({loc6[0]:+.6f}); pressure parts {pre6[0]:+.6f} vs "
        f"{pre6[1]:+.6f} -- the obstruction is genuinely outside the gradient tensor")
report["M5_M6"] = {
    "M5": {"zeta": float(cp["zeta"][0]), "A": float(cp["A"][0]),
           "s_I": float(cp["s"][0]), "s_II": float(cp["s"][1]),
           "dxi_I": float(locp[0]), "dxi_II": float(locp[1])},
    "M6": {"local": float(loc6[0]), "pressure_I": float(pre6[0]),
           "pressure_II": float(pre6[1])}}

print()
print("=" * 76)
print("M7  finite-difference sanity check (secondary)")
print("=" * 76)
# central differences along the frozen-H gradient flow: second-order accurate,
# so the residual is a truncation floor rather than a disagreement
dt = 1e-5
Gs = G[:2000]
Hs = H[:2000]
cs = state_coords(Gs)
DA = -(Gs @ Gs) - Hs
xi_p = state_coords(Gs + dt * DA)["xi"]
xi_m = state_coords(Gs - dt * DA)["xi"]
dxi_fd = (xi_p - xi_m) / (2 * dt)
lo2, pr2, _ = dxi_formula(cs, Hs)
abs_err = np.abs(dxi_fd - (lo2 + pr2))
med_rel = float(np.median(abs_err / (np.abs(lo2 + pr2) + 1e-12)))
require("M7 central finite differences agree with the closed-form law",
        float(np.max(abs_err)) < 1e-6 and med_rel < 1e-9,
        f"max absolute error {np.max(abs_err):.2e}, median relative error "
        f"{med_rel:.2e} at dt = {dt} (O(dt^2) truncation floor)")
report["M7"] = {"dt": dt, "max_abs_error": float(np.max(abs_err)),
                "median_rel_error": med_rel}

print()
print("=" * 76)
print("K  SPECIAL STATES")
print("=" * 76)
table = []


def row(name, zeta, A, phi, note, Hdev_mode="zero"):
    G_ = build_state(zeta, A, phi)
    if G_ is None:
        table.append(dict(state=name, status="not admissible", note=note))
        return
    cx = state_coords(G_[None])
    Hx = np.zeros((1, 3, 3)) + np.eye(3) * ((cx["E_W"] - cx["E_S"]) / 3.0)[:, None, None]
    if Hdev_mode == "aligned":
        Hx = Hx + cx["S"] / np.sqrt(cx["E_S"])[:, None, None]      # B = +1
    elif Hdev_mode == "misaligned":
        lam_, vec_ = np.linalg.eigh(cx["S"][0])
        Hd_ = vec_ @ np.diag(lam_[[1, 2, 0]]) @ vec_.T
        Hd_ = Hd_ - np.eye(3) * np.trace(Hd_) / 3.0
        Hx = Hx + Hd_[None] / np.linalg.norm(Hd_)
    lo_, pr_, _ = dxi_formula(cx, Hx)
    nHd = float(np.linalg.norm(Hx[0] - np.eye(3) * np.trace(Hx[0]) / 3.0))
    B = (float(np.sum(cx["S"][0] * (Hx[0] - np.eye(3) * np.trace(Hx[0]) / 3.0))) /
         (np.sqrt(cx["E_S"][0]) * nHd)) if nHd > 1e-14 else 0.0
    table.append(dict(state=name, zeta=float(cx["zeta"][0]), xi=float(cx["xi"][0]),
                      h=float(cx["h"][0]), A=float(cx["A"][0]), s=float(cx["s"][0]),
                      local=float(lo_[0]), pressure=float(pr_[0]), B=float(B),
                      total=float(lo_[0] + pr_[0]), note=note))
    print(f"   {name:34s} zeta={cx['zeta'][0]:+.3f} xi={cx['xi'][0]:+.3f} "
          f"A={cx['A'][0]:+.4f} s={cx['s'][0]:+.3f} | local={lo_[0]:+.4f} "
          f"press={pr_[0]:+.4f} (B={B:+.3f})")


# a NON-degenerate strain state (phi = pi/6 gives lambda ~ (1,0,-1)/sqrt2, s = 0)
PHI_GEN = np.pi / 6
lam_gen = np.sqrt(2 / 3) * np.cos(PHI_GEN + 2 * np.pi * np.arange(3) / 3)
print(f"   generic strain state used for the alignment rows: "
      f"lambda/||S|| = {np.round(np.sort(lam_gen)[::-1], 4)}, s = "
      f"{sv.lund_rogers_s(np.diag(lam_gen)):+.3f}")
row("balanced, omega on e_max", 0.0, float(np.max(lam_gen)), PHI_GEN,
    "most extensional eigenvector; A attains lambda_max/||S||")
row("balanced, omega on e_int", 0.0, float(np.sort(lam_gen)[1]), PHI_GEN,
    "intermediate eigenvector; here lambda_int = 0 so A = 0")
row("balanced, omega on e_min", 0.0, float(np.min(lam_gen)), PHI_GEN,
    "most compressive eigenvector; A attains lambda_min/||S||")
row("strain-rich zeta = -0.8", -0.8, 0.5, 0.0, "approaching the pure-strain endpoint")
row("rotation-rich zeta = +0.8", 0.8, 0.5, 0.0,
    "approaching the pure-rotation endpoint; note the e^{2xi} amplification")
row("balanced, H_dev = 0", 0.0, 0.5, 0.0, "restricted-Euler state", "zero")
row("balanced, H_dev aligned with S", 0.0, 0.5, 0.0, "B = +1, maximal coupling",
    "aligned")
row("balanced, H_dev eigen-permuted", 0.0, 0.5, 0.0,
    "same spectrum, permuted eigenframe: |B| < 1", "misaligned")
report["K_special_states"] = table

# endpoint limits, treated analytically
print("\n   endpoint limits (analytic):")
print("     zeta -> -1 (pure strain):  xi -> -inf, e^{2xi} -> 0,  "
      "D xi/Dt -> ||S||[A - s/sqrt6] + S:H_dev/E_S   (finite; xi itself is singular)")
print("     zeta -> +1 (pure rotation): e^{2xi} = E_W/E_S -> inf and "
      "||S|| e^{2xi} = E_W/||S|| -> inf, so D xi/Dt DIVERGES unless A -> 0")
report["endpoint_limits"] = {
    "zeta_to_minus_1": "coordinate singularity of xi; D xi/Dt has a finite limit "
                       "||S||[A - s/sqrt6] + S:H_dev/E_S along any fixed omega direction",
    "zeta_to_plus_1": "genuine divergence: ||S|| e^{2xi} A = (E_W/||S||) A -> infinity "
                      "as the strain sector is extinguished at fixed enstrophy"}

print()
print("=" * 76)
print("L  CONDITIONING")
print("=" * 76)
for zt in (0.0, 0.5, 0.9, 0.99):
    dxi_dz = 1 / (1 - zt ** 2)
    dh_dz = -zt / (2 * np.sqrt(1 - zt ** 2))
    print(f"   zeta = {zt:5.2f}:  d xi/d zeta = {dxi_dz:10.4f} (= R = 1/(4h^2))    "
          f"d h/d zeta = {dh_dz:+8.4f}")
require("L xi is regular at the fold where h is degenerate",
        abs(1 / (1 - 0.0 ** 2) - 1.0) < 1e-12 and abs(-0.0 / 2) < 1e-12,
        "d xi/d zeta = 1 at zeta = 0, while dh/d zeta = 0 there")
require("L xi stretches the endpoints",
        1 / (1 - 0.99 ** 2) > 50,
        f"d xi/d zeta = {1 / (1 - 0.99 ** 2):.1f} at zeta = 0.99")
report["L_conditioning"] = {str(zt): {"dxi_dzeta": 1 / (1 - zt ** 2),
                                      "dh_dzeta": -zt / (2 * np.sqrt(1 - zt ** 2))}
                            for zt in (0.0, 0.5, 0.9, 0.99)}

report["failures"] = FAILURES
(OUT / "xi_numeric.json").write_text(json.dumps(report, indent=2, default=float))
print(f"\n{'ALL CHECKS PASSED' if not FAILURES else 'FAILURES: ' + str(FAILURES)}")
print(f"written: {OUT / 'xi_numeric.json'}")
sys.exit(0 if not FAILURES else 1)
