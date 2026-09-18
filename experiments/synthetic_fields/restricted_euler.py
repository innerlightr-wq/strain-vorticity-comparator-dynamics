"""PHASE 7a -- dynamical test on an analytically tractable system.

The restricted Euler (Vieillefosse) equation for the velocity-gradient tensor
G = grad u along a fluid trajectory,

    dG/dt = -G^2 + (1/3) tr(G^2) I,

is the exact Euler velocity-gradient equation with the ANISOTROPIC part of the
pressure Hessian and the viscous term dropped.  It is a closed 8-dimensional ODE,
so it is exactly integrable numerically and, crucially, the enstrophy balance
along a trajectory reduces to

    d/dt (|omega|^2 / 2) = P = omega . S omega       (no dissipation, no transport),

which removes the confound that the brief warns about: here P > 0 really is the
whole of the enstrophy budget.  What restricted Euler does NOT contain is the
nonlocal pressure coupling, so nothing below is a claim about Navier-Stokes; the
DNS in src/dns_taylor_green.py carries that burden.

Tests
  T1  matched-zeta families: fix zeta_0 and the scale, sweep A_0, measure the
      finite-time enstrophy amplification.
  T2  variance decomposition of the amplification over (zeta_0, A_0).
  T3  the phase-2 C5 pair: identical (scale, zeta, A, P) but opposite strain
      state s -- do they evolve identically?

Run:  .venv/bin/python src/restricted_euler.py
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
import statlib as st  # noqa: E402
import svcore as sv  # noqa: E402

OUT = _ROOT / "results"
OUT.mkdir(exist_ok=True)
FAILURES: list[str] = []


def require(name, ok, detail=""):
    print(f"[{'PASS' if ok else 'FAIL':4s}] {name}" + (f"   {detail}" if detail else ""))
    if not ok:
        FAILURES.append(name)


I3 = np.eye(3)


def re_rhs(G):
    G2 = G @ G
    tr = np.trace(G2, axis1=-2, axis2=-1)
    return -G2 + (tr / 3.0)[..., None, None] * I3


def integrate(G0, t_end, dt=2e-4, cap=1e4, record_every=50):
    """RK4 integration with a magnitude cap (restricted Euler blows up)."""
    G = G0.copy()
    n = len(G)
    alive = np.ones(n, bool)
    t = 0.0
    times, ens, norms = [0.0], [sv.coordinates_from_gradient(G)["w2"].copy()], \
                        [np.linalg.norm(G, axis=(-2, -1)).copy()]
    t_blow = np.full(n, np.nan)
    step = 0
    while t < t_end - 1e-12:
        h = min(dt, t_end - t)
        k1 = re_rhs(G)
        k2 = re_rhs(G + 0.5 * h * k1)
        k3 = re_rhs(G + 0.5 * h * k2)
        k4 = re_rhs(G + h * k3)
        Gn = G + (h / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
        nrm = np.linalg.norm(Gn, axis=(-2, -1))
        newly = alive & (nrm > cap)
        t_blow[newly] = t + h
        alive &= ~newly
        G = np.where(alive[:, None, None], Gn, G)
        t += h
        step += 1
        if step % record_every == 0:
            times.append(t)
            ens.append(sv.coordinates_from_gradient(G)["w2"].copy())
            norms.append(np.linalg.norm(G, axis=(-2, -1)).copy())
    return dict(G=G, alive=alive, t_blow=t_blow, times=np.array(times),
                enstrophy=np.array(ens), norms=np.array(norms))


def strain_eigenvalues(phi):
    """All traceless unit-Frobenius strain states, one parameter.

    lambda_i = sqrt(2/3) cos(phi + 2 pi i / 3): traceless for every phi, with
    sum lambda_i^2 = 1.  phi in [0, 2 pi / 3) covers every strain state once up to
    relabelling (phi = 0 gives (2,-1,-1)/sqrt 6, phi = pi/3 gives (1,1,-2)/sqrt 6).
    """
    return np.sqrt(2.0 / 3.0) * np.cos(phi + 2.0 * np.pi * np.arange(3) / 3.0)


def sample_cos2_with_A(lam, A, rng):
    """cos^2 weights (c1,c2,c3) on the simplex with sum_i lambda_i c_i = A exactly.

    The feasible set is a segment; we sample the free coordinate uniformly on it,
    so the two invariants that (zeta, A) does not fix are genuinely randomized.
    """
    order = np.argsort(lam)                       # ascending: lo, mid, hi
    lo, mid, hi = order
    llo, lmid, lhi = lam[lo], lam[mid], lam[hi]
    if not (llo - 1e-12 <= A <= lhi + 1e-12):
        return None
    # c_mid free; c_hi = ((A - lmid c_mid) - llo (1 - c_mid)) / (lhi - llo)
    lows, highs = [], []
    for cm in (0.0, 1.0):                         # feasibility is affine in c_mid
        pass
    # solve 0 <= c_hi <= 1 - c_mid  for c_mid in [0, 1]
    den = lhi - llo
    cand = []
    for cm in np.linspace(0.0, 1.0, 2001):
        chi = ((A - lmid * cm) - llo * (1.0 - cm)) / den
        if -1e-12 <= chi <= 1.0 - cm + 1e-12:
            cand.append(cm)
    if not cand:
        return None
    cm = rng.uniform(min(cand), max(cand))
    chi = np.clip(((A - lmid * cm) - llo * (1.0 - cm)) / den, 0.0, 1.0)
    clo = max(0.0, 1.0 - cm - chi)
    c = np.zeros(3)
    c[lo], c[mid], c[hi] = clo, cm, chi
    return c / c.sum()


def build_state(zeta, A, seed=0, n=1, rng=None, random_frame=True, phi=None):
    """States with prescribed (zeta, A) and ||grad u||_F = 1, residual invariants random.

    ||S||_F^2 = a = (1-zeta)/2 and |omega|^2 = 2 b = 1 + zeta, so every state
    shares the same total gradient magnitude and hence the same time scale.
    """
    rng = rng or np.random.default_rng(seed)
    a = 0.5 * (1.0 - zeta)
    w2 = 1.0 + zeta
    out = []
    tries = 0
    while len(out) < n and tries < 200 * n:
        tries += 1
        ph = rng.uniform(0.0, 2.0 * np.pi / 3.0) if phi is None else phi
        lam_u = strain_eigenvalues(ph)             # unit Frobenius norm
        c = sample_cos2_with_A(lam_u, A, rng)
        if c is None:
            continue
        S = np.diag(lam_u * np.sqrt(a))
        w = np.sqrt(c * w2)
        if random_frame:
            Qr = np.linalg.qr(rng.normal(size=(3, 3)))[0]
            if np.linalg.det(Qr) < 0:
                Qr[:, 0] *= -1.0
            S = Qr @ S @ Qr.T
            w = Qr @ w
        out.append(S + sv.skew_from_vorticity(w))
    return np.stack(out) if out else None


print("=" * 72)
print("PHASE 7a: restricted Euler (Vieillefosse) dynamics")
print("=" * 72)

# sanity: the ODE reproduces d/dt(|omega|^2/2) = P
G_test = build_state(0.2, 0.4, seed=5, n=64)
c0 = sv.coordinates_from_gradient(G_test)
h = 1e-6
G_step = G_test + h * re_rhs(G_test)
c1 = sv.coordinates_from_gradient(G_step)
dedt = (0.5 * c1["w2"] - 0.5 * c0["w2"]) / h
require("restricted Euler reproduces d/dt(|omega|^2/2) = omega.S omega",
        np.max(np.abs(dedt - c0["P"]) / (np.abs(c0["P"]) + 1e-9)) < 1e-4,
        f"max relative error {np.max(np.abs(dedt - c0['P']) / (np.abs(c0['P']) + 1e-9)):.2e}")

# ---------------------------------------------------------------- T1 / T2
print("\n--- T1: matched-zeta families, sweeping A ---")
zeta_grid = np.array([-0.8, -0.5, -0.2, 0.0, 1 / 3, 0.6, 0.8])
A_grid = np.linspace(-0.40, 0.80, 13)
n_per = 80                      # residual invariants (strain state, alignment split)
Gs, meta = [], []
rng_states = np.random.default_rng(4242)
for iz, z0 in enumerate(zeta_grid):
    for ia, A0 in enumerate(A_grid):
        G = build_state(z0, A0, rng=rng_states, n=n_per)
        if G is None:
            continue
        Gs.append(G)
        meta += [(z0, A0, 0.0)] * len(G)
Gs = np.concatenate(Gs)
meta = np.array(meta)
c_init = sv.coordinates_from_gradient(Gs)
require("built states realize the prescribed (zeta, A)",
        np.max(np.abs(c_init["zeta"] - meta[:, 0])) < 1e-9 and
        np.max(np.abs(c_init["A"] - meta[:, 1])) < 1e-9,
        f"n = {len(Gs)} states, max errors "
        f"{np.max(np.abs(c_init['zeta'] - meta[:, 0])):.1e} / "
        f"{np.max(np.abs(c_init['A'] - meta[:, 1])):.1e}")

T_END = 0.6            # in units of 1/||grad u||_F(0); restricted Euler blows up later
res = integrate(Gs, T_END, dt=2e-4, cap=1e3)
e0 = c_init["w2"] / 2.0
eT = res["enstrophy"][-1] / 2.0
amp = np.log(np.maximum(eT, 1e-300) / np.maximum(e0, 1e-300))
print(f"  integrated {len(Gs)} restricted-Euler trajectories to t = {T_END} "
      f"({np.mean(~res['alive']):.1%} hit the magnitude cap)")

zeta0, A0s, ssign = meta[:, 0], meta[:, 1], meta[:, 2]
good = np.isfinite(amp)
Dz, _, _ = st.bin_dummies(zeta0[good], len(zeta_grid))
zeta_only = st.nested_r2(amp[good], [Dz], [np.zeros_like(amp[good])])
with_A = st.nested_r2(amp[good], [Dz], [A0s[good], A0s[good] ** 2])
# how much variance remains once BOTH zeta and A are fixed?  (the residual
# invariants that the (scale, zeta, A) triple discards)
cell = (np.round(zeta0[good], 6) * 1e6).astype(np.int64) * 10_000 + \
    (np.round(A0s[good], 6) * 1000).astype(np.int64)
_, cell_idx = np.unique(cell, return_inverse=True)
res_cell = st.group_mean_residual(amp[good], cell_idx)
frac_unexplained = float(np.sum(res_cell ** 2) / np.sum((amp[good] - amp[good].mean()) ** 2))
A_only = st.nested_r2(amp[good], [A0s[good], A0s[good] ** 2], [np.zeros_like(amp[good])])
print(f"  R^2 of ln(e_omega(T)/e_omega(0)) on zeta alone (7 bins): "
      f"{zeta_only['r2_base']:.4f}")
print(f"  R^2 on A alone (quadratic):                              {A_only['r2_base']:.4f}")
print(f"  R^2 on zeta bins + A:                                    {with_A['r2_full']:.4f}"
      f"   (delta over zeta alone: {with_A['delta_r2']:+.4f})")
print(f"  variance of the amplification left UNEXPLAINED by (zeta, A) jointly: "
      f"{frac_unexplained:.4%}  (residual strain state / alignment split)")
require("T1 at matched zeta, A controls the finite-time enstrophy amplification",
        with_A["delta_r2"] > 0.3 and with_A["r2_full"] > zeta_only["r2_base"] + 0.3,
        f"delta R^2 = {with_A['delta_r2']:+.4f} adding A to a saturated zeta model")

# per-(zeta) spread of the amplification at fixed zeta
print("\n  amplification ln(e(T)/e(0)) at fixed zeta, by A:")
table = []
for z0 in zeta_grid:
    sel = (zeta0 == z0) & good
    lo = amp[sel & (A0s < -0.3)]
    mid = amp[sel & (np.abs(A0s) < 0.05)]
    hi = amp[sel & (A0s > 0.69)]
    row = dict(zeta=float(z0),
               amp_A_neg=float(lo.mean()) if len(lo) else float("nan"),
               amp_A_zero=float(mid.mean()) if len(mid) else float("nan"),
               amp_A_high=float(hi.mean()) if len(hi) else float("nan"),
               spread=float(amp[sel].max() - amp[sel].min()))
    table.append(row)
    print(f"    zeta = {z0:+.3f}:  A<-0.2: {row['amp_A_neg']:+.4f}   "
          f"A~0: {row['amp_A_zero']:+.4f}   A>0.7: {row['amp_A_high']:+.4f}   "
          f"(full spread {row['spread']:.4f})")
require("T1 the sign of the amplification at fixed zeta follows the sign of A",
        all(r["amp_A_neg"] < r["amp_A_zero"] < r["amp_A_high"] for r in table),
        "monotone in A at every zeta tested")

# EXACT prediction to test: the production contribution to the RELATIVE growth rate
# is  P / (||grad u||_F e_omega) = A sqrt(2(1-zeta)),  so the A-sensitivity of the
# amplification should fall off as sqrt(1-zeta).
sens = []
for z0 in zeta_grid:
    sel = (zeta0 == z0) & good
    slope = np.polyfit(A0s[sel], amp[sel], 1)[0]
    sens.append(dict(zeta=float(z0), slope_dAmp_dA=float(slope),
                     predicted_shape=float(np.sqrt(2.0 * (1.0 - z0)))))
ratio = np.array([r["slope_dAmp_dA"] for r in sens]) / \
    np.array([r["predicted_shape"] for r in sens])
print("\n  d(amplification)/dA at fixed zeta vs the exact sqrt(2(1-zeta)) prediction:")
for r, q in zip(sens, ratio):
    print(f"    zeta = {r['zeta']:+.3f}: slope = {r['slope_dAmp_dA']:.4f}, "
          f"sqrt(2(1-zeta)) = {r['predicted_shape']:.4f}, ratio = {q:.4f}")
require("T2 the A-sensitivity of the amplification follows the exact "
        "sqrt(2(1-zeta)) allocation law",
        float(np.std(ratio) / np.mean(ratio)) < 0.12,
        f"slope / sqrt(2(1-zeta)) is constant to "
        f"{np.std(ratio) / np.mean(ratio):.2%} across zeta "
        f"(mean {np.mean(ratio):.4f})")

# size of the sufficiency gap: spread of the amplification WITHIN (zeta, A) cells
cell_spread = []
for b in range(cell_idx.max() + 1):
    m = cell_idx == b
    if m.sum() > 5:
        cell_spread.append(float(np.ptp(amp[good][m])))
cell_spread = np.array(cell_spread)
print(f"\n  within-(zeta, A)-cell spread of the amplification: median "
      f"{np.median(cell_spread):.4f}, max {cell_spread.max():.4f} nats "
      f"(vs a total range of {np.ptp(amp[good]):.4f})")

# ---------------------------------------------------------------- T3
print("\n--- T3: the phase-2 C5 pair (identical zeta, A, P; opposite strain state) ---")
cex = json.loads((OUT / "counterexamples.json").read_text())["C5_residual_invariants"]
S_I, w_I = np.array(cex["S_I"]), np.array(cex["w_I"])
S_II, w_II = np.array(cex["S_II"]), np.array(cex["w_II"])
G_pair = np.stack([S_I + sv.skew_from_vorticity(w_I),
                   S_II + sv.skew_from_vorticity(w_II)])
cp = sv.coordinates_from_gradient(G_pair)
require("T3 the pair is identical in (scale, zeta, A, P)",
        abs(cp["zeta"][0] - cp["zeta"][1]) < 1e-12 and
        abs(cp["A"][0] - cp["A"][1]) < 1e-12 and
        abs(cp["P"][0] - cp["P"][1]) < 1e-12 and
        abs(cp["Qtot"][0] - cp["Qtot"][1]) < 1e-12,
        f"zeta = {cp['zeta'][0]:.6f}, A = {cp['A'][0]:.6f}, P = {cp['P'][0]:.6f}")
res_pair = integrate(G_pair, 0.8, dt=1e-4, cap=1e4, record_every=200)
amp_pair = np.log(res_pair["enstrophy"][-1] / res_pair["enstrophy"][0])
traj = {"times": res_pair["times"].tolist(),
        "enstrophy_I": res_pair["enstrophy"][:, 0].tolist(),
        "enstrophy_II": res_pair["enstrophy"][:, 1].tolist(),
        "s_I": cex["s_I"], "s_II": cex["s_II"]}
print(f"  state I  (s = {cex['s_I']:+.1f}): ln amplification = {amp_pair[0]:+.6f}")
print(f"  state II (s = {cex['s_II']:+.1f}): ln amplification = {amp_pair[1]:+.6f}")
require("T3 FALSIFIER for sufficiency: identical (scale, zeta, A, P) states diverge "
        "dynamically",
        abs(amp_pair[0] - amp_pair[1]) > 0.01,
        f"amplification differs by {abs(amp_pair[0] - amp_pair[1]):.4f} nats -- "
        "(scale, zeta, A) is exactly sufficient for INSTANTANEOUS P, not for evolution")

# ---------------------------------------------------------------- persist
out = {
    "model": "restricted Euler (Vieillefosse): dG/dt = -G^2 + (1/3)tr(G^2) I",
    "caveat": "no nonlocal pressure Hessian, no viscosity; not a Navier-Stokes claim",
    "T_end": T_END, "n_trajectories": int(len(Gs)),
    "fraction_capped": float(np.mean(~res["alive"])),
    "r2_zeta_only": zeta_only["r2_base"], "r2_A_only": A_only["r2_base"],
    "r2_zeta_plus_A": with_A["r2_full"], "delta_r2_adding_A": with_A["delta_r2"],
    "variance_unexplained_by_zeta_and_A": frac_unexplained,
    "amplification_table": table,
    "A_sensitivity_vs_exact_law": sens,
    "A_sensitivity_ratio_cv": float(np.std(ratio) / np.mean(ratio)),
    "within_cell_spread_median": float(np.median(cell_spread)),
    "within_cell_spread_max": float(cell_spread.max()),
    "total_amplification_range": float(np.ptp(amp[good])),
    "C5_pair": {"amp_I": float(amp_pair[0]), "amp_II": float(amp_pair[1]),
                "difference": float(abs(amp_pair[0] - amp_pair[1])), "trajectory": traj},
    "failures": FAILURES,
}
(OUT / "restricted_euler.json").write_text(json.dumps(out, indent=2, default=float))
np.savez_compressed(OUT / "restricted_euler_samples.npz",
                    zeta0=zeta0, A0=A0s, s_sign=ssign, amplification=amp,
                    t_blow=res["t_blow"])
print(f"\n{'ALL CHECKS PASSED' if not FAILURES else 'FAILURES: ' + str(FAILURES)}")
print(f"written: {OUT / 'restricted_euler.json'}")
sys.exit(0 if not FAILURES else 1)
