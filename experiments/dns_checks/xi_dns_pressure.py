"""XI EVOLUTION AUDIT -- how large is the nonlocal obstruction in a real flow?

Marches the parent audit's pseudo-spectral Taylor-Green solver (imported read-only
from ../../../src/dns_taylor_green.py) to a snapshot, solves the pressure Poisson
equation spectrally, forms H = Hess(p), and then:

  P1  verifies tr H = Lap p = E_W - E_S pointwise (a check on the solve)
  P2  verifies the derived law against the direct assembly from DA/Dt, pointwise
  P3  measures |pressure term| vs |local term| in D xi/Dt
  P4  measures how often dropping H_dev (i.e. restricted Euler) flips the SIGN of
      D xi/Dt -- the practical cost of the restricted-Euler closure
  P5  tests whether the pressure term is predictable from the local coordinates
      (Q, zeta, A, s) -- the empirical closure question

Run:  ../../../.venv/bin/python src/xi_dns_pressure.py [--n 96] [--re 400] [--t 8.0]
"""

from __future__ import annotations

import argparse
import json
import pathlib
import sys
import time

import numpy as np

# repo-root finder: works from src/ and from any experiments/ subdirectory
_ROOT = next(p for p in pathlib.Path(__file__).resolve().parents
             if (p / "src").is_dir())
sys.path.insert(0, str(_ROOT / "src"))
sys.path.insert(0, str(_ROOT / "experiments" / "dns_checks"))
import dns_taylor_green as dns  # noqa: E402  (module-level definitions only)
import statlib as st  # noqa: E402
import svcore as sv  # noqa: E402

OUT = _ROOT / "results"
OUT.mkdir(exist_ok=True)
FAILURES: list[str] = []


def require(name, ok, detail=""):
    print(f"[{'PASS' if ok else 'FAIL':4s}] {name}" + (f"   {detail}" if detail else ""))
    if not ok:
        FAILURES.append(name)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=96)
    ap.add_argument("--re", type=float, default=400.0)
    ap.add_argument("--t", type=float, default=8.0)
    ap.add_argument("--dt", type=float, default=0.02)
    args = ap.parse_args()

    nu = 1.0 / args.re
    sp = dns.Spectral(args.n, nu)
    uh = dns.taylor_green(sp)
    nsteps = int(round(args.t / args.dt))
    print(f"DNS: N = {args.n}^3, nu = {nu:g}, marching {nsteps} steps to t = {args.t}")
    t0 = time.time()
    for step in range(nsteps):
        uh = sp.step_rk4(uh, args.dt)
        if step % 100 == 0:
            print(f"   step {step:4d}  ({time.time() - t0:.0f}s)")

    # ---- fields at the snapshot
    G = sp.velocity_gradient(uh)                       # (N,N,N,3,3), (grad u)_ij
    S = sv.sym(G)
    W = sv.skew(G)
    w = sv.vorticity_from_skew(W)
    E_S = np.sum(S * S, axis=(-2, -1))
    E_W = np.sum(W * W, axis=(-2, -1))
    Q = E_S + E_W

    # ---- pressure Hessian, spectrally:  Lap p = -tr(A^2) = E_W - E_S
    f = E_W - E_S
    fh = sp.fft(f)
    K2 = np.where(sp.K2 == 0, 1.0, sp.K2)
    H = np.empty_like(G)
    for i in range(3):
        for j in range(3):
            # p_hat = -f_hat/k^2 ;  H_ij_hat = -k_i k_j p_hat = +k_i k_j f_hat / k^2
            H[..., i, j] = sp.ifft(sp.K[i] * sp.K[j] * fh / K2)
    trH = np.trace(H, axis1=-2, axis2=-1)
    rel_tr = float(np.max(np.abs(trH - f)) / np.max(np.abs(f)))
    require("P1 tr H = Lap p = E_W - E_S pointwise (pressure solve verified)",
            rel_tr < 1e-10, f"max relative residual {rel_tr:.2e}")

    # ---- viscous pieces, spectrally
    LapG = np.empty_like(G)
    for i in range(3):
        for j in range(3):
            LapG[..., i, j] = sp.ifft(-sp.K2 * sp.fft(G[..., i, j]))
    LapS, LapW = sv.sym(LapG), sv.skew(LapG)

    # ---- the three terms of D xi/Dt
    P = np.einsum("...i,...ij,...j->...", w, S, w)
    s_lr = sv.lund_rogers_s(S)
    nS = np.sqrt(E_S)
    xi = 0.5 * np.log(E_W / E_S)
    A = P / (nS * np.sum(w * w, axis=-1))
    H_dev = H - np.eye(3) * (trH / 3.0)[..., None, None]
    local = nS * ((1.0 + 0.5 * np.exp(2.0 * xi)) * A - s_lr / np.sqrt(6.0))
    press = np.sum(S * H_dev, axis=(-2, -1)) / E_S
    visc = nu * (np.sum(W * LapW, axis=(-2, -1)) / E_W -
                 np.sum(S * LapS, axis=(-2, -1)) / E_S)

    # ---- direct assembly from the gradient equation itself
    DA = -(G @ G) - H + nu * LapG
    dES = 2.0 * np.sum(S * sv.sym(DA), axis=(-2, -1))
    dEW = 2.0 * np.sum(W * sv.skew(DA), axis=(-2, -1))
    direct = 0.5 * (dEW / E_W - dES / E_S)
    tot = local + press + visc
    m = np.isfinite(direct) & np.isfinite(tot) & (E_S > 1e-12) & (E_W > 1e-12)
    rel = float(np.max(np.abs(tot[m] - direct[m]) / (np.abs(direct[m]) + 1e-9)))
    require("P2 the derived law reproduces D xi/Dt pointwise in the DNS field",
            rel < 1e-8,
            f"max relative error {rel:.2e} over {int(m.sum())} grid points")

    # ---- P3: how big is the nonlocal term?
    L, Pr, V = local[m], press[m], visc[m]
    ratio = np.abs(Pr) / (np.abs(L) + 1e-300)
    stats = {q: float(np.quantile(ratio, q)) for q in (0.1, 0.25, 0.5, 0.75, 0.9)}
    frac_dom = float(np.mean(np.abs(Pr) > np.abs(L)))
    print(f"\n   |pressure| / |local|  quantiles: " +
          ", ".join(f"{k:.0%}: {v:.3f}" for k, v in stats.items()))
    print(f"   pressure term is larger than the local term at {frac_dom:.1%} of points")
    print(f"   std(local) = {np.std(L):.4f}, std(pressure) = {np.std(Pr):.4f}, "
          f"std(viscous) = {np.std(V):.4f}")
    require("P3 the nonlocal term is of the same order as the local term",
            0.05 < np.median(ratio) < 20.0,
            f"median |pressure|/|local| = {np.median(ratio):.3f}")

    # ---- P4: does dropping H_dev flip the sign of D xi/Dt?
    flip = float(np.mean(np.sign(L + V) != np.sign(L + Pr + V)))
    require("P4 restricted Euler flips the sign of D xi/Dt at a substantial "
            "fraction of points",
            flip > 0.05,
            f"sign of D xi/Dt differs with and without H_dev at {flip:.1%} of points")

    # ---- P5: is the pressure term predictable from the local coordinates?
    sub = np.random.default_rng(0).choice(len(L), min(300_000, len(L)), replace=False)
    zt = ((E_W - E_S) / Q)[m][sub]
    Av, sv_, Qv = A[m][sub], s_lr[m][sub], Q[m][sub]
    Pr_s = Pr[sub]
    qb = lambda x, n: st.quantile_bins(x, n)[0]        # statlib returns (idx, edges)
    cell = ((qb(zt, 8) * 8 + qb(Av, 8)) * 8 + qb(sv_, 8)) * 8 + qb(np.log(Qv), 8)
    _, cidx = np.unique(cell, return_inverse=True)
    resid = st.group_mean_residual(Pr_s, cidx)
    r2 = 1.0 - float(np.sum(resid ** 2)) / float(np.sum((Pr_s - Pr_s.mean()) ** 2))
    print(f"\n   R^2 of the pressure term on 8^4 bins in (zeta, A, s, log Q): {r2:.4f}")
    require("P5 the pressure term is NOT determined by the local coordinates",
            r2 < 0.5,
            f"a saturated 4-way binning of the local state explains only "
            f"{r2:.1%} of the variance of S:H_dev/E_S")

    # ---- P6: the UNWEIGHTED volume integral of the pressure term vanishes
    # d_i d_j S_ij = Lap(div u) = 0, so  int S:H dx = int p d_i d_j S_ij dx = 0
    SH = np.sum(S * H_dev, axis=(-2, -1))
    mean_SH = float(np.mean(SH))
    scale_SH = float(np.mean(np.abs(SH)))
    require("P6 the volume mean of S:H_dev vanishes (integration by parts: "
            "d_i d_j S_ij = Lap div u = 0)",
            abs(mean_SH) / scale_SH < 1e-10,
            f"<S:H_dev> = {mean_SH:+.3e} against <|S:H_dev|> = {scale_SH:.4f} "
            f"(ratio {abs(mean_SH)/scale_SH:.2e})")
    mean_ES, mean_EW = float(np.mean(E_S)), float(np.mean(E_W))
    require("P6a <E_S> = <E_W> (the parent audit's homogeneity identity, re-verified)",
            abs(mean_EW - mean_ES) / mean_ES < 1e-10,
            f"<E_S> = {mean_ES:.6f}, <E_W> = {mean_EW:.6f}")
    weighted = float(np.mean(press[m]))
    require("P6b BUT the WEIGHTED mean <S:H_dev/E_S> does NOT vanish",
            abs(weighted) / scale_SH > 1e-6,
            f"<S:H_dev/E_S> = {weighted:+.5f}: the obstruction to the comparator is "
            "created entirely by the pointwise weight 1/E_S, not by the pressure "
            "term's own mean")

    out = {
        "parameters": {"n": args.n, "re": args.re, "t": args.t, "dt": args.dt},
        "P1_trace_residual": rel_tr,
        "P2_max_rel_error": rel,
        "P3": {"ratio_quantiles": stats, "frac_pressure_dominates": frac_dom,
               "std_local": float(np.std(L)), "std_pressure": float(np.std(Pr)),
               "std_viscous": float(np.std(V)),
               "mean_local": float(np.mean(L)), "mean_pressure": float(np.mean(Pr)),
               "mean_viscous": float(np.mean(V))},
        "P4_sign_flip_fraction": flip,
        "P5_r2_pressure_on_local_bins": r2,
        "P6": {"mean_S_Hdev": mean_SH, "mean_abs_S_Hdev": scale_SH,
               "mean_E_S": mean_ES, "mean_E_W": mean_EW,
               "mean_weighted_pressure_term": weighted},
        "failures": FAILURES,
    }
    (OUT / "xi_dns_pressure.json").write_text(json.dumps(out, indent=2, default=float))
    np.savez_compressed(OUT / "xi_dns_terms.npz",
                        local=L[sub], pressure=Pr[sub], viscous=V[sub],
                        zeta=zt, A=Av, s=sv_, Q=Qv)
    print(f"\n{'ALL CHECKS PASSED' if not FAILURES else 'FAILURES: ' + str(FAILURES)}")
    print(f"written: {OUT / 'xi_dns_pressure.json'}")
    return 0 if not FAILURES else 1


if __name__ == "__main__":
    sys.exit(main())
