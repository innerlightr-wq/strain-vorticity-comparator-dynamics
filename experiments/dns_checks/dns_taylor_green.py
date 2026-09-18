"""PHASE 7b -- Navier-Stokes DNS: does A predict enstrophy evolution beyond zeta?

A standard dealiased pseudo-spectral solver for the incompressible Navier-Stokes
equations on the 2pi-periodic box, started from the Taylor-Green vortex.  At
selected times the FULL pointwise enstrophy balance is evaluated spectrally,

    D_t (|omega|^2 / 2) = P - nu |grad omega|^2 + nu Laplacian(|omega|^2 / 2),
                          ^P    ^eps_omega        ^T_nu  (transport, sign-indefinite)

with no forcing (f = 0), so the production term is never confused with net growth.
Two tests are then run:

  D1  conditional test: inside joint bins of (zeta, log ||grad u||_F) -- i.e. at
      matched allocation AND matched scale -- does A predict the normalized
      material rate of change of enstrophy?  Is the answer changed by whether the
      viscous terms offset the production?
  D2  Lagrangian test: particles seeded at t_seed are advected with the flow and
      their enstrophy amplification over a finite window is compared across
      A-quartiles inside matched (zeta, scale) cells.

Run:  .venv/bin/python src/dns_taylor_green.py  [--n 64] [--re 400] [--tend 12]
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
import statlib as st  # noqa: E402
import svcore as sv  # noqa: E402

OUT = _ROOT / "results"
OUT.mkdir(exist_ok=True)
FAILURES: list[str] = []


def require(name, ok, detail=""):
    print(f"[{'PASS' if ok else 'FAIL':4s}] {name}" + (f"   {detail}" if detail else ""))
    if not ok:
        FAILURES.append(name)


# ----------------------------------------------------------------- solver setup

class Spectral:
    def __init__(self, N, nu):
        self.N, self.nu = N, nu
        kx = np.fft.fftfreq(N, 1.0 / N)
        kz = np.fft.rfftfreq(N, 1.0 / N)
        self.KX, self.KY, self.KZ = np.meshgrid(kx, kx, kz, indexing="ij")
        self.K = np.stack([self.KX, self.KY, self.KZ])
        self.K2 = np.sum(self.K ** 2, axis=0)
        self.K2inv = np.where(self.K2 == 0, 1.0, self.K2)
        self.kmax = N // 3                          # 2/3 rule, spherical
        self.mask = (np.sqrt(self.K2) <= self.kmax).astype(float)
        self.x = np.linspace(0, 2 * np.pi, N, endpoint=False)

    def fft(self, f):
        return np.fft.rfftn(f, axes=(0, 1, 2))

    def ifft(self, fh):
        return np.fft.irfftn(fh, s=(self.N,) * 3, axes=(0, 1, 2))

    def vfft(self, f):
        return np.stack([self.fft(f[i]) for i in range(3)])

    def vifft(self, fh):
        return np.stack([self.ifft(fh[i]) for i in range(3)])

    def curl_h(self, uh):
        K = self.K
        return 1j * np.stack([K[1] * uh[2] - K[2] * uh[1],
                              K[2] * uh[0] - K[0] * uh[2],
                              K[0] * uh[1] - K[1] * uh[0]])

    def project(self, fh):
        kdotf = np.sum(self.K * fh, axis=0)
        return fh - self.K * (kdotf / self.K2inv)

    def rhs(self, uh):
        u = self.vifft(uh)
        w = self.vifft(self.curl_h(uh))
        nl = np.stack([u[1] * w[2] - u[2] * w[1],       # u x omega (rotational form)
                       u[2] * w[0] - u[0] * w[2],
                       u[0] * w[1] - u[1] * w[0]])
        nlh = self.vfft(nl) * self.mask
        return self.project(nlh) - self.nu * self.K2 * uh

    def step_rk4(self, uh, dt):
        k1 = self.rhs(uh)
        k2 = self.rhs(uh + 0.5 * dt * k1)
        k3 = self.rhs(uh + 0.5 * dt * k2)
        k4 = self.rhs(uh + dt * k3)
        return uh + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)

    def gradient(self, fh):
        """grad of a scalar field given its transform: returns 3 real fields."""
        return np.stack([self.ifft(1j * self.K[j] * fh) for j in range(3)])

    def velocity_gradient(self, uh):
        """(grad u)_ij = d u_i / d x_j as a (N,N,N,3,3) real array."""
        N = self.N
        G = np.empty((N, N, N, 3, 3))
        for i in range(3):
            for j in range(3):
                G[..., i, j] = self.ifft(1j * self.K[j] * uh[i])
        return G

    def laplacian(self, f):
        return self.ifft(-self.K2 * self.fft(f))


def taylor_green(sp):
    X, Y, Z = np.meshgrid(sp.x, sp.x, sp.x, indexing="ij")
    u = np.stack([np.sin(X) * np.cos(Y) * np.cos(Z),
                  -np.cos(X) * np.sin(Y) * np.cos(Z),
                  np.zeros_like(X)])
    return sp.project(sp.vfft(u)) * sp.mask


def interp_trilinear(field, pos, N):
    """Periodic trilinear interpolation of a scalar field at positions pos (n,3)."""
    g = pos / (2 * np.pi / N)
    i0 = np.floor(g).astype(int)
    f = g - i0
    i0 %= N
    i1 = (i0 + 1) % N
    out = np.zeros(len(pos))
    for dx in (0, 1):
        for dy in (0, 1):
            for dz in (0, 1):
                wx = f[:, 0] if dx else 1 - f[:, 0]
                wy = f[:, 1] if dy else 1 - f[:, 1]
                wz = f[:, 2] if dz else 1 - f[:, 2]
                ix = i1[:, 0] if dx else i0[:, 0]
                iy = i1[:, 1] if dy else i0[:, 1]
                iz = i1[:, 2] if dz else i0[:, 2]
                out += wx * wy * wz * field[ix, iy, iz]
    return out


def interp_vector(fields, pos, N):
    return np.stack([interp_trilinear(fields[i], pos, N) for i in range(3)], axis=-1)


# -------------------------------------------------------------- balance fields

def balance_fields(sp, uh):
    """Pointwise enstrophy balance and audit coordinates for the whole box."""
    G = sp.velocity_gradient(uh)
    S = sv.sym(G)
    w = sv.vorticity_from_skew(sv.skew(G))
    co = sv.coordinates(S, w)
    e_w = 0.5 * co["w2"]
    # eps_omega = nu |grad omega|^2 ; T_nu = nu Laplacian(e_omega)
    wh = sp.curl_h(uh)
    gsum = np.zeros_like(e_w)
    for i in range(3):
        gi = sp.gradient(wh[i])
        gsum += np.sum(gi ** 2, axis=0)
    eps_w = sp.nu * gsum
    T_nu = sp.nu * sp.laplacian(e_w)
    Dt_e = co["P"] - eps_w + T_nu
    return dict(G=G, S=S, w=w, co=co, e_w=e_w, eps_w=eps_w, T_nu=T_nu, Dt_e=Dt_e)


def snapshot_analysis(sp, uh, tag):
    """Conditional test D1 on a single snapshot."""
    B = balance_fields(sp, uh)
    co, e_w = B["co"], B["e_w"]
    flat = lambda a: np.asarray(a).reshape(-1)
    zeta, A, P = flat(co["zeta"]), flat(co["A"]), flat(co["P"])
    Q = flat(co["Qtot"])
    ew, eps, Tn, Dte = flat(e_w), flat(B["eps_w"]), flat(B["T_nu"]), flat(B["Dt_e"])
    m = np.isfinite(A) & (ew > 1e-12) & (Q > 1e-12)
    zeta, A, P, Q, ew, eps, Tn, Dte = (v[m] for v in (zeta, A, P, Q, ew, eps, Tn, Dte))
    nrm = np.sqrt(Q) * ew                       # dimensional scale of a rate * enstrophy
    g_tot = Dte / nrm                           # dimensionless material growth rate
    g_P = P / nrm                               # exactly A sqrt(2(1-zeta))
    g_visc = (-eps + Tn) / nrm
    # exact identity check
    id_err = float(np.max(np.abs(g_P - A * np.sqrt(2.0 * (1.0 - zeta)))))
    # joint conditioning on (zeta, log Q)
    iz, _ = st.quantile_bins(zeta, 10)
    iq, _ = st.quantile_bins(np.log(Q), 10)
    cell = iz * 10 + iq
    _, cidx = np.unique(cell, return_inverse=True)
    res = lambda v: st.group_mean_residual(v, cidx)
    rA, rgt, rgP, rgv = res(A), res(g_tot), res(g_P), res(g_visc)
    corr = lambda a, b: float(np.sum(a * b) / np.sqrt(np.sum(a ** 2) * np.sum(b ** 2)))
    ss = lambda v: float(np.sum(v ** 2))
    out = {
        "tag": tag, "n_points": int(len(zeta)),
        "identity_max_error_gP_minus_A_sqrt": id_err,
        "mean_P": float(P.mean()), "mean_eps": float(eps.mean()),
        "mean_T_nu": float(Tn.mean()), "mean_Dt_e": float(Dte.mean()),
        "frac_P_positive": float((P > 0).mean()),
        "frac_Dt_e_positive": float((Dte > 0).mean()),
        "frac_P_pos_but_Dt_e_neg": float(((P > 0) & (Dte < 0)).mean()),
        "mean_A": float(A.mean()), "mean_zeta": float(zeta.mean()),
        "partial_corr_A_with_g_tot_given_zeta_and_scale": corr(rA, rgt),
        "partial_corr_A_with_g_P_given_zeta_and_scale": corr(rA, rgP),
        "partial_corr_A_with_g_visc_given_zeta_and_scale": corr(rA, rgv),
        "delta_r2_A_on_g_tot": corr(rA, rgt) ** 2 * ss(rgt) /
                               float(np.sum((g_tot - g_tot.mean()) ** 2)),
        "r2_cells_only_g_tot": 1.0 - ss(rgt) /
                               float(np.sum((g_tot - g_tot.mean()) ** 2)),
        "var_ratio_visc_to_P": float(np.var(g_visc) / np.var(g_P)),
    }
    # conditional means of the growth rate by A-quartile inside matched cells
    qa = np.zeros(len(A), int)
    for b in range(cidx.max() + 1):
        s_ = cidx == b
        if s_.sum() >= 20:
            qa[s_] = np.clip(np.searchsorted(
                np.quantile(A[s_], [0.25, 0.5, 0.75]), A[s_], side="right"), 0, 3)
    out["by_A_quartile_within_cells"] = [
        {"quartile": int(q), "n": int((qa == q).sum()),
         "mean_A": float(A[qa == q].mean()),
         "mean_g_tot": float(g_tot[qa == q].mean()),
         "mean_g_P": float(g_P[qa == q].mean()),
         "mean_g_visc": float(g_visc[qa == q].mean()),
         "frac_Dt_e_positive": float((Dte[qa == q] > 0).mean())}
        for q in range(4)]
    sample = np.random.default_rng(0).choice(len(zeta), min(80_000, len(zeta)),
                                             replace=False)
    arrays = {f"{tag}_{k}": v[sample] for k, v in
              dict(zeta=zeta, A=A, P=P, Q=Q, e_w=ew, eps=eps, T_nu=Tn, Dt_e=Dte,
                   g_tot=g_tot, g_P=g_P, g_visc=g_visc).items()}
    return out, arrays, B


# ------------------------------------------------------------------------ main

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=64)
    ap.add_argument("--re", type=float, default=400.0)
    ap.add_argument("--tend", type=float, default=12.0)
    ap.add_argument("--dt", type=float, default=0.02)
    ap.add_argument("--seed-time", type=float, default=7.0)
    ap.add_argument("--window", type=float, default=1.5)
    ap.add_argument("--n-part", type=int, default=40_000)
    ap.add_argument("--snaps", type=float, nargs="*", default=[4.0, 6.0, 8.0, 10.0])
    args = ap.parse_args()

    nu = 1.0 / args.re
    sp = Spectral(args.n, nu)
    uh = taylor_green(sp)
    print(f"DNS: N = {args.n}^3, nu = {nu:g} (Re = {args.re:g}), dt = {args.dt}, "
          f"t_end = {args.tend}, k_max = {sp.kmax}")

    snap_times = list(args.snaps)
    snaps, arrays, history = [], {}, []
    rng = np.random.default_rng(7)
    parts = None
    t = 0.0
    t0 = time.time()
    nsteps = int(round(args.tend / args.dt))
    for step in range(nsteps + 1):
        # -------- diagnostics
        u = sp.vifft(uh)
        w = sp.vifft(sp.curl_h(uh))
        E = 0.5 * float(np.mean(np.sum(u ** 2, axis=0)))
        Zt = 0.5 * float(np.mean(np.sum(w ** 2, axis=0)))
        eps_tot = 2.0 * nu * Zt                      # = nu <|omega|^2> for periodic
        history.append(dict(t=t, energy=E, enstrophy=Zt, dissipation=eps_tot))

        for ts in list(snap_times):
            if abs(t - ts) < 0.5 * args.dt:
                s_out, s_arr, _ = snapshot_analysis(sp, uh, f"t{ts:g}")
                snaps.append(s_out)
                arrays.update(s_arr)
                snap_times.remove(ts)
                print(f"  t = {t:5.2f}: snapshot analysed "
                      f"(<P> = {s_out['mean_P']:.4f}, <eps_w> = {s_out['mean_eps']:.4f}, "
                      f"partial corr(A, g_tot | zeta, scale) = "
                      f"{s_out['partial_corr_A_with_g_tot_given_zeta_and_scale']:+.4f})")

        # -------- Lagrangian seeding / tracking
        if parts is None and t >= args.seed_time - 0.5 * args.dt:
            B = balance_fields(sp, uh)
            pos = rng.uniform(0, 2 * np.pi, size=(args.n_part, 3))
            co = B["co"]
            flat3 = lambda f: f
            init = {}
            for k in ("zeta", "A", "Qtot", "w2", "P"):
                init[k] = interp_trilinear(np.ascontiguousarray(co[k]), pos, sp.N)
            init["e_w"] = 0.5 * init.pop("w2")
            init["Dt_e"] = interp_trilinear(np.ascontiguousarray(B["Dt_e"]), pos, sp.N)
            parts = dict(pos=pos, init=init, t_seed=t)
            print(f"  t = {t:5.2f}: seeded {args.n_part} Lagrangian particles")
        elif parts is not None and t < parts["t_seed"] + args.window - 0.5 * args.dt:
            u_at = interp_vector(u, parts["pos"], sp.N)
            mid = (parts["pos"] + 0.5 * args.dt * u_at) % (2 * np.pi)
            u_mid = interp_vector(u, mid, sp.N)           # RK2 (u frozen within step)
            parts["pos"] = (parts["pos"] + args.dt * u_mid) % (2 * np.pi)
        elif parts is not None and "final" not in parts:
            w2f = np.sum(sp.vifft(sp.curl_h(uh)) ** 2, axis=0)
            parts["final"] = {"e_w": 0.5 * interp_trilinear(np.ascontiguousarray(w2f),
                                                            parts["pos"], sp.N)}
            parts["t_final"] = t
            print(f"  t = {t:5.2f}: Lagrangian window closed")

        if step == nsteps:
            break
        uh = sp.step_rk4(uh, args.dt)
        t += args.dt
        if step % 50 == 0:
            print(f"    step {step:4d}  t = {t:6.3f}  E = {E:.6f}  Z = {Zt:.4f}  "
                  f"({time.time() - t0:.0f}s)")

    # ---------------------------------------------------------- resolution check
    hist = history
    eps_peak = max(h["dissipation"] for h in hist)
    eta = (nu ** 3 / eps_peak) ** 0.25
    kmax_eta = sp.kmax * eta
    u_rms = np.sqrt(2.0 / 3.0 * hist[0]["energy"] * 2)
    print(f"\n  resolution: peak dissipation {eps_peak:.5f} at "
          f"t = {max(hist, key=lambda h: h['dissipation'])['t']:.2f}, "
          f"eta = {eta:.4f}, k_max * eta = {kmax_eta:.3f}")
    require("DNS is adequately resolved (k_max * eta >= 1)", kmax_eta >= 1.0,
            f"k_max * eta = {kmax_eta:.3f} (>= 1 is the usual adequacy criterion)")

    # ------------------------------------------------------------------ D1 test
    print("\n--- D1: at matched (zeta, scale), does A predict the material rate? ---")
    for s_ in snaps:
        print(f"  {s_['tag']}: partial corr(A, g_tot) = "
              f"{s_['partial_corr_A_with_g_tot_given_zeta_and_scale']:+.4f}, "
              f"partial corr(A, g_P) = "
              f"{s_['partial_corr_A_with_g_P_given_zeta_and_scale']:+.4f}, "
              f"partial corr(A, g_visc) = "
              f"{s_['partial_corr_A_with_g_visc_given_zeta_and_scale']:+.4f}, "
              f"delta R^2 = {s_['delta_r2_A_on_g_tot']:.4f}")
        print(f"      <P> = {s_['mean_P']:+.5f}, <eps_w> = {s_['mean_eps']:.5f}, "
              f"<T_nu> = {s_['mean_T_nu']:+.2e}, <Dt e> = {s_['mean_Dt_e']:+.5f}; "
              f"P > 0 at {s_['frac_P_positive']:.1%} of points, "
              f"P > 0 but Dt e < 0 at {s_['frac_P_pos_but_Dt_e_neg']:.1%}")
        for q in s_["by_A_quartile_within_cells"]:
            print(f"      A-quartile {q['quartile']}: <A> = {q['mean_A']:+.4f}  "
                  f"<g_tot> = {q['mean_g_tot']:+.5f}  <g_P> = {q['mean_g_P']:+.5f}  "
                  f"<g_visc> = {q['mean_g_visc']:+.5f}  "
                  f"Dt e > 0 for {q['frac_Dt_e_positive']:.1%}")
    require("D1 the exact identity g_P = A sqrt(2(1-zeta)) holds pointwise in the DNS",
            bool(snaps) and max(s_["identity_max_error_gP_minus_A_sqrt"] for s_ in snaps) < 1e-9,
            f"max error {max(s_['identity_max_error_gP_minus_A_sqrt'] for s_ in snaps):.2e}")
    require("D1 the PROBABILITY of net enstrophy growth increases monotonically with A "
            "inside matched (zeta, scale) cells",
            all([q["frac_Dt_e_positive"] for q in s_["by_A_quartile_within_cells"]] ==
                sorted([q["frac_Dt_e_positive"] for q in s_["by_A_quartile_within_cells"]])
                for s_ in snaps),
            "; ".join(s_["tag"] + ": " +
                      "->".join(f"{q['frac_Dt_e_positive']:.0%}"
                                for q in s_["by_A_quartile_within_cells"])
                      for s_ in snaps))
    print("  NOTE: the Pearson partial correlation of A with the PER-ENSTROPHY rate is "
          "near zero\n        at late times because that normalization is singular as "
          "e_omega -> 0.\n        src/dns_reanalysis.py repeats the test with "
          "D_t e / ||grad u||_F^3 and with\n        rank statistics; see RESULTS.md.")
    require("D1 CONTROL: production alone does not determine net enstrophy growth",
            all(s_["frac_P_pos_but_Dt_e_neg"] > 0.05 for s_ in snaps),
            "a nontrivial fraction of points have P > 0 and D_t e_omega < 0")

    # ------------------------------------------------------------------ D2 test
    lag = {}
    if parts is not None and "final" in parts:
        print("\n--- D2: Lagrangian enstrophy amplification over a finite window ---")
        i0 = parts["init"]
        amp = np.log(np.maximum(parts["final"]["e_w"], 1e-300) /
                     np.maximum(i0["e_w"], 1e-300))
        m = np.isfinite(amp) & np.isfinite(i0["A"]) & (i0["e_w"] > 1e-10)
        amp, zA, AA, QQ = amp[m], i0["zeta"][m], i0["A"][m], i0["Qtot"][m]
        iz, _ = st.quantile_bins(zA, 6)
        iq, _ = st.quantile_bins(np.log(QQ), 6)
        _, cidx = np.unique(iz * 6 + iq, return_inverse=True)
        rA = st.group_mean_residual(AA, cidx)
        ra = st.group_mean_residual(amp, cidx)
        pc = float(np.sum(rA * ra) / np.sqrt(np.sum(rA ** 2) * np.sum(ra ** 2)))
        ss_tot = float(np.sum((amp - amp.mean()) ** 2))
        r2_cells = 1.0 - float(np.sum(ra ** 2)) / ss_tot
        dr2 = pc ** 2 * float(np.sum(ra ** 2)) / ss_tot
        quart = np.zeros(len(AA), int)
        for b in range(cidx.max() + 1):
            s_ = cidx == b
            if s_.sum() >= 20:
                quart[s_] = np.clip(np.searchsorted(
                    np.quantile(AA[s_], [0.25, 0.5, 0.75]), AA[s_], side="right"), 0, 3)
        rows = [{"quartile": int(q), "n": int((quart == q).sum()),
                 "mean_A": float(AA[quart == q].mean()),
                 "mean_amplification": float(amp[quart == q].mean()),
                 "median_amplification": float(np.median(amp[quart == q])),
                 "frac_growing": float((amp[quart == q] > 0).mean())} for q in range(4)]
        lag = {"t_seed": parts["t_seed"], "t_final": parts["t_final"],
               "n_particles": int(m.sum()),
               "partial_corr_A_amplification_given_zeta_scale": pc,
               "r2_cells_only": r2_cells, "delta_r2_adding_A": dr2,
               "by_A_quartile": rows}
        print(f"  window t = {parts['t_seed']:.2f} -> {parts['t_final']:.2f}, "
              f"{m.sum()} particles")
        for r in rows:
            print(f"    A-quartile {r['quartile']}: <A> = {r['mean_A']:+.4f}  "
                  f"<ln e(T)/e(0)> = {r['mean_amplification']:+.5f}  "
                  f"(median {r['median_amplification']:+.5f}, "
                  f"growing: {r['frac_growing']:.1%})")
        print(f"  partial corr(A, amplification | zeta, scale) = {pc:+.4f}, "
              f"delta R^2 = {dr2:.4f}")
        require("D2 at matched (zeta, scale), initial A predicts finite-time "
                "Lagrangian enstrophy amplification",
                pc > 0.1 and rows[3]["mean_amplification"] > rows[0]["mean_amplification"],
                f"partial corr {pc:+.4f}; quartile means "
                f"{[round(r['mean_amplification'], 4) for r in rows]}")

    out = {"parameters": vars(args), "nu": nu, "kmax": int(sp.kmax),
           "resolution": {"eta": eta, "kmax_eta": kmax_eta, "eps_peak": eps_peak},
           "history": hist, "snapshots": snaps, "lagrangian": lag,
           "failures": FAILURES}
    (OUT / "dns_taylor_green.json").write_text(json.dumps(out, indent=2, default=float))
    np.savez_compressed(OUT / "dns_samples.npz", **arrays)
    if lag:
        np.savez_compressed(OUT / "dns_lagrangian.npz",
                            zeta0=parts["init"]["zeta"], A0=parts["init"]["A"],
                            Q0=parts["init"]["Qtot"], e0=parts["init"]["e_w"],
                            eT=parts["final"]["e_w"])
    print(f"\n{'ALL CHECKS PASSED' if not FAILURES else 'FAILURES: ' + str(FAILURES)}")
    print(f"written: {OUT / 'dns_taylor_green.json'}")
    return 0 if not FAILURES else 1


if __name__ == "__main__":
    sys.exit(main())
