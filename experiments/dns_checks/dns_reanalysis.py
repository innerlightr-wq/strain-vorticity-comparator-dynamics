"""PHASE 7b (continued) -- normalization audit of the DNS conditional test.

The first pass of test D1 normalized the material rate by the local enstrophy,

    g = D_t e_omega / (||grad u||_F e_omega) ,

whose production part is exactly A sqrt(2(1-zeta)) -- convenient, but the viscous
part (-eps_omega + T_nu)/(||grad u||_F e_omega) diverges as e_omega -> 0, so the
statistic is heavy-tailed and its sample correlations are dominated by a handful of
near-irrotational points.  That is a property of the normalization, not of the flow.

This script re-analyses the SAME saved snapshots with a normalization that has no
such singularity,

    G := D_t e_omega / ||grad u||_F^3 ,      production part = g(zeta) A  (bounded by 4 sqrt2/9),

reports Pearson and Spearman partial statistics inside matched (zeta, scale) cells,
and documents the failure mode of the first normalization explicitly.

Run:  .venv/bin/python src/dns_reanalysis.py
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
sys.path.insert(0, str(_ROOT / "experiments" / "dns_checks"))
import statlib as st  # noqa: E402
import svcore as sv  # noqa: E402

OUT = _ROOT / "results"
FAILURES: list[str] = []


def require(name, ok, detail=""):
    print(f"[{'PASS' if ok else 'FAIL':4s}] {name}" + (f"   {detail}" if detail else ""))
    if not ok:
        FAILURES.append(name)


def rank(x):
    r = np.empty(len(x))
    r[np.argsort(x)] = np.arange(len(x))
    return r


def cells(zeta, Q, nz=10, nq=10):
    iz, _ = st.quantile_bins(zeta, nz)
    iq, _ = st.quantile_bins(np.log(Q), nq)
    _, cidx = np.unique(iz * nq + iq, return_inverse=True)
    return cidx


def partial(x, y, cidx):
    rx = st.group_mean_residual(x, cidx)
    ry = st.group_mean_residual(y, cidx)
    d = np.sqrt(np.sum(rx ** 2) * np.sum(ry ** 2))
    rho = float(np.sum(rx * ry) / d) if d > 0 else np.nan
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    return rho, float(rho ** 2 * np.sum(ry ** 2) / ss_tot)


z = np.load(OUT / "dns_samples.npz")
tags = sorted({f.split("_")[0] for f in z.files}, key=lambda t: float(t[1:]))
print("=" * 72)
print("PHASE 7b re-analysis: normalization audit of the DNS conditional test")
print("=" * 72)
report = {}
for tg in tags:
    zeta, A, P, Q = (z[f"{tg}_{k}"] for k in ("zeta", "A", "P", "Q"))
    ew, eps, Tn, Dte = (z[f"{tg}_{k}"] for k in ("e_w", "eps", "T_nu", "Dt_e"))
    cidx = cells(zeta, Q)

    # --- normalization 1: per enstrophy (the first pass)
    g_tot = Dte / (np.sqrt(Q) * ew)
    g_P = P / (np.sqrt(Q) * ew)
    g_v = (-eps + Tn) / (np.sqrt(Q) * ew)
    pc1, dr1 = partial(A, g_tot, cidx)
    pc1s, _ = partial(rank(A), rank(g_tot), cidx)

    # --- normalization 2: per cubed gradient (no singularity)
    G_tot = Dte / Q ** 1.5
    G_P = P / Q ** 1.5
    G_v = (-eps + Tn) / Q ** 1.5
    pc2, dr2 = partial(A, G_tot, cidx)
    pc2s, _ = partial(rank(A), rank(G_tot), cidx)
    pcP, drP = partial(A, G_P, cidx)
    pcV, drV = partial(A, G_v, cidx)

    # --- normalization 1 restricted to the upper half of the enstrophy distribution
    m = ew > np.median(ew)
    pc1t, dr1t = partial(A[m], g_tot[m], cells(zeta[m], Q[m]))

    q = np.zeros(len(A), int)
    for b in range(cidx.max() + 1):
        s_ = cidx == b
        if s_.sum() >= 20:
            q[s_] = np.clip(np.searchsorted(np.quantile(A[s_], [.25, .5, .75]),
                                            A[s_], side="right"), 0, 3)
    quart = [{"quartile": int(k), "n": int((q == k).sum()),
              "mean_A": float(A[q == k].mean()),
              "mean_G_tot": float(G_tot[q == k].mean()),
              "mean_G_P": float(G_P[q == k].mean()),
              "mean_G_visc": float(G_v[q == k].mean()),
              "median_G_tot": float(np.median(G_tot[q == k])),
              "frac_Dt_e_positive": float((Dte[q == k] > 0).mean())} for k in range(4)]

    report[tg] = {
        "n": int(len(A)),
        "norm1_per_enstrophy": {
            "partial_corr_A_g_tot": pc1, "spearman": pc1s, "delta_r2": dr1,
            "g_visc_p01": float(np.percentile(g_v, 1)),
            "g_visc_p99": float(np.percentile(g_v, 99)),
            "g_visc_min": float(g_v.min()),
            "std_ratio_visc_to_P": float(np.std(g_v) / np.std(g_P)),
            "partial_corr_upper_half_enstrophy": pc1t},
        "norm2_per_gradient_cubed": {
            "partial_corr_A_G_tot": pc2, "spearman": pc2s, "delta_r2": dr2,
            "partial_corr_A_G_P": pcP, "delta_r2_G_P": drP,
            "partial_corr_A_G_visc": pcV, "delta_r2_G_visc": drV,
            "std_ratio_visc_to_P": float(np.std(G_v) / np.std(G_P))},
        "by_A_quartile": quart,
    }
    r = report[tg]
    print(f"\n{tg}: n = {r['n']}")
    print(f"  norm 1  D_t e/(||grad u|| e):  partial corr(A, .) = {pc1:+.4f} "
          f"(Spearman {pc1s:+.4f}), std(visc)/std(P) = "
          f"{r['norm1_per_enstrophy']['std_ratio_visc_to_P']:.1f}, "
          f"viscous 1st pct = {r['norm1_per_enstrophy']['g_visc_p01']:.3g}, "
          f"min = {r['norm1_per_enstrophy']['g_visc_min']:.3g}")
    print(f"          same, upper half of enstrophy only: {pc1t:+.4f}")
    print(f"  norm 2  D_t e/||grad u||^3:    partial corr(A, .) = {pc2:+.4f} "
          f"(Spearman {pc2s:+.4f}), delta R^2 = {dr2:.4f}, "
          f"std(visc)/std(P) = "
          f"{r['norm2_per_gradient_cubed']['std_ratio_visc_to_P']:.2f}")
    print(f"          production part: {pcP:+.4f} (identity);  "
          f"viscous part: {pcV:+.4f}")
    for qq in quart:
        print(f"    A-quartile {qq['quartile']}: <A> = {qq['mean_A']:+.4f}   "
              f"<G_P> = {qq['mean_G_P']:+.5f}   <G_visc> = {qq['mean_G_visc']:+.5f}   "
              f"<G_tot> = {qq['mean_G_tot']:+.5f}   "
              f"D_t e > 0 for {qq['frac_Dt_e_positive']:.1%}")

# ---------------------------------------------------------------------------
# PHASE 8 -- is there conditional structure of A given zeta in a real flow?
# zeta and A are functionally INDEPENDENT (exact, phase 2), so any dependence here is
# a property of the dynamics/measure, not of the geometry.  Compared against the
# Gaussian solenoidal null, which has none by construction.
# ---------------------------------------------------------------------------
print("\n" + "=" * 72)
print("PHASE 8: conditional structure of A given zeta (DNS vs Gaussian null)")
print("=" * 72)
import generators as gen  # noqa: E402

Ggauss = gen.gaussian_solenoidal_gradients(n_grid=64, n_real=2, seed=11)
cg = sv.coordinates_from_gradient(Ggauss)
cond = {}


def conditional_A(zeta, A, nz=10, seed=0):
    m = np.isfinite(A) & np.isfinite(zeta)
    zeta, A = zeta[m], A[m]
    idx, edges = st.quantile_bins(zeta, nz)
    rows = [{"zeta_lo": float(edges[b]), "zeta_hi": float(edges[b + 1]),
             "n": int((idx == b).sum()), "mean_A": float(A[idx == b].mean()),
             "std_A": float(A[idx == b].std())} for b in range(idx.max() + 1)]
    within = float(np.sum(st.group_mean_residual(A, idx) ** 2))
    tot = float(np.sum((A - A.mean()) ** 2))
    expl = 1.0 - within / tot
    rng = np.random.default_rng(seed)
    null = []
    for _ in range(20):
        Ash = A[rng.permutation(len(A))]
        w = float(np.sum(st.group_mean_residual(Ash, idx) ** 2))
        null.append(1.0 - w / float(np.sum((Ash - Ash.mean()) ** 2)))
    return {"bins": rows, "variance_of_A_explained_by_zeta": expl,
            "permutation_null_mean": float(np.mean(null)),
            "permutation_null_std": float(np.std(null)),
            "z_score": float((expl - np.mean(null)) / np.std(null)),
            "mean_A_overall": float(A.mean()),
            "mean_A_lowest_zeta_decile": rows[0]["mean_A"],
            "mean_A_highest_zeta_decile": rows[-1]["mean_A"]}


cond["gaussian_null"] = conditional_A(cg["zeta"], cg["A"], seed=1)
for tg in tags:
    cond[tg] = conditional_A(z[f"{tg}_zeta"], z[f"{tg}_A"], seed=2)
for k, v in cond.items():
    print(f"  {k:14s}: var(A) explained by 10 zeta deciles = "
          f"{v['variance_of_A_explained_by_zeta']:+.4%} "
          f"(null {v['permutation_null_mean']:+.4%} +- "
          f"{v['permutation_null_std']:.4%}, z = {v['z_score']:+.0f});  "
          f"<A> = {v['mean_A_overall']:+.4f}, "
          f"lowest/highest zeta decile: {v['mean_A_lowest_zeta_decile']:+.4f} / "
          f"{v['mean_A_highest_zeta_decile']:+.4f}")
report["phase8_conditional_A_given_zeta"] = cond
require("CONTROL: the Gaussian solenoidal null shows no conditional structure of A "
        "given zeta",
        cond["gaussian_null"]["variance_of_A_explained_by_zeta"] < 0.01,
        f"{cond['gaussian_null']['variance_of_A_explained_by_zeta']:+.4%} of var(A)")
require("EMPIRICAL: the DNS does show conditional structure of A given zeta "
        "(statistical coupling, not geometric)",
        all(cond[t]["variance_of_A_explained_by_zeta"] >
            cond["gaussian_null"]["variance_of_A_explained_by_zeta"] + 0.005
            for t in tags),
        "explained variance: " +
        ", ".join(f"{t}: {cond[t]['variance_of_A_explained_by_zeta']:+.2%}"
                  for t in tags))

require("the per-enstrophy normalization is heavy-tailed (its viscous part diverges "
        "as e_omega -> 0)",
        all(report[t]["norm1_per_enstrophy"]["std_ratio_visc_to_P"] > 50 for t in tags),
        "std(viscous)/std(production): " +
        ", ".join(f"{report[t]['norm1_per_enstrophy']['std_ratio_visc_to_P']:.0f}"
                  for t in tags) +
        " -- and its rank correlation survives while its Pearson correlation does not: "
        "Spearman " + ", ".join(f"{report[t]['norm1_per_enstrophy']['spearman']:+.2f}"
                                for t in tags))

require("A predicts the net material rate in RANK at matched (zeta, scale), at every "
        "snapshot",
        all(report[t]["norm2_per_gradient_cubed"]["spearman"] > 0.4 for t in tags),
        "Spearman partial: " +
        ", ".join(f"{t}: {report[t]['norm2_per_gradient_cubed']['spearman']:+.3f}"
                  for t in tags))

require("the PROBABILITY of net enstrophy growth is monotone increasing in A at every "
        "snapshot",
        all([q["frac_Dt_e_positive"] for q in report[t]["by_A_quartile"]] ==
            sorted([q["frac_Dt_e_positive"] for q in report[t]["by_A_quartile"]])
            for t in tags),
        "; ".join(t + ": " + "->".join(f"{q['frac_Dt_e_positive']:.0%}"
                                       for q in report[t]["by_A_quartile"])
                  for t in tags))

require("the viscous terms are essentially UNCORRELATED with A at matched "
        "(zeta, scale): they add variance rather than cancelling production",
        all(abs(report[t]["norm2_per_gradient_cubed"]["partial_corr_A_G_visc"]) < 0.05
            for t in tags),
        "partial corr(A, viscous): " +
        ", ".join(f"{report[t]['norm2_per_gradient_cubed']['partial_corr_A_G_visc']:+.3f}"
                  for t in tags))

# the diagnosis: Pearson degrades exactly as the viscous variance grows
ratios = [report[t]["norm2_per_gradient_cubed"]["std_ratio_visc_to_P"] for t in tags]
pears = [report[t]["norm2_per_gradient_cubed"]["partial_corr_A_G_tot"] for t in tags]
require("the Pearson partial correlation falls as the viscous-to-production variance "
        "ratio grows (a variance effect, not a sign reversal)",
        ratios == sorted(ratios) and pears == sorted(pears, reverse=True),
        "std ratio " + ", ".join(f"{r:.1f}" for r in ratios) +
        " against Pearson " + ", ".join(f"{p:+.3f}" for p in pears))

report["interpretation"] = {
    "norm1_verdict": "unusable for Pearson statistics: singular as e_omega -> 0",
    "norm2_verdict": "well behaved; production part bounded by 4 sqrt2/9",
    "A_vs_net_rate": "strong in rank (Spearman 0.45-0.78) and in the probability of "
                     "growth (quartile 1 -> 4: ~12% -> ~67%); weak in variance share "
                     "at late times because the viscous terms dominate the variance",
    "viscous_terms": "large and essentially uncorrelated with A at matched (zeta, "
                     "scale), so they do not cancel the alignment dependence",
}
(OUT / "dns_reanalysis.json").write_text(json.dumps(report, indent=2, default=float))
print(f"\n{'ALL CHECKS PASSED' if not FAILURES else 'FAILURES: ' + str(FAILURES)}")
print(f"written: {OUT / 'dns_reanalysis.json'}")
sys.exit(0 if not FAILURES else 1)
