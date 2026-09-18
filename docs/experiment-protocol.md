# EXPERIMENT — protocol, as pre-registered and as executed

> *Reproduced verbatim from the audit run that produced it; only file paths have
> been updated to this repository's layout. Statements of the form "nothing outside
> this directory was modified" refer to that original run.*

The protocol was fixed by the brief before any computation, and is reproduced here
with what was actually run against each phase. Everything is deterministic: every
generator is seeded and every seed is recorded below, so all numbers in
`RESULTS.md` are reproducible with the commands in `README.md`.

## Design decisions taken before looking at any output

1. **Exact first, numerics second.** Phases 1–3 are pure algebra and were settled
   symbolically (`sympy`) before any ensemble was generated, so no numerical result
   could shape a "derivation".
2. **Undefined stays undefined.** `A` is `NaN` wherever `||S||_F = 0` or `|omega| = 0`
   (`svcore.coordinates` takes `eps = 0` by default). No regularized denominators
   anywhere — the regularized Omega-method style `a + b + eps` is not used, because
   it would silently invent values at exactly the degenerate states that matter.
3. **Normalization declared in advance.** Because `P/(||S||_F |omega|^2) = A`
   identically, that ratio is never presented as an empirical correlation. The
   scale-free production variable used throughout is
   `p := P/||grad u||_F^3`, which is a function of *all three* coordinates and
   therefore not circular. Where a rate is wanted, the dimensionless material rate
   `g := D_t e_omega / (||grad u||_F e_omega)` is used, whose production part is
   exactly `A sqrt(2(1-zeta))`.
4. **Deterministic dependence must be quotiented out, not "controlled for".**
   `h`, `D`, `R`, `eta_T` are exact functions of `zeta`, so the correct null for "does
   `h` add information" is *exactly zero*, and any nonzero estimate is an artifact.
   Two artifacts were anticipated and are measured rather than assumed away:
   * **binning leakage** — conditioning on finite `zeta` strata lets `h` resolve
     `zeta` inside a stratum. Diagnostic: refine the conditioning and check that the
     signal decays (it must, if it is leakage) while a genuine predictor's does not.
   * **plug-in MI bias** — conditional mutual information with many strata has a
     positive bias of order (cells)/(2·samples) per stratum. Diagnostic: subtract a
     within-stratum-shuffled null, which measures exactly that floor.
5. **No constant hunting.** No search for `1/e`, `pi/8`, `phi`, or any other
   constant was performed. Two numbers arose from derivations (`zeta = 1/3` and
   `4 sqrt 2/9`); both are stated with their proofs, and the near-coincidence between
   the derived `b = 2/3` and the framework's cubic landmark `b_* = 0.68233` is
   reported as a **non**-identity, with the exact residual `-1/27`.

## Phase-by-phase execution

| phase | script | what it does |
|---|---|---|
| 1 | `src/exact_algebra.py` | symbolic decomposition, Thales map, information content |
| 2 | `src/counterexamples.py` | exact witnesses, surjectivity construction, dof count |
| 3 | `src/exact_algebra.py` | three-factor factorization, sharp bounds, both allocation laws |
| 4 | `src/canonical_flows.py` | five canonical flows + Burgers sweeps, all symbolic |
| 5 | `experiments/synthetic_fields/ensembles.py` | four synthetic ensembles, conditional spreads, sign prediction |
| 6 | `experiments/synthetic_fields/ensembles.py` | Thales identities, determinism, CMI, residual-altitude trap |
| 7a | `experiments/synthetic_fields/restricted_euler.py` | restricted-Euler trajectories, matched-`zeta` families |
| 7b | `experiments/dns_checks/dns_taylor_green.py` | Navier–Stokes DNS, full enstrophy balance, Lagrangian test |
| 7b' | `experiments/dns_checks/dns_reanalysis.py` | normalization audit of the conditional test (see below) |
| 8 | `experiments/dns_checks/homogeneity_check.py` | is apex occupancy dynamical or kinematic? |
| 8 | `experiments/dns_checks/dns_reanalysis.py` | conditional structure of `A` given `zeta`, DNS vs Gaussian null |
| 8 | `src/exact_algebra.py` | envelope, bounded efficiencies, shared degeneracy at `zeta = ±1` |
| 9 | `experiments/synthetic_fields/ensembles.py` | six adversarial controls |

## Ensembles (Phase 5) and seeds

| tag | construction | n | seed |
|---|---|---|---|
| `E1_gaussian_solenoidal` | gradients of a divergence-free isotropic Gaussian field, spectrum `E(k) ∝ k^4 exp(-2k^2/k_0^2)`, `k_0 = 4`, `64^3` grid × 4 realizations | 1 048 576 | 11 |
| `E2_naive_independent` | random traceless symmetric `S`, independent isotropic `omega`, lognormal magnitudes (the construction suggested in the brief) | 400 000 | 12 |
| `E3_zeta_designed` | `zeta ~ U(-0.995, 0.995)`, `||S||_F = 1`, isotropic orientation — the ensemble for "at fixed `zeta`, what remains?" | 400 000 | 13 |
| `E4_aligned_max` / `E4_aligned_mid` | as E3 but `omega` tilted (`kappa = 4`) toward the most-extensional / intermediate strain eigenvector | 400 000 each | 14 / 15 |

E1 is the physically meaningful measure (a real random incompressible field);
E2 is the naive one; E3 and E4 are *designs*, used to answer conditional
questions without the confound of an ensemble-specific `zeta` distribution. All
four are reported separately and never pooled.

## Phase 7a — restricted Euler

`dG/dt = -G^2 + (1/3) tr(G^2) I`, RK4, `dt = 2e-4`, integrated to `t = 0.6` in units
of the initial `1/||grad u||_F` (all initial states normalized to
`||grad u||_F = 1`, so times are comparable across `zeta`). Enstrophy balance in this
model is exactly `d e_omega/dt = P` — no dissipation, no transport, no forcing — which
is why it is used for the matched-family test. Its omissions (the anisotropic
pressure Hessian, viscosity) are the reason nothing from it is reported as a
Navier–Stokes result.

Grid: `zeta_0 ∈ {-0.8,-0.5,-0.2,0,1/3,0.6,0.8}` × `A_0 ∈ linspace(-0.40, 0.80, 13)`,
80 states per cell, with the **two residual invariants randomized** inside each cell
(strain state `phi` uniform on `[0, 2 pi/3)`, and the alignment split sampled
uniformly along the feasible segment of the `cos^2` simplex at fixed `A`). Random
global rotations are also applied, but note they are dynamically trivial — the
equation is rotation-equivariant — so they are not what makes the cells nontrivial.
Seed 4242.

## Phase 7b — Navier–Stokes DNS

Dealiased pseudo-spectral solver, `2 pi`-periodic box, rotational form of the
nonlinear term, spherical `2/3`-rule truncation, RK4, **unforced**, Taylor–Green
initial condition. Production parameters: `96^3`, `Re = 1/nu = 400`, `dt = 0.02`,
`t_end = 12`. Resolution is reported as `k_max * eta_K` with
`eta_K = (nu^3/eps_peak)^{1/4}`, and the run is only used if `k_max eta_K >= 1`.

At `t = 4, 6, 8, 10` the **full pointwise balance** is evaluated spectrally:

```
D_t e_omega = P - nu |grad omega|^2 + nu Laplacian(e_omega) ,     e_omega = |omega|^2/2
```

The transport term `T_nu` is retained and reported; the brief's prohibition on
equating `P > 0` with net growth is enforced by reporting `P`, `-eps_omega + T_nu`
and `D_t e_omega` separately in every table, plus the fraction of points where
`P > 0` while `D_t e_omega < 0`.

* **D1 (conditional):** points are binned jointly into 10 `zeta` quantiles × 10
  `log ||grad u||_F^2` quantiles — matched allocation *and* matched scale — and the
  partial correlation of `A` with the normalized material rate is computed by exact
  group-mean residualization inside those cells.
* **D2 (Lagrangian):** 40 000 particles seeded uniformly at `t = 7.0`, advected with
  the flow (RK2, trilinear interpolation) to `t = 8.5`, and their enstrophy
  amplification `ln[e_omega(t_1)/e_omega(t_0)]` compared across `A`-quartiles inside
  matched `(zeta, scale)` cells. This is a finite-time test, not an instantaneous one.

## Phase 9 — the six controls

1. rotate `omega` isotropically in a fixed strain eigenframe, all magnitudes fixed
   (200 000 draws);
2. same, read as "preserve `zeta`, change alignment";
3. scale `S` and `omega` together — preserves `A`, `zeta`, and `p`, moves `P` over
   `2e16`;
4. randomize strain eigenvectors at fixed eigenvalues (50 000 random `SO(3)` frames);
5. isotropic vs preferentially aligned ensembles, including the Gaussian-field null;
6. partial correlation of `h` with `|p|` under conditioning on `zeta` refined from 10
   to 2560 strata.

Plus one control that the brief did not ask for and that turned out to matter: the
**residual-altitude test** of Appendix B of the Archimedean note was applied to a
target constructed to be an exact function of the allocation alone, so that the
correct answer is known to be "no independent information".

## The normalization audit (added after the first DNS pass)

The first pass of test D1 normalized the material rate by the local enstrophy,
`g = D_t e_omega/(||grad u||_F e_omega)`, because its production part is then exactly
`A sqrt(2(1-zeta))`. That statistic is singular as `e_omega -> 0`: its viscous part
diverges, `std(viscous)/std(production)` reaches `7e5`, and the Pearson correlation is
then set by a handful of near-irrotational points. `experiments/dns_checks/dns_reanalysis.py` repeats the
test on the same saved fields with `D_t e_omega/||grad u||_F^3` (no singularity, and
the production part is bounded by `4 sqrt 2/9`), reports Spearman alongside Pearson,
and reports the probability of net growth per `A`-quartile — a rank statistic that no
tail can dominate. Both the artifact and the corrected result are in `RESULTS.md`;
the first-pass number is reported, not quietly replaced.

## Statistical methods

Deliberately elementary (`src/statlib.py`), so every reported number is traceable:
OLS with HC3 robust errors; exact `O(n)` group-mean residualization for bin models
(identical to OLS on bin dummies); partial correlations and exact `Delta R^2`
increments from those residuals; quantile-binned mutual information with a
within-stratum permutation null; bootstrap CIs where used. No `statsmodels`, no
`sklearn`.

## Known limitations of the experiment

* The DNS is a single unforced Taylor–Green run at `Re = 400`: transitional, with a
  well-defined dissipation peak but no inertial range. Alignment statistics in
  developed high-`Re` turbulence (Ashurst et al. and successors) are not reproduced
  here and are cited, not re-measured.
* Restricted Euler omits the nonlocal pressure Hessian; its quantitative
  amplification numbers are model results.
* The Gaussian-field ensemble is a *null*, not a model of turbulence — by
  construction it has zero mean production.
* Nothing here tests the Archimedean framework's own empirical domain (the
  interplanetary-shock analysis). The audit is confined to the strain/rotation
  partition of an incompressible velocity gradient.
