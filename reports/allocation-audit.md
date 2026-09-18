# RESULTS

> *Reproduced verbatim from the audit run that produced it; only file paths have
> been updated to this repository's layout. Statements of the form "nothing outside
> this directory was modified" refer to that original run. Two later additions are
> marked as such and dated: bracketed editorial notes, and a prior-work addendum at the
> end of the file.*

# Classification: **A — pure reparameterization**

> **A. Pure reparameterization.** *The Thales geometry contributes no information
> beyond `zeta`, and the useful dynamics are entirely captured by established
> strain–vorticity interaction quantities.*

Both clauses are established here, the first **exactly** and the second by
reproducing the standard interaction term as the unique object that carries the
extra content. Outcome **B** was considered and rejected on evidence: the Thales
altitude is not merely redundant but **strictly lossy** (it is two-to-one in `zeta`,
mapping solid-body rotation and pure strain onto the same point), so it cannot be
defended as a better state-space *picture* either. Outcomes **C** and **D** were not
reached. Three candidate routes to **C** were examined and all three are reported:
two were artifacts that died under controls (the residual-altitude test, F7; our own
binning-leakage result, F8), and the third is a real but small statistical coupling
between allocation and alignment in the DNS (§2.6) which we argue explicitly does not
meet the bar. **D** would need a genuinely new invariant or theorem; what emerged
instead is a short list of exact but elementary sharpenings (§5), each of which is
one line from material already in the two papers.

The central hypothesis itself — *does `(zeta, A)` carry materially more dynamical
information than `zeta` alone* — is **true, and provably so** (`zeta` and `A` are
functionally independent; `zeta` determines neither the sign nor the magnitude of
production). But the additional content is exactly the content of the classical
vortex-stretching term `P = omega . S omega`, which Paper 2 already puts at the
centre. The allocation half of the hypothesis is the half that fails.

---

## 1. What was mathematically proved — **EXACT**

All 44 checks in `src/exact_algebra.py` pass on symbolic input; the witnesses in
`src/counterexamples.py` are exact matrices.

**1.1 The Thales layer is a reparameterization, and a lossy one.**
`a = (1-zeta)/2`, `b = (1+zeta)/2`, `L = zeta/2`, `h = (1/2)sqrt(1-zeta^2)`,
`D = (1-sqrt(1-zeta^2))/2`, `R = 1/(1-zeta^2)`. The lateral coordinate is *equivalent*
to `zeta`; the altitude, deficit, response and efficiency are **even in `zeta`**, hence
two-to-one, discarding `sign(zeta)` — the rotation- vs strain-dominated distinction
itself. `L^2 + h^2 = 1/4` and `L^2 = D(1-D)` are the same identity written twice.

**1.2 The exact three-factor decomposition the brief asked for.**

```
P = ||grad u||_F^3 · g(zeta) · A ,     g(zeta) = (1+zeta) sqrt((1-zeta)/2) = 2 b sqrt(a)
```

`production = scale × allocation × alignment`, exactly, with no residual.

**1.3 The allocation factor is not the Thales altitude.** `g = h sqrt(2(1+zeta))`, and
`g/h` is non-constant. The partition enters production as the asymmetric monomial
`a^{1/2} b`, which is neither the symmetric geometric mean `(ab)^{1/2}` nor linear in
`a` — i.e. neither case of the amplitude/additive classification in Appendix B of
Paper 1. Relatedly, the amplitude cross term of this partition, `S:Omega`, is
**identically zero**: the decomposition that makes the partition Pythagorean is the
one that annihilates the interference term an altitude would measure.

**1.4 A sharp global pointwise bound.**

```
omega . S omega  <=  (4 sqrt 2 / 9) ||grad u||_F^3  =  0.628539361... ||grad u||_F^3
```

with equality iff `||Omega||_F^2 = 2||S||_F^2` (`zeta = 1/3`), `S` axisymmetric with
eigenvalues `||S||_F(2,-1,-1)/sqrt 6`, and `omega` along the distinct eigenvector.
Elementary, and **no priority is claimed** — the `sqrt(2/3)` ingredient is classical
(Wolkowicz–Styan; RVP Thm 4.1) and the allocation optimization is one line; we did
not locate the combined constant in the literature searched, which is weak evidence
either way.

**1.5 A derived allocation landmark, and a warning about it.** The production-optimal
allocation is `zeta = 1/3`, i.e. `||Omega||_F^2 = 2||S||_F^2` — **not** the Thales apex.
But the *relative* (per-enstrophy) law is
`P/(||grad u||_F e_omega) = A sqrt(2(1-zeta))`, monotone decreasing, maximal as
`zeta -> -1`. Two exact laws, two different optima, neither at the apex: any claim
that a partition landmark is dynamically optimal must first name its normalization.

**1.6 Complete non-identifiability.** The map `(S, omega) -> (zeta, A)` is **onto the
rectangle** `(-1,1) x [-sqrt(2/3), sqrt(2/3)]` (explicit inverse construction, verified
to `3.6e-16`). Hence no conditional bound, no forbidden region, no symmetry relation
links allocation to alignment, and `zeta` carries *exactly zero* information about the
sign of `P`. Exact witnesses: three states with identical
`||S||_F, |omega|, zeta, h, L, D` and `P = +2, 0, -2`; the same at the apex with
`P = +4, 0, -4`.

**1.7 Invariant accounting.** `(S, omega)` has 8 degrees of freedom, 5 modulo `SO(3)`:
`(||S||_F, s, |omega|, cos^2_1, cos^2_2)`. `zeta` is a function of 2 of them (a ratio);
`A` collapses 3 into 1; `(||grad u||_F, zeta, A)` is exactly **sufficient for `P`** and
each of the three is **necessary**; 2 invariants are discarded.

**1.8 An exact null for mean production.** For isotropically oriented vorticity
`E[A | S, |omega|] = (sum lambda_i)/(3||S||_F) = 0`, so `<P> = 0` for any traceless `S`.
Nonzero mean production is therefore *entirely* an alignment-correlation effect, and
`zeta` is invariant under the rotations generating this null.

**1.9 Apex occupancy is forced by homogeneity.** `||S||_F^2 - ||Omega||_F^2 = tr((grad u)^2)`
exactly, and for any homogeneous incompressible field `<tr((grad u)^2)> = 0`. So

```
<||S||_F^2> = <||Omega||_F^2>   =>   <a> = <b> = 1/2 :   the mean allocation of EVERY
homogeneous incompressible flow sits exactly at the Thales apex.
```

Confirmed numerically in a Gaussian solenoidal field (ratio `1.000000`,
`<tr(grad u^2)>/<||grad u||^2> = 7e-20`) and in the DNS at all four snapshots
(ratio-of-means `a = 0.498–0.504`), while the *pointwise* mean of `zeta` is **not** zero
(`+0.027` to `-0.187`). In Paper 1's own vocabulary this is **occupancy without
meaning**: finding a turbulent flow "at the apex on average" is a kinematic identity,
not evidence of a selection principle.

## 2. What was computationally demonstrated — **COMPUTATIONAL / EMPIRICAL**

**2.1 Bounds hold, sharply, in 2.25M random states** (`experiments/synthetic_fields/ensembles.py`). Across four
ensembles: the largest `|A|` observed is `0.815131–0.816176` against the exact ceiling
`0.8164966`, and the largest `|P|/||grad u||_F^3` is `0.623705–0.628216` against the
exact ceiling `0.6285394` — approached but never exceeded. The `zeta`-dependent
envelope `|p| <= sqrt(2/3) g(zeta)` is respected at every one of the 2.25M samples.

**2.2 At fixed `zeta`, essentially all of the variation in `A` and `P` remains.**
Conditioning on 20 `zeta`-strata removes **0.00–0.01 %** of the variance of `A` in every
ensemble, including the physically generated Gaussian solenoidal field (0.0019 %).
Predicting `sign(P)` from 20 `zeta`-strata beats the majority baseline by at most
**+0.3 %**. (`sign(P) = sign(A)` is an identity and is not reported as a finding.)

**2.3 The Thales coordinates add nothing, measured three ways** (`experiments/synthetic_fields/ensembles.py`):
`R^2(h | exact formula in zeta) = 1.0000000000`; bias-corrected
`I(p ; h | zeta) -> 0.00017` nats as the conditioning is refined to 512 strata while
`I(p ; A | zeta)` holds at `~2.0` nats; partial `corr(h, |p|)` decays as `1/n_strata`
(`+0.096 -> +0.0017` for 10 -> 2560 strata) while partial `corr(|A|, |p|)` holds at
`+0.95`.

**2.4 Restricted Euler: alignment controls finite-time amplification at matched
allocation and scale.** 7280 trajectories, residual invariants randomized inside each
`(zeta, A)` cell. `R^2` of `ln[e_omega(T)/e_omega(0)]` on a saturated `zeta` model:
`0.167`; adding `A`: `0.935` (`Delta R^2 = +0.768`). The amplification is monotone in
`A` at every `zeta`. The measured `A`-sensitivity follows the **exact** law
`sqrt(2(1-zeta))` to `4.9 %` (coefficient of variation; a systematic `15 %` drift across
the `zeta` range remains, attributable to the finite integration window). Variance left
unexplained by `(zeta, A)` jointly: `0.44 %`.

**2.5 Navier–Stokes DNS** (`96^3`, `Re = 400`, unforced Taylor–Green, `k_max eta_K = 1.10`,
peak dissipation `0.01103` at `t = 8.76`; `experiments/dns_checks/dns_taylor_green.py`,
`experiments/dns_checks/dns_reanalysis.py`). The exact identity
`P/(||grad u||_F e_omega) = A sqrt(2(1-zeta))` holds pointwise to `4.3e-15`. With the
**full** balance `D_t e_omega = P - eps_omega + T_nu` evaluated spectrally:

* `P > 0` at `52.6–62.5 %` of points, but `P > 0` while `D_t e_omega < 0` at
  `12.0–30.2 %` — production is emphatically *not* net growth, exactly as the brief
  insists.
* Inside matched `(zeta, scale)` cells, the **probability of net enstrophy growth rises
  monotonically with `A` at every snapshot**:
  `9 % -> 34 % -> 57 % -> 74 %` (`t=4`), `15 -> 31 -> 47 -> 68` (`t=6`),
  `13 -> 28 -> 43 -> 64` (`t=8`), `12 -> 28 -> 46 -> 67` (`t=10`) across `A`-quartiles.
  (These are the 80 000-point saved subsamples; the solver's own full-field pass over
  all 884 736 points gives the same ordering to within one percentage point.)
* Rank (Spearman) partial correlation of `A` with the net material rate at matched
  `(zeta, scale)`: `+0.78, +0.50, +0.45, +0.46`.
* The viscous terms are **essentially uncorrelated with `A`** at matched `(zeta, scale)`
  (`|partial corr| < 0.04`): they do not cancel the alignment dependence, they add
  variance to it — and that variance grows from `0.8x` to `34x` the production
  variance between `t=4` and `t=10`, which is why the *Pearson* correlation with the
  net rate falls from `+0.62` to `+0.02` while the rank association persists. This is
  a variance effect, not a sign reversal, and it is the reason a normalization audit
  was necessary (see §4).
* **Lagrangian test** (40 000 particles, `t = 7.0 -> 8.5`): at matched `(zeta, scale)`,
  mean amplification by initial-`A` quartile is
  `-0.431, -0.199, -0.054, +0.118`; only the top `A`-quartile grows on average
  (`57.5 %` of its particles growing vs `23.9 %` in the bottom quartile); partial
  `corr(A, amplification | zeta, scale) = +0.298`, `Delta R^2 = 0.081`.

**2.6 Phase-8 finding: allocation and alignment are statistically coupled in a real
flow, though geometrically independent — EMPIRICAL.** `zeta` and `A` are functionally
independent (§1.6), so any dependence between them is a property of the flow's
measure. Measured against a permutation null: the Gaussian solenoidal field shows
`+0.005 %` of `var(A)` explained by 10 `zeta`-deciles (i.e. none), while the DNS shows
`+1.24 %` (`t=4`), `+1.24 %` (`t=6`), `+0.83 %` (`t=8`), `+0.52 %` (`t=10`), at
`z = +132` to `+207` above null. Conditional means of `A` move by `~0.16` between the
extreme `zeta` deciles at `t=4` (`-0.119 -> +0.038`).

This is real conditional structure and it is exactly what Phase 8 asked to look for.
It is **not** enough for outcome **C**, for three stated reasons: (i) it is a property
of this flow's measure, not of the proposed representation — the same coupling is
visible in `(Q, alignment-cosine)` statistics already standard in the literature;
(ii) it is `~1 %` of variance in a single transitional `Re = 400` run, and the
*direction* of the trend is not stable across snapshots (`t=8` and `t=10` are much
weaker and not monotone in the same sense); (iii) it says nothing about the Thales
layer, which remains an exact function of `zeta`. Whether it survives at higher
Reynolds number is a well-posed follow-up (§7).

**2.7 Controls** (`experiments/synthetic_fields/ensembles.py` §9). Rotating `omega` isotropically leaves every
magnitude coordinate fixed to `5.6e-16` while `A` sweeps `[-0.408, +0.816]` and `P`
changes sign for `57.8 %` of orientations; scaling `S` and `omega` together fixes
`(zeta, A, p)` exactly while `P` moves by a factor `2.1e16`; randomizing strain
eigenvectors at fixed eigenvalues leaves `zeta` fixed and moves `A` over its full
range; mean `A` is `+0.0007` for isotropic orientation vs `+0.646` when tilted toward
the most extensional eigenvector; and a Gaussian solenoidal field has
`<P/||grad u||^3> = +4.6e-5`, i.e. zero.

## 3. What failed

Full list in `FAILURES_AND_COUNTEREXAMPLES.md`. The load-bearing failures:

| # | failure | label |
|---|---|---|
| F1 | `h`, `D`, `R` cannot separate solid-body rotation from pure strain (same point) | EXACT |
| F2 | the Thales apex carries no production information (`h = 1/2`, `P = 0` for simple shear; `P = +4,0,-4` at the apex) | EXACT |
| F3 | `zeta` fixes nothing about `P`: the `(zeta,A)` image is a rectangle. Our own Phase-8 search for forbidden regions and conditional bounds found none | EXACT |
| F4 | `zeta` sweeps its entire range inside the Lamb–Oseen vortex, where `P ≡ 0` | EXACT |
| F5 | in self-similar Burgers coordinates `(zeta, A)` are constant while `P ∝ a^3` (12 decades) — allocation *and* alignment without scale is also insufficient | EXACT |
| F6 | the geometric mean is the wrong allocation factor (`a^{1/2}b`, not `(ab)^{1/2}`) | EXACT |
| F7 | the residual-altitude test of Paper 1's Appendix B is a **nonlinearity detector**: applied to a target with provably zero independent `h` content it reports `Delta R^2 = 0.56`, `t = 1094`, `p < 1e-300` | COMPUTATIONAL |
| F8 | our own first pass found "`h` adds information beyond `zeta`" — an artifact of finite-stratum conditioning plus plug-in MI bias, killed by refinement | COMPUTATIONAL |
| F9 | the derived optimum `b = 2/3` is **not** the framework's cubic landmark `b_* = 0.68233` (cubic residual exactly `-1/27`, gap `2.30 %`) | EXACT |
| F10 | `(scale, zeta, A)` is sufficient for `P` but not for its evolution: identical-coordinate states diverge (`0.0285` nats) | COMPUTATIONAL |
| F12 | no threshold and no corridor exists in `(zeta, A)` — reproduces Paper 2's own negative result | EXACT (negative) |

## 4. What was merely a reparameterization

* **Everything the Thales construction contributes to this partition.** `a, b, L` are
  affine in `zeta`; `h, D, R, eta_T` are functions of `|zeta|`. `zeta` in turn is an affine
  rescaling of the second invariant `Q` (and of the unregularized Omega method), as
  Paper 2 already states. So the chain
  `Q -> zeta -> (a,b) -> (h,L,D,R)` adds no content at any link, and loses `sign(zeta)`
  at the last one.
* **`P/(||S||_F |omega|^2) = A`.** Identically `A`; never reported here as a
  correlation. Likewise `sign(P) = sign(A)`, `R^2 = 1` for `P` on `(scale, zeta, A)`,
  and the partial correlation of `A` with the *production part* of the rate
  (`+0.86–0.97` in the DNS): all identities, listed as such.
* **One normalization choice that behaved like a result and was not.** The
  per-enstrophy rate `D_t e_omega/(||grad u||_F e_omega)` has its production part equal
  to `A sqrt(2(1-zeta))` (attractive) but its viscous part diverges as `e_omega -> 0`,
  giving `std(viscous)/std(production)` up to `7e5` and a near-zero Pearson
  correlation that is pure tail behaviour. Re-normalizing by `||grad u||_F^3` and using
  rank statistics recovers the association. The first-pass number was a normalization
  artifact, and is reported as one.

## 5. What genuinely new insight, if any, emerged

> *Editorial note added September 2026:* "new" here means **new relative to the two
> audited manuscripts and the earlier audits**, which is what this section was written to
> assess. For the comparison against the published literature — where most of these items
> are classical or reparameterizations — see the prior-work addendum at the end of this
> report and [`../docs/LITERATURE_REVIEW.md`](../docs/LITERATURE_REVIEW.md).

Nothing that overturns either paper. Five statements that are, as far as we can tell,
not in either manuscript, all elementary and all exact:

1. **The exact factorization** `P = ||grad u||_F^3 g(zeta) A` with
   `g(zeta) = 2 b sqrt(a)`, which is the precise form of the
   "scale × allocation × alignment" decomposition the brief hypothesized, and which
   shows the allocation enters as `a^{1/2} b` — a functional form outside Paper 1's
   amplitude/additive dichotomy.
2. **The sharp constant** `sup P/||grad u||_F^3 = 4 sqrt 2/9` with its complete equality
   case, and the derived allocation landmark `zeta = 1/3` (`||Omega||^2 = 2||S||^2`).
3. **The two-normalization theorem** (§1.5): absolute and relative optima are
   `zeta = 1/3` and `zeta -> -1`; the apex is neither.
4. **The Burgers vortex core saturates the sharp alignment bound exactly**
   (`A = sqrt(2/3)` on the axis, for every `a`, `Gamma`, `nu`), and in self-similar
   coordinates the whole `(zeta, A)` portrait of the Burgers family is invariant while
   `P ∝ a^3` — a sharper form of Paper 2's decisive counterexample, and one that
   identifies the trend it reports as a core-radius effect.
5. **Apex occupancy is a kinematic identity** for any homogeneous incompressible flow
   (§1.9) — which removes "the flow sits at the apex" from the list of things that
   could ever be evidence.

Methodologically, one finding is likely to be more useful to the author than any of
the above: **the residual-altitude test of Appendix B cannot serve as the general
criterion it is offered as** (F7). It is a nonlinearity detector, so it will report
independent altitude information for any observable that is merely nonlinear in the
partition fraction. Neither published verdict is overturned by this — the shock case
is a null, and the two-slit case is true by construction because the observable *is*
the geometric mean — but the test should not be applied to a new system as it stands.
The fix is cheap: control for a flexible function of `a`, not a linear one.

## 6. The strongest counterexample found

**The Burgers vortex in self-similar coordinates** (F5). At fixed
`eta = r/sqrt(4 nu/a)`, `zeta` and `A` are *exactly* independent of the axial strain `a`
while `P ∝ a^3`: over `a = 10^-2 -> 10^2` at `eta = 1`, `zeta ≡ +0.0718` and
`A ≡ +0.5191` while `P` runs `8.6e-6 -> 8.6e+6`. It is exact, it lives in the canonical
flow Paper 2 identifies as decisive, and it defeats both candidate coordinate systems
at once — the allocation layer *and* the joint `(zeta, A)` pair. The only way to
restore a production diagnostic is to admit the scale coordinate, at which point the
triple is the factorization of `P` itself: the established interaction term, rewritten.

Runner-up, and the sharpest strike against the Thales layer specifically: the three
exact states at the **apex** `zeta = 0, h = 1/2, D = 0` with `P = +4, 0, -4`.

## 7. Is further investigation justified?

**For the Thales/allocation layer applied to the strain–rotation partition: no.** The
question is closed by exact algebra, not by statistics — `h`, `D`, `R` are functions of
`|zeta|`, so no data set can ever give them independent content here, and the one
place where the partition does enter the physics (the factor `a^{1/2}b`) is not the
altitude. Further empirical work on this pairing would be measuring an identity.

**Two narrow follow-ups are well posed** (both **CONJECTURE** / open, neither claimed):

1. **Who sits near the sharp bound?** The bound `P <= (4 sqrt 2/9)||grad u||_F^3` is
   saturated exactly on the Burgers axis. Whether concentrating or
   blowup-candidate flows approach it — and in particular the sign and magnitude of
   `sigma_eff` along the core of the forced-blowup construction that Paper 2 §8.2
   explicitly leaves open — is a sharply posed, answerable question, and the
   saturation deficit `sqrt(2/3) g(zeta) - P/||grad u||^3` is a concrete scalar to test
   it with. This is a question about the *interaction* geometry, not about the
   allocation geometry.
2. **Does the DNS's `zeta`–`A` coupling survive at higher Reynolds number?** §2.6
   found `0.5–1.24 %` of `var(A)` explained by `zeta` in a transitional run, with an
   unstable sign — either a genuine feature of developed turbulence's alignment
   statistics (in which case it belongs to the Ashurst programme and should be
   reported as alignment cosines conditioned on `Q`, which is how that literature
   already frames it) or a transitional artifact. One higher-`Re`, larger-box run
   would decide it.
3. **The two discarded invariants.** `(scale, zeta, A)` is sufficient for `P` but leaves
   `0.44 %` of restricted-Euler amplification variance unexplained, and in the DNS the
   viscous terms — uncorrelated with `A` — dominate the variance of the net rate at
   late times. The natural next coordinates are the three alignment cosines
   themselves, i.e. the programme Ashurst et al. (1987) opened; that is a
   continuation of the existing literature rather than a new framework.

For the author's purposes, the actionable items are the Appendix-B methodological
correction (F7), the retraction-strength fact that apex occupancy is kinematically
forced (§1.9), and the note that `2/3` is not `0.68233` (F9).

---

## Appendix: verification status of every artifact

| script | checks | status |
|---|---|---|
| `src/exact_algebra.py` | 44 symbolic | all pass |
| `src/counterexamples.py` | 14 exact/constructive | all pass |
| `src/canonical_flows.py` | 21 symbolic | all pass |
| `experiments/synthetic_fields/ensembles.py` | 19 (ensembles, Thales test, 6 controls) | all pass |
| `experiments/synthetic_fields/restricted_euler.py` | 7 | all pass |
| `experiments/dns_checks/dns_taylor_green.py` | 5 | all pass |
| `experiments/dns_checks/dns_reanalysis.py` | 5 | all pass |
| `experiments/dns_checks/homogeneity_check.py` | 3 | all pass |

Each script exits non-zero if any of its own checks fail. Raw output is in
`results/*.json`; sampled fields in `results/*.npz`; figures in `figures/`.

---

## Prior-work addendum (added September 2026, after a literature review)

Added after the report above, which is otherwise unchanged. See
[`../docs/LITERATURE_REVIEW.md`](../docs/LITERATURE_REVIEW.md) for the full audit and
[`../references.bib`](../references.bib) for verified records.

* `ζ` is an affine rescaling of the second invariant `Q` (Hunt, Wray & Moin 1988) and is
  **exactly** `2Ω − 1` for the Ω vortex-identification measure of Liu et al. (2016); in
  two dimensions it reparameterizes the Okubo–Weiss discriminant (Okubo 1970;
  Weiss 1991). Verdict **A** is therefore reinforced, not weakened: the coordinate the
  Thales layer reduces to was already published.
* The independence of alignment from the magnitude coordinates — §2 of this report — is
  classical, the founding observation of the alignment literature (Ashurst et al. 1987;
  reviews: Meneveau 2011; Johnson & Wilczek 2024). The exact witnesses here demonstrate
  it inside these coordinates; they do not establish it.
* The sharp constant `4√2/9` was not located in the reviewed literature, but it is an
  elementary optimization over two classical ingredients (the alignment bound of
  Wolkowicz & Styan 1980 and the normalization identity). **No priority is claimed**, and
  a prior appearance would not be surprising.
* `⟨E_S⟩ = ⟨E_W⟩` for homogeneous flow, used in §2, is a classical homogeneity identity
  (Betchov 1956), now completely classified by Carbone & Wilczek (2022).

No result in the report is contradicted by the reviewed literature, and no number
changed.
