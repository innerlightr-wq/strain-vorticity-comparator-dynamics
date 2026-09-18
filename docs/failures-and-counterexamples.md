# FAILURES AND COUNTEREXAMPLES

> *Reproduced verbatim from the audit run that produced it; only file paths have
> been updated to this repository's layout. Statements of the form "nothing outside
> this directory was modified" refer to that original run.*

What did **not** survive. Includes falsifications of the audit's own ideas, not only
of the hypothesis. Each entry names the artifact that produces it.

---

## F1 — The Thales altitude cannot distinguish rotation from strain — **EXACT**

Solid-body rotation (`zeta = +1`) and pure extensional strain (`zeta = -1`) both have

```
h = 0 ,   D = 1/2 ,   R = 1/(4h^2) = infinity ,   eta_T = 2h = 0 .
```

`h`, `D`, `R`, `eta_T` are **even** functions of `zeta`, hence two-to-one: they discard
`sign(zeta)`. The entire purpose of the `Q`-criterion / Omega method / Okubo–Weiss
lineage is that sign. *Artifact:* `src/exact_algebra.py`, `src/canonical_flows.py`.

This also falsifies outcome **B** (useful visualization): a coordinate that maps the
two opposite extremes of the state space onto the same point is not a better picture
of that space, it is a worse one.

## F2 — The Thales apex carries no production information — **EXACT**

Simple shear sits exactly at the apex (`zeta = 0`, `h = 1/2`, `D = 0`) and has
`P = 0`, `A = 0`. And at the same apex, three explicit states with identical
`||S||_F`, `|omega|`, `h`, `L`, `D` have `P = +4`, `0`, `-4`. Maximal Thales "coherence"
is compatible with stretching, nothing, or compression.
*Artifact:* `src/counterexamples.py` (C1'), `src/canonical_flows.py`.

## F3 — `zeta` determines nothing about production: the admissible set is a rectangle — **EXACT**

The map `(S, omega) -> (zeta, A)` is **onto** `(-1,1) x [-sqrt(2/3), sqrt(2/3)]`
(explicit inverse construction, verified on a 41×41 grid to `3.6e-16`). There is no
forbidden region, no conditional bound, no symmetry relation coupling allocation to
alignment.

This is the audit's own search for structure failing: Phase 8 explicitly looked for
forbidden regions and conditional bounds in the `(zeta, A)` plane. **There are none** —
and that is exactly why `A` is not redundant. *Artifact:* `src/counterexamples.py` (C3).

Numerically, in every ensemble, conditioning on `zeta` removes **0.00–0.01 %** of the
variance of `A` (E1 Gaussian field: 0.0019 %; E3 by construction: 0.0045 %), and 20
`zeta`-strata buy at most `+0.3 %` accuracy over the majority baseline for predicting
`sign(P)`. *Artifact:* `experiments/synthetic_fields/ensembles.py`.

## F4 — `zeta` sweeps its whole range in a flow with identically zero production — **EXACT**

The unstretched Lamb–Oseen vortex has `P = 0` and `A = 0` exactly at every radius and
time, while `zeta(eta)` runs from `+1` in the core through `0` at `eta = 1.1141` to `-1`
outside. A `zeta`-based coherence reading labels the core of a purely diffusing
vortex "rotationally favourable" and its skirt "strain-dominated", with production
zero in both. *Artifact:* `src/canonical_flows.py`.

## F5 — `(zeta, A)` without scale is also insufficient: the Burgers family — **EXACT**

At fixed similarity radius `eta = r/sqrt(4 nu/a)`, every velocity-gradient component
of the Burgers vortex scales linearly in the axial strain `a`, so

```
zeta(eta) and A(eta) are exactly independent of a  ,  while  P ∝ a^3 .
```

Over `a = 10^-2 -> 10^2` at `eta = 1` (`Re_Gamma = 100`): `zeta ≡ +0.0718`,
`A ≡ +0.5191`, `P = 8.6e-6 -> 8.6e+6`. Twelve decades of production at a fixed point
of the `(zeta, A)` plane.

This **refines rather than contradicts** Paper 2's Burgers counterexample: the
reported fall of `zeta` with increasing `a` is real at fixed *physical* radius, but it
is a core-radius effect (the core shrinks as `sqrt(4 nu/a)`), not a dynamical trend.
*Artifact:* `src/canonical_flows.py`.

## F6 — The geometric mean is the wrong allocation factor — **EXACT**

Exactly, `P = ||grad u||_F^3 · 2 b sqrt(a) · A`. The allocation monomial is
`a^{1/2} b^1`, so

```
g(zeta) = h · sqrt(2(1 + zeta)) ,
```

and `g/h` is non-constant: no rescaling of the Thales altitude is the allocation
factor. The partition enters enstrophy production through a **mixed-power,
asymmetric** combination that is neither the symmetric geometric mean `(ab)^{1/2}`
(the "amplitude" case of the Archimedean note's Appendix B) nor linear in `a` (its
"additive" case). The two-case classification does not contain the case that actually
occurs here. *Artifact:* `src/exact_algebra.py`.

Consistently with this, the amplitude cross term for this partition, `S:Omega`, is
**identically zero** (Prop. 1.1) — the decomposition that makes the partition
Pythagorean is the one that annihilates the interference term the altitude would
have measured.

## F7 — The residual-altitude test of Appendix B is a nonlinearity detector — **COMPUTATIONAL**

Appendix B of Paper 1 proposes: fit `h` linearly on `a`, take the residual
`h_resid`, and ask whether `h_resid` predicts the observable beyond `a` (criteria
`p < 0.01`, `Delta R^2 > 0.01`, bootstrap CI excluding zero).

We applied that exact procedure to a target **constructed to be an exact function of
the allocation alone** (`target = g(zeta) + independent noise`), where the true
independent content of `h` is provably zero. With the prescribed *linear* control in
`a`, the test reports

```
Delta R^2 = 0.5613 ,   t = 1094 ,   p < 1e-300     -> "h carries independent information"
```

which is false by construction. Refining the control in `a` collapses it:

| control in `a` | linear | 5 bins | 20 | 40 | 80 | 160 | 320 | 640 |
|---|---|---|---|---|---|---|---|---|
| `Delta R^2` | 0.561 | 0.143 | 1.2e-2 | 3.2e-3 | 8.3e-4 | 2.2e-4 | 6.1e-5 | 1.5e-5 |

The statistic measures *curvature of the target in `a`*, which `h` — being a
nonlinear function of `a` — inevitably proxies. Any target that is nonlinear in the
partition fraction will "pass".

Being precise about what this does and does not overturn:

* It threatens **positive** verdicts, not null ones. Appendix B's Case II (shocks) is
  a null result, and a nonlinearity detector that finds nothing is, if anything,
  stronger evidence of no `h`-dependence — that conclusion stands.
* Its Case I (two-slit) is a positive verdict, but there the observable *is* the
  geometric mean (`V = 2 sqrt(p_1 p_2)`), so the conclusion is true by construction
  and needs no test at all.
* What does not survive is the **general criterion**: Appendix B closes by offering
  the residual-altitude test as "an empirical method for distinguishing these cases
  in new systems", and for that purpose it is unsound — on a new system it will
  report independent altitude information whenever the observable is merely nonlinear
  in the partition fraction, which is the common case.

The fix is cheap: control for a flexible function of `a` (bins or splines) rather than
a linear one, or compare `h` against a genuinely independent regressor.
*Artifact:* `experiments/synthetic_fields/ensembles.py`.

## F8 — A trap the audit walked into and had to climb out of — **COMPUTATIONAL**

Our own first pass reported `I(p ; h | zeta) = 0.03–0.06` nats with `z`-scores of
several hundred against a shuffled null — apparently "`h` adds information after
conditioning on `zeta`". It does not: `h` is an exact function of `zeta`, so the true
value is zero. Two artifacts were responsible, and both are diagnosable:

* **binning leakage.** Conditioning on finite `zeta` strata leaves `h` able to resolve
  `zeta` *within* a stratum. Refining the conditioning kills it:
  partial `corr(h, |p|)` = `+0.096 -> +0.036 -> +0.0133 -> +0.0046 -> +0.0017` for
  `10, 40, 160, 640, 2560` strata (`~1/n_strata`), while partial `corr(|A|, |p|)`
  holds at `+0.948 -> +0.954`. Bias-corrected CMI: `I(p;h|zeta)` falls
  `0.049 -> 0.00017` nats from 8 to 512 strata while `I(p;A|zeta)` holds at `~2.0`.
* **plug-in MI bias.** With many strata the naive estimator's own floor (of order
  cells/2n per stratum) *rises* and can masquerade as signal; subtracting the
  within-stratum-shuffled null removes it.

Reported here in full because the same two artifacts would manufacture exactly the
"nontrivial conditional structure" that outcome **C** describes. *Artifact:*
`experiments/synthetic_fields/ensembles.py`, figure `fig7_leakage_scan.png`.

## F9 — A near-coincidence that is not an identity — **EXACT**

The derived production-optimal allocation is `b = 2/3` exactly. The
constraint-geometry framework's cubic landmark is the root of `b^3 + b - 1 = 0`,
`b_* = 0.6823278038...`. Substituting `b = 2/3` into the cubic leaves exactly
`-1/27`; the relative gap is `2.30 %`. They are different numbers.

No attempt was made to reconcile them, and none should be: `0.667` vs `0.682` is
precisely the kind of proximity that invites a retrospective fit. *Artifact:*
`src/exact_algebra.py`.

## F10 — `(scale, zeta, A)` is sufficient for `P` but **not** for the dynamics — **COMPUTATIONAL**

Two states with **identical** `||grad u||_F`, `zeta`, `h`, `L`, `D`, `A` and `P` but
opposite strain state (`s = -1` vs `s = +1`) diverge under restricted-Euler
evolution: `ln` enstrophy amplification `+0.2841` vs `+0.2556` over `t = 0.8`
(a `0.0285`-nat gap). Across the whole matched-cell design the amplification variance
left unexplained by `(zeta, A)` jointly is `0.44 %` (median within-cell spread
`0.036` nats against a total range of `1.22`).

So the sufficiency proved in `DERIVATIONS.md` §4 is sufficiency for *instantaneous*
production only — small but nonzero dynamical information lives in the two invariants
that `A` discards, exactly as Paper 2's Remark 4.4 warns. *Artifact:*
`experiments/synthetic_fields/restricted_euler.py` (T1, T3).

## F11 — Mean production is invisible to `zeta`, exactly — **EXACT**

For isotropically oriented vorticity, `E[cos^2 theta_i] = 1/3`, so
`E[A | S, |omega|] = (sum lambda_i)/(3||S||_F) = 0` for any traceless `S`, and `<P> = 0`.
A Gaussian solenoidal field confirms it numerically
(`<P/||grad u||_F^3> = +4.6e-5`, i.e. zero at `10^6` samples). Nonzero mean production
is therefore *entirely* an alignment-correlation effect — a property of the dynamics,
not of the kinematic partition — and `zeta` is invariant under the rotations that
generate the null. *Artifact:* `src/exact_algebra.py`, `experiments/synthetic_fields/ensembles.py` (9.5).

## F12 — No threshold, no corridor — **EXACT (negative)**

The production envelope `|P|/||grad u||_F^3 <= sqrt(2/3) g(zeta)` is a *bound*, not a
threshold: nothing in it separates dynamical regimes at a critical value of `zeta`,
`A`, `h`, or `D`. We searched (Phase 8) for invariants, forbidden regions, symmetry
relations and bifurcation-like structure in `(zeta, A)` and found a rectangle with a
smooth envelope. This reproduces Paper 2's own "no corridor found" conclusion and
extends it to the joint coordinate.

## F13 — Not tested, therefore not claimed

* Nothing here bears on the Archimedean framework's other domains (the
  interplanetary-shock analysis, the two-slit case, the `1/e` operating point): the
  audit is confined to the strain/rotation partition of an incompressible velocity
  gradient. F6 and F7 bear on the *method* used there, not on its conclusions.
* The alignment statistics of developed high-Reynolds turbulence (preferential
  alignment with the intermediate eigenvector) are cited, not re-measured; the DNS
  here is a single transitional `Re = 400` Taylor–Green run.
* The sharp constant `4 sqrt 2/9` is elementary. We did not find it stated in the
  literature we searched, but we make **no priority claim** — a bound of this kind is
  exactly what a careful reader of Wolkowicz–Styan plus the Pythagorean split would
  write down, and it may well be known.

---

## The strongest counterexample found

**F5 (Burgers vortex in self-similar coordinates).** It is exact, it lives in the
canonical flow that Paper 2 already identifies as decisive, and it is the only one
that defeats *both* candidate coordinate systems at once: the Thales/`zeta`
allocation layer **and** the joint `(zeta, A)` pair are simultaneously **constant**
along a family in which enstrophy production varies by twelve orders of magnitude.
Whatever else is true, no coordinate pair built from normalized allocation and
alignment can be a production diagnostic without the scale coordinate — and once
scale, allocation and alignment are all admitted, the triple is exactly the
factorization of `P` itself, which is to say: the established interaction term,
rewritten.
