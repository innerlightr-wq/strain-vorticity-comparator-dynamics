# The material evolution of the comparator ξ

> *Reproduced verbatim from the audit run that produced it; only file paths have
> been updated to this repository's layout. Statements of the form "nothing outside
> this directory was modified" refer to that original run.*

Continuation of `reports/xi-evolution-audit.md` (classification **B**) and of the parent audit in
`reports/allocation-audit.md` (classification **A**). Nothing outside this directory was created or
modified; the parent's `svcore.py`, `statlib.py`, `generators.py` and
`dns_taylor_green.py` are imported read-only.

**Verification status.** `src/xi_exact.py` — 49/49 symbolic checks;
`src/xi_numeric.py` — 13/13 numerical checks (200 000 random admissible states);
`experiments/dns_checks/xi_dns_pressure.py` — 8/8 checks on a 96³ Navier–Stokes field (883 592 points).

**Central answer, up front.** ξ's evolution *does* isolate the local nonlinear
dynamics from the nonlocal pressure obstruction, exactly and completely: the entire
pressure coupling reduces to one scalar, and the locally-computable part of the
pressure Hessian cancels identically. But it does not *control* that scalar, and no
combination of the framework's quantities can remove it. **Classification: B.**

---

## 1. Conventions (re-locked, not re-chosen)

Taken verbatim from `src/svcore.py` and `docs/comparator-results.md`, and re-verified
symbolically in §A of `src/xi_exact.py`:

```
A := grad u ,  (grad u)_ij = d u_i/d x_j        S = (A+A^T)/2     W = (A-A^T)/2
omega_i = eps_ijk d_j u_k        =>        W_ij = -(1/2) eps_ijk omega_k
||W||_F^2 = |omega|^2 / 2         <-- PROJECT CONVENTION, re-verified
E_S := ||S||_F^2    E_W := ||W||_F^2 = |omega|^2/2 = enstrophy density
Q := ||grad u||_F^2 = E_S + E_W          a = E_S/Q   b = E_W/Q   a + b = 1
zeta = b - a        h = sqrt(ab)         xi = artanh(zeta)
P = omega.S omega        A_align = P/(||S||_F |omega|^2)        |A_align| <= sqrt(2/3)
s = -3 sqrt6 det S / ||S||_F^3      (Lund-Rogers strain state, |s| <= 1)
momentum:  Du/Dt = -grad p + nu Lap u        H := Hess(p),  symmetric
```

## 2. Exact definition of ξ

```
xi = artanh(zeta) = (1/2) log(E_W / E_S)
```

verified equal (the shared normalization `Q` cancels). ξ is a pure *log-ratio of the
two sector energies*; this is the representation the whole derivation uses.

## 3. Full derivation

Differentiating the momentum equation, `D A/Dt = -A^2 - H + nu Lap A`. Its trace is
the pressure Poisson equation. Splitting by symmetry (`sym(A²) = S² + W²`,
`skew(A²) = SW + WS`):

```
DS/Dt = -(S^2 + W^2) - H   + nu Lap S
DW/Dt = -(SW + WS)         + nu Lap W          <-- no H: Hess(p) is symmetric
```

Contracting and using the derived (not assumed) coefficients
`S:W² = ¼[ω·Sω − |ω|² tr S]` and `W:(SW+WS) = −2 S:W²`:

```
D E_W/Dt = P + 2 nu W:Lap W
D E_S/Dt = -2 tr(S^3) - (1/2) P - 2 S:H_dev + 2 nu S:Lap S
```

Normalizing (`P/E_W = 2‖S‖A`, `−½P/E_S = −‖S‖A e^{2ξ}`,
`−2tr(S³)/E_S = +2s‖S‖/√6`) and taking half the difference of the logarithmic rates:

## 4. The exact Navier–Stokes evolution equation

```
                ┌──────────────── LOCAL ────────────────┐   ┌─── PRESSURE ───┐
D xi/Dt  =  ||S||_F [ (1 + e^{2 xi}/2) A  -  s/sqrt(6) ]  +  S:H_dev/||S||_F^2
                                                          +  nu[ W:Lap W/E_W - S:Lap S/E_S ]
                                                             └───── VISCOUS ─────┘
```

Equivalently, in `(‖∇u‖_F, ζ)` form,

```
D xi/Dt|_local = ||grad u||_F [ A (3-zeta) / (2 sqrt(2(1-zeta)))
                                - (s/sqrt6) sqrt((1-zeta)/2) ]
```

since `‖S‖_F = ‖∇u‖_F √a` and `(1+a)/(2a) = (3−ζ)/(2(1−ζ)) = 1 + e^{2ξ}/2`.

## 5. Local / pressure / viscous decomposition

| group | term | what it needs |
|---|---|---|
| vortex stretching | `‖S‖ A` (from `P/E_W`) **and** `‖S‖ A e^{2ξ}/2` (from `−P/2E_S`) | `‖S‖, ξ, A` |
| strain self-amplification | `−‖S‖ s/√6` (from `−2 tr S³`) | `‖S‖, s` |
| allocation weighting | the factor `1 + e^{2ξ}/2 = 1 + b/(2a)` | `ξ` |
| **pressure** | `S:H_dev/E_S` | the **nonlocal** Poisson solution |
| viscous | `ν[W:ΔW/E_W − S:ΔS/E_S]` | **second** derivatives of `∇u` |

The same production term `P` enters **both** sectors — with `+1` in the rotation
sector and `−1/2` in the strain sector — which is why vortex stretching appears twice
in `Dξ/Dt` with the combined weight `1 + e^{2ξ}/2`. Stretching moves the state toward
rotation dominance twice over: it grows enstrophy *and* drains strain energy.

**Encoded by the parent's coordinates:** scale, allocation, comparator and alignment
capture the stretching group entirely. **Not encoded:** the strain-state invariant
`s` (one further local invariant), the viscous group (second derivatives), and the
pressure scalar (nonlocal).

## 6. Pressure-Hessian audit — **EXACT**

1. **The vorticity sector is exactly pressure-blind.** `H = Hess(p)` is symmetric, so
   its antisymmetric part is identically zero and it cannot appear in `DW/Dt`.
   (Equivalently: `curl grad p = 0`.) Verified numerically by perturbing `H`
   arbitrarily and observing `D E_W/Dt` unchanged to `1e-10`.
2. **Only the deviatoric part acts.** `S:(tr H/3)I = (tr H/3) tr S = 0` for
   incompressible flow, so `S:H = S:H_dev` exactly.
3. **The part that cancels is exactly the part that is locally known.**
   `tr H = Δp = E_W − E_S = Q ζ = Q tanh ξ` — a pointwise function of the gradient.
   So `ξ` discards precisely the computable part of the pressure Hessian and retains
   precisely the part that requires solving the nonlocal Poisson problem. This is the
   sharpest positive structural result of the audit.
4. **The surviving scalar.** `S:H_dev/E_S = B ‖H_dev‖_F/‖S‖_F` with
   `B := S:H_dev/(‖S‖_F‖H_dev‖_F) ∈ [−1,1]`.
5. **Contrast.** The same exercise on `P` is *not* clean: `DP/Dt` contains
   `−ω·Hω = −(ωω)_dev:H_dev − (tr H/3)|ω|²`, and the second piece does **not** vanish.
   ξ is special in this respect; `P` is not.

**The pressure source is the comparator itself.** `Δp = Q tanh ξ`. The coordinate
whose evolution is being studied is, multiplied by the scale, exactly the source of
the field that obstructs its own closure.

## 7. Euler limit

Setting `ν = 0` removes the viscous pair exactly:

```
D xi/Dt|_Euler = ||S||_F[(1 + e^{2xi}/2) A - s/sqrt6] + S:H_dev/||S||_F^2
```

This is a genuine simplification — the viscous group is the only one requiring
*second* derivatives, so in Euler the equation involves only the gradient tensor plus
one nonlocal scalar. **The pressure Hessian is untouched by the inviscid limit and
remains the sole obstruction.** No regularity claim is made or implied.

## 8. Restricted-Euler comparison

```
FULL NS:            D xi/Dt = ||S||[(1 + e^{2xi}/2)A - s/sqrt6]
                              + S:H_dev/||S||^2
                              + nu[W:Lap W/E_W - S:Lap S/E_S]

RESTRICTED EULER:   D xi/Dt = ||S||[(1 + e^{2xi}/2)A - s/sqrt6]

DIFFERENCE:         S:H_dev/||S||^2      (+ the viscous pair)
```

Restricted Euler replaces `H` by its isotropic part `(tr H/3)I`, which contracts to
zero against `S`. **The entire nonlocal obstruction is exactly the term restricted
Euler discards** — no more and no less.

**Consequence for the parent audit.** In restricted Euler, `Dξ/Dt` closes *exactly* on
`(‖S‖, ξ, A, s)`. That retro-explains the parent audit's finding that `0.44 %` of the
restricted-Euler amplification variance was unexplained by `(ζ, A)`: the missing
variable is `s`, and in that model there is nothing else. Every apparent simplicity
of the parent's restricted-Euler results is attributable to suppressing `H_dev`.

## 9. Closure test — **PARTIAL (answer: NO for the parent's coordinate set)**

Two independent failures, both exhibited exactly (`src/xi_numeric.py` M5, M6):

**(a) The local part already fails to close on `(Q, h, ξ, A)`.** With
`ζ = 0.2`, `A = 0.3`, `‖∇u‖ = 1` and `H_dev = 0`:

| state | `ζ` | `h` | `ξ` | `A` | `s` | `Dξ/Dt` |
|---|---|---|---|---|---|---|
| I | +0.2 | 0.489898 | +0.202733 | +0.3 | **−1** | **+0.590238** |
| II | +0.2 | 0.489898 | +0.202733 | +0.3 | **+1** | **+0.073840** |

Identical in every coordinate the framework carried, an eightfold difference in the
comparator's drift, with the pressure Hessian identically zero. The missing variable
is the strain state `s`. Adding it *does* close the local part.

**(b) With `(Q, h, ξ, A, s)` all fixed, the pressure still moves it.** Same state,
two environments: `H_dev = 0` gives `Dξ/Dt = +0.5902`; a unit misaligned `H_dev` gives
`+1.2030`. The local parts are identical to machine precision; the difference is
entirely `S:H_dev/E_S`.

**Verdict.** `Dξ/Dt` closes on `(‖S‖, ξ, A, s)` **plus one nonlocal scalar** (plus
second derivatives when `ν > 0`). It does not close on the gradient tensor alone.

## 10. Pressure-Hessian geometry

`S:H_dev` is a contraction of two traceless symmetric tensors, so it depends on
**both** spectra and their **relative eigenframe orientation** — it is a full
5-dimensional inner product, not a smaller invariant. Writing
`S:H_dev = ‖S‖‖H_dev‖ B`:

* `B` is structurally parallel to the alignment coordinate `A`: both are normalized
  direction cosines, `|A| ≤ √(2/3)`, `|B| ≤ 1`.
* The parallel **breaks** on the magnitude factor. `A` multiplies `‖S‖`, a local
  quantity bounded by the state; `B` multiplies `‖H_dev‖/‖S‖`, which is **nonlocal and
  unbounded**. A large `‖H_dev‖` with `B = 0` contributes nothing; a modest
  `‖H_dev‖` perfectly aligned dominates. Magnitude alone does not control the term.
* Special states 6–8 (`SPECIAL_STATES.md`) make this concrete: identical local state,
  pressure contribution `0`, `+1.4142`, `−0.7071` for `B = 0, +1, −0.5`.

**We do not add `B` to the framework.** It is not a state variable of the fluid
element: it depends on the environment through the Poisson solve, and it is
independent of the local invariants (§13, `R² = 0.37` on a saturated local binning).
The honest hierarchy is therefore *not* `scale × allocation × comparator × alignment +
pressure coupling` as a closed list of state variables, but

```
    local state (scale, comparator, alignment, strain state)  ⊕  an external field
```

where the second summand is not a coordinate of the element at all.

## 11. Special states

Full table in `SPECIAL_STATES.md`. The two structural readings:

* **The fold `ζ = 0` is a regular point of the comparator equation.** `Dξ/Dt` is
  finite and generic there. The degeneracy that made `h` two-to-one and
  ill-conditioned is invisible to `ξ`.
* **The two endpoints are not symmetric.** `ζ → −1` (pure strain) is a *coordinate*
  singularity: `ξ → −∞` but `Dξ/Dt → ‖S‖[A − s/√6] + S:H_dev/E_S`, finite.
  `ζ → +1` (pure rotation) is a *genuine* divergence: `‖S‖e^{2ξ}A = (E_W/‖S‖)A → ∞`.
  The reason is structural — the rotation sector is forced only by `P`, which vanishes
  with `E_W`, whereas the strain sector is forced by `−P/2` and `−2tr(S³)`, which do
  not vanish as `E_S → 0`. A vanishing strain sector is regenerated at a finite rate
  from an infinitesimal base.

## 12. Conditioning

```
d xi/d zeta = 1/(1 - zeta^2) = R = 1/(4 h^2)
```

— the Jacobian of the comparator with respect to the allocation is exactly the
framework's own *response* coordinate. Consequences:

| `ζ` | 0 | 0.5 | 0.9 | 0.99 |
|---|---|---|---|---|
| `dξ/dζ` | 1.000 | 1.333 | 5.263 | 50.25 |
| `dh/dζ` | 0 (fold) | −0.289 | −1.032 | −3.509 |

`ξ` **regularizes the fold** (where the previous continuation found the branch
reconstruction ill-conditioned, amplification `4h/|ζ| → ∞`) and **stretches the
endpoints**. For dynamics this is the right trade: motion through balance is
well-resolved, and sector extinction is pushed to infinite coordinate distance — but
the drift itself diverges there (§11), so `ξ` does not tame the endpoints, it only
relabels them.

## 13. Attempted cancellation — **IMPOSSIBLE, with proof**

The pressure enters `D E_W/Dt` not at all and `D E_S/Dt` only as `−2 S:H_dev`.
Therefore, for **any** scalar `F(E_S, E_W)`:

```
   pressure contribution to DF/Dt  =  (dF/dE_S) · (-2 S:H_dev)
```

> **Theorem.** `DF/Dt` is free of the pressure Hessian ⟺ `∂F/∂E_S = 0` ⟺ `F` is a
> function of the **enstrophy alone**.

This disposes of the entire family the brief proposed, at once:
`ξ + α log h`, `log g(ζ)`, `log(2h^{3/2}e^{ξ/2})` and every other function of the two
sector energies. Explicitly, `h = 1/(2cosh ξ)` gives `D log h/Dt = −ζ Dξ/Dt`, so
`D/Dt[ξ + α log h] = (1 − αζ) Dξ/Dt`: the combination merely **rescales** the pressure
term by `(1 − αζ)`, which cannot vanish identically on `(−1,1)`. Any `F` that avoids
the pressure must ignore the strain sector — and then it is not a comparator.

**Enlarging the family does not help.** Admitting `F(E_S, E_W, P)` brings in `DP/Dt`,
which contains the *different* pressure functional `ω·H_dev ω = (ωω)_dev:H_dev`.
Cancellation would require `c₁S + c₂(ωω)_dev = 0` as traceless symmetric tensors;
these are generically linearly independent (verified). The exceptional set is
`S ∝ (ωω)_dev`, i.e. an axisymmetric `S` with `ω` on its distinct axis — which is
exactly the Burgers-vortex core state that saturates the parent audit's alignment
bound. A codimension-4 set, not a mechanism.

**What *does* cancel — and why it does not help here.** By parts on a periodic or
decaying domain, `∂_i∂_j S_ij = Δ(div u) = 0`, hence

```
   ∫ S:H_dev dx = ∫ p ∂_i∂_j S_ij dx = 0        (exactly)
```

The pressure term vanishes in the **unweighted volume integral** — which is the
classical reason the integrated strain-energy balance is pressure-free and Betchov's
relation exists. The obstruction in `Dξ/Dt` is therefore created *entirely by the
pointwise weight `1/E_S`*, not by the pressure term's own mean. Confirmed in the DNS
(§15, P6): the unweighted mean is zero to a relative `2.7e-18`, while the
`1/E_S`-weighted mean is `−0.0156`. Two caveats, stated so the point is not oversold:
that weighted mean is small compared with the term's typical magnitude
(`median |S:H_dev/E_S| = 0.28`), so the weighted term is mostly fluctuation with a
weak negative bias; and the vanishing of the unweighted integral is the classical
fact behind Betchov's relation, not a discovery of this audit.

But the unweighted route is empty for the comparator: for homogeneous incompressible
flow `∫E_S = ∫E_W` exactly (the parent audit's §1.9 identity), so the global
comparator `Ξ = ½log(∫E_W/∫E_S)` is **identically zero at all times** and its
evolution is trivially `0`. The useful objects therefore lie strictly between the
pointwise and the unweighted-global: weighted or conditional averages, where the
pressure term's fate is genuinely open. That is the recommended next test (§19).

## 14. Symbolic verification

`src/xi_exact.py`, 49/49 checks, all on symbolic `3×3` matrices with free entries —
no coefficient assumed. Covers: the convention re-lock; `S:W² = ¼[ω·Sω − |ω|²trS]`;
`W:(SW+WS) = −2S:W²`; `W²ω = 0`; `tr(S³) = 3det S = −s‖S‖³/√6`; both sector equations
coefficient-by-coefficient against the textbook forms; the Poisson trace; the
deviatoric-only property; the assembled law; both `ζ`- and `ξ`-forms of the
coefficients; the Euler and restricted-Euler reductions; the cancellation theorem and
its extension; `Δp = Q tanh ξ`; `dξ/dζ = R`; and the additive/multiplicative duality
`D log P/Dt = 3 D log‖∇u‖/Dt + (3/2) D log h/Dt + ½ Dξ/Dt + D log A/Dt`.

## 15. Numerical verification

**Random admissible states** (`src/xi_numeric.py`, 200 000 samples, random symmetric
`H` satisfying `tr H = −tr(A²)`, random `ΔS`, `ΔW`, `ν = 0.37`):

| check | result |
|---|---|
| full law vs. direct assembly from `DA/Dt` | max relative error `1.7e-11` |
| `D E_W/Dt = P + 2νW:ΔW` | `< 1e-12` relative |
| `D E_S/Dt = −2trS³ − P/2 − 2S:H + 2νS:ΔS` | `< 1e-12` relative |
| rotation sector pressure-blind (arbitrary `δH`) | max change `< 1e-10` |
| `H → H + cI` leaves `Dξ/Dt` fixed | max change `2.7e-14` |
| restricted-Euler reduction | exact |
| central finite differences | max abs error `2.5e-08` at `dt = 1e-5` (`O(dt²)` floor) |

**Navier–Stokes field** (`experiments/dns_checks/xi_dns_pressure.py`; unforced Taylor–Green, `96³`,
`Re = 400`, `t = 8.0`, 883 592 points; pressure solved spectrally):

| check | result |
|---|---|
| P1 `tr H = Δp = E_W − E_S` pointwise | max relative residual `8.3e-16` |
| P2 derived law vs. direct assembly, pointwise | max relative error `1.2e-11` |
| P3 `median |pressure|/|local|` | **1.26** (quartiles 0.45, 3.65) |
| P3 pressure term larger than local term | at **56.3 %** of points |
| P3 typical magnitudes (medians) | local `0.223`, pressure `0.282`, viscous `0.283` |
| P4 sign of `Dξ/Dt` changed by dropping `H_dev` | at **20.7 %** of points |
| P5 `R²` of the pressure term on `8⁴` bins in `(ζ, A, s, log Q)` | **0.368** (noise floor `0.011`) |
| P6 `⟨S:H_dev⟩` volume mean | `+4.1e-18` against `⟨\|S:H_dev\|⟩ = 1.50` (ratio `2.7e-18`) |
| P6a `⟨E_S⟩ = ⟨E_W⟩` (parent's homogeneity identity) | `2.137477` vs `2.137477` |
| P6b `⟨S:H_dev/E_S⟩` weighted mean | `−0.0156`, **not** zero |

Two further readings: `corr(local, pressure) = −0.33`, so the nonlocal term partially
*opposes* the local drive rather than reinforcing it; and the viscous group is heavily
tailed (median `0.283`, std `37`), the same `1/E_S` artifact the parent audit
documented — its typical size is comparable to the other two, its variance is not a
meaningful statistic.

## 16. Classification: **B**

> **B — ξ cleanly separates all local dynamics, while a single pressure-Hessian
> scalar remains as the nonlocal obstruction.**

With one amendment that must be stated: "all local dynamics" closes on
`(‖S‖, ξ, A, s)` — **four** local invariants, one more than the framework carried.
The separation itself is clean and complete:

* every term that is a function of `∇u` at the point is in the local group, in closed
  form;
* exactly **one** scalar is nonlocal, and it is the minimal possible such object (a
  single contraction, not a tensor);
* the locally-computable part of the pressure Hessian cancels **identically**.

**C is excluded** by the theorem of §13 — no cancellation exists, within the sector
family or its natural enlargement. **D is not claimed**: no new estimate, inequality,
monotone quantity, or closure theorem was obtained, and none is implied by anything
above.

## 17. What is genuinely new

Exact, and not present in the parent audit, the previous continuation, or either
source manuscript:

1. **The closed form of `Dξ/Dt`** under full incompressible Navier–Stokes, with every
   coefficient verified symbolically, and its `ζ`-form.
2. **The pressure reduction**: exactly one surviving scalar `S:H_dev/E_S`; the
   vorticity sector exactly pressure-blind; and — the sharp part — the cancelling
   piece is precisely the locally-known isotropic part `tr H = Q tanh ξ`.
3. **The cancellation impossibility theorem** (`∂F/∂E_S = 0 ⟺ F` is a function of the
   enstrophy alone), plus its extension to `F(E_S,E_W,P)` with the exceptional set
   identified as the Burgers-core configuration.
4. **`Δp = Q tanh ξ`**: the comparator, times the scale, *is* the pressure source.
5. **`∫S:H_dev = 0` versus `⟨S:H_dev/E_S⟩ ≠ 0`**: the obstruction is created by the
   pointwise weight, not by the pressure term's mean — and the unweighted global
   comparator is kinematically trivial. This locates precisely where a useful
   statement could still live.
6. **The retro-explanation of the parent's restricted-Euler residual** as the
   `s`-dependence, since RE closes exactly on `(‖S‖, ξ, A, s)`.
7. **`dξ/dζ = R = 1/(4h²)`** — the conditioning Jacobian is the framework's own
   response coordinate.

## 18. What is merely coordinate change

* **The two sector equations are textbook.** The enstrophy equation and the
  strain-norm equation (with `−2tr S³`, `−¼ω·Sω` in the half-norm convention, `−S:∇∇p`)
  are standard. `Dξ/Dt` is their difference of logarithmic rates: a *recombination*,
  not a new dynamical law. Nothing in §4 could contradict anything already known.
* **`A`, `s`, `P` are classical** (Betchov 1956; Ashurst et al. 1987; Lund–Rogers
  1994), as the parent audit already recorded.
* **`ξ` adds no information over `ζ`** — established in the previous continuation and
  unchanged here. The value of `ξ` in this audit is that it makes the *rate* structure
  additive (a difference of log-rates), which is what exposes the pressure term
  cleanly; that is a presentational gain, not an informational one.
* **The three-layer language is bookkeeping.** Calling `S:H_dev/E_S` a "pressure-
  environment coupling layer" would be relabelling: it is not a coordinate of the
  fluid element.

## 19. Consequences and the strongest justified next test

**Consequences for the framework.** The comparator is the natural coordinate in which
to *state* the strain–rotation balance dynamics, because it makes the pressure
coupling minimal and purely deviatoric. It is not a coordinate in which that balance
*closes*. Any diagnostic built on `ξ` alone will be wrong at the ~56 % of points where
the nonlocal term exceeds the local one, and will have the wrong sign at ~21 %.

**Explicitly, on Navier–Stokes.** No new estimate, inequality, monotonicity result,
coercive quantity, or closure theorem was obtained, and no regularity progress is
claimed. The one classical fact worth recording as context — *not* as a result of this
audit — is that `‖H‖_{L^q} ≤ C_q‖Δp‖_{L^q} = C_q‖Q tanh ξ‖_{L^q}` for `1 < q < ∞` by
Calderón–Zygmund: the obstruction is uncontrolled pointwise but its `L^q` norms are
controlled by the scale times the comparator. This is exactly why an integrated
statement is the only plausible route, and why it fails at `q = 1, ∞`.

**Recommended next test** (concrete, falsifiable, and with a definite negative outcome
available). §13 shows that the pressure term vanishes under the *unweighted* integral
and that the unweighted global comparator is trivial. The open interval is weighted
averages. Take a weight `w ≥ 0` and study

```
   d/dt ∫ w xi dx    or    Xi_w := (1/2) log( ∫ w E_W / ∫ w E_S )
```

and ask, for which `w` does the pressure contribution `∫ w S:H_dev/E_S` admit a bound
by local quantities? Two specific candidates, in order:

1. **`w = E_S`** (strain-energy weighting). Then the pressure contribution is
   `∫S:H_dev = 0` exactly — the weight cancels the `1/E_S`. So the `E_S`-weighted
   comparator has a **pressure-free** evolution. The question is whether
   `Ξ_{E_S} := ½log(∫E_S E_W/∫E_S²)` or the corresponding weighted average of `ξ`
   retains any diagnostic content, or whether the weighting has destroyed exactly what
   the comparator was for. This is a short calculation and should be done first: it is
   the only weight in sight for which the obstruction provably disappears.
2. If (1) is vacuous, test intermediate weights `w = E_S^α`, `0 < α < 1`, and use
   Calderón–Zygmund with Hölder to see whether the residual pressure term can be
   bounded by `‖Q tanh ξ‖_{L^q}` times a local factor, i.e. whether a *conditional*
   (not pointwise) comparator inequality exists. Expect this to fail for `α < 1`; a
   clean demonstration of the failure would close the line.

The prior from this audit is that (1) is vacuous and (2) fails — the comparator's
pressure coupling is structurally tied to the `1/E_S` weight that makes it a
comparator in the first place. Demonstrating that cleanly would be a legitimate
negative result and the right place to stop.

---

## Artifacts in this repository

| kind | path |
|---|---|
| exact derivation | `src/xi_exact.py` (49 symbolic checks) |
| numerical suite | `src/xi_numeric.py` (13 checks, counterexamples, special states) |
| Navier-Stokes check | `experiments/dns_checks/xi_dns_pressure.py` (8 checks, 96^3 field) |
| figure script | `src/figures_xi.py` |
| special-states table | `docs/special-states.md` |
| results | `results/xi_exact.json`, `results/xi_numeric.json`, `results/xi_dns_pressure.json` |
| figure | `figures/fig_xi_terms.png` |

The 96^3 sampled term fields (`xi_dns_terms.npz`, ~14 MB) are regenerated by the DNS
script rather than stored; see `../README.md`.
