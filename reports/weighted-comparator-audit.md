# Weighted comparator functionals: does exact pressure-freedom survive the product rule?

> *Reproduced verbatim from the audit run that produced it; only file paths have
> been updated to this repository's layout. Statements of the form "nothing outside
> this directory was modified" refer to that original run.*

Continuation of `reports/xi-evolution-audit.md` (classification **B**), itself continuing
`docs/comparator-results.md` (**B**) and the parent audit `reports/allocation-audit.md` (**A**). Nothing outside
this directory was created or modified; the parent's `svcore.py` and
`dns_taylor_green.py` are imported read-only.

**Verification status.** `src/weighted_exact.py` — 30/30 symbolic checks;
`src/weighted_numeric.py` — 11/11 checks on genuine periodic incompressible fields
and a Navier–Stokes snapshot; `experiments/synthetic_fields/synthetic_redistribution.py` — 7/7 checks.

> **Headline, and a retraction.** The previous audit's §19 recommendation — that
> weighting by `E_S` removes the nonlocal pressure term after integration — is
> **wrong**, and this audit retracts it. It ignored the product-rule term
> `ξ D E_S/Dt`. The correct pressure contribution to `d/dt ∫E_S ξ dx` is
> `−2∫ξ S:H_dev dx`, which does not vanish. Worse, a general theorem shows the
> failure is structural: **within the tested local two-sector class, and under the
> stated nondegeneracy condition (§7 (v)), every exactly pressure-free functional
> reduces effectively to a vorticity/enstrophy-only class, and therefore carries no
> comparator information at all.** Classification: **E**.

---

## 1. Prior result being tested

From `docs/comparator-results.md`: the pointwise law

```
Dξ/Dt = ‖S‖_F[(1 + e^{2ξ}/2)A − s/√6] + (S:H_dev)/E_S + ν[W:ΔW/E_W − S:ΔS/E_S]
```

with `D E_W/Dt = P + 2νW:ΔW` pressure-blind, `D E_S/Dt = −2tr(S³) − P/2 − 2S:H_dev +
2νS:ΔS`, and the integral identity `∫S:H_dev dx = 0`. The proposal under test: since
the pointwise pressure term carries a `1/E_S`, multiplying by `E_S` and integrating
should cancel it.

## 2. Exact definitions

```
Xi_α(t) := ∫ E_S^α ξ dx ,        ⟨ξ⟩_w := (∫ w ξ dx)/(∫ w dx) ,
J_Φ(t)  := ∫ Φ(E_S,E_W) dx ,     ξ = ½ log(E_W/E_S) .
```

Periodic box, `div u = 0`, unit density. All integrals are over the box; `∫` and `⟨·⟩`
differ only by the (constant) volume, which never matters below.

## 3. Full product-rule derivation

First, advection integrates away: for `div u = 0`, `u·∇Φ = ∇·(uΦ)`, so

```
d/dt ∫Φ dx = ∫ ∂_t Φ dx = ∫ DΦ/Dt dx        (verified symbolically)
```

Then, by the chain rule on `Φ(E_S,E_W)`,

```
d/dt J_Φ = ∫ [ Φ_{E_S} · D E_S/Dt  +  Φ_{E_W} · D E_W/Dt ] dx
```

and since only `D E_S/Dt` carries the pressure, with coefficient `−2`:

```
   ┌─────────────────────────────────────────────────────────────┐
   │  pressure contribution to d J_Φ/dt  =  −2 ∫ Φ_{E_S} S:H_dev dx │
   └─────────────────────────────────────────────────────────────┘
```

**This is the whole audit in one line.** Everything below is the consequence.

For `Φ = E_S^α ξ` the product rule gives two terms — the weight differentiates, *and*
`ξ` itself depends on `E_S` through `∂ξ/∂E_S = −1/(2E_S)`:

```
Φ_{E_S} = α E_S^{α−1} ξ  +  E_S^α ·(−1/(2E_S))  =  E_S^{α−1} ( α ξ − ½ )
```

## 4. Pressure contribution for general α

```
   pressure part of dXi_α/dt  =  −2 ∫ E_S^{α−1} ( α ξ − ½ ) S:H_dev dx
```

## 5. Does α = 1 cancel? — **NO**

At `α = 1`, `Φ_{E_S} = ξ − ½`, so

```
   pressure part = −2∫(ξ − ½) S:H_dev = −2∫ξ S:H_dev + ∫S:H_dev = −2∫ξ S:H_dev dx
```

The constant half cancels (that was the heuristic); the `ξ`-weighted half does not.
This is outcome **(3)** of the brief: *a new weighted pressure term `ξ S:H_dev`
survives*.

**Numerically** (`weighted_numeric.py` N1, N2, M), with
`rel := |∫X| / ∫|X|` measuring how far an integral is from cancelling:

| field | `∫S:H_dev` (rel) | `∫ξ S:H_dev` (rel) |
|---|---|---|
| random solenoidal ×3 | `9.6e-18`, `…`, `…` | **`0.0155`, `0.0210`, `0.0107`** |
| structured asymmetric (ABC k=2 + 2×cellular k=1) | `0` | **`0.0713`** |
| **Navier–Stokes snapshot (64³, Re=200, t=4)** | `1.2e-17` | **`0.263`** |

The unweighted identity holds to round-off everywhere; the `ξ`-weighted one misses by
1–26 %. In the actual Navier–Stokes field the failure is decisive
(`∫ξ S:H_dev = −0.173`).

**Caveat, stated because it cuts against the headline.** For sufficiently *symmetric*
fields the `ξ`-weighted integral can vanish accidentally: pure ABC (`rel = 1.1e-5`),
pure Taylor–Green cellular (`8.3e-5`), plane shear (`ξ ≡ 0` identically). These are
properties of those fields, not of the functional — mixing two different wavenumbers
destroys the cancellation immediately (`rel = 0.071`).

**No α works.** `Φ_{E_S} = E_S^{α−1}(αξ − ½)` is constant in the independent variables
`(E_S, ξ)` only if `∂/∂ξ = α E_S^{α−1} = 0`, i.e. `α = 0`; and `α = 0` gives
`Φ_{E_S} = −1/(2E_S)`, still not constant. Verified numerically over
`α ∈ {−1, −½, 0, ½, 1, 3/2, 2}`: the pressure part is nonzero at every α
(relative sizes `0.031, 0.032, 0.012, 0.008, 0.010, 0.005, 0.003`).

## 6. The general weight condition

The question "which `w` preserve cancellation?" has an exact answer. Integrating by
parts twice on the periodic box and using `∂_i∂_j S_ij = Δ(div u) = 0` and
`∂_j S_ij = ½Δu_i`:

```
   ∫ w S:H dx  =  ∫ p [ S : ∇∇w  +  ∇w · Δu ] dx
```

(verified symbolically on a `u = curl Ψ` test field and numerically for
`w = 1, ξ, E_S`). **Every surviving term carries at least one derivative of `w`.**
Hence:

> **Proposition (I).** The pressure cancellation `∫w S:H_dev dx = 0` holds for every
> incompressible field precisely when `w` is spatially constant. Any state-dependent
> weight `w = W(E_S,E_W)` has `∇w = W_S∇E_S + W_W∇E_W ≠ 0` wherever the sector
> energies vary, and reintroduces the pressure.

Numerically (N6): constant weight `rel = 9.6e-18`; every non-constant weight tested is
nonzero — `ξ`: `1.6e-2`, `E_W`: `2.6e-2`, `sigmoid(E_S)`: `8.3e-3`, `log E_S`:
`8.4e-4`, `E_S`: `5.4e-4`.

**Corollary (closing a door left open in the previous audit).** An indicator function
is a non-constant weight, so *conditional and regional* comparator statistics are
obstructed too, with the obstruction a pure surface term on the region boundary.
Verified (N9): `∫_Ω S:H_dev` over a half-box, over `{ξ > 0}`, and over
`{E_S > median}` gives relative `1.4e-2`, `1.3e-2`, `5.8e-3`. The previous audit's
hope that "the useful objects lie between the pointwise and the unweighted-global" is
therefore also closed.

## 7. The general Φ theorem

> **Theorem.** Let `J_Φ = ∫Φ(E_S,E_W)dx` on the periodic incompressible domain.
>
> (i) *Sufficiency.* If `Φ_{E_S} ≡ c` (a constant) then `dJ_Φ/dt` is exactly
> independent of `H_dev`, because the pressure part is `−2c∫S:H_dev dx = 0`.
>
> (ii) `Φ_{E_S} ≡ c` ⟺ `Φ(E_S,E_W) = c E_S + ψ(E_W)` for some `ψ`.
>
> (iii) *Collapse.* Since `∫E_S dx = ∫E_W dx` exactly on this domain,
> `J_Φ = c∫E_W + ∫ψ(E_W) = ∫ψ̃(E_W) dx` — **a functional of the enstrophy density
> alone.**
>
> (iv) *Consequence.* No exactly pressure-free member of the class depends on the
> pointwise pairing of the two sectors. `Φ = cE_S + ψ(E_W)` is additively separated:
> it sees `∫E_S` and the distribution of `E_W`, never `ξ`.
>
> (v) *Necessity.* If `Φ_{E_S}` is non-constant, the pressure part
> `−2∫Φ_{E_S}S:H_dev dx` does not vanish identically. **Status: established under an
> explicit nondegeneracy condition, not proved in full generality.** The mechanism is
> Proposition (I) plus a two-region argument: choose a weight that is nearly constant
> on each of two sets on which `S:H_dev` integrates to equal and opposite nonzero
> amounts (guaranteed, since the total is zero); the weighted integral is then
> `(w_1 − w_2)·I_1 ≠ 0`. Demonstrated numerically with a sigmoid weight (N6) and for
> every weight and α tested. Individual symmetric fields can still give zero (§5).

**Why this is the right answer rather than an accident.** The vorticity equation
`Dω/Dt = Sω` has no pressure term, so *any* functional of `ω` alone is automatically
pressure-free. The theorem is the converse: within this class, there is nothing else.
Exact pressure-freedom is not a property one can engineer — it is the signature of
having quietly discarded the strain sector.

## 8. Comparator-information test

Exact blindness, demonstrated by rearranging the `E_S` field at fixed `E_W` field
(N8; a test of functional form, not a dynamical evolution):

| state field | pressure-free `J = ⟨E_S + E_W²⟩` | comparator `Ξ_1 = ⟨E_S ξ⟩` |
|---|---|---|
| original | `0.0009949782` | `−0.00026751` |
| random rearrangement of `E_S` | `0.0009949782` | `−0.00028042` |
| comonotone (extremal) rearrangement | `0.0009949782` | `−0.00000852` |

Relative change: `J` moves by `2.2e-16` (i.e. not at all); `Ξ_1` moves by **97 %**.

## 9. Collapse to known invariants

Yes, completely. `J = ∫ψ̃(E_W)dx` with `E_W = |ω|²/2` is exactly the family of
enstrophy moments / vorticity `L^p` norms (`ψ̃(E) = E^{p/2}` gives `‖ω‖_p^p` up to a
constant). These are standard objects; the audit produces no new one. Reporting
`∫E_S` as an additional ingredient would be a repackaging, since `∫E_S = ∫E_W`.

## 10. Normalized averages

Classified separately, as required:

* **Denominator** `∫E_S dx` **is** pressure-free (`Φ = E_S`, `Φ_{E_S} = 1`) — and it is
  exactly `∫E_W`, i.e. the enstrophy. Consistent with the theorem.
* **Numerator** `∫E_S ξ dx` is **not** (§5).
* Therefore `d/dt ⟨ξ⟩_{E_S}` has pressure part `(−2∫ξ S:H_dev dx)/∫E_S dx ≠ 0`.
  Normalization cannot repair a pressure-carrying numerator.

## 11. Viscous terms of the surviving class

For `J = ∫ψ̃(E_W)`, using `W:ΔW = ½ΔE_W − ‖∇W‖²` (verified as a differential identity)
and one integration by parts:

```
   viscous part  =  −ν ∫ ψ̃''(E_W) |∇E_W|² dx  −  2ν ∫ ψ̃'(E_W) ‖∇W‖² dx
```

**Sign-definite (dissipative) iff `ψ̃` is nondecreasing and convex.** For `ψ̃(E) = E`
this is exactly `−2ν∫‖∇W‖²`, the classical palinstrophy dissipation. Nothing new; the
sign structure is inherited from the standard enstrophy balance.

## 12. Euler limit

`ν = 0`: `d/dt ∫ψ̃(E_W)dx = ∫ψ̃'(E_W) P dx`. Not conserved, not monotone, and
sign-indefinite (`P` is not sign-definite; `⟨P⟩ > 0` is an empirical regularity of
developed turbulence, not a theorem — parent audit F12). It collapses to the standard
evolution of vorticity `L^p` norms, i.e. exactly the objects of Beale–Kato–Majda-type
analysis. No new Euler invariant appears.

## 13. Numerical diagnostics

On the Navier–Stokes snapshot (64³, `Re = 200`, `t = 4`, adequately resolved by the
parent audit's criterion):

| α | `Xi_α` | pressure part of `dXi_α/dt` |
|---|---|---|
| −1.0 | `+3.140` | `−6.941` |
| −0.5 | `+0.478` | `−0.712` |
| 0.0 | `−0.123` | `−0.112` |
| +0.5 | `−0.410` | `+0.077` |
| +1.0 | `−0.785` | `+0.346` |
| +1.5 | `−1.671` | `+1.379` |
| +2.0 | `−4.400` | `+6.701` |

The pressure part never vanishes and is comparable to or larger than the functional
itself. Diagnostically, `⟨ξ⟩_{E_S} = −0.849` versus plain `⟨ξ⟩ = −0.123`: the
strain-weighted comparator is far more strain-dominated than the unweighted one, which
is sensible (strain-heavy points are strain-dominated) and is a genuine diagnostic
contrast — but it is not pressure-free, so it is not a candidate here.

## 14. Synthetic redistribution test

Five genuine periodic incompressible fields, all rescaled to
`∫E_S = ∫E_W = 1` exactly:

| field | `std(ξ)` | `frac(ξ>0)` | `⟨ξ⟩_{E_S}` | `⟨E_W²⟩` | pressure part of `dXi_1/dt` |
|---|---|---|---|---|---|
| random broadband (k 1–8) | 0.597 | 0.465 | `−0.277` | 1.640 | `+1.2e-4` |
| random small-scale (k 6–12) | 0.601 | 0.464 | `−0.281` | 1.672 | `+1.5e-3` |
| ABC (Beltrami, symmetric) | 1.141 | 0.500 | `−0.490` | 1.333 | `−5.1e-6` |
| Taylor–Green cellular (symmetric) | 5.192 | 0.494 | `−1.804` | 1.875 | `−1.1e-4` |
| ABC(k=2)+2×cellular(k=1) (asymmetric) | 0.768 | 0.489 | `−0.398` | 1.462 | `−2.5e-2` |

The `ξ` distributions differ by a factor **8.7** in width and the comparator
functional by a factor **6.5**, while **the entire strain-sector content permitted to
a pressure-free functional (`c∫E_S`) is identical to `1.1e-16` across all five**. The
remaining freedom, `∫ψ̃(E_W)`, can only read the enstrophy distribution (e.g.
`⟨E_W²⟩` above), which is not comparator information.

Incidental observation: the two random fields have nearly identical `ξ` statistics
despite occupying disjoint wavenumber bands — for Gaussian solenoidal fields the
single-point `ξ` distribution appears essentially spectrum-independent. Noted, not
pursued.

## 15. Classification: **E**

> **E — a broader impossibility theorem shows that no local two-sector comparator
> functional can eliminate pressure nontrivially.**

Outcome **A** is also established, as the special case that motivated the audit: the
`E_S` weighting does not cancel the pressure once the evolving weight is treated
correctly. **E** is reported as the classification because it explains *why*: the
failure at `α = 1` is not a bad choice of exponent but an instance of a structural
obstruction that admits no cure inside the class. **B** is excluded (there is no
unique pressure-free weighting; there is an infinite family, and all of it is
comparator-blind). **C** and **D** are excluded a fortiori: no pressure-free
comparator-sensitive observable exists, so neither its evolution structure nor a
coercive property can be discussed.

## 16. Strongest exact result

```
   pressure contribution to d/dt ∫Φ(E_S,E_W) dx  =  −2 ∫ Φ_{E_S} S:H_dev dx ,

   and       ∫ w S:H dx  =  ∫ p [ S:∇∇w + ∇w·Δu ] dx ,

   so exact pressure-freedom ⟺ Φ_{E_S} = const ⟺ Φ = c E_S + ψ(E_W)
   ⟹ (using ∫E_S = ∫E_W)  J = ∫ ψ̃(E_W) dx :  a functional of the vorticity alone.
```

## 17. What is genuinely new

1. **The product-rule correction and the retraction it forces.** `α = 1` does not
   cancel the pressure; the surviving term is `−2∫ξ S:H_dev dx`. The previous audit's
   §19 recommendation is withdrawn.
2. **The general theorem** (§7): pressure-freedom ⟺ `Φ = cE_S + ψ(E_W)`, and the
   collapse to an enstrophy-only functional via `∫E_S = ∫E_W`.
3. **The weight identity** `∫w S:H = ∫p[S:∇∇w + ∇w·Δu]`, which localizes the
   obstruction in derivatives of the weight and yields Proposition (I) — cancellation
   holds exactly for spatially constant weights and nothing else.
4. **The corollary for conditional/regional statistics**: indicators are non-constant
   weights, so regional comparator averages are obstructed by a boundary term. This
   closes the specific direction the previous audit recommended.
5. **The blindness demonstrations**: exact invariance of pressure-free functionals
   under rearrangement (97 % change in the comparator, `2e-16` in `J`), and five
   genuine fields with matched sector totals and factor-8.7 different `ξ` widths.

## 18. What is merely a reformulation

* `∫S:H_dev dx = 0` is classical — it is why the integrated strain balance is
  pressure-free and why Betchov's relation exists. Not a discovery here.
* `∫ψ̃(E_W)dx` are the standard enstrophy moments / vorticity `L^p` norms, and their
  viscous and Euler evolutions above are textbook.
* The pressure-freedom of the vorticity equation is classical. The only non-classical
  content is the **converse** (§7): within this class nothing *else* is pressure-free.
* `ξ` remains informationally equivalent to `ζ` (established two audits ago); nothing
  here changes that.

**No Navier–Stokes progress.** No new estimate, inequality, coercive quantity,
monotonicity, or closure theorem was obtained, and none is implied.

## 19. Should this branch continue? — **No. Stop.**

The stopping criterion set out in the brief (§O) is met, in both of its clauses at
once: every exactly pressure-free functional in the local two-sector class (a) loses
all comparator dependence and (b) reduces to already-known global quantities
(vorticity norms). The three natural escape routes are closed by the same identity:
power weights (§5), general state-dependent weights (§6), and conditional/regional
restriction (§6, corollary).

**Answer to the audit's central question.** *Exact pressure-freedom forces the
comparator information to disappear.* The two properties are not merely hard to
combine — within this class they are mutually exclusive, and the reason is structural:
the pressure enters the dynamics only through the strain sector, so the only way to
be blind to it is to be blind to the strain sector, and the comparator is precisely
the quantity that is not.

## 20. If continued anyway — the single strongest next test

The theorem is a statement about the class `∫Φ(E_S,E_W)dx`. The only honest way out
is to leave that class, and exactly one exit is not already closed by the weight
identity: functionals that are **not integrals of a pointwise function of the state**
— specifically, ones built from the pressure field itself, so that the obstruction is
absorbed rather than avoided. The concrete candidate is

```
   K(t) := ∫ [ Φ(E_S,E_W) + λ p Ψ(E_S,E_W) ] dx
```

for which `S:H_dev` can in principle be cancelled against the `Dp/Dt` terms rather
than annihilated. The test is short and has a definite negative outcome available:
`Dp/Dt` is itself nonlocal (it requires `Δ^{-1}` of the field), so the candidate
exchanges one nonlocal object for another, and the calculation should be done only to
confirm that the exchange is not favourable. We expect it is not, and we do not
recommend it; it is listed because the brief asks for the strongest remaining test,
not because we believe it will work.

---

## Artifacts in this repository

| kind | path |
|---|---|
| exact derivation | `src/weighted_exact.py` (30 symbolic checks) |
| numerical suite | `src/weighted_numeric.py` (11 checks, periodic fields + NS snapshot) |
| redistribution test | `experiments/synthetic_fields/synthetic_redistribution.py` (7 checks) |
| results | `results/weighted_exact.json`, `results/weighted_numeric.json`, `results/synthetic_redistribution.json`, `results/ns_snapshot_fields.npz` |

No figures: every result here is a table or an identity, and no chart would add to
them. See `../README.md` for the full repository layout.
