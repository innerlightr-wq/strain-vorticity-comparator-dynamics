# The pressure Hessian: consolidated results

A single technical thread pulled out of audits 3 and 4
([`reports/xi-evolution-audit.md`](../reports/xi-evolution-audit.md),
[`reports/weighted-comparator-audit.md`](../reports/weighted-comparator-audit.md)),
where it is developed with full derivations and check logs. This file states the
results in one place, in the order in which they constrain each other.

Conventions as in [`../README.md`](../README.md). `H := Hess(p)`,
`H_dev := H − (tr H/3)I`.

---

## 1. Where the pressure can and cannot enter

From `DA/Dt = −A² − H + νΔA` with `A = ∇u`, split by symmetry:

```
DS/Dt = −(S² + W²) − H   + νΔS
DW/Dt = −(SW + WS)       + νΔW          ← no H
```

**The vorticity sector is exactly pressure-blind**, because `Hess(p)` is symmetric and
the antisymmetric sector sees only antisymmetric parts — equivalently `curl grad p = 0`.
Verified numerically by perturbing `H` arbitrarily: `D E_W/Dt` does not move
(`< 1e-10`).

Contracting gives the two sector balances, coefficient-for-coefficient against the
textbook forms:

```
D E_W/Dt = P + 2ν W:ΔW
D E_S/Dt = −2 tr(S³) − P/2 − 2 S:H_dev + 2ν S:ΔS
```

**EXACT.** So the pressure reaches any two-sector quantity through exactly one channel,
the strain sector, with coefficient `−2`.

## 2. Only the deviatoric part acts — and that is the nonlocal part

`S:(tr H/3)I = (tr H/3) tr S = 0` for incompressible flow, so `S:H = S:H_dev`
identically. What drops out is precisely the part that is *locally computable*:

```
tr H = Δp = E_W − E_S = Q ζ = Q tanh ξ
```

**The comparator, times the scale, is the pressure source.** The coordinate whose
evolution one wants to study is, up to the scale factor, exactly the source of the
field that obstructs it. And the piece of `H` that this source determines pointwise is
exactly the piece that cancels; what survives, `S:H_dev`, is the piece that requires
solving the nonlocal Poisson problem.

This is special to `ξ`. The same exercise on the production term is not clean:
`DP/Dt` contains `−ω·Hω = −(ωω)_dev:H_dev − (tr H/3)|ω|²`, and the second piece does
**not** vanish.

## 3. The comparator evolution

```
Dξ/Dt = ‖S‖_F[(1 + e^{2ξ}/2)A − s/√6]  +  (S:H_dev)/E_S  +  ν[W:ΔW/E_W − S:ΔS/E_S]
        └──────────── local ─────────┘     └─ nonlocal ─┘     └──── viscous ────┘
```

**EXACT**, verified symbolically (49/49 checks) and pointwise on a 96³ Navier–Stokes
field to a relative `1.2e-11` over 883 592 grid points.

Reading of the three groups:

* **Local.** Closes on four invariants `(‖S‖_F, ξ, A, s)` — one more than the
  framework carried, namely the Lund–Rogers strain state `s`. Vortex stretching
  appears twice, with combined weight `1 + e^{2ξ}/2`, because the same `P` feeds `+1`
  into the rotation sector and `−1/2` into the strain sector.
* **Nonlocal.** One scalar. Writing `S:H_dev = B‖S‖_F‖H_dev‖_F`, the term is
  `B‖H_dev‖_F/‖S‖_F`: structurally parallel to the alignment term `A`, but its
  magnitude factor is unbounded and nonlocal, whereas `|A| ≤ √(2/3)`.
* **Viscous.** Requires second derivatives, so it is not a function of the
  velocity-gradient state at all. It vanishes in the Euler limit; the pressure term
  does not.

**Restricted Euler** replaces `H` by its isotropic part, which contracts to zero
against `S`. So restricted Euler discards *exactly* the nonlocal term and nothing
else, and in that model `Dξ/Dt` closes exactly on `(‖S‖, ξ, A, s)` — which
retro-explains the `0.44 %` of restricted-Euler amplification variance that audit 1
found unexplained by `(ζ, A)`.

**Size in a real flow** (96³, `Re = 400`, `t = 8`): median `|pressure|/|local| = 1.26`;
the pressure term is the larger of the two at `56 %` of points; dropping it flips the
*sign* of `Dξ/Dt` at `21 %`; and a saturated `8⁴` binning of the local state
`(ζ, A, s, log Q)` explains only `37 %` of its variance (noise floor `1.1 %`). It is
also *anti*-correlated with the local term (`−0.33`), i.e. it partially brakes rather
than reinforces.

## 4. Can the nonlocal scalar be cancelled? — No

### 4.1 Pointwise: no

For any scalar `F(E_S,E_W)`, the pressure contribution to `DF/Dt` is
`(∂F/∂E_S)(−2 S:H_dev)`. It vanishes pointwise iff `∂F/∂E_S = 0`, i.e. iff `F`
ignores the strain sector — and then it is not a comparator. Concretely,
`D(log h)/Dt = −ζ Dξ/Dt`, so every combination `ξ + α log h` merely rescales the
pressure term by `(1 − αζ)`, which cannot vanish identically.

Enlarging the family to `F(E_S, E_W, P)` does not help: `DP/Dt` brings the *different*
functional `ω·H_dev ω`, and `S` and `(ωω)_dev` are generically linearly independent.
The exceptional set `S ∝ (ωω)_dev` is an axisymmetric `S` with `ω` on its distinct
axis — the Burgers-core state — a codimension-4 set, not a mechanism.

### 4.2 Integrated: the classical cancellation, and why it does not rescue anything

On a periodic (or decaying) domain, `∂_i∂_j S_ij = Δ(div u) = 0`, hence

```
∫ S:H_dev dx = ∫ p ∂_i∂_j S_ij dx = 0        (exactly; verified to 1e-17–1e-18)
```

This is the identity behind Betchov's relation and the pressure-free integrated strain
balance. It tempts the conclusion that weighting the `ξ` equation by `E_S` and
integrating removes the obstruction. **That conclusion is false**, and audit 4 retracts
it: the weight evolves, and `ξ` itself depends on `E_S`. Carrying the product rule,

```
Φ = E_S^α ξ   ⟹   Φ_{E_S} = E_S^{α−1}(α ξ − ½)
pressure part of dJ/dt = −2 ∫ Φ_{E_S} S:H_dev dx
α = 1  ⟹  −2∫ξ S:H_dev dx  ≠ 0
```

Measured relative miss `|∫X|/∫|X|`: the unweighted integral misses by `1e-17`; the
`ξ`-weighted one by `1.1–2.1 %` on random solenoidal fields, `7.1 %` on a structured
asymmetric field, and `26 %` on a Navier–Stokes snapshot. (For sufficiently symmetric
fields — ABC, Taylor–Green, plane shear — it can vanish accidentally; that is a
property of those fields, not of the functional.) No `α` works.

### 4.3 Which weights survive at all

Integrating by parts twice:

```
∫ w S:H dx = ∫ p [ S:∇∇w + ∇w·Δu ] dx
```

Every surviving term carries at least one derivative of `w`. **Cancellation holds for
spatially constant weights and nothing else.** Two corollaries:

* any state-dependent weight `w(E_S,E_W)` varies in space wherever the state does, so
  it reintroduces the pressure;
* an indicator is a non-constant weight, so **conditional and regional comparator
  statistics are obstructed too**, by a pure boundary term. Measured over a half-box,
  over `{ξ > 0}`, and over `{E_S > median}`: `1.4 %`, `1.3 %`, `0.6 %`.

## 5. The impossibility result (scoped)

> For `J = ∫Φ(E_S,E_W)dx` on the periodic incompressible domain, `dJ/dt` is exactly
> independent of `H_dev` **if** `Φ_{E_S}` is constant, i.e. `Φ = c E_S + ψ(E_W)`
> (exact); and **only if**, under the two-region nondegeneracy condition stated
> below (not proved unconditionally for all incompressible flows).
> Since `∫E_S = ∫E_W` exactly, every such `J` equals `∫ψ̃(E_W)dx` — a functional of
> the **enstrophy density alone**, and therefore blind to the comparator.

Sufficiency is exact. Necessity is established under an explicit two-region
nondegeneracy and verified for every weight and exponent tested. The theorem's content
is really a converse: `Dω/Dt = Sω` has no pressure term, so any functional of `ω`
alone is automatically pressure-free — the result is that within this class there is
**nothing else**.

Demonstrated blindness: rearranging the `E_S` field at fixed `E_W` field leaves a
pressure-free `J` unchanged to `2e-16` while the comparator functional moves by `97 %`;
and across five genuine incompressible fields with `∫E_S = ∫E_W = 1` and `std(ξ)`
differing by a factor `8.7`, the entire strain-sector content a pressure-free
functional may carry is identical to `1.1e-16`.

The surviving class is the standard one: `∫ψ̃(E_W)` are the enstrophy moments /
vorticity `L^p` norms, dissipative under viscosity iff `ψ̃` is nondecreasing and
convex, and in the Euler limit evolving by `∫ψ̃′(E_W)P` — not conserved, not monotone,
sign-indefinite.

## 6. Status

**Within the tested local two-sector class, and under the stated nondegeneracy
condition, exact pressure-freedom and comparator sensitivity are mutually
exclusive.** The reason is structural rather than technical: the pressure enters only
through the strain sector, so being blind to it means being blind to strain, and the
comparator is precisely the quantity that is not.

No estimate, inequality, coercive quantity, monotonicity result, or closure theorem is
claimed. The one classical fact worth keeping in view — recorded as context, not as a
result of this work — is that `‖H‖_{L^q} ≤ C_q‖Δp‖_{L^q} = C_q‖Q tanh ξ‖_{L^q}` for
`1 < q < ∞` by Calderón–Zygmund: the obstruction is uncontrolled pointwise but its
`L^q` norms are controlled by the scale times the comparator, and the bound fails
exactly at `q = 1, ∞`.
