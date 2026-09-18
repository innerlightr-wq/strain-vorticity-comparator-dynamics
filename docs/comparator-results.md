# The comparator layer: magnitude, orientation, alignment

> *Reproduced verbatim from the audit run that produced it; only file paths have
> been updated to this repository's layout. Statements of the form "nothing outside
> this directory was modified" refer to that original run.*

A continuation of the allocation audit (`reports/allocation-audit.md`) (classification **A**, 44/44 exact checks).
Nothing in the parent directory, in either source manuscript, or in any other
repository was modified. All new work lives in `continuation_comparator/`.

**Verification status.** `src/comparator_exact.py` — 56/56 symbolic checks;
`src/comparator_numeric.py` — 13/13 numerical checks over 1 048 576 random
incompressible velocity-gradient states plus the parent's four DNS snapshots.

---

## 1. Previous result being extended

The parent audit established, exactly:

* `∇u = S + Ω`, `S:Ω = 0`, so `‖∇u‖²_F = ‖S‖²_F + ‖Ω‖²_F` and the two sector shares
  `a = ‖S‖²/‖∇u‖²`, `b = ‖Ω‖²/‖∇u‖²` satisfy `a + b = 1`;
* `ζ = b − a`, and the Thales coordinates are exact functions of it:
  `L = ζ/2`, `h = ½√(1−ζ²)`, `D = ½ − h`, `R = 1/(1−ζ²)`, `η = 2h`;
* `h, D, R, η` are **even** in `ζ` — they identify `+ζ` with `−ζ`;
* `P = ω·Sω = ‖∇u‖³_F · g(ζ) · A` with `g(ζ) = (1+ζ)√((1−ζ)/2)` and
  `A = P/(‖S‖_F|ω|²)`, `|A| ≤ √(2/3)`;
* `ζ` and `A` are functionally independent (joint image = full rectangle), so `ζ`
  alone determines neither the sign nor the size of `P`.

## 2. New comparator hypothesis

A symmetric magnitude coordinate loses orientation under exchange of the two
competing components; an antisymmetric comparator restores it; magnitude +
comparator may still be insufficient when the dynamics depends on relative
alignment. Proposed hierarchy: **magnitude → comparator → alignment.**

## 3. Definitions

```
exchange involution   ι : (‖S‖²_F, ‖Ω‖²_F) ↦ (‖Ω‖²_F, ‖S‖²_F)      ⇒  ζ ↦ −ζ
magnitude layer       M ∈ {h, D, R, η, |ζ|}      (ι-invariant)
comparator layer      C ∈ {ζ, L, χ, ξ}, minimal σ = sign ζ   (ι-odd / ι-flipping)
rapidity              ξ = artanh ζ = ½ log(b/a),   ζ = tanh ξ,   h = 1/(2 cosh ξ)
alignment layer       A = ω·Sω / (‖S‖_F|ω|²) = Σᵢ (λᵢ/‖S‖_F) cos²θᵢ
scale                 Q = ‖∇u‖²_F
```

## 4. Exact symmetry / information-loss result — **EXACT**

**Theorem 4.1 (what is lost).** Every member of the altitude family satisfies
`G(ζ) = G(−ζ)`. Its level sets are exactly the orbits `{ζ, −ζ}` of `ι`; the quotient
map is `ζ ↦ |ζ|`, and `h = ½√(1−|ζ|²)` factors through it. For `h ∈ (0, ½)` every
fibre has **exactly two** points; the fibre over `h = ½` (`ζ = 0`) is a single point.

**Corollary 4.2 (size of the loss).** The information discarded by the altitude
family is **exactly one bit per state**, and no less: any `F` making `(h, F)`
injective must separate each two-point fibre, so `|range F| ≥ 2`, and `σ = sign ζ`
attains that bound.

**Measured loss (EMPIRICAL).** That bit is very nearly a *full* bit in practice —
`h` carries almost no statistical information about its own branch:

| ensemble | `P(σ=+1)` | `H(σ)` | `H(σ | h)` |
|---|---|---|---|
| Gaussian solenoidal field | 0.4642 | 0.9963 bits | **0.9845 bits** |
| DNS `t=4` | 0.5083 | 0.9998 | **0.9221** |
| DNS `t=6` | 0.4219 | 0.9824 | **0.9412** |
| DNS `t=8` | 0.3683 | 0.9494 | **0.8873** |
| DNS `t=10` | 0.3805 | 0.9584 | **0.8866** |

## 5. Comparator reconstruction result — **EXACT**

**Theorem 5.1 (exact inverse).** With `σ(0) := +1`,

```
ζ = σ · √(1 − 4h²) ,     and conversely   h = ½√(1−ζ²),  σ = sign ζ.
```

`(h, σ) ↔ ζ` is a bijection. Adjoining scale,

```
‖S‖²_F = Q(1−ζ)/2 ,    ‖Ω‖²_F = Q(1+ζ)/2 ,    |ω|² = Q(1+ζ),
```

so `(Q, h, σ)` reconstructs **both component magnitudes exactly** — verified to
`3.6e-15` on 10⁶ random states. It does **not** reconstruct `S` and `ω`: those carry
5 `SO(3)`-invariants and `(Q, ζ)` fixes 2.

**Counterexample (the comparator must be odd).** `C = |ζ|` or `C = ζ²` repairs
nothing: `(h, C)` is constant on each class.

**Theorem 5.2 (conditioning of the repair) — new.** `ζ² = 1 − 4h²` gives
`dζ/dh = −4h/ζ`. The inverse is exact but its condition number diverges at the fold
`ζ → 0`. Numerically: over 10⁶ states the reconstruction is exact to `2.4e-14` for
`|ζ| > 10⁻²`, while the worst residual `4.5e-11` occurs at `|ζ| = 4.2e-6`, matching
the predicted amplification `4h/|ζ| = 4.8e5` exactly. **Restoring the sign is free;
restoring the magnitude through the quotient is not numerically free near balance.**

## 6. Dynamic sufficiency test — **the answer is NO**

**Theorem 6.1 (insufficiency, by symmetry).** `Q`, `|ζ|`, `σ` — indeed `‖S‖_F`,
`‖Ω‖_F`, `|ω|` — are invariant under the **independent** action
`(S, ω) ↦ (R₁SR₁ᵀ, R₂ω)`, `R₁,R₂ ∈ SO(3)`. `P = ω·Sω` is invariant only under the
**diagonal** action `R₁ = R₂`. A function of independent-action invariants is constant
on independent-action orbits, and `P` is not: holding `S` and `|ω|` fixed and rotating
`ω` alone sweeps `A` over `[λ_min/‖S‖_F, λ_max/‖S‖_F]`. Hence `P` is **not** a
function of (scale, magnitude, comparator). ∎

This is stronger than a counterexample: it identifies *why* no such function can
exist, and it yields the general criterion of §H below.

**Exact witnesses.** `S = diag(1, 0, −1)`, `|ω|² = 2`, so `Q = 3`, `h = √2/3`,
`σ = −1` for all three:

| `ω` | `Q` | `h` | `σ` | `cos²(ω, e_asc)` | `A` | `P` |
|---|---|---|---|---|---|---|
| `(√2,0,0)` | 3 | 0.471405 | −1 | `(0,0,1)` | `+1/√2` | **+2** |
| `(0,√2,0)` | 3 | 0.471405 | −1 | `(0,1,0)` | `0` | **0** |
| `(0,0,√2)` | 3 | 0.471405 | −1 | `(1,0,0)` | `−1/√2` | **−2** |

Identical scale, magnitude and comparator; three different productions with all three
signs. **What differs is exactly the direction cosines of `ω` in the strain eigenframe**
— nothing else: the eigenvalues, `‖S‖_F`, `|ω|`, `ζ`, `h`, `D`, `R` are all identical.

**Statistical version (COMPUTATIONAL).** Binning 10⁶ random states into 200 joint
`(Q, h, σ)` cells: **99.97 %** of the variance of `A` survives inside the cells, and
**100 %** of cells with more than 50 members contain both `P > 0` and `P < 0`.

## 7. Alignment analysis — the minimal third coordinate

**Proposition 7.1.** Given `(Q, ζ)`, the map `A ↦ P = Q^{3/2} g(ζ) A` is affine with
nonzero slope, so a **single real scalar** is both necessary (P varies continuously
over an interval) and sufficient. `A` is canonical: if `P = Φ(Q, ζ, F)` for some `F`,
then at fixed `(Q, ζ)` the map `F ↦ A` is a bijection, i.e. `F` is a fibre-preserving
reparameterization of `A`. No tensor, no eigenvalue triple, and no set of three
direction cosines is needed for `P`.

**Caveat carried forward from the parent audit.** `A` is minimal *for instantaneous
`P`*, not for the dynamics: states identical in `(Q, ζ, A, P)` but with opposite
strain-state parameter `s` diverge under restricted-Euler evolution (0.0285 nats over
`t = 0.8`; 0.44 % of amplification variance unexplained by `(ζ, A)`).

## 8. Archimedean / Thales quotient interpretation — **EXACT**

**Theorem 8.1.** The Thales *point* `(L, h)` is **faithful**, not a quotient: `L = ζ/2`
recovers `ζ`. The quotient arises only when the state is read through the **altitude
family** `{h, D, R, η}`, which is precisely the algebra of functions invariant under
reflection of the semicircle about its vertical axis. So:

```
   allocation line  [−1,1] ∋ ζ         (faithful; = Thales point, = L, = χ, = ξ)
              ↓  quotient by  ι: ζ ↦ −ζ
   altitude line   [0, ½] ∋ h          (the diagnostic scalars η, D, R live here)
```

The correct state object is `(Thales point)` — equivalently `(h, σ)` — and **not** the
altitude alone. The framework's own `χ` and `ξ` already carry the branch; the exposure
is confined to diagnostics built from `η`, `D`, or `R` alone.

**Theorem 8.2 (the fold is tangential).** `dh/dζ = −ζ/(2√(1−ζ²)) → 0` as `ζ → 0`: the
two sheets meet tangentially at the apex. This is the geometric statement behind
Theorem 5.2 — the sign is lost smoothly, which is exactly why recovering the magnitude
from `h` is ill-conditioned there, and why a symmetric coordinate is *insensitive to
motion through balance*.

**Formulation.** *The symmetric Archimedean geometry is not wrong; it is the quotient
geometry obtained after forgetting orientation, and the comparator is the section that
selects the sheet.*

## 9. Exact identities (all symbolically verified)

**The canonical split of the allocation factor — the central new identity.**
`g > 0` on `(−1,1)`, so `log g` decomposes uniquely into even + odd parts;
exponentiating gives a **unique** factorization:

```
g(ζ) g(−ζ) = 4h³                         (symmetric invariant)
g(ζ) / g(−ζ) = e^{ξ}                     (pure comparator)
────────────────────────────────────────────────────────────
g(ζ) = 2 h^{3/2} · e^{ξ/2}               (unique even × odd split)
      = 2 (ab)^{3/4} (b/a)^{1/4}          (same thing as a monomial identity)
```

the last line being `a^p b^q = (ab)^{(p+q)/2}(b/a)^{(q−p)/2}` with `(p,q) = (½,1)`:
**symmetric exponent 3/4, antisymmetric exponent 1/4**. The nonvanishing of that `1/4`
is the exact reason no function of the altitude family can equal `g`.

**Production in Archimedean coordinates** (verified to `8.0e-13` relative on 10⁶
states):

```
P = ‖∇u‖³_F · 2h^{3/2} · e^{ξ/2} · A
     scale     magnitude   comparator  alignment
```

This is exactly the four-factor schematic of the hypothesis, and it is exact rather
than schematic.

**The landmark is a balance of the two layers.**

```
d log g_sym/dζ = −3ζ / [2(1−ζ²)]      (maximal at the apex ζ = 0)
d log g_cmp/dζ = +1  / [2(1−ζ²)]      (strictly increasing)
sum = (1−3ζ)/[2(1−ζ²)] = 0   ⟺   ζ = 1/3
```

The parent audit's landmark `ζ = 1/3` is exactly where the magnitude layer's pull
toward balance cancels the comparator layer's push toward rotation. See
`figures/fig_comparator_split.png`.

**Branch asymmetry of the sharp envelope.** `g(ζ)/g(−ζ) = e^{ξ} = √((1+ζ)/(1−ζ))`,
verified exactly at `|ζ| = 1/3, 0.6, 0.8` → ratios `1.414214, 2.000000, 3.000000`.
A bound written in symmetric coordinates alone must take the larger branch, and so
**overestimates the weaker branch by exactly `e^{|ξ|}`**, which diverges as `|ζ| → 1`.

**Why rapidity is the natural comparator for dynamics.**

```
Dξ/Dt = ½ [ D_t log‖Ω‖²_F − D_t log‖S‖²_F ]
```

— the comparator's material derivative is exactly half the *difference of the two
sectors' logarithmic growth rates*, additively separable in the two sectors. No even
coordinate has this property; indeed `dh/dζ → 0` at the fold.

## 10. Counterexamples

1. **Even "comparators" repair nothing**: `(h, |ζ|)` and `(h, ζ²)` are constant on each
   equivalence class.
2. **Magnitude + comparator ≠ dynamics**: the three exact states of §6 (same
   `Q, h, σ`; `P = +2, 0, −2`).
3. **Branch pairs**: `ζ = ±1/3` share `h = 0.471405` exactly but have production
   envelopes differing by a factor `√2`; `ζ = ±0.8` share `h = 0.3` and differ by a
   factor 3.
4. **Endpoint identification**: pure extensional strain (`ζ = −1`) and solid-body
   rotation (`ζ = +1`) share `h = 0`, `D = ½`, `R = ∞`; only `σ` separates them — and
   in both cases `A` is *undefined* (`0/0`), so the third layer degenerates exactly
   where the first two do.

## 11. Numerical checks

| test | result |
|---|---|
| N1 endpoints / symmetry pairs | exact values reproduced |
| N2 reconstruction `ζ = σ√(1−4h²)` | max error `2.4e-14` (`|ζ|>10⁻²`, n = 1 033 355); worst overall `4.5e-11` explained exactly by the fold conditioning |
| N2 magnitude recovery | max error `3.6e-15` |
| N3 `P = Q^{3/2}·2h^{3/2}·e^{ξ/2}·A` | max relative error `8.0e-13` over 1 048 576 states |
| N4 insufficiency | 99.97 % of `var(A)` survives inside `(Q,h,σ)` cells; 100 % of cells hold both signs of `P` |
| N5 exact witnesses | `P = +2, 0, −2` at identical `(Q,h,σ)` |
| N6 branch asymmetry | ratio `= e^{ξ}` to `<1e-12` |
| N7 discarded information | `H(σ|h) = 0.887–0.985` bits |

## 12. Classification: **B**

> **B — the comparator restores information lost by the symmetric geometry, but
> magnitude + comparator still does not determine the dynamics because alignment
> remains independent.**

Strongest support for each half:

* **Restoration** (Theorem 5.1): `ζ = σ√(1−4h²)`, an exact inverse; the loss is
  exactly one bit (Corollary 4.2) and `σ` attains the minimum.
* **Insufficiency** (Theorem 6.1): `Q, |ζ|, σ` are `SO(3)×SO(3)`-invariant while `P`
  is only diagonally invariant — plus the exact triple with `P = +2, 0, −2`.

Outcome **C** is excluded by that theorem. Outcome **A** is excluded *relative to the
altitude family* — the comparator is not a reparameterization of `h`, it is
strictly additional — but note the honest caveat in §14. Outcome **D** is not claimed.

## 13. What is genuinely new

Structural, exact, and (to our knowledge) not stated in either manuscript or the
parent audit:

1. **The unique even × odd factorization** `g(ζ) = 2h^{3/2} e^{ξ/2}`, hence
   `P = ‖∇u‖³ · 2h^{3/2} · e^{ξ/2} · A` — the hypothesis's four-layer schematic
   realized exactly, in the framework's *own* coordinates (`h` = altitude,
   `ξ` = rapidity).
2. **The landmark as a balance**: `ζ = 1/3` is exactly the zero of
   `d log g_sym/dζ + d log g_cmp/dζ`.
3. **The price of the quotient, quantified twice**: the envelope branch asymmetry is
   exactly `e^{|ξ|}` (an estimate written symmetrically is that far from sharp), and
   the magnitude reconstruction is ill-conditioned at the fold with amplification
   `4h/|ζ|`.
4. **Exactly one bit**, with the minimality proof and its measured value in real
   fields (`H(σ|h) ≈ 0.89–0.98` bits).
5. **A criterion for when the third layer is forced** (§H): the observable's
   invariance group. Diagonal-only ⇒ alignment layer required; full product group ⇒
   two layers suffice.

## 14. What is only reparameterization

Stated plainly, because it is the easiest thing to get wrong here:

* **The comparator is not a new coordinate for this system.** `C(a,b) = b − a` *is*
  `ζ`, the coordinate the parent audit already used, and `σ = sign ζ = sign Q_HWM` is
  the sign of the classical `Q`-criterion. Nothing in §5 adds information to `ζ`;
  it restores information that the *altitude family* discarded. Relative to `ζ`,
  the comparator layer is empty.
* **The factorization adds no information.** `(h, ξ) ↔ ζ` is a bijection, so
  `P = ‖∇u‖³ 2h^{3/2} e^{ξ/2} A` and `P = ‖∇u‖³ g(ζ) A` are the same identity in
  different letters. Its value is organizational: it shows *where* the Archimedean
  altitude sits (the even half, with exponent 3/2) and *exactly what it omits*.
* **`A` is classical.** It is the vortex-stretching alignment of Betchov (1956) and
  Ashurst et al. (1987), normalized; the parent audit already said so.
* **No new mechanism.** Nothing here is a dynamical mechanism; every statement is
  kinematic algebra plus measurements on existing fields.

## 15. Implications for Navier–Stokes work

**No new estimate is claimed, and none was obtained.** No inequality, monotonicity
result, coercive quantity, or a-priori bound beyond the parent audit's sharp pointwise
constant has been produced here.

The one analytically relevant observation is negative and quantitative: any a-priori
estimate of the vortex-stretching term written in *symmetric* allocation variables
(`h`, `D`, `R`, `η`, `|ζ|`, or equivalently `|Q_HWM|`) cannot be sharp on both
branches, and is off by exactly `e^{|ξ|} = √((1+|ζ|)/(1−|ζ|))` on the weaker one —
a factor that diverges at the endpoints, i.e. precisely in the strain-dominated and
rotation-dominated regimes where such estimates matter. That is a statement about the
cost of symmetric bookkeeping, not a bound.

## 16. Recommended next mathematical test

Anchored on the exact identity of §9:

```
Dξ/Dt = ½ [ D_t log‖Ω‖²_F − D_t log‖S‖²_F ]
```

Substitute the exact evolution equations for the two sectors — the enstrophy equation
(`D_t‖Ω‖² = P − νε_ω + νΔ(...)`) and the strain equation (which contains
`−tr(S³) − ¼ ω·Sω − S:∇∇p + viscous`) — and ask whether the *comparator* has a
cleaner balance than either sector alone: specifically, whether the nonlocal pressure
Hessian enters `Dξ/Dt` only through the strain sector, and whether `−tr(S³)` and
`ω·Sω` appear in the combination fixed by the Betchov relation
`4⟨tr S³⟩ = −3⟨ω·Sω⟩`. If the pressure term does *not* cancel, the comparator has no
better closure than its parts and the line should be closed. This is a concrete,
falsifiable calculation with a definite negative outcome available, and it is the only
route we see from the present kinematic algebra toward anything dynamical.

A cheaper preliminary: measure `D_t log‖Ω‖²` and `D_t log‖S‖²` separately in the
existing DNS fields and test whether `Dξ/Dt` is better predicted by `(ζ, A)` than
either sector's rate is.

---

## Artifacts in this repository

| kind | path |
|---|---|
| exact derivation | `src/comparator_exact.py` (56 symbolic checks) |
| numerical suite | `src/comparator_numeric.py` (13 checks) |
| figure script | `src/figures_comparator.py` |
| results | `results/comparator_exact.json`, `results/comparator_numeric.json` |
| figure | `figures/fig_comparator_split.png` |

Reuses `src/svcore.py` and `src/generators.py`; see `../README.md` for the full
repository layout.
