# CANONICAL FLOWS (Phase 4)

> *Reproduced verbatim from the audit run that produced it; only file paths have
> been updated to this repository's layout. Statements of the form "nothing outside
> this directory was modified" refer to that original run.*

All values below are **EXACT**: each velocity field is written in Cartesian
components and differentiated symbolically in `src/canonical_flows.py`
(sympy), so no cylindrical-frame convention can leak in. `div u = 0` is verified
for every flow. Undefined quantities are reported as undefined, never
regularized. Raw output: `results/canonical_flows.json`.

## Summary table

| flow | `zeta` | `a` | `b` | `h` | `L` | `D` | `A` | `P` |
|---|---|---|---|---|---|---|---|---|
| solid-body rotation `u=(-Omega y, Omega x, 0)` | `+1` | `0` | `1` | `0` | `+1/2` | `1/2` | **undefined** (`||S||=0`) | `0` |
| pure extensional strain `u=(a x, -a y, 0)` | `-1` | `1` | `0` | `0` | `-1/2` | `1/2` | **undefined** (`|omega|=0`) | `0` |
| simple shear `u=(gamma y,0,0)` | `0` | `1/2` | `1/2` | `1/2` | `0` | `0` | `0` | `0` |
| Lamb–Oseen (unstretched) | `(r,t)`-dependent, sweeps `(-1,+1)` | — | — | — | — | — | `0` exactly | `0` exactly |
| Burgers, on the axis | `(c^2-3)/(c^2+3)`, `c = Gamma/(4 pi nu)` | | | | | | `+sqrt(2/3)` **exactly** | `a |omega|^2 > 0` |
| Burgers, at `eta = 1` (`Re_Gamma=100`) | `+0.0718` | | | | | | `+0.5191` | `> 0` |
| Burgers, at `eta = 2` (`Re_Gamma=100`) | `-0.9932` | | | | | | `+0.5650` | `> 0` |

Three separate failures of magnitude-only reasoning are visible in this table
before any dynamics is invoked.

## 1. Solid-body rotation vs 2. pure extensional strain — the sign-loss failure

Solid rotation has `S = 0` identically: `zeta = +1`, `P = 0`, and `A` is **undefined**
(`0/0`), not zero — there is no strain direction for the vorticity to align with.
Pure extensional strain has `omega = 0`: `zeta = -1`, `P = 0` trivially, and `A` is again
undefined, for the opposite reason.

These two flows are physically as different as two incompressible flows can be.
In the Thales coordinates they are **the same point**:

```
h = 0 ,   D = 1/2 ,   R = infinity ,   eta_T = 0     for both
```

Only the signed coordinates `zeta` (or equivalently `L = zeta/2`) separate them. This is
Corollary 2.2 of `DERIVATIONS.md` realized in the two most elementary flows in fluid
mechanics. **FAILED** for the altitude/deficit family as a state coordinate here.

## 3. Simple shear — the apex is not a production landmark

`||S||_F^2 = ||W||_F^2 = gamma^2/2`, so simple shear sits **exactly at the Thales apex**:
`zeta = 0`, `h = 1/2` (maximal "coherence"), `D = 0` (zero "deficit"), `R = 1`. And yet
`A = 0` and `P = 0` exactly: the vorticity `(0,0,-gamma)` is orthogonal, in the
quadratic-form sense, to the only nonzero entries of `S`.

Combined with the exact witnesses of `DERIVATIONS.md` §4 — three states at `zeta = 0`
with `P = +4, 0, -4` — this settles the question the brief asks: the apex carries no
production meaning. The maximum of `h` is a critical point of `h`, nothing more
(`dh/dzeta = 0` at `zeta = 0` is an algebraic fact about the semicircle).

## 4. Lamb–Oseen vortex — `zeta` sweeps its entire range at `P = 0`

For the unstretched Lamb–Oseen swirl `u_theta = (Gamma/2 pi r)(1 - e^{-r^2/4 nu t})`,
vorticity is purely axial and `S` acts only in the `r`–`theta` plane, so

```
P = 0   and   A = 0   exactly, for every r and every t.
```

Meanwhile, in the similarity variable `eta = r/sqrt(4 nu t)`,

```
zeta(eta) = (eta^4 - (eta^2 - e^{eta^2} + 1)^2) / (eta^4 + (eta^2 - e^{eta^2} + 1)^2)
```

which runs from `+1` (core) through `0` at `eta = 1.1141` to `-1` (outer region):
**`zeta` takes essentially every value in its admissible range inside a single flow
whose enstrophy production is identically zero.** Any diagnostic that reads
`zeta > 0` as "rotationally favourable" will label the core of a purely decaying
vortex as favourable, and the same flow's outer region as strain-dominated, while
the exact production term is zero everywhere. Correspondingly the Thales
coordinates sweep from `h = 0` to `h = 1/2` and back to `h = 0` across the profile.

The dynamics is `D_t omega = nu Laplacian(omega)`: pure diffusion, no production
term at all. This is also the cleanest possible illustration of the brief's
warning: `P = 0` here coexists with *decaying* enstrophy — production is not the
enstrophy budget.

## 5. Burgers vortex — the principal counterexample, sharpened

With `u_r = -a r/2`, `u_theta = (Gamma/2 pi r)(1 - e^{-a r^2/4 nu})`, `u_z = a z`:

```
P = (Gamma^2 a^3 / 16 pi^2 nu^2) e^{-a r^2 / 2 nu}  >  0        (exact, all r)
```

### 5a. On the axis, the sharp alignment bound is saturated — **EXACT**

At `r = 0` the swirl shear vanishes and `S` is the imposed background strain
`diag(-a/2, -a/2, a)`, so `||S||_F^2 = (3/2) a^2` and, since `omega` is axial,

```
A(axis) = a / (a sqrt(3/2)) = sqrt(2/3) = 0.816496581...      for every a > 0, Gamma, nu.
```

**The Burgers vortex core is exactly the equality case of Theorem 3.1**: an
axisymmetric strain state with `omega` along its distinct eigenvector. Off axis,
the swirl contributes an `r`–`theta` shear entry that enters `||S||_F` but not
`omega . S omega`, so `A(r) = a/||S||_F(r) < sqrt(2/3)` strictly, and the deficit from
the bound measures exactly the swirl shear. We have not seen this stated in either
source paper; it is elementary once `A` is defined, and it explains *why* the
Burgers vortex is the extremal sustained structure in this coordinate system.

At `Re_Gamma = Gamma/nu = 100`: `A(eta)` falls from `0.8165` on the axis to `0.519` at
`eta = 1`, with a shallow minimum outside the core before rising again toward
`sqrt(2/3)` in the irrotational far field (where `||S||_F` is again dominated by the
background strain).

### 5b. The axis value of `zeta` does not depend on the axial strain — **EXACT**

```
zeta(axis) = (c^2 - 3)/(c^2 + 3) ,      c = Gamma/(4 pi nu)
```

— a function of the circulation Reynolds number alone. `d zeta/da = 0`. More
generally, at any fixed **similarity** radius `eta = r/sqrt(4 nu/a)`, every component
of `grad u` scales linearly in `a`, so

```
zeta(eta)  and  A(eta)  are exactly independent of a ,     while   P ∝ a^3 .
```

(Verified symbolically: `d zeta/da = d A/da = 0` at fixed `eta`; `a d(ln P)/da = 3`.)

**This refines the counterexample as stated in the RVP note.** That note reports
that increasing `a` drives `zeta` toward its strain-dominated extreme while
production rises — true at fixed *physical* radius, and we reproduce it
(`r = 1`, `Re_Gamma = 100`: `zeta` falls monotonically from `+0.909` at `a = 0.01` to
`-1.000` at `a = 100`). The refinement is that this trend is **a core-radius effect,
not a dynamical trend**: the core radius `sqrt(4 nu/a)` shrinks as `a` grows, so a
fixed `r` slides outward through the profile. In the self-similar frame the entire
`(zeta, A)` portrait of the Burgers family is invariant under changing `a`, and the
whole `a`-dependence of production sits in the **scale** coordinate
`||grad u||_F ∝ a`. So the correct statement is not "`zeta` moves the wrong way" but
"**`zeta` and `A` do not move at all, and `P` varies by six orders of magnitude**"
(`a = 0.01 -> 100` at `eta = 1`: `P = 8.6e-6 -> 8.6e+6`, `zeta ≡ +0.0718`,
`A ≡ +0.5191`). This is the sharpest available demonstration that the scale
coordinate is not optional.

### 5c. Does movement in `zeta` predict production? — **FAILED**

Directly answering the Phase 4 question, for the Burgers family:

* at fixed `eta`: `zeta` is constant while `P` varies over `10^12`. Movement in `zeta`
  predicts nothing because there is none.
* at fixed `r`: `zeta` decreases monotonically in `a` while `P` first rises and then
  collapses (the core withdraws inside `r`). `P` peaks at `a = 6.31` for `r = 1`,
  where `zeta = -0.521` — so along this path the extremum of production occurs at
  neither `zeta = 0` nor `zeta = 1/3`. *(This is not a violation of Theorem 3.4: that
  theorem optimizes the allocation at fixed `||grad u||_F`, whereas this path varies
  scale, allocation and alignment simultaneously. A landmark derived under one
  constraint does not transfer to a different family — a point worth keeping in view
  whenever a partition landmark is claimed to be physically located.)*
* `A > 0` throughout, and `sign(P) = sign(A)` identically. Alignment resolves the
  ambiguity that magnitude leaves open — as it must, by Theorem 3.2, and this is
  algebra rather than evidence.

## 6. What the canonical flows establish

| claim | status |
|---|---|
| `h`, `D`, `R` cannot distinguish pure rotation from pure strain | **EXACT**, realized in flows 1–2 |
| the Thales apex carries no production information | **EXACT**, realized in flow 3 (`h = 1/2`, `P = 0`) |
| `zeta` can sweep its whole range at `P ≡ 0` | **EXACT**, flow 4 |
| the Burgers core saturates `|A| <= sqrt(2/3)` | **EXACT**, flow 5 |
| in self-similar coordinates `zeta`, `A` are `a`-independent while `P ∝ a^3` | **EXACT**, flow 5 |
| movement in `zeta` predicts neither direction nor strength of production | **FAILED** for `zeta` alone |
| `A` fixes the sign and, with scale and allocation, the value of `P` | **EXACT** (identity, not evidence) |
