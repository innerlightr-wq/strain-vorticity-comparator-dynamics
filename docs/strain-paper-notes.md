# Notes on the two audited manuscripts

Reading notes taken before any computation, kept here so the audit's starting
assumptions are inspectable. Both manuscripts are by the same author as this
repository; neither was modified by the audit. Page-level claims are paraphrased,
not quoted at length.

---

## Paper 1 — *Archimedean Compensation on the Thales Semicircle: Variance Redistribution, Harmonic Dynamics, and the Two Structural Parameters* (Revision 3, February 2026)

### What it sets up

A diagnostic vocabulary for constrained binary partitions. For a partition `(a,b)`
with `a + b = 1` on the unit-diameter Thales semicircle:

```
h = √(ab)   (altitude / geometric mean)      L = (b−a)/2   (lateral displacement)
D = R − h = ½ − h   (coherence deficit)      R = ½         (radius)
exchange identity:  L² + h² = ¼ ,  equivalently  L² = D(1−D)
```

Two parameters are claimed to organize the construction: the radius `R = 1/2`
(constraint surface, contraction rate `f′ = 1/2`, information unit `ln 2`,
exchange rate) and `χ = 1/e` (Shannon-optimal asymmetry under an exponential noise
model). Derived quantities include the random-partition baseline `⟨h⟩ = π/8`, the
"decision gap" `≈14 %`, and an angle hierarchy near `20–22.5°`.

### What the paper is careful about

The manuscript is explicit about its own limits, and the audit took those at face
value: a scope-and-non-claims section; a remark bounding the `1/e` derivation to
binary partitions with exponential noise and a Shannon criterion; the decision gap
reframed as a geometric difference plus a testable hypothesis; the previously
hypothesized cubic scaling `Δ ∝ h³` explicitly **retracted** on 17 tests; and an
epistemic table separating theorem from derived from hypothesis from empirical.

### What the audit used

Only the partition geometry — `a + b = 1`, `h`, `L`, `D`, `R`, `η = 2h` — applied to
the strain/rotation split. None of `1/e`, `π/8`, the angle hierarchy, the plasma
application, or the two-parameter claim was tested, and nothing here bears on them.

### What the audit found about it

* The Thales coordinates of the strain/rotation partition are **exact functions of
  `ζ`**, and the altitude family `{h, D, R, η}` is **even** in `ζ` — two-to-one,
  discarding exactly one bit per state (the sign). `L` is equivalent to `ζ`.
  So the framework's *state* coordinate is fine; the exposure is confined to
  diagnostics built from `η`, `D`, or `R` alone. (Audit 1, §1.1; Audit 2, §4–5.)
* Appendix B's **residual-altitude test** — fit `h` linearly on `a`, ask whether the
  residual predicts the observable — is a *nonlinearity detector*, not an
  independent-information test. Applied to a target constructed to have exactly zero
  independent `h` content it reports `ΔR² = 0.56`, `t = 1094`, `p < 10⁻³⁰⁰`; refining
  the control in `a` collapses it to `1.5e-5`. Neither published verdict is
  overturned (the shock case is a null; the two-slit case is true by construction
  because the observable *is* the geometric mean), but the test should not be used on
  a new system as it stands. The fix is one line: control for a flexible function of
  `a`, not a linear one. (Audit 1, F7.)
* Apex occupancy carries no evidential weight in this domain: for **any** homogeneous
  incompressible flow, `⟨E_S⟩ = ⟨E_W⟩` exactly, so the mean allocation sits at the
  apex by kinematics. In the manuscript's own vocabulary this is *occupancy without
  meaning*. (Audit 1, §1.9.)
* The framework's cubic landmark `b_* = 0.6823278…` is **not** the production-optimal
  allocation `b = 2/3` derived here; the cubic residual at `2/3` is exactly `−1/27`,
  a `2.30 %` gap. Reported as a non-identity precisely because the two are close
  enough to invite a retrospective fit. (Audit 1, F9.)

---

## Paper 2 — *Strain–Vorticity Interaction and Rotational Coherence: From Magnitude Comparison to Enstrophy Production* (interaction-response revision, September 2026)

### What it sets up

The exact kinematic decomposition `∇u = S + Ω`, `S:Ω = 0`,
`‖∇u‖²_F = ‖S‖²_F + ‖Ω‖²_F`, and the normalized imbalance

```
ζ = (‖Ω‖²_F − ‖S‖²_F)/(‖Ω‖²_F + ‖S‖²_F) = 2Q/‖∇u‖²_F
```

with `Q` the Hunt–Wray–Moin second invariant. The paper then **replaces** its own
earlier magnitude-comparison framing with the exact interaction term
`P = ω·Sω = |ω|²σ_eff`, re-derives the sharp bound
`|ω·Sω| ≤ √(2/3)‖S‖_F|ω|²` (a specialization of Wolkowicz–Styan), works five
canonical flows, states the exact production/dissipation balance, and argues that no
universal scalar threshold or corridor follows from it.

### What the paper is careful about

Unusually so, and this shaped the audit's standard. It states plainly that `ζ` is an
affine rescaling of `Q` (Hunt, Wray & Moin 1988) and of the unregularized Omega method
(Liu et al. 2016) — indeed `ζ = 2Ω − 1` exactly — and in 2D a reparametrization of
Okubo–Weiss (Okubo 1970; Weiss 1991); that the alignment bound is a kinematic identity,
not a threshold; that `⟨P⟩ > 0` is an empirical regularity of developed turbulence,
not a theorem; that a local `P/ε_ω` ratio has no invariant meaning without the
transport term; that the OpenAI forced-blowup construction is used as structural
context only, with the sign of `σ_eff` along its core explicitly left open; and, in
Remark 4.4, that collapsing three alignment cosines into the single number `A`
discards exactly the "which eigenvector" question that has been the interesting one
since Ashurst et al. (1987); see also the reviews of Meneveau (2011) and Johnson &
Wilczek (2024).

It also records a prior negative result: §7.2 argues that the author's binary-partition
("Two-Face Bridge") geometry is insufficient for the **production/dissipation**
partition `a = P₊/(P₊+D)`, `b = D/(P₊+D)`.

### What the audit added

The audit's target is a *different* partition from the one §7.2 rejects — the
strain/rotation energy split — so §7.2 does not already settle it. On that partition:

* The three-factor decomposition is exact:
  `P = ‖∇u‖³_F · g(ζ) · A` with `g(ζ) = 2b√a`, and `g` splits uniquely as
  `2h^{3/2}·e^{ξ/2}` (even × odd). The allocation enters as `a^{1/2}b`, which is
  neither the symmetric geometric mean nor linear in `a`. (Audit 1 §1.2–1.3;
  Audit 2 §9.)
* Optimizing the allocation gives the sharp pointwise bound
  `ω·Sω ≤ (4√2/9)‖∇u‖³_F`, attained iff `ζ = 1/3`, `S` axisymmetric, `ω` on its
  distinct eigenvector. Elementary, possibly known; the subsequent prior-art search did
  not locate it, and no priority is claimed
  ([`LITERATURE_REVIEW.md`](LITERATURE_REVIEW.md) §4).
* The **Burgers vortex core saturates the paper's own alignment bound exactly**:
  `A = √(2/3)` on the axis for every `a`, `Γ`, `ν`. And in self-similar coordinates
  `ζ` and `A` are exactly `a`-independent while `P ∝ a³` — which sharpens the paper's
  decisive counterexample by identifying its `ζ` trend as a core-radius effect rather
  than a dynamical one. (Audit 1, F5; `docs/canonical-flows.md`.)
* Remark 4.4's caveat is quantified: `(scale, ζ, A)` is exactly sufficient for
  instantaneous `P` but leaves `0.44 %` of restricted-Euler amplification variance
  unexplained — and Audit 3 identifies that residual exactly as the strain-state
  invariant `s`.
* The paper's "no corridor found" conclusion is reproduced and extended to the joint
  coordinate: the admissible set of `(ζ, A)` is a full rectangle, so `ζ` imposes no
  constraint on alignment whatsoever. The independence itself is classical — it is why
  the alignment literature exists — and what is added is the exact demonstration in
  these coordinates.

A claim-by-claim comparison of all of the above against the published literature,
including which items are classical and which are reparameterizations, is in
[`LITERATURE_REVIEW.md`](LITERATURE_REVIEW.md); the verified bibliography is
[`../references.bib`](../references.bib).

---

## What neither paper claims, and neither does this repository

No modification of established physics; no new vortex-identification method; no
universal threshold; and no progress on Navier–Stokes regularity. Where this audit
produces an exact statement, it is elementary algebra applied to standard equations,
and it is labelled as such.
