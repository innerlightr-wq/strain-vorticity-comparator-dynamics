# Manuscript plan

**Status: not written.** This directory is a placeholder holding the plan and the
scope decisions, so that the repository's claims and a future manuscript's claims
cannot drift apart. Nothing here is a draft.

## Is there a paper here?

Honest answer: at most a short note, and only if it is framed as a negative /
clarifying result. The audit produced no new estimate, inequality, monotonicity
result, coercive quantity, or closure theorem, and every exact statement in it is
elementary algebra applied to standard equations. What it does have is a coherent
chain of three sharp negative results and one clean structural identity, all
machine-checked — which is publishable as a note, and is not publishable as a
mechanism.

## The defensible claim set

In descending order of strength:

1. **The comparator's evolution equation and its pressure reduction.** `Dξ/Dt` in
   closed form; the vorticity sector exactly pressure-blind; exactly one surviving
   nonlocal scalar `S:H_dev/E_S`; and the point that the piece which cancels is
   precisely the locally computable piece `tr H = Δp = Q tanh ξ`.
2. **The impossibility theorem.** Every exactly pressure-free functional
   `∫Φ(E_S,E_W)dx` reduces to a functional of the enstrophy density alone, hence
   carries no comparator information. With the by-parts identity
   `∫w S:H = ∫p[S:∇∇w + ∇w·Δu]` as the mechanism, and the corollary that conditional
   and regional averages are obstructed too.
3. **The exact factorization** `P = ‖∇u‖³ · 2h^{3/2} · e^{ξ/2} · A` with its unique
   even × odd split, the sharp constant `4√2/9`, and the derived landmark `ζ = 1/3`
   as the balance of the two layers.
4. **The information-loss result** for symmetric partition coordinates: the altitude
   family is two-to-one in `ζ`, discarding exactly one bit, with the exact inverse
   `ζ = σ√(1−4h²)` and its conditioning `dζ/dh = −4h/ζ`.

Items 3 and 4 are elementary and may well be known; they should be presented as
observations with no priority claim. Item 1 is a recombination of two textbook sector
equations. Item 2 is the only statement that would carry a paper on its own, and it
is a negative one.

## What must not appear

* Any suggestion of progress on regularity. The Calderón–Zygmund remark is context
  and is classical; it is not a result of this work.
* Any framing of `ξ`, `h`, `L`, `D` or `R` as new physical quantities. They are exact
  reparameterizations of the second invariant `Q`; this must be stated early, not in
  a limitations section.
* Any reuse of the `1/e`, `π/8`, or angle-hierarchy material. None of it was tested
  and none of it bears on the fluid problem.
* The near-coincidence `2/3 ≈ 0.68233` as anything but a non-identity (exact cubic
  residual `−1/27`).

## Suggested shape, if written

A short note, ~8 pages, titled for the negative result rather than the coordinates —
something in the register of *"Exact pressure-freedom forces loss of strain–rotation
comparator information"*. Sections: conventions; the two sector equations; `Dξ/Dt`
and its three groups; the weight identity; the theorem and its corollaries; the
numerical confirmations; an explicit statement of what is classical. The
factorization and the sharp bound belong in a preliminaries section, not in the
abstract.

Before drafting, two things should be checked that this audit did not:

1. a literature search for the constant `4√2/9` and for the closed form of `Dξ/Dt` —
   both are elementary enough that prior appearances are likely;
2. whether the impossibility theorem, in the form "pressure-free ⟹ vorticity-only",
   is already standard folklore in the velocity-gradient-dynamics literature. It
   follows quickly from the pressure-free vorticity equation, so it may be.

If either check finds prior art, the note should cite it and shrink accordingly —
possibly to nothing, which would be an acceptable outcome.

## Source manuscripts

The two audited papers are listed in [`../CITATION.cff`](../CITATION.cff) under
`references`, with notes in
[`../docs/strain-paper-notes.md`](../docs/strain-paper-notes.md). Any manuscript
arising from this repository must cite both, and must state plainly which of their
claims it supports, which it sharpens, and which it corrects (the Appendix-B residual
test).
