# Manuscript plan

**Status: not written.** This directory is a placeholder holding the plan and the
scope decisions, so that the repository's claims and a future manuscript's claims
cannot drift apart. Nothing here is a draft.

## Is there a paper here?

**Assessment after the literature review: outcome C — repository documentation only.**
(A = standalone paper justified; B = short note; C = repository note/documentation only;
D = not determinable.) Recorded September 2026, after the prior-art audit in
[`../docs/LITERATURE_REVIEW.md`](../docs/LITERATURE_REVIEW.md); this supersedes the
earlier reading of "at most a short note", which was written before that audit.

The reasoning is the relationship to prior work, not a judgement about effort:

* Items 3 and 4 below are **classical or reparameterizations**. `ζ` is an affine
  rescaling of `Q` (Hunt, Wray & Moin 1988) and exactly `2Ω − 1` for the published Ω
  measure (Liu et al. 2016), so `ξ` is its logit; the sign-blindness of a symmetric
  coordinate is the reason every standard discriminant is signed; and the independence of
  alignment from the magnitudes is the founding observation of the alignment literature
  (Ashurst et al. 1987).
* Item 1 is an **immediate combination of two standard sector equations**
  (Nomura & Post 1998), and its most striking-sounding component — the reduction of the
  pressure coupling to one scalar — is coordinate compression. The genuine reduction of
  the anisotropic pressure Hessian is already published, and is stronger:
  Carbone, Iovieno & Bragg (2020).
* Item 2 is the only surviving candidate, and it is small: `∫S:H_dev = 0` is classical
  and now completely classified (Carbone & Wilczek 2022; Zhou & Yang 2023); the
  classification itself follows from standard integration by parts; necessity is
  conditional on a stated nondegeneracy rather than proved; and it sits beside a
  published completeness theorem that is more general in its own domain.

One conditional classification in a narrow functional class, provable by standard
machinery and adjacent to a stronger published theorem, is a documented repository
result — not a paper. It is fully stated in
[`../reports/weighted-comparator-audit.md`](../reports/weighted-comparator-audit.md) §7
and [`../docs/pressure-hessian-results.md`](../docs/pressure-hessian-results.md) §5,
which is the appropriate venue.

**What would move this to B.** Removing the nondegeneracy condition so that necessity
holds unconditionally for all incompressible flows; *or* establishing that the
classification does not follow from Carbone & Wilczek (2022) by a reading this review
missed, and positioning it explicitly against that theorem; *or* extending the class
beyond `Φ(E_S,E_W)` far enough that the obstruction becomes a statement about
functionals of the velocity gradient generally. Nothing weaker justifies a submission.

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

Items 3 and 4 are elementary and are classified as known / reparameterized or classical
by the literature review; they carry no priority claim. Item 1 is a recombination of two
textbook sector equations. Item 2 is the only statement for which no equivalent was
identified in the reviewed literature — and, per the assessment above, it is not enough
for a standalone paper.

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

## The two pre-drafting checks: run, and their outcome

Both were carried out in September 2026
([`../docs/LITERATURE_REVIEW.md`](../docs/LITERATURE_REVIEW.md)).

1. **`4√2/9` and the closed form of `Dξ/Dt`.** Neither was located in the reviewed
   sources. Both remain elementary — a two-variable constrained optimization over
   classical ingredients, and half the difference of two standard logarithmic rates — so
   no priority is claimed for either, and `4√2/9` is classified *uncertain — more search
   needed* rather than new.
2. **Whether "pressure-free ⟹ vorticity-only" is folklore.** No equivalent statement was
   identified. The closest prior art is a *stronger* theorem about a *different* object:
   Carbone & Wilczek (2022) prove completeness of the Betchov homogeneity constraints and
   remark that their method extends to relations involving the pressure Hessian, and
   Zhou & Yang (2023) derive mixed pressure-Hessian moment constraints. Those classify
   static moment identities; the statement here concerns a functional's time evolution.
   The questions are adjacent, and the review could not derive ours from theirs — which is
   *not* evidence of significance, only of non-identity.

The checks therefore found no prior art that contradicts the results, and also found
that the framing shrank: see the outcome-C assessment above. Required public wording
throughout: *"No equivalent classification result was identified in the literature
reviewed."* Not "the first proof", and not "new" unqualified.

## Source manuscripts

The two audited papers are listed in [`../CITATION.cff`](../CITATION.cff) under
`references`, with notes in
[`../docs/strain-paper-notes.md`](../docs/strain-paper-notes.md). Any manuscript
arising from this repository must cite both, and must state plainly which of their
claims it supports, which it sharpens, and which it corrects (the Appendix-B residual
test). It must also cite the prior art identified in
[`../docs/LITERATURE_REVIEW.md`](../docs/LITERATURE_REVIEW.md) — at minimum
Liu et al. (2016), Nomura & Post (1998), Buaria & Pumir (2023),
Carbone, Iovieno & Bragg (2020), Carbone & Wilczek (2022) and Zhou & Yang (2023) — and
use the wording recommended in its §10 novelty matrix. Records:
[`../references.bib`](../references.bib).
