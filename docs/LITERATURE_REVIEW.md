# Literature review and prior-art audit

*Conducted September 2026, after the computational work was complete and published.
Purpose: to establish which results are classical, which are reparameterizations of
known quantities, which are new derivations of known ingredients, and which remain
apparently distinct after comparison with the literature. The review was carried out
to **reduce** novelty claims wherever prior work warrants it.*

**Method and its limits.** Bibliographic records were verified individually by DOI
content negotiation against publisher metadata (Crossref / DOI Foundation), and
claim-level comparison was done from abstracts and, where needed, full-text
statements. This is a targeted prior-art search in the velocity-gradient-dynamics
literature, not a systematic review: absence of a prior result here means *not found
in the sources reviewed*, never *does not exist*.

**Zotero.** The review itself was carried out without Zotero, which was not installed at
the time; the bibliography was verified instead against publisher metadata by DOI content
negotiation. Zotero 10.0.3 was installed afterwards, on 18 September 2026, and the library
was then built from this same verified bibliography: 25 items in a
`Strain–Vorticity Comparator Dynamics` collection with ten subcollections. Zotero's stored
metadata agrees with `../references.bib` on all 25 records — DOI, year, volume, pages and
journal, 0/25 discrepancies — so nothing in this review was revised as a result. Details,
including the local-API write flow, are in [`ZOTERO_SETUP.md`](ZOTERO_SETUP.md).

---

## 1. Verified bibliography

25 entries in [`../references.bib`](../references.bib), 23 with publisher DOIs. Three
records I initially recalled from memory had **wrong DOIs** that resolved to
unrelated papers (Vlaykov & Wilczek, Lund & Rogers, Meneveau 2011); all were corrected
against Crossref. This is recorded because it is the reason every entry was verified
rather than typed.

| key | work | role in this repository |
|---|---|---|
| `Betchov1956` | Betchov, *JFM* **1**, 497–504 | the homogeneity relation `4⟨tr S³⟩ = −3⟨ω·Sω⟩`; the reason integrated strain balances are pressure-free |
| `AshurstEtAl1987` | Ashurst, Kerstein, Kerr & Gibson, *Phys. Fluids* **30**, 2343 | vorticity–strain eigenframe alignment; the classical content of our `A` |
| `Vieillefosse1982`, `Vieillefosse1984` | *J. Physique* **43**, 837; *Physica A* **125**, 150 | origin of the restricted-Euler / local velocity-gradient system |
| `Cantwell1992` | *Phys. Fluids A* **4**, 782 | exact solution of the restricted Euler equation |
| `ChongPerryCantwell1990` | *Phys. Fluids A* **2**, 765 | `Q–R` topology classification of velocity-gradient states |
| `HuntWrayMoin1988` | CTR-S88, 193–208 | the `Q` criterion; `ζ = 2Q/‖∇u‖²_F` |
| `LiuEtAl2016` | *Sci. China Phys. Mech. Astron.* **59**, 684711 | the Ω vortex-identification method; `ζ = 2Ω_Liu − 1` exactly |
| `Okubo1970`, `Weiss1991` | *Deep-Sea Res.* **17**, 445; *Physica D* **48**, 273 | the two-dimensional strain/rotation discriminant |
| `JeongHussain1995` | *JFM* **285**, 69 | `λ₂`; the other standard strain/rotation-based criterion |
| `LundRogers1994` | *Phys. Fluids* **6**, 1838 | the strain-state parameter `s` used throughout |
| `WolkowiczStyan1980` | *Lin. Alg. Appl.* **29**, 471 | eigenvalue–trace bound specialized to `|A| ≤ √(2/3)` |
| `NomuraPost1998` | *JFM* **377**, 65 | structure and dynamics of vorticity **and** rate of strain; sector equations |
| `Meneveau2011` | *Annu. Rev. Fluid Mech.* **43**, 219 | review: Lagrangian velocity-gradient dynamics and models |
| `WilczekMeneveau2014` | *JFM* **756**, 191 | pressure-Hessian and viscous contributions to gradient statistics |
| `VlaykovWilczek2019` | *JFM* **861**, 422 | small-scale structure and its impact on the pressure field |
| `CarboneIovienoBragg2020` | *JFM* **900**, A38 | gauge symmetry and **rank reduction** of the anisotropic pressure Hessian |
| `CarboneWilczek2022` | *JFM* **948**, R2 | **completeness**: only two Betchov homogeneity constraints exist |
| `ZhouYang2023` | *Phys. Rev. Fluids* **8**, 024601 | homogeneity constraints on **mixed moments of velocity gradient and pressure Hessian** |
| `BuariaPumir2023` | *JFM* **973**, A23 | decomposition of `H` into **local isotropic** and **nonlocal deviatoric** parts |
| `YangEtAl2024` | *JFM* **983**, A19 | structure of the pressure Hessian in strong-vorticity regions |
| `JohnsonWilczek2024` | *Annu. Rev. Fluid Mech.* **56**, 463 | current review of multiscale velocity gradients |
| `BealeKatoMajda1984`, `Burgers1948` | | continuation criterion; the Burgers vortex |

## 2. Prior-art audit: the comparator coordinate `ζ`, `ξ`

The decisive question is **information content**, not notation. A monotone bijection of
a known coordinate is the same coordinate.

| prior coordinate | definition | relation to `ζ` | relation to `ξ` | equivalent information? | citation |
|---|---|---|---|---|---|
| `Q` (second invariant) | `½(‖Ω‖² − ‖S‖²)` | `ζ = 2Q/‖∇u‖²_F` | `ξ = artanh(2Q/‖∇u‖²)` | **yes**, once normalized by `‖∇u‖²` | `HuntWrayMoin1988` |
| `Ω` (Omega method) | `‖Ω‖²/(‖S‖²+‖Ω‖²)` | `ζ = 2Ω − 1` | `ξ = ½log(Ω/(1−Ω))` — the **logit** of `Ω` | **yes**, exactly | `LiuEtAl2016` |
| Okubo–Weiss `W_OW` | `s_n²+s_s²−ω²` (2-D) | `W_OW = −2Q_{2D}` | monotone in `ξ` at fixed scale | **yes** in 2-D | `Okubo1970`, `Weiss1991` |
| `Q–R` plane | `(Q,R)` invariants | `ζ` is `Q` normalized | — | `ζ` is the `Q` axis only | `ChongPerryCantwell1990` |
| `λ₂` | 2nd eigenvalue of `S²+Ω²` | not a function of `ζ` | — | no (different object) | `JeongHussain1995` |
| ratio `r = E_W/E_S` | sector energy ratio | `ζ = (r−1)/(r+1)` | `ξ = ½log r` | **yes** — `ξ` is literally `½log r` | used implicitly throughout |

**Finding.** `ζ` is an affine rescaling of the normalized second invariant and is
*exactly* the affine rescaling of the unregularized Ω method — a point the source
paper already states. `ξ = artanh ζ = ½log(E_W/E_S)` is a **monotone bijective
reparameterization** of that same quantity: it is the logit of `Ω_Liu`, or equivalently
half the log of the sector-energy ratio. It therefore carries **identical information**.

**Classification: KNOWN / REPARAMETERIZED.** The `ξ` coordinate must not be presented
as a new state variable. Its only defensible advantages are presentational: it is
additive under sector-ratio multiplication, its Jacobian `dξ/dζ = 1/(1−ζ²)` is regular
at the balance point where the symmetric altitude coordinate degenerates, and it makes
the rate equation a difference of logarithmic growth rates. No prior use of the
*logarithmic* form as a named coordinate was found in the reviewed sources, but that is
a notational observation, not an information-theoretic one.

## 3. Prior-art audit: the symmetric-quotient (parity) observation

The claim audited: `h = √(ab)`, `D`, `R`, `η` are even in `ζ`, so they identify the
rotation-dominated and strain-dominated branches; the sign restores them.

**Finding.** No turbulence coordinate system in the reviewed literature uses `√(ab)`;
the standard strain/rotation measures (`Q`, `Ω_Liu`, `Δ`, Okubo–Weiss) are all *signed*
by construction, precisely because the sign is the physically decisive part. The parity
property is therefore a statement about the **Archimedean/Thales coordinates imported
from the companion manuscript**, not a discovery about the turbulence literature, and
the "restoration" is the observation that one should use the signed invariant — which
is what the field already does.

**Classification: KNOWN / REPARAMETERIZED** (trivial parity property of a
non-standard coordinate). It should be presented as an internal consistency check on
the imported geometry, not as a result about turbulence.

## 4. Prior-art audit: the production factorization

`P = ‖∇u‖³_F · g(ζ) · A = ‖∇u‖³_F · 2h^{3/2} · e^{ξ/2} · A`.

1. **Has this exact algebraic factorization appeared?** Not found in the reviewed
   literature.
2. **Has an equivalent factorization appeared in other variables?** Yes in substance.
   `P = |ω|²‖S‖_F A` is definitional once `A` is defined as the normalized
   alignment; writing `|ω|² = 2bQ` and `‖S‖_F = √(aQ)` is a one-line substitution. The
   decomposition of production into a magnitude factor and an alignment factor, with
   the alignment expanded in strain-eigenframe direction cosines, is the standard
   framing of `AshurstEtAl1987` and of the reviews `Meneveau2011`,
   `JohnsonWilczek2024`.
3. **Is "scale × allocation × comparator × alignment" interpretive repackaging?**
   **Largely yes.** The four-way reading is a labelling of an algebraic identity that
   follows from the definitions in one substitution. It is useful as organisation; it
   is not a decomposition theorem.
4. **Which component is definitely classical?** `A` (alignment), the sharp bound
   `|A| ≤ √(2/3)` (`WolkowiczStyan1980` specialized), and the eigenframe expansion
   `A = Σᵢ(λᵢ/‖S‖_F)cos²θᵢ` (`AshurstEtAl1987`).

**The sharp constant `sup P/‖∇u‖³_F = 4√2/9`** with equality at `ζ = 1/3`,
axisymmetric `S`, `ω` on the distinct eigenvector: not located in the reviewed
literature, but it is an elementary constrained optimization over two classical
ingredients. **Classification: UNCERTAIN — MORE SEARCH NEEDED**, presumption of prior
existence. No priority should be claimed.

**Overall classification: NEW DERIVATION OF KNOWN INGREDIENTS.**

## 5. Prior-art audit: alignment independence

The claim audited: scale, allocation magnitude and comparator do not determine
`P = ω·Sω`, because the relative orientation of `ω` and the strain eigenframe is
independent.

**Finding — this is classical.** It is the founding observation of the alignment
literature: `AshurstEtAl1987` measured the alignment precisely because production is
not determined by magnitudes; `Meneveau2011` and `JohnsonWilczek2024` state the
magnitude/alignment separation as standard background. The whole vortex-identification
debate (`HuntWrayMoin1988`, `ChongPerryCantwell1990`, `JeongHussain1995`,
`LiuEtAl2016`) exists because `Q`-type magnitude measures do not determine the
stretching dynamics. The companion manuscript at
[10.5281/zenodo.22675817](https://doi.org/10.5281/zenodo.22675817) already makes this
its central point.

**Classification: CLASSICAL.** Our contribution is an **exact demonstration inside the
comparator coordinate system** — explicit witness matrices with identical
`(‖S‖_F, |ω|, ζ, h, L, D)` and `P = +2, 0, −2`, and the statement that the joint image
of `(ζ, A)` is a full rectangle. The counterexamples retain value as verification
artifacts; the independence itself must be attributed.

## 6. Prior-art audit: the `ξ` evolution equation

`Dξ/Dt = ‖S‖_F[(1 + e^{2ξ}/2)A − s/√6] + (S:H_dev)/E_S + ν[·]`.

**Finding.** The two ingredients are textbook: the enstrophy equation and the
strain-norm equation (with `−2tr S³`, `−½ω·Sω`, `−2S:H`, viscous) appear in
`NomuraPost1998` and in the reviews. `NomuraPost1998` in particular studies the
coupled dynamics of vorticity **and** rate of strain in exactly this decomposition.
Taking half the difference of their logarithmic rates is one line of calculus. No
prior paper in the reviewed set writes the evolution equation for `log(E_W/E_S)` (or
for `Q/‖∇u‖²`) explicitly in closed form.

**Classification: NEW DERIVATION OF KNOWN INGREDIENTS** — *an immediate combination of
standard equations that is not usually written explicitly*. It is not structurally
distinct, and it should not be described as a new equation of motion.

## 7. Prior-art audit: the pressure-Hessian interpretation

| our statement | status | closest prior work |
|---|---|---|
| the vorticity equation is pressure-blind | **CLASSICAL** (`curl grad p = 0`) | textbook; `Meneveau2011` |
| the strain equation carries `−2S:H` | **CLASSICAL** | `NomuraPost1998`, `Meneveau2011` |
| `H` splits into a **local isotropic** part (fixed by `Δp = E_W−E_S`) and a **nonlocal deviatoric** part | **CLASSICAL** — stated in exactly these terms | `BuariaPumir2023`: "Decomposing **H** into local isotropic (**H**^I) and non-local deviatoric (**H**^D) components"; restricted Euler retains only `H^I` (`Vieillefosse1982`, `Cantwell1992`) |
| the isotropic part drops out of `Dξ/Dt` because `tr S = 0` | elementary consequence of the above | — |
| the obstruction reduces to a **single scalar** `S:H_dev/E_S` | **mostly coordinate compression** | see below |

**On "single scalar obstruction".** Any *scalar* observable can only ever see one
contraction of `H`; reducing to one scalar is a property of having chosen a scalar
observable, not a reduction of the pressure Hessian's own complexity. The genuine
dimensional-reduction result is `CarboneIovienoBragg2020`, which uses a gauge symmetry
to perform a **rank reduction of the anisotropic pressure Hessian itself**, confining
its dynamical activity to two-dimensional manifolds everywhere in the flow — a
substantially stronger statement about the tensor. `YangEtAl2024` further gives an
approximate closed form for `H` in strong-vorticity regions.

**Repository action taken:** the wording "the sharpest positive structural result" was
removed; the claim is now stated as coordinate compression with the classical
local/nonlocal split cited.

## 8. Prior-art audit: the pressure-free functional classification

The audited result: for `J = ∫Φ(E_S,E_W)dx` on the periodic incompressible domain, the
pressure contribution is `−2∫Φ_{E_S} S:H_dev dx`; it vanishes for all fields iff
`Φ_{E_S}` is constant (necessity under a stated two-region nondegeneracy), whence
`Φ = cE_S + ψ(E_W)` and, using `∫E_S = ∫E_W`, `J = ∫ψ̃(E_W)dx`.

**What is classical.** `∫S:H_dev dx = 0` is a first-order homogeneity relation of
exactly the kind this literature systematizes. `CarboneWilczek2022` prove that the
Betchov constraints are the **only** homogeneity constraints for incompressible
isotropic velocity-gradient fields, and state that their analysis "allows the
derivation of homogeneity relations involving the velocity gradient and other
dynamically relevant quantities, **such as the pressure Hessian and viscous
stresses**". `ZhouYang2023` derive and test the mixed-moment relation
`⟨tr(m h^p m)⟩ = −½⟨(tr m²)²⟩`. Our vanishing integral belongs to this family and is
**CLASSICAL**; the audit already described it as the identity behind Betchov's
relation, and now cites the modern classification work.

**What was not found.** A classification of which *functionals of the two sector
energies* have **time evolution** exactly free of the deviatoric pressure Hessian. The
prior results above classify *static moment identities*; ours asks a dynamical
question about a functional's evolution. These are adjacent but not the same
statement, and `CarboneWilczek2022` does not imply ours as far as this review could
determine.

**Classification: NEGATIVE RESULT NOT FOUND IN REVIEW**, with three caveats that must
travel with it: (i) it follows from standard integration-by-parts machinery;
(ii) necessity is conditional on the stated nondegeneracy and is not proved
unconditionally; (iii) it sits beside a published completeness theorem
(`CarboneWilczek2022`) that is more general in its own domain. Required wording:
*"No equivalent classification result was identified in the literature reviewed."*
Never "the first proof".

## 9. Conditional / regional averages

The corollary that non-constant spatial weighting destroys the cancellation, leaving
`∫w S:H dx = ∫p[S:∇∇w + ∇w·Δu]dx` and, for an indicator weight, a boundary term.

**Finding: standard integration-by-parts machinery.** Conditional pressure-Hessian
statistics are routine in this literature (`BuariaPumir2023` conditions on intense
vorticity; `YangEtAl2024` on strong-vorticity regions), and the loss of homogeneity
identities under conditioning is well understood — it is why those identities are
stated for homogeneous averages. **Classification: CLASSICAL machinery; the explicit
weight identity is a useful corollary inside our framework, attributable to standard
calculus rather than to this work.**

---

## 10. Novelty matrix

| # | claim / result | our statement | closest prior literature | relationship | classification | recommended public wording | keys |
|---|---|---|---|---|---|---|---|
| 1 | `ζ` as strain–rotation coordinate | `ζ = (E_W−E_S)/‖∇u‖²` | `Q` criterion; Ω method | affine rescaling of both | **KNOWN / REPARAMETERIZED** | "an affine rescaling of the normalized second invariant, exactly `2Ω_Liu − 1`" | `HuntWrayMoin1988`, `LiuEtAl2016` |
| 2 | `ξ = artanh ζ` | `½log(E_W/E_S)` | same as above | monotone bijection ⇒ same information | **KNOWN / REPARAMETERIZED** | "a logarithmic reparameterization; the logit of the Ω measure. Presentational, not informational" | `LiuEtAl2016` |
| 3 | altitude family is sign-blind; comparator restores | parity in `ζ` | signed standard criteria | property of an imported non-standard coordinate | **KNOWN / REPARAMETERIZED** | "an internal consistency check: the standard invariants are signed for this reason" | `HuntWrayMoin1988`, `ChongPerryCantwell1990` |
| 4 | production factorization | `P = ‖∇u‖³ g(ζ) A` | magnitude × alignment framing | one-line substitution into definitions | **NEW DERIVATION OF KNOWN INGREDIENTS** | "an exact but elementary rearrangement; the alignment content is classical" | `AshurstEtAl1987`, `Meneveau2011` |
| 5 | sharp bound `4√2/9` | `ω·Sω ≤ (4√2/9)‖∇u‖³` | Wolkowicz–Styan specialization | elementary optimization over classical ingredients | **UNCERTAIN — MORE SEARCH NEEDED** | "elementary; we did not locate it, and claim no priority" | `WolkowiczStyan1980` |
| 6 | `|A| ≤ √(2/3)` | sharp alignment bound | specialization of an eigenvalue–trace bound | identical | **CLASSICAL** | cite Wolkowicz–Styan | `WolkowiczStyan1980` |
| 7 | alignment independence | `(ζ,A)` image is a rectangle | founding observation of alignment studies | exact demonstration of a classical fact | **CLASSICAL** (+ new witnesses) | "an exact demonstration, in these coordinates, of a classical independence" | `AshurstEtAl1987`, `JohnsonWilczek2024` |
| 8 | `Dξ/Dt` closed form | see §6 | enstrophy + strain-norm equations | half the difference of two standard log-rates | **NEW DERIVATION OF KNOWN INGREDIENTS** | "an immediate combination of standard sector equations, not usually written explicitly" | `NomuraPost1998`, `Meneveau2011` |
| 9 | vorticity sector pressure-blind | `curl grad p = 0` | textbook | identical | **CLASSICAL** | state as classical | `Meneveau2011` |
| 10 | isotropic `H` local, deviatoric nonlocal | `tr H = Δp = Qζ` | stated in these terms | identical | **CLASSICAL** | cite Buaria & Pumir | `BuariaPumir2023` |
| 11 | "single scalar obstruction" | `S:H_dev/E_S` | rank reduction of anisotropic `H` | ours is coordinate compression; theirs reduces the tensor | **KNOWN / REPARAMETERIZED** | "the contraction a scalar observable can see; not a reduction of the pressure Hessian" | `CarboneIovienoBragg2020` |
| 12 | `∫S:H_dev dx = 0` | integrated cancellation | homogeneity-constraint family | member of a classified family | **CLASSICAL** | cite Betchov and the modern classification | `Betchov1956`, `CarboneWilczek2022`, `ZhouYang2023` |
| 13 | pressure-free ⇒ enstrophy-only | `Φ = cE_S + ψ(E_W)` | completeness of Betchov constraints | adjacent, different question (dynamical vs. moment identity) | **NEGATIVE RESULT NOT FOUND IN REVIEW** | "no equivalent classification result was identified in the literature reviewed; conditional on the stated nondegeneracy" | `CarboneWilczek2022`, `ZhouYang2023` |
| 14 | weighted/regional obstruction | `∫wS:H = ∫p[S:∇∇w+∇w·Δu]` | standard by-parts | standard machinery | **CLASSICAL** | "standard integration by parts, applied in our coordinates" | — |
| 15 | apex occupancy is kinematic | `⟨E_S⟩ = ⟨E_W⟩` for homogeneous flow | classical homogeneity identity | identical | **CLASSICAL** | cite as a homogeneity identity | `Betchov1956`, `CarboneWilczek2022` |

## 11. Top-five audit

**1. Symmetric geometry loses branch information; comparator restores it.**
*Classical:* that the strain/rotation discriminant must be signed — every standard
criterion is. *Coordinate change:* the entire statement, since `h` is an imported
non-standard coordinate and the "comparator" is `sign Q`. *Distinct:* nothing.
*Closest prior:* `HuntWrayMoin1988`, `LiuEtAl2016`. *Repository says:* an internal
consistency check on the imported geometry.

**2. Production separates into scale × magnitude × comparator × alignment.**
*Classical:* the magnitude/alignment separation and `A` itself. *Coordinate change:*
the four-way labelling. *Distinct:* the specific algebraic grouping and the constant
`4√2/9` — elementary, no priority claimed. *Closest prior:* `AshurstEtAl1987`,
`WolkowiczStyan1980`. *Repository says:* an exact but elementary rearrangement.

**3. Alignment remains independent of allocation/comparator.**
*Classical:* the independence itself. *Coordinate change:* expressing it in `(ζ,A)`.
*Distinct:* the exact witnesses and the rectangle statement, as verification.
*Closest prior:* `AshurstEtAl1987`. *Repository says:* an exact demonstration of a
classical independence, not a new degree of freedom.

**4. `ξ` evolution isolates pressure into `S:H_dev/E_S`.**
*Classical:* both sector equations; the pressure-blind vorticity sector; the
local/nonlocal split of `H`. *Coordinate change:* the isolation into one scalar.
*Distinct:* the closed form itself, as an explicit combination. *Closest prior:*
`NomuraPost1998`, `BuariaPumir2023`, `CarboneIovienoBragg2020`. *Repository says:* an
immediate combination of standard equations; the "single scalar" is coordinate
compression.

**5. Exact pressure-freedom removes comparator dependence.**
*Classical:* `∫S:H_dev = 0` and its homogeneity-constraint context. *Coordinate
change:* none essential. *Distinct:* the classification of the functional class — the
only result surviving this review as apparently distinct. *Closest prior:*
`CarboneWilczek2022` (completeness for moment constraints), `ZhouYang2023`.
*Repository says:* no equivalent classification result was identified in the
literature reviewed; conditional on the stated nondegeneracy; adjacent to a stronger
published completeness theorem in a different domain.

## 12. Terminology audit

| our term | collides with standard usage? | standard equivalent | action |
|---|---|---|---|
| **comparator** | no collision found | the signed strain–rotation discriminant; `sign Q` | keep, but define on first use as "the signed allocation, i.e. the sign/value of the normalized second invariant" |
| **rapidity** | no turbulence meaning; borrowed from relativity | *log strain–rotation ratio*, `½log(E_W/E_S)`; the logit of `Ω_Liu` | keep as an internal geometric analogy; **always give the standard name alongside** |
| **allocation** | no collision found | relative strain/rotation contribution to `‖∇u‖²_F`; `Ω_Liu` | keep with definition |
| **Thales quotient** | none — internal to the companion framework | — | mark explicitly as internal terminology, not turbulence usage |
| **pressure obstruction** | no collision, but risks sounding like a theorem | "the nonlocal (deviatoric) pressure-Hessian contribution" | keep informally; use the standard phrase in claim statements |

No renaming of the repository is warranted; no serious scientific ambiguity was found.

## 13. Sources that would have strengthened the deposited paper

Identified **during this subsequent repository literature review**, not present in the
deposited manuscript at [10.5281/zenodo.22675817](https://doi.org/10.5281/zenodo.22675817)
as far as its reference list shows. These are recorded here for future work; the
Zenodo record is not altered, and no claim is made that they were cited there.

* `CarboneIovienoBragg2020` — the rank reduction of the anisotropic pressure Hessian;
  directly relevant to any claim about compressing the pressure contribution.
* `CarboneWilczek2022` — completeness of the Betchov homogeneity constraints; the
  natural home for any classification statement about pressure-related identities.
* `ZhouYang2023` — mixed-moment homogeneity constraints involving the pressure Hessian.
* `BuariaPumir2023` — the local-isotropic / nonlocal-deviatoric decomposition of `H`
  and its opposite effects on vortex stretching.
* `YangEtAl2024` — pressure-Hessian structure in strong-vorticity regions.
* `JohnsonWilczek2024` — the current review, which supersedes `Meneveau2011` as the
  orientation reference.

## 14. Mathematical corrections triggered by the literature

**None.** No reviewed source contradicts a computed result in this repository. The
changes made are to *framing and attribution*, not to mathematics: the 179 symbolic
checks are unaffected, and no equation was altered.

## 15. Standalone-paper assessment

Recorded in [`../paper/README.md`](../paper/README.md). Summary: results 1–4 are
classical, reparameterizations, or elementary recombinations of standard equations;
result 5 is the only candidate and is conditional, derivable by standard machinery,
and adjacent to a stronger published theorem. **Recommendation: C — repository note
only.**
