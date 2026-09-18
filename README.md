# Strain–Vorticity Comparator Dynamics

Symbolic verification, numerical validation and adversarial audit of a coordinate
framework for the local strain–rotation balance of incompressible three-dimensional
flow, carried out in four stages. Each stage was run as a falsification attempt on the
stage before it. Two of the four returned negative results, and one retracts a
recommendation made by an earlier stage of the same audit.

**Author:** Elias De Jesús — Independent Researcher
([ORCID 0009-0007-0190-9143](https://orcid.org/0009-0007-0190-9143))

> **This repository does not claim a proof of Navier–Stokes global regularity.**
> It contains no new a-priori estimate, inequality, coercive quantity, monotonicity
> result, or closure theorem, and nothing here should be read as progress toward one.

## Relationship to the published research

This repository **extends the strain–vorticity research program** associated with the
published paper

> *Strain–Vorticity Interaction and Rotational Coherence: From Magnitude Comparison to
> Enstrophy Production* — [10.5281/zenodo.22675817](https://doi.org/10.5281/zenodo.22675817)

That DOI is the scholarly anchor for the **paper**. It is **not** a DOI for this
repository; this repository has no DOI of its own. The paper supplies the exact
interaction term `P = ω·Sω`, the alignment coordinate `A`, the sharp bound
`|A| ≤ √(2/3)`, the canonical-flow analysis, and the conclusion that no universal
scalar threshold follows from the enstrophy balance. A second manuscript — *Archimedean
Compensation on the Thales Semicircle* (Revision 3) — supplies the binary-partition
geometry. Reading notes on both, including which of their claims this audit supports,
sharpens, or corrects, are in [`docs/strain-paper-notes.md`](docs/strain-paper-notes.md).

## How to read the claims in this repository

Results are kept in separate categories on purpose, and are labelled individually
rather than by section:

| category | what it means | where |
|---|---|---|
| **published-paper results** | established in the Zenodo paper above; used here as input, not re-derived as new | `docs/strain-paper-notes.md` |
| **later exact identities** | algebra derived in this audit, elementary, no priority claimed | `docs/derivations.md`, `docs/pressure-hessian-results.md` |
| **symbolic verification** | machine-checked on symbolic input with free entries (179 checks) | `src/*_exact.py`, `results/*_exact.json` |
| **numerical validation** | finite-precision confirmation on random admissible states | `src/*_numeric.py`, `experiments/synthetic_fields/` |
| **DNS observations** | measured in a 96³ unforced Taylor–Green simulation; empirical, single flow | `experiments/dns_checks/`, `reports/` |
| **negative / impossibility results** | the audit's strongest output; stated with their scope conditions | `reports/weighted-comparator-audit.md`, `docs/failures-and-counterexamples.md` |
| **open questions** | explicitly not settled here | end of each report; `paper/README.md` |

Epistemic labels used throughout: **EXACT** · **COMPUTATIONAL** · **EMPIRICAL** ·
**CONJECTURE** · **FAILED**.

## The structural hierarchy

The audit's organising picture, in the order in which information is added:

```
        scale
      + allocation
      + comparator
      + alignment
      + nonlocal pressure environment
```

* **scale** — `‖∇u‖_F`;
* **allocation** — the symmetric magnitude split `ζ` (equivalently `h`, `D`, `R`, `η`);
* **comparator** — the signed orientation of that split, `σ = sign ζ` or the rapidity
  `ξ = artanh ζ`, which is what the symmetric altitude family discards;
* **alignment** — `A = ω·Sω/(‖S‖_F|ω|²)`, the classical vortex-stretching alignment;
* **nonlocal pressure environment** — the single deviatoric scalar `S:H_dev/E_S`,
  which is not a coordinate of the fluid element at all.

The first four are exactly sufficient for the *instantaneous* production `P`. The
fifth is required for its *evolution*.

**Interpretation, stated as an interpretation and not as a theorem:** *a local state
description is not the same thing as closed local dynamics.* What is proved is
narrower and specific — `(‖∇u‖_F, ζ, A)` determines `P` exactly, while `Dξ/Dt`
additionally requires the strain-state invariant `s` and a nonlocal pressure scalar
that no member of the tested functional class can remove. The general reading is a way
of organising those particular results, not a separate claim.

## Conventions (fixed once, in [`src/svcore.py`](src/svcore.py))

```
(∇u)_ij = ∂u_i/∂x_j        S = sym(∇u)        Ω = W = skew(∇u)
ω_i = ε_ijk ∂_j u_k        ‖W‖²_F = |ω|²/2 = E_W       E_S = ‖S‖²_F
Q = ‖∇u‖²_F = E_S + E_W    a = E_S/Q   b = E_W/Q   a + b = 1
ζ = b − a                  h = √(ab)   L = ζ/2   D = ½ − h   ξ = artanh ζ
P = ω·Sω                   A = P/(‖S‖_F|ω|²),  |A| ≤ √(2/3)
s = −3√6 det S/‖S‖³_F      (Lund–Rogers strain state)
```

## The four audits and their verdicts

| # | audit | question | verdict |
|---|---|---|---|
| 1 | [allocation](reports/allocation-audit.md) | does the Thales/Archimedean layer add information to `ζ`? | **A — pure reparameterization** |
| 2 | [comparator](docs/comparator-results.md) | does an antisymmetric comparator repair what the symmetric geometry loses? | **B — restores it, but magnitude+comparator still does not determine the dynamics** |
| 3 | [ξ evolution](reports/xi-evolution-audit.md) | does `Dξ/Dt` isolate local dynamics from the nonlocal pressure Hessian? | **B — clean separation, one irreducible nonlocal scalar** |
| 4 | [weighted comparator](reports/weighted-comparator-audit.md) | can a weighted global functional cancel that scalar exactly? | **E — impossibility within the tested class; branch closed** |

### Later exact identities

```
P = ‖∇u‖³_F · g(ζ) · A ,            g(ζ) = (1+ζ)√((1−ζ)/2) = 2 b √a
  = ‖∇u‖³_F · 2h^{3/2} · e^{ξ/2} · A      (unique even × odd split of g)

ω·Sω ≤ (4√2/9) ‖∇u‖³_F ,  sharp, attained iff ζ = 1/3, S axisymmetric, ω on its
                           distinct eigenvector

Dξ/Dt = ‖S‖_F[(1 + e^{2ξ}/2)A − s/√6] + (S:H_dev)/E_S + ν[W:ΔW/E_W − S:ΔS/E_S]

pressure part of d/dt ∫Φ(E_S,E_W) dx = −2 ∫ Φ_{E_S} S:H_dev dx
```

### The negative result, stated at its safe strength

> Within the **tested local two-sector functional class** `J = ∫Φ(E_S,E_W)dx` on the
> periodic incompressible domain, and **under the stated nondegeneracy condition**,
> exact elimination of the deviatoric pressure-Hessian contribution forces the
> functional to reduce effectively to a **vorticity/enstrophy-only class** — and
> therefore to carry no strain–rotation comparator information.

Sufficiency (`Φ_{E_S}` constant ⟹ pressure-free ⟹ `J = ∫ψ̃(E_W)dx`, using
`∫E_S = ∫E_W`) is exact. Necessity is established under an explicit two-region
nondegeneracy condition and verified for every weight and exponent tested; it is
**not** proved unconditionally for all incompressible flows, and individual symmetric
fields can make the relevant integral vanish accidentally. See
[`reports/weighted-comparator-audit.md`](reports/weighted-comparator-audit.md) §7 and
[`docs/pressure-hessian-results.md`](docs/pressure-hessian-results.md) §5 for the full
statement and its scope.

### What is not new

`ζ` is an affine rescaling of the second invariant `Q_HWM`; `A` is the classical
vortex-stretching alignment (Betchov 1956; Ashurst et al. 1987); `|A| ≤ √(2/3)`
specializes Wolkowicz–Styan; `∫S:H_dev = 0` is the identity behind Betchov's relation;
`∫ψ̃(E_W)` are the standard vorticity `L^p` norms. The Thales coordinates are exact
functions of `ζ`, and `ξ` is informationally equivalent to it. The derived identities
above are elementary, and **no priority is claimed** for any of them.

### Open questions

Listed at the end of each report. The principal ones: whether the sharp constant
`4√2/9` and the closed form of `Dξ/Dt` have prior appearances in the
velocity-gradient literature; whether the `ζ`–`A` statistical coupling seen in the DNS
(0.5–1.2 % of `var(A)`) survives at higher Reynolds number; and whether the
impossibility result is already folklore, since it follows quickly from the
pressure-free vorticity equation.

## Layout

```
docs/      source-paper notes, derivations, canonical flows, counterexamples,
           protocol, environment, special states, comparator and pressure results
src/       coordinate algebra + every exact/symbolic and ensemble-numeric script
experiments/dns_checks/        Navier-Stokes solver and DNS-based measurements
experiments/synthetic_fields/  random-field ensembles, restricted Euler, redistribution
figures/   9 figures, regenerated from results/ by the src/figures_*.py scripts
results/   JSON check logs and summary numbers (large .npz are regenerated, not stored)
reports/   the three main audit reports
paper/     manuscript plan (not written)
```

## Reproducing

```bash
python3 -m venv .venv && .venv/bin/pip install -r requirements.txt

# exact algebra (fast; 179 symbolic checks in total)
.venv/bin/python src/exact_algebra.py        # 44/44   allocation audit
.venv/bin/python src/comparator_exact.py     # 56/56   comparator audit
.venv/bin/python src/xi_exact.py             # 49/49   xi evolution
.venv/bin/python src/weighted_exact.py       # 30/30   weighted comparator

# numerics
.venv/bin/python src/counterexamples.py
.venv/bin/python src/canonical_flows.py
.venv/bin/python src/comparator_numeric.py
.venv/bin/python src/xi_numeric.py
.venv/bin/python src/weighted_numeric.py

# experiments
.venv/bin/python experiments/synthetic_fields/ensembles.py
.venv/bin/python experiments/synthetic_fields/restricted_euler.py
.venv/bin/python experiments/synthetic_fields/synthetic_redistribution.py
.venv/bin/python experiments/dns_checks/dns_taylor_green.py --n 96 --re 400 --tend 12 \
    --snaps 4 6 8 10 --seed-time 7.0 --window 1.5 --n-part 40000
.venv/bin/python experiments/dns_checks/dns_reanalysis.py
.venv/bin/python experiments/dns_checks/homogeneity_check.py
.venv/bin/python experiments/dns_checks/xi_dns_pressure.py --n 96 --re 400 --t 8.0

# figures (need the .npz sample files produced by the runs above)
.venv/bin/python src/figures_allocation.py
.venv/bin/python src/figures_comparator.py
.venv/bin/python src/figures_xi.py
```

Every script exits non-zero if any of its own checks fail and writes its numbers to
`results/`. Total runtime is roughly 45 minutes, dominated by the two 96³ DNS runs.
Machine, versions and timings: [`docs/environment.md`](docs/environment.md).

## Status and non-claims

A verification and falsification record, not a physics result. It introduces no new
mechanism and modifies no established theory. Where a quantity is a reparameterization
of a known one, it is labelled as such; where a near-coincidence with a framework
landmark occurs (`b = 2/3` versus the cubic root `0.68233`), it is reported as a
**non**-identity with the exact residual `−1/27`. And, to repeat the point at the top:
**no proof of Navier–Stokes global regularity is claimed here, in whole or in part.**

## License

Dual-licensed, by material type:

* **Code** (`src/`, `experiments/`, and any scripts) — **Apache License 2.0**, see
  [`LICENSE`](LICENSE).
* **Original research prose, documentation, diagrams and figures** (`docs/`,
  `reports/`, `paper/`, `figures/`, this README) — **CC BY 4.0**, see
  [`LICENSE-CONTENT.md`](LICENSE-CONTENT.md).

Third-party material retains its original license, and these terms do not change the
license of any separately published Zenodo deposit. Details in
[`LICENSE-CONTENT.md`](LICENSE-CONTENT.md).

## Citation

See [`CITATION.cff`](CITATION.cff). If you use the results, cite the published paper
([10.5281/zenodo.22675817](https://doi.org/10.5281/zenodo.22675817)); if you refer to
the audit itself, cite this repository.
