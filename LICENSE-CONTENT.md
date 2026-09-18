# Content license

This repository is licensed in two parts.

## Software and code — Apache License 2.0

All source code in this repository — everything under `src/` and `experiments/`, and
any scripts elsewhere — is licensed under the **Apache License, Version 2.0**. The
full text is in [`LICENSE`](LICENSE).

```
SPDX-License-Identifier: Apache-2.0
```

## Original research prose, documentation, diagrams and figures — CC BY 4.0

Unless otherwise noted in the file itself, all original written content in this
repository — the research prose, documentation, derivations, tables, diagrams and
figures, specifically everything under `docs/`, `reports/`, `paper/`, `figures/`, and
`README.md` — is licensed under the
**Creative Commons Attribution 4.0 International License (CC BY 4.0)**.

```
SPDX-License-Identifier: CC-BY-4.0
```

Full text: <https://creativecommons.org/licenses/by/4.0/legalcode>

To comply, give appropriate credit (author, title, link to the license), indicate
whether changes were made, and do not suggest the licensor endorses you or your use.

## Numerical results

The machine-generated outputs in `results/` (JSON check logs, summary numbers,
`.npz` sample arrays) are released under the same terms as the code
(**Apache-2.0**), since they are direct outputs of the scripts rather than authored
prose.

## Third-party material

**Third-party material retains its original license.** Nothing in this repository
relicenses work that is not the author's own. In particular:

* the Python dependencies listed in [`requirements.txt`](requirements.txt) (NumPy,
  SciPy, SymPy, Matplotlib) are separately licensed by their respective projects and
  are neither vendored nor redistributed here;
* the classical results relied upon and cited throughout — the Betchov relation, the
  Ashurst et al. alignment statistics, the Wolkowicz–Styan eigenvalue bound, the
  Hunt–Wray–Moin second invariant, the Lund–Rogers strain-state parameter, the
  Beale–Kato–Majda criterion — belong to their authors and publishers and are cited,
  not reproduced;
* any quoted phrasing from the source manuscripts is quoted for scholarly comment
  under the ordinary conventions of citation.

If any file in this repository is subject to a different license, that license is
stated in or alongside the file, and it takes precedence over this document.

## Relationship to the published paper

This repository accompanies, and audits, published research deposited separately on
Zenodo (see [`CITATION.cff`](CITATION.cff)), including the strain–vorticity paper
[10.5281/zenodo.22675817](https://doi.org/10.5281/zenodo.22675817).

**The licenses stated here apply to the contents of this repository only. They do not
change, restate, or retroactively affect the license under which any Zenodo deposit
was published.** The license of each deposited paper is whatever that deposit
records, and that record is authoritative for the paper. Where this repository
reproduces or paraphrases material from those deposits, it does so as a citing work.

## Attribution

Elias De Jesús — Independent Researcher —
[ORCID 0009-0007-0190-9143](https://orcid.org/0009-0007-0190-9143)

Suggested attribution for the content portion:

> De Jesús, E. (2026). *Strain–Vorticity Comparator Dynamics.*
> https://github.com/innerlightr-wq/strain-vorticity-comparator-dynamics
> Licensed under CC BY 4.0.
