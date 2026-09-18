# Zotero: status and reconstruction instructions

## Status: Zotero is not available on this machine

The literature review was specified to be carried out inside Zotero — a library
collection, subcollections, tags, per-item notes, and a Better BibTeX export. **It
could not be:**

| checked | result |
|---|---|
| `zotero` / `zotero6` / `zotero7` on `PATH` | not found |
| `~/Zotero`, `~/.zotero`, `~/snap/zotero*`, `~/.var/app/org.zotero.Zotero` | do not exist |
| `zotero.sqlite` anywhere under `$HOME` | not found |
| local Zotero HTTP API / connector on `127.0.0.1:23119` | no listener |
| Zotero MCP server or Zotero tooling in this session | none available |

Nothing was installed, and no Zotero account, web API key, or cloud library was
configured, so the following specified steps **were not executed**: creating the
library collection and subcollections, adding items through Zotero, retrieving PDFs
through Zotero, tagging items in Zotero, writing Zotero item notes, and exporting
`references.bib` via Better BibTeX.

**What was done instead.** Each source was identified and verified individually
against publisher metadata by DOI content negotiation
(`curl -LH 'Accept: application/x-bibtex' https://doi.org/<DOI>`) and, for
title-based lookups, the Crossref REST API. The resulting
[`../references.bib`](../references.bib) is standard BibTeX and imports into Zotero
directly. The organisation that would have lived in Zotero collections and tags lives
in [`LITERATURE_REVIEW.md`](LITERATURE_REVIEW.md) instead: §1 states what each item is
cited for, and §§2–9 hold the per-claim comparison that Zotero item notes were meant
to carry.

This is reported plainly rather than simulated. No Zotero export was fabricated.

## Reconstructing the intended library in one import

1. Zotero → *File* → *Import…* → *A file* → select `references.bib`. All 25 items
   import with DOIs; Zotero will resolve remaining metadata from the DOIs.
2. Create the collection **Strain–Vorticity Comparator Dynamics** with these
   subcollections, and file the keys as listed:

   | subcollection | items |
   |---|---|
   | *Velocity-gradient dynamics* | `Vieillefosse1982`, `Vieillefosse1984`, `Cantwell1992`, `Meneveau2011`, `JohnsonWilczek2024` |
   | *Strain–rotation criteria* | `HuntWrayMoin1988`, `ChongPerryCantwell1990`, `JeongHussain1995`, `LiuEtAl2016`, `Okubo1970`, `Weiss1991` |
   | *Alignment and production* | `AshurstEtAl1987`, `NomuraPost1998`, `LundRogers1994`, `WolkowiczStyan1980` |
   | *Pressure Hessian* | `WilczekMeneveau2014`, `VlaykovWilczek2019`, `CarboneIovienoBragg2020`, `BuariaPumir2023`, `YangEtAl2024` |
   | *Homogeneity constraints* | `Betchov1956`, `CarboneWilczek2022`, `ZhouYang2023` |
   | *Regularity and exact solutions* | `BealeKatoMajda1984`, `Burgers1948` |

3. Apply these tags, which are the ones the review actually uses:
   `classical`, `known-reparameterized`, `new-derivation-known-ingredients`,
   `apparently-distinct`, `negative-result-not-found`, `uncertain-more-search`,
   plus `closest-prior-art` on `CarboneIovienoBragg2020`, `CarboneWilczek2022`,
   `ZhouYang2023`, `BuariaPumir2023`, `LiuEtAl2016`, `AshurstEtAl1987`.
4. For item notes, paste the relevant row of the novelty matrix
   ([`LITERATURE_REVIEW.md`](LITERATURE_REVIEW.md) §10) onto each item.

## PDFs

No full texts were downloaded, attached, or redistributed. All claim-level comparison
was done from publisher-hosted abstracts and openly readable statements, and the
repository contains no copies of any third-party paper.
