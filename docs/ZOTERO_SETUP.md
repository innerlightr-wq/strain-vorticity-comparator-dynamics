# Zotero: installation, library structure, and how to rebuild it

**Status as of 18 September 2026: Zotero is installed, running, and the library for this
project is built.** This file previously recorded that Zotero was unavailable; that was
accurate when the literature review was carried out and is now superseded. The review
itself was not redone — the bibliography in [`../references.bib`](../references.bib) was
already verified against publisher metadata, and Zotero confirmed it (see §4).

## 1. Installation

| | |
|---|---|
| executable | `~/Downloads/Zotero-10.0.3_linux-x86_64 (1)/Zotero_linux-x86_64/zotero` (Bash launcher → `zotero-bin`, ELF x86-64) |
| installation type | extracted upstream Linux **archive** — not a system package, AppImage, Flatpak or Snap |
| version | Zotero **10.0.3**, API version 3, schema version 44 |
| profile | `~/.zotero/zotero/5lik04lk.default` (preferences in `prefs.js`) |
| data directory | `~/Zotero/` — the database is `~/Zotero/zotero.sqlite` |

`zotero.sqlite` was never edited directly. Every change below went through Zotero's own
HTTP interface on `127.0.0.1:23119`.

## 2. Enabling the local API

The local API is off by default and returns `403` until the preference is set:

```
extensions.zotero.httpServer.localAPI.enabled = true
```

The Zotero helper tooling sets it and restarts the application:

```bash
Z=~/.codex/plugins/cache/openai-curated-remote/zotero/0.1.2/skills/zotero
python3 $Z/scripts/zotero.py enable --restart
python3 $Z/scripts/zotero.py status --json      # api_running: true
```

Two practical notes. The helper keeps a timestamped backup of `prefs.js` before editing
it. Its restart step launches a bare `zotero`, so with an extracted archive the binary
must be on `PATH` or Zotero will shut down without coming back — either add the directory
to `PATH` or relaunch it yourself:

```bash
setsid "$HOME/Downloads/Zotero-10.0.3_linux-x86_64 (1)/Zotero_linux-x86_64/zotero" &
```

## 3. Reads and writes

**Reads** need no credentials: `GET http://127.0.0.1:23119/api/users/0/...` implements
Zotero Web API v3 for the logged-in desktop user.

**Writes** are supported by the local API in Zotero 10, in three steps. They are *not*
available through the helper script, which can only list collections and import through
the Connector into whichever collection is selected in the UI.

1. Read the server ID from the response headers of `GET /api/`
   (`Zotero-Server-ID: <id>`). Without it a write returns `428`.
2. `POST /api/local/authorize` with `{"appName":"<name>"}` and that header. **This blocks
   on a "Local API Authorization" dialog in the Zotero window until a human approves it.**
   The response is `{"key":"…","remember":true}`; because the approval is remembered,
   later calls for the same `appName` return immediately without a dialog.
3. Send the key as `Zotero-API-Key: <key>` (or `Authorization: Bearer <key>`) on
   `POST`/`PATCH`/`DELETE`. Keys go stale: fetch one and use it in the same session, or a
   write returns `401 Invalid or expired API key`. Item edits take
   `If-Unmodified-Since-Version: <version>` so a concurrent change fails the write rather
   than overwriting it.

## 4. What the library contains

25 items, matching [`../references.bib`](../references.bib) one-for-one, imported once
through the Connector (`import-bibtex --file references.bib`) into an empty library:
23 journal articles, 1 conference paper (`HuntWrayMoin1988`), 1 book section
(`Burgers1948`). Trash is empty.

**Metadata cross-check: no discrepancies.** All 25 items matched by title, and DOI, year,
volume, pages and journal agree with the verified BibTeX in every case — 0/25 differences.
Zotero confirmed the records rather than correcting them, so no citation in the repository
required revision.

Collection `Strain–Vorticity Comparator Dynamics` (key `LMCGVKDK`) holds ten
subcollections. Items live in the subcollections, not in the parent:

| subcollection | items | keys |
|---|---|---|
| 00 — Reviews & Orientation | 2 | `Meneveau2011`, `JohnsonWilczek2024` |
| 01 — Betchov & Exact Identities | 4 | `Betchov1956`, `WolkowiczStyan1980`, `CarboneWilczek2022`, `ZhouYang2023` |
| 02 — Vorticity–Strain Alignment | 3 | `AshurstEtAl1987`, `NomuraPost1998`, `LundRogers1994` |
| 03 — Restricted Euler & Velocity-Gradient Dynamics | 7 | `Vieillefosse1982`, `Vieillefosse1984`, `Cantwell1992`, `Burgers1948`, `BealeKatoMajda1984`, `Meneveau2011`, `JohnsonWilczek2024` |
| 04 — Pressure Hessian & Nonlocality | 5 | `WilczekMeneveau2014`, `VlaykovWilczek2019`, `CarboneIovienoBragg2020`, `BuariaPumir2023`, `YangEtAl2024` |
| 05 — Strain / Rotation Coordinates | 6 | `HuntWrayMoin1988`, `ChongPerryCantwell1990`, `JeongHussain1995`, `LiuEtAl2016`, `Okubo1970`, `Weiss1991` |
| 06 — Comparator / Rapidity Prior-Art Search | 5 | `LiuEtAl2016`, `HuntWrayMoin1988`, `ChongPerryCantwell1990`, `Okubo1970`, `Weiss1991` |
| 07 — Integral & Pressure-Cancellation Identities | 3 | `Betchov1956`, `CarboneWilczek2022`, `ZhouYang2023` |
| 08 — Directly Cited in Repository | 25 | all |
| 09 — Novelty Uncertain / Read Closely | 5 | `CarboneWilczek2022`, `CarboneIovienoBragg2020`, `ZhouYang2023`, `WolkowiczStyan1980`, `LiuEtAl2016` |

**65 memberships across 25 distinct items, with no duplicates.** An item belongs to every
subcollection that applies — `LiuEtAl2016` is in four — which is why the counts sum past
25. Every item is in at least one subcollection.

Subcollection 09 marks the items the novelty matrix depends on most: the two that
establish `ζ`/`ξ` as a published measure and a classical bound, and the three whose
relationship to the §7 classification determines how that result must be worded. Read
those closely before restating any novelty claim.

## 5. Rebuilding this from scratch

With Zotero running and the local API enabled:

1. Import the bibliography once, into an empty library or a selected collection:
   `python3 $Z/scripts/zotero.py import-bibtex --file references.bib --yes`.
2. Obtain a write key as in §3.
3. `POST /api/users/0/collections` with `[{"name":"Strain–Vorticity Comparator Dynamics"}]`,
   then again with `[{"name":"00 — …","parentCollection":"<parent key>"}, …]` for the ten
   subcollections (one request takes the whole batch).
4. For each item, `PATCH /api/users/0/items/<key>` with `{"collections":[…]}` and its
   `If-Unmodified-Since-Version`. Patching membership **adds** the item to further
   collections; it does not copy the item, which is what keeps the count at 25.

Do not re-import per subcollection: each Connector import creates new items, so importing
the ten subsets would produce 65 items instead of 25 memberships of the same 25.

## 6. Tags and notes

Not applied. The per-claim commentary that Zotero item notes would carry lives in
[`LITERATURE_REVIEW.md`](LITERATURE_REVIEW.md) — §1 states what each item is cited for,
§§2–9 hold the prior-art comparison, and §10 is the novelty matrix. To mirror the matrix
in Zotero, the useful tags are `classical`, `known-reparameterized`,
`new-derivation-known-ingredients`, `apparently-distinct`, `negative-result-not-found`,
`uncertain-more-search`, plus `closest-prior-art` on the five items in subcollection 09
and on `BuariaPumir2023` and `AshurstEtAl1987`.

## 7. PDFs

No full texts were downloaded, attached, or redistributed, and the repository contains no
copies of any third-party paper. Claim-level comparison was done from publisher-hosted
abstracts and openly readable statements.
