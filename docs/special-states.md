# Special states of the comparator equation

> *Reproduced verbatim from the audit run that produced it; only file paths have
> been updated to this repository's layout. Statements of the form "nothing outside
> this directory was modified" refer to that original run.*

All values from `src/xi_numeric.py` (section K), exact to machine precision.
The law being evaluated is

```
D xi/Dt = ||S||_F [ (1 + e^{2 xi}/2) A - s/sqrt(6) ]   +   S:H_dev / ||S||_F^2   +   nu[...]
          \___________________ local ___________________/   \____ nonlocal ____/
```

with every state normalized to `||grad u||_F = 1` and `nu = 0` (the viscous pair needs
second derivatives and is not a function of the state).

## Interior states

| # | state | `zeta` | `xi` | `A` | `s` | local | pressure | `B` |
|---|---|---|---|---|---|---|---|---|
| 1 | balanced, `omega` on `e_max` | 0 | 0 | +0.7071 | 0 | **+0.7500** | 0 | – |
| 2 | balanced, `omega` on `e_int` | 0 | 0 | 0 | 0 | **0** | 0 | – |
| 3 | balanced, `omega` on `e_min` | 0 | 0 | −0.7071 | 0 | **−0.7500** | 0 | – |
| 4 | strain-rich, `zeta = −0.8` | −0.8 | −1.0986 | +0.5 | −1 | **+0.8880** | 0 | – |
| 5 | rotation-rich, `zeta = +0.8` | +0.8 | +1.0986 | +0.5 | −1 | **+0.9987** | 0 | – |
| 6 | balanced, `H_dev = 0` (restricted Euler) | 0 | 0 | +0.5 | −1 | +0.8190 | **0** | – |
| 7 | balanced, `H_dev` aligned with `S` | 0 | 0 | +0.5 | −1 | +0.8190 | **+1.4142** | +1 |
| 8 | balanced, `H_dev` eigen-permuted | 0 | 0 | +0.5 | −1 | +0.8190 | **−0.7071** | −0.5 |

Rows 1–3 use the non-degenerate strain state `lambda/||S||_F = (0.7071, 0, −0.7071)`
(`s = 0`), so "intermediate eigenvector" is meaningful; rows 4–8 use the axisymmetric
state `(2,−1,−1)/sqrt6` (`s = −1`).

Readings:

* **1 vs 3** — reversing the alignment reverses the sign of the comparator's motion:
  vorticity on the extensional axis drives the state toward rotation dominance,
  on the compressive axis toward strain dominance.
* **2** — vorticity on the intermediate axis of a state with `lambda_int = 0` and
  `s = 0` gives `D xi/Dt = 0` exactly: a stationary point of the comparator that is
  *not* an equilibrium of the flow.
* **4 vs 5** — at equal `|zeta|` and equal `A`, the rotation-rich branch moves faster
  (`+0.9987` vs `+0.8880`), because the stretching coefficient `1 + e^{2 xi}/2` grows
  with `xi`. The comparator accelerates its own drift toward rotation dominance.
* **6 vs 7 vs 8** — identical local state, three different nonlocal environments. The
  pressure contribution ranges over `[−0.7071, +1.4142]` while the local part is
  frozen at `+0.8190`: in row 8 the *sign* of `D xi/Dt` is still positive, but a
  moderately larger misaligned `H_dev` would reverse it.

## Endpoint limits (analytic, not evaluated numerically)

| limit | behaviour of `xi` | behaviour of `D xi/Dt` | status |
|---|---|---|---|
| `zeta -> −1` (pure strain, `E_W -> 0`) | `xi -> −infinity` | `-> \|\|S\|\|[A − s/sqrt6] + S:H_dev/E_S`, **finite** along any fixed `omega` direction | coordinate singularity only |
| `zeta -> +1` (pure rotation, `E_S -> 0`) | `xi -> +infinity` | `\|\|S\|\| e^{2 xi} A = (E_W/\|\|S\|\|) A -> infinity` unless `A -> 0` | **genuine divergence** |
| `zeta = 0` (balance) | `xi = 0` | regular; `d xi/d zeta = 1` | regular point (the fold of `h` is *not* a fold of `xi`) |

The asymmetry between the two endpoints is structural, not a convention: the
rotation sector is forced only by `P` (which vanishes with `E_W`), whereas the strain
sector is forced by `−P/2` and `−2 tr(S^3)`, which do **not** vanish as `E_S -> 0`.
A vanishing strain sector is therefore regenerated at a finite rate from an
infinitesimal base, and `log E_S` moves infinitely fast.

## Conditioning (section L)

| `zeta` | `d xi/d zeta = R = 1/(4h^2)` | `d h/d zeta` |
|---|---|---|
| 0.00 | 1.0000 | 0 (fold) |
| 0.50 | 1.3333 | −0.2887 |
| 0.90 | 5.2632 | −1.0324 |
| 0.99 | 50.2513 | −3.5090 |

`xi` regularizes exactly where `h` degenerates, and stretches exactly where `h` is
well-behaved. The Jacobian is the framework's own response coordinate `R`.
