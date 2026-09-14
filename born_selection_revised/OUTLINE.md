# Born selection, revised: outline of the specification paper

**Opened 2026-09-14 at JB's instruction, after the internal review round of the Dirac–Kuramoto
record's Paper 1 v0.10 and the consistency pass (v0.11).** The decision, JB's: a clean draft in
MathX, of the second of three shapes offered — a *specification* paper. Not a model paper (the
model is the photodetection unravelling with JB's vocabulary, which the sibling paper already
states) and not a bare negative-result paper (the record already carries that). The claim is what
any microscopic selection mechanism must supply to reproduce Born's rule and one record per
quantum, each clause with a registered run that shows what happens without it, and the model of
`MODEL_SECTION.md` §12 as the worked example that meets every clause by inheriting it.

## Provenance

The three ideas of JB's that shaped the runs (README, 2026-09-14): the vacuum interacting with
the individual atoms of the bulk; the non-commutative order of vacuum, bound electron and free
particle; many atoms in the aperture. The formalization, the runs and the prose are Claude's. The
Dirac–Kuramoto record (`rayolddog/DiracKuramotoFramework`) is cited as the record, never edited.

## The clauses (one section each, one run each)

| clause | what must hold | the run that shows what happens otherwise | where |
|---|---|---|---|
| 1 intensity | the per-site hazard is proportional to the local intensity \|ψ\|², and a term that enters by interference with the amplitude gives √Born | vacuum race V1 against V2 | `vacuum_race/` Table 1 |
| 2 proportion, not maximum | selection is a first-passage race, not a maximum over sites; a maximum gives the brighter port every record unless the noise is Gumbel at unit scale, which is the race in disguise | vacuum race V0, V3–V5 against V1; the β sweep | `vacuum_race/` Tables 1, 3 |
| 3 many sites | the phase factor of the deposit averages out only from about ten sites per port; at one site the model predicts a 12 % departure | second study Table 2; vacuum race N = 1 | `second_spinor_phase_race/`, `vacuum_race/` |
| 4 memorylessness | the clock is renewed per arrival; a frozen factor is harmless only if light-tailed; the Fano factor of the split reads a correlation time in the model but has no power in a device | frozen/fresh addenda 1 and 3 | `vacuum_race/` Tables 6, 7, 9 |
| 5 the conditional state | under staggered arrival the port share is the branch weight; rates on the survival of the whole quantum give the nearer port; a local repair gives 0.952; only the normalized conditional state gives Born; P5 acts twice | staggered addenda 2, 2b | `vacuum_race/` Table 8 |
| 6 fairness in the exchange class | in a closed pot any symmetric zero-mean kernel with absorbing zero is fair; the √e scaling is neither necessary nor unique; a reflecting zero biases toward the poor from three sites up; a threshold biases toward the rich; a rate-linear commit does not | exchange fairness A–D | `exchange_fairness/` |
| 7 the bound site's phase | a bound Dirac electron's lower spinor is slaved; the phase mechanism has nothing to select in a passive detector; a driven absorber shows coherent control, which quantum mechanics predicts and the deposit rule as written does not | SPAD study Parts D, E, F; L13, L14 | `spad_surface_contact/` |
| 8 the record | downstream of the vertex the record copies the selection and adds an efficiency; the SPAD toy shows the placement and establishes nothing about a device | SPAD study Parts A–C | `spad_surface_contact/` |
| 9 energy | a divisible deposit with escape registers nothing at a balanced split; the indivisible quantum keeps Born, exclusivity, whole records and closure at once, and is the postulate | Part G | `spad_surface_contact/` L21–L23 |
| 10 the pair | a local rule gives Bell's foil; the detection gate rides Garg–Mermin to the closed loophole; the conditional update at the second wing must be rate-linear | pair case | `pair_case/` |

## What the paper says the specification is

Stated in one paragraph, then clause by clause: a per-site hazard linear in intensity (1), a race
and not a maximum (2), many sites per port for the phase factor to average (3), a memoryless
clock, renewed per arrival or light-tailed if frozen (4), rates on the normalized conditional
state with the first fire ending every clock (5), and, if an exchange stage precedes commitment,
a closed pot with symmetric zero-mean exchange, an absorbing zero and a Born-neutral termination
(6). Clauses 1 and 5 are Born's rule and the projection postulate; the specification does not
derive them, it locates them. Clauses 2, 3, 4 and 6 are what the runs added: they are the ways a
mechanism that has inherited Born can still lose it.

## What it does not claim

No derivation of the square; no mechanism beyond the unravelling; no prediction beyond quantum
mechanics except the two the model makes and the record excludes (the single-atom phase departure
of clause 3, and the driven-detector phase dependence of clause 7 in its registered form).

## Owed before drafting

- [x] `vacuum_race/` addendum 3, frozen light-tailed multipliers (clause 4).
- [x] `exchange_fairness/` (clause 6): run; two registered misses kept (termination at the cap; the two-site parcel model is fair because the reflecting zero is a three-body effect).
- [x] `PAPER_DRAFT.md` v0.2.

## Owed before circulation

- The references marked [verify] checked against the sources.
- JB's reading of the draft, in particular §10's naming of the worked example and §11's list of what remains.
- Whether the two model predictions that experiment excludes (§5, §9) should be stated as the model's or dropped with the phase weight.
