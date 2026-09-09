# Results: two SPAD pixels, surface contact by least rotation, the record as an induced current

**Date:** 2026-09-08. Definitions and predictions: `PREREGISTRATION.md` (unchanged by this file).
Script: `spad_contact.py` (seed 20260908, 3 s). Side ledger of departures from quantum mechanics
or experiment: `LEDGER.md`, continuing the second study's numbering.

## Read this before quoting any number

- The model is JB's, built as he described it on 2026-09-07 and 2026-09-08 (pre-registration §0).
  Nothing here is a measurement. The pixel model's parameters (lognormal disorder, 60 rows, loss
  0.1, ionization 0.5, threshold 100) are my defaults with no physical scale except the gap ratio
  θ = 0.45 for silicon at 500 nm, which is the record's number.
- **The pixel statistics and the single record per quantum are supplied by the trajectory reading
  (A1, A7), not by the phase mechanism or by the induced current.** The particle is at one place;
  the second pixel is silent by A7. The same least-rotation rule in the field reading gives the
  atom count or the brighter pixel, never Born (Table 4). This is the one thing to take from the
  study, and it was written before the run (pre-registration §0, C1, C5).
- One registered verification failed: at n = 101 the mean rotation is 1.89° against π/(n + 1) =
  1.76° with a standard error of 0.03°. Cause: D4 excludes atoms outside the pixel, so arrivals
  within 0.5 of a pixel edge see fewer than n atoms and use more rotation. Arrivals more than 0.5
  from an edge give 1.80° (a diagnostic added after the run and labelled as such in
  `verification.txt`). `run_info.txt` records the FAIL; no definition was changed.
- One prediction missed: the diffusive charge's offset from the contact has a standard deviation
  of 0.76 cells against about 0.4 predicted. Contacts ten or more columns from a lattice edge give
  0.05, contacts within ten columns give 2.36: a reflection at the lattice edge that the scratch
  model's contact distribution did not sample. Recorded under "Predictions against outcome".
- **Two addenda of the same day at the end of this file:** the empty pixel (afternoon) and the
  crystal as a dephasing transition coherence (evening, Part D). The afternoon one: the empty pixel given the same
  contact and current (the third step reopens; my prediction P7 was wrong), the SPAD's field as a
  lateral penalty on the channel, and the bulk's phase distribution.
- 4000 quanta per intensity ratio. A share's standard error is about 0.008; a miss rate's about
  0.006. At r = 1 the arrivals fell 0.4355 / 0.4145 against 0.4254 each (1.3 standard errors), so
  the r = 1 shares of 0.50 to 0.52 are that fluctuation, not a bias.

## The answer in one paragraph

With the packet carried by its own Dirac current to two SPAD pixels, one record forms per
quantum and pixel 1 takes 0.803 ± 0.008 of the records at r = 4 (Born 0.800; the density's own
share within the two apertures 0.803), and 0.50 at r = 1, under least rotation and under R2, at
every n from 1 to 1001, under both deposits. That statistic is Bohm's equivariance restricted to
the apertures. What the mechanism decides is the rest: the rotation used is π/(n + 1) under least
rotation (8.3° at n = 21, 0.19° at n = 1001) against 54° under R2; under the weighted deposit the
least-rotation atom fails to make a pair at the rate 0.468ⁿ, negligible beyond about ten atoms in
reach, while the sampled atom fails at 0.155 whatever n is, a phase-limited efficiency ceiling of
0.85. The efficiency of the primary configuration is 0.656, the fill factor 0.847 times the
breakdown probability 0.775. The induced current carries all the charge down one channel under
the greedy least-resistance rule (electrode width 0), which wanders 5.2 cells from the contact on
a fixed landscape and lands in 10 of the 121 columns of the device; the dielectric-breakdown split
and diffusion spread the charge over 6.1 cells and keep its centroid at the contact. The greedy
path crosses a third of the resistance the diffusive charge does (0.58 against 1.65 per row).

## Table 1. Arrivals by region (D2, D3; `verification.txt`)

| region | ∫ρ, r = 4 | arrivals, r = 4 | ∫ρ, r = 1 | arrivals, r = 1 |
|---|---|---|---|---|
| pixel 1, [1, 7] | 0.6865 | 0.6787 | 0.4254 | 0.4355 |
| pixel 2, [−7, −1] | 0.1680 | 0.1678 | 0.4254 | 0.4145 |
| gap, (−1, 1) | 0.0790 | 0.0860 | 0.0828 | 0.0777 |
| outside, \|x\| > 7 | 0.0665 | 0.0675 | 0.0665 | 0.0722 |

All within 3 standard errors (0.004 to 0.008). Share of ρ in pixel 1 among the two pixels: 0.8033
at r = 4 (Born 0.800; the lobes spread into the gap and outside), 0.5000 at r = 1.

## Table 2. Records with the whole quantum deposited (D4a), r = 4 (`summary.csv`)

| n in reach | least rotation: records | efficiency | share pixel 1 | rotation used | R2: records | efficiency | share pixel 1 | rotation used |
|---|---|---|---|---|---|---|---|---|
| 1 | 2651 | 0.663 | 0.806 | 90.1° | 2688 | 0.672 | 0.798 | 91.1° |
| 3 | 2651 | 0.663 | 0.804 | 46.1° | 2638 | 0.659 | 0.796 | 62.5° |
| 11 | 2638 | 0.659 | 0.801 | 15.7° | 2632 | 0.658 | 0.801 | 55.6° |
| 21 (primary) | 2625 | 0.656 | 0.803 | 8.3° | 2576 | 0.644 | 0.800 | 54.4° |
| 101 | 2660 | 0.665 | 0.800 | 1.9° | 2596 | 0.649 | 0.796 | 53.8° |
| 1001 | 2582 | 0.645 | 0.801 | 0.19° | 2668 | 0.667 | 0.805 | 54.6° |

Standard error of a share 0.008. Breakdown probability per configuration 0.76 to 0.79 (closed form
0.778). Doubles 0 in every configuration (A7). At r = 1, n = 21: share 0.505 ± 0.010 (least
rotation), 0.518 ± 0.010 (R2), against the arrivals' own 0.512; efficiency 0.659 and 0.658.

## Table 3. The weighted deposit (D4b): no pair when w < θ = 0.45, r = 4

| n in reach | least rotation: miss | predicted 0.468ⁿ | efficiency | R2: miss | predicted | efficiency |
|---|---|---|---|---|---|---|
| 1 | 0.471 | 0.468 | 0.350 | 0.487 | 0.470 | 0.347 |
| 3 | 0.115 | 0.103 | 0.588 | 0.232 | 0.234 | 0.511 |
| 11 | 0.0003 | 0.0002 | 0.659 | 0.170 | 0.167 | 0.547 |
| 21 | 0 | 10⁻⁷ | 0.656 | 0.156 | 0.158 | 0.547 |
| 101 | 0 | 0 | 0.665 | 0.156 | 0.153 | 0.548 |
| 1001 | 0 | 0 | 0.645 | 0.155 | 0.152 | 0.565 |

The shares in pixel 1 under D4b are 0.79 to 0.81 in every row: the miss is blind to the pixel.

## Table 4. The field reading of the same rules (C1; `field_reading.csv`)

Every surface atom of both pixels in play, the density under each atom as its available energy,
4000 trials.

| atoms per pixel | pure least rotation (intensity ignored) | least rotation weighted by deposit (R3) | R2 |
|---|---|---|---|
| 121 | 0.518 ± 0.008 | 1.000 | 0.804 ± 0.006 |
| 5 | 0.502 ± 0.008 | 0.995 ± 0.001 | 0.779 ± 0.007 |

Born 0.800; the atom-count share 0.500. In the field reading the least-action rule returns either
the count of atoms or the brighter pixel; Born appears only under R2, from the |ψ|² weight in the
deposit (second study, Table 2), or under the trajectory reading, from equivariance (Table 2 here).

## Table 5. The induced current (D7; `current.csv`), primary configuration, 3386 pixel arrivals

| transfer rule | landscape | electrode width (cells) | offset from contact, sd | corr(contact, centroid) | resistance crossed per row | distinct end columns, pixel 1 |
|---|---|---|---|---|---|---|
| greedy (least resistance) | fixed | 0.00 | 5.15 | 0.987 | 0.579 | 10 of 121 |
| greedy | fresh | 0.00 | 6.27 | 0.981 | 0.561 | 121 |
| DBM split ∝ 1/ρ | fixed | 6.11 | 1.04 | 1.000 | 0.920 | |
| DBM | fresh | 6.07 | 1.12 | 0.999 | 0.912 | |
| diffusion | fixed | 6.14 | 0.76 | 1.000 | 1.650 | |
| diffusion | fresh | 6.14 | 0.76 | 1.000 | 1.649 | |

Predicted: greedy width 0, offset sd 5.3 (fixed) and 6.1 (fresh), about 12 end columns; DBM 6.1
and 1.0; diffusion 6.32 and 0.4; resistance per row 0.56, between, 1.65. The greedy channel is
the lightning reading: one path, a third of the resistance, and on one device the same ten
channels for every quantum (`fig3_current.png`, middle panel). The charge conservation of the
DBM and diffusive transfers is 4 × 10⁻¹⁵.

## Verification (`results/verification.txt`)

| Check | Result |
|---|---|
| Arrivals in four regions against ∫ρ, both r | all within 3 s.e., pass |
| Galton–Watson extinction against 0.2222 | 0.2255 (s.e. 0.0066), pass |
| Diffusion width against √(2H/3) = 6.325 within 5 % | 6.140 (2.9 % low, the lattice edges), pass |
| Charge conservation, DBM and diffusion, both landscapes | ≤ 4 × 10⁻¹⁵, pass |
| Least-rotation mean against π/(n + 1), n = 1, 3, 11, 21, 1001 (r = 4) and 21 (r = 1) | pass |
| Least-rotation mean at n = 101 | 1.89° against 1.76°, s.e. 0.03°: **FAIL** (pixel-edge exclusion; interior arrivals 1.80°) |
| D4b miss against 0.468ⁿ at n = 1, 3 | 0.471, 0.115 against 0.468, 0.103, pass |

## The cases named in §3, resolved

- **C1, the field reading.** Computed (Table 4): pure least rotation 0.518 and 0.502, R3 1.000
  and 0.995, R2 0.804 and 0.779. The Born share of this study comes from the trajectory reading.
- **C2, the efficiency ceiling under R2 with D4b.** 0.845 to 0.85 at every n ≥ 11 (Table 3).
  Ledger L8.
- **C3, the undeposited remainder under D4b.** Not modelled; ledger L8.
- **C4, the record's location.** The greedy channel misplaces the contact by 5.2 cells (0.26
  Compton lengths) and lands in ten device-fixed columns; the correlation with the contact stays
  0.987 because the contacts spread over the whole pixel. Not testable in a SPAD; recorded.
- **C5, the third step closes by A7 alone.** Confirmed in the sense that nothing else in the run
  acts on the second pixel: doubles are 0 because the rule says so. Ledger L7.

## Predictions against outcome

- P1 held: arrivals within 3 s.e.; share 0.803 ± 0.008 against 0.8033; doubles 0; the share the
  same under both rules, both deposits, every n.
- P2 held: least rotation π/(n + 1) at every n, with the edge excess of 0.1° to 1° at finite n
  (the n = 101 verification failed on it, above); R2 91° at n = 1 tending to 54°.
- P3 held: least rotation 0.471, 0.115, 0.0003, 0, 0, 0; R2 0.487, 0.232, 0.170, 0.156, 0.156,
  0.155 against 0.470, 0.234, 0.167, 0.158, 0.153, 0.152.
- P4 held: breakdown 0.775 (primary) against 0.778; efficiency 0.656 against 0.665 (the arrivals
  fell 0.008 low in the pixels), D4b with R2 0.547 against 0.560.
- P5 held for greedy (0, 5.15, 10 columns) and DBM (6.11, 1.04); diffusion width 6.14 against
  6.32 (within the 5 % allowance); **diffusion offset sd 0.76 against 0.4 missed**, the lattice-edge
  reflection (0.05 in the interior, 2.36 near an edge). Times 0.579, 0.920, 1.650 against 0.56,
  between, 1.65.
- P6 held: correlations 0.981 to 1.000 against ≥ 0.96.

## The three steps, in this study

| Step | State here |
|---|---|
| 1, free flight | computed: the packet's own Dirac current, arrivals equivariant |
| 2, contact | computed: least rotation picks the closest of n, π/(n + 1); the pair forms at the gap |
| 3, the record | computed for the pixel that has the particle: avalanche, channel, location; single per quantum by A7, which is the reading, not a dynamics |

**What this study did not do (as of the morning run).** It did not test whether the second pixel
stays silent when its wave is allowed to deposit. That test was run the same afternoon at JB's
instruction and is the addendum below: the second pixel does not stay silent. JB's listed
extensions, unbound surface electrons and a depth of penetration with scatter, sit on top of that.

---

# Addendum, 2026-09-08 (afternoon): the empty pixel with the same contact and current; the bulk's phase distribution

Definitions and predictions: the dated addendum to `PREREGISTRATION.md` (D10 to D15, P7 to P12,
C6 to C8). Script: `empty_pixel.py` (seed 20260909, 5 s). Outputs: `results/empty_pixel.csv`,
`results/bulk_phase.csv`, `results/verification_addendum.txt`, `fig4_empty_pixel.png`,
`fig5_bulk_phase.png`. New landscapes were drawn for this run, so the device asymmetry differs
from the main run's.

## Read this before quoting any number

- **My registered prediction P7 was wrong as written, and the verification built on it failed.**
  I wrote that at r = 4 the empty pixel's share, 0.168, is below every gap ratio, so the empty
  pixel never pairs. That holds only when the particle lands in the bright pixel. When it lands
  in the dim one (a fifth of quanta) the *bright* pixel is the empty one, with 0.687 of the
  quantum, above the gap at every wavelength run. So under the share reading E1 the empty pixel
  pairs in 0.19 of quanta at r = 4 and 0.33 at r = 2, and the check "E1/H0 at r = 4 reproduces
  the main run" failed: share 0.880 against 0.803. The failure is the result.
- Two further misses, both recorded below: the greedy channel's doubles at short drain times
  came out lower than predicted from the main run's landscapes (the tie rate depends on the
  particular device), and the B0 rotation check failed at 3.4 standard errors on the same
  pixel-edge exclusion as the main run's n = 101 check. The pre-registration's sentence "the
  greedy channel gives the most doubles" also contradicted its own numbers, which had DBM
  highest; the run has DBM highest.
- The drain (D12) is a rule I wrote to give JB's lightning picture a definite form: the first
  channel to connect withdraws the other pixel's deposit unless that pixel connects within τ_d.
  At τ_d = 0 it is an instantaneous stop across the pixels, which is the projection postulate as
  dynamics (second study, Part G); every number at τ_d = 0 should be read with that in mind.
- Everything else in the main run's "read this first" still applies.

## The answer in one paragraph

Given the same contact and the same induced current, the pixel without the particle forms a pair
whenever its wave carries more than the gap fraction of the quantum: the bright pixel does so at
every wavelength when the particle went to the dim one, and both pixels do so at a balanced
split below about 470 nm. Once it pairs, the empty pixel races. With no head start (H0) the
record goes to whichever device channel connects first, so the share of records is pushed above
Born by the bright pixel stealing them (0.880 at r = 4, 0.811 at r = 2 against 0.800 and 0.667)
or, when both pixels pair, to the device's own asymmetry (0.48 at every r): Born is lost, and the
doubles rise with the drain time as the tie rate of the two channels (0.33 per quantum at τ_d = 5
for greedy, 0.52 with no drain). With the channel speed proportional to the deposit (H1) the
particle's pixel connects first by a margin of 1.5 to 3.5 channel times, doubles are 0 up to a
drain time of 10, and the share returns to 0.81 at r = 4 and 0.69 at r = 2; what remains is a
wrong-pixel rate of 0.03 to 0.15, the cases where the particle's own avalanche failed and the
empty pixel's succeeded. The SPAD's field, as a lateral penalty on the channel, straightens the
path (offset sd 4.4 against 6.3 cells) and widens the spread of connection times, so it lowers
the tie rate (0.23 against 0.33 at τ_d = 5) but silences nothing. On the phases: the bulk was
flat, uniform on the circle and independent, the most random an angle can be; an ordered or
locked crystal makes the rotation used and the pair rate depend on the particle's own offset,
a periodic crystal behaves like a random bulk with a hard bound on the rotation.

## Table 6. The empty pixel, greedy channel (D10 to D14; `empty_pixel.csv`)

Share = records at pixel 1 over all records; wrong = records at the pixel the particle did not
reach, per quantum; doubles per quantum; "both" = quanta with a candidate at each pixel.

| r | θ (nm) | empty pixel pairs | reading | head start | τ_d | efficiency | share pixel 1 | Born | wrong pixel | doubles | both |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 4 | any of 0.45 to 0.28 | 0.19 (the bright pixel, when the particle went dim) | E1 share | H0 none | 0 | 0.698 | 0.880 | 0.800 | 0.077 | 0 | 0.105 |
| 4 | | | E1 | H0 | 5 | 0.698 | 0.847 | 0.800 | 0.044 | 0.060 | 0.105 |
| 4 | | | E1 | H1 ∝ deposit | 0 | 0.698 | 0.811 | 0.800 | 0.029 | 0 | 0.105 |
| 4 | | | E1 | H1 | 5 | 0.698 | 0.812 | 0.800 | 0.029 | 0.003 | 0.105 |
| 4 | | 1.00 | E2 whole | H0 = H1 | 0 | 0.817 | 0.486 | 0.800 | 0.413 | 0 | 0.520 |
| 4 | | | E2 | | 5 | 0.817 | 0.483 | 0.800 | 0.252 | 0.328 | 0.520 |
| 2 | 0.45, 0.407 (500, 450) | 0.33 (the bright pixel) | E1 | H0 | 0 | 0.702 | 0.811 | 0.667 | 0.133 | 0 | 0.167 |
| 2 | | | E1 | H0 | 5 | 0.702 | 0.765 | 0.667 | 0.077 | 0.103 | 0.167 |
| 2 | | | E1 | H1 | 0, 5 | 0.702 | 0.694 | 0.667 | 0.051 | 0 | 0.167 |
| 2 | 0.28 (310) | 1.00 (either pixel) | E1 | H0 | 0 | 0.804 | 0.484 | 0.667 | 0.416 | 0 | 0.507 |
| 2 | 0.28 | | E1 | H1 | 0, 5, 10 | 0.804 | 0.606 | 0.667 | 0.153 | 0 | 0.507 |
| 2 | 0.28 | | E1 | H1 | 20 | 0.804 | 0.625 | 0.667 | | 0.041 | 0.507 |
| 2 | any | 1.00 | E2 | | 0 | 0.804 | 0.484 | 0.667 | 0.416 | 0 | 0.507 |
| 1 | 0.45 (500) | 0 (0.425 < 0.45) | E1 | any | any | 0.657 | 0.496 | 0.500 | 0 | 0 | 0 |
| 1 | 0.407 and below (450 and below) | 1.00 | E1 | H0 | 5 | 0.799 | 0.478 | 0.500 | 0.244 | 0.317 | 0.509 |
| 1 | 0.407 and below | | E1 | H1 | 0 to 10 | 0.799 | 0.496 | 0.500 | 0.142 | 0 | 0.509 |
| 1 | 0.407 and below | | E1 | H1 | 20 | 0.799 | 0.492 | 0.500 | | 0.006 | 0.509 |
| 1 | any | 1.00 | E2 | | 5 | 0.799 | 0.478 | 0.500 | 0.244 | 0.317 | 0.509 |

Where both pixels pair (E2, or E1 at r = 2 with θ = 0.28, or r = 1 below 470 nm) the share at
τ_d = 0 is the device's channel asymmetry and changes with the transfer rule: at r = 2,
θ = 0.28, E1/H0 it is 0.484 (greedy), 0.437 (DBM), 0.605 (diffusion), 0.415 (field-biased
greedy), for a Born value of 0.667.

## Table 7. Doubles against the drain time, E2 with no head start, r = 4 (per quantum; 0.520 of quanta have candidates at both pixels)

| τ_d | 0 | 1 | 2 | 5 | 10 | 20 | ∞ |
|---|---|---|---|---|---|---|---|
| greedy | 0 | 0.086 | 0.151 | 0.328 | 0.476 | 0.517 | 0.519 |
| DBM split ∝ 1/ρ | 0 | 0.136 | 0.277 | 0.470 | 0.519 | 0.519 | 0.519 |
| diffusion | 0 | 0.065 | 0.124 | 0.280 | 0.441 | 0.519 | 0.519 |
| greedy with the field, b = 3 (D12b) | 0 | 0.050 | 0.095 | 0.229 | 0.401 | 0.513 | 0.519 |

Predicted as fractions of the both-candidate rate: greedy 0.21, 0.40, 0.81, 0.98, 1.00 (from the
main run's landscapes; got 0.17, 0.29, 0.63, 0.92, 1.00, a miss at τ_d ≤ 5); DBM 0.25, 0.47,
0.89, 1.00 (got 0.26, 0.53, 0.91, 1.00); diffusion 0.11, 0.22, 0.52, 0.86 (got 0.13, 0.24,
0.54, 0.85); field 0.09, 0.19, 0.43, 0.76, 0.99 (got 0.10, 0.18, 0.44, 0.77, 0.99).

## Table 8. The field (D12b): the particle's contacts at r = 4

| rule | offset from contact, sd (cells) | resistance per row | sd of the connection time | corr(contact, end) |
|---|---|---|---|---|
| greedy, b = 0 | 6.31 | 0.574 | 4.44 | 0.981 |
| greedy with the field, b = 3 | 4.37 | 0.865 | 6.23 | 0.991 |

Predicted (P12): the offset sd falls by about 0.6 (got 0.69), the resistance per row rises to
about 0.88 (got 0.865), the spread of T rises by about 1.6 (got 1.40). The field straightens the
channel and makes it costlier and less regular in time.

## Table 9. The bulk's phase distribution (D15; `bulk_phase.csv`), 21 atoms in reach, least rotation, D4b with θ = 0.45, r = 4

| bulk model | δ_p = 0° | 45° | 90° | 135° | 180° |
|---|---|---|---|---|---|
| B0 random, uniform iid (as run): rotation | 8.65° | 8.57° | 8.38° | 8.52° | 8.48° |
| B0: no pair | 0 | 0 | 0 | 0 | 0 |
| B1 ordered, von Mises κ = 1: rotation | 4.07° | 5.27° | 10.80° | 19.82° | 25.85° |
| B1 κ = 1: no pair | 0 | 0 | 0 | 0.0003 | 0.0044 |
| B1 ordered, κ = 4: rotation | 1.81° | 5.57° | 31.96° | 74.98° | 108.97° |
| B1 κ = 4: no pair | 0 | 0 | 0 | 0.114 | 0.791 |
| B2 locked crystal, one phase: rotation | 0° | 45° | 90° | 135° | 180° |
| B2: no pair | 0 | 0 | 0 | 1 | 1 |
| B3 periodic crystal, q = 137.5°: rotation | 5.21° | 4.98° | 5.25° | 4.87° | 4.87° |
| B3 q = 137.5°: no pair | 0 | 0 | 0 | 0 | 0 |
| B3 periodic crystal, q = 120° (δ_p = 0, 30, 60°): rotation | 0° | 30° | 60° | | |
| B3 q = 120°: no pair | 0 | 0 | 0 | | |

Predicted (P11) and held, with the same edge excess on the random and golden models (8.65°
against 8.18°, 5.2° against 4.5°): an ordered or locked bulk turns the pair rate into a function
of the particle's own offset, ½[1 + cos δ_p]-like at κ = 4 and a step at 95.7° for the locked
crystal; the periodic crystal is a random bulk with a hard bound on the rotation. Under D4a
(the winner takes the whole quantum) none of this changes any statistic.

## Verification (`results/verification_addendum.txt`)

| Check | Result |
|---|---|
| Empty-pixel contact points against ρ restricted to the pixel, both pixels, three r | max CDF difference 0.009 to 0.022, pass |
| E1/H0 at r = 4, θ = 0.45 reproduces the main run's share, doubles 0 | share 0.880 against 0.803: **FAIL** (the prediction was wrong; the bright empty pixel pairs) |
| Doubles at τ_d = 0 | 0 in every configuration, pass |
| Doubles at τ_d = ∞ under E2 equal the both-candidate rate | 0.5195 = 0.5195 (predicted 0.515 ± 0.008), pass |
| B0 rotation against π/22 | 8.65° against 8.18°, s.e. 0.14°: **FAIL** (pixel-edge exclusion, as the main run's n = 101) |
| B2 rotation equals \|δ_p\| | to 10⁻¹⁶, pass |
| B3 q = 120° at δ_p = 30° | 30.00000000°, pass |

## Predictions against outcome

- **P7 wrong as written.** The empty pixel pairs at r = 4 and r = 2 whenever the particle went
  to the dim pixel; only its second half (r = 1 pairs below 470 nm; the race is a tie rate; the
  device sets the winner) held. The share under E1/H0 is above Born, not Born.
- **P8 held:** E2 gives the device's asymmetry, 0.48 here at every r, wrong-pixel 0.41, doubles
  0.52 with no drain.
- **P9 held in form:** H1 gives doubles 0 up to τ_d = 10 and the particle's pixel first; but the
  share is not Born, it is Born plus the wrong-pixel leakage (0.811 at r = 4, 0.694 at r = 2,
  0.606 at r = 2 with both pixels pairing), and the leakage is 0.03 to 0.15, not the 0.15 I
  wrote for every case.
- **Table 7:** the greedy tie rate missed at τ_d ≤ 5 (device-dependent); DBM, diffusion and the
  field-biased greedy held. The pre-registration's sentence naming greedy the worst was
  inconsistent with its own numbers; DBM is the worst.
- **P11 held** in every number to within the edge excess. **P12 held** in every ratio.

## The cases named in the addendum, resolved

- **C6 held:** the lightning channel does not silence the second pixel. Under H0 the record's
  pixel is a property of the device, and the doubles are the tie rate the drain time allows.
- **C7 held and widened:** exclusivity under E1 rests on the gap, and the gap fails not only at
  a balanced split below 470 nm but at any split, at any wavelength, whenever the particle goes
  to the dim pixel and the bright one is empty. Only a head start (H1) or an instantaneous drain
  (τ_d = 0, the nonlocal stop) keeps one record per quantum.
- **C8:** energy is not conserved under E1 with D4a. Ledger.

## What the empty-pixel test says about the third step

The main run closed the third step by assuming the empty pixel deposits nothing. Lifting that
assumption, with JB's own contact and current rules at the empty pixel, the step reopens exactly
as in the second study's Part G and the third study: a second pixel with enough of the wave forms
its own record unless something acts across the pixels faster than the channels connect, and the
lightning channel is not that something. What the induced current adds that the earlier sinks did
not is a *timing* variable, the connection time, so that a head start proportional to the deposit
(H1) is enough to keep doubles at zero for drain times up to a third of the channel time. Whether
a SPAD's avalanche is faster for a larger deposit is a physical question the model does not
answer; in silicon the first carrier's energy above the gap is thermalised in about a
picosecond before the avalanche, which argues against it.

**JB's question of 2026-09-08, after the run: does the interference cancel the deposited energy?**
At a mismatched atom it does: the deposit weight is the interference term w = cos²(Δφ/2), zero
for opposite second spinors (A1 of the second study). But least rotation, JB's own rule, picks the
atom in reach where the interference is constructive, and with n atoms in reach one is within
0.45 rad (w ≥ 0.95) with 95 % confidence once n ≥ 20 (second study, Table 1). So the cancellation
protects the empty pixel only when very few atoms are in reach. Closed forms for the empty pixel's
pair rate with the weighted deposit s·w ≥ θ (θ = 0.45), least rotation, n atoms: the bright pixel
empty at r = 4 (s = 0.687) needs |Δφ| ≤ 71.9°, pair rate 1 − 0.601ⁿ = 0.40, 0.78, 0.996, 1.000 at
n = 1, 3, 11, 21; at r = 2 (s = 0.570) needs |Δφ| ≤ 54.6°, 1 − 0.697ⁿ = 0.30, 0.66, 0.98, 1.00; at
r = 1 (s = 0.425) never, since θ/s > 1 at 500 nm; under E2 (s = 1) 1 − 0.468ⁿ = 0.53, 0.90,
1.00, 1.00. A silicon surface holds of order 10⁶ atoms within a wavelength of the arrival, so in a
real SPAD the interference cancels nothing at the selected atom. It would cancel across a whole
pixel only in a locked crystal (B2) whose phase is opposite to the particle's offset, and then the
pixel that records is chosen by the particle's phase and not by the intensity, which is the phase
readout of Table 9, not Born.

---

# Addendum 2, 2026-09-08 (evening): Part D, the crystal as a transition coherence with thermal dephasing

Definitions and predictions: addendum 2 of `PREREGISTRATION.md` (D16 to D18, P13, C9). Script:
`transition_coherence.py` (seed 20260910, 3 s). Outputs: `results/dephasing_single.csv`,
`results/dephasing_pixels.csv`, `results/verification_dephasing.txt`, `fig6_dephasing.png`.

## Read this before quoting any number

- Every verification passed and every registered number held (P13). No definition changed.
- **The cost case C9 comes first.** Part D describes a detector whose transition coherence was
  driven within about one T2 of the particle's arrival. In silicon at room temperature T2 is of
  order 10 to 100 fs, and at equilibrium the coherence is absent, suppressed by exp(−ħω/kT) ≈
  e⁻⁹⁷ at 2.5 eV. A passive SPAD is the random bulk of every earlier run. What Part D shows is
  what the model would show in a *driven* detector, and where its phase mechanism would be
  looked for.
- The MRI mapping is exact in form: the preparation is the 90° pulse, μ_k the pulse phase at
  pixel k, T2 the irreversible random walk, T2* the static site offsets, and R(t) the transverse
  magnetization. No echo was run; the T2* part is refocusable in principle and the T2 part is not.
- One pixel uses pure phase draws with exactly 21 atoms, so the pixel-edge excess of Parts A to
  C is absent there (the rotation at t = 10 is 8.18° against π/22 = 8.18°). Two pixels use the
  geometry.

## The answer in one paragraph

Given a shared transition coherence at t = 0, the bulk is the locked crystal of Part C and
dephases into the random bulk over about one T2, with the mean resultant R(t) = exp(−s²/2)
reproduced to 0.002. Least rotation tracks the dephasing: for a particle in phase with the
preparation the rotation used rises from 0 to 8.2°; for one in antiphase it falls from 180°
through 124° at t = 0.1 T2 and 83° at 0.3 T2 to 8.2° by 3 T2; the pair rate under the weighted
deposit follows, 0 until t ≈ 0.1 T2 and 1 from t ≈ T2. Two pixels pulsed in antiphase take every
record of an in-phase particle at t ≤ 0.03 T2 (share 1.000 against Born 0.800), 0.994 at 0.1 T2,
0.846 at 0.3 T2 and Born from t ≈ T2; for particles of random offset the two pixels' pair rates
are equal by symmetry and the share is Born at every delay, while the efficiency rises from 0.36
at t = 0 to 0.67 by t ≈ T2. So in this model the Born share and a phase-independent efficiency
are what a dephased detector shows, and the phase mechanism leaves a trace only inside one T2 of
a preparation, in the pixel share for a source of definite offset and in the efficiency for any
source.

## Table 10. One pixel, homogeneous dephasing (D17; `dephasing_single.csv`), 21 atoms, least rotation, θ = 0.45

| t/T2 | R(t) | rotation, δ_p = 0 | δ_p = 90° | δ_p = 180° | pair rate, δ_p = 180° | pair rate, δ_p uniform | winner's w, uniform |
|---|---|---|---|---|---|---|---|
| 0 | 1.000 | 0° | 90.0° | 180.0° | 0 | 0.525 | 0.498 |
| 0.03 | 0.971 | 0.8° | 63.4° | 149.3° | 0 | 0.676 | 0.644 |
| 0.1 | 0.905 | 1.5° | 41.6° | 123.9° | 0.021 | 0.798 | 0.750 |
| 0.3 | 0.741 | 2.6° | 16.3° | 83.4° | 0.704 | 0.964 | 0.884 |
| 1 | 0.369 | 4.7° | 8.5° | 23.2° | 0.999 | 1.000 | 0.980 |
| 3 | 0.048 | 7.6° | 8.2° | 9.1° | 1 | 1 | 0.990 |
| 10 | 0.002 | 8.2° | 8.1° | 8.3° | 1 | 1 | 0.990 |

Predicted: R 1, 0.970, 0.905, 0.741, 0.368, 0.050, 0; rotations 0, 0.8, 1.5, 2.6, 4.6, 7.5, 8.2 /
90, 63.5, 41.7, 16.3, 8.5, 8.2, 8.2 / 180, 149, 124, 83, 23, 9.1, 8.2; pair rates 0, 0, 0.021,
0.71, 0.999, 1, 1 and 0.532, 0.679, 0.802, 0.964, 1, 1, 1. All held (the uniform pair rate at
t = 0 is 0.525 against 0.532, s.e. 0.0035, within 3 s.e.). The inhomogeneous and combined cases
give the same numbers at equal s, with s = t/T2* and √(2t/T2 + t²/T2*²): at t = 0.3 the
inhomogeneous bulk still has R = 0.956 and a uniform-offset pair rate of 0.73, the combined
0.708 and 0.98.

## Table 11. Two pixels, r = 4, trajectory reading (D18; `dephasing_pixels.csv`); arrival share 0.801, Born 0.800

| t/T2 | preparation | δ_p | pair rate pixel 1 | pair rate pixel 2 | efficiency | share pixel 1 |
|---|---|---|---|---|---|---|
| 0 | opposite | fixed 0 | 1.000 | 0.000 | 0.534 | **1.000** |
| 0.03 | opposite | fixed 0 | 1.000 | 0.000 | 0.534 | **1.000** |
| 0.1 | opposite | fixed 0 | 1.000 | 0.024 | 0.538 | 0.994 ± 0.002 |
| 0.3 | opposite | fixed 0 | 1.000 | 0.735 | 0.631 | 0.846 ± 0.007 |
| 1 | opposite | fixed 0 | 1.000 | 0.999 | 0.666 | 0.802 ± 0.008 |
| 3, 10 | opposite | fixed 0 | 1 | 1 | 0.667 | 0.802 ± 0.008 |
| 0 | opposite | uniform | 0.544 | 0.527 | 0.355 | 0.808 ± 0.010 |
| 0.1 | opposite | uniform | 0.797 | 0.796 | 0.531 | 0.805 ± 0.009 |
| 0.3 | opposite | uniform | 0.960 | 0.966 | 0.641 | 0.800 ± 0.008 |
| ≥ 1 | opposite | uniform | 1 | 1 | 0.667 | 0.802 ± 0.008 |
| 0 | same | uniform | 0.529 | 0.539 | 0.355 | 0.800 ± 0.011 |
| 0.1 | same | uniform | 0.796 | 0.810 | 0.535 | 0.798 ± 0.009 |
| any | same | fixed 0 | 1 | 1 | 0.667 | 0.802 ± 0.008 |

Inhomogeneous (T2* = T2) at t = 0.3, opposite, fixed 0: pair rate pixel 2 0.000, share 1.000; at
t = 1: 0.956 and 0.809. Combined at t = 0.3: 0.810 and 0.835. Predicted (hom): 1.000, 1.000,
0.995, 0.851, 0.803, 0.803, 0.803 for the fixed offset; 0.803 throughout for the uniform offset
and the same preparation; efficiencies 0.35, 0.45, 0.53, 0.64, 0.66 (uniform) and 0.53, 0.53,
0.53, 0.62, 0.66 (opposite, fixed): all held.

## Verification (`results/verification_dephasing.txt`)

| Check | Result |
|---|---|
| R measured against exp(−s²/2), seven delays | worst difference 0.0016, pass |
| Rotation at t = 0 equals \|δ_p\|, five offsets | to 10⁻¹⁶, pass |
| Rotation at t = 10 against π/22 | 8.18° against 8.18°, pass |
| Uniform-offset pair rate at t = 0 against 0.5325 | 0.5250 (s.e. 0.0035), pass |
| Two-pixel share at t = 10 against the arrival share | 0.8016 against 0.8014, pass |

## What Part D says

- **On JB's question whether a crystal reduces the phase spread between atoms:** a prepared
  coherence gives zero spread at t = 0, and thermal dephasing returns it to the random bulk in
  about one T2. The lattice by itself does not prepare the coherence; a pulse does.
- **Where the mechanism would show:** in a detector driven within T2 of the arrival, as a pixel
  share that follows the particle's offset relative to the pulse phase (1.000 against 0.800 for
  antiphase pixels) and as an efficiency 0.36 against 0.67 for random offsets. Both are gone by
  t ≈ T2, 10 to 100 fs in silicon.
- **On Born:** in this model the intensity share is what every dephased detector shows, under the
  trajectory reading, and the deviation is a coherence effect with the lifetime of the
  coherence. That is a statement the model makes and quantum mechanics does not; ledger L13.

---

# Addendum 3, 2026-09-08 (night): Part E, the bound electron's second spinor, and the deposit weight with unequal amplitudes

Definitions and predictions: addendum 3 of `PREREGISTRATION.md` (D19 to D22, P17, C10). Script:
`bound_spinor.py` (17 s; dense diagonalisation of the 1+1-D Dirac Hamiltonian with a well).
Outputs: `results/bound_spinor.csv`, `results/bound_visibility.csv`,
`results/verification_bound.txt`, `fig8_bound_spinor.png`.

## Read this before quoting any number

- Every verification passed and every registered number held. No definition changed.
- This part answers two statements of JB's, both correct: the second Dirac spinor of a bound
  electron is a complete complex wave with its own amplitude, so |A| ≠ |B| in the interference
  weight; and the ground state of a bound electron is not the ground state of a free electron.
  What the computation adds is *what* that second spinor is: its amplitude, its phase, and whether
  either is free.
- The wells are far deeper than a silicon valence electron's (bindings 0.03 to 0.35 of the rest
  energy against 10⁻⁵); the silicon numbers are extrapolated by the scaling of the two shallowest
  wells and are order-of-magnitude estimates.
- Both readings of "the second spinor" are measured on the same states: JB's, the lower row of
  the Dirac spinor; and the second study's A4, the negative-energy free component.

## The answer in one paragraph

A Dirac electron bound in a well has a second spinor that is slaved to the first: in this
representation χ = −iφ′/(E_b + m − V) exactly, so its amplitude is set by the gradient of the
first spinor and its phase sits at +90° on one side of the well and −90° on the other, fixed to
10⁻⁸ degrees, the same for every identical atom. Its amplitude relative to the first spinor is
0.06 for a binding of 0.03 rest energies and 0.30 for 0.35, scaling as the square root of the
binding, about 10⁻³ for a silicon valence electron. In the free-particle reading the bound state's
negative-energy content is far smaller, 1.5 × 10⁻⁶ of the norm at the shallow end, scaling as the
2.75 power of the binding, about 10⁻¹⁶ for silicon; its offset against the positive-energy part
is not a single number (resultant over momentum 0.42 to 0.10), and both parts rotate together at
E_b, not against each other. With unequal amplitudes the model's weight becomes
[a² + b² + 2ab cos Δφ]/(a + b)², of visibility 2ab/(a² + b²): zero for an incident particle
without a second-spinor amplitude of its own, unity only when the two amplitudes are equal.

## Table 12. The bound electron (D19, D20; `bound_spinor.csv`)

| well U₀, L | binding / m | bound states | lower-row fraction f_L | amplitude a_L | phase of χ against φ (left, right) | −E fraction f₋ | amplitude a₋ | ±E offset: mean, resultant over p |
|---|---|---|---|---|---|---|---|---|
| 0.05, 4 | 0.028 | 1 | 0.0041 | 0.064 | −90.0°, +90.0° | 1.5 × 10⁻⁶ | 0.0012 | 168°, 0.42 |
| 0.2, 2 | 0.111 | 2 | 0.0167 | 0.129 | −90.0°, +90.0° | 6.8 × 10⁻⁵ | 0.0082 | −5°, 0.21 |
| 0.6, 1 | 0.302 | 2 | 0.0537 | 0.232 | −90.0°, +90.0° | 1.3 × 10⁻³ | 0.036 | −89°, 0.12 |
| 1.0, 0.5 | 0.346 | 1 | 0.0904 | 0.301 | −90.0°, +90.0° | 4.7 × 10⁻³ | 0.068 | −84°, 0.10 |

The phase spread of χ against φ across each well is below 10⁻⁸ degrees; the first spinor is real
to 10⁻¹³ (a global phase removed). Scaling through the two shallowest wells: f_L ∝ binding^1.01,
f₋ ∝ binding^2.75. Silicon valence electron (binding 10⁻⁵): f_L ≈ 1.4 × 10⁻⁶ (a_L ≈ 1.2 × 10⁻³),
f₋ ≈ 5 × 10⁻¹⁶ (a₋ ≈ 2 × 10⁻⁸). Eigenvector residuals 3 × 10⁻¹⁴; the identity for χ holds to
3 × 10⁻¹².

## Table 13. The model's phase mechanism with unequal amplitudes (D22; `bound_visibility.csv`)

Visibility V = 2ab/(a² + b²); the weight's minimum is [(a − b)/(a + b)]².

| bound electron's amplitude a | incident b = 0 | b = a | b = 0.1 | b = 0.5 | b = 1 |
|---|---|---|---|---|---|
| lower row, shallow well, 0.064 | 0 | 1 | 0.91 | 0.25 | 0.13 |
| lower row, deepest well, 0.301 | 0 | 1 | 0.61 | 0.88 | 0.55 |
| −E component, shallow well, 0.0012 | 0 | 1 | 0.024 | 0.005 | 0.002 |
| silicon, lower row, ≈ 10⁻³ | 0 | 1 | ≈ 0.02 | ≈ 0.005 | ≈ 0.002 |
| silicon, −E component, ≈ 10⁻⁸ | 0 | 1 | ≈ 10⁻⁷ | ≈ 10⁻⁸ | ≈ 10⁻⁸ |

A physical free electron has b = 0 in the negative-energy reading (a lag θ gives b = sin(θ/2) at a
pair-scale energy cost, L4) and b = p/(E + m) in the lower-row reading (0.1 at p = 0.2 m); a photon
has neither.

## What Part E says

- **JB's two statements hold, and they cut against the model's bulk.** The bound electron's
  second spinor has its own amplitude, but that amplitude and its phase are functions of the first
  spinor: quadrature, fixed, identical across identical atoms. The random, independent, counter-
  rotating second-spinor phases of the second study's D1 and D2 are not properties of a Dirac bound
  state in either reading. A bulk of bound electrons is the locked crystal B2 of Part C with the
  phase at ±90°, not the random bulk B0.
- **With the physical amplitudes the phase mechanism's visibility is small or zero.** In the
  negative-energy reading a physical incident electron has b = 0 and the mechanism is off; a lag
  admixture equal to the bound electron's (10⁻⁸ for silicon) would switch it on at a pair-scale
  energy cost. In the lower-row reading a slow free electron has b ≈ 0.1 against the bound
  electron's 10⁻³, visibility 0.02.
- **The L13 comparison is revised, not reversed.** With a ≠ b the model's cross term is 2ab cos Δφ
  on a floor of a² + b², which has the form of the quantum-mechanical cross term with A and B as
  the two second-spinor amplitudes; the differences in sign, order in the field, and the
  population factor are as in `L13_TEST.md`, and the visibility is now a further suppression.
- **On "not the same ground state as a free electron."** Correct, and in a specific way: the bound
  ground state expanded in free spinors has negative-energy content, with a definite but
  momentum-dependent phase, and that content is 10⁻⁶ of the norm at a binding of 0.03 rest
  energies and about 10⁻¹⁶ for silicon.


---

# Addendum 4, 2026-09-08 (night): Part F, the bulk rebuilt on the slaved second spinor

Definitions and predictions: addendum 4 of `PREREGISTRATION.md` (D23 to D26, P18, C11). Script:
`slaved_bulk.py` (seed 20260911, 5 s). Outputs: `results/slaved_bulk.csv`, `results/slaved_a_of_d.csv`,
`results/verification_slaved.txt`, `fig9_slaved.png`.

## Read this before quoting any number

- Every verification passed and every registered number held (P18). No definition changed.
- The bulk is now what Part E measured: every atom presents a second spinor of one amplitude a and
  a phase of +90° or −90° according to the side on which the particle touches it. Nothing in the
  bulk is random. The incident particle brings its own amplitude b and offset δ_p: 0° for a
  running wave, ±90° for a standing one, uniform as the sensitivity for a partially reflected wave.
- Part E's amplitudes are used: a = 0.064 for the shallow well (0.078 at a contact four Compton
  lengths from the nucleus, 0.119 in the far tail) and a ≈ 0.0012 for a silicon valence electron,
  an order-of-magnitude extrapolation.

## The answer in one paragraph

On the slaved bulk the phase race has nothing to select. Every atom in reach presents one of two
phases and the same amplitude; a running incident wave, offset 0°, is 90° from both, so the
interference term vanishes at every atom and the weight is the same everywhere,
(a² + b²)/(a + b)², which is never below ½. Least rotation then picks a side of an atom, not an
atom, and the record's pixel is the arrival's: share 0.804 ± 0.008 against the arrival share
0.803 in every one of the 48 configurations that record. At the silicon gap ratio 0.45 every
configuration makes a pair, so the efficiency is 0.653 throughout, the fill factor times the
breakdown probability. The weight bites only above a gap ratio of ½ with the two amplitudes
comparable: at 0.6 (660 nm) a running wave meeting a bulk of amplitude 0.064 makes no pair at all
for b = 0.025, 0.098 or 0.064, while silicon's amplitude 0.0012 against an electron's 0.025 or
0.098 gives a weight of 0.91 to 0.98 and a pair every time, and a photon (b = 0) a weight of 1
always. The randomness JB asked about is not in the atom: the rotation used is 90° for every
running wave, 0° for a standing wave at the matching position, and uniform on 0° to 90° only if
the incident wave's own offset is uniform.

## Table 14. The slaved bulk (D23 to D25; `slaved_bulk.csv`), r = 4, least rotation, 21 atoms in reach

| a | b | δ_p | rotation used | chosen weight | pair rate, θ = 0.45 | pair rate, θ = 0.6 | efficiency (0.45) | share pixel 1 |
|---|---|---|---|---|---|---|---|---|
| 0.064 | 0 (photon) | any | 90° / 0° / 44.6° | 1.000 | 1 | 1 | 0.653 | 0.804 |
| 0.064 | 0.025 | 0° | 90° | 0.596 | 1 | **0** | 0.653 | 0.804 |
| 0.064 | 0.098 | 0° | 90° | 0.522 | 1 | **0** | 0.653 | 0.804 |
| 0.064 | 0.064 | 0° | 90° | 0.500 | 1 | **0** | 0.653 | 0.804 |
| 0.064 | 0.025, 0.098, 0.064 | 90° | 0° | 1.000 | 1 | 1 | 0.653 | 0.804 |
| 0.064 | 0.025 | uniform | 44.6° | 0.856 | 1 | 0.994 | 0.653 | 0.804 |
| 0.064 | 0.098 | uniform | 44.6° | 0.829 | 1 | 0.904 | 0.653 | 0.804 |
| 0.064 | 0.064 | uniform | 44.6° | 0.821 | 1 | 0.883 | 0.653 | 0.804 |
| 0.0012 (silicon) | 0.025 | 0° | 90° | 0.913 | 1 | 1 | 0.653 | 0.804 |
| 0.0012 | 0.098 | 0° | 90° | 0.976 | 1 | 1 | 0.653 | 0.804 |
| 0.0012 | 0.0012 | 0° | 90° | 0.500 | 1 | **0** | 0.653 | 0.804 |
| 0.0012 | 0.025, 0.098 | uniform | 44.6° | 0.97, 0.99 | 1 | 1 | 0.653 | 0.804 |

Predicted (P18): weights 0.522, 0.500, 0.596, 1, 0.976 at δ_p = 0; pair rates at θ = 0.6 of 0 for
the three a = 0.064 rows and 0.90, 0.87, 0.99 for the uniform ones; rotation 90°, 0°, mean 45°
(got 44.63° ± 0.44°); share the arrival share throughout. All held. The side chosen at δ_p = 90°
is the atoms left of the arrival in every case; at δ_p = 0° the nearest atom, left in 0.503.

## Table 15. The amplitude the atom presents at the contact (`slaved_a_of_d.csv`, Part E's shallow well)

| contact distance d (Compton lengths) | 0.5 | 1 | 2 | 4 | 6 | 8 | 12 | 16 | tail limit κ/(E_b + m) |
|---|---|---|---|---|---|---|---|---|---|
| a(d) = \|χ\|/\|φ\| | 0.011 | 0.022 | 0.043 | 0.078 | 0.100 | 0.112 | 0.118 | 0.116 | 0.119 |
| phase of χ against φ | +90° throughout d > 0 (−90° for d < 0) |

## What Part F says

- **The phase race is gone.** Once the bulk's second spinor is what a Dirac bound state has, no
  atom differs from another in phase except by the side of contact, and a running incident wave
  is in quadrature with both sides. The which-atom selection reduces to the nearest atom, the
  which-pixel statistics to the arrival's, Born by the trajectory reading as in every earlier run.
- **What is left of the mechanism is a threshold.** For gap ratios above ½ and comparable
  amplitudes a running wave deposits too little to make a pair, a detector blind to running waves
  and open to standing ones. With silicon's amplitude against any electron or photon the weight
  stays above 0.9 and the threshold never acts. Ledger L16.
- **On JB's question.** Within one bound electron the first-to-second-spinor phase difference is
  fixed at ±90° by the side, not random; in three dimensions it is −i times a direction-dependent
  spin factor. The only random phase in the rebuilt model is the incident particle's own, and
  that exists for a partially reflected wave and not for a photon.

---

# Addendum 5, 2026-09-08 (late night): Part G, release and retry at the excited electron

Definitions and predictions: addendum 5 of `PREREGISTRATION.md` (D27 to D31, P19, C12 to C14).
Script: `release_retry.py` (seed 20260914, 3000 quanta per configuration, 21 s). Outputs:
`results/release_retry.csv`, `results/verification_release.txt`, `fig10_release.png`.

## Read this before quoting any number

- Every verification passed and every registered number held (P19). No definition changed.
- The rates are my defaults in units of the re-emission rate; nothing in them carries a physical
  scale. The sensitivities show which conclusions depend on them (the efficiency) and which do
  not (the share, the doubles, the record energy).
- Two readings of the deposit are run side by side and the difference between them is the
  result. Reading A is JB's shares as energies with the three release channels. Reading B is the
  indivisible quantum, the postulate, with the same channels.

## The answer in one paragraph

With the excited electron given three ways to release its holding, ejection above the gap,
re-emission back to the wave for another try, and phonon loss, the energy of the quantum closes
to 10⁻¹⁵ in every configuration, so the ledger item L11 is resolved in both readings. In reading
A the re-emission cycle (about 15 returns per quantum) does what the drain of Part B could not:
above a gap ratio of one half a double record is impossible, because the energy left after a
record of at least θ is less than θ, and at the silicon ratio 0.45 none occurred in 36 000
quanta; doubles return only below θ = 0.36 (0.20 at a balanced split, 0.30 at r = 2 for
θ = 0.28). The price is Born and the record's energy. The share becomes a step in r: the brighter
pixel takes 1.000 at r = 4 and 0.97 at r = 2 for θ = 0.6, 0.996 and 0.93 at 0.45, against 0.80
and 0.67; and each record carries 0.45 to 0.72 of the quantum, 0.65 at r = 4 and 0.45, never the
whole, while calorimetric detectors measure the whole photon energy per absorption. In reading B
everything holds at once, Born to 0.01 at every r and θ, no doubles, the full quantum per record,
closure exact, with the channels adding only an efficiency of 0.58 to 0.60; and nothing is
derived, because the indivisibility of the deposit is the assumption that does all the work.

## Table 16. Reading A, shares as energies, primary rates (`release_retry.csv`)

| r | θ (nm) | share pixel 1 | Born | doubles | efficiency | record energy | cycles | heat | escape |
|---|---|---|---|---|---|---|---|---|---|
| 1 | 0.6 (660) | 0.52 | 0.50 | 0 | 0.32 | 0.63 | 17.8 | 0.44 | 0.37 |
| 1 | 0.45 (500) | 0.49 | 0.50 | 0 | 0.39 | 0.60 | 17.4 | 0.41 | 0.35 |
| 1 | 0.36 (400) | 0.50 | 0.50 | 0.20 | 0.80 | 0.47 | 11.8 | 0.26 | 0.28 |
| 1 | 0.28 (310) | 0.50 | 0.50 | 0.21 | 0.83 | 0.46 | 12.5 | 0.25 | 0.28 |
| 2 | 0.6 | 0.97 | 0.67 | 0 | 0.18 | 0.72 | 17.8 | 0.48 | 0.39 |
| 2 | 0.45 | 0.93 | 0.67 | 0 | 0.60 | 0.59 | 15.7 | 0.33 | 0.31 |
| 2 | 0.36 | 0.83 | 0.67 | 0 | 0.69 | 0.56 | 15.9 | 0.31 | 0.30 |
| 2 | 0.28 | 0.58 | 0.67 | 0.30 | 0.85 | 0.45 | 11.3 | 0.23 | 0.26 |
| 4 | 0.6 | 1.000 | 0.80 | 0 | 0.49 | 0.71 | 16.2 | 0.34 | 0.31 |
| 4 | 0.45 | 0.996 | 0.80 | 0 | 0.70 | 0.65 | 16.0 | 0.27 | 0.28 |
| 4 | 0.36 | 0.97 | 0.80 | 0 | 0.72 | 0.63 | 15.5 | 0.27 | 0.27 |
| 4 | 0.28 | 0.88 | 0.80 | 0.07 | 0.79 | 0.60 | 14.7 | 0.24 | 0.26 |

Standard errors: 0.01 in a share, 0.008 in a rate. Sensitivities at r = 4, θ = 0.45: the share
stays 0.99 to 1.00 and the record energy 0.62 to 0.67 for every rate change; the efficiency
runs from 0.36 (slow ejection) to 0.89 (fast ejection, or no phonon loss), and the number of
cycles from 2 (fast phonon loss or slow re-emission) to 360 (no phonon loss).

## Table 17. Reading B, the indivisible quantum, same channels

| r | share pixel 1 at θ = 0.6, 0.45, 0.36, 0.28 | Born | doubles | efficiency | record energy | cycles |
|---|---|---|---|---|---|---|
| 1 | 0.52, 0.50, 0.50, 0.50 | 0.50 | 0 | 0.57 to 0.59 | 1.000 | 0.6 |
| 2 | 0.67, 0.67, 0.65, 0.67 | 0.67 | 0 | 0.58 to 0.59 | 1.000 | 0.6 |
| 4 | 0.80, 0.78, 0.80, 0.82 | 0.80 | 0 | 0.58 to 0.60 | 1.000 | 0.6 |

All within 3 standard errors of the normalized Born shares. The efficiency is the fill factor
times the ejection branch, thinned by the escape fraction on each retry.

## Verification (`results/verification_release.txt`)

| Check | Result |
|---|---|
| Energy closure, 31 configurations | ≤ 2 × 10⁻¹⁵, pass |
| Doubles exactly 0 at θ = 0.6 (both readings) and at every θ in reading B | pass |
| Reading B's record energy exactly 1 | pass |
| Reading A's shares at r = 4, θ = 0.45 and r = 1 against P19 | 0.996, 0.494, 0.521 against 0.996, 0.48, 0.50, pass |
| Reading B's shares against Born, 12 configurations | pass |

## The cases named in the addendum, resolved

- **C12 held:** reading A's records carry 0.45 to 0.72 of the quantum; a transition-edge sensor
  reports the whole photon energy per absorption. Ledger L21.
- **C13 held:** reading A's share is a step in r, 1.00 at r = 4 for the silicon gap ratio, against
  Born's 0.80. Ledger L21.
- **C14 held:** reading B satisfies every constraint by assuming the indivisible quantum. Ledger
  L22.

## What Part G says about the third step

The release cycle is the first dynamics in this program that keeps one record per quantum
without a stop across pixels and without assuming the empty pixel silent: energy conservation
does it, whenever the gap is more than half the photon energy, and JB's re-emission channel is
what lets a sub-threshold holding go back to the wave instead of latching. But it does so only
in the reading where the deposit is divisible, and in that reading the share is not Born and the
record is not the whole quantum. Both of those are measured facts against it. The reading in
which the share is Born and the record is whole is the reading in which the quantum is deposited
indivisibly, and that is the postulate the program set out to derive. The two readings cannot be
combined: a divisible deposit that always ends in a whole-quantum record is a deposit that draws
the rest from the wave at the moment of ejection, which is the stop across pixels again.

## Addendum 5b, 2026-09-08 (late night): re-emitted photons leave (`release_retry.py escape`)

JB's objection to D30, registered as D30b before the run: a re-emitted photon does not interact
with the detectors again. A sub-threshold holding re-emitted as a photon is below the gap, and
silicon is transparent below its gap; an above-gap photon from a surface atom leaves into the
half-space away from the material. The accumulation of D29 was a two-photon process. Outputs:
`results/release_escape.csv`, `results/verification_release_escape.txt`.

**Read this first.** 23 of 24 checks passed. One failed on my closed form: I wrote that at r = 2
the bright pixel records at every gap ratio, but its share is 0.570, below the 0.6 of 660 nm, so
at that wavelength reading A makes no record at all at r = 2. The run is right; the FAIL stays in
the verification file.

| reading | r | θ = 0.6 | 0.45 | 0.36 | 0.28 |
|---|---|---|---|---|---|
| A: efficiency | 1 | **0** | **0** | 0.68 | 0.68 |
| A: efficiency | 2 | **0** | 0.45 | 0.44 | 0.70 |
| A: efficiency | 4 | 0.44 | 0.43 | 0.43 | 0.43 |
| A: share pixel 1 | 1 | none | none | 0.51 | 0.50 |
| A: share pixel 1 | 2 | none | 1.000 | 1.000 | 0.50 |
| A: share pixel 1 | 4 | 1.000 | 1.000 | 1.000 | 1.000 |
| A: doubles | 1 | 0 | 0 | 0.19 | 0.18 |
| A: doubles | 2 | 0 | 0 | 0 | 0.18 |
| B: share pixel 1 (Born 0.50, 0.67, 0.80) | 1, 2, 4 | 0.52, 0.69, 0.80 | 0.52, 0.66, 0.81 | 0.51, 0.68, 0.79 | 0.48, 0.66, 0.81 |
| B: efficiency, doubles, record energy | any | 0.35 to 0.39, 0, 1.000 | | | |

Energy closes exactly in every configuration (nothing is ever re-deposited).

**What it says.** Without the retry, reading A detects nothing at a balanced split for every
wavelength above 470 nm and nothing at r = 2 above 553 nm, and where it detects, the brighter
pixel takes every record. The retry of D30 was the only thing that gave reading A an efficiency
at a balanced split, and JB's objection removes it on physical grounds. Reading B loses only
efficiency (0.37 against 0.58), keeps Born, exclusivity and the whole quantum per record. The
physical retry that does exist, recycling of whole-quantum photons inside a thick absorber, is
reading B's kind, and silicon does it weakly (radiative lifetime of microseconds against
picosecond thermalization). Ledger L23.
