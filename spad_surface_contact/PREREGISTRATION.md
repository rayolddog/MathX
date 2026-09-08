# Pre-registration: two SPAD pixels, surface contact by least rotation, the record as an induced current

**Date:** 2026-09-08. Written before any run; `RESULTS.md` cites this file and does not change the
definitions below. Predicted numbers were computed in scratch from closed forms, the second study's
propagated packet and a small Monte Carlo before the script was written; the script must reproduce
them. **Exploration:** `MathX/spad_surface_contact/`. **Reuses:** the free-flight code and the
1+1-D Dirac conventions of `../second_spinor_phase_race/bohm_dirac_flight.py` (ħ = c = m = 1,
H = σ_x p + σ_z m, the standard Dirac representation, velocity ψ†σ_x ψ / ψ†ψ).

## 0. Read this first

**The model (JB, 2026-09-07 and 2026-09-08), in his words as far as possible.** The phase
difference between the positive and negative spinor persists in the wave function after the
initial interaction with an electron bound to an atom. The first spinor of both the incident
particle and the interacting electron rotate to a matched state; the second spinor of the incident
particle is coupled to the first, so the phase difference is matched. The incident particle
interacts with the atom that has the most similar phase, minimizing the rotation needed to match
the spinor phases (a least-action reading). The detector is an array, beginning with two
detectors, each with an aperture within which many atoms may interact. Begin with a model that
oversimplifies: the incident particle interacts with an electron bound to an atom on the surface;
no depth of penetration, no scatter before absorption, no unbound electrons (later extensions).
The next step is to model how bumping the energy orbit of the electron results in a record of the
location of interaction. Do not presume the wave function diffuses through the medium of the
detector: the physics of lightning suggests a current following the path of least resistance,
which may be where the Dirac–Bohm trajectory really comes into play. Look at it as an induced
current. Silicon SPAD only; film and the nanowire are set aside.

**What this study can and cannot show, stated before the run.** In the trajectory reading (A1)
the particle is at one place, so the second pixel is silent and the pixel statistics are the
arrival statistics, which are Born's by Bohm's equivariance (second study, C6; ledger L5). One
record per quantum and the Born share are therefore supplied by the reading, not by the phase
mechanism or by the induced current. What the mechanism decides here is which surface atom, how
much rotation, whether a pair forms, and what the record inside the pixel looks like. The field
reading of the same rule is computed beside it as the cost case (C1), so that the reader sees
which reading carries the statistics.

**Choices the description leaves open, fixed here before the run. Audit these first.**

- **A1, transport.** The wave is transported by its own Dirac current: the free Dirac equation to
  the detector plane, and the particle is where its Bohm–Dirac trajectory puts it. No diffusion
  is used anywhere (the third study's diffusive transport is not used). The one in-medium
  transport is the induced current of D7.
- **A2, contact at the surface.** Each pixel is a line of bound electrons on its surface; the
  arriving particle interacts with one of the atoms within its reach. No penetration, no scatter.
- **A3, the selection rule.** Primary, JB's statement of 2026-09-07: *least rotation*, the atom
  in reach whose second-spinor mismatch |Δφ_i| is smallest. Comparator, JB's choice of 2026-09-07
  for the second study: R2, one atom sampled with probability ∝ cos²(Δφ_i/2). Both run.
- **A4, the deposit.** D4a (primary): the winning atom's bound electron takes the whole quantum
  E_p. D4b (sensitivity): it takes the interference share w E_p with w = cos²(Δφ/2); the
  remainder is not deposited and its fate is a ledger item, not a definition.
- **A5, the pair.** A bumped electron is a real electron–hole pair when the deposit reaches the
  gap: E_dep ≥ θ E_p with θ = 0.45, silicon at 500 nm, taken from the record's Paper 1 §6.1
  (E_gap/E_photon at 500 nm). Under D4a a pair always forms.
- **A6, the record as an induced current, in two separable parts.** *Where* the charge goes is a
  transfer of charge down a lattice of random local resistances, cell by cell, by a rule: greedy,
  all charge to the least-resistance neighbour (the lightning leader; primary); a dielectric-
  breakdown split ∝ 1/resistance (η = 1); equal split (diffusion, the comparator JB rejected,
  kept so the difference is a number). *Whether* it becomes a record is a Galton–Watson avalanche
  with a loss and an ionization probability per row, independent of the transfer rule, because in
  a linear transfer the gain is a common factor. The record forms when the avalanche reaches the
  electrode with at least Q_th carriers.
- **A7, the empty pixel.** The wave at a pixel where the particle is not deposits nothing (the
  trajectory reading; ledger L5). This is the single assumption that closes the third step here,
  and it is named as such in C5.
- **A8, the landscape.** The resistance landscape is fixed per pixel for the whole run (one
  device; primary) and drawn fresh per quantum as the sensitivity. The bulk's second-spinor
  phases are drawn fresh per quantum (second study, A5).
- **A9, units.** Lengths in Compton lengths, times in Compton times, as in the second study.
  Nothing in the pixel model carries a physical scale except θ; the lattice cell is the unit.

## 1. Definitions (fixed)

- **D1, particle and atoms.** As the second study's D1 to D3: offsets δ_p and δ_i, contact shifts
  both spinors of the particle, mismatch Δφ_i = δ_p − δ_i, uniform on the circle for a random bulk.
- **D2, flight.** The second study's Part C packet: two Gaussian lobes of width σ = 1 at
  x = ±4, upper component only at t = 0, intensity ratio r ∈ {4 (primary), 1}, propagated exactly
  in momentum space to t_d = 4; M = 4000 trajectories per r, started from ρ(x, 0) and integrated by
  RK4 (dt = 0.01) against the exact field. The detector plane is x at t = t_d.
- **D3, pixels.** Pixel 1 is the interval [1, 7], pixel 2 is [−7, −1]. The gap (−1, 1) and the
  region |x| > 7 hold no atoms: a quantum arriving there makes no record (the fill factor).
- **D4, surface atoms in reach.** The n atoms in reach of an arrival x_a are the n points of a
  lattice of spacing a_n = 1/(n − 1) (a₁ = 1) centred on the lattice point nearest x_a; they span
  x_a ± 0.5 as in Part C. n ∈ {1, 3, 11, 21 (primary), 101, 1001}: "the number of atoms
  available for interaction" (JB). Atoms outside the pixel are excluded.
- **D5, contact.** Least rotation: i* = argmin_i |Δφ_i|. R2: i sampled ∝ w_i = cos²(Δφ_i/2).
  Recorded: the rotation used, |Δφ_{i*}|, and w_{i*}. The contact position is x_{i*}.
- **D6, deposit and pair.** D4a: E_dep = E_p. D4b: E_dep = w_{i*} E_p. Pair iff E_dep ≥ θ E_p.
- **D7, the induced current.** Each pixel is a lattice of W = 121 columns (spacing 0.05, so
  column c sits at x = 1 + 0.05c in pixel 1 and x = −7 + 0.05c in pixel 2) and H = 60 rows; the
  contact column is the column nearest x_{i*}. Cell resistance ρ_c = exp(ξ), ξ ~ N(0, 1),
  independent per cell. The charge starts as 1 in (row 0, contact column). From row k to k + 1 the
  charge in column c moves to c − 1, c, c + 1 (columns outside the lattice excluded) by the rule:
  greedy, all to the least-resistance destination; DBM, split ∝ 1/ρ of the destinations;
  diffusion, equal split. The path time is T = Σ_rows Σ_c n_{k+1}(c) ρ_{k+1}(c), the
  charge-weighted resistance crossed. At the electrode (row H): centroid c̄ and rms width of
  n_H(c); offset = c̄ − contact column.
- **D8, breakdown.** From one carrier at row 0, per row each carrier is lost with p_l = 0.1, else
  ionizes with p_i = 0.5 (two carriers). Record iff the count at row H is ≥ Q_th = 100.
- **D9, what a record is, and the statistics.** A record is a pair (D6) that breaks down (D8). Per
  quantum: pixel of arrival (1, 2, gap, outside), contact atom, rotation, w, pair, breakdown,
  record, contact column, electrode centroid, width, offset, T under each transfer rule.
  Efficiency = records / quanta. Share = records in pixel 1 / records. Doubles = quanta with
  records in both pixels (zero by A7, counted anyway). Miss = pixel arrivals with no pair.

## 2. Predictions (from scratch; to be reproduced)

- **P1, arrivals and share.** r = 4: pixel 1 0.6865, pixel 2 0.1680, gap 0.0790, outside 0.0665
  of ρ(x, t_d); share 0.8033 (Born 0.800; the difference is the lobes' spreading into the gap and
  outside), fill 0.8545. r = 1: 0.4254 each, share 0.500, fill 0.8507. Arrivals match these
  within 3 standard errors of 4000. The share is the same under least rotation and R2, at every n,
  under D4a and D4b, because the pixel is decided before the race. Doubles 0.
- **P2, rotation used.** Least rotation: mean π/(n + 1) = 90°, 45°, 15.0°, 8.2°, 1.8°, 0.18° for
  n = 1, 3, 11, 21, 101, 1001. R2: 90° at n = 1, tending to 53.5°.
- **P3, the miss under D4b.** Least rotation: (1 − φ₀/π)^n with φ₀ = 2 arccos √θ = 95.7°, that
  is 0.468^n: 0.468, 0.103, 2.4 × 10⁻⁴, 1.2 × 10⁻⁷, 0, 0. R2: 0.470, 0.234, 0.167, 0.158,
  0.153, 0.152, tending to 0.152: a phase-limited efficiency ceiling of 0.85 at every n. Under
  D4a the miss is 0 for both.
- **P4, breakdown.** Extinction 0.2222 (Galton–Watson, mean offspring 1.35), P_b = 0.778.
  Efficiency, primary (D4a, r = 4): 0.8545 × 0.778 = 0.665; D4b with R2 at n = 21:
  0.8545 × 0.842 × 0.778 = 0.560.
- **P5, the induced current (fixed landscape, 60 rows).** Greedy: one channel, electrode width 0,
  offset sd ≈ 5.3 cells, and the channels of one device end in about 12 of 121 columns. DBM
  η = 1: width ≈ 6.1, offset sd ≈ 1.0. Diffusion: width √(2H/3) = 6.32, offset sd ≈ 0.4. Time
  per row: greedy 0.56 (E[exp min of three N(0, 1)]), diffusion 1.65 (E[exp N]), DBM between.
  Fresh landscape: greedy offset sd ≈ 6.1, width 0.
- **P6, location.** The correlation between contact column and electrode centroid is ≥ 0.96 under
  every rule, because the contact's spread across the pixel (sd about 25 cells) exceeds the
  channel's wander; the greedy record misplaces the contact by about 5 cells = 0.27 Compton
  lengths, the diffusive one by 0.4 cells but spread over 6.3.

## 3. The cases that cost the model, named so they are not skipped

- **C1, the field reading of least rotation.** With every surface atom of both pixels in play and
  the density under each atom as its available energy, "least rotation" has two readings and
  neither is Born: pure least rotation (argmin |Δφ| over all atoms, intensity ignored) gives the
  atom-count share, 0.500 at r = 4 with 121 atoms per pixel; deposit-weighted least rotation (R3,
  argmax ρ_i w_i) gives the brighter pixel in essentially every trial (second study, Table 2). R2
  in the same reading gives 0.80. All three are computed here (4000 trials, 121 and 5 atoms per
  pixel). This is the computation that shows where Born comes from in this study.
- **C2, D4b's efficiency ceiling under R2.** 0.85 at any n is a claim about the detector's quantum
  efficiency that no material property sets. Ledger.
- **C3, the undeposited remainder under D4b.** (1 − w) E_p goes nowhere in the model. Ledger.
- **C4, the record's location.** The greedy channel decouples the electrode location from the
  contact by about 5 cells and funnels into device-fixed columns; a SPAD reports no position
  inside a pixel, so this is not testable, and is recorded as what the model says.
- **C5, the third step closes by A7 alone.** If the empty pixel's wave is given any ability to
  deposit, the doubles of the second study's Part G and of the third study return (on record).
  This study does not test that; it assumes it away and says so.

## 4. What will be recorded

`results/`: `trials.csv` (the primary configuration per quantum), `summary.csv` (rule × deposit ×
n × r), `field_reading.csv` (C1), `current.csv` (rule × landscape), `verification.txt`,
`run_info.txt`; `fig1_arrivals_pixels.png`, `fig2_rotation_miss.png`, `fig3_current.png`.
`RESULTS.md` states the uncertainty and the assumptions at the top, then the tables. Departures
from quantum mechanics or experiment go to `LEDGER.md`, dated, continuing the second study's
numbering.

## 5. Verification (must pass before the tables are read)

Arrivals in each region against ∫ρ within 3 standard errors; Galton–Watson extinction against
0.2222 within 3 standard errors; diffusion width against √(2H/3) within 5 %; transfer conserves
charge to 10⁻¹⁰; least-rotation mean mismatch against π/(n + 1) within 3 standard errors at every
n; D4b miss against 0.468^n within 3 standard errors at n ≤ 3.

## 6. Method note

Definitions and predictions were fixed in this file before the script was written to disk and
run. Any later change to a definition is to be made as a dated addendum below, not by editing
§0 to §3.

---

## Addendum, 2026-09-08 (afternoon): the empty pixel, and the bulk's phase distribution

Written before `empty_pixel.py` was run, at JB's instruction of 2026-09-08: "run the empty-pixel
test with the same contact and current", and his concern that the atoms' phases were drawn from a
flat distribution and that a formed silicon crystal would change that. §0 to §3 above are
unchanged; A7 is lifted for Part B only. Predicted numbers were computed in scratch from the main
run's `trials.csv`, the propagated packet at r = 2, and a small Monte Carlo.

**What was used for the phases, stated plainly.** In the second study and in the main run above
every atom's second-spinor offset δ_i is drawn uniformly on the circle, independently of every
other atom, fresh for each quantum (second study D2 and A5; D1 above). That is the flat
distribution JB asks about. It is the maximum-entropy distribution of a phase: there is no "more
randomized" distribution of an angle. What a crystal changes is not the randomness but (i) the
marginal, if the bound electrons share a preferred phase, and (ii) the independence, if the phase
advances regularly from atom to atom. Both are run in Part C.

### Part B: the empty pixel with the same contact and current

- **D10, the empty pixel's contact.** At the pixel where the particle is not, a contact point is
  sampled from ρ(x, t_d) restricted to that pixel (the wave arriving as it would if the particle
  were there), the n = 21 atoms in reach are drawn with fresh phases, and the least-rotation atom
  is the contact (D4, D5). R2 is not run in Part B; the rule is not what is tested.
- **D11, the empty pixel's deposit, two readings, both run.** E1 (share): the wave's norm in that
  pixel, s_k = ∫_pixel ρ dx, times E_p. E2 (whole): E_p, the empty lobe acting as if it carried the
  particle. The particle's pixel keeps D4a (whole quantum). A pair forms at either pixel iff its
  deposit ≥ θ E_p. θ ∈ {0.45, 0.407, 0.361, 0.316, 0.28}: silicon (E_gap = 1.12 eV) at 500, 450,
  400, 350 and 310 nm. r ∈ {1, 2, 4}; pixel norms s at t_d: r = 1: 0.4254 each; r = 2: 0.5699
  and 0.2819; r = 4: 0.6865 and 0.1680.
- **D12, the channels and the drain.** Both contacts launch their induced current at t_d on
  their own pixel's lattice (D7; greedy primary, DBM and diffusion beside; fixed landscapes, one
  per pixel). A pixel is a candidate when it has a pair and its avalanche breaks down (D8,
  independent per pixel). Channel k connects at T_k, the resistance crossed. The first candidate
  to connect drains the wave (the lightning return stroke) and the other candidate's deposit is
  withdrawn, unless its own channel connects within τ_d of the first, in which case both record.
  τ_d ∈ {0, 1, 2, 5, 10, 20, ∞} in the unit of the resistance crossed (the greedy channel takes
  about 35). τ_d = 0 is an instantaneous stop across the pixels; τ_d = ∞ is no drain.
- **D13, head start, two cases.** H0 (primary, the cost case): the channel time does not depend on
  the deposit, as in D7. H1 (sensitivity): the channel time scales inversely with the deposit,
  T_k → T_k E_p / E_dep,k, so under E1 the empty pixel's channel is slower by 1/s_k; under E2 both
  deposits are E_p and H1 is H0.
- **D14, statistics.** Per quantum: no record, a record at the particle's pixel, a record at the
  empty pixel (a *wrong-pixel* record), or a double. Share = records at pixel 1 / all records (a
  double counts once at each pixel). Efficiency = quanta with at least one record / quanta.

**Predictions.**

- **P7, E1 with H0.** At r = 4 and r = 2 the empty pixel's share (0.168; 0.282) is below every
  θ except 0.28 at r = 2, so Part B reproduces the main run: share 0.80 (r = 4), doubles 0. At
  r = 1 the empty pixel pairs for θ ≤ 0.407 (450 nm and below), and at r = 2, θ = 0.28, both
  pixels pair whichever holds the particle. Then the two channels race with no head start: the
  doubles per quantum are (in-pixel fraction 0.85) × P_b² × P(|T₁ − T₂| ≤ τ_d) ≈ 0.515 ×
  {0, 0.21, 0.40, 0.81, 0.98, 1.00, 1.00} for greedy at τ_d = 0, 1, 2, 5, 10, 20, ∞; DBM
  {0, 0.25, 0.47, 0.89, 1.00, 1.00, 1.00}; diffusion {0, 0.11, 0.22, 0.52, 0.86, 1.00, 1.00}. The
  greedy channel, with the smallest spread of connection times, gives the most doubles at a given
  drain time. Where both pixels are candidates the record's pixel is the device's faster channel:
  on this run's landscapes P(T₁ < T₂) ≈ 0.34 (greedy; pixel 2's landscape is faster by 1.5),
  0.50 (DBM), 0.53 (diffusion). At r = 2, θ = 0.28, H0, τ_d = 0 the share is therefore about
  0.173 × 0.669 + 0.605 × 0.34 + 0.173 × 0.34 ≈ 0.38 for greedy (Born 0.669): the intensity
  ratio no longer sets the share.
- **P8, E2 with H0.** The empty pixel always pairs, so at every r and θ the record's pixel is
  independent of the particle's pixel: share ≈ 0.40 (greedy), 0.50 (DBM), 0.53 (diffusion) at
  τ_d = 0, the device's channel asymmetry, and the wrong-pixel rate ≈ 0.85 × [0.173 + 0.605 ×
  P(empty first)] ≈ 0.30. Doubles per quantum 0.515 × P(|ΔT| ≤ τ_d), reaching 0.515 at τ_d = ∞.
  Born is destroyed at every r.
- **P9, E1 with H1.** The empty pixel's channel is slower by 1/s: at r = 1, θ ≤ 0.407, by 2.35
  (greedy 82 against 35), at r = 2, θ = 0.28, by 3.5 or 1.75. The particle's pixel always connects
  first; doubles 0 for τ_d ≤ 20 (the gap in connection times is at least 26); share Born (0.5 at
  r = 1; 0.669 at r = 2); the wrong-pixel rate is P_b(1 − P_b) × 0.85 × P(empty pairs) ≈ 0.15
  where the empty pixel pairs (the particle's own avalanche failed and the empty one succeeded).
- **P10, E2 with H1** equals P8.

**The cost cases, named.** C6: under H0 the lightning channel does not silence the second pixel;
it makes the record's pixel a property of the device and produces doubles at the rate the drain
time allows. C7: under E1 the model's exclusivity at a balanced split rests on the gap: it holds
at 500 nm because 0.425 < 0.45 and fails at 450 nm and below (heralded 505 nm silicon data give
g²(0) = 0.0023, which is inside the window; no data below 470 nm were located by the record's
2026-09 search). C8: E1 with D4a gives the two pixels E_p(1 + s_empty) in total: energy is not
conserved in the model as written. Ledger.

### Part C: the bulk's phase distribution

- **D15, bulk models** at the particle's pixel, n = 21, least rotation, the weighted deposit D4b
  (the reading in which the phase leaves a trace; under D4a the winner takes the whole quantum and
  no bulk model changes any statistic), θ = 0.45, r = 4, the main run's arrivals.
  B0, random: δ_i uniform on the circle, iid (as run). B1, ordered: δ_i ~ von Mises(μ = 0, κ),
  κ ∈ {1, 4} (second study, D7). B2, locked crystal: δ_i = 0 for every atom. B3, periodic
  crystal: δ_i = q·i with i the atom's lattice index, q = 2π/φ² = 137.5° (incommensurate) and
  q = 120° (commensurate, three phases). The particle's offset δ_p ∈ {0, 45, 90, 135, 180}°
  (for q = 120°: {0, 30, 60}°).
- **P11.** B0: 8.2°, miss 0, at every δ_p. B1 κ = 1: 3.8°, 5.1°, 10.1°, 19.2°, 25.0° with miss
  ≤ 0.003; κ = 4: 1.7°, 5.3°, 31.9°, 74.6°, 108.4°, miss 0, 0, 0, 0.11, 0.79. B2: rotation
  |δ_p| exactly, miss 0 for |δ_p| ≤ 95.7° and 1 beyond; efficiency over a uniform δ_p 0.53.
  B3 golden: 4.5°, miss 0, independent of δ_p. B3 q = 120°: 0°, 30°, 60° at δ_p = 0, 30, 60°,
  miss 0. So an ordered or locked crystal turns the SPAD's efficiency into a function of the
  particle's own offset, ½[1 + cos δ_p]-like (the second study's C5 in SPAD language), and a
  periodic crystal behaves like a random bulk with a hard bound on the rotation.

**Verification for the addendum.** The empty pixel's contact points follow ρ restricted to the
pixel (max CDF difference ≤ 0.03); E1/H0 at r = 4, θ = 0.45 reproduces the main run's share
within 3 s.e. and gives doubles 0; doubles at τ_d = 0 are 0; doubles at τ_d = ∞ under E2 equal
the both-candidates fraction (0.515 ± 0.008); B0 rotation against π/22; B2 rotation equals
|δ_p| to 10⁻¹²; B3 q = 120° rotation at δ_p = 30° equals 30° to 10⁻⁶.

**What will be recorded.** `results/empty_pixel.csv` (r × θ × E × H × rule × τ_d), `results/
bulk_phase.csv`, `results/verification_addendum.txt`, `fig4_empty_pixel.png`, `fig5_bulk_phase.png`;
a dated section in `RESULTS.md`; ledger entries L10 onward.

**D12b, the field (added before the run, on JB's question of 2026-09-08 whether the SPAD's field,
which amplifies the electron current, straightens the path).** A SPAD's depletion region carries
the applied reverse-bias field, of order 10⁵ V/cm in silicon, perpendicular to the junction; it
accelerates the carriers to impact ionization (the gain) and makes the drift nearly straight, with
lateral motion only by diffusion and scattering. In the lattice this is a fourth transfer rule,
*greedy_field*: the greedy rule with the two lateral destinations' resistances multiplied by
(1 + b), b = 3, for the choice only; the time counts the chosen cell's own resistance. b = 0 is
the greedy rule of D7. Run in Part B beside the other three rules, and on the particle's
contacts at r = 4 for its own width, offset and time.

**P12.** With b = 3 a lateral move is taken in 0.26 of rows against 0.67 (the chance that a
side cell's resistance times four is still the smallest); the offset sd falls to about 3.2 cells
from 5.2; the resistance crossed per row rises to about 0.88 from 0.57, because the straighter
path cannot take the cheapest cell, and its spread rises (sd of T about 5.9 against 3.7); so the
doubles at a given drain time fall: P(|ΔT| ≤ τ_d) ≈ 0.09, 0.19, 0.43, 0.76, 0.99 at τ_d = 1, 2, 5,
10, 20 against 0.21, 0.40, 0.81, 0.98, 1.00 for b = 0 on the main run's landscapes. The field
straightens the channel and, by widening the spread of connection times, makes the two-pixel
race less often a tie; it does not silence the second pixel.
