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

---

## Addendum 2, 2026-09-08 (evening): Part D, the crystal as a transition coherence with thermal dephasing

Written before `transition_coherence.py` was run, at JB's instruction of 2026-09-08 ("run the
crystal as a transition coherence with thermal dephasing"), after the exchange on whether atoms
have a physical phase: a stationary state has none that can be read; a superposition of two
levels has a relative phase that turns at the transition frequency, and that is what MRI reads
as the Larmor phase. Part D gives the bulk that phase. §0 to §3 and the first addendum are
unchanged. Predicted numbers were computed in scratch from the closed forms and a Monte Carlo of
2 × 10⁵ draws.

**The mapping, in MRI language.** A preparation at t = 0 gives every atom of pixel k the same
coherence phase μ_k (the 90° pulse). Afterwards each atom's phase drifts: by a random walk,
irreversible, which is T2 (homogeneous dephasing); and by a static frequency offset of its own
site, refocusable by an echo, which is T2* (inhomogeneous dephasing). The magnitude of the
ensemble's mean phasor, the transverse magnetization, is R(t) = exp(−t/T2) for the first,
exp(−t²/2T2*²) for the second, and their product for both. R is the second study's mean
resultant R_κ, now a function of the delay between the preparation and the particle's arrival.
No echo is run.

**The physical scale, and the cost case named first (C9).** The transition coherence of a
bound electron in silicon at 2.5 eV turns at about 4 × 10¹⁵ rad/s and dephases, at room
temperature, in about 10 to 100 fs (semiconductor interband coherence; the record's Paper 1 puts
silicon thermalisation at 10 fs to 1 ps). At thermal equilibrium the coherence is not merely
dephased but absent: its amplitude carries the Boltzmann factor exp(−ħω/kT) ≈ e⁻⁹⁷. So Part D
describes a detector whose coherence was *driven* within about T2 of the arrival, not a passive
SPAD, and a passive SPAD is the random bulk B0 of every earlier run. Nuclear spins have T2 of
milliseconds to seconds, but zero transverse coherence without a pulse, and the deposit in this
model is electronic. Ledger.

- **D16, the phases.** Atom i of pixel k at delay t has δ_i(t) = μ_k + s(t) Z_i with Z_i standard
  normal, independent per atom and drawn fresh per quantum (a fresh draw of static offsets is the
  same as a fresh set of atoms). s_hom(t) = √(2t/T2); s_inh(t) = t/T2* with T2* = T2;
  s_both = √(s_hom² + s_inh²). t/T2 ∈ {0, 0.03, 0.1, 0.3, 1, 3, 10}. R(t) = exp(−s²/2).
- **D17, one pixel.** n = 21 atoms exactly (pure phase draws, no pixel geometry, so the edge excess
  of Part C is absent), μ = 0, least rotation, the weighted deposit D4b with θ = 0.45. The particle's
  offset δ_p ∈ {0, 45, 90, 135, 180}° and uniform. Recorded: R measured, rotation used, the pair
  rate, the winner's weight.
- **D18, two pixels.** The r = 4 flight (D2, fresh seed), the trajectory reading, the pixel
  geometry of D3 and D4, least rotation, D4b, breakdown by D8. Preparations: *same*, μ₁ = μ₂ = 0;
  *opposite*, μ₁ = 0, μ₂ = π (two devices pulsed in antiphase). δ_p fixed at 0 and uniform.
  Recorded: pair rate at each pixel, efficiency, share of records at pixel 1 against the arrival
  share 0.803 and Born 0.800.

**Predictions (P13).**

- R(t) for hom: 1, 0.970, 0.905, 0.741, 0.368, 0.050, 0; inh: 1, 1.000, 0.995, 0.956, 0.607,
  0.011, 0; both: 1, 0.970, 0.900, 0.708, 0.223, 0.001, 0.
- One pixel, hom, rotation used at δ_p = 0: 0, 0.8°, 1.5°, 2.5°, 4.6°, 7.5°, 8.2°; at 90°: 90°,
  63.5°, 41.7°, 16.3°, 8.5°, 8.2°, 8.2°; at 180°: 180°, 149°, 124°, 83°, 23°, 9.1°, 8.2°. Pair rate
  at 180°: 0, 0, 0.021, 0.71, 0.999, 1, 1; at δ_p = 0 and 90°: 1 throughout. Uniform δ_p: 0.532,
  0.679, 0.802, 0.964, 1.000, 1, 1. inh and both: the same numbers at equal s, so the same curves
  against s and different against t.
- Two pixels, opposite preparation, fixed δ_p = 0: share 1.000, 1.000, 0.995, 0.851, 0.803, 0.803,
  0.803 (the pixel in phase with the particle takes every record until t ≈ 0.3 T2, Born by
  t ≈ T2); uniform δ_p: 0.803 at every t, because the two pixels' pair rates are equal by symmetry
  (I first estimated 0.82 by hand and the Monte Carlo corrected it before the run). Same
  preparation: 0.803 at every t. Efficiency = 0.847 × 0.778 × pair rate: uniform δ_p 0.35, 0.45,
  0.53, 0.64, 0.66, 0.66, 0.66; fixed δ_p = 0 with the same preparation 0.66 throughout; with the
  opposite preparation 0.53, 0.53, 0.53, 0.62, 0.66, 0.66, 0.66.
- So the trace of the phase mechanism is a Born deviation and an efficiency that both depend on
  the delay since preparation, gone after about one T2. In this model Born statistics are what a
  dephased detector shows.

**Verification.** R measured against exp(−s²/2) within 0.005; rotation at t = 0 equals |δ_p| to
10⁻¹²; rotation at t = 10 against π/22 within 3 s.e.; uniform-δ_p pair rate at t = 0 against
0.5325 within 3 s.e.; two-pixel share at t = 10 (opposite, fixed δ_p) against the arrival share
within 3 s.e.

**What will be recorded.** `results/dephasing_single.csv`, `results/dephasing_pixels.csv`,
`results/verification_dephasing.txt`, `fig6_dephasing.png`; a dated section in `RESULTS.md`;
ledger L13.

---

## Addendum 3, 2026-09-08 (night): Part E, the bound electron's second spinor, and the deposit weight with unequal amplitudes

Written before `bound_spinor.py` was run, at JB's instruction of 2026-09-08: the second Dirac
spinor is a complete complex wave, so |A| ≠ |B| in the interference weight; and the ground state
of an electron bound to an atom is not the ground state of a free electron, so revise the
calculation on that basis. Predicted numbers were computed in scratch on the same grid.

**What is computed.** A Dirac electron bound in a well, in the record's 1+1-D conventions
(ħ = c = m = 1, H = σ_x p + σ_z m, the standard representation), so that both readings of "the
second spinor" can be measured on the same state: JB's reading, the lower row of the Dirac spinor
(the "original Dirac spinors"); and A4's reading of the second study, the negative-energy
component in the free-particle expansion. For each: its amplitude relative to the first spinor,
its phase relative to the first spinor, whether that phase is fixed or free, and how it rotates in
time. Then the deposit weight of A1 with the two amplitudes unequal.

- **D19, the bound electron.** V(x) = −U₀ exp(−x²/2L²) added to H as a potential energy. Grid
  1024 points at spacing 0.05. The ground state is the lowest eigenvalue in (−m, m). Four wells,
  (U₀, L) = (0.05, 4), (0.2, 2), (0.6, 1), (1.0, 0.5), with bindings from 0.03 to 0.35 of the rest
  energy. A valence electron in silicon is bound by a few eV, 10⁻⁵ of the rest energy, shallower
  than any well run; the shallow end of the run gives the scaling to extrapolate.
- **D20, measured.** E_b. Lower-row fraction f_L = ∫|χ|²/∫|ψ|² and amplitude a_L = √f_L. The
  pointwise phase of χ relative to φ. The identity χ = −i φ′/(E_b + m − V) (the lower row of the
  Dirac equation in this representation). Free ±E content f₋ = Σ|c₋|²/Σ(|c₊|² + |c₋|²) with c±
  from Part C's projectors, amplitude a₋ = √f₋. The offset arg(c₋/c₊) over p: its mean and its
  resultant. The eigenvector residual |Hψ − E_bψ|, which is the statement that c₊ and c₋ both
  rotate as e^(−iE_b t), together, where a free particle's rotate as e^(∓iE(p)t), apart.
- **D21, the incident particle.** Its second-spinor amplitude in A4's reading is b = sin(θ/2) for
  a lag θ (second study, L4); a physical free electron has θ = 0 and b = 0; a photon has no lower
  spinor (2026-09-06 ledger). In JB's reading a free electron of momentum p has b_L = p/(E + m),
  0.1 at p = 0.2 m.
- **D22, the revised weight.** With amplitudes a (bound electron) and b (incident particle),
  w = [a² + b² + 2ab cos Δφ]/(a + b)², visibility V = 2ab/(a² + b²). Every phase effect of the
  model (the second study's C5, Part D's fringe and pair rates) scales with V; the minimum of w is
  [(a − b)/(a + b)]², so the miss under D4b needs V close to 1. Table of V over the four wells'
  a_L and a₋ against b ∈ {0, a, 0.1, 0.5, 1}.

**Predictions (P17, from scratch on the same grid).**

- E_b = 0.9721, 0.8888, 0.6979, 0.6543; f_L = 0.0041, 0.0167, 0.0537, 0.0904 (a_L = 0.064, 0.129,
  0.232, 0.301); f₋ = 1.5 × 10⁻⁶, 6.8 × 10⁻⁵, 1.3 × 10⁻³, 4.7 × 10⁻³ (a₋ = 0.0012, 0.0082, 0.036,
  0.068). The lower fraction is about 0.15 of the binding in rest-energy units at the shallow end;
  the −E content falls faster, about as the cube of the binding.
- χ = −i φ′/(E_b + m − V) to 10⁻¹¹. Since φ is real for the ground state, χ is purely imaginary
  relative to it: the second spinor's phase is +90° on one side of the well and −90° on the other,
  fixed by the first spinor's gradient, identical for every identical atom. It is not a random
  variable and not a free one.
- The ±E offset is not a single number for a bound state: its resultant over p is 0.42, 0.21,
  0.12, 0.10, with well-dependent means. Both components rotate at E_b together: no counter-
  rotation. The eigenvector residual is below 10⁻⁹.
- Visibility: b = 0 gives V = 0 for every well: for a physical incident electron or a photon the
  phase mechanism is off. V = 1 only at b = a. Shallow well, a_L = 0.064: V = 0.91 at b = 0.1,
  0.25 at b = 0.5, 0.13 at b = 1; a₋ = 0.0012: V = 0.024 at b = 0.1. Extrapolated to silicon
  (binding 10⁻⁵): a_L ≈ 10⁻³, a₋ ≈ 10⁻⁸, an order-of-magnitude estimate from the scaling.

**The cost case, named (C10).** If the second spinor of a bound electron is slaved to the first
and in quadrature with it, the bulk of the second study, random and independent phases, and its
counter-rotation D1, are not what a Dirac bound state does; and if the incident particle carries
no negative-energy admixture, the weight's visibility is zero. Ledger.

**Verification.** Eigenvector residual < 10⁻⁹; χ against −i φ′/(E_b + m − V) < 10⁻¹⁰; E_b and
f_L against the predictions to 10⁻⁴; f₋ to 2 %; V table exact by construction.

**What will be recorded.** `results/bound_spinor.csv`, `results/bound_visibility.csv`,
`results/verification_bound.txt`, `fig8_bound_spinor.png`; a dated section in `RESULTS.md`;
ledger L15; a note in `L13_TEST.md`.

---

## Addendum 4, 2026-09-08 (night): Part F, the bulk rebuilt on the slaved second spinor

Written before `slaved_bulk.py` was run, at JB's instruction of 2026-09-08 ("rebuild the bulk on
Part E's slaved second spinor and run it"), and his question whether the phase difference between
the first and second spinor within one bound electron is variable or random.

**The answer to the question, from Part E, which fixes the bulk.** Within one bound electron the
difference is not random and not free. In the standard representation χ = −iφ′/(E_b + m − V),
and φ is real for the ground state, so χ is in quadrature with φ everywhere: +90° on the side
where φ falls, −90° on the side where it rises, the sign flipping at the nucleus. In three
dimensions the same holds with χ = −i σ·r̂ φ′(r)/(E + m − V): the factor −i is intrinsic, and
the direction of contact decides which spin components carry it, not the phase. So what an
incident particle meets at atom i is a second spinor of amplitude a(d), set by the distance d
of the contact from the nucleus, and phase +90° if it touches the atom on the side facing the
arrival from the left, −90° from the right. The "180° apart" of JB's question is the difference
between the two sides. All randomness of the earlier bulk is gone; what remains variable is the
side of contact and the incident particle's own offset.

**The incident particle's own offset, in the same reading.** A free positive-energy plane wave has
χ = p φ/(E + m) with p real: phase 0° for motion to the right, 180° to the left, amplitude
b = p/(E + m). A standing wave (total reflection at the surface) has χ ∝ −iφ′ with φ real, so
its offset is ±90° by position. A photon has no second spinor.

- **D23, the slaved bulk B5.** Atom i at x_i presents to an arrival at x_a the phase
  δ_i = +90° if x_a > x_i and −90° if x_a < x_i, and the amplitude a. The n = 21 atoms in reach are
  those of D4; the atoms left of the arrival present +90° and those right of it −90°.
- **D24, the weight and the choice.** w_i = [a² + b² + 2ab cos(δ_p − δ_i)]/(a + b)² (D22), w = 1
  for b = 0. Least rotation picks the smaller |δ_p − δ_i|; a tie (δ_p = 0° or 180°, both sides
  90° away) is broken by the nearest atom. Deposit D4b: pair iff w ≥ θ, θ ∈ {0.45 (500 nm),
  0.6 (660 nm)}. Breakdown D8. The r = 4 flight of D2, fresh seed.
- **D25, the grid.** a ∈ {0.064 (Part E's shallow well, the bulk value), 0.0012 (silicon, Part E's
  extrapolation)}; b ∈ {0 (a photon), 0.025 (an electron at p = 0.05 m), 0.098 (p = 0.2 m), a
  (equal amplitudes, the second study's A1)}; δ_p ∈ {0° (a running wave), 90° (a standing wave
  at the matching position), uniform (a partially reflected wave, the sensitivity)}.
- **D26, recorded.** Rotation used, the chosen atom's w, the side chosen, pair rate, efficiency,
  share of records at pixel 1 against the arrival share, and a(d) of the shallow well against the
  contact distance d (from Part E's state, recomputed).

**Predictions (P18, closed forms).**

- The chosen atom's weight is [a² + b² + 2ab |sin δ_p|]/(a + b)²: at δ_p = 0° it is
  (a² + b²)/(a + b)², which is ≥ ½ for every a, b and equals ½ at a = b: 0.522 (a = 0.064,
  b = 0.098), 0.500 (b = a), 0.596 (b = 0.025), 1 (b = 0), 0.976 (a = 0.0012, b = 0.098). At
  δ_p = 90° it is 1 for every a, b. The interference term vanishes identically for a running wave
  meeting a bound electron: cos(±90°) = 0.
- Rotation used: 90° exactly at δ_p = 0°, 0° at 90°, uniform on [0°, 90°] with mean 45° for a
  uniform δ_p. The side chosen at δ_p = 90°: always the atoms left of the arrival.
- Pair rate at θ = 0.45: 1 in every configuration (the weight never falls below ½). At θ = 0.6
  with δ_p = 0°: 0 for a = 0.064 with b = 0.098, b = a, and b = 0.025 (0.596 < 0.6); 1 for b = 0
  and for silicon's a = 0.0012 with b = 0.098. At θ = 0.6 with uniform δ_p: 0.90 (a = 0.064,
  b = 0.098), 0.87 (b = a), 0.99 (b = 0.025), 1 (b = 0 and silicon).
- Share: the arrival share (0.80) wherever pairs form, in every configuration; the choice of atom
  never moves the pixel. Efficiency 0.847 × 0.778 × pair rate = 0.66 × pair rate.
- a(d) of the shallow well: 0.011 at d = 0.5, 0.022 at 1, 0.043 at 2, 0.078 at 4, 0.100 at 6,
  0.112 at 8, tending to κ/(E_b + m) = 0.119 in the tail; phase +90° throughout d > 0.

**What the run will show, stated before it (C11).** On the slaved bulk the phase race has nothing
to select: every atom in reach presents one of two phases and the same amplitude, a running
incident wave is 90° from both, the weight is equal at every atom, and the record's pixel is the
arrival's. Randomness enters only through the incident particle's own offset, and then it selects
a side of the atom, not an atom. The one place the weight bites is a gap ratio above ½ with the
two amplitudes comparable, where a running wave makes no pair at all; with silicon's amplitudes
that never happens. Ledger.

**Verification.** Rotation at δ_p = 0° equals 90° and at 90° equals 0° to 10⁻¹²; the chosen weight
against the closed form to 10⁻¹²; uniform-δ_p rotation mean 45° within 3 s.e.; share against the
arrival share within 3 s.e. in every configuration with pair rate 1; a(d) against P18 to 0.001.

**What will be recorded.** `results/slaved_bulk.csv`, `results/slaved_a_of_d.csv`,
`results/verification_slaved.txt`, `fig9_slaved.png`; a dated section in `RESULTS.md`; ledger L16.
