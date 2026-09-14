# Pre-registration: the vacuum race, "highest interference term wins" with a per-atom vacuum term

**Date:** 2026-09-14. Written before any run, at JB's instruction "write the pre-registration and
run it", after the exchange of the same day on the modulus against the real-axis projection and on
the vacuum term as a property of the individual atom. `RESULTS.md` cites this file and does not
change the definitions below. Predicted numbers were computed in scratch (10⁵ to 2 × 10⁵ trials
per cell) before the script was written. **Exploration:** `MathX/vacuum_race/`. **Ledger:**
`LEDGER.md`, continuing the pair case's numbering (L24 onward).

## 0. Read this first

**The proposal (JB, 2026-09-14).** Within the aperture of a detector there are many atoms with
varying phases. The atom with the highest interference term with the incident particle is the
atom selected for the interaction. Each atom also carries its own vacuum term, unaffected by the
number of atoms in the aperture. Within an aperture the incident wavefunction is then the larger
factor in deciding the site; between apertures the vacuum adds a stochastic element to the
deterministic propagation of the wavefunction.

**What this run tests, stated before it.** The rule is a maximum over atoms. The pixel share it
gives at an intensity ratio r depends on the distribution of the per-atom vacuum term, and on
nothing else once the aperture holds more than a few atoms. The claim to be tested (C1) is that
exactly one family of distributions turns the maximum into a proportion, the Gumbel form in the
logarithm of the deposit, and that with it the rule is identical to the second study's sampled
rule R2, which gives Born and is the golden rule's rate law. Every other form tried is predicted to
give the brighter pixel every time once the aperture is full, because the maximum of N draws from
a bounded or Gaussian distribution concentrates as N grows while the per-atom term does not change.

**What is not modeled (C4).** JB's geometry of the vacuum, four-dimensional against fermions in
three dimensions plus spacetime, does not enter. The vacuum term here is one random number per
atom per trial, and the study tests only its distribution and how it enters the deposit
(multiplying the intensity, or interfering with the amplitude). One winner is selected by
construction, so the single record per quantum is assumed here and not tested.

## 1. Definitions (fixed)

- **D1, the pixels.** Two pixels A and B with intensities I_A = r and I_B = 1, r = 4 (primary,
  Born share 0.800) and r = 2 (sensitivity, 0.667). N atoms per pixel, N ∈ {1, 10, 100, 1000}.
- **D2, the phase weight.** w_i = cos²(Δφ_i/2) with Δφ_i uniform on the circle, independent per
  atom (the second study's D1 to D3 and A1). Control: w_i = 1, no phase mechanism.
- **D3, the deposit.** d_i = I_pixel · w_i, the interference-weighted energy available at atom i.
- **D4, the vacuum term.** One independent draw per atom per trial, from a distribution that does
  not depend on N. Six forms:
  - **V0**, none: E_i = d_i (the second study's R3, the deterministic limit).
  - **V1**, Gumbel on the deposit: E_i = d_i · exp(G_i/β), G_i standard Gumbel; β = 1 (primary),
    β ∈ {½, 2} (sensitivity). Equivalent forms: argmax(log d_i + G_i/β), or the first to fire among
    exponential clocks with rates d_i^β.
  - **V2**, Gumbel on the amplitude: E_i = √d_i · exp(G_i). The vacuum interferes with the
    particle's amplitude, so the random part of the deposit is linear in |ψ|, as the cross term
    2|ψ_p||ψ_v| cos(φ_p − φ_v) is.
  - **V3**, multiplicative Gaussian: E_i = d_i (1 + 0.5 g_i), g_i standard normal.
  - **V4**, additive Gaussian: E_i = d_i + 2 g_i, in units of I_B.
  - **V5**, additive uniform: E_i = d_i + U_i, U_i uniform on [0, 4].
- **D5, the selection.** The atom with the largest E_i among all 2N atoms takes the quantum.
  One winner per trial by construction.
- **D6, what is recorded.** The pixel share P(A); the exponent γ from P(A)/P(B) = (√r)^γ, Born
  giving 2; the fraction of trials in which the winner is the best phase-matched atom of its
  pixel; and the scaling test: the same rule at (I_A, I_B) = (40, 10), N = 1 and 10, for V1, V3,
  V4, V5, with and without the phase weight. Born depends on the ratio only, so a form whose
  share changes under the scaling reads the intensity difference, not the ratio.
- **D7, statistics.** T = 10⁵ trials per cell for N ≤ 100, 2 × 10⁴ for N = 1000; standard error
  of a share about 0.0013 and 0.003. Seed 20260914.

## 2. Predictions (from scratch)

- **P1, V0 (no vacuum).** P(A) = 0.79 at N = 1 and 1.000 from N = 10, the second study's R3 row
  (Table 2: 0.790, 1.000). The best phase-matched atom of the winning pixel wins in every trial.
- **P2, V1 at β = 1 is R2.** With the phase weight: P(A) = 0.705, 0.795, 0.800, 0.800 at N = 1,
  10, 100, 1000, the second study's R2 row (0.704, 0.794, 0.800, 0.800) within 3 s.e. Without the
  phase weight: 0.800 at every N, exactly Born, and 2/3 at r = 2. At r = 2 with the phase weight,
  N = 10: 0.663. Reason: argmax(log d_i + G_i) selects atom i with probability d_i/Σd_j, the
  Gumbel-max identity, and that is R2's sampling rule.
- **P3, the exponent is β.** At N = 100: P(A) = 0.667 at β = ½, 0.800 at β = 1, 0.941 at β = 2,
  i.e. r^β/(r^β + 1), exponent γ = 2β. Only β = 1 gives Born; the scale of the vacuum term is a
  free parameter that has to be set to Born's value.
- **P4, V2 (the vacuum on the amplitude).** P(A) = √r/(√r + 1) = 0.667 from N = 10 (0.634 at
  N = 1), exponent 1: the rule selects in proportion to |ψ|, not |ψ|².
- **P5, V3, V4, V5 concentrate.** With the phase weight: at N = 1 between 0.68 and 0.79 (V3 0.75,
  V4 0.68, V5 0.74); at N = 10 between 0.88 and 1.00 (V3 0.999, V4 0.88, V5 0.99); at N ≥ 100,
  1.000 for all three. The brighter pixel takes every record once the aperture holds a hundred
  atoms, whatever the per-atom vacuum term, because the maximum of N draws climbs to the ceiling.
- **P6, the scaling test.** V1 is unchanged at (40, 10) within 3 s.e. (0.80 without the phase
  weight at N = 1 and 10; 0.705 and 0.795 with it). V3 is unchanged under scaling (multiplicative)
  but still concentrates with N (0.93 at N = 1 without the phase weight, 1.000 at N = 10). V4 and
  V5 change under scaling: without the phase weight, V4 goes from 0.855 to 1.000 at N = 1 and V5
  from 0.97 to 1.000; with it, V4 from 0.68 to 0.79 and V5 from 0.74 to 0.79 (the V0 value, the
  noise having become negligible). An additive vacuum term reads the intensity difference.
- **P7, within the pixel.** Under V1 the best phase-matched atom of the winning pixel wins in a
  fraction about 2/N of trials (0.20 at N = 10, 0.02 at N = 100, 0.002 at N = 1000): most likely
  when N is small, never certain. Under V0 it wins always.

## 3. The cases that cost the model, named

- **C1**, V1 at β = 1 is R2 and R2 is a rate proportional to |ψ|²: the square is written into the
  vacuum's distribution, not derived from the race. It is Fermi's golden rule in the form of a
  first-passage race. It buys a dynamics with one record per quantum and the Born share at once,
  which no other dynamics in this program has done without a stop across pixels, and it derives
  nothing.
- **C2**, the distribution is the physics. Nothing in the model selects the Gumbel form or β = 1;
  a Gaussian or bounded vacuum term, the natural first guesses, gives the brighter pixel every
  time. A model that adopts V1 has to say why the vacuum acts on the deposit as an exponential
  clock, and why with that scale.
- **C3**, V2 shows that the vacuum has to multiply the intensity, not add to the amplitude. The
  interference term between particle and vacuum is linear in |ψ_p|, so a vacuum that enters by
  interference selects in proportion to |ψ|, the square root of Born.
- **C4**, above: the geometry of the vacuum is not modeled; exclusivity is by construction.

## 4. Verification

V1 without the phase weight: P(A) within 3 s.e. of r/(r + 1) at every N and both r. V1 with the
phase weight against the second study's R2 row (0.704, 0.794, 0.800, 0.800) within 3 s.e., taking
that row's own s.e. of 0.0014 into account (combined 0.002). V0 at N = 10: P(A) ≥ 0.999. V2
without the phase weight at N = 1: P(A) within 3 s.e. of √r/(√r + 1). V1 at β = 2, N = 100,
without the phase weight: within 3 s.e. of r²/(r² + 1) = 0.941.

## 5. What will be recorded

`results/share.csv` (every form, N, r, with and without the phase weight), `results/beta.csv`,
`results/scaling.csv`, `results/best_atom.csv`, `results/verification.txt`, `results/run_info.txt`,
`results/fig1_vacuum_race.png`; `RESULTS.md` with the uncertainty first; `LEDGER.md`.

---

## Addendum, 2026-09-14 (later): frozen against fresh vacuum draws

Written before `frozen_fresh.py` was run, at JB's instruction, after his proposal that the
vacuum acts on the fermion first, so that the atom is already dressed when the free particle
arrives. The main run drew the vacuum term fresh at every event. This addendum asks what happens
when the atom's vacuum term is a fixed property of the atom (frozen for the life of the detector),
or refreshes on a correlation time of τ events, and what a single detector's counting statistics
would show in each case.

- **D8, the detector realization.** N = 10 atoms per pixel (primary; N = 100 as sensitivity),
  r = 4, V1 with β = 1. A realization is one detector: its atoms, their phase weights, and their
  vacuum terms. M = 1000 events per realization, K = 300 realizations (200 at N = 100).
- **D9, the vacuum's correlation time.** The Gumbel term of each atom is redrawn every τ events:
  τ = 1 (fresh, the main run), 10, 100, and frozen (drawn once per realization, never redrawn).
- **D10, the phase weight.** Fresh per event (primary, as the main run); frozen per realization
  (sensitivity: the bulk's second-spinor phases as a fixed property of the detector); and absent,
  w = 1 (the clean case).
- **D11, what is recorded.** For each realization the share s_k of events at pixel A. Across
  realizations: the mean, the standard deviation, the Fano factor F = Var(s_k)/[p(1 − p)/M] with p
  the mean (F = 1 is the binomial split of independent events, what a Bernoulli trial per photon
  gives), and the step fraction, the fraction of realizations with s_k < 0.05 or > 0.95.
- **D12, the physical reading.** A single detector's split of M photons between two pixels is
  binomial with F = 1 in every SPAD experiment at every intensity ratio. The model must give F ≈ 1
  for a single detector, so the vacuum term must decorrelate between arrivals. At single-photon
  counting rates the arrivals are about a microsecond apart; any vacuum timescale in silicon
  (10 to 100 fs) is eight orders of magnitude shorter, so τ ≪ 1 event is the physical case, and
  "the atom is dressed before the particle arrives" has to mean dressed afresh for each arrival.
  Frozen terms are the signature of static disorder (the SPAD study's fixed resistances), not of a
  vacuum.

**Predictions (from scratch).**

- **P8, fresh (τ = 1).** Mean 0.795 (R2 at N = 10), F = 1.0 ± 0.1, step fraction 0. Without the
  phase weight: mean 0.800, F = 1.
- **P9, frozen, no phase weight.** Each realization's share is exactly 0 or 1 (step fraction
  1.000, F = M = 1000), and the mean across realizations is 0.800 to within 3 s.e. (0.023): the
  difference of the two pixels' largest Gumbel terms is logistic, so the ensemble of detectors
  still gives Born while every single detector gives all or nothing.
- **P10, frozen, fresh phase weight.** The frozen term acts as a fixed gain per atom; the fresh
  phase weight smooths the step. Mean 0.79, standard deviation about 0.28, F about 470, step
  fraction about 0.49. At N = 100 the same (0.80, 0.28, 470, 0.48): the frozen vacuum does not
  average out with N, because the maximum picks one atom per pixel whatever N.
- **P11, the correlation time sets the Fano factor.** Without the phase weight F = τ (blocks of τ
  identical outcomes): 9 ± 2 at τ = 10, 100 ± 20 at τ = 100. With the fresh phase weight about
  τ/2: 5 ± 1 and 48 ± 10. Step fraction 0 at τ = 10, 0.02 to 0.11 at τ = 100.
- **P12, frozen phase weight, fresh vacuum.** The detector's fixed phase sums bias each
  realization by W_A/W_B: mean 0.794, F = 18 ± 4 at N = 10, 2.5 ± 0.7 at N = 100, 1.2 ± 0.3 at
  N = 1000. A frozen bulk phase distribution is a second source of excess variance that a
  single detector would show at small N and lose at large N; the frozen vacuum's excess (P10)
  does not go away with N.

**The cases that cost the model.**

- **C5**, a frozen vacuum term gives every single detector an all-or-nothing split and Born only
  across an ensemble of detectors (P9, P10). Every SPAD pair shows the Born split within one
  detector with binomial statistics (D12). So the vacuum term of the model must be redrawn between
  arrivals, which is a fluctuation, not a dressing fixed in advance; the physical timescales make
  this automatic (τ about 10⁻⁸ events) and leave no room for "dressed first, then the particle"
  to mean a standing property of the atom.
- **C6**, the Fano factor of a single detector's split is the model's own observable for the
  vacuum's correlation time in units of the arrival interval (P11). It is F = 1 in experiment at
  every counting rate used, which bounds τ below one arrival interval at the highest rates
  (about 10⁸ per second in fast SPADs), i.e. the vacuum term decorrelates in under 10 ns, which
  is consistent with, and far weaker than, any vacuum timescale.

**Verification.** τ = 1 without the phase weight: mean within 3 s.e. of 0.800, F within
[0.7, 1.3]. Frozen without the phase weight: step fraction 1.000 exactly, F = 1000 within 1 %,
mean within 3 s.e. of 0.800 (s.e. 0.4/√300 = 0.023). τ = 10 without the phase weight: F within
[7, 13].

**Recorded.** `results/frozen_fresh.csv`, `results/frozen_fresh_shares.csv` (every realization's
share), `results/verification_frozen.txt`, `results/fig2_frozen_fresh.png`; a section in
`RESULTS.md`; ledger L26.

---

## Correction note, 2026-09-14 (later still)

Written after both runs, on reading the Dirac–Kuramoto record's `NEGATIVE_RESULT.md` and Paper 1
§2 against this file. C1 above says the race gives one record per quantum "without a stop across
pixels". That clause is withdrawn; the text is left as written so the error stays on record. The
first-passage race has one winner because every other clock is cancelled when the first fires, the
losing pixel's included, and that cancellation is the one-quantum stop, the record's premise P5,
which independent site clocks do not supply (they give coincidences at a balanced split) and which
must act faster than light across separated arms. D5's "one winner per trial by construction" was
already the same admission; C1's clause contradicted it. No definition, prediction or number
changes. Ledger L27.

---

## Addendum 2, 2026-09-14 (night): staggered arrival, the survival hazard against the conditional-state hazard

Written before `staggered_arrival.py` was run, after the anticipated-findings ledger for the
Paper 1 v0.10 round (Dirac–Kuramoto record, `born_selection_AI_review_2026-09/`) predicted that
the fresh-vacuum race inherits the failure that record found for its selection postulate in
September: a hazard that runs while the wave is present gives the click to the nearer detector
whenever the arms are unequal. JB's instruction: the failure can be removed. This addendum tests
the two readings of the race's rate under staggered arrival, with the fix the record adopted (the
hazard acts on the conditional state, in the quantum-jump unravelling) as the second reading.

- **D13, the arrival.** The packet occupies pixel A during [0, T] and pixel B during [Δ, Δ + T],
  with Δ/T ∈ {0, 0.1, 0.5, 1, 2}. Two pixels, N = 10 atoms each, intensities I_A = r, I_B = 1,
  r = 4 (shares s_A = 0.8, s_B = 0.2) and r = 2 (sensitivity). Phase weights w_i = cos²(Δφ_i/2)
  fresh per event (primary) and absent, w = 1 (the clean case). W_P = Σ_i w_i over pixel P.
- **D14, reading S, the survival hazard.** While the packet is at pixel P and nothing has fired
  anywhere, atom i of P fires at rate λ I_P w_i. The first atom to fire anywhere takes the quantum
  and every other clock stops (P5). This is the race of the main run, extended in time; at Δ = 0
  it is that run. The coupling is set by κ = λ N T, the integrated per-atom hazard scale over one
  passage at unit intensity: κ ∈ {1, 10, 100}.
- **D15, reading C, the conditional-state hazard.** Atom i of P fires at rate λ w_i a_P²(t), where
  a_P² is the unnormalized amplitude-squared of the wave's P-branch: a_P²(t_P) = s_P at the packet's
  arrival at P, and da_P²/dt = −λ W_P a_P² while the packet is at P, the branch being consumed by
  the port's own no-fire absorption. The integrated hazard of port P is then s_P (1 − e^{−λ W_P T}):
  the port's share times its efficiency, whatever the other port does and whenever it is reached.
  The first fire anywhere still stops every clock (P5). This is the quantum-jump unravelling's
  structure for a single-photon packet at two absorbers, imported as a premise, and it is reading
  B's bookkeeping: nothing is deposited unless a fire occurs. Sampling is exact in both readings
  (exponential first-passage within each window; no time stepping).
- **D16, recorded.** P(A fires), P(B fires), P(no fire); the share among fired events
  P(A)/(P(A) + P(B)), against Δ/T and κ; T = 2 × 10⁵ trials per cell (s.e. ≤ 0.001 on a share).

**Predictions (closed forms; the script computes them beside the samples).**

- **P13, reading S gives the nearer pixel.** Without the phase weight, with λ_A = κr/T and
  λ_B = κ/T: at Δ = 0 the share is r/(r + 1) = 0.800 for every κ (competing exponentials over a
  common window). At Δ ≥ T, P(A) = 1 − e^{−κr}: the share is 1.000 at κ = 10 and 100, and 0.988 at
  κ = 1. For 0 < Δ < T, P(A) = (1 − e^{−λ_A Δ}) + e^{−λ_A Δ} [λ_A/(λ_A + λ_B)] (1 − e^{−(λ_A+λ_B)(T−Δ)}).
  At κ = 10 the share is 0.800, 0.996, 1.000, 1.000, 1.000 at Δ/T = 0, 0.1, 0.5, 1, 2; at κ = 1 it is
  0.800, 0.866, 0.971, 0.988, 0.988. A tenth of a packet length of stagger at physical coupling
  (κ ≥ 10) gives the nearer pixel every record. With the phase weight the same with κ W_P/N in
  place of κ, so slightly slower. This is the September round's finding N1, reproduced in the race.
- **P14, reading C removes it.** The share is s_A = 0.800 (r = 4) and 0.667 (r = 2) at every Δ/T
  and every κ, within 3 s.e.; the efficiency P(fired) is 1 − e^{−κ} without the phase weight
  (0.632, 1.000, 1.000 at κ = 1, 10, 100) and about 1 − e^{−κ/2} with it (0.39, 0.99, 1.00). At
  κ = 1 with the phase weight the share may deviate from 0.800 by up to 0.01, because the two
  ports' efficiencies differ trial by trial.
- **P15, what it costs.** Reading C's integrated hazard equals the port's share by construction:
  the amplitude-squared, put in as the initial condition of the no-fire decay, is the Born weight
  once more, now as the *total* probability that a port ever fires. The race decides which atom
  within a port and, through P5, that only one port records; it does not decide the share between
  ports under staggered arrival, and could not, since the ports are never in the race at the same
  time. So the fresh-vacuum reading's premise 2 (MODEL_SECTION §12.1) has to read: the rate is
  proportional to the conditional-state deposit, and a port's integrated hazard is its share. The
  Born share is then inherited twice over, and the ledger entry (L28) says so.

**Verification.** Reading S, no phase weight, κ = 10: share at Δ = 0 within 3 s.e. of 0.800; P(A) at
Δ = T ≥ 0.999; every sampled share within 3 s.e. of its closed form at all five Δ/T and three κ.
Reading C, no phase weight, κ = 10: share within 3 s.e. of 0.800 at all five Δ/T; P(fired) ≥ 0.999.
Reading C, κ = 1, no phase weight: P(fired) within 3 s.e. of 1 − e⁻¹ = 0.632.

**Recorded.** `results/staggered.csv`, `results/verification_staggered.txt`,
`results/fig3_staggered.png`; a section in `RESULTS.md`; ledger L28.

---

## Addendum 2b, 2026-09-14 (night, after addendum 2 was run): the normalized conditional state

Written after `staggered_arrival.py` and before `staggered_conditional.py`. **P14 failed as
registered** (`results/verification_staggered.txt`): reading C gave a share of 0.857 at Δ = 0 and
0.952 from Δ = T/2 at κ ≥ 10, against the registered 0.800, with P(fired) = 0.84 against 1.000.
The cause is in D15, not in the code. D15 let each port's branch decay on its own and coupled the
ports only by the stop on firing. Then P(B fires) = (1 − s_A η_A) s_B η_B: the second port fires
with its share times the probability that the first did not, and the share among records is
s_A η_A / (s_A η_A + (1 − s_A η_A) s_B η_B) = 0.8/(0.8 + 0.04) = 0.952 at η = 1, which is the number
the run gave. D15's sentence "whatever the other port does" was wrong: it holds for the
unconditional probabilities only if a null at one port raises the rate at the other. That is
what the quantum-jump unravelling does through the normalization of the conditional state, and
it is the second face of P5: the stop cancels the losing clocks when one fires, and the null
renormalizes the survivors when none does. A reading that is local port by port does not remove
the staggered-arrival failure; it replaces one wrong share (the nearer port's 1.000) with another
(0.952).

- **D17, reading C′.** As D15, with the rates taken on the *normalized* conditional state: at time
  t the unnormalized branch weights are a_P²(t) = s_P e^{−λ W_P τ_P(t)}, τ_P the exposure so far at
  P, and n(t) = a_A²(t) + a_B²(t) the survival norm. The probability that the first fire anywhere
  occurs in [t, t + dt] given none before is −dn/n, and the port is P with probability
  proportional to λ W_P a_P²(t) [P exposed at t]. Sampled exactly by inverting n(t) on a grid of
  4000 times (n is monotone) and drawing the port at the crossing.
- **D18, everything else as D13, D14, D16.** The same κ, Δ/T, r, phase settings, T = 2 × 10⁵
  trials. Reading S is not rerun. One labelled diagnostic, not a registered check: the reading-S
  cell (κ = 100, Δ = 0, no phase weight), which failed its closed form by 3.0 s.e. in addendum 2
  (0.7973 against 0.8000), re-drawn once with a fresh seed.

**Predictions.**

- **P16, reading C′ gives the share at every stagger.** Without the phase weight: share = s_A
  = 0.800 (r = 4) and 0.667 (r = 2) at all five Δ/T and all three κ, within 3 s.e.; P(fired) =
  1 − e^{−κ} = 0.632, 1.000, 1.000. With the phase weight: the same shares within 0.01 at κ ≥ 10
  (efficiencies about 1 − e^{−κ/2} = 0.993, 1.000); at κ = 1 the share is s_A η_A/(s_A η_A + s_B η_B)
  averaged over the ports' unequal efficiencies, predicted 0.80 ± 0.02 with P(fired) about 0.39.
- **P17, the closed forms of reading C, for the record.** P(A) = s_A η_A; P(B) = (1 − s_A η_A) s_B η_B
  for Δ ≥ T; the run's 0.952 and 0.841 are these at η = 1 − e^{−10}. Not rerun.
- **P18, the diagnostic.** The κ = 100, Δ = 0 reading-S cell lands within 3 s.e. of 0.800 on a
  fresh seed; if it does not, the sampling of very short exponential times is suspect and is to
  be examined.

**What it means, stated before the run.** If P16 holds, the staggered-arrival failure is removed
from the fresh-vacuum race at a stated price: the race's rates are those of the normalized
conditional state, so P5 acts twice, cancelling on a fire and renormalizing on a null, and the
share between ports under staggered arrival is not decided by any race but by the branch weights
put in as initial conditions. MODEL_SECTION §12.1 premise 2 and the Paper 1 v0.10 Addendum A.7 are
to be amended accordingly, and the review ledger's item 3 marked removed before the freeze.

**Verification.** Reading C′, no phase weight, κ = 10: share within 3 s.e. of 0.800 at all five
Δ/T and P(fired) ≥ 0.999; κ = 1: P(fired) within 3 s.e. of 0.632; r = 2, κ = 10, Δ = T: share
within 3 s.e. of 0.667.

**Recorded.** `results/staggered_conditional.csv`, `results/verification_staggered_conditional.txt`,
`results/fig3_staggered.png` (regenerated with reading C′ replacing reading C; the addendum-2
figure is kept as `fig3_staggered_local.png`); a section in `RESULTS.md`; ledger L28.
