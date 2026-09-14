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
