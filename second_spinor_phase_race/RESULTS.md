# Results: the second-spinor phase race at a detector

**Date:** 2026-09-07. Definitions and predictions were fixed in `PREREGISTRATION.md` before the
run and are not restated here; read that file first, including its dated addendum (§6). Run:
`python3 phase_race.py`, `python3 bohm_dirac_flight.py`, `python3 interference_visibility.py`,
`python3 energy_cost.py`, `python3 sink_flow.py`, `python3 sink_race.py`, then
`python3 splitter_asymmetry.py` and `python3 splitter_asymmetry.py slow`; every number below is in `results/`.
The dated side ledger of departures from quantum mechanics or experiment is `LEDGER.md`.

## Read this before quoting any number

- The model is JB's, built as he described it (pre-registration §0). Nothing here is a
  measurement, and nothing here is a result of the Dirac–Kuramoto record. What the model departs
  from is in the ledger and is not argued here.
- Parts A, B and the ordered bulk are exact properties of the definitions: closed forms checked by
  Monte Carlo (worst difference 1.9 × 10⁻³ in a probability and 1.0 × 10⁻³ in a mean deposit,
  against up to 4 × 10⁵ trials). Part C is a 1+1-D free Dirac packet propagated exactly in
  momentum space with RK4 trajectories against the exact field; 4000 arrivals match ρ(x, t_d)
  to 0.013 in the cumulative distribution.
- Every number depends on the choices A1 to A5 of the pre-registration, above all the cos² form
  of "closest in phase" (A1) and the three readings of "based on the deposited energy" (A2). The
  band form is shown beside the cos² form where it changes the answer.
- The record rule is R2, sampling by deposit (JB, 2026-09-07); R3 stays in the tables as the
  deterministic limit. The interference reading of the weight is Part D at the end.
- One registered prediction missed: the zitterbewegung frequency 2 (in units of mc²/ħ) holds for a
  particle at rest, and a localized packet beats faster. It is recorded under "Predictions
  against outcome", with the addendum's wider packet.

## The answer in one paragraph

With a random bulk the phase race decides *which* atom fires and *how much* energy it takes, and
nothing about *where* the record lands or about the incident particle's own second-spinor phase.
Among N atoms in reach there is one within 0.45 rad of the particle's second-spinor phase (energy
weight ≥ 0.95) with 95 % confidence once N ≥ 20, and within 0.20 rad (weight ≥ 0.99) once
N ≥ 46. Under the primary rule R2 the atom that fires is typically 54° off in phase whatever N is;
under R3 it is the closest one, π/(N + 1). Across two lobes of unequal intensity, R2 returns the
intensity ratio (Born's r/(1 + r), exponent 2) once each lobe holds more than about ten atoms,
and R3 returns the brighter lobe every time. Only an ordered bulk reads the particle's phase: the
deposit and the firing probability then follow ½[1 + R_κ cos(δ_p − μ)]. In free flight the
internal offset between the positive- and negative-energy components is a constant of the motion
to machine precision, and the Bohm–Dirac trajectory carries it as the phase of its
zitterbewegung at the instant of contact.

## Table 1. Part A: uniform wave, random bulk, N atoms in reach (Monte Carlo; closed forms in the CSV)

Source: `results/part_a_uniform.csv`. Band ε = 0.35 rad. A random atom is 90° off on average.

| N atoms | best mismatch | P(best takes ≥ 95 %) | P(best takes ≥ 99 %) | R1: best atom's share | R2: fired atom's mismatch | band: some atom in reach |
|---|---|---|---|---|---|---|
| 1 | 89.9° | 0.144 | 0.064 | 1.0000 | 89.9° | 0.112 |
| 2 | 60.0° | 0.267 | 0.124 | 0.7500 | 69.7° | 0.210 |
| 3 | 45.0° | 0.371 | 0.179 | 0.5817 | 62.7° | 0.298 |
| 5 | 30.0° | 0.539 | 0.280 | 0.3878 | 58.2° | 0.446 |
| 10 | 16.3° | 0.788 | 0.484 | 0.2024 | 55.5° | 0.694 |
| 20 | 8.6° | 0.955 | 0.732 | 0.1015 | 54.4° | 0.906 |
| 26 | 6.7° | 0.982 | 0.820 | 0.0779 | 54.2° | 0.954 |
| 46 | 3.8° | 0.999 | 0.952 | 0.0439 | 53.8° | 0.995 |
| 100 | 1.8° | 1.000 | 0.998 | 0.0201 | 53.8° | 1.000 |
| 1000 | 0.2° | 1.000 | 1.000 | 0.0020 | 53.9° | 1.000 |
| 10000 | 0.0° | 1.000 | 1.000 | 0.0002 | 53.9° | 1.000 |

The closest atom's mismatch falls as π/(N + 1); the R2-fired atom's mismatch stays at the closed
form π/2 − 2/π = 53.5°; the R1 share falls as 2/N.

## Table 2. Part B: two lobes, intensity ratio r = 4 (Born share 0.800)

Source: `results/part_b_two_lobes.csv`. γ from P₁/P₂ = (√r)^γ; Born gives 2.

| n atoms per lobe | 1 | 2 | 5 | 10 | 100 | 1000 |
|---|---|---|---|---|---|---|
| R2: P(lobe 1) | 0.704 | 0.759 | 0.788 | 0.794 | 0.800 | 0.800 |
| R2: exponent γ | 1.25 | 1.66 | 1.89 | 1.95 | 2.00 | 2.00 |
| R3: P(lobe 1) | 0.790 | 0.924 | 0.997 | 1.000 | 1.000 | 1.000 |
| R3: exponent γ | 1.92 | 3.61 | 8.21 | 16 | ∞ | ∞ |
| R2, band form: P(lobe 1) | 0.521 | 0.533 | 0.581 | 0.651 | 0.791 | 0.799 |
| band form: miss rate | 0.791 | 0.624 | 0.304 | 0.094 | 0.000 | 0.000 |

Exponent against the amplitude ratio with one atom per lobe (n = 1):

| amplitude ratio √r | 1.41 | 2 | 3 | 4 |
|---|---|---|---|---|
| Born P(lobe 1) | 0.667 | 0.800 | 0.900 | 0.941 |
| R2: P(lobe 1), γ | 0.609, 1.27 | 0.704, 1.25 | 0.795, 1.23 | 0.843, 1.21 |
| R3: P(lobe 1), γ | 0.695, 2.37 | 0.790, 1.92 | 0.863, 1.67 | 0.897, 1.56 |

R2 tends to Born's 2 from below as the aperture fills, because the phase factor averages to ½ in
both lobes; R3 tends to a step. The record's own race (Paper 1) gave 1.56, interval 1.44 to 1.69.

## Table 3. Ordered bulk: the deposit against the particle's own phase

Source: `results/ordered_bulk.csv` (κ = 0, 1, 4 at 15° steps). Two halves of 100 atoms, μ = 0
beside μ = π. Closed form ½[1 + R_κ cos δ_p] for the deposit and for the R2 two-half probability.

| κ | R_κ | δ_p | mean deposit per atom (half μ = 0) | closed form | two halves, R2: P(μ = 0 half) | two halves, R3 |
|---|---|---|---|---|---|---|
| 1 | 0.446 | 0° | 0.723 | 0.723 | 0.724 | 0.884 |
| 1 | 0.446 | 45° | 0.658 | 0.658 | 0.658 | 0.806 |
| 1 | 0.446 | 90° | 0.500 | 0.500 | 0.500 | 0.499 |
| 1 | 0.446 | 135° | 0.342 | 0.342 | 0.342 | 0.194 |
| 1 | 0.446 | 180° | 0.277 | 0.277 | 0.277 | 0.116 |
| 4 | 0.864 | 0° | 0.932 | 0.932 | 0.932 | 1.000 |
| 4 | 0.864 | 45° | 0.805 | 0.805 | 0.805 | 0.996 |
| 4 | 0.864 | 90° | 0.500 | 0.500 | 0.500 | 0.501 |
| 4 | 0.864 | 135° | 0.195 | 0.195 | 0.194 | 0.004 |
| 4 | 0.864 | 180° | 0.068 | 0.068 | 0.068 | 0.000 |

At κ = 0 (the random bulk) every entry is 0.500. The two-half R2 curve is the deposit curve itself;
R3 steepens it toward a step of width about 1/κ.

## Part C: free flight by the Bohm–Dirac velocity

Source: `results/flight_arrivals.csv`, `results/flight_zitter.csv`, `fig4_flight.png`,
`fig5_zitter.png`. Two lobes r = 4 at x = ±4, σ = 1, t_d = 4, 4000 trajectories, 21 atoms in reach.

- norm at t_d: |∫ρ − 1| = 2.2e-16;  arrivals at x > 0: 0.8050, ∫_(x>0) ρ(t_d) = 0.7922 (s.e. 0.0064), |diff| = 0.0128   -> PASS
- max |CDF(arrivals) − CDF(ρ(t_d))| = 0.0133   (tolerance 0.03)   -> PASS
- r/(1+r) = 0.800; the lobes' spreading across x = 0 by t_d moves ∫_(x>0) ρ(t_d) to 0.7922
- winner's energy weight, mean: R2 0.7360 (predicted 0.75), R3 0.9905 (predicted ≈ 0.99); blur x_record − x_arrival: R2 mean +0.0019, sd 0.3037; R3 mean +0.0033, sd 0.3026 (uniform on ±0.5: sd 0.2887) The discrete uniform on the 21 lattice sites has sd 0.303, which is what both rules give: the race adds a position blur of ±L_c and nothing else.

The internal offset in free flight (D9; θ = 60°):

| case | σ | δ | v at t = 0 | sin θ cos δ | ±E offset over the packet: mean, spread | offset drift to t = 12 | beat, early / late | 2⟨E⟩ | beat amplitude | position half-range |
|---|---|---|---|---|---|---|---|---|---|---|
| registered | 1 | θ = 0 | max \|v\| 0.741 | 0 | none | | none | | | 2.14 (drift only) |
| registered | 1 | 0° | 0.8660 | 0.8660 | 0.0°, 3.0° | 5e-12° | 2.42 / 2.36 | 2.22 | 0.820 (sin θ = 0.866) | 1.09 (0.43 if no drift) |
| registered | 1 | 90° | 0.0000 | 0.0000 | 90.0°, 12.7° | 2e-12° | 2.30 / 2.01 | 2.22 | 0.858 (sin θ = 0.866) | 0.43 (0.43 if no drift) |
| addendum 2026-09-07 | 3 | 0° | 0.8660 | 0.8660 | -0.0°, 0.0° | 6e-13° | 2.11 / 2.07 | 2.03 | 0.850 (sin θ = 0.866) | 0.51 (0.43 if no drift) |
| addendum 2026-09-07 | 3 | 90° | 0.0000 | 0.0000 | 90.0°, 5.3° | 4e-13° | 2.06 / 2.05 | 2.03 | 0.868 (sin θ = 0.866) | 0.42 (0.43 if no drift) |

Reading: at the instant of contact the trajectory's velocity is sin θ cos δ exactly, so the
internal offset is the phase of the zitterbewegung. The ±E relative phase at each momentum, minus
2E(p)t, does not move at all over twelve Compton times (drift below 10⁻¹¹ degrees): the free
equation maintains the offset. What does move is the beat in position space, because the packet's
momentum components beat at different 2E(p): the fitted frequency sits above 2⟨E⟩ early and
falls toward 2 as the faster components leave the centre, and the offset itself varies with p
across a narrow packet (spread 3° to 13° at σ = 1, 0° to 5° at σ = 3), since (1, 0) and (0, 1) are
the ±E spinors only at p = 0. The δ = 0 trajectory starts at its maximum velocity and drifts with
the spreading packet, so its position range exceeds the stationary prediction 0.43; the velocity
oscillation itself has the predicted amplitude.

## Verification (`results/verification.txt`)

| Check | Result |
|---|---|
| Part A: Monte Carlo against closed forms (mean best mismatch, P(w_max ≥ 0.95), P(w_max ≥ 0.99), band hit), all N | worst difference 1.9 × 10⁻³, pass |
| Part B: R2 at r = 4, n = 1000 against Born 0.800 | 0.8000, pass |
| Part B: pre-registered r = 4 table | reproduced to 0.001 in P₁ and 0.01 in γ |
| Ordered bulk: mean deposit against ½[1 + R_κ cos δ_p] | worst difference 1.0 × 10⁻³, pass |
| Part C: norm, arrivals against ∫ρ, cumulative distribution | 2 × 10⁻¹⁶; within 2 s.e.; 0.013, pass |
| Part C: v(0, 0) = sin θ cos δ; ±E offset maintained; beat within 10 % of 2⟨E⟩ | exact; drift < 10⁻¹¹ °; pass |

## The consequences named in §3, resolved

- **C1, R1 spreads the quantum thin.** The best atom's share is 0.75 at N = 2, 0.20 at N = 10,
  0.020 at N = 100 and 2 × 10⁻⁴ at N = 10⁴. R1 alone puts no quantum at one site.
- **C2, under R2 a random bulk hides the particle's phase.** Confirmed exactly: every Part A
  statistic is independent of δ_p (Δφ_i is uniform), the fired atom is 54° off on average for
  every N, and the two-lobe probability is the intensity share once n ≳ 10 (Table 2).
- **C3, under R3 the brighter lobe always fires once the aperture is full.** P₁ = 0.997 at
  n = 5 and 1.000 from n = 10. At n = 1 the exponent is 1.92 for r = 4 but 2.37 at r = 2 and
  1.56 at r = 16: it is not a constant, so it is not Born's law at any n.
- **C4, the band form has a miss rate** 0.79 at one atom per lobe and 0.09 at ten, 0 from a
  hundred; with the aperture full its two-lobe answer is the same as the cos² form's.
- **C5, the distinctive prediction.** In an ordered bulk the deposit and the R2 firing probability
  follow ½[1 + R_κ cos(δ_p − μ)] to 10⁻³ (Table 3): visibility 0.45 at κ = 1, 0.86 at κ = 4.
  This is the one place the model's phase mechanism leaves a trace, and where it would be tested.
- **C6, in the trajectory reading the site statistics are Bohm's.** Arrivals follow ρ(x, t_d);
  the race moves the record by a uniform ±0.5 and changes nothing else.

## Predictions against outcome

- Parts A, B and the ordered bulk: every closed form and the scratch Monte Carlo table of §2 were
  reproduced (Part B to 0.001 in P₁).
- Part C (i), equivariance: held (0.805 against 0.792, two standard errors). (ii), the race as a
  uniform blur with winner weights 0.75 and 0.99: held (0.736 with 21 atoms, 0.990).
- Part C (iii), the zitterbewegung: the phase at t = 0 equals δ exactly and the offset is
  maintained exactly, as predicted. **The frequency prediction 2 missed**: the registered
  packet beats at 2.30 to 2.42 early, against 2⟨E⟩ = 2.22; the addendum's σ = 3 packet at 2.06
  to 2.11 against 2.03. The prediction was written for a particle at rest and a packet is not at
  rest. The velocity amplitude held (0.82 to 0.87 against 0.866); the position range for δ = 0
  exceeded 0.43 because the trajectory drifts with the spreading packet.

## Part D (added 2026-09-07): the interference term, with the record rule fixed at R2

Source: `results/part_d_visibility.csv`, `part_d_visibility_ordered.csv`, `part_d_bound_amplitude.csv`,
`fig6_interference.png`. Registered in the pre-registration's second addendum (§6) before the run.
JB's clarification: the deposit is highest at the smallest mismatch because the incident particle's
second spinor interferes with the bound electron's; sampling by deposit is the chosen rule. The
cos² weight of A1 is that interference term at full visibility, so the tables above already
assume it. What the clarification adds is the size of the bound electron's amplitude relative to
the incident one, in two readings.

**D10, the interference gates the particle's energy** (visibility V; V = 1 is the study as run).
The place of the record does not move; every phase effect scales with V.

| D10, visibility V | 1 | 0.5 | 0.2 | 0 |
|---|---|---|---|---|
| fired atom's mean mismatch (closed) | 53.8° (53.5°) | 72.0° (71.8°) | 82.8° (82.7°) | 90.1° (90.0°) |
| fired atom within the 0.35 rad band (closed) | 0.219 (0.221) | 0.165 (0.166) | 0.134 (0.133) | 0.112 (0.111) |
| two lobes r = 4, P₁ at n = 1 | 0.706 | 0.787 | 0.798 | 0.800 |
| two lobes r = 4, P₁ at n = 2 | 0.759 | 0.793 | 0.799 | 0.800 |
| two lobes r = 4, P₁ at n = 5 | 0.788 | 0.797 | 0.800 | 0.800 |
| two lobes r = 4, P₁ at n = 10 | 0.794 | 0.799 | 0.800 | 0.800 |
| two lobes r = 4, P₁ at n = 100 | 0.799 | 0.800 | 0.800 | 0.800 |
| two lobes r = 4, P₁ at n = 1000 | 0.800 | 0.800 | 0.800 | 0.800 |
| ordered bulk κ = 4, two halves, P₁ at δ_p = 0 (closed) | 0.932 (0.932) | 0.716 (0.716) | 0.586 (0.586) | 0.500 (0.500) |

**D11, the bound electron's own energy joins the deposit** (β = bound amplitude in units of the
incident amplitude in the weaker lobe). The intensity share is diluted toward ½ as β grows, and
the dilution is already complete at one atom per lobe.

| D11, bound amplitude β | 0 | 0.5 | 1 | 2 |
|---|---|---|---|---|
| two lobes r = 4, P₁ at n = 1 | 0.800 | 0.771 | 0.698 | 0.543 |
| two lobes r = 4, P₁ at n = 2 | 0.800 | 0.772 | 0.704 | 0.585 |
| two lobes r = 4, P₁ at n = 5 | 0.800 | 0.772 | 0.710 | 0.605 |
| two lobes r = 4, P₁ at n = 10 | 0.800 | 0.772 | 0.713 | 0.611 |
| two lobes r = 4, P₁ at n = 100 | 0.800 | 0.773 | 0.714 | 0.615 |
| two lobes r = 4, P₁ at n = 1000 | 0.800 | 0.773 | 0.714 | 0.616 |
| closed form as n → ∞, (r + β²)/(r + 1 + 2β²) | 0.800 | 0.773 | 0.714 | 0.615 |

Verification: worst |MC − closed| over the fired atom's mean mismatch, its band fraction, P₁ at n = 1000 against Born, and the ordered two-half law: 5.29e-03   (tolerance 0.01)   -> PASS worst |P₁(n = 1000) − (r + β²)/(r + 1 + 2β²)| = 1.60e-04   (tolerance 0.01)   -> PASS

**Does the interference term change the result?** Under D10, no: Born's share is reached at every
visibility, sooner for lower V, and the random bulk still hides the particle's phase; the one
thing that changes is the visibility of the ordered-bulk law, V R_κ. Under D11, yes, in one
direction: the bound electron's self-term adds a phase-independent deposit to every atom the
wave reaches, so the brighter lobe's share falls from 0.800 toward ½ as the bound amplitude grows
past the incident one. Which reading JB means is the definition to fix next, together with the
size of β. The consequence stated in the addendum for atoms the wave does not reach, a record
outside the wave, follows from D11 alone and was not run.

## Part E (added 2026-09-07): a small energy cost of the phase change, shared equally

Source: `results/part_e_cost.csv`, `part_e_cost_aligned.csv`, `fig7_cost.png`. Registered in the
pre-registration's third addendum (§6) before the run. JB's addition: the shift of the particle's
spinors at contact costs a small energy, shared equally by the two particles. As fixed there, the
particle pays half of c(Δ) = η(1 − cos Δ)/2 out of what it can deposit, with Δ the first-spinor
shift, and an atom whose deposit cannot pay does not fire.

**Random first spinors in the bulk (primary).**

| η (cost of a half-turn, in units of the weaker lobe's energy) | 0 | 0.05 | 0.1 | 0.2 | 0.5 |
|---|---|---|---|---|---|
| two lobes r = 4, P₁ at n = 1 | 0.704 | 0.706 | 0.708 | 0.712 | 0.721 |
| two lobes r = 4, P₁ at n = 2 | 0.759 | 0.761 | 0.764 | 0.767 | 0.780 |
| two lobes r = 4, P₁ at n = 5 | 0.788 | 0.791 | 0.794 | 0.799 | 0.813 |
| two lobes r = 4, P₁ at n = 10 | 0.795 | 0.797 | 0.800 | 0.805 | 0.820 |
| two lobes r = 4, P₁ at n = 100 | 0.799 | 0.802 | 0.805 | 0.810 | 0.825 |
| two lobes r = 4, P₁ at n = 1000 | 0.800 | 0.803 | 0.805 | 0.811 | 0.826 |
| P₁ as n → ∞ (quadrature; no-floor closed form) | 0.800 (0.800) | 0.803 (0.803) | 0.805 (0.806) | 0.811 (0.812) | 0.826 (0.833) |
| exponent γ as n → ∞ (Born 2) | 2.00 | 2.03 | 2.05 | 2.10 | 2.25 |
| miss rate at n = 1 | 0.000 | 0.002 | 0.004 | 0.008 | 0.021 |
| fired atom's mean mismatch, 50 atoms | 53.8° | 53.2° | 52.6° | 51.4° | 48.6° |

**Aligned first spinors in the bulk.** The cost is then one number set by the particle's own
first-spinor phase against the bulk's, and that phase enters the site statistics for the first time.

| aligned first spinors, n = 1000: particle's shift Δ | 0° | 45° | 90° | 135° | 180° |
|---|---|---|---|---|---|
| η = 0.1: P₁ (no-floor closed form) | 0.800 (0.800) | 0.802 (0.802) | 0.806 (0.806) | 0.809 (0.811) | 0.811 (0.812) |
| η = 0.5: P₁ (no-floor closed form) | 0.800 (0.800) | 0.808 (0.809) | 0.827 (0.833) | 0.846 (0.862) | 0.854 (0.875) |

Verification: worst |P₁(n = 1000) − quadrature large-n| over η: 1.15e-04   (tolerance 0.01)   -> PASS worst |P₁(n = 1000) − quadrature large-n| over η and Δ: 1.58e-04   (tolerance 0.01)   -> PASS

**Does the cost change the result much?** No, while it is small. It is a threshold: the atoms whose
deposit cannot pay it are cut, and they are more often in the weaker lobe, so the brighter lobe's
share rises above Born by about η/(r + 1): half a percent at η = 0.05, a percent at η = 0.1,
three percent at η = 0.5, where the exponent reaches 2.25. The fired atom is also a little
closer in phase, because the weakest deposits are the ones cut. With aligned first spinors the
particle's own first-spinor phase enters by the same threshold and the same small amount, and
nothing at all changes when the particle needs no shift. The random-bulk conclusion for the
second-spinor phase (C2) and the form of the ordered-bulk law (C5) are unchanged. A cost taken as
a fixed fraction of the deposit, instead of subtracted from it, would cancel in the ratio and
change nothing.

## Part F (added 2026-09-07): the winning atom becomes a sink

Source: `results/part_f_sink_summary.csv`, `part_f_sink_series.csv`, `fig8_sink.png`. Registered
in the pre-registration's Part F addendum (§6) before the run, with the defaults JB asked for:
sink strength from the deposit, growth in proportion to the energy collected, width one lattice
spacing, each swept one at a time. Time is measured from the deposit, in Compton times; the state
is Part C's packet at the detector time. "Near" and "far" are the norms on the sink's side of
x = 0 and on the other side. Verification: the split-operator scheme reproduces the exact free
propagator to 2 × 10⁻¹⁴, a well alone conserves the norm to 10⁻¹³, the initial collection rate
equals the registered Γ₀ ρ(x_s) √(2π) L₀ to three decimals, and near + far + collected = 1 to
machine precision in every run. The Klein diagnostic stayed at or below 0.07 while the
remainder was above 0.05 of the quantum, so no well entered pair-production territory.

| run | near side's initial share | initial rate | collected at t = 5 | t = 10 | t = 20 | t = 40 | near side half-time | far side at t = 40 (free reference) | sink width reached |
|---|---|---|---|---|---|---|---|---|---|
| base (γ₀ = 1, u₀ = 0.5, L₀ = 0.05, κ = 10) | 0.792 | 0.0130 | 0.124 | 0.178 | 0.287 | 0.382 | 17.0 | 0.362 (0.458) | 0.05 |
| absorber ×0.1 | 0.792 | 0.0013 | 0.011 | 0.015 | 0.024 | 0.047 | never | 0.406 (0.458) | 0.05 |
| absorber ×10 | 0.792 | 0.1301 | 0.263 | 0.277 | 0.306 | 0.322 | 5.5 | 0.414 (0.458) | 0.05 |
| no well | 0.792 | 0.0130 | 0.113 | 0.151 | 0.227 | 0.311 | 19.5 | 0.396 (0.458) | 0.05 |
| well 1.5 mc² from the start | 0.792 | 0.0130 | 0.147 | 0.218 | 0.347 | 0.432 | 14.5 | 0.338 (0.458) | 0.05 |
| width 0.5 | 0.792 | 0.1271 | 0.419 | 0.446 | 0.481 | 0.498 | 2.0 | 0.324 (0.458) | 0.50 |
| width 2 | 0.792 | 0.4099 | 0.862 | 0.878 | 0.887 | 0.890 | 0.5 | 0.106 (0.458) | 2.00 |
| no growth (κ = 0) | 0.792 | 0.0130 | 0.086 | 0.112 | 0.166 | 0.251 | 34.0 | 0.371 (0.458) | 0.05 |
| channel extends, λ = 10 | 0.792 | 0.0130 | 1.000 | 1.000 | 1.000 | 1.000 | 1.0 | 0.000 (0.458) | 10.05 |
| winner in the weak lobe (x_s = −4) | 0.207 | 0.0032 | 0.023 | 0.031 | 0.103 | 0.233 | never | 0.554 (0.542) | 0.05 |
| wide lobes ±10, σ = 3, base | 0.799 | 0.0094 | 0.118 | 0.267 | 0.452 | 0.522 | 15.5 | 0.276 (0.259) | 0.05 |
| wide lobes, channel λ = 10 | 0.799 | 0.0094 | 0.994 | 0.999 | 1.000 | 1.000 | 1.5 | 0.000 (0.259) | 10.05 |

**Reading.** A sink at one atom does not draw the packet in. It eats what arrives, and what
arrives is set by the packet's own spreading, not by any pull: the far side's norm follows the
free run to within the amount that spreads across x = 0 and is then eaten (base 0.36 against
0.46). The base sink has collected 0.38 of the quantum after 40 Compton times and has not yet
halved the near lobe's norm until t = 17. The well helps little (0.31 without it, 0.43 with a
well at the cap from the start), because a well of width 0.05 binds weakly; growth helps more
(0.25 without growth). A stronger absorber is not a faster one: at ten times the rate the sink
collects more in the first five Compton times and less by forty, the reflecting side of an
imaginary potential. A wider sink collects faster in proportion to what it covers: at width 2 it
takes the near lobe in half a Compton time and then eats the far lobe's spreading tail, ending
above the near lobe's share. The channel run collects the whole quantum within five Compton
times, and the way it does so is the finding: as the sink widens with what it has collected, its
absorbing region grows over the far lobe and consumes it in place. In one dimension a channel
from the atom to the far lobe is the interval between them, so "the channel reaches the
reservoir" and "the sink covers the reservoir" are the same thing. Nothing flowed. The wide-lobe
runs say the same with the two lobes kept apart: 0.52 collected by the fixed sink, all of it by
the channel.

**Predictions against outcome.** F1 held exactly. F2 held: collection is sublinear in the
absorber strength and lower by t = 40 at ten times the rate. F3 held for widths 0.05 and 0.5 and
**missed at width 2**, where the sink ends above the near share by eating the far lobe's tail. F4
held for its first half, nothing flows toward a fixed sink, and for its second half in the letter
but not the mechanism: the channel run collected the far lobe, but by covering it, not after
draining the near lobe first. F5 held (0.23 against 0.38). F6 held once the diagnostic was
restricted to a remainder above 0.05; on a vanished remainder the fraction is noise, and that
masking was added after the first run and is recorded here.

## Part G (added 2026-09-07): two sinks racing, with reversible charging, a latch, and vacuum noise

Source: `results/part_g_race_summary.csv`, `part_g_race_trials.csv`, `fig9_race.png`. Registered
in the pre-registration's Part G addendum (§6) before the run, with the defaults stated there:
each site charges reversibly from the norm in its well and loses charge with a ten-Compton-time
relaxation; at a charge of 0.03 the site latches and its absorber switches on; the well, the
absorber and the channel grow as in Part F; an Ornstein–Uhlenbeck noise of strength s on each
charge stands for the vacuum. Time from the deposit, 20 Compton times, 168 trials in all.

| r | site 2 offset | noise s | channel | trials | site 1 first | site 2 first | tie | both latch | loser's reversal ratio | latch t₁ | latch t₂ |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 0 | 0 | no channel | 1 | 0.00 | 0.00 | 1.00 | 1.00 | – | 4.12 | 4.12 |
| 1 | 0 | 0 | channel | 1 | 0.00 | 0.00 | 1.00 | 1.00 | – | 4.12 | 4.12 |
| 1 | 0 | 0.005 | no channel | 10 | 0.70 | 0.30 | 0.00 | 1.00 | 0.882 | 3.69 | 3.89 |
| 1 | 0 | 0.005 | channel | 10 | 0.50 | 0.50 | 0.00 | 1.00 | 0.264 | 3.81 | 3.64 |
| 1 | 0 | 0.02 | no channel | 10 | 0.50 | 0.50 | 0.00 | 1.00 | 0.874 | 2.39 | 2.37 |
| 1 | 0 | 0.02 | channel | 10 | 0.30 | 0.70 | 0.00 | 1.00 | 0.191 | 3.50 | 1.75 |
| 1 | 0.5 | 0 | no channel | 1 | 1.00 | 0.00 | 0.00 | 1.00 | 0.998 | 4.12 | 11.76 |
| 1 | 0.5 | 0 | channel | 1 | 1.00 | 0.00 | 0.00 | 0.00 | 0.280 | 4.12 | never |
| 1 | 0.5 | 0.005 | no channel | 10 | 0.60 | 0.40 | 0.00 | 1.00 | 0.954 | 4.10 | 5.35 |
| 1 | 0.5 | 0.005 | channel | 10 | 0.80 | 0.20 | 0.00 | 1.00 | 0.291 | 4.14 | 4.70 |
| 1 | 0.5 | 0.02 | no channel | 10 | 0.50 | 0.50 | 0.00 | 1.00 | 0.928 | 2.36 | 2.16 |
| 1 | 0.5 | 0.02 | channel | 10 | 0.40 | 0.60 | 0.00 | 1.00 | 0.207 | 2.86 | 2.34 |
| 1 | 1 | 0 | no channel | 1 | 1.00 | 0.00 | 0.00 | 1.00 | 0.985 | 4.12 | 9.85 |
| 1 | 1 | 0 | channel | 1 | 1.00 | 0.00 | 0.00 | 0.00 | 0.295 | 4.12 | never |
| 1 | 1 | 0.005 | no channel | 10 | 0.80 | 0.20 | 0.00 | 1.00 | 0.968 | 3.99 | 5.94 |
| 1 | 1 | 0.005 | channel | 10 | 0.90 | 0.10 | 0.00 | 0.50 | 0.269 | 3.61 | 4.80 |
| 1 | 1 | 0.02 | no channel | 10 | 0.80 | 0.20 | 0.00 | 1.00 | 0.942 | 1.56 | 2.65 |
| 1 | 1 | 0.02 | channel | 10 | 0.70 | 0.30 | 0.00 | 1.00 | 0.206 | 2.09 | 3.26 |
| 4 | 0 | 0 | no channel | 1 | 1.00 | 0.00 | 0.00 | 1.00 | 1.000 | 1.92 | 19.43 |
| 4 | 0 | 0 | channel | 1 | 1.00 | 0.00 | 0.00 | 0.00 | 0.182 | 1.92 | never |
| 4 | 0 | 0.005 | no channel | 10 | 1.00 | 0.00 | 0.00 | 1.00 | 1.000 | 1.84 | 17.41 |
| 4 | 0 | 0.005 | channel | 10 | 1.00 | 0.00 | 0.00 | 0.00 | 0.184 | 1.99 | never |
| 4 | 0 | 0.02 | no channel | 10 | 0.90 | 0.10 | 0.00 | 1.00 | 0.962 | 2.15 | 5.69 |
| 4 | 0 | 0.02 | channel | 10 | 0.90 | 0.10 | 0.00 | 0.90 | 0.179 | 1.75 | 6.94 |

Reversal ratio: the loser's charge at t = 20 over its peak. "Never" means the site did not
latch within 20 Compton times in any trial.

**Reading.**

- **The first sink drains the second only through the channel.** Without it, both sites latch
  in every trial, at both lobe ratios, and the loser keeps its charge (ratio 0.87 to 1.00). With
  it, the winner's region grows over the loser's site, empties its well, and the loser's charge
  decays with the relaxation time (ratio 0.18 to 0.30). That is the reversal JB described, and
  it needs the sinks to interact; fixed-width sinks do not.
- **Vacuum noise decides ties, and it also makes double records.** At r = 1 with no offset the
  tie is exact without noise and is broken by noise at close to one half. But noise kicks the
  loser over the threshold before the winner's region reaches it: with the channel, the second
  latch is prevented in every deterministic run with an offset, and in almost none of the noisy
  ones at the stronger noise. At r = 4 the weak-lobe site never latches by itself with the
  channel on and little noise, and latches in nine trials of ten at the stronger noise.
- **Offset and noise trade off as a sigmoid.** At r = 1 the probability that the favoured site
  latches first rises with its head start and falls with the noise, which is the shape
  registered, with the rise steeper than the rough numbers given.
- **Double records are the cost case at this level.** For one quantum the model produces two
  latched sites in every no-channel run and in most channel runs with noise. A single record
  needs the sinks to interact faster than they charge, and the only interaction the model has
  is the channel covering the loser.

**Predictions against outcome.** G1: site 1 first at t = 2.9 (predicted 3 to 4) and the channel
reversal 0.18 (predicted 0.3 to 0.5, so lower); **but site 2 did latch by itself without the
channel, at t = 19, because the bright lobe's spreading feeds it**, which the prediction
missed. G2 held exactly (tie, both latch). G3 held in shape; the rise with offset was steeper
than the rough Φ estimate. G4 held: double records in every no-channel run and in most noisy
channel runs. **G5 missed narrowly**: Klein diagnostic max over all runs while the remainder exceeds 0.05: 0.1219; trials above the 0.10 flag: 3 of 168   -> FLAGGED (prediction G5 missed); see the cap sensitivity. The cap sensitivity shows the well is not the
source, since the fraction and every outcome are unchanged with the well capped at 1.0 mc²:

- r = 1, offset 0: cap 1.5 -> winner 2, both 1, reversal nan, Klein 0.122; cap 1.0 -> winner 2, both 1, reversal nan, Klein 0.122   -> same outcome
- r = 1, offset 0.5: cap 1.5 -> winner 0, both 0, reversal 0.28, Klein 0.085; cap 1.0 -> winner 0, both 0, reversal 0.28, Klein 0.085   -> same outcome
- r = 1, offset 1: cap 1.5 -> winner 0, both 0, reversal 0.30, Klein 0.083; cap 1.0 -> winner 0, both 0, reversal 0.30, Klein 0.083   -> same outcome
- r = 4, offset 0: cap 1.5 -> winner 0, both 0, reversal 0.18, Klein 0.042; cap 1.0 -> winner 0, both 0, reversal 0.18, Klein 0.042   -> same outcome

The admixture comes with the strong absorber and its growing width, and it does not change
which site wins, whether both latch, or the reversal.

## Part H (added 2026-09-07): an uneven split at the beam splitter from a spinor asymmetry

Source: `results/part_h_splitter.csv`, `part_h_splitter_slow.csv`, `part_h_targeted_trials.csv`,
`fig10_splitter.png`, `fig10_splitter_slow.png`. Registered in the pre-registration's Part H
addendum before the run, with two later addenda: the slow-charging run, added because the top of
the registered threshold sweep produced no record at all rather than a slow one; and a
fifty-event check at the edge of the window. JB's proposal: the two spinors carry unequal energy
shares, the symmetry breaking ε, and at a beam splitter the second spinor's mismatch with the
splitter's random phases decides which output gets the larger share, so the split is uneven by ε
on every event and even on average. Two forms of the dependence: a switch and a smooth cosine.
Each event's lobes then enter the two-sink race of Part G.

**Registered run, τ_c = 10, twenty Compton times.** Entries: double-record rate (larger share's
side recorded first), ten events per cell.

| form | C_th | ε = 0 | ε = 0.05 | ε = 0.1 | ε = 0.2 |
|---|---|---|---|---|---|
| switch | 0.015 | 1.00 (0.40) | 1.00 (0.40) | 1.00 (0.70) | 1.00 (0.90) |
| switch | 0.03 | 0.90 (0.50) | 0.90 (1.00) | 1.00 (1.00) | 0.30 (1.00) |
| switch | 0.06 | no record | no record | no record | no record |
| smooth | 0.015 | 1.00 (0.60) | 1.00 (0.75) | 1.00 (0.70) | 1.00 (0.80) |
| smooth | 0.03 | 1.00 (0.50) | 1.00 (0.60) | 1.00 (0.65) | 0.60 (1.00) |
| smooth | 0.06 | no record | no record | no record | no record |

**Slow-charging addendum, τ_c = 40, forty Compton times.**

| form | C_th | ε = 0 | ε = 0.05 | ε = 0.1 | ε = 0.2 |
|---|---|---|---|---|---|
| switch | 0.06 | 0.90 (0.30) | 1.00 (1.00) | 0.60 (1.00) | 0.00 (1.00) |
| switch | 0.12 | 1.00 (0.60) | 0.70 (1.00) | 0.40 (1.00) | 0.00 (1.00) |
| smooth | 0.06 | 1.00 (0.50) | 1.00 (0.70) | 1.00 (1.00) | 0.50 (1.00) |
| smooth | 0.12 | 1.00 (0.50) | 1.00 (0.80) | 0.40 (0.90) | 0.30 (1.00) |

**Targeted check, switch form, τ_c = 40, C_th = 0.12, fifty events each:** ε = 0.15: double records 1 of 50 (rate 0.02; 95 % upper bound nan if zero), larger share first 1.00, visibility 0.9887; ε = 0.2: double records 0 of 50 (rate 0.00; 95 % upper bound 0.058 if zero), larger share first 1.00, visibility 0.9798

The per-event split's fringe visibility, 2√(f(1 − f)), is √(1 − ε²) for the switch form (0.9987,
0.995, 0.980 at ε = 0.05, 0.1, 0.2) and about 1 − ε²/4 for the smooth form (0.9994, 0.997, 0.990),
against the 0.98 measured with single photons in 1986.

**Reading.** The asymmetry does what it was brought in to do, and only in one corner. With the
registered charging rate the double-record rate stays near one for every ε below 0.2 and falls
to 0.3 at ε = 0.2 in the switch form, 0.6 in the smooth form. With slow charging the switch form
reaches zero doubles in ten events at ε = 0.2 and none in fifty at the check, which bounds the
rate at 0.06 with 95 % confidence; the smooth form stays at 0.3 to 0.5 there, and at ε = 0.1 both
forms still give doubles in 0.4 to 0.6 of events. So the model's window is the switch form at
ε ≈ 0.2 with slow charging, and ε = 0.2 is exactly where the switch form's visibility meets the
1986 value of 0.98. The window has no width to spare on either side: at ε = 0.15 the check gave
one double in fifty, and the visibility bound admits nothing above 0.2. The anticorrelation bound
of 0.01 was not reached by this check and would need three hundred events; that is stated, not
rounded.

**Predictions against outcome.** P1 held (0.9 to 1.0 at ε = 0). P2 held in the slow run and
missed in the registered run's top threshold, where no site latched at all (the addendum records
why). P3 held everywhere a record formed: the smooth form's rate was at or above the switch
form's. P4 held, steeper than expected: the larger share's side is recorded first in every event
from ε = 0.05 in the switch form once charging is slow. P5: a window exists, at ε = 0.2 with slow
charging, and only there.

## Files

`results/part_a_uniform.csv`, `part_b_two_lobes.csv`, `ordered_bulk.csv`, `flight_arrivals.csv`,
`flight_zitter.csv`, `verification.txt`, `run_info.txt`; figures `fig1_part_a.png`,
`fig2_two_lobes.png`, `fig3_ordered_bulk.png`, `fig4_flight.png`, `fig5_zitter.png`; Part D:
`part_d_visibility.csv`, `part_d_visibility_ordered.csv`, `part_d_bound_amplitude.csv`, `fig6_interference.png`; Part E: `part_e_cost.csv`,
`part_e_cost_aligned.csv`, `fig7_cost.png`; Part F: `part_f_sink_summary.csv`,
`part_f_sink_series.csv`, `fig8_sink.png`; Part G: `part_g_race_summary.csv`,
`part_g_race_trials.csv`, `fig9_race.png`; Part H: `part_h_splitter.csv`, `part_h_splitter_slow.csv`,
`part_h_targeted_trials.csv`, `fig10_splitter.png`, `fig10_splitter_slow.png`.
