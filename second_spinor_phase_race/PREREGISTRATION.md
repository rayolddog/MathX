# Pre-registration: the second-spinor phase race at a detector

**Date:** 2026-09-07. Written before any run; `RESULTS.md` cites this file and does not change
the definitions below. Predicted numbers were computed in scratch from the closed forms and a
small Monte Carlo before the scripts were written, and the scripts must reproduce them.
**Exploration:** `MathX/second_spinor_phase_race/`.
**Related:** `../kuramoto_one_two_rotations/` (the first study: the first-spinor match followed
by an Adler pull on the second spinor; that pull is switched off here). The 1+1-D Dirac
conventions of the free-flight part are those of the record's `../dk_origin/dirac_arnold_tongue.py`
(ħ = c = m = 1, H = σ_x p + σ_z m, velocity ψ†σ_x ψ, zitterbewegung frequency 2).

## 0. Read this first

**The model (JB, 2026-09-06 and 2026-09-07), in his words as far as possible.** Each particle
carries two coupled helical spinors rotating in opposite directions. In an interaction, one spinor
of the incident particle couples at once with the first spinor of the other particle; as a result
the second spinor's phase is rotated as well, and the original phase difference between the
particle's first and second spinors is maintained. In this first simulation there is no
synchronization of the second spinors. The atoms of the bulk have a random distribution of
second-spinor phases, hence a random distribution of first–second phase differences. A large number
of atoms sit in the aperture of the detector. The interaction of the incident particle with that
distribution is based on the distribution of deposited energy, and the deposited energy is highest
where the second spinors are closest in phase. In free flight the Bohm–Dirac motion of the Weyl
spinors is what is happening, and it is used for the free-flight part.

**Contrast with the record.** The record's closed-system script matched the zitterbewegung
*frequency* of two confined particles (an Arnold tongue in coupling and detuning). Here the first
spinor's *phase* is matched, at once, and each particle's internal offset is carried through.

**Mode.** The model is built as described. Departures from quantum mechanics or from experiment
are kept in `LEDGER.md`, dated, and are not argued in this file or in `RESULTS.md`. What this file
does name are the consequences of the model itself, including the ones that cost it, so that they
are computed and not skipped.

**Choices the description leaves open, fixed here before the run. Audit these first.**

- **A1, the form of "closest in phase".** w(Δφ) = cos²(Δφ/2): 1 when the second spinors coincide,
  0 when they are opposite, the interference form of two amplitudes at relative phase Δφ.
  Sensitivity: a band, w = 1 for |Δφ| ≤ ε with ε = 0.35 rad (the first study's band), else 0.
- **A2, what "based on the distribution of deposited energy" means for the record.** Three
  readings are run (D5). Primary: one atom takes the quantum, with probability proportional to
  its share of the deposited energy (the closest atom is the most likely to fire).
- **A3, how the incident wave's amplitude enters.** Two readings, both run: the field reading,
  in which the energy available at an atom is the wave's intensity there (Part B); and the
  trajectory reading, in which the particle is where its Bohm–Dirac trajectory puts it and the
  race is among the atoms within its reach (Part C).
- **A4, which components are the two spinors in free flight.** The positive- and
  negative-energy components of the Dirac spinor, which rotate as e^(−iEt) and e^(+iEt): opposite
  directions, coupled by the mass, with a relative phase that the free equation does not change.
- **A5, the bulk's phases are drawn fresh for each incident particle.** The model gives the
  atoms' second spinors no dynamics, so a fixed bulk and a fresh one differ only in whether the
  same atom wins twice; the fresh draw is the simpler statement of "random".

## 1. Definitions (fixed)

- **D1, particle.** Two helical spinors with opposite rotation, φ₁(t) = +ωt + a and
  φ₂(t) = −ωt + b. The internal offset is δ = b − a, the second spinor's phase minus the first's
  at the instant of contact. The "second-spinor phase" is b.
- **D2, atom.** Atom i of the bulk has the same structure with a_i, b_i and δ_i = b_i − a_i.
  In a random bulk δ_i is uniform on (−π, π], independent across atoms. The first-spinor phases
  a_i drop out after D3 and are not sampled.
- **D3, contact.** The particle's first spinor is shifted to the atom's, and the same shift
  Δ_i = a_i − a_p is applied to both spinors of the particle, so δ_p is maintained. No pull acts
  on the second spinor. The second-spinor mismatch at atom i is therefore
  Δφ_i = (b_p + Δ_i) − b_i = δ_p − δ_i: the difference of the two internal offsets. For a
  random bulk Δφ_i is uniform on the circle whatever δ_p is.
- **D4, deposited energy.** E_i = E_p · A_i · w(Δφ_i) / Σ_j A_j w(Δφ_j), with w from A1 and A_i
  the energy available at atom i: A_i = 1 for a uniform wave (Part A), A_i = |ψ(x_i)|² in the
  field reading (Part B), A_i = 1 for the atoms within reach and 0 otherwise in the trajectory
  reading (Part C).
- **D5, the record.** R1, shared: every atom keeps E_i and the record is the map.
  R2, sampled (primary): one atom takes E_p, atom i with probability E_i/E_p.
  R3, largest: the atom with the largest E_i takes E_p. R3 is the deterministic limit of R2.
- **D6, aperture.** N atoms. Part A: uniform wave, N from 1 to 10⁴. Part B: two lobes of
  intensity r and 1 (amplitude ratio √r), n atoms in each, r ∈ {1, 2, 4, 9, 16}, n from 1 to 10³.
  Born share of lobe 1: r/(1 + r). Exponent γ defined by P₁/P₂ = (√r)^γ; Born gives γ = 2.
- **D7, ordered bulk.** δ_i drawn from a von Mises distribution of mean μ and concentration κ,
  κ ∈ {0, 1, 4} (κ = 0 is the random bulk), with mean resultant R_κ = I₁(κ)/I₀(κ). The two-half
  detector: N atoms with μ = 0 beside N atoms with μ = π, uniform wave, and the incident δ_p swept.
- **D8, free flight.** 1+1-D Dirac in the record's conventions, x the coordinate across the
  aperture. Bohm–Dirac velocity v = ψ†σ_x ψ / ψ†ψ. Initial positions drawn from ρ(x, 0); each
  trajectory is integrated to the detector time t_d and its arrival x(t_d) recorded. Two-lobe
  state: ψ(x, 0) ∝ [√r g(x − x₁) + g(x + x₁)] u₊ with g a Gaussian of width σ = 1, x₁ = 4, r = 4,
  u₊ = (1, 0) the rest-frame positive-energy spinor, t_d = 4. Atoms on a lattice of spacing
  d = 0.05; reach L_c = 0.5, so 21 atoms within reach of an arrival point.
- **D9, the internal offset in free flight.** ψ(x, 0) = g(x)[cos(θ/2) u₊ + sin(θ/2) e^(iδ) u₋],
  u₋ = (0, 1): first spinor the positive-energy part, second the negative-energy part, internal
  offset δ. θ = 60°, δ ∈ {0°, 90°}; also θ = 0 as the reference with no second spinor.

## 2. Predictions (fixed before the run)

**Part A, uniform wave, N atoms, random bulk.** The mismatches |Δφ_i| are uniform on [0, π].
Best mismatch m* = min_i |Δφ_i|: P(m* > a) = (1 − a/π)^N, mean π/(N + 1). Energy weight of the
best atom w_max = cos²(m*/2).

| | closed form | N for 95 % confidence |
|---|---|---|
| w_max ≥ 0.95 | 1 − 0.8565^N | 20 |
| w_max ≥ 0.99 | 1 − 0.9362^N | 46 |
| band: some atom within ε = 0.35 | 1 − 0.8886^N | 26 |

R1: the best atom's share of the quantum is w_max/Σw ≈ 2/N (Σw has mean N/2).
R2: the fired atom's mismatch has density cos²(Δφ/2)/π on the circle, whatever N: mean
|Δφ| = π/2 − 2/π = 0.934 rad (53.5°); probability of being within the 0.35 rad band 0.221,
against 0.111 for a random atom. R3: the fired atom is the best, mean mismatch π/(N + 1).

**Part B, two lobes, r = 4 (Born 0.800).** Scratch Monte Carlo, 2 × 10⁵ trials:

| n per lobe | 1 | 2 | 5 | 10 | 100 | 1000 |
|---|---|---|---|---|---|---|
| R2: P₁ | 0.705 | 0.759 | 0.788 | 0.795 | 0.800 | 0.800 |
| R2: γ | 1.26 | 1.66 | 1.89 | 1.95 | 2.00 | 2.00 |
| R3: P₁ | 0.791 | 0.923 | 0.997 | 1.000 | 1.000 | 1.000 |
| R3: γ | 1.92 | 3.6 | 8.3 | 16 | ∞ | ∞ |

Exponent at n = 1 against r: R2 gives 1.27, 1.25, 1.23, 1.22 for r = 2, 4, 9, 16; R3 gives
2.36, 1.92, 1.67, 1.57. So R2 tends to Born's 2 from below as the aperture fills, because the
phase factor averages to ½ in both lobes; R3 tends to a step. The record's own race (Paper 1)
gave 1.56 with interval 1.44 to 1.69.

**Ordered bulk.** Mean deposit per atom ½[1 + R_κ cos(δ_p − μ)], R₁ = 0.446, R₄ = 0.864. The
two-half detector under R2, large N: P₁(δ_p) = ½[1 + R_κ cos δ_p], a Malus-type law in the
particle's own internal offset. Under R3 the same curve steepens toward a step of width ~1/κ.

**Part C, free flight.** (i) The fraction of arrivals at x > 0 equals ∫_{x>0} ρ(x, t_d) dx to
sampling error (Bohm equivariance); that integral equals r/(1 + r) = 0.800 up to the overlap of
the two lobes at t_d, which is reported. (ii) The record position under either rule is the
arrival point plus a uniform blur on [−L_c, L_c] and nothing else, since the atoms' phases are
independent of their positions; the winner's energy weight has mean 0.75 under R2
(E[w²]/E[w]) and about 0.99 under R3 (21 atoms). (iii) For D9 the Bohm–Dirac velocity at the
packet centre is v(t) = sin θ · cos(2t + δ): zitterbewegung at frequency 2 with the internal
offset δ as its phase and sin θ = 0.866 as its amplitude; position amplitude sin θ / 2 = 0.433
Compton lengths. δ is maintained through the flight up to the dispersion of the packet, whose
drift of the fitted phase over t = 0 to 12 is recorded. θ = 0 gives no zitterbewegung.

## 3. Consequences to be recorded, including the ones that cost the model

- **C1, R1 spreads the quantum thin.** The best atom's share is about 2/N; at N = 10⁴ it is
  2 × 10⁻⁴ of the quantum. A record needs the quantum at one site, so R1 alone makes no record.
- **C2, under R2 a random bulk hides the particle's phase.** Δφ_i is uniform whatever δ_p, so
  no statistic depends on δ_p; the fired atom is typically 54° off in phase; and the site
  statistics are the intensity weighting A_i with the phase factor averaged out.
- **C3, under R3 the brighter lobe always fires once the aperture is full.** γ → ∞. At n = 1,
  r = 4 the exponent is 1.92, near Born's 2, but it falls with r (1.57 at r = 16).
- **C4, the band form has a miss rate** (1 − ε/π)^N; 95 % efficiency needs N ≥ 26 atoms in reach.
- **C5, the distinctive prediction.** In an ordered bulk the deposit and the firing probability
  follow ½[1 + R_κ cos(δ_p − μ)]: a phase-ordered detector reads the incident particle's internal
  offset. This is where the model says something the golden rule does not, and where it would be
  tested. What is known about phase-ordered detectors is in the ledger, not here.
- **C6, in the trajectory reading the site statistics are Bohm's, not the race's.** The race
  decides only which atom within reach fires and how much energy it takes.

## 4. What will be recorded

`results/`: `part_a_uniform.csv`, `part_b_two_lobes.csv`, `ordered_bulk.csv`,
`flight_arrivals.csv`, `flight_zitter.csv`, `verification.txt`, `run_info.txt`;
figures `fig1_part_a.png`, `fig2_two_lobes.png`, `fig3_ordered_bulk.png`,
`fig4_flight.png`, `fig5_zitter.png`. `RESULTS.md` states the uncertainty and the assumptions
before the tables. `LEDGER.md` holds the dated side ledger.

## 5. Method note

Definitions and predictions were fixed in this file before the scripts were written and run.
Any later change to a definition is to be made as a dated addendum below, not by editing
§1 to §3.

## 6. Addenda

**2026-09-07, after the first run of `bohm_dirac_flight.py`.** Prediction §2 (iii) gave the
zitterbewegung frequency of the trajectory as 2 (units mc²/ħ). That is the value for a particle at
rest; the registered packet (σ = 1) has a momentum spread of 1/(2σ) and beats at about 2⟨E⟩,
with ⟨E⟩ the packet's mean energy, and its beat in position space dephases as the momentum
components separate. The first run showed the miss (fitted 2.30 to 2.42 against 2). The
definitions are unchanged; two things are added. (1) The claim that the internal offset is
maintained is tested directly: the relative phase of the ±E components at each momentum, minus
2E(p)t, must equal its value at t = 0 at all later times. (2) A wider packet, σ = 3, is run as a
sensitivity case, with the expectation that the beat approaches 2 as the packet widens
(2⟨E⟩ = 2.03 at σ = 3). The σ = 1 case remains the registered one and its miss is recorded in
`RESULTS.md` under "Predictions against outcome".

**2026-09-07, later, after JB's clarification of the interference term and his choice of the
record rule.** JB: the deposited energy is highest at the smallest mismatch *because* the second
spinors interfere, the incident particle's with the bound electron's; and sampling by deposit is
the chosen rule. Two things follow, written before Part D was run.

- **The record rule is fixed: R2.** R3 stays in the tables above as the deterministic limit but
  is no longer a candidate reading.
- **The cos² weight of A1 is the interference term for two amplitudes of equal size:**
  |1 + e^(iΔφ)|² = 4 cos²(Δφ/2). What the description leaves open is the size of the bound
  electron's second-spinor amplitude b relative to the incident amplitude a at the atom. Two
  readings, both run as Part D (`interference_visibility.py`):
  - **D10, gating reading.** The energy deposited is the particle's, and the interference gates
    how much: E_i ∝ A_i · (1 + V cos Δφ_i)/2 with V = 2ab/(a² + b²) ≤ 1 the visibility of the
    interference. V = 1 is A1 as run. V ∈ {1, 0.5, 0.2, 0}.
  - **D11, energy-density reading.** The deposit is the energy density of the superposition at
    the atom, E_i ∝ |a_i e^(iφ_p) + b e^(iφ_i)|² = A_i + β² + 2β√A_i cos Δφ_i, with β = b/a₀ the
    bound amplitude in units of the incident amplitude in the weaker lobe. β ∈ {0, 0.5, 1, 2}.
    Only atoms the wave reaches (A_i > 0) are in the aperture. If atoms with A_i = 0 were
    included, each would carry β² on its own and a record could land outside the wave with
    probability n_out β² / [n(r + 1 + 2β²) + n_out β²]; that consequence is stated, not run.

Predictions (closed forms and scratch Monte Carlo, r = 4, R2):

| D10, visibility V | 1 | 0.5 | 0.2 | 0 |
|---|---|---|---|---|
| fired atom's mean mismatch, π/2 − 2V/π | 53.5° | 71.8° | 82.7° | 90° |
| fired atom within the 0.35 rad band, (ε + V sin ε)/π | 0.221 | 0.166 | 0.133 | 0.111 |
| two lobes, P₁ at n = 1 | 0.705 | 0.787 | 0.798 | 0.800 |
| two lobes, P₁ at n = 10 | 0.795 | 0.799 | 0.800 | 0.800 |
| two lobes, P₁ as n → ∞ | 0.800 | 0.800 | 0.800 | 0.800 |
| ordered bulk, two halves: ½(1 + V R_κ cos δ_p), visibility at κ = 4 | 0.86 | 0.43 | 0.17 | 0 |

| D11, bound amplitude β | 0 | 0.5 | 1 | 2 |
|---|---|---|---|---|
| two lobes, P₁ at n = 1 | 0.800 | 0.771 | 0.698 | 0.543 |
| two lobes, P₁ as n → ∞, (r + β²)/(r + 1 + 2β²) | 0.800 | 0.773 | 0.714 | 0.615 |

So under D10 the visibility changes nothing about where the record lands and scales every phase
effect by V; under D11 the bound electron's own energy dilutes the intensity share toward ½ as
β grows. Neither reading moves C2 (a random bulk hides the particle's phase) or the form of C5
(the ordered-bulk law), whose visibility becomes V R_κ.

**2026-09-07, later again, after JB's addition of an energy cost for the phase change.** JB: add a
small energy cost of the phase change to the mechanism, shared equally by both particles, and ask
whether it changes the result much. Written before Part E was run.

- **D12, the cost.** At contact the particle's spinors are shifted together by Δ_i = a_i − a_p
  (D3). The shift costs c(Δ) = η (1 − cos Δ)/2 in units of the incident energy in the weaker
  lobe: nothing for no shift, η for a half-turn. Shared equally: the particle pays half out of
  what it can deposit and the atom pays the other half from its own energy. The record rule R2
  samples on the particle's deposit, E_i = max(0, A_i w(Δφ_i) − c(Δ_i)/2). An atom whose
  deposit would be negative cannot fire; a particle with no positive deposit anywhere is a miss.
  η ∈ {0, 0.05, 0.1, 0.2, 0.5}; "small" is η ≤ 0.1.
- **D13, the bulk's first-spinor phases return.** D2 dropped them because the match at contact
  made them irrelevant; the cost depends on the shift, so they matter again. Primary: random,
  Δ_i uniform on the circle. Sensitivity: aligned, all a_i equal, so Δ is a single angle set by
  the particle's own first-spinor phase against the bulk's; swept over 0°, 45°, 90°, 135°, 180°.

Predictions (r = 4, R2). Closed form ignoring the floor at zero, large n, random first spinors:
P₁ = (r − η/2)/(r + 1 − η) = 0.800, 0.803, 0.806, 0.813, 0.833 for η = 0, 0.05, 0.1, 0.2, 0.5;
with the floor the shift is a little smaller. The exact large-n value is the ratio of the two
lobes' mean truncated deposits, computed by quadrature in the script. Scratch Monte Carlo:

| η | 0 | 0.05 | 0.1 | 0.2 | 0.5 |
|---|---|---|---|---|---|
| P₁ at n = 1 (miss rate) | 0.705 (0) | 0.706 (0.002) | 0.709 (0.004) | 0.711 (0.008) | 0.721 (0.021) |
| P₁ at n = 10 | 0.795 | 0.798 | 0.800 | 0.806 | 0.820 |
| P₁ at n = 1000 | 0.800 | 0.803 | 0.806 | 0.811 | 0.826 |

Aligned first spinors, n = 1000: P₁ = (r − ηc)/(r + 1 − 2ηc) without the floor, c = (1 − cos Δ)/2;
at η = 0.1: 0.800, 0.806, 0.811 for Δ = 0°, 90°, 180°; at η = 0.5: 0.800, 0.827, 0.854 (scratch
Monte Carlo; no-floor 0.800, 0.833, 0.875).

So the cost is a threshold. It removes the atoms whose deposit cannot pay it, which are more
often in the weaker lobe, and so pushes the brighter lobe's share *above* Born, by about
η/(r + 1) at small η: about a percent at η = 0.1, three percent at η = 0.5. With aligned first
spinors the particle's own first-spinor phase enters the site statistics for the first time,
through the same threshold, by the same small amount. The random-bulk conclusion for the
second-spinor phase (C2) is unchanged; the ordered-bulk law (C5) is unchanged in form.
Recorded: `results/part_e_cost.csv`, `part_e_cost_aligned.csv`, `fig7_cost.png`.

**2026-09-07, Part F: the winning atom becomes a sink (JB's step of 2026-09-07, lightning reading;
run with my defaults at JB's instruction).** Written before `sink_flow.py` was run.

- **D14, the sink.** After the deposit, the winning atom at x_s becomes a well plus an absorber
  at its own site, both with the profile f(x) = exp[−(x − x_s)²/2L²]. The well is a scalar
  potential V = −U f added to the Dirac Hamiltonian (same sign on both components); the absorber
  removes norm at the rate Γ f. The collected fraction of the quantum is the absorbed norm,
  E_coll(t) = 1 − ∫ρ dx: each unit of norm carries the quantum's energy.
- **D15, strength from the deposit.** The initial strength is the deposit itself: Γ₀ = γ₀ w_win
  and U₀ = u₀ w_win, with w_win = 0.75 (the R2 winner's mean weight, Part C) and γ₀, u₀ the
  scales. Defaults γ₀ = 1 (in mc²/ħ), u₀ = 0.5 (in mc²).
- **D16, growth (the channel heats).** Both grow with the energy already collected:
  Γ(t) = Γ₀ [1 + κ E_coll(t)], U(t) = U₀ [1 + κ E_coll(t)], default κ = 10, with the well capped at
  1.5 mc² so it stays below the Klein threshold 2 mc². Channel extension, off by default:
  L(t) = L₀ + λ E_coll(t), λ ∈ {0, 10} Compton lengths.
- **D17, width.** L₀ = one lattice spacing, 0.05 Compton lengths (D8). Sweep {0.05, 0.5, 2}.
- **D18, state and site.** Primary: Part C's packet at t_d (r = 4, lobes at ±4, σ = 1), winner at
  the bright lobe's centre x_s = +4; the weak-lobe winner x_s = −4 as a case. Sensitivity: wide
  lobes (σ = 3 at ±10) so that spreading and flow can be told apart. Run for 40 Compton times
  after t_d, dt = 0.01, split-operator with the exact free propagator. "Near" and "far" are the
  norms on the sink's side of x = 0 and on the other side. A free run (Γ = U = 0) of the same
  state is the reference for spreading across x = 0.
- **Diagnostic.** The negative-energy fraction of the remaining wave, to flag a well that has
  entered Klein territory (flag above 0.10).

Runs (one parameter at a time around the base): base; γ₀ = 0.1, 10; u₀ = 0, 1.5; L₀ = 0.5, 2;
κ = 0; λ = 10; x_s = −4; wide lobes with the base; wide lobes with λ = 10.

Predictions.
- **F1, the start is absorption-limited.** Initial rate dE_coll/dt = Γ₀ ∫ρ f dx, about
  Γ₀ ρ(x_s) √(2π) L₀ ≈ 0.013 per Compton time at the base (ρ(4, t_d) ≈ 0.14).
- **F2, a strong narrow absorber reflects.** Collection is sublinear in γ₀ and may fall again at
  γ₀ = 10 (the absorber acts as a barrier): the quantum-Zeno side of an imaginary potential.
- **F3, a point sink does not swallow the quantum on this timescale.** For every fixed-width run
  E_coll(40) stays below the near lobe's initial share (0.79 for the bright lobe), and much of
  what is collected has arrived by the packet's own spreading rather than by any pull, because a
  well of width 0.05 and depth ≤ 1.5 mc² binds weakly (strength U L ≤ 0.075) and its reach is
  its width.
- **F4, the far lobe.** With λ = 0 the far side's norm follows the free reference to within the
  amount that spreads across x = 0 and is then eaten; nothing flows toward the sink from afar.
  With λ = 10 the channel reaches the far lobe only after collecting most of the near one, and
  then collects it: E_coll(40) rises above the near share in that run and in no other.
- **F5, the weak-lobe winner** collects less: its near share is 0.21.
- **F6, Klein diagnostic** stays below 0.10 in every run with the cap.

Recorded: `results/part_f_sink_summary.csv`, `part_f_sink_series.csv`, `fig8_sink.png`, and
verification lines (norm bookkeeping to 10⁻⁸; U = 0, Γ = 0 reproduces the exact propagator;
Γ ≠ 0, U ≠ 0 with κ = 0 conserves nothing but ∫ρ + E_coll; the initial rate against F1 to 2 %).

**2026-09-07, Part G: several sinks, the first draining the others, with vacuum noise (JB's
thoughts of 2026-09-07; run with my defaults at JB's instruction).** Written before `sink_race.py`
was run.

- **D19, two sites.** Site 1 at the bright lobe's centre x = +4, site 2 at the weak lobe's
  centre x = −4, in Part C's packet at t_d (σ = 1). Both deposits w = 0.75. The equal-lobe case
  r = 1 is run as well, with site 2 displaced toward the middle by 0, 0.5 or 1.0 Compton lengths
  to give a controlled imbalance (the density there is 1, 0.975 and 0.905 of site 1's).
- **D20, charging below threshold (the Josephson side).** Each site gathers a charge C_i from
  the norm inside its well, N_i = ∫ρ f_i dx, reversibly: dC_i/dt = N_i − C_i/τ_c with
  τ_c = 10 Compton times, so an unfed site loses its charge. Below threshold the site is a well
  only, U_i = min(1.5, u₀ w [1 + κ (C_i + E_i)]) of width L_i = L₀ + λ E_i, with u₀ = 0.5,
  κ = 10, L₀ = 0.05 and λ ∈ {0, 10} as in Part F. No absorber.
- **D21, the latch (the lightning side).** When C_i + η_i first reaches C_th = 0.03 the site
  latches: its absorber switches on, Γ_i = γ₀ w [1 + κ (C_i + E_i)] with γ₀ = 1, and E_i, the
  norm it has absorbed, grows from then on. A latched site stays latched. Two sites latching on
  the same step is a tie.
- **D22, vacuum noise.** η_i is an Ornstein–Uhlenbeck noise on each site's charge with
  correlation time 1 Compton time and strength s ∈ {0, 0.005, 0.02}, stationary sd s/√2:
  0.0035 and 0.014 against a threshold of 0.03. Ten trials per noisy setting, one otherwise.
- **D23, recorded.** Which site latches first and when; whether both latch (two records for one
  quantum); the loser's reversal ratio, its charge at t = 20 over its peak charge; the norm each
  site has absorbed by t = 20; the Klein diagnostic. Run for 20 Compton times after the deposit.

Predictions.
- **G1, r = 4, no noise.** Site 1 latches at about t = 3 to 4. Site 2 never latches by itself:
  its charge saturates near τ_c N₂ ≈ 0.02, below threshold. With the channel (λ = 10) site 1's
  region covers site 2 within a few Compton times of latching, site 2's well is emptied and its
  charge decays with τ_c: reversal ratio 0.3 to 0.5 at t = 20. Without the channel site 2 keeps
  its charge (ratio near 1). This is the "reversal" JB described, and it needs the channel.
- **G2, r = 1, no noise, no offset.** An exact tie: both latch on the same step, and both stay
  latched whether or not the channel is on. Two records for one quantum.
- **G3, r = 1 with noise.** P(site 1 latches first) = 0.5 at zero offset and rises with the
  offset along a sigmoid Φ(ΔC/(√2 σ_C)), with ΔC the imbalance in charge at the latch (about 5 %
  and 10 % of C_th for offsets 0.5 and 1.0) and σ_C the noise sd: roughly 0.62 and 0.73 at
  s = 0.005, 0.53 and 0.56 at s = 0.02.
- **G4, double records.** Without the channel both sites latch in nearly every r = 1 trial, since
  fixed-width sinks do not interact (Part F). With the channel the second latch is prevented only
  when the winner's region reaches the loser before the loser's own charge crosses threshold;
  at r = 1 the head starts are fractions of a Compton time, so most trials still give two
  records. This is the cost case at this level: two records for one quantum unless the sinks
  interact faster than they charge. The fraction of double records is recorded per setting.
- **G5.** The Klein diagnostic stays below 0.10 while the remainder exceeds 0.05.

Recorded: `results/part_g_race_summary.csv`, `part_g_race_trials.csv`, `fig9_race.png`.

**2026-09-07, Part H: an uneven split at the beam splitter from a spinor energy asymmetry, fed
into the two-sink race (JB's proposal of 2026-09-07; registered with my defaults at his
instruction).** Written before `splitter_asymmetry.py` was run.

- **D24, the splitter as a bulk.** The splitter's atoms carry random second-spinor phases. At the
  splitter the contact rule of D3 applies: the first spinor matches at once and the second
  spinor's mismatch Δφ_s against the atom met is uniform on the circle, fresh per event.
- **D25, the asymmetry.** The two spinors carry energy shares (1 + ε)/2 and (1 − ε)/2, ε the
  symmetry breaking, ε ∈ {0, 0.05, 0.1, 0.2}. The mismatch decides which output the larger share
  goes to, so the per-event split is f = ½ + (ε/2) g(Δφ_s) with two forms of g:
  the **switch**, g = sign(cos Δφ_s), and the **smooth** form, g = cos Δφ_s. The split is even on
  average in both.
- **D26, the race.** Each event's lobes carry f and 1 − f and enter the two-sink race of Part G
  unchanged (channel λ = 10, noise s = 0.005, twenty Compton times, sites at the lobe centres),
  with the charging threshold swept, C_th ∈ {0.015, 0.03, 0.06}, to move the ratio of charging
  time to suppression time. Ten events per cell of the sweep.
- **D27, recorded.** The double-record rate, the rate at which the larger share's side is
  recorded first, the mean per-event fringe visibility 2√(f(1 − f)) the same split would give in
  an interferometer, and the latch-time gap.

Predictions (closed form for the visibility; the race numbers are expectations from Part G).

| | ε = 0 | 0.05 | 0.1 | 0.2 |
|---|---|---|---|---|
| visibility, switch, √(1 − ε²) | 1 | 0.9987 | 0.995 | 0.980 |
| visibility, smooth, ≈ 1 − ε²/4 | 1 | 0.9994 | 0.9975 | 0.990 |
| share of events within 55/45 of even, switch | all | all | none | none |
| share of events within 55/45 of even, smooth | all | all | 0.33 | 0.16 |

- **P1.** At ε = 0 the double-record rate is that of Part G at equal lobes, near 1.
- **P2, the switch form.** The double-record rate falls with ε and with C_th, since a larger
  imbalance and slower charging both lengthen the gap the winner has to suppress the loser. Part G
  gave a rate of 0.5 at a ten percent density imbalance with this noise; the expectation is below
  0.5 at ε = 0.2 for the two larger thresholds and near zero at C_th = 0.06.
- **P3, the smooth form.** The rate stays above the switch form's at every ε and every C_th, by
  at least the share of near-even events in the table, because those events have no head start.
- **P4.** The larger share's side is recorded first with probability rising from ½ at ε = 0
  toward one at ε = 0.2 in the switch form, following the same sigmoid as Part G.
- **P5, the window.** With the double-record rate required below 0.01 and the visibility above
  the best measured single-photon value, the switch form has a window in ε only if the rate at
  ε = 0.2 reaches 0.01 at some threshold in the sweep; whether it does is the result.

Recorded: `results/part_h_splitter.csv`, `part_h_splitter_trials.csv`, `fig10_splitter.png`.

**2026-09-07, Part H addendum, after the first run.** At C_th = 0.06 no site latched in any
event: with τ_c = 10 the charge saturates near τ_c ρ w, below that threshold, so the top of the
sweep gave no record rather than a slow one, and "double record 0.00" there is not a result. The
slow-charging regime the registration meant needs the relaxation time raised with the threshold.
Added, as a second run of the same script: τ_c = 40 with C_th ∈ {0.06, 0.12}, forty Compton times,
both forms, the same ε list and event count. Prediction P2 is read against this run: the switch
form's double-record rate at ε = 0.2 should fall below the first run's 0.30, and the question of
P5 is whether it reaches 0.01.

**2026-09-07, Part H, targeted check, after the slow-charging run.** The switch form at ε = 0.2 gave
no double record in ten events at either slow threshold, which bounds the rate only at about 0.2.
To read P5 against the anticorrelation bound, fifty events are run at ε = 0.2 and fifty at
ε = 0.15 (switch form, τ_c = 40, C_th = 0.12, forty Compton times). Zero doubles in fifty bounds
the rate at 0.06 with 95 % confidence; the 0.01 bound itself would need three hundred events and
is not reached by this check, which is stated.
