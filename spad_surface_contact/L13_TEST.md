# The test of L13: the model's deposit rule against quantum mechanics in a driven detector

**Date:** 2026-09-08 (evening), at JB's instruction "proceed with the test of L13". The closed forms
in §2 were written before `l13_qm_comparison.py` was run and the script must reproduce them.

## 0. Read this first

- **L13 was wrong in one clause.** It said quantum mechanics predicts no dependence of a
  detector's response on any phase of the incident particle, driven or not. For *coherent* light
  on a coherently prepared absorber that is false: the second pulse's absorption depends on its
  phase relative to the polarization the first pulse left, decaying with T2. This is coherent
  control, measured since the early 1990s; in semiconductor quantum wells an antiphase second pulse
  removes the excitons the first created (Heberle, Baumberg and Köhler, Phys. Rev. Lett. 75, 2598,
  1995, checked today). The clause holds only for single quanta in Fock states, where the field has
  no phase. Ledger L14 records the correction.
- **So the test splits in two.** For coherent light both the model and quantum mechanics predict a
  fringe, and they differ in its *form*: sign, order in the field, and the population factor
  (§3). For single quanta, quantum mechanics predicts no fringe and an absorption proportional to
  the population difference; the model predicts a fringe only for a quantum whose offset δ_p is
  definite relative to the pump, a variable no source controls and which the model has not defined
  for a photon (2026-09-06 ledger: photons have no lower spinor).
- **Verdict, stated before the tables.** Where the test can be done, with coherent light, existing
  experiments already decide it against the model's deposit rule as written: the rule cannot give
  net stimulated emission in antiphase, and its fringe is second order in the probe field where the
  measured one is first order. Where the model's prediction is distinctive, single quanta of
  definite offset, no experiment exists and none can be designed until the model says what sets a
  quantum's offset. For a passive detector the two agree, as Part D found.
- No literature search located a phase-locked pulse-pair experiment at the single-photon level
  (a Fock-state probe on a coherently prepared absorber within T2); heralded single-photon
  absorption by a single atom has been done without a prepared coherence. The absence is an
  absence in my search of today, not a result.

## 1. What L13 predicted

Part D: the pixel share follows the particle's offset relative to the pulse phase (1.000 against
Born's 0.800 for antiphase pixels, an in-phase particle), and the offset-averaged efficiency is
0.36 against 0.67, both decaying with R(t) = e^(−t/T2). The deposit weight behind both numbers is
the second study's A1, w = cos²(Δφ/2) = ½[1 + cos Δφ], with the bulk's phase spread giving
½[1 + R(t) cos Δφ] on average (second study, C5; Part D, Table 10).

## 2. Quantum mechanics for the same situation (closed forms, written before the run)

**Two-level absorber, semiclassical.** Resonant pulses in the rotating frame; a pulse of area θ
and phase φ rotates the Bloch vector about the axis (cos φ, sin φ, 0). From the ground state,
pulse 1 (θ₁, φ₁) gives u = −sin θ₁ sin φ₁, v = sin θ₁ cos φ₁, w = −cos θ₁. Dephasing for a delay t
multiplies u and v by e^(−t/T2) (T1 neglected). Pulse 2 (θ₂, φ₂) then changes the excited
population by

  Δρ_ee = ½ [ cos θ₁ (1 − cos θ₂) + sin θ₁ sin θ₂ e^(−t/T2) cos(φ₂ − φ₁) ].

The first term is ordinary absorption, weighted by the population difference cos θ₁ = ρ_gg − ρ_ee.
The second is the interference of pulse 2's field with the polarization pulse 1 left: signed,
first order in each field, decaying with T2. Checks: θ₁ = θ₂ = π/2 gives the Ramsey fringe
ρ_ee = ½[1 + e^(−t/T2) cos Δφ]; θ₁ = π/2 with a weak probe gives Δρ_ee = ½ θ₂ e^(−t/T2) cos Δφ,
pure interference, negative (net emission) in antiphase, zero in quadrature; θ₁ = 0 gives θ₂²/4.
The ratio of the fringe amplitude to the phase-averaged absorption is
sin θ₁ sin θ₂ e^(−t/T2) / [cos θ₁ (1 − cos θ₂)] ≈ 2 tan θ₁ e^(−t/T2) / θ₂ for a weak probe: it
*grows* as the probe weakens.

**Single quantum, Jaynes–Cummings.** Atom prepared as cos(θ₁/2)|g⟩ − i e^(−iφ₁) sin(θ₁/2)|e⟩, one
mode with coupling g, resonant. For a Fock state |1⟩ the amplitudes |g,1⟩ → |e,0⟩ (absorption) and
|e,1⟩ → |g,2⟩ (stimulated emission) end in different states and do not interfere, so the net
photon change is independent of φ₁ and, to order (gt)², equals (gt)² [2 sin²(θ₁/2) − cos²(θ₁/2)]:
absorption proportional to the population difference, with the factor 2 from stimulated plus
spontaneous emission into the mode. For a coherent state |α⟩ of mean photon number 1 and phase
φ₂ the two paths interfere and the net change carries a term ∝ |α| sin θ₁ cos(φ₂ − φ₁).

**The model's rule for comparison.** Deposit = (ordinary deposit) × ½[1 + R(t) cos Δφ]: never
negative, second order in the probe field (it multiplies the ordinary deposit), no population
factor, and a fringe-to-average ratio equal to R(t) whatever the probe strength.

## 3. Predictions for the script

- P14: the numerical two-level evolution reproduces the closed form for Δρ_ee to 10⁻¹² at
  θ₁ ∈ {0.3, π/2}, θ₂ = 0.2, t/T2 ∈ {0, 0.3, 1, 3}, over Δφ.
- P15: Jaynes–Cummings, gt = 0.3, θ₁ = π/2: the Fock-state net photon change is flat in φ to
  10⁻¹² and equals about (gt)² × ½ = +0.045 (net emission, since 2 sin² − cos² = ½ at π/2);
  the coherent-state change shows a fringe of amplitude about gt |α| sin θ₁ = 0.3.
- P16: fringe-to-average ratio at θ₁ = 0.3, t = 0: quantum mechanics 2 tan(0.3)/θ₂ = 6.2 at
  θ₂ = 0.1, 0.62 at θ₂ = 1 (exact form in the script); the model 1.0 at every θ₂.

## 4. Sources

- Heberle, Baumberg, Köhler, "Ultrafast coherent control and destruction of excitons in quantum
  wells", Phys. Rev. Lett. 75, 2598 (1995): https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.75.2598
- Piro et al., "Heralded single-photon absorption by a single atom" (arXiv:1004.4158): absorption
  without a prepared coherence; located today, not read in full.
- The record's Paper 1 §6.1 for the silicon timescales; this folder's Part D for the model's
  numbers.

## 5. Results (`results/verification_l13.txt`, `l13_two_level.csv`, `l13_jaynes_cummings.csv`, `fig7_l13.png`)

| Check | Result |
|---|---|
| P14, two-level numeric against the closed form of §2 | worst difference 2 × 10⁻¹⁶, pass |
| P15, Fock \|1⟩ on a π/2-prepared atom, gt = 0.3 | net photon change +0.041 (predicted about +0.045), spread over phase 0: no fringe, net emission, pass |
| P15, coherent state of mean photon number 1, same atom | fringe amplitude 0.27 (predicted about 0.3), pass |
| P16, fringe / phase-averaged absorption, θ₁ = 0.3, t = 0 | quantum mechanics 20.6, 6.2, 2.1, 0.57, 0.20 at θ₂ = 0.03, 0.1, 0.3, 1, 2; the model 1.0 throughout, pass |

In the middle panel of `fig7_l13.png` the coherent-state fringe is shifted by a quarter turn
relative to the model's curve: the phase origin of a coherent state against an atomic coherence is
a convention (the −i of the preparation pulse), and only the presence, amplitude and sign of the
fringe are compared, not its origin.

## 6. The comparison, in one table

| | quantum mechanics (§2) | the model's rule (A1, Part D) |
|---|---|---|
| coherent probe, in phase | absorption raised by a term first order in the probe field | deposit × 1 |
| coherent probe, antiphase | net stimulated emission (deposit negative); measured as exciton destruction (Heberle 1995) | deposit × 0, never negative |
| coherent probe, quadrature | ordinary absorption × (ρ_gg − ρ_ee), zero for a π/2 preparation | deposit × ½ |
| fringe / average against probe strength | grows as 1/θ₂ for a weak probe | constant R(t) |
| single quantum, Fock state | no fringe; absorption ∝ ρ_gg − ρ_ee; net emission at inversion ½ | fringe of visibility R(t) if the quantum's δ_p is definite; ½ on average |
| decay of the fringe | e^(−t/T2), the same | R(t) = e^(−t/T2), the same |
| passive detector | ordinary absorption | ordinary deposit, random atom in phase (Table 1 of the second study) |

## 7. Verdict

The test of L13 can be carried out today only with coherent light, and there it is already
decided: the phase-locked pulse-pair experiments show a signed, first-order fringe with net
emission in antiphase, and the model's deposit rule gives an unsigned, second-order one with no
population factor. The two share only the T2 decay. To pass, the rule would have to become
(ρ_gg − ρ_ee) × ordinary + sin θ₁ sin θ₂ R(t) cos Δφ, which is the golden rule's coherent term; the
phase mechanism would then add nothing for coherent light, and for single quanta it would keep
only the fringe that a definite per-quantum offset δ_p would give, a hidden variable no source
controls and which the model has not defined for a photon. For a passive detector, the case of
every SPAD that has ever counted a photon, the model and quantum mechanics agree, as Part D found.
L13 is corrected in the ledger as L14; the Part D numbers stand as what the model says.

## 8. Revision after Part E (2026-09-08, night)

JB's objection that |A| ≠ |B| is right, and Part E measures the amplitudes on a bound Dirac
electron: a_L = 0.06 to 0.30 for the lower row (about 10⁻³ for silicon), a₋ = 10⁻³ to 0.07 for the
negative-energy component (about 10⁻⁸ for silicon). With those amplitudes the model's weight is
[a² + b² + 2ab cos Δφ]/(a + b)², whose cross term 2ab cos Δφ has the form of §2's interference
term with the two second-spinor amplitudes as A and B. The verdict of §7 stands and is
strengthened: the sign, the order in the field and the population factor still differ, and the
visibility 2ab/(a² + b²) is a further suppression, zero for a physical incident electron in the
negative-energy reading and for any photon, and about 0.02 for a slow electron in the lower-row
reading against a silicon bound electron. The bound electron's second spinor is not a free phase:
it is in quadrature with the first, fixed, and the same for every identical atom (ledger L15).
