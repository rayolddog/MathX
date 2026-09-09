# Pre-registration: the pair case, JB's model at two analyzers

**Date:** 2026-09-08 (night). Written before any run, at JB's instruction "proceed with the
microscopic model simulation of the pair case", after the exchange in which single-photon Bell
experiments were seen to use one-photon absorption at each detector, so that the correlations
must be carried by the pair's amplitudes. `RESULTS.md` cites this file and does not change the
definitions below. Predicted numbers were computed in scratch (4 × 10⁵ pairs) before the script
was written. **Exploration:** `MathX/pair_case/`. **Ledger:** `LEDGER.md`, continuing the SPAD
study's numbering (L17 onward).

## 0. Read this first

**The model (JB, 2026-09-06, as recorded that day).** Entanglement is equal phases: the two
particles of a pair carry the same second-spinor phase. At contact the first spinor locks to the
bulk's axis and the second keeps the pair value; the outcome is the projection of the second
spinor's phase on the bulk's axis. In the photon reading established on 2026-09-08 the
second-spinor phase is the polarization, the bulk with an axis is the analyzer, and the projection
rule is Malus's law or its deterministic limit.

**What this run can and cannot show, stated before it (C1).** A shared phase and a rule that
acts at each analyzer on that phase and the local setting is a local hidden-variable model in
Bell's sense, whatever the rule. Bell's theorem then bounds the CHSH combination by 2 for every
rule when every photon is detected. The run cannot give more, and it is registered as such: what
it records is the model's own correlation function against the quantum one, and the single route
the model has above 2, which is to let its phase-dependent efficiency (Part D of the SPAD study)
gate the detection, the detection loophole, closed by the experiments of 2015.

## 1. Definitions (fixed)

- **D1, the pair.** Two photons with the same hidden polarization angle φ, uniform on [0, π)
  (JB: entanglement = equal phases). Sensitivities: the second photon at φ + π/2 (the same
  magnitudes with sign changes; not tabulated); the two angles independent (the unentangled
  control, E = 0).
- **D2, the analyzer.** A polarizer at angle a with two ports: a phase-ordered bulk with axis a
  (second study, D7, κ → ∞). The photon's offset from the axis is Δ = φ − a.
- **D3, the outcome rule, three readings.** R-least (JB's projection, deterministic): port +
  if cos 2Δ ≥ 0, the port needing the smaller rotation, else −. R-Malus (the second study's R2):
  port + with probability cos² Δ. QM, the comparator: the joint distribution of the state
  (|HH⟩ + |VV⟩)/√2, P(s_a, s_b) = ¼[1 + s_a s_b cos 2(a − b)], sampled.
- **D4, detection.** Ideal (primary): every photon detected. Gated (the model's own
  phase-dependent efficiency, SPAD study Part D, Table 10): a photon that went to port ± is
  detected only if its weight, cos² Δ for + and sin² Δ for −, is at least g, with g ∈ {0, 0.5,
  0.6, 0.63, 0.75, 0.85, 0.9}. Recorded: the single-side efficiency η₁, the conditional efficiency
  η_c = P(coincidence)/P(single), and the CHSH value among coincidences S_c.
- **D5, settings and statistics.** CHSH at (a, a′, b, b′) = (0°, 45°, 22.5°, 67.5°), S = E(a, b) −
  E(a, b′) + E(a′, b) + E(a′, b′). The correlation E(θ) at θ = a − b from 0° to 90° in 7.5° steps.
  N = 4 × 10⁵ pairs per setting: standard error about 0.0016 in E and 0.003 in S.
- **D6, the spin version.** The same with period 360° and E_QM = −cos θ for the singlet gives the
  same CHSH values with θ halved; stated, not run.

## 2. Predictions (from scratch)

- **P1, ideal detection.** R-least: E(θ) = 1 − 4θ/π on [0°, 90°]: 1, 0.667, 0.5, 0.333, 0, −0.333,
  −0.5, −0.667, −1 at 0°, 15°, 22.5°, 30°, 45°, 60°, 67.5°, 75°, 90° (Bell's own sawtooth); S = 2.000,
  Bell's bound met exactly. R-Malus: E(θ) = ½ cos 2θ: 0.5, 0.433, 0.354, 0.25, 0, …; S = 1.414.
  QM: cos 2θ; S = 2.828. Independent angles: E = 0 at every θ, S = 0.
- **P2, gated detection.** R-least with g = 0.5 changes nothing (the + condition is cos² Δ ≥ ½).
  Above it the model saturates the Garg–Mermin bound S_c = 4/η_c − 2: g = 0.6 gives S_c = 2.69 at
  η₁ = 0.873, η_c = 0.853; g = 0.63 about 2.83 at η_c about 0.83; g = 0.75 gives S_c = 4.00 (the
  algebraic maximum, perfect correlations among coincidences) at η₁ = 0.666, η_c = 0.625; g = 0.85
  and 0.9 give 4.00 at η₁ = 0.51 and 0.41. R-Malus does not saturate the bound: g = 0.5 gives
  S_c = 2.62 at η₁ = 0.82; g = 0.6 gives 3.16 at 0.75; g ≥ 0.75 gives 4.00.
- **P3, the physical reading of P2.** The quantum value 2.83 is reached by the gated
  least-rotation model at a conditional efficiency of 2/(1 + √2) = 0.828, the CHSH detection-
  loophole threshold. The loophole-free experiments of 2015 (Hensen et al., event-ready; Giustina
  et al. and Shalm et al., photon pairs with efficiencies above the threshold of their
  inequalities) exclude this route. The model therefore has no way to the measured correlations.

## 3. The cases that cost the model, named

- **C1**, above: Bell's theorem fixes S ≤ 2 at full detection for every local rule.
- **C2**, the only route above 2 is the gate, which is the detection loophole, closed.
- **C3**, "entanglement = equal phases" is, in quantum mechanics, a product state with a shared
  classical variable (2026-09-06 ledger); its correlations are those of Bell's local foil.
- **C4**, what a rule would have to do: act on both photons at once, conditioning the outcome at
  one analyzer on the setting at the other, which no rule local to each analyzer can do
  (MODEL_SECTION.md §6).

## 4. Verification

E_least(θ) and E_Malus(θ) against their closed forms within 3 s.e. at every θ; S values within
3 s.e. of 2.000, 1.414, 2.828, 0; the Garg–Mermin identity S_c = 4/η_c − 2 for R-least at g = 0.6
within 0.02; g = 0.5 leaves R-least unchanged to 10⁻¹².

## 5. What will be recorded

`results/correlation.csv`, `results/chsh.csv`, `results/gated.csv`, `results/verification.txt`,
`results/run_info.txt`, `fig1_pair.png`; `RESULTS.md` with the uncertainty first; `LEDGER.md`.

---

## Addendum, 2026-09-08 (late night): JB's two options for the violation

Written before `pair_options.py` was run. JB (2026-09-08): two ways a violation could occur.
**Option 1**, nonlocality of the wave function: the rotation of the spinor phases in the
entangled particles occurs with the first interaction, so the second particle's spinors are
rotated with the first's. **Option 2**, no nonlocal communication: the first spinors of the bulk
at both detectors are synchronized (the second spinors are not); the entangled particle becomes
unentangled at the first interaction but keeps the same (or opposite) phases; when the second
particle's first spinor rotates to match the second detector's first spinor, which is
synchronized to the same bulk phase as the first detector's, the rotation is equal, and the
coupled rotation of the second spinor (equal to that of the first particle's second spinor)
produces the violation.

- **D7, option 2 (V1, local).** Both bulks share one first-spinor phase μ (fixed for the run, and
  drawn per pair as the sensitivity). Particle k has first-spinor phase a_k and second-spinor
  phase b_k with a shared offset δ = b_k − a_k (opposite offsets as the sensitivity). At contact
  the first spinor rotates to μ, by R_k = μ − a_k, and the second spinor by the same R_k, so its
  phase becomes μ + δ. The outcome is the projection of the second spinor on the local analyzer
  axis: least rotation, + if cos 2(μ + δ − a) ≥ 0; Malus, + with probability cos²(μ + δ − a). The
  rotations R₁ and R₂ are equal when a₁ = a₂ and differ otherwise; either way each outcome is a
  function of (μ, δ, local axis) only.
- **D8, option 1 (V2, nonlocal update).** Particle 1 is read at A as in D3. Its outcome s_A then
  sets particle 2's second-spinor phase to A's axis, a if s_A = +, a + 90° if s_A = −, before
  particle 2 reaches B, where it is read by Malus (primary) or by least rotation (sensitivity).
  This is the projection postulate applied to the pair: the state after A's outcome is
  particle 2 polarized along A's result.
- **D9, statistics.** As D5: E(θ) at 13 angles, CHSH at the standard settings, 4 × 10⁵ pairs.

**Predictions.**

- **P4, option 2 gives exactly the runs already made.** With Λ = μ + δ the outcomes are
  cos²(Λ − a) and cos²(Λ − b), the model of D1 to D3 with φ renamed. E(θ) = 1 − 4θ/π and
  S = 2.000 for least rotation; ½ cos 2θ and 1.414 for Malus; μ fixed or drawn changes nothing;
  opposite offsets flip the sign of E and leave |S| unchanged. The synchronized bulk phase is a
  common cause fixed before the settings are chosen, and a common cause is what Bell's bound
  already allows for. No violation.
- **P5, option 1 gives the experiment or more.** With Malus at B: P(+ at B | s_A) = cos²(a − b) or
  sin²(a − b), so E(θ) = cos 2θ and S = 2.828, the quantum values. With least rotation at B the
  outcome at B is deterministic given s_A: E(θ) = sign(cos 2θ), a step, and S = 4.000, above the
  quantum value. The marginal at B stays ½ in both, so nothing is signalled; the update of
  particle 2 by A's outcome is the nonlocality.

**Verification.** V1 against Tables 1 and 2 within 3 s.e.; V2 with Malus against cos 2θ and 2.828
within 3 s.e.; V2 with least rotation: E(22.5°) = 1.000 exactly, S = 4.000 exactly.

**Recorded.** `results/options.csv`, `results/verification_options.txt`, `fig2_options.png`;
a section in `RESULTS.md`; ledger L20.
