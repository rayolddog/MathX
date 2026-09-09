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
