# Results: the pair case, JB's model at two analyzers

**Date:** 2026-09-08 (night). Definitions and predictions: `PREREGISTRATION.md` (unchanged by this
file). Script: `pair_case.py` (seed 20260912, 4 × 10⁵ pairs per setting, 2 s). Ledger: `LEDGER.md`,
L17 onward.

## Read this before quoting any number

- Bell's theorem fixed the main result before the run, and the pre-registration says so (C1).
  With a shared hidden phase and a rule local to each analyzer, the CHSH value at full detection
  is at most 2 for every rule. The run confirms it and records the model's own correlation
  function, which is Bell's sawtooth for JB's projection rule and a half-amplitude cosine for the
  sampled rule.
- One registered verification failed on its own wording, not on the physics: the check that a
  gate at 0.5 "leaves least rotation unchanged to 10⁻¹²" compared two independent Monte Carlo draws
  of S, which differ by one standard error (1.998 against 2.001). The single-side efficiency was
  1.000000 exactly, which is what the check meant. A labelled diagnostic added after the run gates
  and ungates the same draw: on one draw, gate 0.5 detects 1.000000 of pairs and changes E(a, b) by 0.0e+00. The FAIL stays in `verification.txt` and
  `run_info.txt`.
- One prediction was loosely stated: at gate 0.63 I wrote "about 2.83 at η_c about 0.83"; the gate
  binds a little harder (S_c = 3.01 at η_c = 0.80). The Garg–Mermin identity, which is the exact
  statement, holds to 0.001 at every gate. The quantum value 2.83 is reached at η_c = 0.828.
- Every other number held to within 3 standard errors (E to 0.0016, S to 0.003).

## The answer in one paragraph

Two photons sharing a hidden polarization angle and read by JB's projection rule at each polarizer
give the correlation E(θ) = 1 − 4θ/π, the sawtooth of Bell's own local model, meeting the quantum
cosine at 0°, 45° and 90° and falling short between: 0.50 against 0.71 at 22.5°. Its CHSH value is
2.002 ± 0.003, Bell's bound met exactly. The sampled rule gives ½ cos 2θ and 1.414. Quantum
mechanics gives cos 2θ and 2.827. Independent angles give 0. The model rises above 2 in one way
only: when its own phase-dependent efficiency (SPAD study, Part D) gates detection, the
correlations among coincidences follow the Garg–Mermin bound 4/η_c − 2 exactly, reaching the
quantum 2.83 at a conditional efficiency of 0.828 and the algebraic 4 at 0.625. That is the
detection loophole, and the loophole-free experiments of 2015 exclude it. The model as it stands has
no route to the measured correlations, and a rule that could reach them would have to condition the
outcome at one analyzer on the setting at the other.

## Table 1. The correlation function (`correlation.csv`), 4 × 10⁵ pairs per point, s.e. 0.0016

| θ = a − b | least rotation (JB's projection) | closed form 1 − 4θ/π | Malus, sampled | closed form ½ cos 2θ | quantum | cos 2θ | independent angles |
|---|---|---|---|---|---|---|---|
| 0° | 1.000 | 1.000 | 0.501 | 0.500 | 1.000 | 1.000 | 0.001 |
| 15° | 0.667 | 0.667 | 0.434 | 0.433 | 0.866 | 0.866 | 0.001 |
| 22.5° | 0.501 | 0.500 | 0.353 | 0.354 | 0.708 | 0.707 | −0.001 |
| 30° | 0.332 | 0.333 | 0.249 | 0.250 | 0.499 | 0.500 | 0.001 |
| 45° | −0.002 | 0 | −0.002 | 0 | 0.000 | 0 | −0.001 |
| 60° | −0.336 | −0.333 | −0.249 | −0.250 | −0.499 | −0.500 | 0.001 |
| 67.5° | −0.499 | −0.500 | −0.355 | −0.354 | −0.706 | −0.707 | −0.001 |
| 90° | −1.000 | −1.000 | −0.498 | −0.500 | −1.000 | −1.000 | −0.001 |

## Table 2. CHSH at (0°, 45°, 22.5°, 67.5°), full detection (`chsh.csv`), s.e. 0.003

| rule | S | predicted |
|---|---|---|
| least rotation (JB's projection) | 2.002 | 2.000, Bell's bound |
| Malus, sampled (R2) | 1.414 | 1.414 |
| quantum mechanics | 2.827 | 2.828 |
| least rotation, independent angles | 0.001 | 0 |

## Table 3. Gated detection, the model's phase-dependent efficiency (`gated.csv`)

| rule | gate g | single-side efficiency η₁ | conditional efficiency η_c | S among coincidences | Garg–Mermin 4/η_c − 2 |
|---|---|---|---|---|---|
| least | 0.5 | 1.000 | 1.000 | 1.998 | 2.000 |
| least | 0.6 | 0.872 | 0.853 | 2.689 | 2.690 |
| least | 0.63 | 0.833 | 0.799 | 3.009 | 3.009 |
| least | 0.75 | 0.667 | 0.625 | 4.000 | 4.000 |
| least | 0.85 | 0.505 | 0.506 | 4.000 | 4.000 |
| least | 0.9 | 0.409 | 0.390 | 4.000 | 4.000 |
| Malus | 0.5 | 0.818 | 0.817 | 2.623 | 2.894 |
| Malus | 0.6 | 0.748 | 0.734 | 3.166 | 3.450 |
| Malus | 0.63 | 0.723 | 0.702 | 3.379 | 3.697 |
| Malus | 0.75 | 0.608 | 0.588 | 4.000 | 4.000 |

The least-rotation gate saturates the Garg–Mermin bound at every gate that binds; the Malus gate
stays below it. The quantum value among coincidences is reached by the least-rotation gate at
η_c = 2/(1 + √2) = 0.828, between g = 0.6 and 0.63.

## Verification (`results/verification.txt`)

| Check | Result |
|---|---|
| E(θ) against closed forms, all four rules, 13 angles | all within 3 s.e., pass |
| S for least, Malus, quantum, independent | 2.002, 1.414, 2.827, 0.001 against 2, 1.414, 2.828, 0, pass |
| gate 0.5 leaves least unchanged to 10⁻¹² | **FAIL** as worded (two independent draws, 1.998 against 2.001; η₁ = 1.000000); the same-draw diagnostic shows the gate changes nothing |
| Garg–Mermin identity at g = 0.6 | 2.6893 against 2.6898, pass |

## The cases named in §3, resolved

- **C1 held:** S ≤ 2 at full detection for every local rule; JB's rule meets the bound exactly.
- **C2 held:** above 2 only by gating, on the Garg–Mermin line; closed by the 2015 experiments.
- **C3:** the model's "entanglement = equal phases" gives exactly the correlations of a shared
  classical variable, Bell's local foil. Ledger L17.
- **C4:** the missing rule is one that acts across the pair. Ledger L18.

## The three steps, for the pair

| Step | State |
|---|---|
| 1, free flight | the shared phase carried unchanged (second study, Part C) |
| 2, contact | the local projection at each analyzer, computed here |
| 3, the record | one-photon absorption at each detector, local (SPAD study) |
| the correlations | Bell's sawtooth; the quantum cosine is not reachable by any local rule |

What the model would need is stated in `MODEL_SECTION.md` §6 and is now on record as a number:
the gap between 2.00 and 2.83 at the standard settings, and the one thing that closes it in the
model's own terms, a phase-gated efficiency of 0.83, is what the loophole-free experiments were
built to exclude.
