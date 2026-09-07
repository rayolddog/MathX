# Results: sink field and channel growth in two dimensions

**Date:** 2026-09-07. Definitions and predictions were fixed in `PREREGISTRATION.md` before the
run; its §6 records an erratum found after the first run, a typo, and two added sensitivities.
Run: `python3 sink_field.py` (about 100 s). Every number below is in `results/`. The first run,
made with the field's sign reversed, is kept in `results/first_run_wrong_sign/` and not quoted.

## Read this before quoting any number

- This is a transport model with no physics labels: a density spread by diffusion, atoms with
  random phases, the second study's deposit race as the seed, a charging threshold at every atom,
  absorption at latched atoms, a screened field sourced by the charges, and a conductance tensor
  that grows with the flow it carries. The choices A1 to A4 and the parameter values of §1 are
  what every number depends on; none of them is a measurement.
- Two parameter sets are reported. The registered set is the one the predictions were written
  for. The addendum set (stronger conductance growth, lower threshold) was added after the first
  run showed that the registered values could not put the channel question or the double-record
  question at r = 4 to the test.
- Most registered predictions missed, and the misses are the result. They are listed under
  "Predictions against outcome" with the reason each one missed.
- Verification: Poisson solve, point source at the centre: max residual of (∇² − 1/ℓ²)Φ + q = 4.3e-16; Φ(48,48) = 0.7176, Φ(64,48) = 0.0681, Φ(48,64) = 0.0681; Diffusion only for t = 60 against the Gaussian solution (σ² → σ² + 2Dt): max relative error 4.92e-03; mass conservation per step (ρ + absorbed), worst over the set: 4.4e-16; mass conservation per step (ρ + absorbed), worst over the set: 4.4e-16. All pass.

## The answer in one paragraph

Once every atom can charge and latch, the record is not one atom but a cluster: sixty to two
hundred latched cells in the lobe that reaches threshold first, which absorb the density where it
lies. Nothing then has to travel, so no flow concentrates, no conductance heats, and the channel
question never arises in this regime. The field, with the sign corrected, raises collection
(0.65 with it against 0.51 without, registered set) but does not stop the other lobe from
latching once that lobe can reach threshold on its own: the double-record rate is 1.00 at every
setting of the addendum set and at every equal-lobe setting of the registered set, including the
longest range. And the charging race overrides the seed: the bright lobe is recorded first in
about 0.95 of trials, not the 0.80 the seed was sampled with, so the site statistics move from
Born toward a step as soon as the record has dynamics of its own.

## Registered set (μ = 200, C_th = 0.01), 20 trials per setting

| setting | seed in bright lobe | bright lobe first | double record | winner collects | total collected | latched cells L / R | heated cells |
|---|---|---|---|---|---|---|---|
| base | 0.80 | 0.95 | 0.05 | 0.649 | 0.683 | 0 / 65 | 0 |
| range 4 | 0.95 | 1.00 | 0.00 | 0.577 | 0.577 | 0 / 62 | 0 |
| range 64 | 0.85 | 0.95 | 0.05 | 0.683 | 0.719 | 0 / 68 | 0 |
| no field | 0.90 | 0.95 | 0.05 | 0.505 | 0.532 | 0 / 66 | 0 |
| scalar M | 0.80 | 0.95 | 0.05 | 0.648 | 0.682 | 0 / 66 | 0 |
| no growth | 0.75 | 0.95 | 0.05 | 0.648 | 0.682 | 0 / 67 | 0 |
| r 1 | 0.60 | 0.30 | 1.00 | 0.343 | 0.685 | 32 / 34 | 0 |
| r 1 range 64 | 0.45 | 0.55 | 1.00 | 0.367 | 0.738 | 35 / 35 | 0 |
| no noise | 0.90 | 1.00 | 0.00 | 0.667 | 0.667 | 0 / 50 | 0 |
| r 1 no noise | 0.75 | 0.65 | 1.00 | 0.317 | 0.637 | 22 / 22 | 0 |

## Addendum set (μ = 2000, C_th = 0.004), 20 trials per setting

| setting | seed in bright lobe | bright lobe first | double record | winner collects | total collected | latched cells L / R | heated cells |
|---|---|---|---|---|---|---|---|
| base | 0.85 | 0.95 | 1.00 | 0.731 | 0.904 | 66 / 220 | 0 |
| range 4 | 0.75 | 0.80 | 1.00 | 0.617 | 0.873 | 65 / 229 | 0 |
| range 64 | 0.85 | 0.90 | 1.00 | 0.710 | 0.917 | 66 / 215 | 0 |
| no field | 0.60 | 0.85 | 1.00 | 0.641 | 0.865 | 63 / 237 | 0 |
| scalar M | 0.70 | 0.95 | 1.00 | 0.730 | 0.903 | 63 / 219 | 0 |
| no growth | 0.90 | 1.00 | 1.00 | 0.760 | 0.902 | 64 / 221 | 0 |
| r 1 | 0.60 | 0.50 | 1.00 | 0.454 | 0.909 | 163 / 164 | 0 |
| r 1 range 64 | 0.50 | 0.45 | 1.00 | 0.461 | 0.923 | 163 / 164 | 0 |
| no noise | 0.85 | 1.00 | 1.00 | 0.734 | 0.823 | 17 / 141 | 0 |
| r 1 no noise | 0.45 | 0.45 | 1.00 | 0.416 | 0.833 | 90 / 89 | 0 |

"Winner collects" is the density absorbed in the lobe recorded first by t = 60; "total" adds the
other lobe's. "Heated cells" are those whose conductance trace passed 1.5 times its initial value:
none did, in any run of either set. Mass is conserved to 4 × 10⁻¹⁶ per step.

## Predictions against outcome

- **H1 (no field gives double records above 0.9): missed in the registered set, held in the
  addendum set.** At C_th = 0.01 the weak lobe cannot reach threshold on its own at r = 4 (its
  charge saturates near 0.009), so the double-record rate was 0.05 with or without the field.
  With C_th = 0.004 it was 1.00 with or without the field.
- **H2 (short range behaves like no field): held** in both sets, in the sense that the range
  changed no rate.
- **H3 (long range brings the double-record rate below 0.5): missed** in both sets. At ℓ = 64
  the rate was 0.05 at r = 4 in the registered set, for the threshold reason above, and 1.00 in
  every case where the far lobe could latch. The field moves density but does not stop a lobe
  from charging.
- **H4 (tensor growth makes a line, scalar a patch): not testable.** No cell reached the heated
  criterion in either set, because the flow never concentrated: the density was absorbed by the
  latched cluster where it lay, and the field's gradient inside a spread cluster of charges is
  small. The comparison of tensor and scalar growth therefore stands untested.
- **H5 (bright lobe first in 0.80 ± 0.10 of trials): missed.** The seed landed in the bright
  lobe in 0.75 to 0.95 of trials, as sampled, but the bright lobe was recorded first in 0.95 to
  1.00, because its peak cells charge to threshold in about two time units whether or not they
  hold the seed. At r = 1 the left-first rate was 0.30 to 0.70 across settings, consistent with
  one half at twenty trials.
- **H6 (noise changes nothing beyond the spread): held** for the rates; the number of latched
  cells fell without noise (50 against 65 at r = 4, registered), since noise pushes marginal
  cells over threshold.
- **H7 (collection above 0.9 with the field at long range): missed in the registered set** (0.72
  at ℓ = 64), **held in the addendum set** (0.92), where the lower threshold latches more cells.

## The consequences named in §3

- **C1, double records.** Whenever the far lobe can reach threshold on its own, it does, at every
  range and with or without growth. The field never prevented one.
- **C2, site statistics.** The seed is Born by construction (0.80); the record is not (0.95). The
  charging threshold turns sampling by deposit into a race won by the brighter lobe.
- **C3, the channel.** None formed. In this regime the energy does not travel.
- **C4, what the field costs.** It raised collection by about a seventh at the registered
  threshold and moved density toward the recorded lobe, without preventing the other lobe's
  record wherever that record was possible.

## What this says about the model

The second study ended with the double record as the cost case and the channel as the proposed
fix. This study shows that, once atoms charge and latch by a local rule, the channel is never
needed and never forms: a cluster of atoms eats the wave in place, and each lobe that can reach
threshold makes its own record. One record per quantum would need either a rule that suppresses
latching everywhere else once one record has begun, which is a nonlocal rule, or a threshold that
only the winning lobe can reach, which reproduces the brighter lobe always winning. Neither gives
Born's statistics from the dynamics; Born enters through the seed and leaves through the race.

## Files

`results/summary_registered.csv`, `trials_registered.csv`, `fig1_fields_registered.png`,
`fig2_outcomes_registered.png`, and the same four with `_addendum`; `verification.txt`,
`run_info.txt`; `first_run_wrong_sign/` (the reversed-sign run, kept, not quoted).
