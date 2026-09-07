# Pre-registration: sink field and channel growth in two dimensions

**Date:** 2026-09-07. Written before any run; `RESULTS.md` cites this file and does not change the
definitions below. Registered at JB's instruction with my defaults, after his statement of the
same day that the model should carry as little physics as possible.
**Exploration:** `MathX/sink_field_2d/`.
**Related:** `../second_spinor_phase_race/`, whose Parts A to E (the phase race and the deposit)
stand unchanged and whose Parts F and G (a Dirac packet with wells and absorbers) this study
replaces with a transport model. Nothing here uses the Dirac equation, an energy sign, or a mass.

## 0. Read this first

**The question.** When the wave's energy is transported across a field of atoms by a flow that
the developing sinks themselves shape, does one record form per quantum, and is the path the
energy takes a channel or a patch?

**The model (JB, 2026-09-07).** The wave is a density to be transported. Each atom carries a
second-spinor phase, random across the bulk, and the deposit race of the second study seeds the
first sink. Every atom charges from the density passing over it, weighted by its phase match,
and latches into an absorbing sink at a threshold. The charged sites make a field; the density
flows along that field with a conductance that grows where flow has passed, and the conductance
is a tensor, so it can grow along the direction of the flow and not across it. That is JB's
"tensor in which the flow varies across the field with the presence or absence of sinks". It is
the structure of the dielectric breakdown model (Niemeyer, Pietronero and Wiesmann, 1984) with
the phase race as its seed.

**Choices fixed here for the run. Audit these first.**
- **A1, free spreading is diffusion.** Any transport that spreads the density gives the same
  arrivals (second study, Part C); diffusion is the simplest. It is a stand-in, not a claim.
- **A2, the field is screened.** The charged sites source a potential through
  (∇² − 1/ℓ²)Φ = −q, so a sink's reach is the range ℓ. The unscreened case is ℓ → ∞.
- **A3, one cell per atom.** The lattice spacing is the unit of length.
- **A4, the seed is the second study's R2.** One atom sampled with probability proportional to
  density × phase weight starts with a charge equal to its weight times the threshold. All other
  atoms charge from zero.

## 1. Definitions (fixed)

- **D1, lattice and density.** A 96 × 96 grid, periodic, spacing 1. Initial density: two
  Gaussian lobes of width σ = 6 centred at x = 32 and x = 64 on the row y = 48, with intensity
  ratio r (bright lobe on the right), normalised to 1. Free spreading ∂ρ/∂t = D ∇²ρ, D = 0.5.
- **D2, phases and weights.** Atom i has δ_i uniform on the circle; the incident δ_p = 0
  (the value is immaterial for a random bulk, second study C2). Weight w_i = cos²((δ_p − δ_i)/2).
- **D3, charging.** dC_i/dt = ρ_i w_i − C_i/τ_c + noise, τ_c = 10. The noise is an
  Ornstein–Uhlenbeck term on C_i with correlation time 1 and strength s (default 10⁻³). Threshold
  C_th = 0.01: the site latches when C_i + η_i first reaches it, and stays latched.
- **D4, absorption.** A latched site removes density at the rate Γ_i = γ₀ [1 + κ (C_i + E_i)],
  γ₀ = 1, κ = 10; E_i is the density it has absorbed. The collected fraction is Σ E_i.
- **D5, the field.** Charge q_i = C_i + E_i. (∇² − 1/ℓ²)Φ = −q solved on the grid by FFT,
  F = −∇Φ; default range ℓ = 16, sweep {4, 16, 64}.
- **D6, the flow and the conductance.** J = −D∇ρ + ρ v with v = M·F, |v| capped at 0.5 per unit
  time. M is a symmetric 2 × 2 tensor per cell, initially M₀ times the identity, M₀ = 10, and it
  grows with the flow that has passed: dM/dt = μ |J| ĵ ĵᵀ (tensor, default) or μ |J| I (scalar
  sensitivity), μ = 200, so that a passed flow of 0.05 doubles M₀ along its direction.
  "No field" means M₀ = μ = 0.
- **D7, records.** A lobe is "recorded" when at least one of its cells (x < 48 for the left
  lobe, x > 48 for the right) has latched. A double record is both lobes latched by t = 60.
- **D8, channel index.** Heated cells are those whose M trace exceeds 3 M₀ (6 M₀ in total). The
  index is the bounding box's length along x over its width along y; a line has index above 3,
  a patch below 2. Also the mean anisotropy (λ₁ − λ₂)/(λ₁ + λ₂) of M over heated cells.
- **D9, run.** Forward Euler with upwind drift, dt = 0.1, 600 steps (t = 60). Twenty trials per
  setting with fresh phases and noise. Settings, one at a time around the base
  (r = 4, ℓ = 16, tensor, μ = 200, s = 10⁻³): ℓ = 4; ℓ = 64; no field; scalar M; μ = 0; r = 1;
  s = 0 at r = 4; s = 0 at r = 1.

## 2. Predictions (fixed before the run)

- **H1, no field.** Both lobes latch in nearly every trial: the double-record rate is above 0.9,
  as in the second study's Part G without a channel. The winning lobe collects less than its own
  share of the quantum by t = 60.
- **H2, short range.** ℓ = 4 behaves like no field: the pull does not reach the other lobe.
- **H3, long range.** At ℓ = 64 the seeded lobe's sinks draw the other lobe's density before its
  own cells latch, and the double-record rate falls below 0.5; ℓ = 16 lies between. This is the
  expectation and the least certain number here.
- **H4, tensor against scalar.** With tensor growth the heated region is a line along the axis
  joining the lobes: channel index above 3 and anisotropy above 0.5. With scalar growth it is a
  patch: index below 2, anisotropy 0 by construction. Without growth there is no heated region.
- **H5, which lobe first.** The bright lobe is seeded with probability r/(1 + r) = 0.8 (R2),
  so the bright lobe is recorded first in 0.8 ± 0.1 of trials at r = 4, and the left lobe in
  0.5 ± 0.1 at r = 1.
- **H6, noise.** Removing the noise changes none of the rates by more than the trial-to-trial
  spread, since the phases already randomise the seed; noise matters only in near ties.
- **H7, collection.** With the field and growth at ℓ = 64 the winning lobe collects more than
  0.9 of the quantum by t = 60; with no field, less than the lobe's own share.

## 3. Consequences to be recorded, including the ones that cost the model

- **C1, double records** at each setting: the count of quanta with two recorded lobes.
- **C2, the record's site statistics** against Born's r/(1 + r), which the seed inherits (second
  study C2) and which the race can only distort.
- **C3, the channel** as line or patch, and whether it forms before or after the second lobe
  latches, since a channel that forms afterwards cannot prevent the double record.
- **C4, what the field costs.** The far lobe's density is moved toward the winner before any
  record forms there. If the far lobe's own record forms anyway, the field has moved density
  without preventing a record; that case is counted, not rounded away.

## 4. What will be recorded

`results/summary.csv` (per setting: rates and means), `results/trials.csv` (per trial),
`results/fig1_fields.png` (density, conductance and latched sites at four times, base setting),
`results/fig2_outcomes.png` (double-record rate, bright-first rate, collection and channel index
by setting), `results/verification.txt` (mass conservation to 10⁻⁶ per step, the Poisson solve
against a point source, the diffusion-only run against the Gaussian solution), `run_info.txt`.

## 5. Method note

Definitions and predictions were fixed in this file before the script was written and run. Any
later change is a dated addendum below, not an edit of §1 to §3.

## 6. Addenda

**2026-09-07, after the first run, before the second.** Three things, recorded in the order they
were found.

1. **Erratum in D5, corrected in the code.** With (∇² − 1/ℓ²)Φ = −q the potential is largest at
   the charges, so −∇Φ points away from them. The words of D5 and §0 say the flow is toward the
   sinks; the formula's sign contradicted them and the code followed the formula, so in the first
   run the sinks repelled the density (collection with the field was lower than without it). The
   flow is now along +∇Φ, which is what was meant. The first run's numbers are kept in
   `results/first_run_wrong_sign/` for the record and are not quoted.
2. **Typo in D8.** "M trace exceeds 3 M₀ (6 M₀ in total)" was self-contradictory; the intended
   criterion is trace > 3 M₀, one and a half times the initial trace 2 M₀. Used as such.
3. **Two sensitivities added, the registered values kept.** The first run showed (a) the flow
   magnitudes ρ|v| are at most a few 10⁻³ per unit time, so at μ = 200 no cell can reach the
   heated criterion within t = 60 whatever the sign; and (b) at C_th = 0.01 the weak lobe at
   r = 4 cannot latch on its own (its charge saturates near 0.009), so the double-record question
   cannot be asked at r = 4 with the registered threshold. Added, as a second parameter set run
   in full: μ = 2000 (a passed flow of 0.005 doubles M₀) and C_th = 0.004 (the weak lobe can latch
   by itself in about ten time units). Predictions H1 to H7 are read against both sets; the
   registered set is the one they were written for, and its misses stand as misses.
