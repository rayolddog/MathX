# Pre-registration: fairness in the conserving exchange class

**Date:** 2026-09-14 (night). Written before any run, for the specification paper
(`born_selection_revised/`), after the Born Selection paper's internal review round found (R1 M5)
that Theorem 1 of that paper analyses an independent-noise model, de_i = σ√e_i dW_i, that does not
conserve the pot: the survivor's energy at absorption is not one quantum. The question here is
what fairness requires in the class that *does* conserve the pot, pairwise exchange between sites,
and whether the √e scaling that Theorem 1 calls the unique fair point is needed there.
**Exploration:** `MathX/exchange_fairness/`. Predictions from the review round's scratch numbers
(R1) and from Paper 1's own closed form for stopping at a threshold. Every number in `results/`.

## 0. Read this first

A share s_i = e_i/E is a martingale, and the site's probability of ending with the whole pot equals
its initial share, whenever the increments of e_i have zero mean and E does not change (the pot is
closed), provided the process terminates with one site holding everything. In the conserving
pairwise class E is constant by construction, so the drift of s_i is the drift of e_i, and any
symmetric zero-mean exchange kernel is fair. Theorem 1's drift computation, E[ds_i] ≠ 0 for
additive or multiplicative noise, comes from E fluctuating under independent noise; it is a
statement about the non-conserving class. What can still spoil fairness in the conserving class is
a boundary: a depleted site that can only gain (a reflecting zero) has an upward drift. And what can
still spoil Born is the termination: stopping at a threshold below the full share biases the
winner. This study runs the four cases.

## 1. Definitions (fixed)

- **D1, the sites.** Two sites with initial energies (0.8, 0.2) (Born 0.800, 0.200); three sites
  with (0.64, 0.29, 0.07) as Paper 1 §4(b)'s configuration. E(0) = 1.
- **D2, model A, independent noise (Paper 1 Theorem 1).** de_i = σ√e_i dW_i with independent
  Wiener processes; Euler–Maruyama, σ = 0.3, dt = 2 × 10⁻³ (R1's parameters); a site whose energy
  would go negative is set to zero and frozen (absorbing); the game ends when at most one site is
  above zero (last survivor). Recorded: the winner, and the survivor's energy E_surv at that time.
- **D3, model B, conserving exchange, geometric-mean kernel.** For every unordered pair (i, j) at
  each step, a transfer δ_ij = σ√(e_i e_j) √dt ξ_ij with ξ_ij standard normal, added to e_i and
  subtracted from e_j; a transfer larger than the giver's energy is capped at it (the cap keeps the
  pot closed and zero absorbing). σ = 0.3, dt = 2 × 10⁻³. Ends when one site holds at least
  1 − 10⁻⁶. Recorded: the winner; E at the end (conservation check).
- **D4, model C, conserving exchange, multiplicative kernel.** As D3 with δ_ij = σ_c e_i e_j √dt ξ_ij,
  σ_c = 1.2 (chosen so that the exchange time is comparable to model B's). This kernel's diffusion
  coefficient vanishes quadratically at zero, so a site is never absorbed and the game never ends
  by itself. Two terminations: (C1) rate-linear commit, each site commits with probability κ s_i dt
  per step, κ = 1, the first commit winning (Paper 1's Theorem 4 structure); (C2) threshold
  commit, the first site to reach s_i ≥ t = 0.95 wins, with argmax at a cap of T = 60 time units
  if none does. Recorded: the winner; the fraction reaching the cap.
- **D5, model D, fixed parcels with clipping (Paper 1 §4(b)).** At each step a random pair and a
  random direction; a parcel δ = 0.02 moves from giver to receiver, capped at the giver's energy
  (so zero is reflecting: a depleted site can only gain). Ends when one site holds at least
  1 − 10⁻⁹. Recorded: the winner.
- **D6, statistics.** 2 × 10⁴ trials per cell (s.e. 0.0028 on a probability near 0.8);
  seed 20260919.

## 2. Predictions

- **P1, model A is fair but does not conserve.** P(site 1) = 0.800 within 3 s.e. The survivor's
  energy: median about 0.7, below 0.5 in about a third of runs and above 1.5 in about a sixth
  (R1: 0.71, 34 %, 15 %; tolerance ±0.08 on each). The closed-pot language of Paper 1 P4(a) does
  not describe this model.
- **P2, model B is fair and conserves.** P(site 1) = 0.800 within 3 s.e.; E = 1 to 10⁻⁹ at every
  step; every trial terminates.
- **P3, model C is fair under a rate-linear commit and biased under a threshold.** C1: 0.800
  within 3 s.e. C2: 0.833 ± 0.010, Paper 1's closed form (s₀ − (1 − t))/(2t − 1) for a continuous
  martingale stopped at t or 1 − t, which holds for any kernel that reaches the boundaries;
  fewer than 1 % of trials reach the cap.
- **P4, model D is unfair, toward the poor.** Three sites: P(site 1) below Born by more than 0.10;
  Paper 1's quoted (0.46, 0.27, 0.26) within ±0.08 each. Two sites: P(site 1) below 0.800 by more
  than 0.05.
- **P5, the clause.** In the conserving class, fairness needs a symmetric zero-mean kernel and an
  absorbing zero, and nothing else; the √e scaling is neither necessary nor unique there. Born
  then needs a termination that does not bias the winner: absorption at the full share, or a
  commit rate linear in the share. This is what a specification should say in place of
  Theorem 1's "unique fair point".

## 3. Verification

Model B: P(site 1) within 3 s.e. of 0.800 and max |E − 1| < 10⁻⁹. Model A: P(site 1) within 3 s.e.
of 0.800. Model C1: within 3 s.e. of 0.800. Model D, three sites: P(site 1) < 0.54.

## 4. What will be recorded

`results/fairness.csv`, `results/survivor_energy.csv`, `results/verification.txt`,
`results/run_info.txt`, `results/fig1_fairness.png`; `RESULTS.md`.
