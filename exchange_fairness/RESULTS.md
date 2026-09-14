# Results: fairness in the conserving exchange class

**Date:** 2026-09-14 (night). Definitions and predictions: `PREREGISTRATION.md` (unchanged by
this file; committed ahead of these results in `1407d2a`). Script: `exchange_fairness.py` (seed
20260919, 2 × 10⁴ trials per cell, 1158 s). Every number below is in `results/`.

## Read this before quoting any number

- These are classical stochastic models of an energy pot shared among sites, run to test what
  fairness requires once the pot is closed. Nothing here is a measurement, and nothing here bears
  on quantum mechanics except through what the Dirac–Kuramoto record's Paper 1 claimed about its
  Theorem 1 (R1 M5 of that paper's review round). The initial shares are the intensities, put in.
- **Two registered predictions missed.** P2 said every model-B trial terminates; 0.44 % had not
  reached the absorbing boundary at the time cap of 60 (the √(e_i e_j) kernel's approach to a
  depleted site is slow near zero, and the cap on the transfer keeps a few trials hovering).
  Nothing else in P2 missed: the pot is closed to 2 × 10⁻¹⁴ and the share among decided trials is
  0.799. P4 said the two-site parcel model would be biased below 0.75; it gave 0.798. The reason
  is clean and is the finding: with two sites a depleted site means the other holds everything,
  so the game ends and the reflecting zero never acts. **The reflecting-zero bias is a three-body
  effect**, and the three-site prediction held (0.425 against Born's 0.64). Both misses kept.
- One registered verification therefore FAILs on the wording "every trial terminates"
  (`results/verification.txt`, run_info FAIL); the numerical checks all pass.

## Table 1. Two sites (0.8, 0.2), Born 0.800; s.e. 0.003

| model | termination | P(site 1) | pot closed | note |
|---|---|---|---|---|
| A, independent noise de_i = σ√e_i dW_i | last survivor | 0.801 | **no** | survivor's energy: median 0.72; below ½ in 0.34; above 1.5 in 0.16 |
| B, conserving √(e_i e_j) kernel | absorption at the full share | 0.799 | to 2 × 10⁻¹⁴ | 0.44 % undecided at the cap (P2 miss) |
| C, conserving e_i e_j kernel | rate-linear commit, κ = 1 | 0.801 | to 7 × 10⁻¹⁵ | sites never absorbed |
| C, conserving e_i e_j kernel | threshold t = 0.95 | 0.830 | | closed form (s₀ − (1 − t))/(2t − 1) = 0.833 |
| D, fixed parcels, reflecting zero | absorption | 0.798 | exact | registered below 0.75: **missed**; two sites cannot reflect |

## Table 2. Three sites (0.64, 0.29, 0.07)

| model | P(site 1) | P(site 2) | P(site 3) |
|---|---|---|---|
| Born | 0.640 | 0.290 | 0.070 |
| A, independent noise | 0.638 | 0.293 | 0.067 |
| B, conserving √(e_i e_j) | 0.635 | 0.289 | 0.070 |
| D, parcels, reflecting zero | **0.425** | 0.298 | **0.278** |
| Paper 1 §4(b)'s quoted values for D | 0.46 | 0.27 | 0.26 |

## Predictions against outcome

- **P1**, model A fair but non-conserving: as registered (0.801; median 0.72, 0.34, 0.16 against
  about 0.7, 0.34, 0.15).
- **P2**, model B fair and conserving: fair (0.799) and conserving (10⁻¹⁴); **termination missed**
  by 0.44 % at the cap.
- **P3**, model C: as registered (0.801 under the rate-linear commit; 0.830 against 0.833 under the
  threshold; no trial reached the cap).
- **P4**, model D: three sites as registered (0.425, within ±0.08 of Paper 1's 0.46; bright site
  0.21 below Born); **two sites missed** (0.798, not below 0.75).
- **P5**, the clause: supported by A to D, with the two-site correction.

## What it says, in plain form

1. **In a closed pot, fairness is automatic.** The share is e_i/E with E constant, so its drift is
   the drift of e_i, and any symmetric zero-mean exchange has none. Model B's geometric-mean kernel
   and model C's multiplicative kernel, which under independent noise would drift in opposite
   directions by Paper 1's Theorem 1, are both exactly fair here. The √e scaling is neither
   necessary nor unique in this class.
2. **Theorem 1's model is not this class.** Independent noise on each site lets the pot wander:
   the survivor ends with a median 0.72 of a quantum and holds less than half of it in a third of
   runs. The share is still a martingale, so the site statistics are Born, but the closed-pot
   arguments that Paper 1 rests exclusivity on (P4(a)) do not describe this model.
3. **What spoils fairness is a boundary, and only with three or more sites.** A depleted site that
   can only gain has an upward drift, and with three sites it takes 0.21 from the bright site's
   share. With two sites the boundary is never reached without the game ending, so the same rule
   is fair. The bias is a property of the many-site game.
4. **What spoils Born is the termination.** Stopping at a threshold below the full share hands the
   leader more than its share, by exactly the closed-form amount Paper 1 gave; absorption at the
   full share, or a commit at a rate linear in the share, does not.
5. **The clause for the specification.** If an exchange stage precedes commitment: a closed pot,
   a symmetric zero-mean kernel, an absorbing zero, and a termination that is either absorption at
   the full share or a rate linear in the share. Nothing about the kernel's scaling.
