# Results: the vacuum race, "highest interference term wins" with a per-atom vacuum term

**Date:** 2026-09-14. Definitions and predictions: `PREREGISTRATION.md` (unchanged by this file).
Script: `vacuum_race.py` (seed 20260914, 53 s). Side ledger: `LEDGER.md`, continuing the pair
case's numbering (L24 onward). Every number below is in `results/`.

## Read this before quoting any number

- The rule is JB's of 2026-09-14 (pre-registration §0): among the atoms of both pixels, the one
  with the highest interference-weighted deposit takes the quantum, and each atom carries its own
  vacuum term, one random number per atom per trial, independent of the number of atoms. Nothing
  here is a measurement. The vacuum's geometry is not modeled (C4); only the distribution of the
  per-atom term and the way it enters the deposit are varied.
- **One winner per trial is by construction (D5).** The single record per quantum is assumed in
  this study, not obtained. What the study decides is the pixel share.
- **Every registered verification passed** (`results/verification.txt`, 15 checks).
- **One prediction missed in its wording.** P5 registered a share of 1.000 at N ≥ 100 for all
  three Gaussian and bounded forms. The additive Gaussian form (V4) gave 0.966 at N = 100 and
  0.990 at N = 1000. Cause: the maximum of N Gaussian draws grows as √(2 ln N), so with a noise
  of 2 against an intensity difference of 3 the concentration is slow; the direction and the
  limit are as registered, the rate is not. The multiplicative Gaussian and the bounded uniform
  forms gave 1.000 as registered. The miss is kept.
- The Born values quoted (0.800 at r = 4, 0.667 at r = 2) are the intensity share, put in by hand
  through D3. Nothing here derives them (C1).

## Table 1. The pixel share at r = 4 with the phase weight (Born 0.800)

Source: `results/share.csv`. γ from P(A)/P(B) = (√r)^γ; Born gives 2. s.e. 0.0013 (N ≤ 100),
0.003 (N = 1000).

| N atoms per pixel | 1 | 10 | 100 | 1000 | registered |
|---|---|---|---|---|---|
| V0, no vacuum (R3) | 0.790 | 1.000 | 1.000 | 1.000 | 0.79, 1.000 |
| V1, Gumbel on the deposit | 0.702 | 0.794 | 0.801 | 0.802 | 0.705, 0.795, 0.800, 0.800 |
| V1: exponent γ | 1.24 | 1.95 | 2.01 | 2.02 | R2 row: 1.25, 1.95, 2.00, 2.00 |
| V2, Gumbel on the amplitude | 0.641 | 0.663 | 0.668 | 0.668 | 0.634, then 0.667 |
| V2: exponent γ | 0.84 | 0.98 | 1.01 | 1.01 | 1 |
| V3, multiplicative Gaussian | 0.751 | 0.999 | 1.000 | 1.000 | 0.75, 0.999, 1.000 |
| V4, additive Gaussian | 0.678 | 0.882 | 0.966 | 0.990 | 0.68, 0.88, **1.000** |
| V5, additive uniform | 0.742 | 0.994 | 1.000 | 1.000 | 0.74, 0.99, 1.000 |

Without the phase weight, V1 gives 0.800, 0.801, 0.799, 0.797 at the four N (verification), and
V2 gives 0.667 at N = 1. The second study's R2 row (Table 2 there: 0.704, 0.794, 0.800, 0.800)
is reproduced by V1 within 3 combined s.e. at every N: **the maximum with a Gumbel vacuum term is
the sampled rule R2**, atom by atom.

## Table 2. The pixel share at r = 2 with the phase weight (Born 0.667)

| N | 1 | 10 | 100 | 1000 |
|---|---|---|---|---|
| V0 | 0.692 | 0.999 | 1.000 | 1.000 |
| V1 | 0.607 | 0.661 | 0.666 | 0.669 |
| V2 | 0.570 | 0.586 | 0.585 | 0.588 |
| V3 | 0.646 | 0.960 | 1.000 | 1.000 |
| V4 | 0.567 | 0.647 | 0.723 | 0.784 |
| V5 | 0.602 | 0.852 | 1.000 | 1.000 |

V1 at N = 10: 0.661 against the registered 0.663. V2 tends to √2/(√2 + 1) = 0.586. The
additive Gaussian form crosses Born's value on its way to 1 (0.647, 0.723, 0.784), which is a
coincidence of N and not a law.

## Table 3. The scale of the vacuum term sets the exponent (V1, N = 100, r = 4)

Source: `results/beta.csv`. Registered: P(A) = r^β/(r^β + 1), γ = 2β.

| β | ½ | 1 | 2 |
|---|---|---|---|
| P(A), phase weight | 0.664 | 0.798 | 0.941 |
| P(A), no phase weight | 0.668 | 0.800 | 0.941 |
| predicted | 0.667 | 0.800 | 0.941 |
| γ | 0.98, 1.01 | 1.98, 2.00 | 4.00, 4.01 |

Only β = 1 is Born. The scale of the vacuum term is a free parameter of the rule, and Born's
value is one point of it.

## Table 4. The scaling test: (I_A, I_B) = (4, 1) against (40, 10)

Source: `results/scaling.csv`. Born depends on the ratio only, so the two columns should agree.

| form | phase weight | N | at (4, 1) | at (40, 10) |
|---|---|---|---|---|
| V1 Gumbel | yes | 1 | 0.705 | 0.706 |
| V1 Gumbel | yes | 10 | 0.794 | 0.796 |
| V1 Gumbel | no | 1 | 0.798 | 0.799 |
| V3 mult Gaussian | no | 1 | 0.927 | 0.928 |
| V3 mult Gaussian | no | 10 | 1.000 | 1.000 |
| V4 add Gaussian | no | 1 | 0.857 | 1.000 |
| V4 add Gaussian | yes | 1 | 0.679 | 0.791 |
| V4 add Gaussian | yes | 10 | 0.882 | 1.000 |
| V5 add uniform | no | 1 | 0.968 | 1.000 |
| V5 add uniform | yes | 1 | 0.742 | 0.790 |

The multiplicative forms (V1, V3) are unchanged under scaling; V3 still concentrates with N.
The additive forms read the intensity difference: ten times the light at the same ratio moves
their share to the no-vacuum value (0.790 at N = 1) or to 1.000. All as registered (P6).

## Table 5. Within the pixel: does the best phase-matched atom win?

Source: `results/best_atom.csv`. Fraction of trials in which the winning atom is the one with the
highest phase weight in its pixel; registered 2/N for V1, 1 for V0.

| N | 1 | 10 | 100 | 1000 | 2/N |
|---|---|---|---|---|---|
| V0 | 1.000 | 1.000 | 1.000 | 1.000 | |
| V1 | 1.000 | 0.201 | 0.021 | 0.0023 | 1, 0.2, 0.02, 0.002 |
| V2 | 1.000 | 0.157 | 0.016 | 0.0015 | |
| V3 | 1.000 | 0.353 | 0.052 | 0.0075 | |
| V4 | 1.000 | 0.260 | 0.035 | 0.0040 | |
| V5 | 1.000 | 0.370 | 0.078 | 0.0174 | |

Under V1 the best-matched atom is the single most likely winner and wins in 2/N of trials, as
registered (P7). Under V0 it wins always. JB's "the wavefunction becomes a larger factor within
the aperture" holds for every form in the sense that the winner is drawn toward the best-matched
atom; it is certain only without a vacuum term.

## Predictions against outcome

- **P1**, V0: 0.790, then 1.000. As registered.
- **P2**, V1 is R2: 0.702, 0.794, 0.801, 0.802 against the R2 row, all within 3 s.e.; 0.800 at every
  N without the phase weight; 0.661 at r = 2, N = 10 against 0.663. As registered.
- **P3**, exponent 2β: 0.98, 2.00, 4.00. As registered.
- **P4**, V2: 0.641 at N = 1, 0.663 to 0.668 after, exponent 1. As registered.
- **P5**, V3, V4, V5 concentrate: V3 and V5 as registered (1.000 from N = 100); **V4 missed**,
  0.966 and 0.990 against 1.000, the slow Gaussian maximum. The direction is as registered.
- **P6**, scaling: as registered for all four forms.
- **P7**, within pixel: 0.201, 0.021, 0.0023 against 2/N. As registered.

## What the run says, in plain form

1. **"Highest wins" with a per-atom vacuum term gives Born for exactly one family of vacuum
   distributions, and with it the rule is the second study's R2.** A Gumbel term on the logarithm
   of the deposit, or equivalently an exponential waiting time with rate proportional to the
   deposit, turns the maximum into a proportion. The pixel share is then the intensity share
   times the phase-sum ratio, which is Born once a pixel holds about ten atoms, and it is
   independent of N and of the absolute intensity, as JB stipulated the vacuum term should be.
2. **Every other form tried gives the brighter pixel every time, or nearly.** Gaussian and
   bounded vacuum terms concentrate under the maximum, quickly for multiplicative and bounded
   forms, slowly for an additive Gaussian, and additive forms also fail the scaling test, reading
   the intensity difference where Born reads the ratio. The per-atom term does not depend on N;
   the maximum of N such terms does.
3. **The vacuum has to multiply the intensity.** A vacuum that enters by interfering with the
   particle's amplitude, so that its random part is linear in |ψ|, selects in proportion to |ψ|
   and gives √r/(√r + 1), the square root of Born. That is the interference term's own scaling,
   and it is the wrong one.
4. **The scale of the vacuum term is a free parameter, and Born is one point of it.** The share is
   r^β/(r^β + 1) for a Gumbel term of scale 1/β. β = 1 is Born; nothing in the rule fixes it.
5. **What this buys the model, and what it costs.** It is the first selection dynamics in this
   program that keeps one record per quantum (by construction here), the whole quantum, and the
   Born share at once, without a stop across pixels and without assuming the empty pixel silent.
   It costs the derivation: a rate proportional to |ψ|² is the golden rule, and the square enters
   in the vacuum's distribution. The model can adopt it as a mechanism for the detector given
   Born; it cannot present it as a source of Born (L24).
