# Results: Kuramoto lock within one or two rotations

**Date:** 2026-09-06. Definitions and predictions were fixed in `PREREGISTRATION.md` before
the run and are not restated here; read that file first. Run: `python3 kuramoto_rotations.py`
then `python3 weyl_split_toy.py`; every number below is in `results/`.

## Read this before quoting any number

- These are properties of the Adler equation dθ/dt = Δ − K sin θ, exact to the stated
  tolerance (closed form against RK45: 9 × 10⁻¹³ rotations; edge inversion against a grid
  measure of the locked set: 9 × 10⁻⁶ in likelihood). Nothing here is a measurement, and
  nothing here is a result of the Dirac–Kuramoto record.
- Every likelihood depends on the lock band ε, which is a criterion. At ε = 0.35 rad the
  floor is 0.111 (starting angles already inside the band). The one-rotation half-point
  moves from x = 0.14 (ε = π/4) to 0.28 (ε = 0.35) to 0.48 (ε = 0.10).
- The starting angle is uniform on the circle (assumption A2). A rotation is one 2π turn of
  the clock itself (A1). The two-Weyl-spinor mapping is JB's posit of 2026-09-06, not a
  result of the record (A3); the toy assumes the Adler pull and does not derive it (C3).

## The answer in one paragraph

A clock pulled toward a bulk locks within one of its own turns, for most starting angles,
only when the pull per radian x = K/ω is about 0.3 or more (K·T ≈ 2). Within two turns the
same holds from x ≈ 0.14 (K·T ≈ 0.9). A near-certain lock (99 %) within one turn needs
x ≈ 0.94; within two turns x ≈ 0.47. The starting angle that locks is bounded by
θ₀,max = 2 arctan[tan(ε/2)·e^(2πnx)], so at x = 0.3 a start up to 99° from the bulk phase
locks in one turn and up to 165° in two. The antipodal start never locks in finite time,
so no finite pull makes the lock certain.

## Table 1. Boundary angle and likelihood, Δ = 0, ε = 0.35 rad

Floor P = 0.111. Source: `results/boundary_angles.csv` (also at ε = 0.10 and π/4).

| x = K/ω | K·T | θ₀,max, one rotation | P₁ | θ₀,max, two rotations | P₂ |
|---|---|---|---|---|---|
| 0.02 | 0.13 | 22.7° | 0.126 | 25.6° | 0.142 |
| 0.05 | 0.31 | 27.2° | 0.151 | 36.7° | 0.204 |
| 0.10 | 0.63 | 36.7° | 0.204 | 63.7° | 0.354 |
| 0.15 | 0.94 | 48.8° | 0.271 | 98.7° | 0.548 |
| 0.20 | 1.26 | 63.7° | 0.354 | 130.8° | 0.727 |
| 0.25 | 1.57 | 80.8° | 0.449 | 152.5° | 0.847 |
| 0.30 | 1.88 | 98.7° | 0.548 | 165.1° | 0.917 |
| 0.40 | 2.51 | 130.8° | 0.727 | 175.7° | 0.976 |
| 0.50 | 3.14 | 152.5° | 0.847 | 178.8° | 0.993 |
| 0.60 | 3.77 | 165.1° | 0.917 | 179.7° | 0.998 |
| 0.70 | 4.40 | 172.0° | 0.956 | 179.9° | 0.999 |
| 0.80 | 5.03 | 175.7° | 0.976 | 180.0° | 1.000 |
| 1.00 | 6.28 | 178.8° | 0.993 | 180.0° | 1.000 |
| 1.25 | 7.85 | 179.7° | 0.999 | 180.0° | 1.000 |
| 1.50 | 9.42 | 179.9° | 1.000 | 180.0° | 1.000 |

The two-rotation column at x is the one-rotation column at 2x: doubling the pull and
doubling the time are the same thing in this equation (the closed form depends on n·x only).

## Table 2. Pull needed for a given likelihood (Δ = 0)

Source: `results/thresholds.csv`. Entries are x = K/ω; multiply by 2π for K·T.

| ε (rad) | rotations | P = 50 % | 90 % | 95 % | 99 % |
|---|---|---|---|---|---|
| 0.10 | 1 | 0.477 | 0.770 | 0.881 | 1.138 |
| 0.10 | 2 | 0.238 | 0.385 | 0.441 | 0.569 |
| 0.35 | 1 | 0.276 | 0.569 | 0.680 | 0.937 |
| 0.35 | 2 | 0.138 | 0.285 | 0.340 | 0.468 |
| π/4 | 1 | 0.140 | 0.434 | 0.545 | 0.801 |
| π/4 | 2 | 0.070 | 0.217 | 0.272 | 0.401 |

## Table 3. Detuning, ε = 0.35 rad

Source: `results/detuned_edges.csv`. No lock at all for x ≤ δ (outside the Arnold tongue).

| δ = Δ/ω | rotations | x = 0.15 | 0.30 | 0.50 | 0.70 | 1.00 |
|---|---|---|---|---|---|---|
| 0 | 1 | 0.271 | 0.548 | 0.847 | 0.956 | 0.993 |
| 0 | 2 | 0.548 | 0.917 | 0.993 | 1.000 | 1.000 |
| 0.1 | 1 | 0.223 | 0.528 | 0.843 | 0.955 | 0.993 |
| 0.1 | 2 | 0.447 | 0.908 | 0.993 | 0.999 | 1.000 |
| 0.25 | 1 | 0 | 0.350 | 0.817 | 0.948 | 0.992 |
| 0.25 | 2 | 0 | 0.822 | 0.988 | 0.999 | 1.000 |
| 0.5 | 1 | 0 | 0 | 0 | 0.914 | 0.988 |
| 0.5 | 2 | 0 | 0 | 0 | 0.996 | 1.000 |

With detuning the locked set of starting angles shifts toward the locked phase θ* and
narrows on the far side. At δ = 0.25, x = 0.3: θ* = 56°, and the one-rotation set runs
from −32° to +94° measured from the bulk phase, a width of 126° against 197° at δ = 0.
Once x is a few times δ the detuning costs little: at x = 0.7 the one-rotation likelihood
is 0.956 (δ = 0), 0.948 (δ = 0.25), 0.914 (δ = 0.5).

## Verification (`results/verification.txt`)

| Check | Result |
|---|---|
| Closed form against RK45, 80 samples, δ = 0 and 0.25, x in [0.05, 1.5] | worst difference 8.6 × 10⁻¹³ rotations, pass |
| Edge inversion against a 720 000-angle grid measure, 4 δ × 3 x × 2 n | worst difference in P 8.7 × 10⁻⁶, pass |
| Mutual clocks K₁ = 0.3, K₂ = 0.2, Δω = 0.1 against Adler with K₁ + K₂ | 6 × 10⁻¹⁵ rotations, pass |
| Pre-registered thresholds (§3) recomputed | identical to three decimals |

So the mutual case (two clocks pulling each other) is the one-sided case with the pulls
added: read Table 1 with x = (K₁ + K₂)/ω.

## The Weyl-split toy (definition D8)

Source: `results/weyl_split_lock_times.csv`, `results/fig4_weyl_split_toy.png`. Upper
component locked at t = 0 by construction; the lower component's lock rotation:

| x = K/ω | θ₀ | lock (rotations) | within 1 | within 2 |
|---|---|---|---|---|
| 0.3 | 30° | 0.22 | yes | yes |
| 0.3 | 90° | 0.92 | yes | yes |
| 0.3 | 150° | 1.62 | no | yes |
| 0.3 | 175° | 2.58 | no | no |
| 0.7 | 30° | 0.09 | yes | yes |
| 0.7 | 90° | 0.39 | yes | yes |
| 0.7 | 150° | 0.69 | yes | yes |
| 0.7 | 175° | 1.11 | no | yes |

The RK45 trajectory and the closed form agree to the trajectory's sampling step (0.0005
rotation). The reading is plain: with the upper component slaved at the first instant, the
L–R relative phase is nothing but the lower component's Adler error, so the "one or two
rotations" question for the spinor is the one-clock question of Table 1. The toy adds no
number that Table 1 does not already contain. Limitation L1 stands: the closed Dirac mass
coupling between the components, which the record found to be a cosine with no attractor,
is switched off. Putting it back is the natural next step and is where the posit would be
tested rather than restated.

## The cases that cost the posit, resolved (pre-registration §4)

- **C1, the antipodal start.** The miss 1 − P at ε = 0.35: one rotation, 0.15 at x = 0.5,
  0.0067 at x = 1, 0.00029 at x = 1.5; two rotations, 0.0067 at x = 0.5, 1.3 × 10⁻⁵ at x = 1.
  The lock is never certain; at x ≥ 1 the miss is below one part in a hundred.
- **C2, the record's own pull.** At the record's Stage-2 ratios (x = 1/64 to 1/8) the
  one-rotation likelihood is 0.123 to 0.236 against a floor of 0.111, and the two-rotation
  likelihood 0.135 to 0.449. Most starting angles do not lock in one turn under the
  record's parameters. For the electron's zitterbewegung clock with any detector-scale pull
  (10⁹ to 10¹⁵ s⁻¹), x ≤ 6 × 10⁻⁷ and the likelihood sits at the floor to six decimals: a
  one-turn lock of that clock is out of reach by six orders of magnitude. A one-turn lock
  needs a pull comparable to the clock rate. Whether a slower hand exists whose rate a
  physical pull can match is JB's question; it is recorded, not answered.
- **C3, no attractor in the closed equation.** Unchanged by the run: the Adler form for the
  lower component is assumed, not derived.

## Predictions against outcome

All closed-form predictions of §3 were reproduced by the independent RK45 route. The
detuned predictions (no lock for x ≤ δ; likelihood below the δ = 0 curve at every x,
steeply near the tongue edge; locked set shifted toward θ*) are confirmed by Table 3.
Nothing in the run contradicted the pre-registration, which is expected for an equation
with a closed-form solution; the run's value is the check and the detuned numbers.

## Files

`results/boundary_angles.csv`, `thresholds.csv`, `lock_map_delta0.csv` (x step 0.01,
θ₀ step 1°), `detuned_edges.csv` (x step 0.01), `verification.txt`,
`weyl_split_lock_times.csv`, `run_info.txt`; figures `fig1_lock_map.png` (rotations to
lock over x and θ₀, with the one- and two-rotation boundaries and the record's Stage-2
range marked), `fig2_likelihood.png` (P₁ and P₂ against x at three ε),
`fig3_detuned.png` (P₁ and P₂ at four δ), `fig4_weyl_split_toy.png`.
