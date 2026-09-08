# Pre-registration: Kuramoto lock within one or two rotations

**Date:** 2026-09-06. Written before any sweep was run; `RESULTS.md` cites this file and
does not change the definitions below.
**Exploration:** `MathX/kuramoto_one_two_rotations/`.
**Related record:** `Projects/Physics/DiracKuramotoFramework` at commit `c9bbf12` (closed).
The two scripts copied from it are in `../dk_origin/` with `PROVENANCE.md`.

## 0. Read this first

**The question (JB, 2026-09-06).** For the Kuramoto/Adler pull of one clock toward a bulk,
record the coupling K and the starting phase difference at which the clock locks within
one or two rotations, and the likelihood of that lock.

**Assumptions stated before the run.**

- **A1, what a rotation is.** One full 2π turn of the clock's own phase, period T = 2π/ω.
  "One or two passes" is read as one or two rotations. Two rotations (4π) is also the
  period of a spinor, so n = 1 and n = 2 are both reported exactly. If JB meant something
  else by "pass", this assumption is the thing to change.
- **A2, the starting angle is uniform.** The clock's phase at contact is taken unrelated
  to the bulk's phase, so the starting difference is uniform on the circle. This is the
  case that costs a lock-based story; a story that pre-aligns the phases would do better.
- **A3, the DK mapping is a posit, not a result of the record.** JB's statement: the
  first (upper) Weyl spinor synchronizes instantaneously; the lower Weyl spinor
  synchronizes as the Kuramoto equation. The closed record says (One World, One Cut v1.5,
  §5.1) that a detector site's phase is slaved from the first instant and that gradual
  Kuramoto/Adler locking exists only for a self-sustaining oscillator. The record's own
  reduction of the Dirac equation gives a cosine, not a sine, in the L–R phase equation:
  no attractor. So this exploration asks a conditional question: *if* the lower component
  locks by the Adler equation, what K and what starting angle let it lock inside one or two
  turns of its own clock. It records the parameters the posit needs. It does not supply a
  mechanism, and it does not reopen the papers.
- **A4, the lock band is a criterion, not physics.** ε = 0.35 rad is the repo's lock
  tolerance (GLOSSARY). Results are also given at ε = 0.10 and ε = π/4 so the reader can
  see how much the criterion moves the numbers.

## 1. The equation (fixed)

The DK repo's form (EQUATIONS.md §2, `tests/kuramoto_sync.py`):

    dφ/dt = ω + K sin(Φ_bulk − φ),      Φ_bulk(t) = ω_bulk t

Phase error θ = φ − Φ_bulk obeys the Adler equation

    dθ/dt = Δ − K sin θ,                Δ = ω − ω_bulk

Stable fixed point θ* = arcsin(Δ/K), unstable θ_u = π − θ*, both requiring |Δ| < K
(the Arnold tongue). Primary case: one clock pulled by a bulk of infinite inertia (the
record's form). The mutual case, two clocks each pulling the other with K₁ and K₂, gives
the same Adler equation for the difference with K_eff = K₁ + K₂; this is checked
numerically, not assumed.

## 2. Definitions (fixed)

- **D1, dimensionless coupling.** x = K/ω, the pull per radian of the clock's own turn.
  Equivalently the pull per rotation K·T = 2πx.
- **D2, dimensionless detuning.** δ = Δ/ω. Only δ ≥ 0 is run; δ < 0 mirrors by θ → −θ.
- **D3, starting angle.** θ₀ = φ(0) − Φ_bulk(0) in (−π, π], measured from the bulk phase,
  not from θ*. For δ = 0 the two coincide.
- **D4, locked.** |θ(t) − θ*| ≤ ε. In a one-dimensional flow with a single stable point
  the band is absorbing as long as it excludes the unstable point, so first entry is
  enough; no dwell is applied (the GLOSSARY dwell was a race criterion and adds nothing
  here).
- **D5, locked within n rotations.** t_ε ≤ nT with n ∈ {1, 2}, t_ε the first entry time.
- **D6, likelihood.** P_n(x, δ, ε) = the fraction of θ₀ uniform on the circle with
  t_ε ≤ nT. Its floor at K = 0 is ε/π (points already inside the band); that floor is the
  criterion's, not physics, and is stated wherever P_n is quoted.
- **D7, boundary angle.** The edges of the set of θ₀ that lock within n rotations.
  For δ = 0 the set is symmetric, [−θ₀,max, +θ₀,max].
- **D8, the Weyl-split toy.** Upper component phase φ_L: set equal to Φ_bulk at t = 0 and
  held there (slaved). Lower component phase φ_R: starts at Φ_bulk(0) + θ₀ and follows the
  equation of §1. Recorded: the L–R relative phase φ_R − φ_L over two rotations, and the
  rotation at which the lower component enters the band. Limitation L1: the closed Dirac
  mass coupling between L and R is switched off; the components are treated as
  independently coupled to the bulk. This is exactly JB's stated split and nothing more.

## 3. Predictions (closed form; to be checked by two independent numerical routes)

For δ = 0 the Adler equation integrates to tan(θ/2) = tan(θ₀/2)·e^(−Kt), so

    t_ε = (1/K) ln[ tan(|θ₀|/2) / tan(ε/2) ]
    θ₀,max(n) = 2 arctan[ tan(ε/2) · e^(2πn·x) ]
    P_n(x) = (2/π) arctan[ tan(ε/2) · e^(2πn·x) ]

Predicted numbers at ε = 0.35 rad (floor 0.111):

| n | x for P = 50 % | 90 % | 95 % | 99 % |
|---|---|---|---|---|
| 1 | 0.276 | 0.569 | 0.680 | 0.937 |
| 2 | 0.138 | 0.285 | 0.340 | 0.468 |

| x = K/ω | 0.05 | 0.1 | 0.2 | 0.3 | 0.5 | 1.0 |
|---|---|---|---|---|---|---|
| θ₀,max (deg), n = 1 | 27 | 37 | 64 | 99 | 152 | 179 |
| θ₀,max (deg), n = 2 | 37 | 64 | 131 | 165 | 179 | 180 |

At ε = 0.10 the one-rotation thresholds roughly double (P = 50 % needs x = 0.48); at
ε = π/4 they halve (x = 0.14). P_n → 1 only as x → ∞; the antipodal start never locks.

For δ > 0: no lock at all unless x > δ; inside the tongue the approach rate is
√(K² − Δ²), so t_ε diverges at the tongue edge and P_n falls below the δ = 0 curve at
every x, most steeply near x = δ. Predicted qualitatively; the numbers come from the run.

Numerical routes: (1) the flow is one-dimensional, so t_ε is a difference of the
antiderivative of 1/(Δ − K sin θ), in closed form; (2) direct RK45 integration with an
event at the band edge on a random sample. They must agree to better than 10⁻⁶ of a
rotation. The mutual two-clock check must match K_eff = K₁ + K₂ to the same tolerance.

## 4. The cases that cost the posit, named so they are not skipped

- **C1, the antipodal start.** For any finite K a starting angle near π takes arbitrarily
  long. P_n < 1 always; the size of the miss at each x is recorded, not rounded away.
- **C2, the record's own coupling ratio.** The record's Stage 2 ran the carrier-to-pull
  ratio ω/K ∈ {8, 16, 32, 64}, that is x = 0.016 to 0.125. Predicted P₁ at ε = 0.35:
  0.123 at x = 1/64 and 0.236 at x = 1/8, against a floor of 0.111; P₂ at x = 1/8: 0.449.
  Under the record's parameters most starting angles do not lock inside one rotation.
  Taking the clock to be the electron's zitterbewegung beat (2mc²/ℏ = 1.55 × 10²¹ rad/s)
  with any detector-scale pull (10⁹ to 10¹⁵ s⁻¹) gives x ≤ 10⁻⁶ and P₁ at the floor to
  seven decimal places. A lock inside one rotation needs K of order 0.3 to 1 ω: a pull
  comparable to the clock rate itself. Whether any slower hand (a Larmor-type
  precession) could supply such a ratio is a question for JB, recorded, not answered here.
- **C3, the closed Dirac coupling has no attractor.** The record's central result is that
  the mass coupling between the Weyl components enters the phase equation as a cosine
  (normal-mode precession, the zitterbewegung beat), not the Adler sine. The Adler form
  for the lower component therefore needs an open, dissipative pull whose K is not m.
  The toy of D8 assumes that pull and does not derive it.

## 5. What will be recorded

`results/`: `boundary_angles.csv` (D7 and D6 for δ = 0 at the three ε and n = 1, 2),
`thresholds.csv` (x needed for P_n = 50, 90, 95, 99 %), `lock_map_delta0.csv` (rotations
to lock over the (x, θ₀) grid at ε = 0.35), `detuned_edges.csv` (edges and P_n for δ = 0,
0.1, 0.25, 0.5), `verification.txt` (closed form against RK45; mutual check),
`weyl_split_lock_times.csv`, three figures and the toy's figure, and `run_info.txt`.
`RESULTS.md` states the uncertainty and the assumptions at the top, then the tables.

## 6. Method note

Definitions and predictions were fixed in this file before the sweep was written to disk
and run. Any later change to a definition is to be made as a dated addendum below, not
by editing §1 to §4.

**Addendum, 2026-09-08 (the clock; recorded on JB's question of 2026-09-08 whether the
zitterbewegung frequency was used in the synchronization runs).** C2 above evaluated the posit
with the electron's zitterbewegung beat as the clock (x ≤ 10⁻⁶). The session of 2026-09-06,
after this study ran, identified the physical clock of the posit as the site's transition
frequency, not the zitterbewegung beat and not the mc² carrier; with any detector-scale pull
that gives x ≈ 10⁻⁵, and "rotations to lock" dissolves into "phase set at contact, amplitude
over 1/Γ". C2 is therefore overstated as written and its zitterbewegung line is superseded;
the table values stand as the record of the run. Noted with it: a lower-component lag θ₀ is a
negative-energy admixture of fraction sin²(θ₀/2) (second study, ledger L4). The
zitterbewegung frequency is used in this folder only as C2's cost case; the record's
`../dk_origin/dirac_arnold_tongue.py` uses it as the clock of two confined Dirac particles, and
the second study's Part C measures it as the beat of the ±E components of a packet (2⟨E⟩).
