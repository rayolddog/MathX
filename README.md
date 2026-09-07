# MathX

Mathematical explorations that sit beside the Dirac–Kuramoto record
(`Projects/Physics/DiracKuramotoFramework`, closed at commit `c9bbf12`). The record is not
reopened by anything here. Each exploration is one folder with a `PREREGISTRATION.md`
written before any run, the scripts, a `results/` folder, and a `RESULTS.md` that states
the uncertainty and the assumptions before the tables.

| Folder | Question | Status |
|---|---|---|
| `kuramoto_one_two_rotations/` | For the Kuramoto/Adler pull of a clock toward a bulk, which coupling K and which starting phase difference lock within one or two turns of the clock, and how likely is that lock over a uniform starting angle | Run 2026-09-06; see its `RESULTS.md` |
| `second_spinor_phase_race/` | JB's two-spinor model with no pull on the second spinor: first spinors matched at contact, each particle's internal offset carried through, a bulk with random second-spinor phases, and the deposited energy highest where the second spinors are closest in phase; free flight by the Bohm–Dirac velocity | Run 2026-09-07; see its `RESULTS.md` and `LEDGER.md` |
| `sink_field_2d/` | A transport model with no physics labels: density spread by diffusion, atoms with random phases, the deposit race as the seed, charging and latching at every atom, a screened field from the charges, a conductance tensor that grows with flow. Does one record form per quantum, and is the path a channel or a patch | Run 2026-09-07; see its `RESULTS.md` (most registered predictions missed; the misses are the result) |
| `dk_origin/` | The two record scripts this work reuses, copied unchanged, with `PROVENANCE.md` | Copies, not run here |

## Relation to the Dirac–Kuramoto framework

JB's posit of 2026-09-06: compared with the framework, the first (upper) Weyl spinor
synchronizes instantaneously, and the lower Weyl spinor synchronizes as the Kuramoto
equation. The record's closing statement is that a site's phase is slaved from the first
instant and that gradual Kuramoto/Adler locking exists only for a self-sustaining
oscillator; its reduction of the Dirac equation gives a cosine in the L–R phase equation,
so no attractor in the closed system. The exploration therefore asks a conditional
question: if the lower component locks by the Adler equation, what pull and what starting
angle let it lock inside one or two turns of its own clock. It records what the posit
needs; it does not supply the mechanism.

The short answer: a lock inside one turn, for most starting angles, needs a pull per radian
K/ω of about 0.3 or more, and a near-certain lock about 0.9; inside two turns, half of
those. The record's own Stage-2 ratios (ω/K = 8 to 64) sit below that, and the electron's
zitterbewegung clock is out of reach by six orders of magnitude with any detector-scale
pull. Whether a slower hand exists whose rate a physical pull can match is the open
question.

## Running

    cd kuramoto_one_two_rotations
    python3 kuramoto_rotations.py     # tables, three figures, verification; about 2 s
    python3 weyl_split_toy.py         # the two-component toy; imports the first script

    cd ../second_spinor_phase_race
    python3 phase_race.py             # Parts A, B and the ordered bulk; about 15 s
    python3 bohm_dirac_flight.py      # free flight, arrival, local race, zitter; about 10 s
    python3 interference_visibility.py # Part D, the interference term with rule R2; about 8 s
    python3 energy_cost.py            # Part E, the energy cost of the shift; about 10 s
    python3 sink_flow.py              # Part F, the winning atom as a sink; about 10 s
    python3 sink_race.py              # Part G, two sinks racing with noise; about 80 s
    python3 splitter_asymmetry.py     # Part H, uneven split from a spinor asymmetry; about 100 s
    python3 splitter_asymmetry.py slow # Part H addendum, slow charging; about 140 s

    cd ../sink_field_2d
    python3 sink_field.py             # the 2-D transport model, two parameter sets; about 100 s

Python 3.13 with numpy, scipy, matplotlib. Outputs go to each folder's `results/`.
