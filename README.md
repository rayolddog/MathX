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
| `spad_surface_contact/` | Two silicon SPAD pixels: the packet carried by its own Dirac current to the surface, contact with the bound electron of least rotation among the atoms in reach, the record as an induced current down a lattice of random resistances (greedy least resistance against a dielectric-breakdown split and diffusion), a Galton–Watson avalanche. Which pixel, how much rotation, whether a pair forms, what the record looks like; the field reading beside it as the cost case. Addendum the same day: the empty pixel given the same contact and current (it does not stay silent), the SPAD's field as a lateral penalty on the channel, and ordered, locked and periodic phase distributions for the bulk. Addendum 2 (Part D): the bulk as a driven transition coherence dephasing by T2 and T2*, in MRI language; the phase trace lives one T2. `L13_TEST.md`: that trace tested against quantum mechanics for a driven absorber (coherent control decides against the deposit rule as written; ledger L14). Part E: a Dirac electron bound in a well, its second spinor measured in both readings (slaved, in quadrature, identical across atoms; negative-energy content 10⁻⁶ to 10⁻¹⁶), and the deposit weight with unequal amplitudes (visibility 2ab/(a² + b²); ledger L15). Part F: the bulk rebuilt on that slaved spinor, ±90° by side of contact: the phase race has nothing to select, a running wave is in quadrature with every atom, Born by the arrival; a threshold survives only above a gap ratio of ½ with comparable amplitudes (ledger L16). Part G: release and retry at the excited electron, the energy tracked to closure: with divisible shares, one record per quantum by energy conservation above θ = ½ but a step share and sub-quantum records (L21); with the indivisible quantum, everything at once and nothing derived (L22); with re-emitted photons leaving, as they physically do, the divisible reading registers nothing at a balanced split (L23) | Run 2026-09-08; see its `RESULTS.md` and `LEDGER.md` |
| `pair_case/` | JB's model at two analyzers: two photons sharing a hidden polarization, read by his projection rule (least rotation) or by Malus's law at each polarizer; the quantum comparator; the model's phase-dependent efficiency as a detection gate. What correlation function, what CHSH value, and the one route above Bell's bound | Run 2026-09-08; see its `RESULTS.md` and `LEDGER.md` (Bell's sawtooth, CHSH 2.00; the gate rides the Garg–Mermin line, the detection loophole). Addendum: JB's two options, a synchronized bulk phase shared by both detectors (a common cause, S = 2.00) against particle 2 set by A's outcome (the projection postulate, S = 2.83) |
| `vacuum_race/` | JB's rule of 2026-09-14: among the atoms of both pixels the one with the highest interference-weighted deposit takes the quantum, each atom carrying its own vacuum term, one random number per atom independent of the number of atoms. Which distribution of the vacuum term gives Born between the pixels, whether the term must multiply the intensity or may add to the amplitude, and whether the best phase-matched atom wins within a pixel | Run 2026-09-14; see its `RESULTS.md` and `LEDGER.md` (only a Gumbel term on the deposit gives Born, and with it the rule is R2, the golden rule as a race; Gaussian and bounded terms give the brighter pixel every time; a term on the amplitude gives the square root of Born) |
| `dk_origin/` | The two record scripts this work reuses, copied unchanged, with `PROVENANCE.md` | Copies, not run here |

`MODEL_SECTION.md` (2026-09-07) is a draft of the model section a rewritten Born Selection paper
would need: the mechanism as computed here, with Born's rule assumed and the two open items named
first. It is a record, not part of the closed DK repository.

`WICK_ROTATION.md` (2026-09-08) is a note stored for JB to revisit: his 4D-to-3D+spacetime tangent
and the Wick rotation, the signature change that makes the light cone. A note, not a study.

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
    python3 sink_field.py seeded      # seeded charging: one, one per lobe, five, all; about 45 s

    cd ../spad_surface_contact
    python3 spad_contact.py           # two SPAD pixels, least rotation, induced current; about 3 s
    python3 empty_pixel.py            # addendum: the empty pixel with the same contact and current, the field, the bulk's phases; about 5 s
    python3 transition_coherence.py   # addendum 2, Part D: the bulk as a transition coherence with T2 and T2* dephasing; about 3 s
    python3 l13_qm_comparison.py      # the test of L13: quantum mechanics for the driven absorber against the model's deposit rule; about 2 s
    python3 bound_spinor.py           # addendum 3, Part E: the bound electron's second spinor (amplitude, phase, rotation) and the weight with unequal amplitudes; about 20 s
    python3 slaved_bulk.py            # addendum 4, Part F: the bulk rebuilt on the slaved second spinor (±90° by side); about 5 s
    python3 release_retry.py          # addendum 5, Part G: release and retry at the excited electron (ejection, re-emission, phonon loss), readings A and B; about 20 s
    python3 release_retry.py escape   # addendum 5b: re-emitted photons leave instead of retrying; about 10 s

    cd ../pair_case
    python3 pair_case.py              # the pair case: shared hidden polarization, two analyzers, CHSH, gated detection; about 2 s
    python3 pair_options.py           # addendum: JB's two options, the synchronized bulk (local, S = 2) and the update by A's outcome (S = 2.83); about 3 s

    cd ../vacuum_race
    python3 vacuum_race.py            # highest deposit wins with six forms of the per-atom vacuum term; share, exponent, scaling, best atom; about 55 s

Python 3.13 with numpy, scipy, matplotlib. Outputs go to each folder's `results/`.

## The session of 2026-09-08, in plain language

JB asked for a microscopic simulation of one particle meeting a detector, built on his
two-spinor picture: every particle carries two rotating phases, the first locks to the
detector's electron at contact and the second turns with it, and the electron whose second
phase is closest to the particle's takes the energy. The detector was a pair of silicon
single-photon avalanche diodes, simplified to a row of bound electrons on the surface, with
the record made by a current running down the device along the path of least resistance,
as lightning does. Everything was written down before it was run, with the numbers it was
expected to give, and every miss is kept.

**What came out, in order.**

- *The main run.* One record per quantum and the Born share (the brighter pixel getting its
  share of the counts in proportion to intensity) appeared exactly, but both came from the
  particle following its guided path to one pixel, not from the phase mechanism. The mechanism
  only decided which electron, and by how much rotation. Run the same rule on the wave instead
  of the particle and Born vanishes.
- *The empty pixel.* Give the pixel without the particle the same rules and it fights back:
  whenever its share of the wave is above the gap it forms its own record, and the only cures
  are an instantaneous stop across the two pixels or a head start for the particle's pixel. My
  own prediction that the empty pixel never pairs was wrong, and that error is on record.
- *A driven detector.* If the detector's electrons are given a shared oscillation, as an MRI
  pulse gives spins one, the phase mechanism leaves a trace for about one coherence time, ten to
  a hundred femtoseconds in silicon, and then the detector is the random one again. A passive
  detector has no such oscillation at all.
- *Testing that trace.* Quantum mechanics also predicts a phase-dependent response to coherent
  light in a driven absorber, and experiments since 1995 confirm it, so my ledger claim that it
  does not was corrected. Where the two theories differ, in the sign and strength of the effect,
  the existing experiments side with quantum mechanics.
- *The bound electron's second spinor.* Computed directly from the Dirac equation for an
  electron in a well: it is not free. It is fixed at ninety degrees to the first spinor, the
  same in every atom, and small. Rebuilding the detector on that fact leaves the phase race
  nothing to choose between.
- *Release and retry.* Letting the excited electron either eject, re-emit, or lose its energy
  to heat, with the energy of the quantum tracked to the last digit, gives one record per
  quantum by energy conservation alone, but only if the quantum can be split between pixels,
  and then the brighter pixel takes every count and each record carries less than a whole
  photon, both against experiment. Once a re-emitted photon is allowed to leave, as JB pointed
  out it must, the split reading detects nothing at a balanced beam splitter. The reading in
  which the quantum is never split gives everything right and derives nothing.
- *Two particles.* JB's model of an entangled pair gives Bell's own classical correlation and
  exactly the classical bound. It rises above that bound only by letting detection depend on
  the hidden phase, which is the detection loophole that experiments closed in 2015. His idea
  of a common bulk phase shared by both detectors changes nothing, because a shared phase set
  in advance is a common cause and the bound already allows for common causes. The accepted
  picture, the second particle's state fixed by the first result, reproduces the experiment.

**The three answers.** Born's rule is reproduced wherever it is put in, through the square of an
amplitude in the deposit or through the guided path, and it is never derived; every dynamics
laid on top pushes the statistics toward the brighter site. The model's hidden phase is either
fixed by the Dirac equation, set by the particle's direction of motion, or simply the photon's
polarization, and as a hidden variable for pairs it is excluded by the Bell experiments. The
model violates detector independence only through its phase gate, the rule that no pair forms
unless the two second spinors interfere constructively enough; that gate is the closed loophole,
and in silicon it never closes, because a bound electron's second spinor is too small to cancel
the deposit.

**Also from today.** The zitterbewegung frequency was never used as the clock of the phase
race, only tested and ruled out as one; all the Dirac code uses the standard representation,
not the Weyl one; and JB's tangent on the light cone is stored as `WICK_ROTATION.md`. Section 10
of `MODEL_SECTION.md` has the same conclusions with the numbers. Several of my registered
predictions failed during the day, among them the empty pixel never pairing, the rate of ties
between two channels, the sign of a mirrored correlation, and one gap threshold, and three
registered checks failed on their own wording or on a pixel-edge effect. Every one is kept in the
verification files and named in the results.
