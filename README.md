# MathX

Mathematical explorations that sit beside the Dirac–Kuramoto record
(`Projects/Physics/DiracKuramotoFramework`, closed at commit `c9bbf12`). The record itself
is not edited by anything here; it was reopened on 2026-09-14 by Paper 1 v0.10's Addendum A, at
JB's instruction, which cites these folders. Each exploration is one folder with a `PREREGISTRATION.md`
written before any run, the scripts, a `results/` folder, and a `RESULTS.md` that states
the uncertainty and the assumptions before the tables.

| Folder | Question | Status |
|---|---|---|
| `kuramoto_one_two_rotations/` | For the Kuramoto/Adler pull of a clock toward a bulk, which coupling K and which starting phase difference lock within one or two turns of the clock, and how likely is that lock over a uniform starting angle | Run 2026-09-06; see its `RESULTS.md` |
| `second_spinor_phase_race/` | JB's two-spinor model with no pull on the second spinor: first spinors matched at contact, each particle's internal offset carried through, a bulk with random second-spinor phases, and the deposited energy highest where the second spinors are closest in phase; free flight by the Bohm–Dirac velocity | Run 2026-09-07; see its `RESULTS.md` and `LEDGER.md` |
| `sink_field_2d/` | A transport model with no physics labels: density spread by diffusion, atoms with random phases, the deposit race as the seed, charging and latching at every atom, a screened field from the charges, a conductance tensor that grows with flow. Does one record form per quantum, and is the path a channel or a patch | Run 2026-09-07; see its `RESULTS.md` (most registered predictions missed; the misses are the result) |
| `spad_surface_contact/` | Two silicon SPAD pixels: the packet carried by its own Dirac current to the surface, contact with the bound electron of least rotation among the atoms in reach, the record as an induced current down a lattice of random resistances (greedy least resistance against a dielectric-breakdown split and diffusion), a Galton–Watson avalanche. Which pixel, how much rotation, whether a pair forms, what the record looks like; the field reading beside it as the cost case. Addendum the same day: the empty pixel given the same contact and current (it does not stay silent), the SPAD's field as a lateral penalty on the channel, and ordered, locked and periodic phase distributions for the bulk. Addendum 2 (Part D): the bulk as a driven transition coherence dephasing by T2 and T2*, in MRI language; the phase trace lives one T2. `L13_TEST.md`: that trace tested against quantum mechanics for a driven absorber (coherent control decides against the deposit rule as written; ledger L14). Part E: a Dirac electron bound in a well, its second spinor measured in both readings (slaved, in quadrature, identical across atoms; negative-energy content 10⁻⁶ to 10⁻¹⁶), and the deposit weight with unequal amplitudes (visibility 2ab/(a² + b²); ledger L15). Part F: the bulk rebuilt on that slaved spinor, ±90° by side of contact: the phase race has nothing to select, a running wave is in quadrature with every atom, Born by the arrival; a threshold survives only above a gap ratio of ½ with comparable amplitudes (ledger L16). Part G: release and retry at the excited electron, the energy tracked to closure: with divisible shares, one record per quantum by energy conservation above θ = ½ but a step share and sub-quantum records (L21); with the indivisible quantum, everything at once and nothing derived (L22); with re-emitted photons leaving, as they physically do, the divisible reading registers nothing at a balanced split (L23) | Run 2026-09-08; see its `RESULTS.md` and `LEDGER.md` |
| `pair_case/` | JB's model at two analyzers: two photons sharing a hidden polarization, read by his projection rule (least rotation) or by Malus's law at each polarizer; the quantum comparator; the model's phase-dependent efficiency as a detection gate. What correlation function, what CHSH value, and the one route above Bell's bound | Run 2026-09-08; see its `RESULTS.md` and `LEDGER.md` (Bell's sawtooth, CHSH 2.00; the gate rides the Garg–Mermin line, the detection loophole). Addendum: JB's two options, a synchronized bulk phase shared by both detectors (a common cause, S = 2.00) against particle 2 set by A's outcome (the projection postulate, S = 2.83) |
| `vacuum_race/` | JB's rule of 2026-09-14: among the atoms of both pixels the one with the highest interference-weighted deposit takes the quantum, each atom carrying its own vacuum term, one random number per atom independent of the number of atoms. Which distribution of the vacuum term gives Born between the pixels, whether the term must multiply the intensity or may add to the amplitude, and whether the best phase-matched atom wins within a pixel | Run 2026-09-14; see its `RESULTS.md` and `LEDGER.md` (only a Gumbel term on the deposit gives Born, and with it the rule is R2, the golden rule as a race; Gaussian and bounded terms give the brighter pixel every time; a term on the amplitude gives the square root of Born). Addendum: the term frozen per atom against redrawn per arrival (a frozen term gives each detector all or nothing and Born only across detectors; the Fano factor of one detector's split reads the vacuum's correlation time, F = τ, against F = 1 in experiment) |
| `exchange_fairness/` | Fairness in the conserving exchange class: whether the √e noise scaling that the record's Paper 1 called the unique fair point is needed once the pot is closed; independent noise against two conserving kernels, with absorption, rate-linear and threshold terminations, and fixed parcels with a reflecting zero | Run 2026-09-14; see its `RESULTS.md` (in a closed pot any symmetric zero-mean exchange is fair; Paper 1's √e model does not close the pot; a reflecting zero biases toward the poor from three sites up; a threshold biases toward the rich; two registered misses kept) |
| `born_selection_revised/` | The specification paper: what a microscopic selection mechanism must supply to reproduce Born and one record per quantum, six clauses each with a registered run, and the fresh vacuum reading as the worked example under its right name | Draft v0.2, 2026-09-14; `OUTLINE.md`, `PAPER_DRAFT.md` |
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
    python3 frozen_fresh.py           # addendum: the vacuum term frozen per detector, refreshed every τ events, or fresh; one detector's counting statistics; about 10 s
    python3 staggered_arrival.py      # addendum 2: the packet reaching the two pixels at different times; the survival hazard against a local conditional reading; about 15 s
    python3 staggered_conditional.py  # addendum 2b: the rates on the normalized conditional state; about 200 s
    python3 frozen_multipliers.py     # addendum 3: a frozen per-atom multiplier with light tails against the infinite-mean e^Gumbel; about 65 s

    cd ../exchange_fairness
    python3 exchange_fairness.py      # independent noise, two conserving kernels, three terminations, fixed parcels; about 20 min

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

## The session of 2026-09-14, in plain language

The day began with a question of arithmetic. JB asked whether the projection of an amplitude onto
the real axis is the same thing as Born's rule. It is not: the projection is one leg of the
triangle and Born's rule is the square of the hypotenuse, and the two agree only after the angle
has been averaged away. That averaging is exactly why the earlier studies reproduced Born without
deriving it. What JB had in mind, the square root of the amplitude times its conjugate, is the
modulus, and its formula is a two-line identity. Born's rule, that the modulus squared is a
probability, is the postulate, and no algebra reaches it.

From there JB proposed a mechanism. Within a detector's aperture there are many atoms with
different phases; the atom whose interference with the incoming particle is strongest takes the
quantum; and each atom carries its own vacuum term, a random contribution that belongs to the atom
and knows nothing about how many atoms are in the aperture. The question put to the computer,
written down before it was run, was which kind of vacuum term makes that rule give Born's share
between two pixels.

**What came out.**

- *Only one kind of vacuum term works, and it is the golden rule in disguise.* If each atom's term
  is drawn from the one distribution that mathematicians call Gumbel, applied to the logarithm of
  the deposit, then "highest wins" is exactly the same as "first to fire," with each atom firing at
  a rate proportional to its share of the light. The share between two pixels comes out at Born's
  value, whatever the number of atoms and whatever the brightness. This is the first mechanism in
  the whole program that gives one record per quantum and the Born share at once. But two things
  were put in. The square of the amplitude was written into the rate, which is Fermi's golden
  rule, so Born is put in, not derived. And the single winner exists only because every other
  atom's clock stops the instant the first one fires, in the other pixel too, which is the
  one-quantum stop across pixels that the Dirac–Kuramoto record names as its premise P5. An
  earlier version of this summary said the race needed no such stop; that was wrong and is
  corrected. The scale of the vacuum term is also a
  free knob, and Born is one setting of it.
- *Every ordinary kind of noise fails.* Gaussian or bounded vacuum terms, whether they multiply the
  deposit or add to it, give the brighter pixel every single time once the aperture holds a hundred
  atoms. The reason is not the vacuum but the word "highest": the largest of many draws crowds up
  against its ceiling, and the brighter pixel's ceiling is higher. Additive noise fails a second
  way too, by reading the difference in brightness where Born reads the ratio.
- *The vacuum must multiply the intensity, not add to the amplitude.* A vacuum contribution that
  interferes with the particle's wave, the natural picture of an added amplitude, selects in
  proportion to the amplitude rather than its square and gives the square root of Born.
- *The vacuum term must be fresh for every arrival.* JB then asked whether the vacuum acts on the
  bound electron first, before the free particle arrives, so that the atom is already dressed. The
  order matters and the runs pick his order. But a term that is a fixed property of the atom makes
  every single detector give all its counts to one pixel or the other, with Born appearing only as
  an average over many detectors. Real detectors show Born's split inside one device, with the
  plain binomial scatter of independent coin tosses. That scatter is the model's own observable:
  it reads how long the vacuum term stays put between arrivals, and experiment says under ten
  nanoseconds. Any vacuum timescale in silicon is a hundred thousand times shorter, so the model
  passes this test, and it is the one test in the program it passes for a reason it supplies
  itself. "Dressed first" has to mean dressed afresh each time.

**Where this leaves the model.** Section 12 of `MODEL_SECTION.md` states it as the fresh vacuum
reading: three premises, of which two are the golden rule's rate law and its lack of memory; a
detector mechanism that then computes which pixel, which atom, one whole record per quantum, and
binomial counting statistics; and the same disclaimer as before, that Born is inherited. The open
question is now sharper than it was. Whatever the vacuum is, in JB's picture a four-dimensional
thing projected onto three dimensions and time, it has to act on a bound electron as an
exponential clock at one particular scale, and not as an amplitude in the electron's state. The
Wick rotation stored in `WICK_ROTATION.md` is where standard theory gets exponential clocks from,
and it delivers their shape but not their rate.

**Where the ideas came from.** Three of JB's ideas drove today's results; the formalization and
the runs were Claude's. First, that the quantum vacuum interacts with the individual atoms of the
detector's bulk, one term per atom. That turned selection into a race with a clock at every atom
and made the distribution of the per-atom term the question, which is what found the one family
that gives Born and, when JB asked whether the atom is dressed before the particle arrives, gave
the frozen-versus-fresh test and the counting-statistics bound. Second, that the order of the
three interactions, vacuum, bound electron, free particle, does not commute. That separated the
vacuum multiplying the atom's rate, which gives Born, from the vacuum adding to the particle's
amplitude, which gives the square root of Born, and fixed the ordering the numbers pick. Third,
from the session of 2026-09-07, that the aperture holds many atoms with varying phases. That is
what makes the difference between "highest wins" and "first to fire" visible at all: with one
atom per pixel the rule without a vacuum term sits within a percent of Born and the others
within twenty, and from ten atoms the largest of
many draws climbs to its ceiling and the brighter pixel takes everything unless the vacuum term
is exponential.

**Later the same night: when the wave reaches the two pixels at different times.** The review
ledger written for the Born Selection paper predicted that the race would give the click to
whichever pixel the wave reaches first, a failure the Dirac–Kuramoto record had already found for
its own postulate. It does: a tenth of a packet length of stagger hands the nearer pixel nearly
every record. JB said the failure could be removed, and it can, at a price the runs made exact.
Letting each pixel consume its own share of the wave, with the pixels tied only by the rule that
the first to fire stops the other, gives 0.95 instead of 0.80; my registered prediction that it
would give 0.80 was wrong and is kept. What gives 0.80 at every stagger is the rule quantum
mechanics uses for photodetection: when one pixel fires the other's clocks stop, and when one
pixel finishes without firing the other pixel's rate is raised to compensate. Both halves act
across the gap between the pixels. So under staggered arrival the split between pixels is not
decided by any race at all; it is the share of the wave each pixel was given, and the race only
picks the atom within the pixel.

**Consistency pass, later the same night.** The internal review round of the Born Selection paper
(its `born_selection_AI_review_2026-09/findings_classification.md`) found statements here that
disagreed with each other or with the record. Corrected in this pass, with the sessions of
2026-09-08 and 2026-09-14 deciding: the "within a thousandth of the commitment time" figure is
withdrawn as an established number in `MODEL_SECTION.md` §12.3 and `vacuum_race/RESULTS.md`, since
the staggered-arrival addendum itself showed the stop to be the conditional-state update, with no
delay and no speed, which is what *One World, One Cut* v1.5 had said; the claim that Gaussian and
bounded vacuum terms give the brighter pixel every time from a hundred atoms is restated with the
additive-Gaussian miss (0.966 and 0.990) in `MODEL_SECTION.md` §11 and §12.4; §12.3's "a proportion,
not a maximum, decides" is qualified to simultaneous arrival; and this README's ledger range and
provenance sentence are corrected. Ledger L29.

**Also from today.** Three registered predictions missed and are kept: the local conditional reading above, and two in their wording: the additive
Gaussian noise approaches the brighter-pixel-always limit more slowly than registered, and the
excess scatter from a detector with frozen bulk phases at ten atoms per pixel came out above the
interval registered for it, which was too narrow for a heavy-tailed quantity. Every other registered
check passed. Ledger entries L24 to L28 record the departures.

## The session of 2026-09-14, the evening and night, in plain language

The day's work was reviewed and then rebuilt. JB asked for a multi-agent critical review of the
Born Selection paper with its evidence spread across this repository and the Dirac–Kuramoto
record. Five reviewers were run under the program's own instrument, each with a different
specialism and none able to see the others or the list of findings written in advance: a
probability theorist, a detector experimentalist, a philosopher of physics, a scholarship checker,
and an auditor who re-ran every script here and checked 294 numbers. The reviewers did not pass
the paper. Every number in it reproduced, and the failures were in what the numbers were said to
show: the paper had come to hold two different selection mechanisms at once; the vacuum race, once
amended, is the standard quantum description of photodetection and had to be called that; the
claim that a frozen per-atom factor "does not average out" was an artefact of a factor with
infinite mean; the wavelength lines in the energy study were an artefact of the toy's pixels
catching only 85 % of the wave; a citation for a beam-splitter experiment applied silicon's gap to
photomultipliers; and the earlier record's own corrections had not been carried through. Two of
those were errors I introduced that same day while correcting other errors. All of it is in the
review folder of the record, and the corrections that were matters of consistency were applied
the same night, with the sessions here as the guide where documents disagreed.

Then JB made three decisions in a row. That the vacuum does not interfere with the incident
packet at a photon-counting detector, only with the atom, which is the textbook position and
settles on physical grounds why the vacuum has to multiply the atom's rate rather than add to the
wave. That the new version of the paper should be written on a clean slate here, in MathX, rather
than patched onto the old one. And that it should be a specification: not a mechanism and not a
bare negative result, but the list of what any microscopic account of detection must supply to
reproduce Born's rule and one record per quantum, each item with a run showing what goes wrong
without it.

**Two runs the specification needed, both registered before they were run.**

- *A frozen per-atom factor with an ordinary distribution.* Uniform, lognormal and exponential
  factors all wash out by a thousand atoms per pixel, so a single detector shows the plain
  binomial split; only the infinite-mean factor of the earlier run does not, and even its average
  over many detectors is not Born. The clause from the afternoon, that the factor must be redrawn
  for every arrival, was wrong and is replaced: a frozen factor is harmless if its tail is light,
  which any bounded physical factor is. My prediction that the infinite-mean case would still
  average to Born across detectors missed and is kept.
- *Fairness when energy is exchanged among sites before one commits.* In a pot closed to one
  quantum, any symmetric zero-mean exchange is fair, whatever its dependence on the energies, so
  the square-root noise scaling the old paper called the unique fair point is neither necessary
  nor unique there, and the old paper's own noise model turns out not to close the pot at all: its
  winner ends with a median 0.72 of a quantum. What does spoil fairness is a floor a depleted site
  can bounce off, and what spoils Born is stopping the exchange short of the full share. My
  prediction that the floor would bias a two-site game missed: with two sites a depleted site
  means the other has everything and the game is over, so the bias needs three sites. Kept.

**The draft.** `born_selection_revised/PAPER_DRAFT.md`, six clauses each tied to a run: the rate
is proportional to the intensity, the modulus squared of the amplitude, never the modulus, and
interference lives inside that square, upstream of selection; selection is a race and not a
maximum; a port needs about ten sites for the phase factor to average; the clock is memoryless and
a frozen factor must be light-tailed; the rates are those of the normalized conditional state, so
that under staggered arrival the share between ports is the branch weight and no race decides it;
and an exchange stage, if there is one, needs a closed pot, a symmetric zero-mean kernel, an
absorbing floor and a termination at the full share or by a rate linear in the share. The worked
example is JB's model under its right name, the quantum-jump description of photodetection with
the vacuum as the name of each atom's clock. JB is first on the byline; the three ideas were his
and the formalization and runs were mine.

**Three questions of JB's, answered in the draft's second section.** Whether the phase of the
amplitude could be what locks, or what projects onto the real axis, to give two values instead of
a continuum: no, on the record's evidence, since a bound site has no free phase to lock, a
real-axis projection gives fifty-fifty whatever the amplitude, and two values come from the
apparatus resolving two modes and the exchange coming in units. Whether it is better to say the
amplitude projects onto the axis of the measurement setup: yes, in the Hilbert-space sense, one
complex component per mode, its modulus squared the intensity, with the absorption at one port as
the event; the word "projection" is now used in that sense only. And whether in Stern–Gerlach the
magnet selects: no, the gradient sets the axis and separates the components reversibly, as
recombination shows, and the screen selects by absorbing the atom in one of the two regions; the
detector registers a position that means spin only because the gradient tied the two together.

**What is left for JB.** A reading of the draft's last two sections; a check of the references,
which are from memory; whether the two predictions of his model that experiment excludes stay in
as the model's or go out with the phase weight; the three decisions in the authors' response in the
record, whether the old paper is retired or rebuilt, whether its clause on frozen factors is
reduced or deleted, and which ontology it holds; and whether this repository becomes public, since
the draft cannot be checked from its pointers otherwise.
