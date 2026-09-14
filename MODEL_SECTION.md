# The microscopic model: a draft section for the record

**Date:** 2026-09-07. Written at JB's instruction as a draft of the model section a rewritten
Born Selection paper would need. It is a record of the mechanism as computed in MathX, with
Born's rule assumed, and it is not part of the closed Dirac–Kuramoto record. Every number is
traceable to a results file named in §9; nothing here is a measurement.

## 0. Read this first

**Addendum of 2026-09-08 at the end (§10):** the plain-form consequences of the SPAD study, the
test of L13 and the pair case, answering whether Born selection is reproduced, whether the local
hidden variable is ruled out, and whether the model implies a violation of detector independence.

**What this section claims.** A single-particle mechanism for the event at a detector: which atom
takes the quantum, how much of the quantum's energy it takes, how the record forms, and where the
mechanism differs from the golden rule. The projection of complex vectors enters as a premise in
two places (§2), and everything after it is computed.

**What it does not claim, stated first so that no reader takes it from the tables.**

1. It does not derive Born's rule. The site statistics inherit it (§3). Every dynamics added on
   top of the selection rule moved the statistics away from Born, toward the brighter site
   always winning (§4.2).
2. It does not complete the third step of the measurement account, the dissipation of the
   packet into one record. Every dynamics tried for that step either spread the record over a
   cluster, left most of the quantum uncollected, or made two records (§5).
3. It has not been run for a pair. In its present form the selection rule is local, and a local
   rule gives the correlations of Bell's own foil, not the measured ones (§6).

## 1. The particle and the bulk

A particle carries two phases rotating in opposite senses, φ₁ = +ωt + a and φ₂ = −ωt + b, its
first and second spinor. The quantity carried through an interaction is the internal offset
δ = b − a, the second spinor's phase minus the first's at the instant of contact. An atom of the
bulk carries the same pair. In a random bulk the atoms' second-spinor phases are uniform on the
circle and independent, so their offsets δ_i are too. At contact the particle's first spinor is
shifted to the atom's and the same shift is applied to its second spinor, so δ is maintained;
no pull acts on the second spinor. The mismatch of the second spinors at atom i is therefore
Δφ_i = δ − δ_i, the difference of two offsets, uniform on the circle whatever the particle's own
offset is. (Second study, D1 to D3.)

Nothing in this paragraph names an energy, a mass, or a sign of energy. The identification of the
two spinors with the positive- and negative-energy components of a Dirac spinor was used in one
place only, the free-flight check of §4.4, and it is separable from the rest.

## 2. Where the projection of complex vectors enters

Twice, and as premises.

- **The deposit weight.** The energy the particle can deposit at atom i is highest where the
  second spinors are closest in phase, because they interfere: w(Δφ_i) = cos²(Δφ_i/2), which is
  |1 + e^(iΔφ)|²/4, the squared overlap of two unit phasors. With unequal amplitudes the same
  form carries a visibility V ≤ 1 (second study, A1 and Part D).
- **The energy available at an atom.** It is the squared amplitude of the wave there, |ψ(x_i)|²
  (second study, D4). In the trajectory reading the wave's density is transported by the
  Bohm–Dirac velocity and the record starts where the particle arrives; in the field reading the
  density is what each atom has to draw on. Both readings were run (§4.4, §4.5).

So the deposit at atom i is E_i ∝ |ψ(x_i)|² cos²(Δφ_i/2): a projection onto position times a
projection between two phases.

## 3. The selection rule and what it inherits

One atom takes the quantum, atom i with probability E_i / Σ_j E_j (rule R2, JB's choice of
2026-09-07). Because Δφ_i is uniform and independent of position, the phase factor averages to
one half over any aperture with more than about ten atoms per lobe, and the probability that a
lobe is recorded is its share of |ψ|²: Born's rule, inherited from the second premise of §2. At
intensity ratio 4 between two lobes the recorded share is 0.704 with one atom per lobe, 0.794
with ten, 0.800 with a thousand, against 0.800 (second study, Table 2). Nothing in the phase
mechanism changes it, at any interference visibility (Part D). A cost of the phase shift lifts it
above Born by about the cost divided by the total, a threshold effect (Part E).

## 4. What the mechanism adds

**4.1 Which atom, and how close in phase.** Among N atoms in reach there is one within 0.45 rad of
the particle's second-spinor phase, taking at least 0.95 of the energy weight, with 95 %
confidence once N ≥ 20, and one within 0.20 rad once N ≥ 46. Under the selection rule the atom
that fires is 53.5° off in phase on average, whatever N; only the deterministic limit of the rule
picks the closest atom, at π/(N + 1). (Second study, Table 1.)

**4.2 The site statistics under the rule and under its dynamics.** Under the rule itself, Born.
Under any dynamics that lets the record form by a race the brighter lobe wins more than its
share: the deterministic limit gives the brighter lobe every time once the aperture holds a few
atoms; a charging threshold at every atom gives the brighter lobe first in 0.95 of trials at a
Born share of 0.80 (third study); a cost of the shift gives an exponent of 2.05 to 2.25 against
Born's 2 (second study, Part E). A single seed sampled by the rule restores Born by construction,
0.85 against 0.80 in twenty trials (third study, seeded charging).

**4.3 The phase-ordered detector.** If the bulk's second-spinor phases cluster about a mean μ with
mean resultant R, the deposit and the recorded share follow ½[1 + V R cos(δ − μ)]: the detector
reads the incident particle's internal offset with visibility V R (second study, Table 3 and
Part D). This is the one place the mechanism says something the golden rule does not, and it is
where it would be tested. The one phase-ordered detector the record examined, the superconducting
nanowire, showed no such dependence (ledger, L3).

**4.4 Free flight.** With the two spinors read as the positive- and negative-energy components of
a Dirac spinor, the internal offset is a constant of the free motion to machine precision, and
the Bohm–Dirac trajectory carries it as the phase of its zitterbewegung at the instant of
contact, v(0) = sin θ cos δ exactly for an admixture angle θ. The beat frequency of a localized
packet is about twice its mean energy, not twice the rest energy, and the registered prediction
of the latter missed (second study, Part C). Arrivals follow |ψ|² by equivariance, so the
trajectory reading adds nothing to the site statistics; it decides which atom within reach.

**4.5 The record as a cluster, and collection.** When every atom in reach can charge and latch,
the record is a cluster of sixty to two hundred latched atoms in one lobe, which absorb the
density where it lies. A cluster is what a real record looks like. A single latched atom, by
contrast, collects about 0.05 of the quantum in sixty time units, a fixed-width sink at one atom
collects 0.38 in forty Compton times of Dirac flow, and a sink collects the whole quantum only
when its absorbing region has grown over the whole packet, which in one dimension is a blanket
and not a channel (third study; second study, Part F).

**4.6 The beam splitter.** If the two spinors carry unequal energy shares (1 ± ε)/2 and the
mismatch with the splitter's random phases decides which output takes the larger share, the split
is uneven by ε on every event and even on average. A smooth dependence on the mismatch leaves a
sixth of events near even at ε = 0.2 and fails the single-photon anticorrelation bound; a
switch-like dependence at ε = 0.2 with slow charging gave no double record in fifty events, and
ε = 0.2 is exactly where the switch form's fringe visibility √(1 − ε²) meets the 0.98 measured
with single photons in 1986. The window has no width: at ε = 0.15 one event in fifty doubled, at
ε = 0.1 half did, and any visibility above 0.98 closes it (second study, Part H).

**4.7 Detection efficiency tied to phase.** With one seed per quantum, a seed of low phase weight
cannot reach threshold, and one quantum in twenty made no record (third study, seeded charging).
With the band form of the weight the miss rate is (1 − ε/π)^N (second study, C4).

## 5. What it has not done: the third step

Every dynamics for dissipation was registered before it was run, and each failed the step in a
specific way, on record:

| Dynamics | Result | Where |
|---|---|---|
| a well plus absorber at the winning atom, growing with what it collects | draws nothing in; collects 0.38 of the quantum in 40 Compton times; a stronger absorber reflects | second study, Part F |
| the absorbing region widening with collection | collects everything, by covering the far lobe and eating it in place | Part F, channel run |
| two sinks with reversible charging, a latch, and vacuum noise | one record only when the head start exceeds the suppression time; at equal lobes, two records in every trial | Part G |
| a screened field from the charges with a conductance tensor, in two dimensions | the record forms as a cluster where the density lies; no flow concentrates, no channel forms; a second lobe that can reach threshold always does | third study |
| charging restricted to sampled seeds | near-equal deposits go from 0.95 to 0.00 of trials; doubles fall to 0.50 and 0.85; one seed gives one record and collects 0.05 | third study, seeded |

Two things remain open, cleanly separated by the last row: a rule that suppresses a second seed
once the first has latched, and a way for a single seed to collect the rest of the quantum. The
only suppression that worked, the growing region of Part F, is not a flow but the sink spreading
over the wave; the only complete collection was the same thing. A rule that acts across the
aperture at once would do both, and it would be the projection postulate written as dynamics.

## 6. The pair case, not run

The selection rule of §3 is local: each detector draws from the intensity arriving at it. For an
entangled pair the joint statistics are not the product of the two arriving intensities. To give
the measured correlations, the record at one detector must condition the amplitudes at the other,
which is a rule acting across the pair, and the model does not yet contain one. Without it the
model reduces to its marginals and gives the correlations of Bell's local foil. A rewritten paper
must either restrict the model to single particles and say so, or state the conditioning and
what the model can say about it. The runs so far bear on neither; the ledger's first entry
records the form.

## 7. Departures from quantum mechanics and experiment

Kept in `second_spinor_phase_race/LEDGER.md`, dated, per JB's instruction of 2026-09-06 that the
model be built before it is judged. In one line each: the internal offset as a variable local to
particle and detector (L1); the deposit rule against the golden rule's lineshape (L2); the
phase-ordered nanowire that showed no phase dependence (L3); the energy content of a lag read as
a negative-energy admixture (L4); Bohm's equivariance fixing the site statistics before the race
begins (L5); the beat frequency of a packet (L6).

## 8. The three steps of the framework paper

| Step | Content | State in this model |
|---|---|---|
| 1, free flight | transport of the density; the offset carried unchanged | computed; any spreading transport gives the same arrivals |
| 2, contact with a bound electron | the match, the deposit, the selection | computed; Born inherited, the phase mechanism deciding which atom and how much |
| 3, dissipation into one record | the sink, the cluster, suppression of other records | open; the failure modes quantified in §5 |

## 9. Sources

- `kuramoto_one_two_rotations/RESULTS.md`: the first-spinor match with an Adler pull, superseded
  by the model above for the second spinor and kept for the record.
- `second_spinor_phase_race/PREREGISTRATION.md`, `RESULTS.md`, `LEDGER.md`: §1 to §4, §5 rows 1
  to 3, §7. Parts A to H, with every registered prediction and its outcome.
- `sink_field_2d/PREREGISTRATION.md`, `RESULTS.md`: §4.2, §4.5, §4.7, §5 rows 4 and 5, including
  the reversed-sign first run kept in `results/first_run_wrong_sign/`.
- The framework paper's three steps: JB, 2026-09-07, from the closed record.

## 10. Consequences of the runs of 2026-09-08, in plain form

**Added 2026-09-08 at JB's instruction, after the SPAD study (Parts A to F), the test of L13 and
the pair case.** Three questions, three answers, then the table. Sources: the study folders named
in the last column.

**Is Born selection reproduced?** Reproduced, never derived. In every run where the Born share
appeared, it came from an ingredient that already contained |ψ|². The sampled rule R2 draws an
atom in proportion to intensity times the phase weight, and the intensity is the square put in by
hand. The trajectory reading gives arrivals proportional to |ψ|² by Bohm's equivariance before any
race begins. Every selection dynamics added on top of those two moved the statistics away from
Born, always toward the brighter site: the deterministic rule gives the brighter lobe every time in
the field reading, charging races give it first in 95 % of trials, and the empty pixel with no head
start gives 0.88 against 0.80 or the device's own asymmetry. So the program's original aim, a
selection dynamics that produces the square, is not met. What the model does is inherit Born and
then distort it.

**Is the local hidden variable ruled out?** Yes, at both ends, and by different means. As a model
of the pair it is a local hidden-variable model in Bell's sense, and the run gives Bell's own
sawtooth with CHSH 2.00. Experiment gives 2.83. The only way the model rises above 2 is by gating
detection on the hidden phase, which is the detection loophole, and the loophole-free experiments
of 2015 closed it. That exclusion is the experiments', not the runs'; the runs show the model
belongs to the class that is excluded. At the single-particle end the hidden variable turned out
not to be hidden or not to be free. In a bound electron the second spinor's phase is fixed by the
Dirac equation at ±90°, the same in every atom. In a free electron it is set by the direction of
motion. In a photon it is the polarization, a measured quantity. In a passive silicon detector the
phase mechanism has no observable consequence at all, and in a driven one it contradicts coherent
control in its form.

**Does the model suggest a violation of detector independence?** In one specific way, yes, and
that is precisely its problem. The weighted deposit makes whether a photon is detected depend on
its hidden phase relative to the analyzer axis. That is a violation of the fair-sampling
assumption, detection independent of the hidden variable and the setting, and it is what lets the
gated model reach 2.83 among coincidences at 83 % conditional efficiency. Two things stand against
it. The loophole-free experiments were built to work without fair sampling and still found the
quantum value. And in silicon the gate never acts: with the bound electron's real second-spinor
amplitude the weight stays above 0.9 and every photon that reaches the gap makes a pair. In the
other sense of the phrase, the outcome at one detector being independent of the setting at the
other, the model as built respects it, which is why it cannot exceed 2, and the measured
correlations say that something must not. The model has no rule that acts across a pair.

| question | answer | where |
|---|---|---|
| Born reproduced | yes, from \|ψ\|² put in by the deposit or the trajectory | second study Table 2; SPAD study main run |
| Born derived | no; every dynamics distorts it toward the brighter site | second study Parts G, H; third study; SPAD study, the empty pixel |
| one record per quantum | by assumption (the empty wave deposits nothing), by a stop across pixels, or by energy conservation in the divisible reading, which then fails balanced splits and calorimetry | SPAD study L7, L10, L21 to L23 |
| local hidden variable | Bell's sawtooth, CHSH 2.00, also with a bulk phase shared by both detectors; excluded by experiment | pair case, L17, L18, L20 |
| detector independence | violated by the phase gate, which is the closed detection loophole; never in silicon | pair case Table 3; SPAD study Part F |
| the phase itself | fixed in bound electrons, kinematic in free ones, polarization in photons | SPAD study Parts E and F |

**Added later the same night, after Part G and the pair options.**

*The third step, once more.* JB's release cycle at the excited electron, ejection above the gap,
re-emission, phonon loss, with the quantum's energy tracked to closure, is the first dynamics in
this program that keeps one record per quantum without a stop across pixels and without
assuming the empty pixel silent: energy conservation does it whenever the gap exceeds half the
photon energy. It does so only in the reading where the deposit is divisible, and that reading
pays twice, with a share that is a step in the intensity ratio (the brighter pixel takes every
record at r = 4) and with records that carry 0.45 to 0.72 of the quantum where calorimeters
measure the whole. When the re-emitted photon is allowed to leave, as it physically must below
the gap, the divisible reading registers nothing at a balanced split above 470 nm and is
excluded outright. The indivisible reading keeps Born, exclusivity, whole-quantum records and
closure at once, and derives nothing, because the indivisibility of the deposit is the postulate.
The two readings cannot be combined: a divisible deposit that always ends in a whole record
draws the remainder from the wave at the moment of ejection, which is the stop across pixels
under another name. So the third step now has its final form in this model: it closes either by
assuming the quantum or by assuming the stop.

*The pair, once more.* JB's second option, both detectors' first spinors synchronized to one
unchanging bulk phase and the second spinors rotating with them, with no communication, is a
common cause fixed before the settings are chosen. It gives exactly the correlations of the
shared-angle model, Bell's sawtooth and CHSH 2.00, because a common cause is what the bound
already allows. His first option, particle 2's phase set by the outcome at A, reproduces the
measured cos 2θ and 2.83 when the second analyzer reads by Malus, and over-shoots to 4 when it
reads deterministically, so the deterministic projection cannot be the rule at the second
analyzer either. The nonlocality is in the update and signals nothing. A bulk that carried A's
outcome to B would be the first option with the bulk as carrier, and any version of it at light
speed or below is excluded by the spacelike-separated experiments.

**What this does to §0.** Its three disclaimers stand and are now numbers: Born is inherited
(§3) and distorted by dynamics (§4.2); the third step is open (§5), closed only by the trajectory
reading's assumption; and the pair case, no longer "not run" (§6), gives Bell's foil. The runs added two things §0 did not foresee. The fate of the phase variable itself: Parts E and
F show it is not free in the detector, and Part D shows where a driven detector would expose it,
inside one T2, where the form of the deposit rule then disagrees with coherent control (L14). And the
shape of the third step's closure: energy conservation can supply exclusivity, but only in the reading
that loses Born and the whole-quantum record, so that the model's choice is between assuming the
quantum and assuming the stop.

## 11. Addendum of 2026-09-14: the vacuum race

**Added at JB's instruction after `vacuum_race/`.** JB's rule: within the aperture the atom with the
highest interference term with the incident particle is selected, each atom carrying its own vacuum
term, independent of the number of atoms. The run varies only the distribution of that term and how
it enters the deposit.

*What it adds to §5.* One dynamics that closes the third step with the Born share and one record per
quantum at once: a per-atom vacuum term that is a Gumbel draw on the logarithm of the deposit, or
equivalently an exponential waiting time with rate proportional to |ψ|² cos²(Δφ/2), the first atom
to fire taking the quantum. Under it "highest wins" is exactly rule R2 of §3, atom by atom (0.702,
0.794, 0.801, 0.802 against R2's 0.704, 0.794, 0.800, 0.800 at N = 1, 10, 100, 1000), independent
of N and of the absolute intensity, with the best phase-matched atom the most likely winner in 2/N
of trials and never certain. It does not assume the empty pixel silent; its exclusivity is by
construction of the maximum, which the run assumes and does not test. (Corrected 2026-09-14,
later: an earlier wording here said the race "needs no stop across pixels". It does. The first
atom to fire takes the quantum only because every other atom's clock, in the other pixel included,
is cancelled at that instant, and that cancellation is the stop, the record's premise P5, acting
faster than light crosses the pixel separation; see §12.3.)

*What it does to §0.* Nothing. The rate proportional to |ψ|² is Fermi's golden rule, and the
square sits in the vacuum's distribution. Every other form of the term tried, Gaussian or bounded,
multiplicative or additive, gives the brighter pixel every time once the aperture holds a hundred
atoms, because the maximum of N draws concentrates while the per-atom term does not; additive forms
also read the intensity difference rather than the ratio. A term that enters by interference with
the amplitude, linear in |ψ|, gives √r/(√r + 1), the square root of Born. And the scale of the
Gumbel term is free, the share being r^β/(r^β + 1), so Born is the point β = 1 and nothing in the
model fixes it (L24, L25). The model can carry this race as its detector mechanism given Born. It
is not a source of Born, and disclaimer 1 of §0 stands.

*Added later the same day, after the frozen-against-fresh addendum.* The vacuum term has to be a
fluctuation faster than the arrivals. Frozen per atom it gives each detector an all-or-nothing split
and Born only across an ensemble of detectors; refreshed every τ arrivals it gives one detector a
Fano factor F = τ for its split, against F = 1 in every SPAD experiment. JB's ordering, the vacuum
acting on the fermion before the free particle arrives, is the ordering the run favors (V1 against
V2), but only with the dressing renewed for each arrival; a standing property of the atom, such as
the SPAD study's fixed resistances, cannot play the part (L26).

## 12. The fresh vacuum reading

**Added 2026-09-14 at JB's instruction, after §11.** This section states the model as it stands
once the vacuum race and its addendum are taken as the detector mechanism: the fresh vacuum
reading. It is written as the model section a rewritten paper would carry for step 2 and step 3 of
§8, in the same form as §1 to §5: premises first, what is computed after, and what the reading does
not do stated with it. Every number is from `vacuum_race/results/`.

**12.1 The premises, stated as premises.** Three, and the reading rests on all three.

1. The deposit at atom i is E_i ∝ |ψ(x_i)|² cos²(Δφ_i/2), the two projections of §2: the
   intensity at the atom and the squared overlap of the two second-spinor phasors.
2. Every atom of the bulk carries a vacuum term that multiplies its deposit, drawn independently
   of every other atom's and of the number of atoms in reach, from the Gumbel distribution on the
   logarithm of the deposit at unit scale. Equivalently, each atom has an exponential waiting
   time with rate proportional to E_i.
3. The term is renewed for every arrival: its correlation time is shorter than the interval
   between quanta.

The first premise was already §2's. The second and third are new and are what §11 found the race
needs. Premise 2 is the golden rule's rate law, and it is where Born sits in this reading; premise 3
is the golden rule's memorylessness. Neither is derived here (L24, L26).

**12.2 The rule.** The atom with the highest vacuum-weighted deposit takes the quantum. Under
premise 2 that is the same statement as: the first atom to fire takes the quantum, and atom i fires
first with probability E_i / Σ_j E_j. So the rule of this reading is rule R2 of §3, now as a
dynamics with a clock rather than a sampling instruction. The identity is exact and was verified
atom by atom (0.702, 0.794, 0.801, 0.802 against R2's 0.704, 0.794, 0.800, 0.800 for N = 1, 10, 100,
1000 at r = 4; 0.800 at every N without the phase weight).

**12.3 What the reading computes.** Given the three premises, without further choices:

- *Which pixel.* The intensity share times the ratio of phase sums, Born once a pixel holds about
  ten atoms, and independent of N and of the absolute intensity (unchanged at ten times the light,
  Table 4 of the study). This is the one selection dynamics in the program that does not distort
  Born toward the brighter site (§4.2), because a proportion, not a maximum, decides.
- *Which atom.* The best phase-matched atom of the winning pixel is the single most likely winner
  and wins in about 2/N of arrivals; it is never certain. The winner's phase mismatch has the
  cos² density of §4.1, 53.5° on average.
- *One record per quantum.* One atom fires first; no second record, and no assumption that the
  empty pixel is silent. The whole quantum goes to that atom. The empty pixel's atoms are in the
  race with their own deposits and lose in proportion, which is the Born share seen from the
  other side. But the single winner is not free: the race has one winner because every other
  clock is cancelled the instant the first fires, and for pixels or beam-splitter arms that are
  separated, that cancellation acts across the separation faster than light. It is the stop
  across pixels, the record's premise P5 (Paper 1 §2, third remark; NEGATIVE_RESULT.md), which the
  record established must act within about a thousandth of the commitment time to match
  single-photon coincidence data; independent clocks that keep running give two records at a
  balanced split. The vacuum race assumes it by taking the maximum over both pixels at once, and
  derives it no more than the readings of §10 did. (An earlier wording of this item, "no stop
  across pixels", was wrong and is corrected 2026-09-14.) With that assumed, the third step
  closes in the sense of §8: the
  cluster, the sink and the channel of §4.5 are then what happens after the winner is chosen, and
  none of them has to choose.
- *One detector's counting statistics.* The split of M quanta between two pixels is binomial:
  Fano factor 1 (0.9 measured against binomial in the run), the same in every detector. This is
  premise 3's content and it is the reading's own observable (12.5).

**12.4 What the reading does not do.**

- It does not derive Born. The square is premise 2's distribution, and the scale of that
  distribution is a free parameter: a Gumbel term at scale 1/β gives the share r^β/(r^β + 1), so
  Born is the point β = 1 and nothing in the model selects it (Table 3, L24). Disclaimer 1 of §0
  stands in this reading exactly as before.
- It does not fix the form of the vacuum term from anything else in the model. Every other form
  tried, Gaussian or bounded, multiplying the deposit or adding to it, gives the brighter pixel
  every time once the aperture is full, and the additive forms read the intensity difference
  rather than the ratio (Tables 1, 2, 4). A term that enters by interference with the particle's
  amplitude, which is what a vacuum amplitude added to the atom's wave packet would do, gives
  √r/(√r + 1), the square root of Born (L25). So the vacuum of this reading is not an amplitude in
  the atom's state; it is a rate. What in the vacuum's geometry makes it a rate at unit Gumbel
  scale is the open question of this reading, and JB's four-dimensional vacuum is recorded as a
  candidate, not modeled (study §0, C4).
- It does not make the phase mechanism observable. The phase weight enters the rate and averages
  to one half over a random bulk, as in §3; the reading inherits §4.3's ordered-detector prediction
  and §10's finding that the bound electron's second spinor is not free, unchanged.
- It does not reach the pair. The rate at each detector is local, so §6 and §10 stand: Bell's
  sawtooth and CHSH 2.00.

**12.5 Where the reading would be tested, and where it already passes.** Premise 3 makes the
correlation time of the vacuum term an observable: refreshed every τ arrivals, the term gives one
detector a Fano factor F = τ for its split (about τ/2 with the phase mechanism), and frozen it gives
each detector all or nothing with Born only across an ensemble of detectors (Tables 6 and 7).
Experiment has F = 1 at every counting rate used, so τ is below one arrival interval, under about
10 ns at the fastest SPAD rates. Any vacuum timescale in silicon, 10 to 100 fs, is five orders of
magnitude inside that bound. The reading passes this test without adjustment, and it is the only
test in this program that the model passes for a reason it supplies itself rather than inherits.
The same bound excludes the alternative: a standing per-atom property, the SPAD study's fixed
resistances or any lattice disorder, cannot be the vacuum term.

**12.6 The order of the three interactions.** JB's proposal that the vacuum acts on the fermion
before the free particle arrives is the order this reading uses, in a definite sense: the term
attaches to the atom's deposit (form V1) and not to the particle's amplitude (form V2), and the
two orders give different statistics, Born against its square root. The order is not commutative
in the model. But "before" means before this arrival and after the last, by premise 3; the atom is
dressed afresh each time. In quantum mechanics this is the dressed-atom picture with the Markov
property of the coupling to the vacuum's continuum of modes, and that is the standard origin of the
exponential clock (L26). The reading therefore has the same content as the golden rule at this
step, with the Born share carried in the rate and the exclusivity in the memorylessness.

**12.7 Summary line for §8.** Step 2 and step 3 in the fresh vacuum reading: computed, given the
rate law; the third step closes with one record per quantum, the whole quantum, the Born share and
binomial counting statistics at once, given the one-quantum stop across pixels (P5); and it
derives nothing that the rate law and the stop do not already contain.
