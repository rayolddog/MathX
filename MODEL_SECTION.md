# The microscopic model: a draft section for the record

**Date:** 2026-09-07. Written at JB's instruction as a draft of the model section a rewritten
Born Selection paper would need. It is a record of the mechanism as computed in MathX, with
Born's rule assumed, and it is not part of the closed Dirac–Kuramoto record. Every number is
traceable to a results file named in §9; nothing here is a measurement.

## 0. Read this first

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
