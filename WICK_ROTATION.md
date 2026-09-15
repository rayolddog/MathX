# Note: the Wick rotation, for JB to revisit

**Date:** 2026-09-08. Stored at JB's request ("store the memory of it in this folder; I will
revisit it in the future"). A note for later, not a study: nothing here was run, and the
references are given from memory, not checked today.

## JB's tangent, in his words (2026-09-08)

"A crazy tangent popped into my head stewing over 4D to 3D + spacetime. What if 4D could mean
any time and location can be reached. 3D + spacetime places a preferred direction of movement in
spacetime? Maybe 4D to 3D+spacetime creates the light cone?"

## What it is

- A Euclidean four-space has the distance s² = x² + y² + z² + τ². Every direction is alike; the
  set of points at a fixed distance from an event is a sphere; there are no null directions and
  no causal order. In that sense any point can be reached from any other and nothing separates
  "time" from "space".
- Minkowski spacetime has s² = x² + y² + z² − t². One direction carries the opposite sign. The
  points at zero distance from an event are no longer a single point but a cone, x² + y² + z² = t²,
  the light cone; points inside it are timelike-separated, points outside spacelike.
- The two are related by t = −iτ (or τ = it). This substitution is the Wick rotation (G. C. Wick,
  Phys. Rev. 96, 1124, 1954): the sphere of equal Euclidean distance becomes the hyperboloids of
  equal Minkowski interval, and the sphere's degenerate limit at zero radius opens into the cone.
  So the answer to the tangent is yes, in a precise sense: giving one direction the opposite
  metric sign is exactly what makes the light cone.

## Where physics uses it

- Quantum field theory: the path-integral weight e^(iS/ħ) becomes e^(−S_E/ħ), a real, positive
  weight, so that Minkowski amplitudes are computed in Euclidean space and continued back. Lattice
  QCD is done entirely in Euclidean signature.
- Gravity: Hawking's imaginary-time treatment of black-hole temperature, and the Hartle–Hawking
  no-boundary proposal (Phys. Rev. D 28, 2960, 1983), in which the early universe is a Euclidean
  four-sphere joined to a Lorentzian spacetime: the "any point reachable" region gives way to the
  region with a light cone.
- The Dirac equation: under the rotation the counter-rotating factors e^(∓iEt) of the positive-
  and negative-energy components become e^(∓Eτ), one decaying and one growing. The
  zitterbewegung beat becomes an exponential. This is the hook for JB's two-spinor picture: what
  reads as two spinors rotating in opposite senses in real time reads as two exponentials in
  imaginary time.

## What it does not give

- The arrow of time. The signature fixes the cone and separates its two halves, but nothing in it
  says which half is the future; that is thermodynamic, or a boundary condition.
- The number of time directions. Why one direction carries the opposite sign, and not none or two,
  is not answered by the rotation; two timelike directions give no well-posed dynamics, and none
  gives no cone.
- A physical process. The rotation is a mathematical continuation. Whether a signature change is
  a physical event (as in the no-boundary picture) or a calculational device is an open question
  that the record's 4D-to-3D+spacetime paper in the closed Dirac–Kuramoto repository would be the
  place to take up. That paper was not reopened for this note.

## Addendum, 2026-09-15: reflection positivity, or the modulus as a fold

Added at JB's request after his question of that night: could Hilbert space be more like a
four-dimensional space, with a measurement the modulus taken in three dimensions plus time.
The literal version fails on counting: Hilbert space is the space of functions on ordinary space
for one particle and on 3N-dimensional configuration space for N, and that growth with N is what
holds entanglement, which no fixed four-space can. But the question has a precise version, and
it is a theorem, in the Euclidean four-space of this note.

- Build the theory in Euclidean signature, with the real positive weight e^(−S_E/ħ) of the
  section above. Choose a three-dimensional hyperplane τ = 0 and call it "now." A state is a
  configuration of the field on one side of it, τ > 0.
- The inner product of a state with itself is the Euclidean weight of that configuration glued
  to its own mirror image across the plane, τ → −τ, with complex conjugation. The modulus squared
  is therefore what a four-dimensional object gives when it is folded onto a three-dimensional
  slice. The condition that this fold gives a non-negative number for every configuration is
  *reflection positivity*, and it is the axiom that lets a Hilbert space, a Hamiltonian bounded
  below, and the Minkowski theory be reconstructed from the Euclidean one (Osterwalder and
  Schrader, Commun. Math. Phys. 31, 83, 1973, and 42, 281, 1975 [verify]).
- So the square of the norm has a geometric origin: it is the reflection. This is more than the
  rotation gave for the exponential above, where the shape came out but not the rate.

What it does not give, in the same spirit as the section above: it delivers the square as the
positive quantity, not the Born rule as a probability of one outcome. It says why the norm is a
modulus squared. It does not say why a measurement yields one result with that frequency, nor
anything about the one-quantum stop across sites (clause 5 of `born_selection_revised/`). The
geometry supplies the form of the square; the event is still the postulate.

Pointers: Osterwalder and Schrader, above; Glimm and Jaffe, *Quantum Physics: A Functional
Integral Point of View*, 2nd ed., Springer 1987, the chapter on reflection positivity [verify];
for a general reader, the "imaginary time" chapters of Hawking already listed below.

## Pointers for the revisit

- Wick 1954, the original; any QFT text's chapter on Euclidean methods.
- Hartle and Hawking 1983; Hawking, *A Brief History of Time*, the "imaginary time" chapters,
  which are this idea for a general reader.
- In this folder: the record's zitterbewegung beat 2⟨E⟩ (second study, Part C) is the quantity
  that turns into a decay rate under the rotation.
