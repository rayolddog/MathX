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

## Pointers for the revisit

- Wick 1954, the original; any QFT text's chapter on Euclidean methods.
- Hartle and Hawking 1983; Hawking, *A Brief History of Time*, the "imaginary time" chapters,
  which are this idea for a general reader.
- In this folder: the record's zitterbewegung beat 2⟨E⟩ (second study, Part C) is the quantity
  that turns into a decay rate under the rotation.
