# Provenance of the copied scripts

Copied unchanged on 2026-09-06 from `Projects/Physics/DiracKuramotoFramework` at commit
`c9bbf12` (the closed record). They are here so the new exploration can be read without
the other repository open. Do not edit them; adapt in `../kuramoto_one_two_rotations/`.

| File | Source path | What it did |
|---|---|---|
| `kuramoto_sync.py` | `tests/kuramoto_sync.py` | Single clock, dφ/dt = ω + K sin(Φ_bulk − φ), free flight then pull to the bulk at detection (K = 15, ω = 2π). The equation and the RK45 pattern reused by `kuramoto_rotations.py`. |
| `dirac_arnold_tongue.py` | `ZBW_synchronization/dirac_arnold_tongue.py` | Two confined 1+1-D Dirac particles with a Breit-like coupling, closed system; order parameter of the zitterbewegung phase difference over a (coupling, detuning) grid. Writes to `/Users/john-bramble/` paths as written; not run here. |
| `arnold_tongue.log` | `ZBW_synchronization/arnold_tongue.log` | The output of that run (R matrix). |

The closed-system script is kept because it is the record's own evidence that the L–R
coupling in the Dirac equation gives coherent precession, not an Adler lock (README,
"The central result"); the new exploration's Adler form for the lower component is a
posit that this script does not support.
