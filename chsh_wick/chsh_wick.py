#!/usr/bin/env python3
"""chsh_wick.py: a CHSH test of "shared complex time" for an entangled electron pair.

The idea under test (JB, 2026-11-17): the two electrons of an entangled pair share
a complex time.  A detector couples its electron to the bath (bulk) time, and that
coupling rotates the complex time.  The spin phase being measured is tied to the
complex time, so the rotation decides the outcome.

The shared complex time is modelled as a unit complex number per pair,

    z = exp(i * phi),      phi drawn uniformly from 0 to 2 pi.

Electron 1 has its spin along the angle phi, electron 2 along phi + pi (opposite,
as in the singlet state).  A detector whose analyzer points along the angle
"setting" rotates the time by that angle, z -> z * exp(-i * setting), and reads
the outcome (+1 or -1) from what is left.

Everything depends on ONE question: when detector A rotates the time, does
electron 2 feel it?

    Model L (local):   no.  Each detector rotates only its own electron's copy.
    Model N (nonlocal): yes.  The time is truly shared, so A's rotation and A's
                        outcome reset electron 2 before B reads it.
    Model Q:           the quantum prediction for the singlet, for comparison.

Run:  python3 chsh_wick.py        (about 2 seconds)
"""
import numpy as np

N_PAIRS = 400_000                 # pairs per pair of settings
SEED = 20261117
rng = np.random.default_rng(SEED)

# Standard CHSH analyzer angles for spin 1/2 (radians).
A1 = 0.0
A2 = np.pi / 2
B1 = np.pi / 4
B2 = 3 * np.pi / 4


# ---------------------------------------------------------------------------
# The detector: rotate the complex time, then read the outcome.
# ---------------------------------------------------------------------------
def rotate_time(z, setting):
    """The detector's coupling: rotate the complex time by the analyzer angle."""
    return z * np.exp(-1j * setting)


def read_deterministic(z):
    """Outcome fixed by the rotated time: +1 if it points 'up' the analyzer."""
    return np.where(z.real >= 0.0, 1, -1)


def read_born(z):
    """Outcome drawn at random with the spin-1/2 Born probability.

    If the spin makes the angle d with the analyzer, P(+1) = cos^2(d / 2).
    Here d is the phase of the rotated time.
    """
    d = np.angle(z)
    p_plus = np.cos(d / 2.0) ** 2
    return np.where(rng.uniform(0.0, 1.0, z.size) < p_plus, 1, -1)


# ---------------------------------------------------------------------------
# Model L: the rotation at each detector stays at that detector.
# ---------------------------------------------------------------------------
def model_local(a, b, read):
    phi = rng.uniform(0.0, 2.0 * np.pi, N_PAIRS)
    z1 = np.exp(1j * phi)                 # electron 1's copy of the shared time
    z2 = np.exp(1j * (phi + np.pi))       # electron 2's copy: spin opposite

    out_a = read(rotate_time(z1, a))      # uses only z1 and the setting a
    out_b = read(rotate_time(z2, b))      # uses only z2 and the setting b
    return out_a, out_b


# ---------------------------------------------------------------------------
# Model N: the time is one shared object, so A's rotation reaches electron 2.
# ---------------------------------------------------------------------------
def model_nonlocal(a, b):
    phi = rng.uniform(0.0, 2.0 * np.pi, N_PAIRS)
    z1 = np.exp(1j * phi)

    out_a = read_born(rotate_time(z1, a))

    # The step that crosses the gap: electron 2 is now set by A's SETTING and
    # A's OUTCOME.  If A read +1 along a, electron 2 points along a + pi.
    # If A read -1, electron 2 points along a.
    angle_2 = np.where(out_a == 1, a + np.pi, a)
    z2 = np.exp(1j * angle_2)

    out_b = read_born(rotate_time(z2, b))
    return out_a, out_b


# ---------------------------------------------------------------------------
# Model Q: sample straight from the quantum probabilities for the singlet.
# ---------------------------------------------------------------------------
def model_quantum(a, b):
    out_a = np.where(rng.uniform(0.0, 1.0, N_PAIRS) < 0.5, 1, -1)
    p_same = np.sin((a - b) / 2.0) ** 2          # P(both outcomes equal)
    same = rng.uniform(0.0, 1.0, N_PAIRS) < p_same
    out_b = np.where(same, out_a, -out_a)
    return out_a, out_b


# ---------------------------------------------------------------------------
# The CHSH quantity.
# ---------------------------------------------------------------------------
def correlation(model, a, b, *args):
    out_a, out_b = model(a, b, *args)
    return np.mean(out_a * out_b)


def chsh(model, *args):
    """S = E(a1,b1) - E(a1,b2) + E(a2,b1) + E(a2,b2)."""
    e11 = correlation(model, A1, B1, *args)
    e12 = correlation(model, A1, B2, *args)
    e21 = correlation(model, A2, B1, *args)
    e22 = correlation(model, A2, B2, *args)
    return e11 - e12 + e21 + e22, e11


def bell_bound_by_counting():
    """The CHSH proof in 16 lines.

    In ANY local model a pair leaves the source with four answers already fixed:
    what A would read at a1 and at a2, and what B would read at b1 and at b2.
    List all 16 possibilities.  Each gives exactly +2 or -2, so no mixture of
    them, whatever the shared variable is, can average to more than 2.
    """
    values = []
    for a1 in (+1, -1):
        for a2 in (+1, -1):
            for b1 in (+1, -1):
                for b2 in (+1, -1):
                    values.append(a1 * b1 - a1 * b2 + a2 * b1 + a2 * b2)
    return values


def main():
    se = 2.0 / np.sqrt(N_PAIRS)           # one standard error on S

    print(f"CHSH test, {N_PAIRS} pairs per setting pair, seed {SEED}")
    print(f"settings: a1 = 0, a2 = 90, b1 = 45, b2 = 135 degrees;  s.e. of S = {se:.3f}")
    print()
    print(f"{'model':58s} {'E(0,45)':>8s} {'|S|':>6s}")
    print("-" * 74)

    rows = [
        ("L  local, outcome fixed by the rotated time", model_local, (read_deterministic,)),
        ("L  local, outcome drawn by the Born rule", model_local, (read_born,)),
        ("N  nonlocal, A's rotation resets electron 2", model_nonlocal, ()),
        ("Q  quantum mechanics, singlet state", model_quantum, ()),
    ]
    for name, model, args in rows:
        s, e11 = chsh(model, *args)
        print(f"{name:58s} {e11:8.3f} {abs(s):6.3f}")

    print("-" * 74)
    print(f"{'Bell bound for every local model':58s} {'':8s} {2.0:6.3f}")
    print(f"{'quantum limit, 2 * sqrt(2)':58s} {'':8s} {2.0 * np.sqrt(2.0):6.3f}")
    print()

    # No-signalling check for model N: B's own statistics must not reveal a.
    _, out_b = model_nonlocal(A1, B1)
    p1 = np.mean(out_b == 1)
    _, out_b = model_nonlocal(A2, B1)
    p2 = np.mean(out_b == 1)
    print("Model N, fraction of +1 at B (b = 45 degrees):")
    print(f"    with A set to  0 degrees: {p1:.4f}")
    print(f"    with A set to 90 degrees: {p2:.4f}     (equal, so B cannot tell what A chose)")
    print()

    values = bell_bound_by_counting()
    print("All 16 sets of pre-fixed answers (a1, a2, b1, b2) give S =")
    print("   ", " ".join(f"{v:+d}" for v in values))


if __name__ == "__main__":
    main()
