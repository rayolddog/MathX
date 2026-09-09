#!/usr/bin/env python3
"""pair_options.py: JB's two options for the Bell violation (addendum to PREREGISTRATION.md)."""
import csv, sys, time
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, str(Path(__file__).resolve().parent))
import pair_case as pc
from pair_case import N, A, AP, B, BP, THETAS, SLOT, TEXT, TEXT2, OUT

rng = np.random.default_rng(20260913); pc.rng = rng

def read(phase, axis, rule):
    d = phase - axis
    if rule == "least":
        return np.where(np.cos(2 * d) >= 0, 1, -1)
    return np.where(rng.uniform(0, 1, phase.size) < np.cos(d) ** 2, 1, -1)

def v1(a, b, rule, mu_mode="fixed", opposite=False):
    mu = np.full(N, 0.7) if mu_mode == "fixed" else rng.uniform(0, np.pi, N)
    a1, a2 = rng.uniform(0, np.pi, N), rng.uniform(0, np.pi, N)          # first-spinor phases, irrelevant after the lock
    delta = rng.uniform(0, np.pi, N)
    d2 = -delta if opposite else delta
    p1 = mu + delta; p2 = mu + d2                                          # second spinors after the equal rotations
    return np.mean(read(p1, a, rule) * read(p2, b, rule))

def v2(a, b, rule_b):
    phi = rng.uniform(0, np.pi, N)
    sA = read(phi, a, "least")                                             # particle 1 read at A (D3, deterministic reading)
    p2 = np.where(sA == 1, a, a + np.pi / 2)                               # particle 2 set to A's result
    return np.mean(sA * read(p2, b, rule_b))

def chsh(fn, *args):
    return fn(A, B, *args) - fn(A, BP, *args) + fn(AP, B, *args) + fn(AP, BP, *args)

def main():
    t0 = time.time(); lines = ["Verification, JB's two options (PREREGISTRATION.md, addendum)", ""]
    se_E = 1 / np.sqrt(N); se_S = 2 * se_E; rows = []
    cases = [("option 2, least, μ fixed", lambda a, b: v1(a, b, "least"), lambda t: 1 - 4 * t / np.pi, 2.0),
             ("option 2, Malus, μ fixed", lambda a, b: v1(a, b, "malus"), lambda t: 0.5 * np.cos(2 * t), np.sqrt(2)),
             ("option 2, least, μ drawn per pair", lambda a, b: v1(a, b, "least", "drawn"), lambda t: 1 - 4 * t / np.pi, 2.0),
             ("option 2, least, opposite offsets", lambda a, b: v1(a, b, "least", "fixed", True), lambda t: -(1 - 4 * t / np.pi), -2.0),
             ("option 1, Malus at B", lambda a, b: v2(a, b, "malus"), lambda t: np.cos(2 * t), 2 * np.sqrt(2)),
             ("option 1, least at B", lambda a, b: v2(a, b, "least"), lambda t: np.sign(np.cos(2 * t)), 4.0)]
    curves = {}
    for name, fn, closed, s_pred in cases:
        Es = [fn(0.0, t) for t in THETAS]; curves[name] = Es
        ok_E = all(abs(e - closed(t)) < 3 * se_E + (1e-12 if "least at B" in name else 0) for e, t in zip(Es, THETAS) if abs(np.cos(2 * t)) > 1e-9)
        S = chsh(fn); ok_S = abs(S - s_pred) < 3 * se_S + 1e-12
        rows.append(dict(case=name, S=S, S_predicted=s_pred, E_22_5=Es[3], E_22_5_predicted=closed(THETAS[3])))
        lines.append(f"{name}: E(22.5°) = {Es[3]:+.4f} against {closed(THETAS[3]):+.4f}; S = {S:+.4f} against {s_pred:+.4f} -> {'PASS' if ok_E and ok_S else 'FAIL'}")
    # diagnostic, not registered (added 2026-09-08 after the run): the "opposite offsets" prediction was written for
    # orthogonal polarizations (particle 2 at particle 1's angle + 90°), not for the mirror image about μ that D7 defines.
    def v1_orth(a, b):
        mu = np.full(N, 0.7); delta = rng.uniform(0, np.pi, N); return np.mean(read(mu + delta, a, "least") * read(mu + delta + np.pi / 2, b, "least"))
    Eo = v1_orth(0.0, THETAS[3]); So = chsh(v1_orth)
    lines.append(f"    diagnostic, not registered (added after the run): option 2 with particle 2 orthogonal to particle 1 (angle + 90°), least: E(22.5°) = {Eo:+.4f} (−0.5 for the sign flip), S = {So:+.4f} (−2)")
    with open(OUT / "options.csv", "w", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=list(rows[0].keys())); wr.writeheader()
        for r in rows: wr.writerow({k: (f"{v:.5g}" if isinstance(v, float) else v) for k, v in r.items()})
    fig, ax = plt.subplots(figsize=(6.2, 3.9)); td = np.degrees(THETAS)
    ax.plot(td, curves["option 2, least, μ fixed"], "o-", color=SLOT[0], label="option 2 (synchronized bulk), least rotation: Bell's sawtooth, S = 2.00")
    ax.plot(td, curves["option 2, Malus, μ fixed"], "s-", color=SLOT[1], label="option 2, Malus: ½ cos 2θ, S = 1.41")
    ax.plot(td, curves["option 1, Malus at B"], "^-", color=TEXT, label="option 1 (particle 2 set by A's outcome), Malus: cos 2θ, S = 2.83")
    ax.plot(td, curves["option 1, least at B"], "v--", color=SLOT[2], label="option 1, least rotation at B: step, S = 4.00")
    ax.set_xlabel("analyzer difference θ (degrees)"); ax.set_ylabel("correlation E(θ)"); ax.set_title("JB's two options for the violation"); ax.legend(fontsize=7, loc="lower left")
    fig.tight_layout(); fig.savefig(OUT / "fig2_options.png", dpi=150); plt.close(fig)
    text = "\n".join(lines) + "\n"; (OUT / "verification_options.txt").write_text(text); print(text)
    ok = all("FAIL" not in l for l in lines)
    with open(OUT / "run_info.txt", "a") as f:
        f.write(f"pair_options.py, run {time.strftime('%Y-%m-%d %H:%M:%S')}, seed 20260913, {time.time() - t0:.1f} s\nverification (options): {'PASS' if ok else 'FAIL'}\n")
    return 0 if ok else 1

if __name__ == "__main__":
    sys.exit(main())
