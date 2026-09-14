#!/usr/bin/env python3
"""
staggered_conditional.py: reading C′ (PREREGISTRATION.md addendum 2b, D17): the race's rates on the
normalized conditional state. Branch weights a_P²(t) = s_P exp(-λ W_P τ_P(t)); the first fire occurs
where the survival norm n(t) crosses a uniform draw; the port is drawn in proportion to the exposed
rates at that time. Reading S is taken from staggered.csv for the figure.
"""
import csv, time
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent; OUT = HERE / "results"
SEED = 20260917
rng = np.random.default_rng(SEED)
N, T, TRIALS, K = 10, 1.0, 200000, 1000  # grid reduced from the registered 4000 for runtime; noted in RESULTS.md
DELTAS = [0.0, 0.1, 0.5, 1.0, 2.0]
KAPPAS = [1.0, 10.0, 100.0]
SURFACE, TEXT, TEXT2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e6e5e1"
SLOT = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#7a5cc7", "#8c8c8c"]
plt.rcParams.update({"font.family": "sans-serif", "font.size": 9, "axes.titlesize": 10, "axes.titleweight": "bold",
    "axes.edgecolor": GRID, "axes.labelcolor": TEXT2, "xtick.color": TEXT2, "ytick.color": TEXT2, "axes.titlecolor": TEXT,
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE, "axes.grid": True, "grid.color": GRID,
    "grid.linewidth": 0.8, "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False, "legend.fontsize": 8, "axes.axisbelow": True})

def weights(phase, n):
    return (np.cos(rng.uniform(0, 2 * np.pi, (n, N)) / 2) ** 2).sum(1) if phase else np.full(n, float(N))

def run(r, kappa, delta, phase, chunk=5000):
    lam = kappa / (N * T); sA, sB = r / (r + 1), 1 / (r + 1)
    t = np.linspace(0, delta + T, K + 1)
    tauA, tauB = np.clip(t, 0, T), np.clip(t - delta, 0, T)
    inA, inB = (t <= T), (t >= delta)
    winA = winB = 0
    for start in range(0, TRIALS, chunk):
        c = min(chunk, TRIALS - start)
        WA, WB = weights(phase, c), weights(phase, c)
        aA = sA * np.exp(-lam * WA[:, None] * tauA[None, :]); aB = sB * np.exp(-lam * WB[:, None] * tauB[None, :])
        n = aA + aB
        u = rng.uniform(0, 1, c)
        crossed = n <= u[:, None]
        fires = crossed.any(1)
        k = np.where(fires, crossed.argmax(1), 0)
        idx = np.arange(c)
        hA = lam * WA * aA[idx, k] * inA[k]; hB = lam * WB * aB[idx, k] * inB[k]
        pA = np.where(hA + hB > 0, hA / np.maximum(hA + hB, 1e-300), 0.0)
        a = fires & (rng.uniform(0, 1, c) < pA); b = fires & ~a
        winA += a.sum(); winB += b.sum()
    PA, PB = winA / TRIALS, winB / TRIALS
    return PA, PB, 1 - PA - PB, PA / (PA + PB)

def main():
    t0 = time.time(); lines = ["Verification (PREREGISTRATION.md addendum 2b)", ""]; rows = []
    for r in (4, 2):
        for phase in (False, True):
            for kappa in KAPPAS:
                for delta in DELTAS:
                    PA, PB, P0, share = run(r, kappa, delta, phase)
                    se = np.sqrt(share * (1 - share) / max((PA + PB) * TRIALS, 1)); born = r / (r + 1)
                    rows.append(dict(reading="Cprime", r=r, phase_weight=phase, kappa=kappa, delta_over_T=delta, P_A=PA, P_B=PB, P_none=P0, share=share, se=se, born=born))
                    if not phase and kappa == 10 and r == 4:
                        ok = abs(share - born) < 3 * se and (PA + PB) >= 0.999
                        lines.append(f"C', no phase, kappa = 10, r = 4, Δ/T = {delta:3.1f}: share {share:.4f} against 0.800 (s.e. {se:.4f}), P(fired) = {PA+PB:.4f} >= 0.999 -> {'PASS' if ok else 'FAIL'}")
                    if not phase and kappa == 1 and r == 4 and delta == 0.0:
                        pf = PA + PB; sef = np.sqrt(pf * (1 - pf) / TRIALS); ok = abs(pf - (1 - np.exp(-1))) < 3 * sef
                        lines.append(f"C', no phase, kappa = 1: P(fired) = {pf:.4f} against 1 - e^-1 = {1-np.exp(-1):.4f} (s.e. {sef:.4f}) -> {'PASS' if ok else 'FAIL'}")
                    if not phase and kappa == 10 and r == 2 and delta == 1.0:
                        ok = abs(share - born) < 3 * se
                        lines.append(f"C', no phase, kappa = 10, r = 2, Δ = T: share {share:.4f} against 0.667 (s.e. {se:.4f}) -> {'PASS' if ok else 'FAIL'}")
    # diagnostic, not registered: the addendum-2 reading-S cell that failed by 3.0 s.e., re-drawn once with a fresh seed
    lam = 100.0 / (N * T); tA = rng.exponential(1 / (lam * 4 * N), TRIALS); tB = rng.exponential(1 / (lam * N), TRIALS)
    tA = np.where(tA < T, tA, np.inf); tB = np.where(tB < T, tB, np.inf); a = tA < tB; b = tB < tA
    sh = a.mean() / (a.mean() + b.mean()); se = np.sqrt(sh * (1 - sh) / TRIALS)
    lines.append(f"    diagnostic, not registered (P18): reading S, kappa = 100, Δ = 0, fresh seed: share {sh:.4f} against 0.8000 (s.e. {se:.4f}) -> {'within 3 s.e.' if abs(sh - 0.8) < 3 * se else 'outside 3 s.e.'}")
    fmt = lambda v: (f"{v:.5g}" if isinstance(v, float) else v)
    with open(OUT / "staggered_conditional.csv", "w", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=list(rows[0].keys())); wr.writeheader()
        for row in rows: wr.writerow({k: fmt(v) for k, v in row.items()})
    S = list(csv.DictReader(open(OUT / "staggered.csv")))
    fig, axes = plt.subplots(1, 2, figsize=(9.8, 3.9))
    for ax, phase, title in ((axes[0], False, "no phase weight"), (axes[1], True, "phase weight cos²(Δφ/2), fresh")):
        for k, kappa in enumerate(KAPPAS):
            sS = [x for x in S if x["reading"] == "S" and x["r"] == "4" and x["phase_weight"] == str(phase) and float(x["kappa"]) == kappa]
            ax.plot([float(x["delta_over_T"]) for x in sS], [float(x["share"]) for x in sS], "o-", color=SLOT[k], label=f"survival hazard, κ = {kappa:g}")
            sC = [x for x in rows if x["r"] == 4 and x["phase_weight"] == phase and x["kappa"] == kappa]
            ax.plot([x["delta_over_T"] for x in sC], [x["share"] for x in sC], "s--", color=SLOT[k], label=f"normalized conditional state, κ = {kappa:g}")
        ax.axhline(0.8, color=TEXT, lw=0.9); ax.text(1.3, 0.81, "Born 0.800", fontsize=8)
        ax.set_xlabel("arrival stagger Δ/T (packet lengths)"); ax.set_ylabel("share of records at the nearer, brighter pixel (r = 4)"); ax.set_ylim(0.75, 1.02)
        ax.set_title(f"Staggered arrival, {title}"); ax.legend(fontsize=7, loc="center right")
    fig.tight_layout(); fig.savefig(OUT / "fig3_staggered.png", dpi=150); plt.close(fig)
    text = "\n".join(lines) + "\n"; (OUT / "verification_staggered_conditional.txt").write_text(text); print(text)
    print("shares, C', r = 4:")
    for phase in (False, True):
        for kappa in KAPPAS:
            sel = [x for x in rows if x["r"] == 4 and x["phase_weight"] == phase and x["kappa"] == kappa]
            print(f"  phase={phase!s:5} kappa={kappa:5.0f}: " + "  ".join(f"Δ/T={x['delta_over_T']:.1f}: {x['share']:.3f} (fired {x['P_A']+x['P_B']:.3f})" for x in sel))
    ok = all("FAIL" not in l for l in lines)
    with open(OUT / "run_info.txt", "a") as f:
        f.write(f"staggered_conditional.py, run {time.strftime('%Y-%m-%d %H:%M:%S')}, seed {SEED}, {time.time() - t0:.1f} s, numpy {np.__version__}\nverification (staggered, conditional): {'PASS' if ok else 'FAIL'}\n")
    print(f"{time.time() - t0:.1f} s")

if __name__ == "__main__":
    main()
