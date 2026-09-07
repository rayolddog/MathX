#!/usr/bin/env python3
"""
energy_cost.py: Part E, a small energy cost of the phase change (PREREGISTRATION.md §6, addendum 3)
===================================================================================================

Rule R2. The particle's deposit at atom i is  E_i = max(0, A_i w(Δφ_i) − c(Δ_i)/2)  with
w = cos²(Δφ/2) the second-spinor interference weight and c(Δ) = η (1 − cos Δ)/2 the cost of the
first-spinor shift Δ, shared equally (the particle pays c/2). Two readings of the bulk's first
spinors: random (Δ_i uniform, primary) and aligned (one Δ for all atoms, swept).
The large-n share is the ratio of the lobes' mean truncated deposits, computed by quadrature and
checked against Monte Carlo in results/verification.txt.
"""
import csv, sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
OUT = HERE / "results"
rng = np.random.default_rng(20260907 + 5)
R_INT, EPS = 4.0, 0.35
ETA_LIST = [0.0, 0.05, 0.1, 0.2, 0.5]
N_LIST = [1, 2, 5, 10, 100, 1000]
DELTA_DEG = [0, 45, 90, 135, 180]

def cost(delta, eta):
    return eta * (1 - np.cos(delta)) / 2

def mean_deposit_quadrature(A, eta, delta_fixed=None, m=1500):
    u = (np.arange(m) + 0.5) / m * 2 * np.pi - np.pi
    w = np.cos(u / 2) ** 2
    if delta_fixed is None:
        c = cost(u, eta)
        return np.mean(np.maximum(0, A * w[:, None] - c[None, :] / 2))
    return np.mean(np.maximum(0, A * w - cost(delta_fixed, eta) / 2))

def share(E1, E2):
    S1, S2 = E1.sum(axis=1), E2.sum(axis=1); live = (S1 + S2) > 0
    return np.mean(S1[live] / (S1[live] + S2[live])), 1 - live.mean()

def run_random(lines):
    rows = []; worst = 0.0
    for eta in ETA_LIST:
        m1, m2 = mean_deposit_quadrature(R_INT, eta), mean_deposit_quadrature(1.0, eta)
        P_inf = m1 / (m1 + m2)
        # the fired atom's mismatch, 50 atoms, uniform wave A = 1
        M, N = 200000, 50
        dphi = rng.uniform(-np.pi, np.pi, (M, N)); dlt = rng.uniform(-np.pi, np.pi, (M, N))
        E = np.maximum(0, np.cos(dphi / 2) ** 2 - cost(dlt, eta) / 2)
        cdf = np.cumsum(E, axis=1); u = rng.uniform(0, 1, M) * cdf[:, -1]
        k = (cdf < u[:, None]).sum(axis=1); mm = np.abs(dphi[np.arange(M), k])
        for n in N_LIST:
            Mn = int(min(400000, 4e6 // n))
            E1 = np.maximum(0, R_INT * np.cos(rng.uniform(-np.pi, np.pi, (Mn, n)) / 2) ** 2 - cost(rng.uniform(-np.pi, np.pi, (Mn, n)), eta) / 2)
            E2 = np.maximum(0, np.cos(rng.uniform(-np.pi, np.pi, (Mn, n)) / 2) ** 2 - cost(rng.uniform(-np.pi, np.pi, (Mn, n)), eta) / 2)
            P1, miss = share(E1, E2)
            if n == 1000:
                worst = max(worst, abs(P1 - P_inf))
            rows.append([eta, n, f"{P1:.4f}", f"{miss:.4f}", f"{P_inf:.4f}", f"{(R_INT - eta/2)/(R_INT + 1 - eta):.4f}",
                         f"{np.log(P_inf/(1-P_inf))/np.log(2):.3f}", f"{np.degrees(mm.mean()):.2f}", f"{(mm <= EPS).mean():.4f}"])
    with open(OUT / "part_e_cost.csv", "w", newline="") as f:
        wr = csv.writer(f)
        wr.writerow(["eta", "n_per_lobe", "R2_P1_r4", "miss_rate", "P1_large_n_quadrature", "P1_large_n_no_floor",
                     "gamma_large_n", "fired_mismatch_mean_deg_N50", "fired_within_band_N50"]); wr.writerows(rows)
    lines.append(f"Part E, random first spinors: worst |P₁(n = 1000) − quadrature large-n| over η: {worst:.2e}   (tolerance 0.01)   -> {'PASS' if worst < 0.01 else 'FAIL'}")
    return rows, worst < 0.01

def run_aligned(lines):
    rows = []; worst = 0.0
    for eta in (0.1, 0.5):
        for D in DELTA_DEG:
            d = np.radians(D)
            m1, m2 = mean_deposit_quadrature(R_INT, eta, d), mean_deposit_quadrature(1.0, eta, d)
            P_inf = m1 / (m1 + m2); c = cost(d, eta)
            for n in (10, 1000):
                Mn = int(min(400000, 4e6 // n))
                E1 = np.maximum(0, R_INT * np.cos(rng.uniform(-np.pi, np.pi, (Mn, n)) / 2) ** 2 - c / 2)
                E2 = np.maximum(0, np.cos(rng.uniform(-np.pi, np.pi, (Mn, n)) / 2) ** 2 - c / 2)
                P1, miss = share(E1, E2)
                if n == 1000:
                    worst = max(worst, abs(P1 - P_inf))
                rows.append([eta, D, n, f"{P1:.4f}", f"{miss:.4f}", f"{P_inf:.4f}", f"{(R_INT - c)/(R_INT + 1 - 2*c):.4f}"])
    with open(OUT / "part_e_cost_aligned.csv", "w", newline="") as f:
        wr = csv.writer(f)
        wr.writerow(["eta", "Delta_deg", "n_per_lobe", "R2_P1_r4", "miss_rate", "P1_large_n_quadrature", "P1_large_n_no_floor"]); wr.writerows(rows)
    lines.append(f"Part E, aligned first spinors: worst |P₁(n = 1000) − quadrature large-n| over η and Δ: {worst:.2e}   (tolerance 0.01)   -> {'PASS' if worst < 0.01 else 'FAIL'}")
    return rows, worst < 0.01

SURFACE, TEXT, TEXT2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e6e5e1"
SLOT = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#7a5cc7"]
plt.rcParams.update({
    "font.family": "sans-serif", "font.size": 9, "axes.titlesize": 10, "axes.titleweight": "bold",
    "axes.edgecolor": GRID, "axes.labelcolor": TEXT2, "xtick.color": TEXT2, "ytick.color": TEXT2,
    "axes.titlecolor": TEXT, "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
    "savefig.facecolor": SURFACE, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.8,
    "grid.linestyle": "-", "axes.spines.top": False, "axes.spines.right": False,
    "legend.frameon": False, "legend.fontsize": 8, "axes.axisbelow": True,
})

def figure(rows_r, rows_a):
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.0))
    ax = axes[0]
    ax.axhline(0.8, color=TEXT2, lw=0.8, alpha=0.6); ax.text(1.05, 0.803, "Born 0.80", color=TEXT2, fontsize=7.5)
    for eta, c in zip(ETA_LIST, SLOT):
        sub = [r for r in rows_r if r[0] == eta]
        ax.plot([r[1] for r in sub], [float(r[2]) for r in sub], "o-", color=c, lw=1.6, ms=4, label=f"η = {eta:g}")
        ax.axhline(float(sub[0][4]), color=c, lw=0.8, ls="--", alpha=0.6)
    ax.set_xscale("log"); ax.set_ylim(0.68, 0.86); ax.set_xlabel("n atoms per lobe"); ax.set_ylabel("P(brighter lobe fires), r = 4, rule R2")
    ax.set_title("Cost of the shift, random first spinors (dashed: n → ∞)"); ax.legend(loc="lower right", ncol=2)
    ax = axes[1]
    ax.axhline(0.8, color=TEXT2, lw=0.8, alpha=0.6)
    for eta, c in zip((0.1, 0.5), SLOT[1:3]):
        sub = [r for r in rows_a if r[0] == eta and r[2] == 1000]
        ax.plot([r[1] for r in sub], [float(r[3]) for r in sub], "o-", color=c, lw=1.6, ms=4, label=f"η = {eta:g}, n = 1000")
        ax.plot([r[1] for r in sub], [float(r[6]) for r in sub], ":", color=c, lw=1.2, label=f"η = {eta:g}, no floor")
    ax.set_xticks(DELTA_DEG); ax.set_ylim(0.78, 0.9); ax.set_xlabel("particle's first-spinor shift Δ against an aligned bulk (deg)")
    ax.set_ylabel("P(brighter lobe fires)"); ax.set_title("Aligned first spinors: the particle's shift Δ enters")
    ax.legend(loc="upper left")
    fig.tight_layout(); fig.savefig(OUT / "fig7_cost.png", dpi=150); plt.close(fig)

def main():
    lines = ["Verification, Part E (PREREGISTRATION.md §6, addendum of 2026-09-07: the energy cost of the shift, rule R2)", ""]
    rows_r, ok_r = run_random(lines)
    rows_a, ok_a = run_aligned(lines)
    figure(rows_r, rows_a)
    text = "\n".join(lines) + "\n"
    with open(OUT / "verification.txt", "a") as f:
        f.write("\n" + text)
    with open(OUT / "run_info.txt", "a") as f:
        f.write(f"verification (energy_cost.py): {'PASS' if ok_r and ok_a else 'FAIL'}\n")
    print(text)
    print("Part E, random first spinors, r = 4, R2:")
    for eta in ETA_LIST:
        sub = [r for r in rows_r if r[0] == eta]
        print(f"  η = {eta:g}: " + ", ".join(f"n={r[1]}: {float(r[2]):.3f}" for r in sub) + f";  n→∞ {float(sub[0][4]):.4f} (no floor {float(sub[0][5]):.4f}), γ {sub[0][6]};"
              f" miss at n=1 {float(sub[0][3]):.3f}; fired mismatch {float(sub[0][7]):.1f}°, in band {float(sub[0][8]):.3f}")
    print("Part E, aligned first spinors, n = 1000:")
    for eta in (0.1, 0.5):
        print(f"  η = {eta:g}: " + ", ".join(f"Δ={r[1]}°: {float(r[3]):.4f}" for r in rows_a if r[0] == eta and r[2] == 1000))
    return 0 if ok_r and ok_a else 1

if __name__ == "__main__":
    sys.exit(main())
