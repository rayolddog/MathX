#!/usr/bin/env python3
"""
interference_visibility.py: Part D, the interference term (PREREGISTRATION.md §6, addendum of 2026-09-07)
=========================================================================================================

Record rule R2 throughout (JB's choice). Two readings of the interference between the incident
particle's second spinor and the bound electron's:

  D10, gating          E_i ∝ A_i (1 + V cos Δφ_i)/2,   V = 2ab/(a² + b²) ≤ 1;  V = 1 is A1 as run
  D11, energy density  E_i ∝ A_i + β² + 2β √A_i cos Δφ_i,   β = b/a₀

Two lobes of intensity r and 1 (n atoms each); the fired atom's mismatch; the ordered-bulk two-half
law with visibility. Closed forms are checked against Monte Carlo in results/verification.txt.
"""
import csv, sys
from pathlib import Path
import numpy as np
from scipy.special import i0, i1
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
OUT = HERE / "results"
rng = np.random.default_rng(20260907 + 4)
EPS, R_INT = 0.35, 4.0
V_LIST = [1.0, 0.5, 0.2, 0.0]
BETA_LIST = [0.0, 0.5, 1.0, 2.0]
N_LIST = [1, 2, 5, 10, 100, 1000]

def p1_r2(E1, E2):
    return np.mean(E1.sum(axis=1) / (E1.sum(axis=1) + E2.sum(axis=1)))

def run_d10(lines):
    rows = []; worst = 0.0
    for V in V_LIST:
        # the fired atom's mismatch, N = 50 atoms, uniform wave
        M, N = 200000, 50
        dphi = rng.uniform(-np.pi, np.pi, (M, N)); w = (1 + V * np.cos(dphi)) / 2
        cdf = np.cumsum(w, axis=1); u = rng.uniform(0, 1, M) * cdf[:, -1]
        k = (cdf < u[:, None]).sum(axis=1); m = np.abs(dphi[np.arange(M), k])
        mm_cf, band_cf = np.pi / 2 - 2 * V / np.pi, (EPS + V * np.sin(EPS)) / np.pi
        worst = max(worst, abs(m.mean() - mm_cf), abs((m <= EPS).mean() - band_cf))
        for n in N_LIST:
            Mn = int(min(400000, 4e6 // n))
            E1 = R_INT * (1 + V * np.cos(rng.uniform(-np.pi, np.pi, (Mn, n)))) / 2
            E2 = (1 + V * np.cos(rng.uniform(-np.pi, np.pi, (Mn, n)))) / 2
            P1 = p1_r2(E1, E2)
            if n == 1000:
                worst = max(worst, abs(P1 - R_INT / (1 + R_INT)))
            rows.append([V, n, f"{P1:.4f}", f"{R_INT/(1+R_INT):.4f}", f"{np.degrees(m.mean()):.2f}", f"{np.degrees(mm_cf):.2f}",
                         f"{(m <= EPS).mean():.4f}", f"{band_cf:.4f}"])
    # ordered bulk, two halves of 100 atoms, κ = 4, δ_p sweep
    kappa = 4.0; Rk = i1(kappa) / i0(kappa); orows = []
    for V in V_LIST:
        for dp_deg in (0, 45, 90, 135, 180):
            dp = np.radians(dp_deg); M, N = 20000, 100
            d1 = rng.vonmises(0.0, kappa, (M, N)); d2 = rng.vonmises(np.pi, kappa, (M, N))
            w1 = (1 + V * np.cos(dp - d1)) / 2; w2 = (1 + V * np.cos(dp - d2)) / 2
            P1 = np.mean(w1.sum(axis=1) / (w1.sum(axis=1) + w2.sum(axis=1)))
            cf = 0.5 * (1 + V * Rk * np.cos(dp)); worst = max(worst, abs(P1 - cf))
            orows.append([V, kappa, f"{V*Rk:.4f}", dp_deg, f"{P1:.4f}", f"{cf:.4f}"])
    with open(OUT / "part_d_visibility.csv", "w", newline="") as f:
        wr = csv.writer(f)
        wr.writerow(["V", "n_per_lobe", "R2_P1_r4", "born_P1", "fired_mismatch_mean_deg_mc", "fired_mismatch_mean_deg_closed",
                     "fired_within_band_mc", "fired_within_band_closed"]); wr.writerows(rows)
    with open(OUT / "part_d_visibility_ordered.csv", "w", newline="") as f:
        wr = csv.writer(f)
        wr.writerow(["V", "kappa", "visibility_V_R_kappa", "delta_p_deg", "two_half_R2_P1_mc", "closed_half_1_plus_VR_cos"]); wr.writerows(orows)
    lines.append(f"Part D, D10 (visibility): worst |MC − closed| over the fired atom's mean mismatch, its band fraction, P₁ at n = 1000 against Born,"
                 f" and the ordered two-half law: {worst:.2e}   (tolerance 0.01)   -> {'PASS' if worst < 0.01 else 'FAIL'}")
    return rows, orows, worst < 0.01

def run_d11(lines):
    rows = []; worst = 0.0
    for b in BETA_LIST:
        cf = (R_INT + b * b) / (R_INT + 1 + 2 * b * b)
        for n in N_LIST:
            Mn = int(min(400000, 4e6 // n))
            E1 = R_INT + b * b + 2 * b * np.sqrt(R_INT) * np.cos(rng.uniform(-np.pi, np.pi, (Mn, n)))
            E2 = 1 + b * b + 2 * b * np.cos(rng.uniform(-np.pi, np.pi, (Mn, n)))
            P1 = p1_r2(E1, E2)
            if n == 1000:
                worst = max(worst, abs(P1 - cf))
            rows.append([b, n, f"{P1:.4f}", f"{cf:.4f}", f"{R_INT/(1+R_INT):.4f}"])
    with open(OUT / "part_d_bound_amplitude.csv", "w", newline="") as f:
        wr = csv.writer(f)
        wr.writerow(["beta", "n_per_lobe", "R2_P1_r4", "closed_large_n", "born_P1"]); wr.writerows(rows)
    lines.append(f"Part D, D11 (bound amplitude): worst |P₁(n = 1000) − (r + β²)/(r + 1 + 2β²)| = {worst:.2e}   (tolerance 0.01)   -> {'PASS' if worst < 0.01 else 'FAIL'}")
    return rows, worst < 0.01

SURFACE, TEXT, TEXT2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e6e5e1"
SLOT = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]
plt.rcParams.update({
    "font.family": "sans-serif", "font.size": 9, "axes.titlesize": 10, "axes.titleweight": "bold",
    "axes.edgecolor": GRID, "axes.labelcolor": TEXT2, "xtick.color": TEXT2, "ytick.color": TEXT2,
    "axes.titlecolor": TEXT, "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
    "savefig.facecolor": SURFACE, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.8,
    "grid.linestyle": "-", "axes.spines.top": False, "axes.spines.right": False,
    "legend.frameon": False, "legend.fontsize": 8, "axes.axisbelow": True,
})

def figure(rows10, rows11):
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.0), sharey=True)
    ax = axes[0]
    ax.axhline(0.8, color=TEXT2, lw=0.8, alpha=0.6); ax.text(1.05, 0.805, "Born 0.80", color=TEXT2, fontsize=7.5)
    for V, c in zip(V_LIST, SLOT):
        sub = [r for r in rows10 if r[0] == V]
        ax.plot([r[1] for r in sub], [float(r[2]) for r in sub], "o-", color=c, lw=1.6, ms=4, label=f"V = {V:g}" + (" (as run)" if V == 1 else ""))
    ax.set_xscale("log"); ax.set_xlabel("n atoms per lobe"); ax.set_ylabel("P(brighter lobe fires), r = 4, rule R2")
    ax.set_title("D10: gating reading, visibility V"); ax.legend(loc="lower right")
    ax = axes[1]
    ax.axhline(0.8, color=TEXT2, lw=0.8, alpha=0.6)
    for b, c in zip(BETA_LIST, SLOT):
        sub = [r for r in rows11 if r[0] == b]
        ax.plot([r[1] for r in sub], [float(r[2]) for r in sub], "o-", color=c, lw=1.6, ms=4, label=f"β = {b:g}")
        ax.axhline(float(sub[0][3]), color=c, lw=0.8, ls="--", alpha=0.7)
    ax.set_xscale("log"); ax.set_xlabel("n atoms per lobe"); ax.set_ylim(0.5, 0.85)
    ax.set_title("D11: bound energy joins the deposit (dashed: n → ∞)")
    ax.legend(loc="lower right")
    fig.tight_layout(); fig.savefig(OUT / "fig6_interference.png", dpi=150); plt.close(fig)

def main():
    lines = ["Verification, Part D (PREREGISTRATION.md §6, addendum of 2026-09-07: the interference term, rule R2)", ""]
    rows10, orows, ok10 = run_d10(lines)
    rows11, ok11 = run_d11(lines)
    figure(rows10, rows11)
    text = "\n".join(lines) + "\n"
    with open(OUT / "verification.txt", "a") as f:
        f.write("\n" + text)
    with open(OUT / "run_info.txt", "a") as f:
        f.write(f"verification (interference_visibility.py): {'PASS' if ok10 and ok11 else 'FAIL'}\n")
    print(text)
    print("D10, r = 4, R2: P₁ by n for each V")
    for V in V_LIST:
        sub = [r for r in rows10 if r[0] == V]
        print(f"  V = {V:g}: " + ", ".join(f"n={r[1]}: {float(r[2]):.3f}" for r in sub) + f";  fired mismatch {float(sub[0][4]):.1f}° (closed {float(sub[0][5]):.1f}°), in band {float(sub[0][6]):.3f} (closed {float(sub[0][7]):.3f})")
    print("D11, r = 4, R2: P₁ by n for each β")
    for b in BETA_LIST:
        sub = [r for r in rows11 if r[0] == b]
        print(f"  β = {b:g}: " + ", ".join(f"n={r[1]}: {float(r[2]):.3f}" for r in sub) + f";  large-n closed {float(sub[0][3]):.3f}")
    print("Ordered bulk, κ = 4, two halves, δ_p = 0: " + ", ".join(f"V={r[0]:g}: {float(r[4]):.3f} (closed {float(r[5]):.3f})" for r in orows if r[3] == 0))
    return 0 if ok10 and ok11 else 1

if __name__ == "__main__":
    sys.exit(main())
