#!/usr/bin/env python3
"""
phase_race.py: the second-spinor phase race at a detector (Parts A, B and the ordered bulk)
=========================================================================================

Definitions are fixed in PREREGISTRATION.md (read it first).

  contact     first spinors matched at once; both spinors of the particle shifted together;
              mismatch of the second spinors at atom i:  Δφ_i = δ_p − δ_i
  weight      w(Δφ) = cos²(Δφ/2)  (A1);  sensitivity: band |Δφ| ≤ ε = 0.35
  deposit     E_i ∝ A_i w(Δφ_i)                                             (D4)
  record      R1 shared map; R2 one atom, probability ∝ E_i (primary); R3 the largest E_i   (D5)

Closed forms are checked against Monte Carlo in results/verification.txt.
"""
import csv, sys, datetime, platform
from pathlib import Path
import numpy as np
from scipy.special import i0, i1
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
OUT = HERE / "results"
OUT.mkdir(exist_ok=True)

EPS = 0.35
rng = np.random.default_rng(20260907)

def wrap(a):
    return (a + np.pi) % (2 * np.pi) - np.pi

def w_cos2(dphi):
    return np.cos(dphi / 2) ** 2

def w_band(dphi):
    return (np.abs(wrap(dphi)) <= EPS).astype(float)

# ----------------------------------------------------------------------------- Part A
N_LIST = [1, 2, 3, 5, 10, 20, 26, 46, 100, 1000, 10000]
A95 = 2 * np.arccos(np.sqrt(0.95)); A99 = 2 * np.arccos(np.sqrt(0.99))

def part_a():
    rows = []; checks = []
    for N in N_LIST:
        M = int(min(400000, 4e6 // N))
        dphi = rng.uniform(-np.pi, np.pi, (M, N))
        w = w_cos2(dphi)
        m_best = np.abs(dphi).min(axis=1)
        w_max = w.max(axis=1)
        share_r1 = w_max / w.sum(axis=1)
        # R2: sample one atom per trial with probability ∝ w
        cdf = np.cumsum(w, axis=1); u = rng.uniform(0, 1, M) * cdf[:, -1]
        k = (cdf < u[:, None]).sum(axis=1)
        m_r2 = np.abs(dphi[np.arange(M), k])
        band_hit = (np.abs(dphi) <= EPS).any(axis=1)
        cf = dict(mean_best=np.pi / (N + 1), p95=1 - (1 - A95 / np.pi) ** N,
                  p99=1 - (1 - A99 / np.pi) ** N, band=1 - (1 - EPS / np.pi) ** N)
        mc = dict(mean_best=m_best.mean(), p95=(w_max >= 0.95).mean(), p99=(w_max >= 0.99).mean(),
                  band=band_hit.mean())
        rows.append([N, f"{cf['mean_best']:.4f}", f"{mc['mean_best']:.4f}", f"{cf['p95']:.4f}", f"{mc['p95']:.4f}",
                     f"{cf['p99']:.4f}", f"{mc['p99']:.4f}", f"{share_r1.mean():.5f}", f"{2/N:.5f}",
                     f"{m_r2.mean():.4f}", f"{(m_r2 <= EPS).mean():.4f}", f"{cf['band']:.4f}", f"{mc['band']:.4f}"])
        for key in cf:
            checks.append(abs(cf[key] - mc[key]))
    with open(OUT / "part_a_uniform.csv", "w", newline="") as f:
        wr = csv.writer(f)
        wr.writerow(["N_atoms", "best_mismatch_mean_closed", "best_mismatch_mean_mc", "P_wmax_ge_0.95_closed", "P_wmax_ge_0.95_mc",
                     "P_wmax_ge_0.99_closed", "P_wmax_ge_0.99_mc", "R1_best_share_mc", "R1_best_share_approx_2_over_N",
                     "R2_fired_mismatch_mean_mc", "R2_fired_within_band_mc", "band_hit_closed", "band_hit_mc"])
        wr.writerows(rows)
    return rows, max(checks)

# ----------------------------------------------------------------------------- Part B
R_LIST = [1, 2, 4, 9, 16]
NB_LIST = [1, 2, 5, 10, 100, 1000]

def gamma_of(P, r):
    if r == 1 or P <= 0 or P >= 1:
        return np.nan
    return np.log(P / (1 - P)) / np.log(np.sqrt(r))

def two_lobes(r, n, wfun, M):
    w1 = wfun(rng.uniform(-np.pi, np.pi, (M, n))); w2 = wfun(rng.uniform(-np.pi, np.pi, (M, n)))
    S1, S2 = w1.sum(axis=1), w2.sum(axis=1)
    tot = r * S1 + S2
    live = tot > 0
    P2 = np.mean(r * S1[live] / tot[live])            # R2: probability lobe 1 fires
    P3 = np.mean(r * w1.max(axis=1) > w2.max(axis=1)) if wfun is w_cos2 else np.nan
    miss = 1 - live.mean()
    return P2, P3, miss

def part_b():
    rows = []
    for r in R_LIST:
        for n in NB_LIST:
            M = int(min(200000, 4e6 // n))
            P2, P3, _ = two_lobes(r, n, w_cos2, M)
            P2b, _, miss = two_lobes(r, n, w_band, M)
            rows.append([r, n, f"{r/(1+r):.4f}", f"{P2:.4f}", f"{gamma_of(P2, r):.3f}", f"{P3:.4f}", f"{gamma_of(P3, r):.3f}",
                         f"{P2b:.4f}", f"{gamma_of(P2b, r):.3f}", f"{miss:.4f}"])
    with open(OUT / "part_b_two_lobes.csv", "w", newline="") as f:
        wr = csv.writer(f)
        wr.writerow(["intensity_ratio_r", "n_per_lobe", "born_P1", "R2_P1", "R2_gamma", "R3_P1", "R3_gamma",
                     "R2_band_P1", "R2_band_gamma", "band_miss_rate"])
        wr.writerows(rows)
    return rows

# ----------------------------------------------------------------------------- ordered bulk
KAPPA_LIST = [0.0, 1.0, 4.0]
DP_GRID = np.radians(np.arange(-180, 181, 15.0))

def R_of(k):
    return 0.0 if k == 0 else i1(k) / i0(k)

def ordered_bulk(N=100, M=20000):
    rows = []; worst = 0.0
    for k in KAPPA_LIST:
        Rk = R_of(k)
        for dp in DP_GRID:
            d1 = rng.vonmises(0.0, k, (M, N)); d2 = rng.vonmises(np.pi, k, (M, N))
            w1 = w_cos2(dp - d1); w2 = w_cos2(dp - d2)
            mean_dep = w1.mean()
            cf = 0.5 * (1 + Rk * np.cos(dp))
            worst = max(worst, abs(mean_dep - cf))
            S1, S2 = w1.sum(axis=1), w2.sum(axis=1)
            P2 = np.mean(S1 / (S1 + S2))
            P3 = np.mean(w1.max(axis=1) > w2.max(axis=1))
            rows.append([k, f"{Rk:.4f}", f"{np.degrees(dp):.0f}", f"{mean_dep:.4f}", f"{cf:.4f}", f"{P2:.4f}", f"{cf:.4f}", f"{P3:.4f}"])
    with open(OUT / "ordered_bulk.csv", "w", newline="") as f:
        wr = csv.writer(f)
        wr.writerow(["kappa", "R_kappa", "delta_p_deg", "mean_deposit_half1_mc", "mean_deposit_closed",
                     "two_half_R2_P1_mc", "two_half_R2_P1_closed_largeN", "two_half_R3_P1_mc"])
        wr.writerows(rows)
    return rows, worst

# ----------------------------------------------------------------------------- figures
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

def fig_part_a(rows):
    N = np.array([r[0] for r in rows], float)
    Ns = np.logspace(0, 4, 200)
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.0))
    ax = axes[0]
    for a, lab, c in [(A95, "best atom takes ≥ 95 % (w_max ≥ 0.95)", SLOT[0]), (A99, "best atom takes ≥ 99 %", SLOT[1])]:
        ax.plot(Ns, 1 - (1 - a / np.pi) ** Ns, color=c, lw=1.6, label=lab)
    ax.plot(Ns, 1 - (1 - EPS / np.pi) ** Ns, color=SLOT[2], lw=1.6, label="band form: some atom within ε = 0.35")
    ax.plot(N, [float(r[4]) for r in rows], "o", color=SLOT[0], ms=4, mfc=SURFACE)
    ax.plot(N, [float(r[6]) for r in rows], "o", color=SLOT[1], ms=4, mfc=SURFACE)
    ax.plot(N, [float(r[12]) for r in rows], "o", color=SLOT[2], ms=4, mfc=SURFACE)
    ax.axhline(0.95, color=TEXT2, lw=0.8, alpha=0.5); ax.text(1.05, 0.955, "95 %", color=TEXT2, fontsize=7.5)
    ax.set_xscale("log"); ax.set_ylim(0, 1.02); ax.set_xlabel("N atoms in the aperture")
    ax.set_ylabel("probability"); ax.set_title("Is there an atom close in phase? (lines closed form, dots Monte Carlo)")
    ax.legend(loc="lower right")
    ax = axes[1]
    ax.plot(Ns, np.degrees(np.pi / (Ns + 1)), color=SLOT[1], lw=1.6, label="R3: the closest atom fires, π/(N+1)")
    ax.axhline(np.degrees(np.pi / 2 - 2 / np.pi), color=SLOT[0], lw=1.6, label="R2: sampled by deposit, 53.5° for every N")
    ax.axhline(90, color=TEXT2, lw=0.8, alpha=0.5); ax.text(1.05, 92, "a random atom: 90°", color=TEXT2, fontsize=7.5)
    ax.plot(N, np.degrees([float(r[9]) for r in rows]), "o", color=SLOT[0], ms=4, mfc=SURFACE)
    ax.plot(N, np.degrees([float(r[2]) for r in rows]), "o", color=SLOT[1], ms=4, mfc=SURFACE)
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_ylim(0.01, 200)
    ax.set_xlabel("N atoms in the aperture"); ax.set_ylabel("mean second-spinor mismatch of the fired atom (deg)")
    ax.set_title("How close in phase is the atom that fires?")
    ax.legend(loc="lower left")
    fig.tight_layout(); fig.savefig(OUT / "fig1_part_a.png", dpi=150); plt.close(fig)

def fig_two_lobes(rows):
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.0))
    ax = axes[0]
    r4 = [r for r in rows if r[0] == 4]
    n = np.array([r[1] for r in r4], float)
    ax.axhline(0.8, color=TEXT2, lw=0.8, alpha=0.6); ax.text(1.05, 0.81, "Born  r/(1+r) = 0.80", color=TEXT2, fontsize=7.5)
    ax.plot(n, [float(r[3]) for r in r4], "o-", color=SLOT[0], lw=1.6, ms=4, label="R2: sampled by deposit (primary)")
    ax.plot(n, [float(r[5]) for r in r4], "o-", color=SLOT[1], lw=1.6, ms=4, label="R3: largest deposit")
    ax.plot(n, [float(r[7]) for r in r4], "o-", color=SLOT[2], lw=1.6, ms=4, label="R2 with the band form")
    ax.set_xscale("log"); ax.set_ylim(0.5, 1.02); ax.set_xlabel("n atoms per lobe"); ax.set_ylabel("P(lobe 1 fires)")
    ax.set_title("Two lobes, intensity ratio r = 4")
    ax.legend(loc="lower right")
    ax = axes[1]
    for n_sel, ls in [(1, "-"), (10, "--"), (1000, ":")]:
        sub = [r for r in rows if r[1] == n_sel and r[0] > 1]
        A = np.sqrt([r[0] for r in sub])
        ax.plot(A, [float(r[4]) for r in sub], ls, color=SLOT[0], lw=1.6, marker="o", ms=4, label=f"R2, n = {n_sel}")
        g3 = np.array([float(r[6]) for r in sub])
        ax.plot(A, np.minimum(g3, 6), ls, color=SLOT[1], lw=1.6, marker="o", ms=4, label=f"R3, n = {n_sel}")
    ax.axhline(2, color=TEXT2, lw=0.8, alpha=0.6); ax.text(1.45, 2.06, "Born: γ = 2", color=TEXT2, fontsize=7.5)
    ax.axhline(1.56, color=TEXT2, lw=0.8, alpha=0.6, ls="--"); ax.text(1.45, 1.4, "record's race: 1.56", color=TEXT2, fontsize=7.5)
    ax.set_ylim(0, 6.2); ax.set_xlabel("amplitude ratio √r"); ax.set_ylabel("exponent γ  (P₁/P₂ = (√r)^γ; clipped at 6)")
    ax.set_title("Exponent against amplitude ratio")
    ax.legend(loc="upper right", ncol=2)
    fig.tight_layout(); fig.savefig(OUT / "fig2_two_lobes.png", dpi=150); plt.close(fig)

def fig_ordered(rows):
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.0), sharex=True)
    dps = np.radians(np.linspace(-180, 180, 361))
    for k, c in zip(KAPPA_LIST, SLOT):
        sub = [r for r in rows if r[0] == k]
        d = np.array([float(r[2]) for r in sub])
        axes[0].plot(np.degrees(dps), 0.5 * (1 + R_of(k) * np.cos(dps)), color=c, lw=1.6, label=f"κ = {k:g}  (R = {R_of(k):.2f})")
        axes[0].plot(d, [float(r[3]) for r in sub], "o", color=c, ms=4, mfc=SURFACE)
        axes[1].plot(np.degrees(dps), 0.5 * (1 + R_of(k) * np.cos(dps)), color=c, lw=1.6, label=f"R2, κ = {k:g}")
        axes[1].plot(d, [float(r[7]) for r in sub], "s--", color=c, lw=1.0, ms=4, mfc=SURFACE, label=f"R3, κ = {k:g}")
    axes[0].set_ylabel("mean deposit per atom, half with μ = 0"); axes[0].set_xlabel("incident internal offset δ_p (deg)")
    axes[0].set_title("Ordered bulk: deposit against the particle's own phase"); axes[0].legend(loc="lower center")
    axes[1].set_ylabel("P(half with μ = 0 fires)"); axes[1].set_xlabel("incident internal offset δ_p (deg)")
    axes[1].set_title("Two-half detector (μ = 0 beside μ = π), 100 atoms each"); axes[1].legend(loc="lower center", ncol=2)
    for ax in axes:
        ax.set_xticks([-180, -90, 0, 90, 180]); ax.set_ylim(-0.02, 1.02)
    fig.tight_layout(); fig.savefig(OUT / "fig3_ordered_bulk.png", dpi=150); plt.close(fig)

# ----------------------------------------------------------------------------- main
def main():
    lines = ["Verification (PREREGISTRATION.md §2): Monte Carlo against closed forms", ""]
    rows_a, worst_a = part_a()
    lines.append(f"Part A: worst |MC − closed| over mean best mismatch, P(w_max ≥ 0.95), P(w_max ≥ 0.99), band hit, all N: {worst_a:.2e}"
                 f"   (tolerance 0.01)   -> {'PASS' if worst_a < 0.01 else 'FAIL'}")
    rows_b = part_b()
    r4 = {r[1]: r for r in rows_b if r[0] == 4}
    d_born = abs(float(r4[1000][3]) - 0.8)
    lines.append(f"Part B: R2 at r = 4, n = 1000 against Born 0.800: P₁ = {float(r4[1000][3]):.4f}, |diff| = {d_born:.2e}   -> {'PASS' if d_born < 0.01 else 'FAIL'}")
    lines.append(f"Part B: pre-registered r = 4 table recomputed: " +
                 ", ".join(f"n={n}: R2 {float(r4[n][3]):.3f}/γ {float(r4[n][4]):.2f}, R3 {float(r4[n][5]):.3f}" for n in NB_LIST))
    rows_o, worst_o = ordered_bulk()
    lines.append(f"Ordered bulk: worst |mean deposit MC − ½(1 + R_κ cos δ_p)| = {worst_o:.2e}   -> {'PASS' if worst_o < 0.01 else 'FAIL'}")
    ok = worst_a < 0.01 and d_born < 0.01 and worst_o < 0.01
    (OUT / "verification.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    fig_part_a(rows_a); fig_two_lobes(rows_b); fig_ordered(rows_o)
    print()
    print("Part A (uniform wave, random bulk):")
    print("   N     best mismatch  P(w_max≥.95)  P(w_max≥.99)  R1 best share  R2 fired mismatch  band hit")
    for r in rows_a:
        print(f"  {r[0]:5d}   {np.degrees(float(r[2])):6.1f}°       {float(r[4]):.3f}         {float(r[6]):.3f}        {float(r[7]):.5f}       {np.degrees(float(r[9])):5.1f}°         {float(r[12]):.3f}")
    print("Part B (two lobes), r = 4, Born 0.800:")
    for n in NB_LIST:
        r = r4[n]; print(f"  n = {n:5d}: R2 P₁ = {float(r[3]):.3f} (γ {float(r[4]):.2f})   R3 P₁ = {float(r[5]):.3f} (γ {r[6]})   band R2 P₁ = {float(r[7]):.3f}, miss {float(r[9]):.3f}")
    info = [f"run: {datetime.datetime.now().isoformat(timespec='seconds')}",
            f"python {platform.python_version()}, numpy {np.__version__}", f"verification (phase_race.py): {'PASS' if ok else 'FAIL'}"]
    (OUT / "run_info.txt").write_text("\n".join(info) + "\n")
    print("\n".join(info))
    return 0 if ok else 1

if __name__ == "__main__":
    sys.exit(main())
