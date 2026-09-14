#!/usr/bin/env python3
"""
staggered_arrival.py: the race under staggered arrival (PREREGISTRATION.md addendum 2, D13 to D16).
Reading S: the survival hazard, rate λ I_P w_i while the packet is at pixel P. Reading C: the
conditional-state hazard, rate λ w_i a_P²(t) with a_P² consumed by the port's own no-fire absorption,
integrating to the port's share. First fire anywhere stops every clock. Exact sampling.
"""
import csv, time
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent; OUT = HERE / "results"; OUT.mkdir(exist_ok=True)
SEED = 20260916
rng = np.random.default_rng(SEED)
N, T, TRIALS = 10, 1.0, 200000
DELTAS = [0.0, 0.1, 0.5, 1.0, 2.0]
KAPPAS = [1.0, 10.0, 100.0]
SURFACE, TEXT, TEXT2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e6e5e1"
SLOT = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#7a5cc7", "#8c8c8c"]
plt.rcParams.update({"font.family": "sans-serif", "font.size": 9, "axes.titlesize": 10, "axes.titleweight": "bold",
    "axes.edgecolor": GRID, "axes.labelcolor": TEXT2, "xtick.color": TEXT2, "ytick.color": TEXT2, "axes.titlecolor": TEXT,
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE, "axes.grid": True, "grid.color": GRID,
    "grid.linewidth": 0.8, "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False, "legend.fontsize": 8, "axes.axisbelow": True})

def weights(phase):
    if phase:
        return np.cos(rng.uniform(0, 2 * np.pi, (TRIALS, N)) / 2).sum(1) * 0 + (np.cos(rng.uniform(0, 2 * np.pi, (TRIALS, N)) / 2) ** 2).sum(1)
    return np.full(TRIALS, float(N))

def fire_time_S(rate, window):
    """First-passage time of a pixel under the survival hazard within its window; inf if none."""
    t = rng.exponential(1.0 / rate)
    return np.where(t < window, t, np.inf)

def fire_time_C(lam_W, share, window):
    """Reading C: P(fire) = share (1 - exp(-lam_W window)); time is the truncated exponential."""
    u = rng.uniform(0, 1, lam_W.size)
    p = share * (1 - np.exp(-lam_W * window))
    fires = u < p
    t = np.full(u.size, np.inf)
    t[fires] = -np.log(1 - u[fires] / share) / lam_W[fires]
    return t

def run(reading, r, kappa, delta, phase):
    lam = kappa / (N * T)
    WA, WB = weights(phase), weights(phase)
    sA, sB = r / (r + 1), 1 / (r + 1)
    if reading == "S":
        tA = fire_time_S(lam * r * WA, T)
        tB = delta + fire_time_S(lam * WB, T)
    else:
        tA = fire_time_C(lam * WA, sA, T)
        tB = delta + fire_time_C(lam * WB, sB, T)
    a = tA < tB; b = tB < tA
    PA, PB = a.mean(), b.mean()
    return PA, PB, 1 - PA - PB, PA / (PA + PB)

def closed_S(r, kappa, delta):
    lA, lB = kappa * r / T, kappa / T
    if delta <= 0:
        f = 1 - np.exp(-(lA + lB) * T); PA, PB = lA / (lA + lB) * f, lB / (lA + lB) * f
    elif delta >= T:
        PA = 1 - np.exp(-lA * T); PB = np.exp(-lA * T) * (1 - np.exp(-lB * T))
    else:
        e1 = np.exp(-lA * delta); f = 1 - np.exp(-(lA + lB) * (T - delta)); e2 = np.exp(-(lA + lB) * (T - delta))
        PA = (1 - e1) + e1 * lA / (lA + lB) * f
        PB = e1 * lB / (lA + lB) * f + e1 * e2 * (1 - np.exp(-lB * delta))
    return PA, PB, PA / (PA + PB)

def main():
    t0 = time.time(); lines = ["Verification (PREREGISTRATION.md addendum 2)", ""]; rows = []
    for r in (4, 2):
        for phase in (False, True):
            for kappa in KAPPAS:
                for delta in DELTAS:
                    for reading in ("S", "C"):
                        PA, PB, P0, share = run(reading, r, kappa, delta, phase)
                        se = np.sqrt(share * (1 - share) / max((PA + PB) * TRIALS, 1))
                        cs = closed_S(r, kappa, delta)[2] if (reading == "S" and not phase) else (r / (r + 1) if reading == "C" else np.nan)
                        rows.append(dict(reading=reading, r=r, phase_weight=phase, kappa=kappa, delta_over_T=delta, P_A=PA, P_B=PB, P_none=P0, share=share, se=se, closed_form=cs))
                        if r == 4 and not phase:
                            if reading == "S":
                                ok = abs(share - cs) < max(3 * se, 1e-4)
                                lines.append(f"S, no phase, kappa = {kappa:5.0f}, Δ/T = {delta:3.1f}: share {share:.4f} against closed form {cs:.4f} (s.e. {se:.4f}) -> {'PASS' if ok else 'FAIL'}")
                                if kappa == 10 and delta == 1.0:
                                    lines.append(f"S, no phase, kappa = 10, Δ = T: P(A) = {PA:.4f} >= 0.999 -> {'PASS' if PA >= 0.999 else 'FAIL'}")
                            else:
                                if kappa == 10:
                                    ok = abs(share - 0.8) < 3 * se and (PA + PB) >= 0.999
                                    lines.append(f"C, no phase, kappa = 10, Δ/T = {delta:3.1f}: share {share:.4f} against 0.800 (s.e. {se:.4f}), P(fired) = {PA+PB:.4f} >= 0.999 -> {'PASS' if ok else 'FAIL'}")
                                if kappa == 1 and delta == 0.0:
                                    pf = PA + PB; sef = np.sqrt(pf * (1 - pf) / TRIALS); ok = abs(pf - (1 - np.exp(-1))) < 3 * sef
                                    lines.append(f"C, no phase, kappa = 1: P(fired) = {pf:.4f} against 1 - e^-1 = {1-np.exp(-1):.4f} (s.e. {sef:.4f}) -> {'PASS' if ok else 'FAIL'}")
    fmt = lambda v: (f"{v:.5g}" if isinstance(v, float) else v)
    with open(OUT / "staggered.csv", "w", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=list(rows[0].keys())); wr.writeheader()
        for row in rows: wr.writerow({k: fmt(v) for k, v in row.items()})
    fig, axes = plt.subplots(1, 2, figsize=(9.8, 3.9))
    for ax, phase, title in ((axes[0], False, "no phase weight"), (axes[1], True, "phase weight cos²(Δφ/2), fresh")):
        for k, kappa in enumerate(KAPPAS):
            sel = [x for x in rows if x["r"] == 4 and x["phase_weight"] == phase and x["kappa"] == kappa]
            ax.plot([x["delta_over_T"] for x in sel if x["reading"] == "S"], [x["share"] for x in sel if x["reading"] == "S"], "o-", color=SLOT[k], label=f"survival hazard, κ = {kappa:g}")
            ax.plot([x["delta_over_T"] for x in sel if x["reading"] == "C"], [x["share"] for x in sel if x["reading"] == "C"], "s--", color=SLOT[k], label=f"conditional-state hazard, κ = {kappa:g}")
        ax.axhline(0.8, color=TEXT, lw=0.9); ax.text(1.3, 0.81, "Born 0.800", fontsize=8)
        ax.set_xlabel("arrival stagger Δ/T (packet lengths)"); ax.set_ylabel("share of records at the nearer, brighter pixel (r = 4)"); ax.set_ylim(0.75, 1.02)
        ax.set_title(f"Staggered arrival, {title}"); ax.legend(fontsize=7, loc="center right")
    fig.tight_layout(); fig.savefig(OUT / "fig3_staggered.png", dpi=150); plt.close(fig)
    text = "\n".join(lines) + "\n"; (OUT / "verification_staggered.txt").write_text(text); print(text)
    print("shares, r = 4:")
    for phase in (False, True):
        for kappa in KAPPAS:
            for reading in ("S", "C"):
                sel = [x for x in rows if x["r"] == 4 and x["phase_weight"] == phase and x["kappa"] == kappa and x["reading"] == reading]
                print(f"  {reading} phase={phase!s:5} kappa={kappa:5.0f}: " + "  ".join(f"Δ/T={x['delta_over_T']:.1f}: {x['share']:.3f} (fired {x['P_A']+x['P_B']:.3f})" for x in sel))
    ok = all("FAIL" not in l for l in lines)
    with open(OUT / "run_info.txt", "a") as f:
        f.write(f"staggered_arrival.py, run {time.strftime('%Y-%m-%d %H:%M:%S')}, seed {SEED}, {time.time() - t0:.1f} s, numpy {np.__version__}\nverification (staggered): {'PASS' if ok else 'FAIL'}\n")
    print(f"{time.time() - t0:.1f} s")

if __name__ == "__main__":
    main()
