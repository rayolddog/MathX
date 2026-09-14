#!/usr/bin/env python3
"""
frozen_multipliers.py: a frozen per-atom multiplier on the hazard, light-tailed against the
infinite-mean e^Gumbel of addendum 1 (PREREGISTRATION.md addendum 3, D19 to D21). First to fire
with fresh exponential clocks; the pixel that fires is sampled exactly from the rate sums.
"""
import csv, time
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent; OUT = HERE / "results"
SEED = 20260918
rng = np.random.default_rng(SEED)
R, M = 4.0, 1000
NS = [10, 100, 1000]
FORMS = ["none", "uniform[0.5,1.5]", "lognormal(0.5)", "exponential(1)", "exp(Gumbel)"]
SURFACE, TEXT, TEXT2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e6e5e1"
SLOT = ["#8c8c8c", "#2a78d6", "#1baf7a", "#eda100", "#eb6834"]
plt.rcParams.update({"font.family": "sans-serif", "font.size": 9, "axes.titlesize": 10, "axes.titleweight": "bold",
    "axes.edgecolor": GRID, "axes.labelcolor": TEXT2, "xtick.color": TEXT2, "ytick.color": TEXT2, "axes.titlecolor": TEXT,
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE, "axes.grid": True, "grid.color": GRID,
    "grid.linewidth": 0.8, "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False, "legend.fontsize": 8, "axes.axisbelow": True})

def multipliers(form, shape):
    if form == "none": return np.ones(shape)
    if form == "uniform[0.5,1.5]": return rng.uniform(0.5, 1.5, shape)
    if form == "lognormal(0.5)": return np.exp(0.5 * rng.standard_normal(shape))
    if form == "exponential(1)": return rng.exponential(1.0, shape)
    if form == "exp(Gumbel)": return np.exp(rng.gumbel(0, 1, shape))
    raise ValueError(form)

def run(form, N, phase, K):
    mA, mB = multipliers(form, (K, N)), multipliers(form, (K, N))
    winA = np.zeros(K)
    for _ in range(M):
        if phase:
            wA = np.cos(rng.uniform(0, 2 * np.pi, (K, N)) / 2) ** 2; wB = np.cos(rng.uniform(0, 2 * np.pi, (K, N)) / 2) ** 2
            LA, LB = R * (wA * mA).sum(1), (wB * mB).sum(1)
        else:
            LA, LB = R * mA.sum(1), mB.sum(1)
        winA += rng.uniform(0, 1, K) < LA / (LA + LB)
    s = winA / M
    p = s.mean(); v = s.var(ddof=1)
    return p, np.sqrt(v / K), v / (p * (1 - p) / M), np.mean((s < 0.05) | (s > 0.95))

def main():
    t0 = time.time(); lines = ["Verification (PREREGISTRATION.md addendum 3)", ""]; rows = []
    for form in FORMS:
        for phase in (False, True):
            for N in NS:
                K = 200 if N == 1000 else 300
                p, se, F, step = run(form, N, phase, K)
                rows.append(dict(form=form, phase_weight=phase, N=N, K=K, M=M, mean=p, se_mean=se, fano=F, step_fraction=step))
                print(f"{form:18s} phase={phase!s:5} N={N:5d}: mean {p:.3f} (s.e. {se:.4f}), F {F:8.1f}, step {step:.3f}")
                if not phase:
                    if form == "none" and N == 10:
                        ok = 0.8 <= F <= 1.2; lines.append(f"no multiplier, N = 10: F {F:.2f} in [0.8, 1.2] -> {'PASS' if ok else 'FAIL'}")
                    if form == "uniform[0.5,1.5]" and N == 1000:
                        ok = 0.7 <= F <= 1.3; lines.append(f"uniform, N = 1000: F {F:.2f} in [0.7, 1.3] -> {'PASS' if ok else 'FAIL'}")
                    if form == "exp(Gumbel)" and N == 1000:
                        ok = F > 50; lines.append(f"exp(Gumbel), N = 1000: F {F:.1f} > 50 -> {'PASS' if ok else 'FAIL'}")
                    if form == "lognormal(0.5)" and N == 100:
                        ok = abs(p - 0.8) < 0.02; lines.append(f"lognormal, N = 100: mean {p:.4f} within 0.02 of 0.800 -> {'PASS' if ok else 'FAIL'}")
    fmt = lambda v: (f"{v:.5g}" if isinstance(v, float) else v)
    with open(OUT / "frozen_multipliers.csv", "w", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=list(rows[0].keys())); wr.writeheader()
        for row in rows: wr.writerow({k: fmt(v) for k, v in row.items()})
    fig, ax = plt.subplots(figsize=(6.2, 3.9))
    for k, form in enumerate(FORMS):
        sel = [x for x in rows if x["form"] == form and not x["phase_weight"]]
        ax.plot([x["N"] for x in sel], [x["fano"] for x in sel], "o-", color=SLOT[k], label=f"frozen multiplier: {form}")
    ax.axhline(1, color=TEXT, lw=0.9); ax.text(11, 1.15, "binomial, F = 1", fontsize=8)
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("atoms per pixel N"); ax.set_ylabel("Fano factor of one detector's split")
    ax.set_title("A frozen multiplier self-averages unless its mean is infinite"); ax.legend(fontsize=7.5)
    fig.tight_layout(); fig.savefig(OUT / "fig4_frozen_multipliers.png", dpi=150); plt.close(fig)
    text = "\n".join(lines) + "\n"; (OUT / "verification_multipliers.txt").write_text(text); print(text)
    ok = all("FAIL" not in l for l in lines)
    with open(OUT / "run_info.txt", "a") as f:
        f.write(f"frozen_multipliers.py, run {time.strftime('%Y-%m-%d %H:%M:%S')}, seed {SEED}, {time.time() - t0:.1f} s, numpy {np.__version__}\nverification (multipliers): {'PASS' if ok else 'FAIL'}\n")
    print(f"{time.time() - t0:.1f} s")

if __name__ == "__main__":
    main()
