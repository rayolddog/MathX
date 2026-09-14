#!/usr/bin/env python3
"""
frozen_fresh.py: the vacuum term frozen per detector, refreshed every τ events, or fresh per event
(PREREGISTRATION.md addendum, D8 to D12). V1 Gumbel on the deposit at β = 1, r = 4. Recorded per
realization: the share of events at pixel A; across realizations: mean, standard deviation, Fano
factor against the binomial split, and the step fraction.
"""
import csv, time
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent; OUT = HERE / "results"; OUT.mkdir(exist_ok=True)
SEED = 20260915
rng = np.random.default_rng(SEED)
R, M = 4.0, 1000
TAUS = [1, 10, 100, None]  # None = frozen
SURFACE, TEXT, TEXT2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e6e5e1"
SLOT = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#7a5cc7", "#8c8c8c"]
plt.rcParams.update({"font.family": "sans-serif", "font.size": 9, "axes.titlesize": 10, "axes.titleweight": "bold",
    "axes.edgecolor": GRID, "axes.labelcolor": TEXT2, "xtick.color": TEXT2, "ytick.color": TEXT2, "axes.titlecolor": TEXT,
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE, "axes.grid": True, "grid.color": GRID,
    "grid.linewidth": 0.8, "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False, "legend.fontsize": 8, "axes.axisbelow": True})

def phase_weights(shape):
    return np.cos(rng.uniform(0, 2 * np.pi, shape) / 2) ** 2

def realizations(N, tau, phase, K):
    """D8 to D10. Returns the K per-realization shares s_k. Vectorized over realizations."""
    winA = np.zeros(K)
    if phase == "frozen":
        wA0, wB0 = phase_weights((K, N)), phase_weights((K, N))
    GA = GB = None
    for m in range(M):
        if tau is None:
            if GA is None:
                GA, GB = rng.gumbel(0, 1, (K, N)), rng.gumbel(0, 1, (K, N))
        elif m % tau == 0:
            GA, GB = rng.gumbel(0, 1, (K, N)), rng.gumbel(0, 1, (K, N))
        if phase == "none":
            wA = wB = 1.0
        elif phase == "fresh":
            wA, wB = phase_weights((K, N)), phase_weights((K, N))
        else:
            wA, wB = wA0, wB0
        winA += (R * wA * np.exp(GA)).max(1) > (wB * np.exp(GB)).max(1)
    return winA / M

def stats(s):
    p = s.mean(); v = s.var(ddof=1)
    fano = v / (p * (1 - p) / M) if 0 < p < 1 else np.nan
    return p, np.sqrt(v), fano, np.mean((s < 0.05) | (s > 0.95)), np.sqrt(v / s.size)

def main():
    t0 = time.time(); lines = ["Verification (PREREGISTRATION.md addendum)", ""]
    rows, shares_rows = [], []
    cells = [(10, tau, ph, 300) for ph in ("none", "fresh") for tau in TAUS]
    cells += [(10, 1, "frozen", 300), (100, 1, "frozen", 200), (100, None, "fresh", 200), (1000, 1, "frozen", 100)]
    for N, tau, phase, K in cells:
        s = realizations(N, tau, phase, K)
        p, sd, fano, step, se = stats(s)
        label = "frozen" if tau is None else str(tau)
        rows.append(dict(N=N, tau=label, phase_weight=phase, K=K, M=M, mean=p, se_mean=se, sd=sd, fano=fano, step_fraction=step))
        for k, v in enumerate(s):
            shares_rows.append(dict(N=N, tau=label, phase_weight=phase, realization=k, share=v))
        print(f"N = {N:4d}, tau = {label:>6}, phase {phase:6s}: mean {p:.3f} (s.e. {se:.4f}), sd {sd:.4f}, Fano {fano:8.1f}, step {step:.3f}")
        if N == 10 and phase == "none":
            if tau == 1:
                ok = abs(p - 0.8) < 3 * se and 0.7 <= fano <= 1.3
                lines.append(f"fresh, no phase weight: mean {p:.4f} against 0.800 (s.e. {se:.4f}), Fano {fano:.2f} in [0.7, 1.3] -> {'PASS' if ok else 'FAIL'}")
            if tau is None:
                ok = step == 1.0 and abs(fano - M) < 0.01 * M and abs(p - 0.8) < 3 * 0.4 / np.sqrt(K)
                lines.append(f"frozen, no phase weight: step {step:.3f} against 1.000, Fano {fano:.1f} against {M} within 1 %, mean {p:.4f} against 0.800 (s.e. {0.4/np.sqrt(K):.4f}) -> {'PASS' if ok else 'FAIL'}")
            if tau == 10:
                ok = 7 <= fano <= 13
                lines.append(f"tau = 10, no phase weight: Fano {fano:.2f} in [7, 13] -> {'PASS' if ok else 'FAIL'}")
    fmt = lambda v: (f"{v:.5g}" if isinstance(v, float) else v)
    for name, rr in (("frozen_fresh.csv", rows), ("frozen_fresh_shares.csv", shares_rows)):
        with open(OUT / name, "w", newline="") as f:
            wr = csv.DictWriter(f, fieldnames=list(rr[0].keys())); wr.writeheader()
            for row in rr: wr.writerow({k: fmt(v) for k, v in row.items()})
    # figure: histograms of per-realization share, N = 10, fresh phase weight, four tau; and Fano against tau
    fig, axes = plt.subplots(1, 2, figsize=(9.8, 3.9))
    ax = axes[0]; bins = np.linspace(0, 1, 41)
    for k, tau in enumerate(TAUS):
        label = "frozen" if tau is None else str(tau)
        s = np.array([r["share"] for r in shares_rows if r["N"] == 10 and r["tau"] == label and r["phase_weight"] == "fresh"])
        ax.hist(s, bins=bins, histtype="step", lw=1.4, color=SLOT[k], label=f"τ = {label}")
    ax.axvline(0.8, color=TEXT, lw=0.9); ax.text(0.81, ax.get_ylim()[1] * 0.9, "Born 0.800", fontsize=8)
    ax.set_xlabel("share of events at the brighter pixel, one detector (M = 1000)"); ax.set_ylabel("detectors (of 300)")
    ax.set_title("What one detector shows, by the vacuum's correlation time"); ax.legend()
    ax = axes[1]
    for ph, c, mk in (("none", SLOT[0], "o"), ("fresh", SLOT[1], "s")):
        sel = [r for r in rows if r["N"] == 10 and r["phase_weight"] == ph]
        xs = [M if r["tau"] == "frozen" else int(r["tau"]) for r in sel]
        ax.plot(xs, [r["fano"] for r in sel], mk + "-", color=c, label=f"phase weight {ph}")
    tt = np.array([1, 10, 100, 1000]); ax.plot(tt, tt, "--", color=TEXT2, label="F = τ"); ax.plot(tt, tt / 2, ":", color=TEXT2, label="F = τ/2")
    ax.axhline(1, color=TEXT, lw=0.9); ax.text(1.1, 1.3, "binomial split, F = 1 (experiment)", fontsize=8)
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("vacuum correlation time τ (events); frozen plotted at M = 1000"); ax.set_ylabel("Fano factor of the split")
    ax.set_title("The Fano factor reads the correlation time"); ax.legend(loc="upper left")
    fig.tight_layout(); fig.savefig(OUT / "fig2_frozen_fresh.png", dpi=150); plt.close(fig)
    text = "\n".join(lines) + "\n"; (OUT / "verification_frozen.txt").write_text(text); print(text)
    ok = all("FAIL" not in l for l in lines)
    with open(OUT / "run_info.txt", "a") as f:
        f.write(f"frozen_fresh.py, run {time.strftime('%Y-%m-%d %H:%M:%S')}, seed {SEED}, {time.time() - t0:.1f} s, numpy {np.__version__}\nverification (frozen): {'PASS' if ok else 'FAIL'}\n")
    print(f"{time.time() - t0:.1f} s")

if __name__ == "__main__":
    main()
