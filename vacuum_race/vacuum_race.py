#!/usr/bin/env python3
"""
vacuum_race.py: "highest interference term wins" with a per-atom vacuum term
(PREREGISTRATION.md of this folder). Two pixels at intensity ratio r, N atoms per pixel, the
second study's phase weight cos²(Δφ/2), and six forms of the vacuum term (D4). The atom with the
largest deposit among all 2N takes the quantum. Recorded: the pixel share, the exponent, the
within-pixel best-atom fraction, and the intensity-scaling test.
"""
import csv, time
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent; OUT = HERE / "results"; OUT.mkdir(exist_ok=True)
SEED = 20260914
rng = np.random.default_rng(SEED)
NS = [1, 10, 100, 1000]
FORMS = ["V0 none", "V1 Gumbel on deposit", "V2 Gumbel on amplitude", "V3 mult Gaussian", "V4 add Gaussian", "V5 add uniform"]
CHUNK = 4000
SURFACE, TEXT, TEXT2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e6e5e1"
SLOT = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#7a5cc7", "#8c8c8c"]
plt.rcParams.update({"font.family": "sans-serif", "font.size": 9, "axes.titlesize": 10, "axes.titleweight": "bold",
    "axes.edgecolor": GRID, "axes.labelcolor": TEXT2, "xtick.color": TEXT2, "ytick.color": TEXT2, "axes.titlecolor": TEXT,
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE, "axes.grid": True, "grid.color": GRID,
    "grid.linewidth": 0.8, "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False, "legend.fontsize": 8, "axes.axisbelow": True})

def trials_for(N):
    return 100000 if N <= 100 else 20000

def vacuum(form, d, beta):
    """D4: the per-atom vacuum term applied to the deposit d (shape trials × N)."""
    if form == "V0 none":
        return d
    if form == "V1 Gumbel on deposit":
        return d * np.exp(rng.gumbel(0.0, 1.0, d.shape) / beta)
    if form == "V2 Gumbel on amplitude":
        return np.sqrt(d) * np.exp(rng.gumbel(0.0, 1.0, d.shape))
    if form == "V3 mult Gaussian":
        return d * (1.0 + 0.5 * rng.standard_normal(d.shape))
    if form == "V4 add Gaussian":
        return d + 2.0 * rng.standard_normal(d.shape)
    if form == "V5 add uniform":
        return d + rng.uniform(0.0, 4.0, d.shape)
    raise ValueError(form)

def run(form, IA, IB, N, phase=True, beta=1.0, T=None):
    """D1 to D6. Returns P(A), s.e., fraction of trials in which the winner is its pixel's best-matched atom."""
    T = trials_for(N) if T is None else T
    winA = 0; best = 0
    for start in range(0, T, CHUNK):
        t = min(CHUNK, T - start)
        wA = np.cos(rng.uniform(0, 2 * np.pi, (t, N)) / 2) ** 2 if phase else np.ones((t, N))
        wB = np.cos(rng.uniform(0, 2 * np.pi, (t, N)) / 2) ** 2 if phase else np.ones((t, N))
        EA, EB = vacuum(form, IA * wA, beta), vacuum(form, IB * wB, beta)
        mA, mB = EA.max(1), EB.max(1)
        a = mA > mB
        winA += a.sum()
        best += np.where(a, EA.argmax(1) == wA.argmax(1), EB.argmax(1) == wB.argmax(1)).sum()
    p = winA / T
    return p, np.sqrt(max(p * (1 - p), 1e-12) / T), best / T

def gamma(p, r):
    return np.log(p / (1 - p)) / np.log(np.sqrt(r)) if 0 < p < 1 else np.inf

def main():
    t0 = time.time(); lines = ["Verification (PREREGISTRATION.md §4)", ""]
    R2_ROW = {1: 0.704, 10: 0.794, 100: 0.800, 1000: 0.800}; SE_R2 = 0.0014
    share, best = [], []
    # Table 1: every form, r = 4 and 2, with and without the phase weight
    for form in FORMS:
        for r in (4, 2):
            for phase in (True, False):
                if r == 2 and not phase and form not in ("V1 Gumbel on deposit",):
                    continue
                for N in NS:
                    p, se, b = run(form, r, 1.0, N, phase)
                    share.append(dict(form=form, r=r, phase_weight=phase, N=N, P_A=p, se=se, born=r / (r + 1), gamma=gamma(p, r)))
                    if phase:
                        best.append(dict(form=form, r=r, N=N, best_atom_wins=b, two_over_N=min(2 / N, 1.0)))
                    if form == "V1 Gumbel on deposit" and not phase:
                        ok = abs(p - r / (r + 1)) < 3 * se
                        lines.append(f"V1, no phase weight, r = {r}, N = {N:4d}: P(A) = {p:.4f} against r/(r+1) = {r/(r+1):.4f} (s.e. {se:.4f}) -> {'PASS' if ok else 'FAIL'}")
                    if form == "V1 Gumbel on deposit" and phase and r == 4:
                        comb = np.sqrt(se ** 2 + SE_R2 ** 2); ok = abs(p - R2_ROW[N]) < 3 * comb
                        lines.append(f"V1, phase weight, r = 4, N = {N:4d}: P(A) = {p:.4f} against R2 row {R2_ROW[N]:.3f} (combined s.e. {comb:.4f}) -> {'PASS' if ok else 'FAIL'}")
                    if form == "V0 none" and phase and r == 4 and N == 10:
                        ok = p >= 0.999
                        lines.append(f"V0, phase weight, r = 4, N = 10: P(A) = {p:.4f} against >= 0.999 -> {'PASS' if ok else 'FAIL'}")
                    if form == "V2 Gumbel on amplitude" and not phase and r == 4 and N == 1:
                        pred = np.sqrt(r) / (np.sqrt(r) + 1); ok = abs(p - pred) < 3 * se
                        lines.append(f"V2, no phase weight, r = 4, N = 1: P(A) = {p:.4f} against sqrt(r)/(sqrt(r)+1) = {pred:.4f} (s.e. {se:.4f}) -> {'PASS' if ok else 'FAIL'}")
    # Table 2: beta
    beta_rows = []
    for beta in (0.5, 1.0, 2.0):
        for phase in (True, False):
            p, se, _ = run("V1 Gumbel on deposit", 4, 1.0, 100, phase, beta)
            pred = 4 ** beta / (4 ** beta + 1)
            beta_rows.append(dict(beta=beta, phase_weight=phase, N=100, r=4, P_A=p, se=se, predicted=pred, gamma=gamma(p, 4), two_beta=2 * beta))
            if beta == 2.0 and not phase:
                ok = abs(p - pred) < 3 * se
                lines.append(f"V1, beta = 2, no phase weight, r = 4, N = 100: P(A) = {p:.4f} against r^2/(r^2+1) = {pred:.4f} (s.e. {se:.4f}) -> {'PASS' if ok else 'FAIL'}")
    # Table 3: scaling
    scal_rows = []
    for form in ("V1 Gumbel on deposit", "V3 mult Gaussian", "V4 add Gaussian", "V5 add uniform"):
        for phase in (True, False):
            for N in (1, 10):
                p1, se1, _ = run(form, 4, 1.0, N, phase); p10, se10, _ = run(form, 40, 10.0, N, phase)
                scal_rows.append(dict(form=form, phase_weight=phase, N=N, P_A_at_4_1=p1, P_A_at_40_10=p10, change=p10 - p1, se_change=np.sqrt(se1 ** 2 + se10 ** 2)))
    fmt = lambda v: (f"{v:.5g}" if isinstance(v, float) else v)
    for name, rows in (("share.csv", share), ("best_atom.csv", best), ("beta.csv", beta_rows), ("scaling.csv", scal_rows)):
        with open(OUT / name, "w", newline="") as f:
            wr = csv.DictWriter(f, fieldnames=list(rows[0].keys())); wr.writeheader()
            for row in rows: wr.writerow({k: fmt(v) for k, v in row.items()})
    # figure
    fig, axes = plt.subplots(1, 2, figsize=(9.8, 3.9))
    ax = axes[0]
    for k, form in enumerate(FORMS):
        sel = [s for s in share if s["form"] == form and s["r"] == 4 and s["phase_weight"]]
        ax.plot([s["N"] for s in sel], [s["P_A"] for s in sel], "o-", color=SLOT[k], label=form)
    ax.axhline(0.8, color=TEXT, lw=0.9); ax.text(1.1, 0.81, "Born 0.800", fontsize=8)
    ax.axhline(2 / 3, color=TEXT, lw=0.9, ls=":"); ax.text(1.1, 0.675, "√r/(√r + 1) = 0.667", fontsize=8)
    ax.set_xscale("log"); ax.set_xlabel("atoms per pixel N"); ax.set_ylabel("P(brighter pixel), r = 4"); ax.set_ylim(0.6, 1.02)
    ax.set_title("Highest deposit wins, six forms of the vacuum term"); ax.legend(fontsize=7)
    ax = axes[1]
    for k, form in enumerate(FORMS):
        sel = [b for b in best if b["form"] == form and b["r"] == 4]
        ax.plot([b["N"] for b in sel], [max(b["best_atom_wins"], 1e-4) for b in sel], "o-", color=SLOT[k], label=form)
    nn = np.array(NS); ax.plot(nn, np.minimum(2 / nn, 1), "--", color=TEXT2, label="2/N")
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("atoms per pixel N"); ax.set_ylabel("P(winner is its pixel's best-matched atom)")
    ax.set_title("Within the pixel: most likely, never certain, under V1"); ax.legend(fontsize=7)
    fig.tight_layout(); fig.savefig(OUT / "fig1_vacuum_race.png", dpi=150); plt.close(fig)
    text = "\n".join(lines) + "\n"; (OUT / "verification.txt").write_text(text); print(text)
    print("share (r = 4, phase weight):")
    for form in FORMS:
        sel = [s for s in share if s["form"] == form and s["r"] == 4 and s["phase_weight"]]
        print(f"  {form:24s} " + "  ".join(f"N={s['N']:4d}: {s['P_A']:.3f} (γ {s['gamma']:.2f})" for s in sel))
    print("scaling:")
    for s in scal_rows: print(f"  {s['form']:24s} phase={s['phase_weight']!s:5} N={s['N']:2d}: {s['P_A_at_4_1']:.3f} -> {s['P_A_at_40_10']:.3f}")
    ok = all("FAIL" not in l for l in lines)
    with open(OUT / "run_info.txt", "a") as f:
        f.write(f"vacuum_race.py, run {time.strftime('%Y-%m-%d %H:%M:%S')}, seed {SEED}, {time.time() - t0:.1f} s, numpy {np.__version__}\nverification: {'PASS' if ok else 'FAIL'}\n")
    print(f"{time.time() - t0:.1f} s")

if __name__ == "__main__":
    main()
