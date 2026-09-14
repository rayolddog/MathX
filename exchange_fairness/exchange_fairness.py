#!/usr/bin/env python3
"""
exchange_fairness.py: fairness in the conserving exchange class (PREREGISTRATION.md D1 to D6).
Model A independent CIR noise (does not conserve); B conserving geometric-mean kernel; C conserving
multiplicative kernel with rate-linear and threshold commits; D fixed parcels with a reflecting zero.
"""
import csv, time
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent; OUT = HERE / "results"; OUT.mkdir(exist_ok=True)
SEED = 20260919
rng = np.random.default_rng(SEED)
TRIALS, DT, SIGMA, SIGC, KAPPA, THR, TCAP = 20000, 2e-3, 0.3, 1.2, 1.0, 0.95, 60.0
TWO, THREE = np.array([0.8, 0.2]), np.array([0.64, 0.29, 0.07])
SURFACE, TEXT, TEXT2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e6e5e1"
SLOT = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#7a5cc7"]
plt.rcParams.update({"font.family": "sans-serif", "font.size": 9, "axes.titlesize": 10, "axes.titleweight": "bold",
    "axes.edgecolor": GRID, "axes.labelcolor": TEXT2, "xtick.color": TEXT2, "ytick.color": TEXT2, "axes.titlecolor": TEXT,
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE, "axes.grid": True, "grid.color": GRID,
    "grid.linewidth": 0.8, "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False, "legend.fontsize": 8, "axes.axisbelow": True})

def model_A(e0):
    e = np.tile(e0, (TRIALS, 1)); alive = np.ones_like(e, bool); done = np.zeros(TRIALS, bool)
    winner = np.full(TRIALS, -1); esurv = np.zeros(TRIALS)
    for _ in range(200000):
        act = ~done
        if not act.any(): break
        de = SIGMA * np.sqrt(np.maximum(e, 0)) * np.sqrt(DT) * rng.standard_normal(e.shape) * alive
        e = e + de * act[:, None]
        dead = e <= 0; e[dead] = 0.0; alive &= ~dead
        n_alive = alive.sum(1)
        fin = act & (n_alive <= 1)
        if fin.any():
            w = np.where(n_alive[fin] == 1, alive[fin].argmax(1), -1)
            winner[fin] = w; esurv[fin] = e[fin].sum(1); done[fin] = True
    return winner, esurv, done

def pair_exchange(e0, kernel, stop):
    n = e0.size; pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    e = np.tile(e0, (TRIALS, 1)); done = np.zeros(TRIALS, bool); winner = np.full(TRIALS, -1); maxdev = 0.0
    steps = 0
    while not done.all() and steps < int(TCAP / DT):
        steps += 1
        for (i, j) in pairs:
            d = kernel(e[:, i], e[:, j]) * np.sqrt(DT) * rng.standard_normal(TRIALS)
            d = np.clip(d, -e[:, i], e[:, j])  # cap at the giver's energy
            d = d * (~done)
            e[:, i] += d; e[:, j] -= d
        maxdev = max(maxdev, np.abs(e.sum(1) - 1).max())
        w = stop(e, steps)
        newly = (~done) & (w >= 0); winner[newly] = w[newly]; done |= newly
    return winner, done, maxdev, steps

def stop_absorb(e, steps):
    m = e.max(1); return np.where(m >= 1 - 1e-6, e.argmax(1), -1)

def make_stop_rate():
    def stop(e, steps):
        u = rng.uniform(0, 1, e.shape); fire = u < KAPPA * e * DT
        any_ = fire.any(1); first = np.where(any_, fire.argmax(1), -1)
        # if two fire in one step, pick proportional to share among them
        multi = fire.sum(1) > 1
        if multi.any():
            for k in np.where(multi)[0]:
                idx = np.where(fire[k])[0]; first[k] = rng.choice(idx, p=e[k, idx] / e[k, idx].sum())
        return first
    return stop

def stop_threshold(e, steps):
    m = e.max(1); return np.where(m >= THR, e.argmax(1), -1)

def model_D(e0, delta=0.02):
    n = e0.size; e = np.tile(e0, (TRIALS, 1)); done = np.zeros(TRIALS, bool); winner = np.full(TRIALS, -1)
    for _ in range(2000000):
        if done.all(): break
        i = rng.integers(0, n, TRIALS); j = (i + rng.integers(1, n, TRIALS)) % n  # giver i, receiver j
        d = np.minimum(delta, e[np.arange(TRIALS), i]) * (~done)
        e[np.arange(TRIALS), i] -= d; e[np.arange(TRIALS), j] += d
        m = e.max(1); fin = (~done) & (m >= 1 - 1e-9); winner[fin] = e[fin].argmax(1); done |= fin
    return winner, done

def main():
    t0 = time.time(); lines = ["Verification (PREREGISTRATION.md §3)", ""]; rows = []; se = np.sqrt(0.8 * 0.2 / TRIALS)
    def rec(model, config, winner, extra=None):
        n = TWO.size if config == "two" else THREE.size
        probs = [(winner == k).mean() for k in range(n)]
        row = dict(model=model, config=config, **{f"P_site{k+1}": probs[k] for k in range(n)}, undecided=(winner < 0).mean())
        if extra: row.update(extra)
        rows.append(row); return probs
    # A
    wA, esurv, doneA = model_A(TWO); pA = rec("A independent CIR", "two", wA, dict(esurv_median=np.median(esurv), esurv_below_half=(esurv < 0.5).mean(), esurv_above_1p5=(esurv > 1.5).mean()))
    ok = abs(pA[0] - 0.8) < 3 * se; lines.append(f"A, two sites: P(site 1) = {pA[0]:.4f} against 0.800 (s.e. {se:.4f}) -> {'PASS' if ok else 'FAIL'}")
    lines.append(f"    A survivor energy: median {np.median(esurv):.3f}, below 0.5 in {(esurv < 0.5).mean():.3f}, above 1.5 in {(esurv > 1.5).mean():.3f} (registered about 0.7, 0.34, 0.15)")
    with open(OUT / "survivor_energy.csv", "w") as f:
        f.write("esurv\n"); f.writelines(f"{v:.5g}\n" for v in esurv)
    wA3, _, _ = model_A(THREE); rec("A independent CIR", "three", wA3)
    # B
    wB, doneB, devB, stB = pair_exchange(TWO, lambda a, b: SIGMA * np.sqrt(np.maximum(a * b, 0)), stop_absorb); pB = rec("B geometric-mean kernel", "two", wB, dict(max_energy_deviation=devB))
    ok = abs(pB[0] - 0.8) < 3 * se and devB < 1e-9 and doneB.all()
    lines.append(f"B, two sites: P(site 1) = {pB[0]:.4f} against 0.800, max |E-1| = {devB:.1e}, terminated {doneB.mean():.4f} -> {'PASS' if ok else 'FAIL'}")
    wB3, _, _, _ = pair_exchange(THREE, lambda a, b: SIGMA * np.sqrt(np.maximum(a * b, 0)), stop_absorb); rec("B geometric-mean kernel", "three", wB3)
    # C
    wC1, dC1, devC1, _ = pair_exchange(TWO, lambda a, b: SIGC * a * b, make_stop_rate()); pC1 = rec("C multiplicative kernel, rate-linear commit", "two", wC1, dict(max_energy_deviation=devC1))
    ok = abs(pC1[0] - 0.8) < 3 * se; lines.append(f"C1, two sites, rate-linear commit: P(site 1) = {pC1[0]:.4f} against 0.800 -> {'PASS' if ok else 'FAIL'}")
    wC2, dC2, devC2, _ = pair_exchange(TWO, lambda a, b: SIGC * a * b, stop_threshold); pC2 = rec("C multiplicative kernel, threshold 0.95", "two", wC2, dict(reached_cap=(wC2 < 0).mean()))
    closed = (0.8 - (1 - THR)) / (2 * THR - 1)
    lines.append(f"    C2, threshold 0.95: P(site 1) = {pC2[0]:.4f} against closed form {closed:.4f} (registered ±0.010), cap reached {(wC2 < 0).mean():.4f}")
    # D
    wD3, _ = model_D(THREE); pD3 = rec("D fixed parcels, reflecting zero", "three", wD3)
    ok = pD3[0] < 0.54; lines.append(f"D, three sites: P = ({pD3[0]:.3f}, {pD3[1]:.3f}, {pD3[2]:.3f}) against Born (0.64, 0.29, 0.07), Paper 1's (0.46, 0.27, 0.26); P(site 1) < 0.54 -> {'PASS' if ok else 'FAIL'}")
    wD2, _ = model_D(TWO); pD2 = rec("D fixed parcels, reflecting zero", "two", wD2)
    lines.append(f"    D, two sites: P(site 1) = {pD2[0]:.4f} (registered below 0.75)")
    keys = sorted({k for r in rows for k in r}, key=lambda k: (k not in ("model", "config"), k))
    with open(OUT / "fairness.csv", "w", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=keys); wr.writeheader()
        for r in rows: wr.writerow({k: (f"{v:.5g}" if isinstance(v, float) else v) for k, v in r.items()})
    fig, axes = plt.subplots(1, 2, figsize=(9.8, 3.9))
    ax = axes[0]; labels = ["A independent\nCIR", "B geometric-\nmean kernel", "C1 multiplicative,\nrate commit", "C2 multiplicative,\nthreshold 0.95", "D parcels,\nreflecting zero"]
    vals = [pA[0], pB[0], pC1[0], pC2[0], pD2[0]]
    ax.bar(range(5), vals, color=SLOT); ax.axhline(0.8, color=TEXT, lw=0.9); ax.text(-0.4, 0.81, "Born 0.800", fontsize=8)
    ax.set_xticks(range(5)); ax.set_xticklabels(labels, fontsize=7.5); ax.set_ylim(0.5, 0.9); ax.set_ylabel("P(site 1), two sites (0.8, 0.2)")
    ax.set_title("Fairness by model and termination")
    ax = axes[1]; ax.hist(np.clip(esurv, 0, 3), bins=60, color=SLOT[0]); ax.axvline(1, color=TEXT, lw=0.9); ax.text(1.05, ax.get_ylim()[1] * 0.9, "one quantum", fontsize=8)
    ax.set_xlabel("survivor's energy at absorption, model A"); ax.set_ylabel("trials"); ax.set_title("Independent noise does not conserve the pot")
    fig.tight_layout(); fig.savefig(OUT / "fig1_fairness.png", dpi=150); plt.close(fig)
    text = "\n".join(lines) + "\n"; (OUT / "verification.txt").write_text(text); print(text)
    ok = all("FAIL" not in l for l in lines)
    (OUT / "run_info.txt").write_text(f"exchange_fairness.py, run {time.strftime('%Y-%m-%d %H:%M:%S')}, seed {SEED}, {TRIALS} trials per cell, {time.time() - t0:.1f} s, numpy {np.__version__}\nverification: {'PASS' if ok else 'FAIL'}\n")
    print(f"{time.time() - t0:.1f} s")

if __name__ == "__main__":
    main()
