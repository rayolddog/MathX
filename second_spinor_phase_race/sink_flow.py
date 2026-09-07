#!/usr/bin/env python3
"""
sink_flow.py: Part F, the winning atom becomes a sink (PREREGISTRATION.md §6, Part F addendum)
==============================================================================================

Split-operator Dirac evolution in the record's conventions (ħ = c = m = 1, H = σ_x p + σ_z m):
exact free propagator over dt between half-steps of the local factor exp(−i V dt/2 − Γ f dt/2),
with V = −U f the well and Γ f the absorber at the winner's site, both growing with the collected
fraction E_coll = 1 − ∫ρ dx. Imports the packet and conventions from bohm_dirac_flight.py.
"""
import csv, sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
from bohm_dirac_flight import X, DX, P, E, M_, gauss, normalise, to_p, psi_at, pm_components, SURFACE, TEXT, TEXT2, GRID

HERE = Path(__file__).resolve().parent
OUT = HERE / "results"
DT, T_F, W_WIN, U_CAP = 0.01, 40.0, 0.75, 1.5
SLOT = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#7a5cc7", "#c2417a", "#3aa7c9", "#8a8a2a"]

def kinetic_factors(dt):
    c, s = np.cos(E * dt), np.sin(E * dt)
    return (c - 1j * s * M_ / E, -1j * s * P / E, -1j * s * P / E, c + 1j * s * M_ / E)

def kinetic(psi, K):
    a, b = np.fft.fft(psi[0]), np.fft.fft(psi[1])
    u = K[0] * a + K[1] * b; d = K[2] * a + K[3] * b
    return np.array([np.fft.ifft(u), np.fft.ifft(d)])

def state(kind):
    """Part C's two-lobe packet at t_d = 4 (narrow) or the wide version."""
    if kind == "narrow":
        g = np.sqrt(4.0) * gauss(X, 4.0, 1.0) + gauss(X, -4.0, 1.0)
    else:
        g = np.sqrt(4.0) * gauss(X, 10.0, 3.0) + gauss(X, -10.0, 3.0)
    psi0 = normalise(np.array([g, np.zeros_like(g)], dtype=complex))
    return psi_at(to_p(psi0), 4.0)

def neg_fraction(psi):
    cp, cm = pm_components(np.fft.fft(psi, axis=1))
    a, b = (np.abs(cp) ** 2).sum(), (np.abs(cm) ** 2).sum()
    return b / (a + b) if a + b > 0 else 0.0

def run(kind, x_s, g0, u0, L0, kappa, lam, record_every=0.5):
    psi = state(kind); rho = np.abs(psi[0]) ** 2 + np.abs(psi[1]) ** 2
    near_side = X > 0 if x_s > 0 else X < 0
    K = kinetic_factors(DT)
    n_steps = int(round(T_F / DT)); every = int(round(record_every / DT))
    E_coll, ts, near, far, negf, Umax, Gmax, Lmax = 0.0, [], [], [], [], 0.0, 0.0, 0.0
    rate0 = None
    for k in range(n_steps + 1):
        if k % every == 0:
            rho = np.abs(psi[0]) ** 2 + np.abs(psi[1]) ** 2
            ts.append(k * DT); near.append(rho[near_side].sum() * DX); far.append(rho[~near_side].sum() * DX)
            negf.append(neg_fraction(psi) if k % (4 * every) == 0 else negf[-1] if negf else neg_fraction(psi))
        if k == n_steps:
            break
        G = g0 * W_WIN * (1 + kappa * E_coll); U = min(U_CAP, u0 * W_WIN * (1 + kappa * E_coll)); L = L0 + lam * E_coll
        Umax, Gmax, Lmax = max(Umax, U), max(Gmax, G), max(Lmax, L)
        f = np.exp(-(X - x_s) ** 2 / (2 * L * L))
        half = np.exp(1j * U * f * DT / 2 - G * f * DT / 2)
        if k == 0:
            rho = np.abs(psi[0]) ** 2 + np.abs(psi[1]) ** 2
            rate0 = G * (rho * f).sum() * DX
        psi = psi * half
        psi = kinetic(psi, K)
        psi = psi * half
        E_coll = 1 - (np.abs(psi[0]) ** 2 + np.abs(psi[1]) ** 2).sum() * DX
    ts, near, far, negf = map(np.array, (ts, near, far, negf))
    coll = 1 - near - far
    negf = np.where(near + far > 0.05, negf, 0.0)      # the diagnostic is meaningless on a vanished remainder
    return dict(t=ts, near=near, far=far, coll=coll, negf=negf, rate0=rate0, Umax=Umax, Gmax=Gmax, Lmax=Lmax)

RUNS = [
    ("base",               "narrow",  4.0, 1.0, 0.5, 0.05, 10.0, 0.0),
    ("gamma0_0.1",         "narrow",  4.0, 0.1, 0.5, 0.05, 10.0, 0.0),
    ("gamma0_10",          "narrow",  4.0, 10.0, 0.5, 0.05, 10.0, 0.0),
    ("no_well",            "narrow",  4.0, 1.0, 0.0, 0.05, 10.0, 0.0),
    ("deep_well_1.5",      "narrow",  4.0, 1.0, 1.5, 0.05, 10.0, 0.0),
    ("width_0.5",          "narrow",  4.0, 1.0, 0.5, 0.5, 10.0, 0.0),
    ("width_2",            "narrow",  4.0, 1.0, 0.5, 2.0, 10.0, 0.0),
    ("no_growth",          "narrow",  4.0, 1.0, 0.5, 0.05, 0.0, 0.0),
    ("channel_lambda_10",  "narrow",  4.0, 1.0, 0.5, 0.05, 10.0, 10.0),
    ("weak_lobe_winner",   "narrow", -4.0, 1.0, 0.5, 0.05, 10.0, 0.0),
    ("wide_lobes_base",    "wide",   10.0, 1.0, 0.5, 0.05, 10.0, 0.0),
    ("wide_lobes_channel", "wide",   10.0, 1.0, 0.5, 0.05, 10.0, 10.0),
]

def t_at(res, key, level, initial):
    """First time the 'key' norm falls to level × its initial value; nan if never."""
    y = res[key]; idx = np.where(y <= level * initial)[0]
    return res["t"][idx[0]] if idx.size else np.nan

def main():
    lines = ["Verification, Part F (PREREGISTRATION.md §6, Part F addendum: the sink)", ""]
    # (a) U = Γ = 0 reproduces the exact propagator; (b) U ≠ 0, Γ = 0 conserves the norm
    psi = state("narrow"); K = kinetic_factors(DT)
    for _ in range(400):
        psi = kinetic(psi, K)
    exact = psi_at(to_p(state("narrow")), 4.0)
    err = np.abs(psi - exact).max()
    lines.append(f"Part F: split-operator with U = Γ = 0 against the exact propagator after 4 Compton times: max |Δψ| = {err:.1e}   -> {'PASS' if err < 1e-9 else 'FAIL'}")
    psi = state("narrow"); f = np.exp(-(X - 4.0) ** 2 / (2 * 0.5 ** 2)); half = np.exp(1j * 1.0 * f * DT / 2)
    for _ in range(400):
        psi = kinetic(psi * half, K) * half
    n = (np.abs(psi) ** 2).sum() * DX
    lines.append(f"Part F: well only (U = 1, L = 0.5), norm after 4 Compton times: |∫ρ − 1| = {abs(n-1):.1e}   -> {'PASS' if abs(n-1) < 1e-9 else 'FAIL'}")
    # the runs
    results = {}; rows = []
    free = {"narrow": run("narrow", 4.0, 0.0, 0.0, 0.05, 0.0, 0.0), "wide": run("wide", 10.0, 0.0, 0.0, 0.05, 0.0, 0.0)}
    for name, kind, xs, g0, u0, L0, kap, lam in RUNS:
        res = run(kind, xs, g0, u0, L0, kap, lam); results[name] = res
        ref = free[kind]; far_ref = ref["far"] if xs > 0 else ref["near"]
        near0 = res["near"][0]
        at = lambda tt: res["coll"][np.searchsorted(res["t"], tt)]
        rows.append([name, kind, xs, g0, u0, L0, kap, lam, f"{near0:.4f}", f"{res['rate0']:.5f}",
                     f"{at(5):.4f}", f"{at(10):.4f}", f"{at(20):.4f}", f"{at(40):.4f}",
                     f"{t_at(res, 'near', 0.5, near0):.1f}", f"{t_at(res, 'near', 0.1, near0):.1f}",
                     f"{res['far'][-1]:.4f}", f"{far_ref[-1]:.4f}", f"{res['negf'].max():.4f}", f"{res['Umax']:.3f}", f"{res['Gmax']:.3f}", f"{res['Lmax']:.3f}"])
        if name == "base":
            rho = np.abs(state("narrow")[0]) ** 2 + np.abs(state("narrow")[1]) ** 2
            pred = 1.0 * W_WIN * rho[np.argmin(np.abs(X - 4.0))] * np.sqrt(2 * np.pi) * 0.05
            lines.append(f"Part F: base initial rate {res['rate0']:.5f} against F1's Γ₀ ρ(x_s) √(2π) L₀ = {pred:.5f}, ratio {res['rate0']/pred:.3f}   -> {'PASS' if abs(res['rate0']/pred - 1) < 0.02 else 'FAIL'}")
        book = np.abs(res["near"] + res["far"] + res["coll"] - 1).max()
        lines.append(f"Part F: {name}: norm bookkeeping |near + far + collected − 1| ≤ {book:.1e}; Klein diagnostic max {res['negf'].max():.4f} {'(flag)' if res['negf'].max() > 0.10 else ''}")
    with open(OUT / "part_f_sink_summary.csv", "w", newline="") as fh:
        wr = csv.writer(fh)
        wr.writerow(["run", "state", "x_s", "gamma0", "u0", "L0", "kappa", "lambda", "near_share_initial", "initial_rate",
                     "collected_t5", "collected_t10", "collected_t20", "collected_t40", "near_half_time", "near_tenth_time",
                     "far_norm_t40", "far_norm_t40_free_reference", "neg_energy_fraction_max_while_norm_gt_0.05", "U_max", "Gamma_max", "L_max"])
        wr.writerows(rows)
    with open(OUT / "part_f_sink_series.csv", "w", newline="") as fh:
        wr = csv.writer(fh)
        names = [r[0] for r in RUNS]
        wr.writerow(["t_after_deposit"] + [f"collected_{n}" for n in names] + ["far_base", "far_free_narrow", "near_base"])
        for i, t in enumerate(results["base"]["t"]):
            wr.writerow([f"{t:.2f}"] + [f"{results[n]['coll'][i]:.5f}" for n in names] +
                        [f"{results['base']['far'][i]:.5f}", f"{free['narrow']['far'][i]:.5f}", f"{results['base']['near'][i]:.5f}"])
    text = "\n".join(lines) + "\n"
    with open(OUT / "verification.txt", "a") as fh:
        fh.write("\n" + text)
    ok = all("FAIL" not in l for l in lines)
    with open(OUT / "run_info.txt", "a") as fh:
        fh.write(f"verification (sink_flow.py): {'PASS' if ok else 'FAIL'}\n")
    print(text)
    print("run                  near0   rate0     E(5)   E(10)  E(20)  E(40)  t½near t10near  far40  far40free  negmax  Umax  Γmax  Lmax")
    for r in rows:
        print(f"{r[0]:20s} {r[8]}  {r[9]}  {r[10]}  {r[11]}  {r[12]}  {r[13]}  {r[14]:>6}  {r[15]:>6}  {r[16]}  {r[17]}     {r[18]}  {r[19]}  {r[20]}  {r[21]}")
    # figure
    plt.rcParams.update({
        "font.family": "sans-serif", "font.size": 9, "axes.titlesize": 10, "axes.titleweight": "bold",
        "axes.edgecolor": GRID, "axes.labelcolor": TEXT2, "xtick.color": TEXT2, "ytick.color": TEXT2,
        "axes.titlecolor": TEXT, "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
        "savefig.facecolor": SURFACE, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.8,
        "grid.linestyle": "-", "axes.spines.top": False, "axes.spines.right": False,
        "legend.frameon": False, "legend.fontsize": 7.5, "axes.axisbelow": True,
    })
    fig, axes = plt.subplots(1, 3, figsize=(13.2, 4.2))
    ax = axes[0]
    for (name, *_), c in zip(RUNS[:9], SLOT + SLOT):
        ax.plot(results[name]["t"], results[name]["coll"], color=c, lw=1.4, label=name.replace("_", " "))
    ax.axhline(results["base"]["near"][0], color=TEXT2, lw=0.8, alpha=0.6); ax.text(39.5, results["base"]["near"][0] - 0.05, "near lobe's share", color=TEXT2, fontsize=7.5, ha="right")
    ax.set_xlabel("time after the deposit (Compton times)"); ax.set_ylabel("collected fraction of the quantum"); ax.set_ylim(0, 1.32)
    ax.set_title("Collection, narrow lobes, winner at +4"); ax.legend(loc="upper center", ncol=3)
    ax = axes[1]
    b = results["base"]; fr = free["narrow"]
    ax.plot(b["t"], b["near"], color=SLOT[0], lw=1.6, label="near side (x > 0), base")
    ax.plot(b["t"], b["far"], color=SLOT[1], lw=1.6, label="far side (x < 0), base")
    ax.plot(fr["t"], fr["far"], color=SLOT[1], lw=1.0, ls="--", label="far side, free reference")
    ax.plot(b["t"], b["coll"], color=SLOT[2], lw=1.6, label="collected")
    ch = results["channel_lambda_10"]
    ax.plot(ch["t"], ch["coll"], color=SLOT[2], lw=1.2, ls=":", label="collected, channel λ = 10")
    ax.plot(ch["t"], ch["far"], color=SLOT[1], lw=1.2, ls=":", label="far side, channel λ = 10")
    ax.set_xlabel("time after the deposit"); ax.set_ylabel("norm"); ax.set_ylim(0, 1.32); ax.set_title("Norm by side: base and the channel run"); ax.legend(loc="upper center", ncol=2)
    ax = axes[2]
    for name, c in zip(("wide_lobes_base", "wide_lobes_channel", "weak_lobe_winner"), SLOT[3:6]):
        ax.plot(results[name]["t"], results[name]["coll"], color=c, lw=1.6, label=name.replace("_", " "))
        ax.axhline(results[name]["near"][0], color=c, lw=0.8, ls="--", alpha=0.6)
    ax.set_xlabel("time after the deposit"); ax.set_ylabel("collected fraction"); ax.set_ylim(0, 1.32)
    ax.set_title("Wide lobes ±10 and the weak-lobe winner (dashed: near share)"); ax.legend(loc="upper center", ncol=3)
    fig.tight_layout(); fig.savefig(OUT / "fig8_sink.png", dpi=150); plt.close(fig)
    return 0 if ok else 1

if __name__ == "__main__":
    sys.exit(main())
