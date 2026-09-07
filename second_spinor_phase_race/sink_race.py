#!/usr/bin/env python3
"""
sink_race.py: Part G, two sinks racing, with reversible charging, a latch, and vacuum noise
==========================================================================================

Definitions D19 to D23 of PREREGISTRATION.md (Part G addendum). Same split-operator Dirac scheme
as sink_flow.py. Each site: a reversible charge C_i fed by the norm in its well (dC/dt = N − C/τ_c),
an Ornstein–Uhlenbeck noise η_i on it, a well U_i of width L_i, and, once C_i + η_i reaches C_th, a
latching absorber Γ_i. E_i is the norm the site has absorbed.
"""
import csv, sys, warnings
from pathlib import Path
import numpy as np
warnings.filterwarnings("ignore", message="Mean of empty slice")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
from bohm_dirac_flight import X, DX, gauss, normalise, to_p, psi_at, SURFACE, TEXT, TEXT2, GRID
from sink_flow import kinetic_factors, kinetic, neg_fraction

HERE = Path(__file__).resolve().parent
OUT = HERE / "results"
DT, T_G = 0.01, 20.0
W, U0, G0, KAPPA, L0, U_CAP = 0.75, 0.5, 1.0, 10.0, 0.05, 1.5
C_TH, TAU_C, TAU_N = 0.03, 10.0, 1.0
SLOT = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#7a5cc7", "#c2417a"]

def state(r):
    g = np.sqrt(r) * gauss(X, 4.0, 1.0) + gauss(X, -4.0, 1.0)
    return psi_at(to_p(normalise(np.array([g, np.zeros_like(g)], dtype=complex))), 4.0)

def run(r, off2, s, lam, rng, record_every=0.1, u_cap=U_CAP):
    psi = state(r); K = kinetic_factors(DT)
    xs = [4.0, -4.0 + off2]
    C = [0.0, 0.0]; eta = [0.0, 0.0]; Eabs = [0.0, 0.0]; latched = [False, False]; t_latch = [np.nan, np.nan]
    n_steps = int(round(T_G / DT)); every = int(round(record_every / DT))
    ts, Cs, Es, Ns, negs = [], [], [], [], []
    for k in range(n_steps + 1):
        t = k * DT
        rho = np.abs(psi[0]) ** 2 + np.abs(psi[1]) ** 2
        f = [np.exp(-(X - xs[i]) ** 2 / (2 * (L0 + lam * Eabs[i]) ** 2)) for i in range(2)]
        N = [(rho * f[i]).sum() * DX for i in range(2)]
        if k % every == 0:
            ts.append(t); Cs.append(C.copy()); Es.append(Eabs.copy()); Ns.append(N.copy())
            negs.append(neg_fraction(psi) if k % (10 * every) == 0 else (negs[-1] if negs else neg_fraction(psi)))
        if k == n_steps:
            break
        for i in range(2):
            eta[i] += -eta[i] * DT / TAU_N + s * np.sqrt(DT) * rng.standard_normal()
            C[i] += (N[i] - C[i] / TAU_C) * DT
            if not latched[i] and C[i] + eta[i] >= C_TH:
                latched[i] = True; t_latch[i] = t
        U = [min(u_cap, U0 * W * (1 + KAPPA * (C[i] + Eabs[i]))) for i in range(2)]
        G = [G0 * W * (1 + KAPPA * (C[i] + Eabs[i])) if latched[i] else 0.0 for i in range(2)]
        phase = sum(U[i] * f[i] for i in range(2)); damp = sum(G[i] * f[i] for i in range(2))
        half = np.exp(1j * phase * DT / 2 - damp * DT / 2)
        n_before = rho.sum() * DX
        psi = kinetic(psi * half, K) * half
        loss = n_before - (np.abs(psi[0]) ** 2 + np.abs(psi[1]) ** 2).sum() * DX
        if loss > 0 and damp.max() > 0:
            wts = [(rho * G[i] * f[i]).sum() for i in range(2)]; tot = sum(wts)
            for i in range(2):
                Eabs[i] += loss * wts[i] / tot
    Cs, Es, Ns, negs = map(np.array, (Cs, Es, Ns, negs))
    remain = 1 - Es.sum(axis=1)
    negs = np.where(remain > 0.05, negs, 0.0)
    return dict(t=np.array(ts), C=Cs, E=Es, N=Ns, neg=negs, t_latch=t_latch, latched=latched)

def outcome(res):
    tl = res["t_latch"]; both = all(res["latched"])
    if not any(res["latched"]):
        winner = -1
    elif both and abs(tl[0] - tl[1]) < DT / 2:
        winner = 2                                   # tie
    else:
        winner = int(np.nanargmin(np.where(np.isnan(tl), np.inf, tl)))
    loser = 1 - winner if winner in (0, 1) else None
    rev = np.nan
    if loser is not None and res["C"][:, loser].max() > 0:
        rev = res["C"][-1, loser] / res["C"][:, loser].max()
    return winner, both, rev

SETTINGS = [(1.0, off, s, lam) for off in (0.0, 0.5, 1.0) for s in (0.0, 0.005, 0.02) for lam in (0.0, 10.0)] + \
           [(4.0, 0.0, s, lam) for s in (0.0, 0.005, 0.02) for lam in (0.0, 10.0)]

def main():
    rng = np.random.default_rng(20260907 + 7)
    lines = ["Verification and outcomes, Part G (PREREGISTRATION.md §6, Part G addendum: two sinks racing)", ""]
    trials = []; summary = []; keep = {}
    for (r, off, s, lam) in SETTINGS:
        n_tr = 1 if s == 0 else 10
        outs = []
        for j in range(n_tr):
            res = run(r, off, s, lam, rng)
            w, both, rev = outcome(res)
            outs.append((w, both, rev, res["t_latch"], res["E"][-1], res["neg"].max()))
            trials.append([r, off, s, lam, j, w, int(both), f"{rev:.3f}", f"{res['t_latch'][0]:.2f}", f"{res['t_latch'][1]:.2f}",
                           f"{res['E'][-1, 0]:.4f}", f"{res['E'][-1, 1]:.4f}", f"{res['neg'].max():.4f}"])
            if j == 0 and (r, off, s, lam) in [(4.0, 0.0, 0.0, 10.0), (4.0, 0.0, 0.0, 0.0), (1.0, 0.0, 0.005, 10.0), (1.0, 1.0, 0.005, 10.0)]:
                keep[(r, off, s, lam)] = res
        p1 = np.mean([o[0] == 0 for o in outs]); p2 = np.mean([o[0] == 1 for o in outs]); pt = np.mean([o[0] == 2 for o in outs])
        pn = np.mean([o[0] == -1 for o in outs]); pboth = np.mean([o[1] for o in outs])
        revs = [o[2] for o in outs if not np.isnan(o[2])]
        tl1 = np.nanmean([o[3][0] for o in outs]); tl2 = np.nanmean([o[3][1] for o in outs])
        summary.append([r, off, s, lam, n_tr, f"{p1:.2f}", f"{p2:.2f}", f"{pt:.2f}", f"{pn:.2f}", f"{pboth:.2f}",
                        f"{np.mean(revs):.3f}" if revs else "", f"{tl1:.2f}" if not np.isnan(tl1) else "never", f"{tl2:.2f}" if not np.isnan(tl2) else "never",
                        f"{np.mean([o[4][0] for o in outs]):.3f}", f"{np.mean([o[4][1] for o in outs]):.3f}", f"{max(o[5] for o in outs):.4f}"])
        lines.append(f"Part G: r = {r:g}, offset {off:g}, noise {s:g}, λ = {lam:g} ({n_tr} trial{'s' if n_tr > 1 else ''}): site 1 first {p1:.2f}, site 2 first {p2:.2f}, tie {pt:.2f}, none {pn:.2f};"
                     f" both latch {pboth:.2f}; loser reversal {np.mean(revs):.2f}" if revs else
                     f"Part G: r = {r:g}, offset {off:g}, noise {s:g}, λ = {lam:g} ({n_tr} trial{'s' if n_tr > 1 else ''}): site 1 first {p1:.2f}, site 2 first {p2:.2f}, tie {pt:.2f}, none {pn:.2f}; both latch {pboth:.2f}")
    kmax = max(float(r[-1]) for r in summary)
    n_flag = sum(1 for t in trials if float(t[-1]) > 0.10)
    lines.append(f"Part G: Klein diagnostic max over all runs while the remainder exceeds 0.05: {kmax:.4f}; trials above the 0.10 flag: {n_flag} of {len(trials)}"
                 f"   -> {'no flag' if n_flag == 0 else 'FLAGGED (prediction G5 missed); see the cap sensitivity'}")
    # cap sensitivity: the deterministic channel settings re-run with the well capped at 1.0 mc²
    lines.append("Part G, cap sensitivity (well capped at 1.0 mc² instead of 1.5), deterministic channel settings:")
    cap_ok = True
    for (r, off, s, lam) in [(1.0, 0.0, 0.0, 10.0), (1.0, 0.5, 0.0, 10.0), (1.0, 1.0, 0.0, 10.0), (4.0, 0.0, 0.0, 10.0)]:
        res15 = run(r, off, s, lam, np.random.default_rng(0)); res10 = run(r, off, s, lam, np.random.default_rng(0), u_cap=1.0)
        o15, o10 = outcome(res15), outcome(res10)
        same = (o15[0] == o10[0]) and (o15[1] == o10[1])
        cap_ok = cap_ok and same
        lines.append(f"  r = {r:g}, offset {off:g}: cap 1.5 -> winner {o15[0]}, both {int(o15[1])}, reversal {o15[2]:.2f}, Klein {res15['neg'].max():.3f};"
                     f" cap 1.0 -> winner {o10[0]}, both {int(o10[1])}, reversal {o10[2]:.2f}, Klein {res10['neg'].max():.3f}   -> {'same outcome' if same else 'OUTCOME CHANGED'}")
    lines.append(f"Part G: outcomes independent of the cap in the flagged regime -> {'PASS' if cap_ok else 'FAIL'}")
    with open(OUT / "part_g_race_summary.csv", "w", newline="") as fh:
        wr = csv.writer(fh)
        wr.writerow(["r", "site2_offset", "noise_s", "lambda", "trials", "P_site1_first", "P_site2_first", "P_tie", "P_none", "P_both_latch",
                     "loser_reversal_ratio_mean", "t_latch_site1_mean", "t_latch_site2_mean", "absorbed_site1_t20", "absorbed_site2_t20", "neg_energy_fraction_max"])
        wr.writerows(summary)
    with open(OUT / "part_g_race_trials.csv", "w", newline="") as fh:
        wr = csv.writer(fh)
        wr.writerow(["r", "site2_offset", "noise_s", "lambda", "trial", "winner_0_1_tie2_none-1", "both_latched", "loser_reversal_ratio",
                     "t_latch_site1", "t_latch_site2", "absorbed_site1_t20", "absorbed_site2_t20", "neg_energy_fraction_max"])
        wr.writerows(trials)
    text = "\n".join(lines) + "\n"
    with open(OUT / "verification.txt", "a") as fh:
        fh.write("\n" + text)
    with open(OUT / "run_info.txt", "a") as fh:
        fh.write(f"verification (sink_race.py): {'PASS' if cap_ok else 'FAIL'}; Klein flag in {n_flag} trials (max {kmax:.3f})\n")
    print(text)
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
    for key, ls, lab in [((4.0, 0.0, 0.0, 10.0), "-", "channel λ = 10"), ((4.0, 0.0, 0.0, 0.0), "--", "no channel")]:
        res = keep[key]
        ax.plot(res["t"], res["C"][:, 0], color=SLOT[0], lw=1.5, ls=ls, label=f"site 1 charge, {lab}")
        ax.plot(res["t"], res["C"][:, 1], color=SLOT[1], lw=1.5, ls=ls, label=f"site 2 charge, {lab}")
        ax.plot(res["t"], res["E"][:, 0], color=SLOT[2], lw=1.2, ls=ls, label=f"site 1 absorbed, {lab}")
    ax.axhline(C_TH, color=TEXT2, lw=0.8, alpha=0.6); ax.text(19.5, C_TH + 0.02, "threshold", color=TEXT2, fontsize=7.5, ha="right")
    ax.set_yscale("log"); ax.set_ylim(1e-3, 1.5); ax.set_xlabel("time after the deposit (Compton times)"); ax.set_ylabel("charge, absorbed norm")
    ax.set_title("r = 4, no noise: latch and reversal"); ax.legend(loc="lower right", ncol=1)
    ax = axes[1]
    for key, c, lab in [((1.0, 0.0, 0.005, 10.0), SLOT[0], "offset 0"), ((1.0, 1.0, 0.005, 10.0), SLOT[1], "offset 1.0")]:
        res = keep[key]
        ax.plot(res["t"], res["C"][:, 0], color=c, lw=1.5, label=f"site 1 charge, {lab}")
        ax.plot(res["t"], res["C"][:, 1], color=c, lw=1.5, ls="--", label=f"site 2 charge, {lab}")
    ax.axhline(C_TH, color=TEXT2, lw=0.8, alpha=0.6)
    ax.set_xlabel("time after the deposit"); ax.set_ylabel("charge"); ax.set_ylim(0, 0.06)
    ax.set_title("r = 1, noise 0.005, channel: one trial each"); ax.legend(loc="lower right")
    ax = axes[2]
    for s, c in zip((0.0, 0.005, 0.02), SLOT[:3]):
        for lam, ls, mk in ((10.0, "-", "o"), (0.0, ":", "s")):
            sub = [row for row in summary if row[0] == 1.0 and row[2] == s and row[3] == lam]
            ax.plot([row[1] for row in sub], [float(row[5]) + 0.5 * float(row[7]) for row in sub], ls, marker=mk, color=c, lw=1.4, ms=5,
                    label=f"noise {s:g}, {'channel' if lam else 'no channel'}")
    ax.set_xlabel("site 2 offset toward the middle (Compton lengths)"); ax.set_ylabel("P(site 1 first) + ½ P(tie)"); ax.set_ylim(0.3, 1.05)
    ax.set_title("r = 1: who latches first"); ax.legend(loc="upper left", ncol=2)
    fig.tight_layout(); fig.savefig(OUT / "fig9_race.png", dpi=150); plt.close(fig)
    return 0 if cap_ok else 1

if __name__ == "__main__":
    sys.exit(main())
