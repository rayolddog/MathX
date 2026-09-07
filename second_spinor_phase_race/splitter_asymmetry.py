#!/usr/bin/env python3
"""
splitter_asymmetry.py: Part H, an uneven per-event split from a spinor energy asymmetry, fed into
the two-sink race of Part G (PREREGISTRATION.md, Part H addendum).
"""
import csv, sys, warnings
from pathlib import Path
import numpy as np
warnings.filterwarnings("ignore", message="Mean of empty slice")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
import sink_race as G
from bohm_dirac_flight import SURFACE, TEXT, TEXT2, GRID

HERE = Path(__file__).resolve().parent
OUT = HERE / "results"
EPS_LIST = [0.0, 0.05, 0.1, 0.2]
FORMS = ["switch", "smooth"]
CTH_LIST = [0.015, 0.03, 0.06]
N_EV, NOISE, LAM = 10, 0.005, 10.0
SLOT = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]

def split_fraction(eps, form, rng):
    d = rng.uniform(-np.pi, np.pi)
    g = np.sign(np.cos(d)) if form == "switch" else np.cos(d)
    return 0.5 + 0.5 * eps * g

def main():
    slow = len(sys.argv) > 1 and sys.argv[1] == "slow"
    rng = np.random.default_rng(20260907 + (17 if slow else 13))
    G.T_G = 40.0 if slow else 20.0
    if slow:
        G.TAU_C = 40.0
    cth_list = [0.06, 0.12] if slow else CTH_LIST
    tag = "_slow" if slow else ""
    trials, summary = [], []
    lines = [f"Part H{' addendum, slow charging (τ_c = 40, 40 Compton times)' if slow else ''} (PREREGISTRATION.md, Part H addendum): uneven split from a spinor asymmetry, into the two-sink race", ""]
    for form in FORMS:
        for eps in EPS_LIST:
            for cth in cth_list:
                G.C_TH = cth
                outs = []
                for j in range(N_EV):
                    f = split_fraction(eps, form, rng)
                    r = f / (1 - f)                                  # bright lobe on the right (site 1) when f > ½
                    res = G.run(r, 0.0, NOISE, LAM, rng)
                    w, both, rev = G.outcome(res)
                    big_side = 0 if f >= 0.5 else 1                  # site index holding the larger share
                    big_first = (w == big_side) if w in (0, 1) else (0.5 if w == 2 else np.nan)
                    tl = res["t_latch"]; gap = abs(tl[0] - tl[1]) if not any(np.isnan(tl)) else np.nan
                    vis = 2 * np.sqrt(f * (1 - f))
                    outs.append((both, big_first, vis, gap, w))
                    trials.append([form, eps, cth, j, f"{f:.4f}", w, int(both), f"{big_first}" if not np.isnan(big_first) else "", f"{vis:.4f}",
                                   f"{gap:.2f}" if not np.isnan(gap) else "", f"{tl[0]:.2f}", f"{tl[1]:.2f}"])
                dbl = np.mean([o[0] for o in outs]); bf = np.nanmean([o[1] for o in outs]); vis = np.mean([o[2] for o in outs])
                gap = np.nanmean([o[3] for o in outs]); none = np.mean([o[4] == -1 for o in outs])
                summary.append([form, eps, cth, N_EV, f"{dbl:.2f}", f"{bf:.2f}", f"{vis:.4f}", f"{gap:.2f}" if not np.isnan(gap) else "", f"{none:.2f}"])
                lines.append(f"  {form:6s} ε = {eps:<4g} C_th = {cth:<5g}: double record {dbl:.2f}, larger share first {bf:.2f}, visibility {vis:.4f}, latch gap {gap:.2f}, none {none:.2f}")
    with open(OUT / f"part_h_splitter{tag}.csv", "w", newline="") as fh:
        wr = csv.writer(fh); wr.writerow(["form", "eps", "C_th", "events", "P_double_record", "P_larger_share_first", "visibility_mean", "latch_gap_mean", "P_none"]); wr.writerows(summary)
    with open(OUT / f"part_h_splitter_trials{tag}.csv", "w", newline="") as fh:
        wr = csv.writer(fh); wr.writerow(["form", "eps", "C_th", "event", "f", "winner_0_1_tie2_none-1", "both_latched", "larger_share_first", "visibility", "latch_gap", "t_latch_site1", "t_latch_site2"]); wr.writerows(trials)
    text = "\n".join(lines) + "\n"
    with open(OUT / "verification.txt", "a") as fh:
        fh.write("\n" + text)
    with open(OUT / "run_info.txt", "a") as fh:
        fh.write(f"Part H{tag} (splitter_asymmetry.py): run; no pass/fail criterion beyond Part G's, which it reuses\n")
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
    fig, axes = plt.subplots(1, 3, figsize=(13.2, 4.0))
    for ax, key, lab in zip(axes[:2], ("P_double_record", "P_larger_share_first"), ("double-record rate", "larger share's side recorded first")):
        for form, ls in (("switch", "-"), ("smooth", "--")):
            for cth, c in zip(cth_list, SLOT):
                sub = [s for s in summary if s[0] == form and s[2] == cth]
                idx = 4 if key == "P_double_record" else 5
                ax.plot([s[1] for s in sub], [float(s[idx]) for s in sub], ls, marker="o", color=c, lw=1.4, ms=4, label=f"{form}, C_th = {cth:g}")
        ax.set_xlabel("asymmetry ε"); ax.set_ylabel(lab); ax.set_ylim(-0.02, 1.05); ax.set_title(lab.capitalize())
        if key == "P_double_record":
            ax.axhline(0.01, color=TEXT2, lw=0.8, alpha=0.6); ax.text(0.2, 0.03, "anticorrelation bound 0.01", color=TEXT2, fontsize=7.5, ha="right")
        ax.legend(loc="best", ncol=2)
    ax = axes[2]
    e = np.linspace(0, 0.3, 200)
    ax.plot(e, np.sqrt(1 - e ** 2), color=SLOT[0], lw=1.6, label="switch form, √(1 − ε²)")
    ax.plot(e, [np.mean(np.sqrt(1 - (x * np.cos(np.linspace(-np.pi, np.pi, 2001))) ** 2)) for x in e], color=SLOT[1], lw=1.6, label="smooth form")
    for form, mk in (("switch", "o"), ("smooth", "s")):
        sub = [s for s in summary if s[0] == form and s[2] == cth_list[1 if slow else 1]]
        ax.plot([s[1] for s in sub], [float(s[6]) for s in sub], mk, color=TEXT, ms=4, mfc=SURFACE)
    ax.axhline(0.98, color=TEXT2, lw=0.8, alpha=0.6); ax.text(0.29, 0.975, "1986 single-photon fringes, 0.98", color=TEXT2, fontsize=7.5, ha="right", va="top")
    ax.set_xlabel("asymmetry ε"); ax.set_ylabel("fringe visibility the split allows"); ax.set_ylim(0.95, 1.002); ax.set_title("Visibility against ε"); ax.legend(loc="lower left")
    fig.suptitle("Part H" + (" addendum: slow charging, τ_c = 40" if slow else ": τ_c = 10"), color=TEXT, fontsize=10, fontweight="bold")
    fig.tight_layout(); fig.savefig(OUT / f"fig10_splitter{tag}.png", dpi=150); plt.close(fig)
    return 0

if __name__ == "__main__":
    sys.exit(main())
