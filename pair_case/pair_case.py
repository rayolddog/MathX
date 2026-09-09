#!/usr/bin/env python3
"""
pair_case.py: JB's model at two analyzers (PREREGISTRATION.md of this folder). A shared hidden
polarization angle, a local outcome rule at each polarizer (least rotation, Malus), the quantum
comparator, and the model's phase-gated detection. CHSH and the correlation function.
"""
import csv, sys, time
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent; OUT = HERE / "results"; OUT.mkdir(exist_ok=True)
rng = np.random.default_rng(20260912)
N = 400000
A, AP, B, BP = 0.0, np.pi / 4, np.pi / 8, 3 * np.pi / 8
THETAS = np.radians(np.arange(0, 90.1, 7.5))
GATES = [0.0, 0.5, 0.6, 0.63, 0.75, 0.85, 0.9]
SURFACE, TEXT, TEXT2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e6e5e1"
SLOT = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#7a5cc7"]
plt.rcParams.update({"font.family": "sans-serif", "font.size": 9, "axes.titlesize": 10, "axes.titleweight": "bold",
    "axes.edgecolor": GRID, "axes.labelcolor": TEXT2, "xtick.color": TEXT2, "ytick.color": TEXT2, "axes.titlecolor": TEXT,
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE, "axes.grid": True, "grid.color": GRID,
    "grid.linewidth": 0.8, "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False, "legend.fontsize": 8, "axes.axisbelow": True})

def outcome(phi, a, rule):
    d = phi - a
    if rule == "least":
        return np.where(np.cos(2 * d) >= 0, 1, -1)
    return np.where(rng.uniform(0, 1, phi.size) < np.cos(d) ** 2, 1, -1)

def pair_run(a, b, rule, gate=0.0, shared=True):
    """Returns E among coincidences, single-side efficiency, coincidence fraction."""
    if rule == "qm":
        c = np.cos(2 * (a - b)); sa = rng.choice([1, -1], N); u = rng.uniform(0, 1, N)
        sb = np.where(u < 0.5 * (1 + c), sa, -sa)
        return np.mean(sa * sb), 1.0, 1.0
    phi = rng.uniform(0, np.pi, N); phi2 = phi if shared else rng.uniform(0, np.pi, N)
    sa, sb = outcome(phi, a, rule), outcome(phi2, b, rule)
    if gate > 0:
        da = np.where(sa == 1, np.cos(phi - a) ** 2, np.sin(phi - a) ** 2) >= gate
        db = np.where(sb == 1, np.cos(phi2 - b) ** 2, np.sin(phi2 - b) ** 2) >= gate
    else:
        da = db = np.ones(N, bool)
    co = da & db
    return np.mean(sa[co] * sb[co]), da.mean(), co.mean()

def chsh(rule, gate=0.0, shared=True):
    e = {k: pair_run(*v, rule, gate, shared) for k, v in dict(ab=(A, B), abp=(A, BP), apb=(AP, B), apbp=(AP, BP)).items()}
    S = e["ab"][0] - e["abp"][0] + e["apb"][0] + e["apbp"][0]
    eta1 = np.mean([x[1] for x in e.values()]); pc = np.mean([x[2] for x in e.values()])
    return S, eta1, pc, e["ab"][0]

def main():
    t0 = time.time(); lines = ["Verification (PREREGISTRATION.md §4)", ""]
    se_E = 1 / np.sqrt(N); se_S = 2 * se_E
    corr = []
    for th in THETAS:
        el = pair_run(0.0, th, "least")[0]; em = pair_run(0.0, th, "malus")[0]; eq = pair_run(0.0, th, "qm")[0]; ei = pair_run(0.0, th, "least", shared=False)[0]
        cl = 1 - 4 * th / np.pi; cm = 0.5 * np.cos(2 * th); cq = np.cos(2 * th)
        corr.append(dict(theta_deg=np.degrees(th), E_least=el, E_least_closed=cl, E_malus=em, E_malus_closed=cm, E_qm=eq, E_qm_closed=cq, E_independent=ei))
        ok = abs(el - cl) < 3 * se_E and abs(em - cm) < 3 * se_E and abs(eq - cq) < 3 * se_E and abs(ei) < 3 * se_E
        lines.append(f"θ = {np.degrees(th):5.1f}°: least {el:+.4f} ({cl:+.4f}), Malus {em:+.4f} ({cm:+.4f}), QM {eq:+.4f} ({cq:+.4f}), independent {ei:+.4f} (0) -> {'PASS' if ok else 'FAIL'}")
    rows_s = []
    for rule, pred in (("least", 2.0), ("malus", np.sqrt(2)), ("qm", 2 * np.sqrt(2))):
        S, eta1, pc, e1 = chsh(rule); ok = abs(S - pred) < 3 * se_S
        rows_s.append(dict(rule=rule, detection="ideal", S=S, S_predicted=pred, E_ab=e1, eta1=eta1, coincidence=pc))
        lines.append(f"CHSH {rule}, ideal: S = {S:.4f} against {pred:.4f} (s.e. {se_S:.4f}) -> {'PASS' if ok else 'FAIL'}")
    S, eta1, pc, e1 = chsh("least", shared=False); ok = abs(S) < 3 * se_S
    rows_s.append(dict(rule="least", detection="ideal, independent angles", S=S, S_predicted=0.0, E_ab=e1, eta1=eta1, coincidence=pc))
    lines.append(f"CHSH least, independent angles: S = {S:+.4f} against 0 -> {'PASS' if ok else 'FAIL'}")
    rows_g = []
    for rule in ("least", "malus"):
        S0 = chsh(rule)[0]
        for g in GATES:
            S, eta1, pc, e1 = chsh(rule, g); eta_c = pc / eta1 if eta1 > 0 else np.nan
            gm = min(4 / eta_c - 2, 4.0) if eta_c > 0 else np.nan
            rows_g.append(dict(rule=rule, gate=g, S_coincidence=S, E_ab=e1, eta_single=eta1, coincidence_fraction=pc, eta_conditional=eta_c, garg_mermin_bound=gm))
            if rule == "least" and g == 0.5:
                ok = abs(S - S0) < 1e-12 and abs(eta1 - 1) < 1e-12
                lines.append(f"gate 0.5 leaves least unchanged: S = {S:.6f} against {S0:.6f}, η₁ = {eta1:.6f} -> {'PASS' if ok else 'FAIL'}")
            if rule == "least" and g == 0.6:
                ok = abs(S - gm) < 0.02
                lines.append(f"Garg–Mermin identity, least at g = 0.6: S_c = {S:.4f} against 4/η_c − 2 = {gm:.4f} (η_c = {eta_c:.4f}) -> {'PASS' if ok else 'FAIL'}")
    # diagnostic, not registered (added 2026-09-08 after the run): the registered gate-0.5 check compared two
    # independent Monte Carlo draws; here the same draw is gated and ungated.
    phi = rng.uniform(0, np.pi, N); sa, sb = outcome(phi, A, "least"), outcome(phi, B, "least")
    da = np.where(sa == 1, np.cos(phi - A) ** 2, np.sin(phi - A) ** 2) >= 0.5; db = np.where(sb == 1, np.cos(phi - B) ** 2, np.sin(phi - B) ** 2) >= 0.5
    lines.append(f"    diagnostic, not registered (added after the run): on one draw, gate 0.5 detects {np.mean(da & db):.6f} of pairs and changes E(a, b) by {abs(np.mean(sa * sb) - np.mean(sa[da & db] * sb[da & db])):.1e}")
    fmt = lambda v: (f"{v:.5g}" if isinstance(v, float) else v)
    for name, rows in (("correlation.csv", corr), ("chsh.csv", rows_s), ("gated.csv", rows_g)):
        with open(OUT / name, "w", newline="") as f:
            wr = csv.DictWriter(f, fieldnames=list(rows[0].keys())); wr.writeheader()
            for r in rows: wr.writerow({k: fmt(v) for k, v in r.items()})
    # figure
    fig, axes = plt.subplots(1, 2, figsize=(9.8, 3.9))
    ax = axes[0]; td = np.degrees(THETAS)
    ax.plot(td, [c["E_least"] for c in corr], "o-", color=SLOT[0], label="least rotation (JB's projection): 1 − 4θ/π")
    ax.plot(td, [c["E_malus"] for c in corr], "s-", color=SLOT[1], label="Malus, sampled (R2): ½ cos 2θ")
    ax.plot(td, [c["E_qm"] for c in corr], "^-", color=TEXT, label="quantum mechanics: cos 2θ")
    ax.plot(td, [c["E_independent"] for c in corr], ":", color=TEXT2, label="independent angles")
    ax.set_xlabel("analyzer difference θ = a − b (degrees)"); ax.set_ylabel("correlation E(θ)"); ax.set_title("A shared hidden polarization at two analyzers"); ax.legend(fontsize=7.5)
    ax = axes[1]
    for rule, c, mk in (("least", SLOT[0], "o"), ("malus", SLOT[1], "s")):
        sel = [r for r in rows_g if r["rule"] == rule and r["gate"] > 0]
        ax.plot([r["eta_conditional"] for r in sel], [r["S_coincidence"] for r in sel], mk + "-", color=c, label=f"{'least rotation' if rule == 'least' else 'Malus'}, gated detection")
    ee = np.linspace(0.55, 1.0, 100); ax.plot(ee, np.minimum(4 / ee - 2, 4), "--", color=TEXT2, label="Garg–Mermin bound 4/η − 2")
    ax.axhline(2, color=TEXT, lw=0.9); ax.text(0.56, 2.05, "Bell's bound 2", fontsize=8)
    ax.axhline(2 * np.sqrt(2), color=TEXT, lw=0.9, ls=":"); ax.text(0.56, 2.88, "quantum 2√2", fontsize=8)
    ax.axvline(2 / (1 + np.sqrt(2)), color=SLOT[2], lw=0.9); ax.text(0.835, 3.5, "η = 0.828", color=SLOT[2], fontsize=8)
    ax.set_xlabel("conditional detection efficiency η_c"); ax.set_ylabel("CHSH among coincidences"); ax.set_ylim(1.2, 4.2)
    ax.set_title("The model's only route above 2 is the detection loophole"); ax.legend(fontsize=7.5, loc="upper right")
    fig.tight_layout(); fig.savefig(OUT / "fig1_pair.png", dpi=150); plt.close(fig)
    text = "\n".join(lines) + "\n"; (OUT / "verification.txt").write_text(text); print(text)
    print("gated:")
    for r in rows_g: print(f"  {r['rule']:<5} g = {r['gate']:<4}: S_c {r['S_coincidence']:.3f}, η₁ {r['eta_single']:.3f}, η_c {r['eta_conditional']:.3f}, bound {r['garg_mermin_bound']:.3f}")
    ok = all("FAIL" not in l for l in lines)
    (OUT / "run_info.txt").write_text(f"pair_case.py, run {time.strftime('%Y-%m-%d %H:%M:%S')}, seed 20260912, N = {N} pairs per setting, {time.time() - t0:.1f} s, numpy {np.__version__}\nverification: {'PASS' if ok else 'FAIL'}\n")
    print(f"{time.time() - t0:.1f} s")
    return 0 if ok else 1

if __name__ == "__main__":
    sys.exit(main())
