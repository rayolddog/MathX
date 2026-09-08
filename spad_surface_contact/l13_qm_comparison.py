#!/usr/bin/env python3
"""
l13_qm_comparison.py: the test of L13 (L13_TEST.md). Quantum mechanics for a coherently prepared
two-level absorber hit by a second pulse (semiclassical Bloch, numerical against closed form) and by
a single quantum (Jaynes–Cummings, Fock against coherent state), beside the model's deposit rule.
"""
import csv, sys, time, math
from pathlib import Path
import numpy as np
from scipy.linalg import expm
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent; OUT = HERE / "results"
SURFACE, TEXT, TEXT2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e6e5e1"
SLOT = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#7a5cc7"]
plt.rcParams.update({"font.family": "sans-serif", "font.size": 9, "axes.titlesize": 10, "axes.titleweight": "bold",
    "axes.edgecolor": GRID, "axes.labelcolor": TEXT2, "xtick.color": TEXT2, "ytick.color": TEXT2, "axes.titlecolor": TEXT,
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE, "axes.grid": True, "grid.color": GRID,
    "grid.linewidth": 0.8, "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False, "legend.fontsize": 8, "axes.axisbelow": True})

# ---------------------------------------------------------------- two-level, semiclassical
g_, e_ = np.array([1, 0], complex), np.array([0, 1], complex)          # basis (g, e)
sp = np.outer(e_, g_); sm = sp.T.copy()
def pulse(theta, phi):
    return expm(-1j * theta / 2 * (np.exp(-1j * phi) * sp + np.exp(1j * phi) * sm))
def dephase(rho, t):
    r = rho.copy(); r[0, 1] *= np.exp(-t); r[1, 0] *= np.exp(-t); return r
def d_rho_ee_numeric(th1, th2, t, dphi, phi1=0.7):
    rho = np.outer(g_, g_.conj())
    U1 = pulse(th1, phi1); rho = U1 @ rho @ U1.conj().T; rho = dephase(rho, t)
    before = rho[1, 1].real
    U2 = pulse(th2, phi1 + dphi); rho = U2 @ rho @ U2.conj().T
    return rho[1, 1].real - before
def d_rho_ee_closed(th1, th2, t, dphi):
    return 0.5 * (np.cos(th1) * (1 - np.cos(th2)) + np.sin(th1) * np.sin(th2) * np.exp(-t) * np.cos(dphi))
def model_weight(t, dphi):
    return 0.5 * (1 + np.exp(-t) * np.cos(dphi))

# ---------------------------------------------------------------- Jaynes–Cummings
NF = 14
a = np.diag(np.sqrt(np.arange(1, NF)), 1); ad = a.T.copy()
I2, IF = np.eye(2), np.eye(NF)
H_jc = np.kron(sp, a) + np.kron(sm, ad)                                   # g = 1
n_op = np.kron(I2, ad @ a)
def atom_state(th1, phi1):
    return np.cos(th1 / 2) * g_ - 1j * np.exp(-1j * phi1) * np.sin(th1 / 2) * e_
def field_state(kind, phi2, alpha_abs=1.0):
    if kind == "fock":
        v = np.zeros(NF, complex); v[1] = 1.0; return v
    al = alpha_abs * np.exp(1j * phi2); k = np.arange(NF)
    v = np.exp(-abs(al) ** 2 / 2) * al ** k / np.sqrt([math.factorial(int(i)) for i in k]); return v / np.linalg.norm(v)
def jc_net_photon_change(th1, kind, dphi, gt=0.3, phi1=0.7):
    psi = np.kron(atom_state(th1, phi1), field_state(kind, phi1 + dphi))
    U = expm(-1j * H_jc * gt); psi_t = U @ psi
    return (psi_t.conj() @ n_op @ psi_t).real - (psi.conj() @ n_op @ psi).real

def main():
    t0 = time.time(); lines = ["Verification, the test of L13 (L13_TEST.md §3)", ""]
    dphis = np.linspace(0, 2 * np.pi, 73)
    rows = []
    worst = 0.0
    for th1 in (0.3, np.pi / 2):
        for t in (0.0, 0.3, 1.0, 3.0):
            for dphi in dphis:
                num = d_rho_ee_numeric(th1, 0.2, t, dphi); cl = d_rho_ee_closed(th1, 0.2, t, dphi)
                worst = max(worst, abs(num - cl))
                rows.append(dict(theta1=th1, theta2=0.2, t_over_T2=t, dphi_deg=np.degrees(dphi), d_rho_ee_numeric=num, d_rho_ee_closed=cl, model_weight=model_weight(t, dphi)))
    lines.append(f"P14: two-level numeric against closed form, worst |difference| = {worst:.1e} (tolerance 1e-12) -> {'PASS' if worst < 1e-12 else 'FAIL'}")
    # Jaynes–Cummings
    jc = {}
    for kind in ("fock", "coherent"):
        jc[kind] = np.array([jc_net_photon_change(np.pi / 2, kind, d) for d in dphis])
    flat = jc["fock"].max() - jc["fock"].min(); mean_fock = jc["fock"].mean()
    ok = flat < 1e-12 and abs(mean_fock - 0.045) < 0.01
    lines.append(f"P15: Fock |1> on a π/2-prepared atom, gt = 0.3: net photon change {mean_fock:+.4f} (predicted about +0.045), spread over phase {flat:.1e} -> {'PASS' if ok else 'FAIL'}")
    amp = (jc["coherent"].max() - jc["coherent"].min()) / 2
    ok = abs(amp - 0.3) < 0.06
    lines.append(f"P15: coherent |α|² = 1 on the same atom: fringe amplitude {amp:.3f} (predicted about 0.3), mean {jc['coherent'].mean():+.4f} -> {'PASS' if ok else 'FAIL'}")
    # fringe-to-average ratio against probe strength
    th2s = np.array([0.03, 0.1, 0.3, 1.0, 2.0]); th1 = 0.3
    ratio_qm = np.sin(th1) * np.sin(th2s) / (np.cos(th1) * (1 - np.cos(th2s)))
    lines.append("P16: fringe / phase-averaged absorption at θ₁ = 0.3, t = 0: quantum mechanics " + ", ".join(f"θ₂ = {a_:g}: {r:.2f}" for a_, r in zip(th2s, ratio_qm)) + "; the model 1.00 at every θ₂"
                 + f" -> {'PASS' if abs(ratio_qm[1] - 6.2) < 0.2 else 'FAIL'}")
    with open(OUT / "l13_two_level.csv", "w", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=list(rows[0].keys())); wr.writeheader()
        for r in rows: wr.writerow({k: (f"{v:.6f}" if isinstance(v, float) else v) for k, v in r.items()})
    with open(OUT / "l13_jaynes_cummings.csv", "w", newline="") as f:
        wr = csv.writer(f); wr.writerow(["dphi_deg", "fock_net_photon_change", "coherent_net_photon_change"])
        for d, x, y in zip(dphis, jc["fock"], jc["coherent"]): wr.writerow([f"{np.degrees(d):.1f}", f"{x:.6f}", f"{y:.6f}"])
    # figure
    fig, axes = plt.subplots(1, 3, figsize=(13, 3.9)); deg = np.degrees(dphis)
    ax = axes[0]
    for (th1, t, c) in ((np.pi / 2, 0.0, SLOT[0]), (np.pi / 2, 1.0, SLOT[1]), (0.3, 0.0, SLOT[2])):
        cl = d_rho_ee_closed(th1, 0.2, t, dphis); ax.plot(deg, cl / (0.2 ** 2 / 4), color=c, lw=1.6, label=f"QM, θ₁ = {th1:.2g}, t = {t:g} T2")
        ax.plot(deg[::6], np.array([d_rho_ee_numeric(th1, 0.2, t, d) for d in dphis[::6]]) / (0.2 ** 2 / 4), "o", color=c, ms=3.5)
    ax.plot(deg, model_weight(0.0, dphis), "--", color=TEXT, lw=1.4, label="the model, ½[1 + cos Δφ], t = 0")
    ax.plot(deg, model_weight(1.0, dphis), ":", color=TEXT, lw=1.4, label="the model, t = T2")
    ax.axhline(0, color=TEXT2, lw=0.8); ax.set_xlabel("phase of the probe relative to the preparation, Δφ (degrees)"); ax.set_ylabel("deposit / ordinary absorption (θ₂ = 0.2)")
    ax.set_title("Coherent probe on a prepared two-level absorber"); ax.legend(fontsize=7)
    ax = axes[1]
    ax.plot(deg, jc["fock"], color=SLOT[0], lw=1.8, label="Fock |1⟩: no fringe")
    ax.plot(deg, jc["coherent"], color=SLOT[1], lw=1.8, label="coherent, |α|² = 1: fringe")
    ax.plot(deg, -0.09 * model_weight(0.0, dphis) * 2, "--", color=TEXT, lw=1.4, label="the model, one quantum of definite offset (scale arbitrary)")
    ax.axhline(0, color=TEXT2, lw=0.8); ax.set_xlabel("Δφ (degrees)"); ax.set_ylabel("net photon change after gt = 0.3 (negative = absorbed)")
    ax.set_title("One quantum on a π/2-prepared atom (Jaynes–Cummings)"); ax.legend(fontsize=7, loc="lower right")
    ax = axes[2]
    th2f = np.logspace(-2, np.log10(2.5), 60)
    ax.plot(th2f, np.sin(0.3) * np.sin(th2f) / (np.cos(0.3) * (1 - np.cos(th2f))), color=SLOT[0], lw=1.8, label="QM, θ₁ = 0.3: ∝ 1/θ₂ for a weak probe")
    ax.axhline(1.0, color=TEXT, ls="--", lw=1.4, label="the model: R(t) = 1, any probe")
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("probe pulse area θ₂ (∝ field)"); ax.set_ylabel("fringe amplitude / phase-averaged deposit")
    ax.set_title("How the fringe scales with the probe, t = 0"); ax.legend(fontsize=7.5, loc="lower left")
    fig.tight_layout(); fig.savefig(OUT / "fig7_l13.png", dpi=150); plt.close(fig)
    text = "\n".join(lines) + "\n"; (OUT / "verification_l13.txt").write_text(text); print(text)
    ok = all("FAIL" not in l for l in lines)
    with open(OUT / "run_info.txt", "a") as f:
        f.write(f"l13_qm_comparison.py, run {time.strftime('%Y-%m-%d %H:%M:%S')}, {time.time() - t0:.1f} s\nverification (L13 test): {'PASS' if ok else 'FAIL'}\n")
    return 0 if ok else 1

if __name__ == "__main__":
    sys.exit(main())
