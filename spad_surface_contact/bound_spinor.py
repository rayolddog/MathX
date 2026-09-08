#!/usr/bin/env python3
"""
bound_spinor.py: Part E (addendum 3 of PREREGISTRATION.md). A Dirac electron bound in a well in the
record's 1+1-D conventions; the amplitude, phase and rotation of its second spinor in both readings
(the lower row; the negative-energy free component); the deposit weight with unequal amplitudes.
"""
import csv, sys, time
from pathlib import Path
import numpy as np
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

NX, DX, M_ = 1024, 0.05, 1.0
X = (np.arange(NX) - NX // 2) * DX
P = 2 * np.pi * np.fft.fftfreq(NX, d=DX)
E = np.sqrt(P ** 2 + M_ ** 2)
WELLS = [(0.05, 4.0), (0.2, 2.0), (0.6, 1.0), (1.0, 0.5)]
PRED = {(0.05, 4.0): (0.9721, 0.0041, 1.5e-6), (0.2, 2.0): (0.8888, 0.0167, 6.8e-5), (0.6, 1.0): (0.6979, 0.0537, 1.3e-3), (1.0, 0.5): (0.6543, 0.0904, 4.7e-3)}

def momentum_operator():
    F = np.fft.fft(np.eye(NX), axis=0); return np.linalg.inv(F) @ (P[:, None] * F)

def free_pm(psi_p):
    nrm = np.sqrt(2 * E * (E + M_)); up = np.array([(E + M_) / nrm, P / nrm]); um = np.array([-P / nrm, (E + M_) / nrm])
    return (up.conj() * psi_p).sum(axis=0), (um.conj() * psi_p).sum(axis=0)

def main():
    t0 = time.time(); lines = ["Verification, Part E (PREREGISTRATION.md, addendum 3)", ""]
    Pop = momentum_operator(); rows = []; states = {}
    for U0, L in WELLS:
        V = -U0 * np.exp(-X ** 2 / (2 * L * L))
        H = np.block([[np.diag(M_ + V), Pop], [Pop, np.diag(-M_ + V)]])
        w, v = np.linalg.eigh(H)
        idx = np.nonzero((w > -M_ + 1e-6) & (w < M_ - 1e-6))[0]
        k = idx[np.argmin(w[idx])]; Eb = w[k]; vec = v[:, k]; psi = vec.reshape(2, NX)
        resid = np.abs(H @ vec - Eb * vec).max()
        # make the upper row real and positive at the centre (a global phase)
        psi = psi * np.exp(-1j * np.angle(psi[0][NX // 2]))
        fL = (np.abs(psi[1]) ** 2).sum() / (np.abs(psi) ** 2).sum()
        dphi = np.fft.ifft(1j * P * np.fft.fft(psi[0])); chi_pred = -1j * dphi / (Eb + M_ - V)
        chi_err = np.abs(psi[1] - chi_pred).max() / np.abs(psi[1]).max()
        ok = (np.abs(psi[0]) > 1e-3 * np.abs(psi[0]).max()) & (np.abs(psi[1]) > 1e-3 * np.abs(psi[1]).max())
        ph = np.degrees(np.angle(psi[1][ok] / psi[0][ok]))
        left = ph[X[ok] < 0]; right = ph[X[ok] > 0]
        cp, cm = free_pm(np.fft.fft(psi, axis=1))
        fm = (np.abs(cm) ** 2).sum() / ((np.abs(cp) ** 2).sum() + (np.abs(cm) ** 2).sum())
        wgt = np.abs(cp) ** 2 + np.abs(cm) ** 2; keep = wgt > 1e-6 * wgt.max()
        ww = wgt[keep] / wgt[keep].sum(); z = (ww * np.exp(1j * np.angle(cm / cp)[keep])).sum()
        imag_upper = np.abs(psi[0].imag).max() / np.abs(psi[0]).max()
        rows.append(dict(U0=U0, L=L, binding=M_ - Eb, E_b=Eb, n_bound=int(idx.size), lower_fraction=fL, a_lower=np.sqrt(fL),
                         chi_phase_left_deg=np.median(left), chi_phase_right_deg=np.median(right), chi_phase_spread_deg=max(np.ptp(left), np.ptp(right)),
                         upper_imag_over_max=imag_upper, negE_fraction=fm, a_negE=np.sqrt(fm), pm_offset_mean_deg=np.degrees(np.angle(z)), pm_offset_resultant=abs(z),
                         eigen_residual=resid, chi_identity_error=chi_err))
        states[(U0, L)] = (psi, V, Eb)
        pe, pf, pm = PRED[(U0, L)]
        lines.append(f"well U0 = {U0:g}, L = {L:g}: eigenvector residual {resid:.1e} -> {'PASS' if resid < 1e-9 else 'FAIL'}; χ against −iφ′/(E_b+m−V): {chi_err:.1e} -> {'PASS' if chi_err < 1e-10 else 'FAIL'};"
                     f" E_b {Eb:.4f} against {pe} -> {'PASS' if abs(Eb - pe) < 1e-4 + 5e-5 else 'FAIL'}; f_L {fL:.4f} against {pf} -> {'PASS' if abs(fL - pf) < 1e-4 + 5e-5 else 'FAIL'};"
                     f" f₋ {fm:.2e} against {pm:.1e} -> {'PASS' if abs(fm - pm) / pm < 0.06 else 'FAIL'}")
    # visibility table
    vis = []
    for r in rows:
        for reading, a in (("lower row", r["a_lower"]), ("negative-energy component", r["a_negE"])):
            for b in (0.0, a, 0.1, 0.5, 1.0):
                V = 2 * a * b / (a * a + b * b) if (a * a + b * b) > 0 else 0.0
                vis.append(dict(U0=r["U0"], L=r["L"], binding=r["binding"], reading=reading, a=a, b=b, visibility=V, w_min=((a - b) / (a + b)) ** 2 if a + b > 0 else np.nan))
    # extrapolation to silicon (binding 1e-5): power laws through the two shallowest wells
    b1, b2 = rows[0]["binding"], rows[1]["binding"]
    sL = np.log(rows[1]["lower_fraction"] / rows[0]["lower_fraction"]) / np.log(b2 / b1)
    sM = np.log(rows[1]["negE_fraction"] / rows[0]["negE_fraction"]) / np.log(b2 / b1)
    fL_si = rows[0]["lower_fraction"] * (1e-5 / b1) ** sL; fm_si = rows[0]["negE_fraction"] * (1e-5 / b1) ** sM
    lines.append(f"scaling through the two shallowest wells: f_L ∝ binding^{sL:.2f}, f₋ ∝ binding^{sM:.2f}; extrapolated to binding 1e-5 (silicon valence): f_L ≈ {fL_si:.1e} (a_L ≈ {np.sqrt(fL_si):.1e}), f₋ ≈ {fm_si:.1e} (a₋ ≈ {np.sqrt(fm_si):.1e}); an order-of-magnitude estimate")
    fmt = lambda v: (f"{v:.6g}" if isinstance(v, float) else v)
    with open(OUT / "bound_spinor.csv", "w", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=list(rows[0].keys())); wr.writeheader()
        for r in rows: wr.writerow({k: fmt(v) for k, v in r.items()})
    with open(OUT / "bound_visibility.csv", "w", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=list(vis[0].keys())); wr.writeheader()
        for r in vis: wr.writerow({k: fmt(v) for k, v in r.items()})
    # figure
    fig, axes = plt.subplots(1, 3, figsize=(13, 3.9))
    ax = axes[0]
    for (U0, L), c in zip((WELLS[0], WELLS[3]), (SLOT[0], SLOT[1])):
        psi, V, Eb = states[(U0, L)]; sc = np.abs(psi[0]).max()
        ax.plot(X, psi[0].real / sc, color=c, lw=1.6, label=f"first spinor φ, U₀ = {U0:g}, L = {L:g}")
        ax.plot(X, psi[1].imag / sc, "--", color=c, lw=1.4, label=f"second spinor, Im χ (Re χ = 0)")
    ax.set_xlim(-12, 12); ax.set_xlabel("x (Compton lengths)"); ax.set_ylabel("amplitude / max φ")
    ax.set_title("A bound Dirac electron: the second spinor is −iφ′/(E+m−V)"); ax.legend(fontsize=7)
    ax = axes[1]
    bind = [r["binding"] for r in rows]
    ax.plot(bind, [r["lower_fraction"] for r in rows], "o-", color=SLOT[0], label="lower-row fraction f_L")
    ax.plot(bind, [r["negE_fraction"] for r in rows], "s-", color=SLOT[1], label="negative-energy fraction f₋")
    ax.plot([1e-5, b1], [fL_si, rows[0]["lower_fraction"]], ":", color=SLOT[0]); ax.plot([1e-5, b1], [fm_si, rows[0]["negE_fraction"]], ":", color=SLOT[1])
    ax.axvline(1e-5, color=TEXT2, lw=0.8, ls="--"); ax.text(1.3e-5, 0.3, "silicon valence\n(binding 10⁻⁵ m)", fontsize=7.5, color=TEXT2)
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("binding energy / m"); ax.set_ylabel("fraction of the norm")
    ax.set_title("How much second spinor a bound electron has"); ax.legend(loc="lower right")
    ax = axes[2]
    bb = np.logspace(-4, 0, 200)
    for r, c in zip(rows, SLOT):
        a = r["a_lower"]; ax.plot(bb, 2 * a * bb / (a * a + bb * bb), color=c, lw=1.6, label=f"lower row, a = {a:.3f} (binding {r['binding']:.2f})")
    a = rows[0]["a_negE"]; ax.plot(bb, 2 * a * bb / (a * a + bb * bb), ":", color=TEXT, lw=1.4, label=f"−E component, shallow well, a = {a:.4f}")
    ax.set_xscale("log"); ax.set_xlabel("incident particle's second-spinor amplitude b"); ax.set_ylabel("visibility 2ab/(a² + b²)")
    ax.set_title("The phase mechanism's visibility with unequal amplitudes"); ax.legend(fontsize=7, loc="upper left")
    fig.tight_layout(); fig.savefig(OUT / "fig8_bound_spinor.png", dpi=150); plt.close(fig)
    text = "\n".join(lines) + "\n"; (OUT / "verification_bound.txt").write_text(text); print(text)
    for r in rows:
        print(f"U0={r['U0']:g} L={r['L']:g}: binding {r['binding']:.4f}, f_L {r['lower_fraction']:.5f} (a {r['a_lower']:.4f}), χ/φ phase left {r['chi_phase_left_deg']:.1f}° right {r['chi_phase_right_deg']:.1f}° spread {r['chi_phase_spread_deg']:.1e}°, Im φ/max {r['upper_imag_over_max']:.1e}, f₋ {r['negE_fraction']:.2e} (a {r['a_negE']:.4f}), ±E offset mean {r['pm_offset_mean_deg']:.1f}° resultant {r['pm_offset_resultant']:.3f}")
    ok = all("FAIL" not in l for l in lines)
    with open(OUT / "run_info.txt", "a") as f:
        f.write(f"bound_spinor.py, run {time.strftime('%Y-%m-%d %H:%M:%S')}, {time.time() - t0:.1f} s\nverification (Part E): {'PASS' if ok else 'FAIL'}\n")
    print(f"{time.time() - t0:.1f} s")
    return 0 if ok else 1

if __name__ == "__main__":
    sys.exit(main())
