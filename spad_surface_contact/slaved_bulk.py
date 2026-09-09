#!/usr/bin/env python3
"""
slaved_bulk.py: Part F (addendum 4 of PREREGISTRATION.md). The bulk rebuilt on Part E's slaved
second spinor: every atom presents ±90° by the side of contact and one amplitude; the incident
particle brings its own offset and amplitude; least rotation, the weight with unequal amplitudes.
"""
import csv, sys, time
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, str(Path(__file__).resolve().parent))
import spad_contact as sc
from spad_contact import PIX, M, N_PRIMARY, SLOT, TEXT, TEXT2, flight, region_of, breakdown

OUT = sc.OUT
rng = np.random.default_rng(20260908 + 3); sc.rng = rng
n = N_PRIMARY
A_LIST = [0.064, 0.0012]; B_LIST = [0.0, 0.025, 0.098, "a"]; DP_LIST = [0.0, 90.0, "uniform"]; THETAS = [0.45, 0.6]

def shallow_well_a_of_d():
    NX, DX, M_ = 1024, 0.05, 1.0
    X = (np.arange(NX) - NX // 2) * DX; P = 2 * np.pi * np.fft.fftfreq(NX, d=DX)
    F = np.fft.fft(np.eye(NX), axis=0); Pop = np.linalg.inv(F) @ (P[:, None] * F)
    V = -0.05 * np.exp(-X ** 2 / (2 * 16.0))
    H = np.block([[np.diag(M_ + V), Pop], [Pop, np.diag(-M_ + V)]]); w, v = np.linalg.eigh(H)
    idx = np.nonzero((w > -M_ + 1e-6) & (w < M_ - 1e-6))[0]; k = idx[np.argmin(w[idx])]; Eb = w[k]
    psi = v[:, k].reshape(2, NX); psi = psi * np.exp(-1j * np.angle(psi[0][NX // 2]))
    rows = []
    for d in (0.5, 1, 2, 4, 6, 8, 12, 16):
        i = NX // 2 + int(round(d / DX))
        rows.append(dict(d=d, a=np.abs(psi[1][i]) / np.abs(psi[0][i]), phase_deg=np.degrees(np.angle(psi[1][i] / psi[0][i])), phi_over_max=np.abs(psi[0][i]) / np.abs(psi[0]).max()))
    return rows, np.sqrt(M_ ** 2 - Eb ** 2) / (Eb + M_)

def main():
    t0 = time.time(); lines = ["Verification, Part F (PREREGISTRATION.md, addendum 4)", ""]
    x0, ts, xs, rho_d = flight(4.0)
    x_arr = xs[-1]; reg = region_of(x_arr); inpix = reg > 0
    xa = x_arr[inpix]; rg = reg[inpix]; m = xa.size
    a_lat = 1.0 / (n - 1)
    base = np.round(xa / a_lat) * a_lat
    atoms = base[:, None] + ((np.arange(n) - n // 2) * a_lat)[None, :]
    lo = np.where(rg == 1, PIX[1][0], PIX[2][0])[:, None]; hi = np.where(rg == 1, PIX[1][1], PIX[2][1])[:, None]
    inside = (atoms >= lo - 1e-9) & (atoms <= hi + 1e-9)
    side = np.where(xa[:, None] > atoms, np.pi / 2, -np.pi / 2)          # D23: +90° if the arrival is to the atom's right
    dist = np.abs(xa[:, None] - atoms)
    bd = breakdown(m) >= sc.Q_TH
    share_arr = np.mean(rg == 1)
    dp_uniform = rng.uniform(-np.pi, np.pi, m)
    rows = []
    for a in A_LIST:
        for bspec in B_LIST:
            b = a if bspec == "a" else bspec
            for dps in DP_LIST:
                dp = dp_uniform if dps == "uniform" else np.full(m, np.radians(dps))
                mis = np.abs(np.angle(np.exp(1j * (dp[:, None] - side))))
                key = np.where(inside, mis + 1e-9 * dist, np.inf)          # D24: ties broken by the nearest atom
                k = key.argmin(axis=1); idx = np.arange(m)
                rot = mis[idx, k]
                w = (a * a + b * b + 2 * a * b * np.cos(dp - side[idx, k])) / (a + b) ** 2 if (a + b) > 0 else np.ones(m)
                w_closed = (a * a + b * b + 2 * a * b * np.abs(np.sin(dp))) / (a + b) ** 2 if (a + b) > 0 else np.ones(m)
                left = np.mean(side[idx, k] > 0)
                for theta in THETAS:
                    pair = w >= theta - 1e-12; rec = pair & bd
                    share = np.mean(rg[rec] == 1) if rec.sum() else np.nan
                    rows.append(dict(a=a, b=b, b_spec=bspec, delta_p=dps, theta=theta, rotation_mean_deg=np.degrees(rot.mean()), w_mean=w.mean(),
                                     w_closed_mean=w_closed.mean(), side_plus90_fraction=left, pair_rate=pair.mean(), efficiency=rec.sum() / M,
                                     share_pixel1=share, share_se=np.sqrt(share * (1 - share) / rec.sum()) if rec.sum() else np.nan, arrival_share=share_arr))
                    if theta == 0.45:
                        if dps == 0.0:
                            ok = abs(rot.max() - np.pi / 2) < 1e-12 and abs(rot.min() - np.pi / 2) < 1e-12
                            lines.append(f"a = {a}, b = {b:.4g}, δ_p = 0: rotation {np.degrees(rot.mean()):.6f}° (90 exactly) -> {'PASS' if ok else 'FAIL'}; w {w.mean():.4f} against closed form {w_closed.mean():.4f} -> {'PASS' if abs(w - w_closed).max() < 1e-12 else 'FAIL'}")
                        elif dps == 90.0:
                            ok = rot.max() < 1e-12
                            lines.append(f"a = {a}, b = {b:.4g}, δ_p = 90: rotation {np.degrees(rot.max()):.2e}° (0 exactly) -> {'PASS' if ok else 'FAIL'}; left side chosen {left:.3f}")
                        else:
                            se = np.degrees(rot.std() / np.sqrt(m)); ok = abs(np.degrees(rot.mean()) - 45) < 3 * se
                            lines.append(f"a = {a}, b = {b:.4g}, δ_p uniform: rotation mean {np.degrees(rot.mean()):.2f}° against 45 (s.e. {se:.2f}) -> {'PASS' if ok else 'FAIL'}")
                        if pair.mean() == 1.0:
                            se = np.sqrt(share_arr * (1 - share_arr) / rec.sum()); ok = abs(share - share_arr) < 3 * se
                            lines.append(f"    share {share:.4f} against the arrival share {share_arr:.4f} (s.e. {se:.4f}) -> {'PASS' if ok else 'FAIL'}")
    aod, tail = shallow_well_a_of_d()
    pred = {0.5: 0.011, 1: 0.022, 2: 0.043, 4: 0.078, 6: 0.100, 8: 0.112}
    for r in aod:
        if r["d"] in pred:
            ok = abs(r["a"] - pred[r["d"]]) < 0.001
            lines.append(f"a(d = {r['d']:g}) = {r['a']:.4f} against {pred[r['d']]:.3f}, phase {r['phase_deg']:.1f}° -> {'PASS' if ok else 'FAIL'}")
    lines.append(f"tail limit κ/(E_b + m) = {tail:.4f} (predicted 0.119)")
    fmt = lambda v: (f"{v:.6g}" if isinstance(v, float) else v)
    with open(OUT / "slaved_bulk.csv", "w", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=list(rows[0].keys())); wr.writeheader()
        for r in rows: wr.writerow({k: fmt(v) for k, v in r.items()})
    with open(OUT / "slaved_a_of_d.csv", "w", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=list(aod[0].keys())); wr.writeheader()
        for r in aod: wr.writerow({k: fmt(v) for k, v in r.items()})
    # figure
    fig, axes = plt.subplots(1, 2, figsize=(9.6, 3.9))
    ax = axes[0]; dps_ = np.linspace(-180, 180, 361); s = np.abs(np.sin(np.radians(dps_)))
    for (a, b, lab), c in zip(((0.064, 0.098, "a = 0.064, b = 0.098 (electron, p = 0.2 m)"), (0.064, 0.064, "a = b (equal amplitudes)"), (0.0012, 0.098, "silicon a = 0.0012, b = 0.098"), (0.064, 0.0, "b = 0 (photon)")), SLOT):
        ax.plot(dps_, (a * a + b * b + 2 * a * b * s) / (a + b) ** 2, color=c, lw=1.6, label=lab)
    ax.axhline(0.45, color=TEXT2, lw=0.8, ls="--"); ax.text(-175, 0.46, "gap ratio 0.45 (500 nm)", fontsize=7.5, color=TEXT2)
    ax.axhline(0.6, color=TEXT2, lw=0.8, ls=":"); ax.text(-175, 0.61, "0.6 (660 nm)", fontsize=7.5, color=TEXT2)
    ax.set_ylim(0.3, 1.05); ax.set_xlabel("incident particle's own offset δ_p (degrees); 0 = running wave"); ax.set_ylabel("weight of the chosen atom")
    ax.set_title("Slaved bulk: the best side gives [a² + b² + 2ab|sin δ_p|]/(a + b)²"); ax.legend(fontsize=7, loc="lower right")
    ax = axes[1]
    ax.plot([r["d"] for r in aod], [r["a"] for r in aod], "o-", color=SLOT[0], label="a(d), shallow well")
    ax.axhline(tail, color=TEXT2, lw=0.8, ls="--"); ax.text(9, tail + 0.004, f"tail limit κ/(E_b + m) = {tail:.3f}", fontsize=8, color=TEXT2)
    ax.set_xlabel("distance of the contact from the nucleus d (Compton lengths)"); ax.set_ylabel("second-spinor amplitude |χ|/|φ|")
    ax.set_title("What amplitude the atom presents at the contact"); ax.legend(loc="lower right")
    fig.tight_layout(); fig.savefig(OUT / "fig9_slaved.png", dpi=150); plt.close(fig)
    text = "\n".join(lines) + "\n"; (OUT / "verification_slaved.txt").write_text(text); print(text)
    print("rows (θ = 0.45 and 0.6):")
    for r in rows:
        print(f"  a={r['a']:<7g} b={r['b']:<7.4g} δ_p={str(r['delta_p']):<8} θ={r['theta']}: rot {r['rotation_mean_deg']:6.2f}° w {r['w_mean']:.4f} left {r['side_plus90_fraction']:.3f} pair {r['pair_rate']:.3f} eff {r['efficiency']:.3f} share {r['share_pixel1']:.4f}")
    ok = all("FAIL" not in l for l in lines)
    with open(OUT / "run_info.txt", "a") as f:
        f.write(f"slaved_bulk.py, run {time.strftime('%Y-%m-%d %H:%M:%S')}, seed 20260911, {time.time() - t0:.1f} s\nverification (Part F): {'PASS' if ok else 'FAIL'}\n")
    print(f"{time.time() - t0:.1f} s")
    return 0 if ok else 1

if __name__ == "__main__":
    sys.exit(main())
