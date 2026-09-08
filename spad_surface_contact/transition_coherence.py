#!/usr/bin/env python3
"""
transition_coherence.py: Part D of the 2026-09-08 addenda to PREREGISTRATION.md. The bulk's
second-spinor phase as a transition coherence prepared at t = 0 and dephased thermally (T2 as a
phase random walk, T2* as static site offsets); one pixel (D17) and two pixels (D18).
"""
import csv, sys, time
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, str(Path(__file__).resolve().parent))
import spad_contact as sc
from spad_contact import PIX, M, THETA, N_PRIMARY, SLOT, TEXT2, flight, region_of, breakdown

OUT = sc.OUT
rng = np.random.default_rng(20260908 + 2)
sc.rng = rng
T_LIST = [0.0, 0.03, 0.1, 0.3, 1.0, 3.0, 10.0]
CASES = {"hom": lambda t: np.sqrt(2 * t), "inh": lambda t: t, "both": lambda t: np.sqrt(2 * t + t * t)}
DPS = [0, 45, 90, 135, 180]
PHI0 = 2 * np.arccos(np.sqrt(THETA))
n = N_PRIMARY
N_SINGLE = 20000

def least(delta, dp, inside=None):
    mis = np.abs(np.angle(np.exp(1j * (dp - delta))))
    if inside is not None:
        mis = np.where(inside, mis, np.inf)
    rot = mis.min(axis=1)
    return rot, np.cos(rot / 2) ** 2

def reach(x_arr, reg):
    a = 1.0 / (n - 1)
    base = np.round(x_arr / a) * a
    atoms = base[:, None] + ((np.arange(n) - n // 2) * a)[None, :]
    lo = np.where(reg == 1, PIX[1][0], PIX[2][0])[:, None]; hi = np.where(reg == 1, PIX[1][1], PIX[2][1])[:, None]
    return (atoms >= lo - 1e-9) & (atoms <= hi + 1e-9)

def main():
    t0 = time.time(); lines = ["Verification, Part D (PREREGISTRATION.md, addendum 2)", ""]
    rows_s, rows_p = [], []
    # ------------------------------------------------------------------ D17: one pixel
    for case, sfun in CASES.items():
        for t in T_LIST:
            s = sfun(t)
            for dp in DPS + ["uniform"]:
                Z = rng.standard_normal((N_SINGLE, n)); delta = s * Z
                dpv = rng.uniform(-np.pi, np.pi, (N_SINGLE, 1)) if dp == "uniform" else np.radians(dp)
                rot, w = least(delta, dpv)
                R = np.abs(np.exp(1j * delta).mean())
                rows_s.append(dict(case=case, t_over_T2=t, s=s, R_measured=R, R_predicted=np.exp(-s * s / 2), delta_p=dp,
                                   rotation_mean_deg=np.degrees(rot.mean()), rotation_se_deg=np.degrees(rot.std() / np.sqrt(N_SINGLE)),
                                   pair_rate=np.mean(w >= THETA), winner_w_mean=w.mean()))
                if case == "hom" and dp == 0:
                    ok = abs(R - np.exp(-s * s / 2)) < 0.005
                    lines.append(f"hom, t = {t:g}: R measured {R:.4f} against exp(−s²/2) = {np.exp(-s * s / 2):.4f} -> {'PASS' if ok else 'FAIL'}")
                if case == "hom" and t == 0.0 and dp != "uniform":
                    err = np.max(np.abs(rot - np.radians(dp))); ok = err < 1e-12
                    lines.append(f"hom, t = 0, δ_p = {dp}°: max |rotation − δ_p| = {err:.1e} -> {'PASS' if ok else 'FAIL'}")
                if case == "hom" and t == 10.0 and dp == 0:
                    pred = np.pi / (n + 1); se = rot.std() / np.sqrt(N_SINGLE); ok = abs(rot.mean() - pred) < 3 * se
                    lines.append(f"hom, t = 10, δ_p = 0: rotation {np.degrees(rot.mean()):.2f}° against π/22 = {np.degrees(pred):.2f}° (s.e. {np.degrees(se):.2f}°) -> {'PASS' if ok else 'FAIL'}")
                if case == "hom" and t == 0.0 and dp == "uniform":
                    pr = np.mean(w >= THETA); se = np.sqrt(0.5325 * 0.4675 / N_SINGLE); ok = abs(pr - 0.5325) < 3 * se
                    lines.append(f"hom, t = 0, uniform δ_p: pair rate {pr:.4f} against 0.5325 (s.e. {se:.4f}) -> {'PASS' if ok else 'FAIL'}")
    # ------------------------------------------------------------------ D18: two pixels
    x0, ts, xs, rho_d = flight(4.0)
    x_arr = xs[-1]; reg = region_of(x_arr); inpix = reg > 0
    ii = np.nonzero(inpix)[0]; m = ii.size; rg = reg[ii]
    inside = reach(x_arr[ii], rg)
    bd = breakdown(m) >= sc.Q_TH
    share_arr = np.mean(rg == 1)
    for case, sfun in CASES.items():
        for t in T_LIST:
            s = sfun(t)
            for prep in ("same", "opposite"):
                mu = np.where(rg == 1, 0.0, 0.0 if prep == "same" else np.pi)
                for dpm in ("fixed0", "uniform"):
                    dpv = np.zeros((m, 1)) if dpm == "fixed0" else rng.uniform(-np.pi, np.pi, (m, 1))
                    delta = mu[:, None] + s * rng.standard_normal((m, n))
                    rot, w = least(delta, dpv, inside)
                    pair = w >= THETA; rec = pair & bd
                    p1 = np.mean(pair[rg == 1]); p2 = np.mean(pair[rg == 2])
                    share = np.mean(rg[rec] == 1) if rec.sum() else np.nan
                    rows_p.append(dict(case=case, t_over_T2=t, s=s, preparation=prep, delta_p=dpm, pair_rate_pixel1=p1, pair_rate_pixel2=p2,
                                       records=int(rec.sum()), efficiency=rec.sum() / M, share_pixel1=share,
                                       share_se=np.sqrt(share * (1 - share) / rec.sum()) if rec.sum() else np.nan, arrival_share=share_arr, born=0.8))
                    if case == "hom" and t == 10.0 and prep == "opposite" and dpm == "fixed0":
                        se = np.sqrt(share_arr * (1 - share_arr) / rec.sum()); ok = abs(share - share_arr) < 3 * se
                        lines.append(f"two pixels, hom, t = 10, opposite, δ_p = 0: share {share:.4f} against the arrival share {share_arr:.4f} (s.e. {se:.4f}) -> {'PASS' if ok else 'FAIL'}")
    fmt = lambda v: (f"{v:.4f}" if isinstance(v, float) else v)
    for name, rows in (("dephasing_single.csv", rows_s), ("dephasing_pixels.csv", rows_p)):
        with open(OUT / name, "w", newline="") as f:
            wr = csv.DictWriter(f, fieldnames=list(rows[0].keys())); wr.writeheader()
            for row in rows: wr.writerow({k: fmt(v) for k, v in row.items()})
    # ------------------------------------------------------------------ figure
    fig, axes = plt.subplots(1, 3, figsize=(13, 3.9))
    tt = np.array(T_LIST); tplot = np.where(tt == 0, 0.01, tt)
    ax = axes[0]
    for dp, c in zip((0, 90, 180), SLOT):
        sel = [x for x in rows_s if x["case"] == "hom" and x["delta_p"] == dp]
        ax.plot(tplot, [max(x["rotation_mean_deg"], 0.05) for x in sel], "o-", color=c, label=f"δ_p = {dp}°")
    ax.axhline(np.degrees(np.pi / 22), color=TEXT2, lw=0.8, ls="--"); ax.text(0.012, 9.5, "random bulk, π/22", color=TEXT2, fontsize=8)
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("delay since preparation, t/T2 (0 drawn at 0.01)"); ax.set_ylabel("rotation used (degrees)")
    ax.set_title("Least rotation as the coherence dephases (T2 only)"); ax.legend(loc="lower right")
    ax = axes[1]
    for case, c in zip(CASES, SLOT):
        sel = [x for x in rows_s if x["case"] == case and x["delta_p"] == "uniform"]
        ax.plot(tplot, [x["pair_rate"] for x in sel], "o-", color=c, label=f"pair rate, {case}")
        ax.plot(tplot, [x["R_measured"] for x in sel], ":", color=c, lw=1.2)
    ax.set_xscale("log"); ax.set_xlabel("t/T2"); ax.set_ylabel("pair rate (solid); R(t) (dotted)"); ax.set_ylim(-0.02, 1.05)
    ax.set_title("Pair rate for a particle of random offset, D4b θ = 0.45"); ax.legend(loc="center right")
    ax = axes[2]
    for case, c in zip(CASES, SLOT):
        sel = [x for x in rows_p if x["case"] == case and x["preparation"] == "opposite" and x["delta_p"] == "fixed0"]
        ax.plot(tplot, [x["share_pixel1"] for x in sel], "o-", color=c, label=f"opposite preparation, δ_p = 0, {case}")
    sel = [x for x in rows_p if x["case"] == "hom" and x["preparation"] == "opposite" and x["delta_p"] == "uniform"]
    ax.plot(tplot, [x["share_pixel1"] for x in sel], "s--", color=TEXT2, label="opposite, δ_p uniform")
    ax.axhline(0.8, color="k", lw=0.8); ax.text(0.012, 0.77, "Born 0.800", fontsize=8)
    ax.set_xscale("log"); ax.set_xlabel("t/T2"); ax.set_ylabel("share of records at pixel 1"); ax.set_ylim(0.7, 1.02)
    ax.set_title("Two pixels pulsed in antiphase: the Born deviation dephases"); ax.legend(loc="upper right", fontsize=7.5)
    fig.tight_layout(); fig.savefig(OUT / "fig6_dephasing.png", dpi=150); plt.close(fig)
    text = "\n".join(lines) + "\n"; (OUT / "verification_dephasing.txt").write_text(text); print(text)
    ok = all("FAIL" not in l for l in lines)
    with open(OUT / "run_info.txt", "a") as f:
        f.write(f"transition_coherence.py, run {time.strftime('%Y-%m-%d %H:%M:%S')}, seed 20260910, {time.time() - t0:.1f} s\nverification (Part D): {'PASS' if ok else 'FAIL'}\n")
    print("one pixel, hom:")
    for x in rows_s:
        if x["case"] == "hom" and x["delta_p"] in (0, 90, 180, "uniform"):
            print(f"  t={x['t_over_T2']:5g} R {x['R_measured']:.3f} δ_p={x['delta_p']!s:>7}: rotation {x['rotation_mean_deg']:7.2f}° pair {x['pair_rate']:.4f} w {x['winner_w_mean']:.3f}")
    print("two pixels:")
    for x in rows_p:
        if x["case"] == "hom" or x["t_over_T2"] in (0.3, 1.0):
            print(f"  {x['case']:<4} t={x['t_over_T2']:5g} {x['preparation']:<8} {x['delta_p']:<7}: pair1 {x['pair_rate_pixel1']:.3f} pair2 {x['pair_rate_pixel2']:.3f} eff {x['efficiency']:.3f} share {x['share_pixel1']:.4f} ± {x['share_se']:.4f}")
    print(f"{time.time() - t0:.1f} s")
    return 0 if ok else 1

if __name__ == "__main__":
    sys.exit(main())
