#!/usr/bin/env python3
"""
empty_pixel.py: Part B (the empty pixel with the same contact and current) and Part C (the bulk's
phase distribution) of the 2026-09-08 addendum to PREREGISTRATION.md. Imports the definitions of
spad_contact.py; its own seed.
"""
import csv, sys, time
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, str(Path(__file__).resolve().parent))
import spad_contact as sc
from spad_contact import (X, DX, PIX, M, W, H, THETA, N_PRIMARY, RULES, SLOT, TEXT, TEXT2,
                          flight, region_of, contact, contact_column, breakdown, landscape, run_current)

OUT = sc.OUT
rng = np.random.default_rng(20260908 + 1)
sc.rng = rng                                   # the imported functions draw from this stream
R_LIST = [4.0, 2.0, 1.0]
THETAS = [(0.45, 500), (0.407, 450), (0.361, 400), (0.316, 350), (0.28, 310)]
TAUS = [0.0, 1.0, 2.0, 5.0, 10.0, 20.0, np.inf]
RULES_B = RULES + ["greedy_field"]                # D12b
BASE = {"greedy_field": ("greedy", 3.0)}
PHI0 = 2 * np.arccos(np.sqrt(THETA))

def pixel_norm(rho_d, pix):
    lo, hi = PIX[pix]
    return rho_d[(X >= lo) & (X <= hi)].sum() * DX

def sample_in_pixel(rho_d, pix, m):
    lo, hi = PIX[pix]
    sel = (X >= lo) & (X <= hi)
    cdf = np.cumsum(rho_d[sel]) * DX; cdf /= cdf[-1]
    return np.interp(rng.uniform(0, 1, m), cdf, X[sel])

# ----------------------------------------------------------------------------- Part B
def part_b(lines, rows):
    res_fixed = {1: landscape((H + 1, W)), 2: landscape((H + 1, W))}
    fig_data = {}
    for r in R_LIST:
        x0, ts, xs, rho_d = flight(r)
        x_arr = xs[-1]; reg = region_of(x_arr); inpix = reg > 0
        s = {1: pixel_norm(rho_d, 1), 2: pixel_norm(rho_d, 2)}
        ii = np.nonzero(inpix)[0]; m = ii.size
        reg_p = reg[ii]; reg_e = 3 - reg_p
        # the particle's pixel: contact as in the main run
        x_cp, rot_p, w_p = contact(x_arr[ii], reg_p, N_PRIMARY, "least")
        col_p = contact_column(x_cp, reg_p)
        # the empty pixel: D10
        x_e = np.empty(m)
        for pix in (1, 2):
            sel = reg_e == pix
            x_e[sel] = sample_in_pixel(rho_d, pix, sel.sum())
            # verification: the sampled contact points follow ρ restricted to the pixel
            lo, hi = PIX[pix]; selx = (X >= lo) & (X <= hi)
            cdf_rho = np.cumsum(rho_d[selx]) * DX; cdf_rho /= cdf_rho[-1]
            cdf_s = np.searchsorted(np.sort(x_e[sel]), X[selx], side="right") / sel.sum()
            ks = np.max(np.abs(cdf_rho - cdf_s))
            lines.append(f"r = {r:g}: empty-pixel contact points in pixel {pix} against ρ: max CDF difference {ks:.4f} (tolerance 0.03) -> {'PASS' if ks < 0.03 else 'FAIL'}")
        x_ce, rot_e, w_e = contact(x_e, reg_e, N_PRIMARY, "least")
        col_e = contact_column(x_ce, reg_e)
        s_e = np.where(reg_e == 1, s[1], s[2])
        bd_p = breakdown(m) >= sc.Q_TH; bd_e = breakdown(m) >= sc.Q_TH
        T = {}
        for rule in RULES_B:
            base, bias = BASE.get(rule, (rule, 0.0))
            Tp = np.empty(m); Te = np.empty(m); Cp = np.empty(m)
            for pix in (1, 2):
                sel = reg_p == pix; out = run_current(col_p[sel], base, True, res_fixed[pix], bias=bias); Tp[sel] = out[2]; Cp[sel] = out[0]
                sel = reg_e == pix; Te[sel] = run_current(col_e[sel], base, True, res_fixed[pix], bias=bias)[2]
            T[rule] = (Tp, Te)
            if r == 4.0 and rule in ("greedy", "greedy_field"):
                off = Cp - col_p
                lines.append(f"channel statistics, r = 4, the particle's contacts, {rule} (b = {bias:g}): offset sd {np.std(off):.2f} cells, resistance per row {Tp.mean() / H:.3f}, sd of T {Tp.std():.2f}, corr(contact, end) {np.corrcoef(col_p, Cp)[0, 1]:.3f}")
        for theta, nm in THETAS:
            for E in ("E1", "E2"):
                pair_e = (s_e >= theta) if E == "E1" else np.ones(m, bool)
                for Hc in ("H0", "H1"):
                    scale = (1.0 / s_e) if (E == "E1" and Hc == "H1") else np.ones(m)
                    for rule in RULES_B:
                        Tp, Te = T[rule]; Te_s = Te * scale
                        cand_p = bd_p; cand_e = pair_e & bd_e
                        for tau in TAUS:
                            both = cand_p & cand_e
                            close = both & (np.abs(Tp - Te_s) <= tau)
                            p_first = both & ~close & (Tp < Te_s); e_first = both & ~close & (Te_s < Tp)
                            rec_p = (cand_p & ~cand_e) | p_first | close
                            rec_e = (cand_e & ~cand_p) | e_first | close
                            n_rec = (rec_p | rec_e).sum()
                            rec1 = (rec_p & (reg_p == 1)).sum() + (rec_e & (reg_e == 1)).sum()
                            share = rec1 / (rec_p.sum() + rec_e.sum()) if n_rec else np.nan
                            rows.append(dict(r=r, theta=theta, nm=nm, E=E, H=Hc, rule=rule, tau=tau,
                                             efficiency=n_rec / M, share_pixel1=share, born=r / (1 + r),
                                             doubles=close.sum() / M, wrong_pixel=(rec_e & ~rec_p).sum() / M,
                                             both_candidates=both.sum() / M, empty_pairs=pair_e.mean(),
                                             empty_first_given_both=(e_first.sum() / max((p_first | e_first).sum(), 1))))
                            if r == 4.0 and theta == 0.45 and E == "E1" and Hc == "H0" and rule == "greedy" and tau == 0.0:
                                se = np.sqrt(0.8033 * 0.1967 / n_rec); ok = abs(share - 0.8027) < 3 * se and close.sum() == 0
                                lines.append(f"E1/H0, r = 4, θ = 0.45, greedy, τ_d = 0: share {share:.4f} against the main run's 0.8027 (s.e. {se:.4f}), doubles {close.sum()} -> {'PASS' if ok else 'FAIL'}")
                            if tau == 0.0 and close.sum() != 0:
                                lines.append(f"doubles at τ_d = 0: {close.sum()} -> FAIL")
                            if r == 4.0 and theta == 0.45 and E == "E2" and Hc == "H0" and rule == "greedy" and tau == np.inf:
                                bc = both.sum() / M; se = np.sqrt(0.515 * 0.485 / M); ok = abs(close.sum() / M - bc) < 1e-12 and abs(bc - 0.515) < 3 * se
                                lines.append(f"E2/H0, τ_d = ∞: doubles {close.sum() / M:.4f} = both candidates {bc:.4f} (predicted 0.515 ± {se:.4f}) -> {'PASS' if ok else 'FAIL'}")
        fig_data[r] = dict(s=s)
    return fig_data

# ----------------------------------------------------------------------------- Part C
def draw_phases(model, shape, idx, kappa=None, q=None):
    if model == "B0":
        return rng.uniform(-np.pi, np.pi, shape)
    if model == "B1":
        return rng.vonmises(0.0, kappa, shape)
    if model == "B2":
        return np.zeros(shape)
    if model == "B3":
        return (q * idx + np.pi) % (2 * np.pi) - np.pi

def part_c(lines, rows):
    x0, ts, xs, rho_d = flight(4.0)
    x_arr = xs[-1]; reg = region_of(x_arr); inpix = reg > 0
    xa = x_arr[inpix]; rg = reg[inpix]; m = xa.size; n = N_PRIMARY
    a = 1.0 / (n - 1)
    base_idx = np.round(xa / a).astype(int)
    idx = base_idx[:, None] + (np.arange(n) - n // 2)[None, :]
    atoms = idx * a
    lo = np.where(rg == 1, PIX[1][0], PIX[2][0])[:, None]; hi = np.where(rg == 1, PIX[1][1], PIX[2][1])[:, None]
    inside = (atoms >= lo - 1e-9) & (atoms <= hi + 1e-9)
    models = [("B0", "random, uniform iid", {}), ("B1", "ordered, von Mises κ = 1", dict(kappa=1.0)),
              ("B1", "ordered, von Mises κ = 4", dict(kappa=4.0)), ("B2", "locked crystal, one phase", {}),
              ("B3", "periodic crystal, q = 137.5°", dict(q=2 * np.pi / ((1 + np.sqrt(5)) / 2) ** 2)),
              ("B3", "periodic crystal, q = 120°", dict(q=2 * np.pi / 3))]
    for model, label, kw in models:
        dps = [0, 30, 60] if (model == "B3" and abs(kw["q"] - 2 * np.pi / 3) < 1e-9) else [0, 45, 90, 135, 180]
        for dp in dps:
            delta = draw_phases(model, atoms.shape, idx, **kw)
            mis = np.abs(np.angle(np.exp(1j * (np.radians(dp) - delta))))
            k = np.where(inside, mis, np.inf).argmin(axis=1)
            rot = mis[np.arange(m), k]; w = np.cos(rot / 2) ** 2
            miss = np.mean(w < THETA)
            rows.append(dict(model=model, label=label, delta_p_deg=dp, rotation_mean_deg=np.degrees(rot.mean()),
                             rotation_se_deg=np.degrees(rot.std() / np.sqrt(m)), miss_D4b=miss, efficiency_D4b=(1 - miss) * 0.85 * 0.778))
            if model == "B0" and dp == 0:
                pred = np.pi / (n + 1); se = rot.std() / np.sqrt(m); ok = abs(rot.mean() - pred) < 3 * se
                lines.append(f"B0 rotation {np.degrees(rot.mean()):.2f}° against π/22 = {np.degrees(pred):.2f}° (s.e. {np.degrees(se):.2f}°) -> {'PASS' if ok else 'FAIL'}")
            if model == "B2":
                err = np.max(np.abs(rot - np.radians(dp))); ok = err < 1e-12
                lines.append(f"B2 rotation at δ_p = {dp}°: max |rotation − δ_p| = {err:.1e} -> {'PASS' if ok else 'FAIL'}")
            if model == "B3" and dps == [0, 30, 60] and dp == 30:
                err = abs(np.degrees(rot.mean()) - 30); ok = err < 1e-6
                lines.append(f"B3 q = 120°, δ_p = 30°: rotation {np.degrees(rot.mean()):.8f}° -> {'PASS' if ok else 'FAIL'}")

# ----------------------------------------------------------------------------- main
def main():
    t0 = time.time(); lines = ["Verification, addendum 2026-09-08 (PREREGISTRATION.md, addendum)", ""]
    rows_b, rows_c = [], []
    fig_data = part_b(lines, rows_b)
    part_c(lines, rows_c)
    fmt = lambda v: (f"{v:.4f}" if isinstance(v, float) else v)
    with open(OUT / "empty_pixel.csv", "w", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=list(rows_b[0].keys())); wr.writeheader()
        for row in rows_b: wr.writerow({k: fmt(v) for k, v in row.items()})
    with open(OUT / "bulk_phase.csv", "w", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=list(rows_c[0].keys())); wr.writeheader()
        for row in rows_c: wr.writerow({k: fmt(v) for k, v in row.items()})
    # figure 4: doubles against the drain time (E2/H0, r = 4, θ = 0.45) and the share under E1 at θ = 0.28 for r = 1, 2, 4
    fig, axes = plt.subplots(1, 2, figsize=(9.4, 3.9))
    ax = axes[0]
    for rule, c in zip(RULES_B, SLOT):
        sel = [x for x in rows_b if x["r"] == 4.0 and x["theta"] == 0.45 and x["E"] == "E2" and x["H"] == "H0" and x["rule"] == rule]
        taus = [x["tau"] if np.isfinite(x["tau"]) else 50.0 for x in sel]
        ax.plot(taus, [x["doubles"] for x in sel], "o-", color=c, label=rule)
    ax.axhline(0.515, color=TEXT2, lw=0.8, ls="--"); ax.text(0.3, 0.53, "both pixels candidates, 0.515", color=TEXT2, fontsize=8)
    ax.set_xscale("symlog", linthresh=1); ax.set_xlabel("drain time τ_d (resistance units; 50 stands for ∞)"); ax.set_ylabel("double records per quantum")
    ax.set_title("E2 (empty lobe acts as if it held the particle), no head start"); ax.legend(loc="center right")
    ax = axes[1]
    xx = np.arange(3); wdt = 0.2
    for j, (E, Hc, lab, c) in enumerate((("E1", "H0", "E1, no head start", SLOT[0]), ("E1", "H1", "E1, channel speed ∝ deposit", SLOT[2]), ("E2", "H0", "E2", SLOT[1]))):
        vals = [next(x["share_pixel1"] for x in rows_b if x["r"] == r and x["theta"] == 0.28 and x["E"] == E and x["H"] == Hc and x["rule"] == "greedy" and x["tau"] == 0.0) for r in (1.0, 2.0, 4.0)]
        ax.bar(xx + (j - 1) * wdt, vals, wdt, color=c, label=lab)
    ax.plot(xx, [0.5, 2 / 3, 0.8], "k_", ms=22, mew=1.6, label="Born r/(1 + r)")
    ax.set_xticks(xx); ax.set_xticklabels(["r = 1", "r = 2", "r = 4"]); ax.set_ylim(0, 1); ax.set_ylabel("share of records at pixel 1")
    ax.set_title("θ = 0.28 (310 nm), greedy channel, τ_d = 0"); ax.legend(loc="upper left", fontsize=7.5)
    fig.tight_layout(); fig.savefig(OUT / "fig4_empty_pixel.png", dpi=150); plt.close(fig)
    # figure 5: the bulk models
    fig, axes = plt.subplots(1, 2, figsize=(9.4, 3.9))
    labels = []
    for row in rows_c:
        if row["label"] not in labels: labels.append(row["label"])
    for lab, c in zip(labels, SLOT):
        sel = [x for x in rows_c if x["label"] == lab]
        axes[0].plot([x["delta_p_deg"] for x in sel], [max(x["rotation_mean_deg"], 0.05) for x in sel], "o-", color=c, label=lab)
        axes[1].plot([x["delta_p_deg"] for x in sel], [1 - x["miss_D4b"] for x in sel], "o-", color=c, label=lab)
    axes[0].set_yscale("log"); axes[0].set_xlabel("particle's offset δ_p (degrees)"); axes[0].set_ylabel("rotation used (degrees)")
    axes[0].set_title("Least rotation, 21 atoms in reach, by bulk model"); axes[0].legend(fontsize=7)
    axes[1].set_xlabel("particle's offset δ_p (degrees)"); axes[1].set_ylabel("fraction making a pair (D4b, θ = 0.45)"); axes[1].set_ylim(-0.02, 1.05)
    axes[1].set_title("The pair rate reads the particle's phase only in an ordered or locked bulk")
    fig.tight_layout(); fig.savefig(OUT / "fig5_bulk_phase.png", dpi=150); plt.close(fig)
    text = "\n".join(lines) + "\n"; (OUT / "verification_addendum.txt").write_text(text); print(text)
    ok = all("FAIL" not in l for l in lines)
    with open(OUT / "run_info.txt", "a") as f:
        f.write(f"empty_pixel.py, run {time.strftime('%Y-%m-%d %H:%M:%S')}, seed 20260909, {time.time() - t0:.1f} s\nverification (addendum): {'PASS' if ok else 'FAIL'}\n")
    # console summary
    print("Part B, greedy, τ_d = 0 and 5, selected rows:")
    for x in rows_b:
        if x["rule"] == "greedy" and x["tau"] in (0.0, 5.0) and (x["theta"] in (0.45, 0.407, 0.28)):
            print(f"  r={x['r']:g} θ={x['theta']:.3f} {x['E']} {x['H']} τ={x['tau']:g}: eff {x['efficiency']:.3f} share {x['share_pixel1']:.3f} (Born {x['born']:.3f}) doubles {x['doubles']:.3f} wrong {x['wrong_pixel']:.3f} both {x['both_candidates']:.3f} emptyfirst|both {x['empty_first_given_both']:.3f}")
    print("Part B, E2/H0, r = 4, θ = 0.45, doubles by rule and τ_d:")
    for rule in RULES_B:
        print("  " + rule + ": " + ", ".join(f"{x['tau']:g}: {x['doubles']:.3f}" for x in rows_b if x["r"] == 4.0 and x["theta"] == 0.45 and x["E"] == "E2" and x["H"] == "H0" and x["rule"] == rule))
    print("Part C:")
    for x in rows_c:
        print(f"  {x['label']:<30} δ_p={x['delta_p_deg']:>3}: rotation {x['rotation_mean_deg']:.2f}° miss {x['miss_D4b']:.4f}")
    print(f"{time.time() - t0:.1f} s")
    return 0 if ok else 1

if __name__ == "__main__":
    sys.exit(main())
