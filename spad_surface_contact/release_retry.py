#!/usr/bin/env python3
"""
release_retry.py: Part G (addendum 5 of PREREGISTRATION.md). The excited electron at each pixel with
three release channels (ejection above the gap, re-emission back to the wave, phonon loss), the
quantum's energy tracked to closure; reading A (shares as energies) and reading B (the indivisible
quantum). A continuous-time Markov chain per quantum.
"""
import csv, sys, time
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent; OUT = HERE / "results"
EMISSION = "escape" if (len(sys.argv) > 1 and sys.argv[1] == "escape") else "retry"     # addendum 5b
rng = np.random.default_rng(20260914 if EMISSION == "retry" else 20260915)
SHARES = {1.0: (0.4254, 0.4254), 2.0: (0.5699, 0.2819), 4.0: (0.6865, 0.1680)}
THETAS = [0.6, 0.45, 0.36, 0.28]; NQ = 3000
PRIMARY = dict(g_dep=10.0, g_em=1.0, g_ej=1.0, g_ph=0.3, g_esc=0.0)
SENS = [dict(g_ej=0.3), dict(g_ej=3.0), dict(g_ph=0.0), dict(g_ph=1.0), dict(g_em=0.3), dict(g_em=3.0), dict(g_esc=0.3)]
SURFACE, TEXT, TEXT2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e6e5e1"
SLOT = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#7a5cc7"]
plt.rcParams.update({"font.family": "sans-serif", "font.size": 9, "axes.titlesize": 10, "axes.titleweight": "bold",
    "axes.edgecolor": GRID, "axes.labelcolor": TEXT2, "xtick.color": TEXT2, "ytick.color": TEXT2, "axes.titlecolor": TEXT,
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE, "axes.grid": True, "grid.color": GRID,
    "grid.linewidth": 0.8, "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False, "legend.fontsize": 8, "axes.axisbelow": True})

def quantum(reading, s, theta, g_dep, g_em, g_ej, g_ph, g_esc, tmax=200.0):
    E = [0.0, 0.0]; Ew = 1.0; heat = esc = 0.0; rec = [0.0, 0.0]; nrec = [0, 0]; t = 0.0; cycles = 0; first = np.nan
    while (E[0] > 0 or E[1] > 0 or Ew > 0) and t < tmax:
        ev, rates = [], []
        if Ew > 0:
            ev.append(("dep", None)); rates.append(g_dep)
            if g_esc > 0: ev.append(("esc", None)); rates.append(g_esc)
        for k in (0, 1):
            if E[k] > 0:
                ev.append(("em", k)); rates.append(g_em); ev.append(("ph", k)); rates.append(g_ph)
                if E[k] >= theta - 1e-12: ev.append(("ej", k)); rates.append(g_ej)
        R = sum(rates); t += rng.exponential(1 / R)
        kind, k = ev[rng.choice(len(ev), p=np.array(rates) / R)]
        if kind == "dep":
            if reading == "A":
                E[0] += Ew * s[0]; E[1] += Ew * s[1]; esc += Ew * (1 - s[0] - s[1])
            else:
                u = rng.uniform()
                if u < s[0]: E[0] += Ew
                elif u < s[0] + s[1]: E[1] += Ew
                else: esc += Ew
            Ew = 0.0
        elif kind == "esc": esc += Ew; Ew = 0.0
        elif kind == "em":
            if EMISSION == "retry": Ew += E[k]
            else: esc += E[k]
            E[k] = 0.0; cycles += 1
        elif kind == "ph": heat += E[k]; E[k] = 0.0
        else:
            rec[k] += E[k]; nrec[k] += 1
            if np.isnan(first): first = E[k]
            E[k] = 0.0
    closure = rec[0] + rec[1] + heat + esc + E[0] + E[1] + Ew
    return nrec[0] > 0, nrec[1] > 0, first, closure, cycles, heat, esc

def config(reading, r, theta, **kw):
    p = dict(PRIMARY); p.update(kw); s = SHARES[r]
    o = np.array([quantum(reading, s, theta, **p) for _ in range(NQ)], dtype=float)
    r1, r2 = o[:, 0], o[:, 1]; tot = r1.sum() + r2.sum()
    return dict(reading=reading, r=r, theta=theta, **p, share_pixel1=(r1.sum() / tot) if tot else np.nan, born=SHARES[r][0] / sum(SHARES[r]),
                doubles=np.mean(r1 * r2), efficiency=np.mean((r1 + r2) > 0), record_energy_mean=np.nanmean(o[:, 2]),
                record_energy_min=np.nanmin(o[:, 2]), cycles_mean=o[:, 4].mean(), heat_mean=o[:, 5].mean(), escape_mean=o[:, 6].mean(),
                closure_max_error=np.abs(o[:, 3] - 1).max(), n=NQ), o[:, 2]

def main_escape():
    t0 = time.time(); lines = ["Verification, Part G, addendum 5b (re-emitted photons leave)", ""]; rows = []
    pred = {("A", 1.0, 0.6): (np.nan, 0.0, 0.0), ("A", 1.0, 0.45): (np.nan, 0.0, 0.0), ("A", 1.0, 0.36): (0.5, 0.68, 0.19), ("A", 1.0, 0.28): (0.5, 0.68, 0.19),
            ("A", 2.0, 0.6): (1.0, 0.435, 0.0), ("A", 2.0, 0.45): (1.0, 0.435, 0.0), ("A", 2.0, 0.36): (1.0, 0.435, 0.0), ("A", 2.0, 0.28): (0.5, 0.68, 0.19),
            ("A", 4.0, 0.6): (1.0, 0.435, 0.0), ("A", 4.0, 0.45): (1.0, 0.435, 0.0), ("A", 4.0, 0.36): (1.0, 0.435, 0.0), ("A", 4.0, 0.28): (1.0, 0.435, 0.0)}
    for reading in ("A", "B"):
        for r in SHARES:
            for th in THETAS:
                d, _ = config(reading, r, th); rows.append(d)
                if reading == "A":
                    ps, pe, pd = pred[(reading, r, th)]
                    ok = d["closure_max_error"] < 1e-12 and abs(d["efficiency"] - pe) < 0.03 and abs(d["doubles"] - pd) < 0.025 and (np.isnan(ps) or abs(d["share_pixel1"] - ps) < 0.03)
                    lines.append(f"A r = {r:g} θ = {th}: share {d['share_pixel1']:.3f} (P20 {ps}), efficiency {d['efficiency']:.3f} ({pe}), doubles {d['doubles']:.3f} ({pd}), closure {d['closure_max_error']:.1e} -> {'PASS' if ok else 'FAIL'}")
                else:
                    se = np.sqrt(d["born"] * (1 - d["born"]) / max(d["efficiency"] * NQ, 1))
                    ok = d["closure_max_error"] < 1e-12 and d["doubles"] == 0.0 and abs(d["record_energy_mean"] - 1) < 1e-12 and abs(d["efficiency"] - 0.372) < 0.03 and abs(d["share_pixel1"] - d["born"]) < 3 * se
                    lines.append(f"B r = {r:g} θ = {th}: share {d['share_pixel1']:.3f} (Born {d['born']:.3f}), efficiency {d['efficiency']:.3f} (0.372), doubles {d['doubles']:.3f}, record energy {d['record_energy_mean']:.4f} -> {'PASS' if ok else 'FAIL'}")
    keys = sorted({k for r in rows for k in r}, key=lambda k: (k != "reading", k))
    with open(OUT / "release_escape.csv", "w", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=keys); wr.writeheader()
        for r in rows: wr.writerow({k: (f"{v:.5g}" if isinstance(v, float) else v) for k, v in r.items()})
    text = "\n".join(lines) + "\n"; (OUT / "verification_release_escape.txt").write_text(text); print(text)
    ok = all("FAIL" not in l for l in lines)
    with open(OUT / "run_info.txt", "a") as f:
        f.write(f"release_retry.py escape, run {time.strftime('%Y-%m-%d %H:%M:%S')}, seed 20260915, {NQ} quanta per configuration, {time.time() - t0:.1f} s\nverification (Part G, 5b): {'PASS' if ok else 'FAIL'}\n")
    return 0 if ok else 1

def main():
    t0 = time.time(); lines = ["Verification, Part G (PREREGISTRATION.md, addendum 5)", ""]; rows = []; energies = {}
    pred_A = {(4.0, 0.45): 0.996, (1.0, 0.45): 0.48, (1.0, 0.6): 0.50}
    for reading in ("A", "B"):
        for r in SHARES:
            for th in THETAS:
                d, en = config(reading, r, th); rows.append(d); energies[(reading, r, th)] = en
                ok = d["closure_max_error"] < 1e-12
                msg = f"{reading} r = {r:g} θ = {th}: closure {d['closure_max_error']:.1e} -> {'PASS' if ok else 'FAIL'}"
                if th == 0.6 or reading == "B":
                    ok2 = d["doubles"] == 0.0; msg += f"; doubles {d['doubles']:.4f} (0 exactly) -> {'PASS' if ok2 else 'FAIL'}"
                if reading == "B":
                    ok3 = abs(d["record_energy_mean"] - 1) < 1e-12 and abs(d["record_energy_min"] - 1) < 1e-12; msg += f"; record energy {d['record_energy_mean']:.6f} (1 exactly) -> {'PASS' if ok3 else 'FAIL'}"
                    se = np.sqrt(d["born"] * (1 - d["born"]) / (d["efficiency"] * NQ)); ok4 = abs(d["share_pixel1"] - d["born"]) < 3 * se
                    msg += f"; share {d['share_pixel1']:.4f} against Born {d['born']:.4f} (s.e. {se:.4f}) -> {'PASS' if ok4 else 'FAIL'}"
                if reading == "A" and (r, th) in pred_A:
                    pv = pred_A[(r, th)]; se = np.sqrt(max(pv * (1 - pv), 1e-4) / max(d["efficiency"] * NQ, 1)); ok5 = abs(d["share_pixel1"] - pv) < 3 * se + 0.005
                    msg += f"; share {d['share_pixel1']:.4f} against P19 {pv} (s.e. {se:.4f}) -> {'PASS' if ok5 else 'FAIL'}"
                lines.append(msg)
    for kw in SENS:
        d, _ = config("A", 4.0, 0.45, **kw); d["sensitivity"] = str(kw); rows.append(d)
        lines.append(f"A sensitivity {kw}: closure {d['closure_max_error']:.1e} -> {'PASS' if d['closure_max_error'] < 1e-12 else 'FAIL'}; share {d['share_pixel1']:.3f}, eff {d['efficiency']:.3f}, E_rec {d['record_energy_mean']:.3f}, cycles {d['cycles_mean']:.1f}")
    keys = sorted({k for r in rows for k in r}, key=lambda k: (k != "reading", k))
    with open(OUT / "release_retry.csv", "w", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=keys); wr.writeheader()
        for r in rows: wr.writerow({k: (f"{v:.5g}" if isinstance(v, float) else v) for k, v in r.items()})
    # figure
    fig, axes = plt.subplots(1, 3, figsize=(13, 3.9))
    ax = axes[0]; rs = [1.0, 2.0, 4.0]
    for th, c in zip(THETAS, SLOT):
        ax.plot(rs, [next(x["share_pixel1"] for x in rows if x["reading"] == "A" and x["r"] == r and x["theta"] == th and "sensitivity" not in x) for r in rs], "o-", color=c, label=f"reading A, θ = {th}")
    ax.plot(rs, [next(x["share_pixel1"] for x in rows if x["reading"] == "B" and x["r"] == r and x["theta"] == 0.45 and "sensitivity" not in x) for r in rs], "s--", color=TEXT, label="reading B (indivisible), any θ")
    ax.plot(rs, [SHARES[r][0] / sum(SHARES[r]) for r in rs], "k_", ms=22, mew=1.6, label="Born r/(1 + r)")
    ax.set_xscale("log"); ax.set_xticks(rs); ax.set_xticklabels(["1", "2", "4"]); ax.set_xlabel("intensity ratio r"); ax.set_ylabel("share of records at pixel 1"); ax.set_ylim(0.4, 1.05)
    ax.set_title("Divisible energy makes the share a step; indivisible gives Born"); ax.legend(fontsize=7, loc="lower right")
    ax = axes[1]
    for r, c in zip(rs, SLOT):
        sel = [x for x in rows if x["reading"] == "A" and x["r"] == r and "sensitivity" not in x]
        ax.plot([x["theta"] for x in sel], [x["doubles"] for x in sel], "o-", color=c, label=f"doubles, r = {r:g}")
        ax.plot([x["theta"] for x in sel], [x["efficiency"] for x in sel], "s:", color=c, label=f"efficiency, r = {r:g}")
    ax.axvline(0.5, color=TEXT2, lw=0.8, ls="--"); ax.text(0.505, 0.9, "θ = ½: doubles impossible above", fontsize=7.5, color=TEXT2)
    ax.set_xlabel("gap ratio θ = E_gap / E_photon"); ax.set_ylabel("per quantum"); ax.set_title("Reading A: doubles and efficiency against the gap ratio"); ax.legend(fontsize=7, loc="center right")
    ax = axes[2]
    en = energies[("A", 4.0, 0.45)]; en = en[~np.isnan(en)]
    ax.hist(en, bins=np.linspace(0.4, 1.0, 31), color=SLOT[0], alpha=0.8, label="reading A, r = 4, θ = 0.45")
    ax.axvline(1.0, color=TEXT, lw=2.0, label="reading B: every record carries the whole quantum")
    ax.axvline(0.45, color=TEXT2, lw=0.8, ls="--"); ax.text(0.455, ax.get_ylim()[1] * 0.9, "gap", fontsize=8, color=TEXT2)
    ax.set_xlabel("energy of the first record (quanta)"); ax.set_ylabel("quanta"); ax.set_title("What a record carries (calorimetry measures 1)"); ax.legend(fontsize=7.5, loc="upper left")
    fig.tight_layout(); fig.savefig(OUT / "fig10_release.png", dpi=150); plt.close(fig)
    text = "\n".join(lines) + "\n"; (OUT / "verification_release.txt").write_text(text); print(text)
    for x in rows:
        if "sensitivity" not in x:
            print(f"  {x['reading']} r={x['r']:g} θ={x['theta']}: share {x['share_pixel1']:.3f} (Born {x['born']:.3f}) doubles {x['doubles']:.3f} eff {x['efficiency']:.3f} E_rec {x['record_energy_mean']:.3f} cycles {x['cycles_mean']:.1f} heat {x['heat_mean']:.3f} esc {x['escape_mean']:.3f}")
    ok = all("FAIL" not in l for l in lines)
    with open(OUT / "run_info.txt", "a") as f:
        f.write(f"release_retry.py, run {time.strftime('%Y-%m-%d %H:%M:%S')}, seed 20260914, {NQ} quanta per configuration, {time.time() - t0:.1f} s\nverification (Part G): {'PASS' if ok else 'FAIL'}\n")
    print(f"{time.time() - t0:.1f} s")
    return 0 if ok else 1

if __name__ == "__main__":
    sys.exit(main_escape() if EMISSION == "escape" else main())
