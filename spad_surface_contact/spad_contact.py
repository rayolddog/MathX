#!/usr/bin/env python3
"""
spad_contact.py: two SPAD pixels, surface contact by least rotation, the record as an induced current
=====================================================================================================

PREREGISTRATION.md of this folder (2026-09-08). The packet, the exact free Dirac propagator and the
Bohm–Dirac trajectories are the second study's (../second_spinor_phase_race/bohm_dirac_flight.py,
ħ = c = m = 1, H = σ_x p + σ_z m, standard Dirac representation). Everything else is defined here.

  flight     D2   two lobes r, σ = 1 at x = ±4, exact propagation to t_d = 4, 4000 RK4 trajectories
  pixels     D3   [1, 7] and [−7, −1]; the gap and |x| > 7 hold no atoms
  reach      D4   n surface atoms on a lattice of spacing 1/(n − 1) centred on the nearest point, ±0.5
  contact    D5   least rotation (argmin |Δφ|) primary; R2 (sampled ∝ cos²(Δφ/2)) comparator
  deposit    D6   D4a whole quantum; D4b the share w; pair iff deposit ≥ θ = 0.45
  current    D7   121 × 60 lattice of lognormal resistances; greedy / DBM η = 1 / diffusion transfer
  breakdown  D8   Galton–Watson, loss 0.1 and ionization 0.5 per row, record iff ≥ 100 at row 60
"""
import csv, sys, time
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "second_spinor_phase_race"))
from bohm_dirac_flight import X, DX, gauss, normalise, to_p, psi_at, rho_v, integrate, SURFACE, TEXT, TEXT2, GRID

OUT = HERE / "results"
OUT.mkdir(exist_ok=True)
rng = np.random.default_rng(20260908)

# ----------------------------------------------------------------------------- fixed definitions
R_LIST = [4.0, 1.0]
X1, SIGMA, T_D, DT, M = 4.0, 1.0, 4.0, 0.01, 4000
PIX = {1: (1.0, 7.0), 2: (-7.0, -1.0)}
N_LIST = [1, 3, 11, 21, 101, 1001]
N_PRIMARY = 21
THETA = 0.45
A_COL, W, H = 0.05, 121, 60
P_LOSS, P_ION, Q_TH = 0.1, 0.5, 100
RULES = ["greedy", "dbm", "diffusion"]

SLOT = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#7a5cc7", "#c2417a"]
plt.rcParams.update({
    "font.family": "sans-serif", "font.size": 9, "axes.titlesize": 10, "axes.titleweight": "bold",
    "axes.edgecolor": GRID, "axes.labelcolor": TEXT2, "xtick.color": TEXT2, "ytick.color": TEXT2,
    "axes.titlecolor": TEXT, "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
    "savefig.facecolor": SURFACE, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.8,
    "grid.linestyle": "-", "axes.spines.top": False, "axes.spines.right": False,
    "legend.frameon": False, "legend.fontsize": 8, "axes.axisbelow": True,
})

# ----------------------------------------------------------------------------- D2: flight
def flight(r):
    g = np.sqrt(r) * gauss(X, X1, SIGMA) + gauss(X, -X1, SIGMA)
    psi0 = normalise(np.array([g, np.zeros_like(g)], dtype=complex))
    psi0_p = to_p(psi0)
    rho0, _ = rho_v(psi0)
    cdf = np.cumsum(rho0) * DX; cdf /= cdf[-1]
    x0 = np.interp(rng.uniform(0, 1, M), cdf, X)
    ts, xs = integrate(psi0_p, x0, T_D, DT)
    rho_d, _ = rho_v(psi_at(psi0_p, T_D))
    return x0, ts, xs, rho_d

def region_of(x):
    """1, 2 for the pixels; 0 for the gap; -1 for outside."""
    reg = np.full(x.shape, -1, dtype=int)
    reg[(x >= PIX[1][0]) & (x <= PIX[1][1])] = 1
    reg[(x >= PIX[2][0]) & (x <= PIX[2][1])] = 2
    reg[np.abs(x) < 1.0] = 0
    return reg

def rho_regions(rho_d):
    inpix = lambda lo, hi: rho_d[(X >= lo) & (X <= hi)].sum() * DX
    p1, p2 = inpix(*PIX[1]), inpix(*PIX[2])
    gap = rho_d[np.abs(X) < 1.0].sum() * DX
    return dict(p1=p1, p2=p2, gap=gap, out=1 - p1 - p2 - gap)

# ----------------------------------------------------------------------------- D4–D6: contact
def contact(x_arr, reg, n, rule):
    """Atoms in reach, the phase draw, the winner. Returns x_contact, rotation |Δφ|, w."""
    a = 1.0 if n == 1 else 1.0 / (n - 1)
    base = np.round(x_arr / a) * a
    atoms = base[:, None] + ((np.arange(n) - n // 2) * a)[None, :]
    lo = np.where(reg == 1, PIX[1][0], PIX[2][0])[:, None]; hi = np.where(reg == 1, PIX[1][1], PIX[2][1])[:, None]
    inside = (atoms >= lo - 1e-9) & (atoms <= hi + 1e-9)
    dphi = rng.uniform(-np.pi, np.pi, atoms.shape)
    w = np.where(inside, np.cos(dphi / 2) ** 2, 0.0)
    idx = np.arange(x_arr.size)
    if rule == "least":
        k = np.where(inside, np.abs(dphi), np.inf).argmin(axis=1)
    else:
        cdf = np.cumsum(w, axis=1); u = rng.uniform(0, 1, x_arr.size) * cdf[:, -1]
        k = np.minimum((cdf < u[:, None]).sum(axis=1), n - 1)
    return atoms[idx, k], np.abs(dphi[idx, k]), w[idx, k]

def contact_column(x_c, reg):
    origin = np.where(reg == 1, PIX[1][0], PIX[2][0])
    return np.clip(np.round((x_c - origin) / A_COL).astype(int), 0, W - 1)

# ----------------------------------------------------------------------------- D8: breakdown
def breakdown(m):
    c = np.ones(m, dtype=np.int64); cap = 10 ** 7
    for _ in range(H):
        live = c < cap
        surv = np.where(live, rng.binomial(np.minimum(c, cap), 1 - P_LOSS), c)
        ion = np.where(live, rng.binomial(np.minimum(surv, cap), P_ION), 0)
        c = np.where(live, surv + ion, c)
    return c

def extinction_closed_form():
    a, b, cc = (1 - P_LOSS) * P_ION, -(1 - (1 - P_LOSS) * (1 - P_ION)), P_LOSS
    return (-b - np.sqrt(b * b - 4 * a * cc)) / (2 * a)

# ----------------------------------------------------------------------------- D7: induced current
def landscape(shape):
    return np.exp(rng.standard_normal(shape))

def transfer_step(n, r_next, rule):
    """One row of charge transfer. n: (m, W). r_next: (W,) or (m, W) resistances of the next row."""
    if r_next.ndim == 1:
        r_next = np.broadcast_to(r_next, n.shape)
    g = 1.0 / r_next
    gm = np.zeros_like(n); gm[:, 1:] = g[:, :-1]      # destination c − 1 seen from source c
    gp = np.zeros_like(n); gp[:, :-1] = g[:, 1:]      # destination c + 1
    g0 = g
    if rule == "diffusion":
        gm, g0, gp = (gm > 0).astype(float), np.ones_like(n), (gp > 0).astype(float)
    tot = gm + g0 + gp
    wm, w0, wp = gm / tot, g0 / tot, gp / tot
    new = n * w0
    new[:, :-1] += (n * wm)[:, 1:]
    new[:, 1:] += (n * wp)[:, :-1]
    return new

def run_current(col0, rule, fixed, res_fixed=None, bias=0.0):
    """Charge from (row 0, col0) to the electrode. Returns centroid, width, time, end column (greedy)."""
    m = col0.size; idx = np.arange(m); cols = np.arange(W)
    if rule == "greedy":
        col = col0.copy(); T = np.zeros(m)
        for k in range(H):
            r_next = res_fixed[k + 1] if fixed else landscape((m, W))
            cands = np.stack([col - 1, col, col + 1], axis=1)
            valid = (cands >= 0) & (cands < W)
            cc = np.clip(cands, 0, W - 1)
            rr = r_next[cc] if r_next.ndim == 1 else r_next[idx[:, None], cc]
            rr = np.where(valid, rr, np.inf)
            j = (rr * np.array([1 + bias, 1.0, 1 + bias])).argmin(axis=1)     # D12b: lateral penalty for the choice only
            col = cc[idx, j]; T += rr[idx, j]
        return col.astype(float), np.zeros(m), T, col
    n = np.zeros((m, W)); n[idx, col0] = 1.0; T = np.zeros(m)
    for k in range(H):
        r_next = res_fixed[k + 1] if fixed else landscape((m, W))
        n = transfer_step(n, r_next, rule)
        T += (n * r_next).sum(axis=1)
    mass = n.sum(axis=1)
    cen = (n * cols).sum(axis=1) / mass
    wid = np.sqrt((n * (cols[None, :] - cen[:, None]) ** 2).sum(axis=1) / mass)
    return cen, wid, T, mass

# ----------------------------------------------------------------------------- main
def main():
    t0 = time.time()
    lines = ["Verification (PREREGISTRATION.md §5)", ""]
    q_pred = extinction_closed_form()
    summary, field_rows, current_rows = [], [], []
    trials_primary = None
    fig1_data = {}
    res_fixed = {1: landscape((H + 1, W)), 2: landscape((H + 1, W))}      # A8: one device, one landscape per pixel

    for r in R_LIST:
        x0, ts, xs, rho_d = flight(r)
        x_arr = xs[-1]; reg = region_of(x_arr)
        rr = rho_regions(rho_d)
        fr = {k: np.mean(reg == v) for k, v in (("p1", 1), ("p2", 2), ("gap", 0), ("out", -1))}
        ok_arr = True
        for k in ("p1", "p2", "gap", "out"):
            se = np.sqrt(rr[k] * (1 - rr[k]) / M); ok = abs(fr[k] - rr[k]) < 3 * se; ok_arr &= ok
            lines.append(f"r = {r:g}: arrivals {k} {fr[k]:.4f} against ∫ρ {rr[k]:.4f} (s.e. {se:.4f}) -> {'PASS' if ok else 'FAIL'}")
        share_rho = rr["p1"] / (rr["p1"] + rr["p2"])
        fig1_data[r] = (xs, x_arr, rho_d, reg)
        n_list = N_LIST if r == 4.0 else [N_PRIMARY]
        for n in n_list:
            for rule in ("least", "R2"):
                inpix = reg > 0
                x_c = np.full(M, np.nan); rot = np.full(M, np.nan); w = np.full(M, np.nan)
                x_c[inpix], rot[inpix], w[inpix] = contact(x_arr[inpix], reg[inpix], n, rule)
                # verification of the rotation and the D4b miss (least rotation only, closed forms)
                if rule == "least":
                    pred = np.pi / (n + 1); se = rot[inpix].std() / np.sqrt(inpix.sum())
                    ok = abs(rot[inpix].mean() - pred) < 3 * se
                    lines.append(f"r = {r:g}, n = {n}: least-rotation mean |Δφ| {np.degrees(rot[inpix].mean()):.2f}° against π/(n+1) = {np.degrees(pred):.2f}° (s.e. {np.degrees(se):.2f}°) -> {'PASS' if ok else 'FAIL'}")
                    edge = np.minimum(np.abs(x_arr - np.where(reg == 1, PIX[1][0], PIX[2][0])), np.abs(x_arr - np.where(reg == 1, PIX[1][1], PIX[2][1])))
                    interior = inpix & (edge > 0.5)
                    lines.append(f"    diagnostic, not registered (added 2026-09-08 after the run): arrivals more than 0.5 from a pixel edge ({interior.sum()}): mean |Δφ| {np.degrees(rot[interior].mean()):.2f}° (s.e. {np.degrees(rot[interior].std() / np.sqrt(interior.sum())):.2f}°)")
                    if n <= 3:
                        phi0 = 2 * np.arccos(np.sqrt(THETA)); pm = (1 - phi0 / np.pi) ** n
                        miss = np.mean(w[inpix] < THETA); se = np.sqrt(pm * (1 - pm) / inpix.sum())
                        ok = abs(miss - pm) < 3 * se
                        lines.append(f"r = {r:g}, n = {n}: D4b miss {miss:.4f} against 0.468^n = {pm:.4f} (s.e. {se:.4f}) -> {'PASS' if ok else 'FAIL'}")
                bd = breakdown(M)
                if r == 4.0 and n == N_PRIMARY and rule == "least":
                    ext = np.mean(bd == 0); se = np.sqrt(q_pred * (1 - q_pred) / M); ok = abs(ext - q_pred) < 3 * se
                    lines.append(f"Galton–Watson extinction {ext:.4f} against closed form {q_pred:.4f} (s.e. {se:.4f}) -> {'PASS' if ok else 'FAIL'}")
                for dep in ("D4a", "D4b"):
                    pair = inpix & ((w >= THETA) if dep == "D4b" else np.ones(M, bool))
                    record = pair & (bd >= Q_TH)
                    n_rec = record.sum()
                    share = np.mean(reg[record] == 1) if n_rec else np.nan
                    se_share = np.sqrt(share * (1 - share) / n_rec) if n_rec else np.nan
                    summary.append(dict(r=r, n=n, rule=rule, deposit=dep, quanta=M, in_pixel=int(inpix.sum()),
                                        pairs=int(pair.sum()), records=int(n_rec), efficiency=n_rec / M,
                                        miss_rate=1 - pair.sum() / inpix.sum(), P_b=np.mean(bd[pair] >= Q_TH) if pair.sum() else np.nan,
                                        share_pixel1=share, share_se=se_share, share_rho=share_rho, born=r / (1 + r),
                                        doubles=0, rotation_mean_deg=np.degrees(np.nanmean(rot)),
                                        rotation_record_mean_deg=np.degrees(rot[record].mean()) if n_rec else np.nan,
                                        w_mean=np.nanmean(w)))
                    if r == 4.0 and n == N_PRIMARY and rule == "least" and dep == "D4a":
                        # the induced current for the primary configuration
                        col0 = contact_column(x_c[inpix], reg[inpix])
                        cur = {}
                        for rule_c in RULES:
                            for fixed in (True, False):
                                cen = np.full(M, np.nan); wid = np.full(M, np.nan); T = np.full(M, np.nan); mass = np.full(M, np.nan)
                                for pix in (1, 2):
                                    sel = reg[inpix] == pix
                                    c_, w_, t_, m_ = run_current(col0[sel], rule_c, fixed, res_fixed[pix])
                                    ii = np.nonzero(inpix)[0][sel]
                                    cen[ii], wid[ii], T[ii], mass[ii] = c_, w_, t_, m_
                                cur[(rule_c, fixed)] = (cen, wid, T, mass)
                                ii = np.nonzero(inpix)[0]
                                off = cen[ii] - col0
                                rec_ii = record[ii]
                                distinct = len(np.unique(cen[ii][(reg[ii] == 1)].round())) if rule_c == "greedy" else np.nan
                                current_rows.append(dict(rule=rule_c, landscape="fixed" if fixed else "fresh",
                                                         width_mean=np.nanmean(wid[ii]), offset_sd=np.std(off), offset_mean=np.mean(off),
                                                         corr_contact_centroid=np.corrcoef(col0, cen[ii])[0, 1],
                                                         time_per_row=np.mean(T[ii]) / H, time_sd=np.std(T[ii]),
                                                         distinct_end_columns_pixel1=distinct,
                                                         charge_conservation_max_err=np.nanmax(np.abs(mass[ii] - 1)) if rule_c != "greedy" else 0.0))
                                if rule_c == "diffusion" and fixed:
                                    pred = np.sqrt(2 * H / 3); got = np.nanmean(wid[ii]); ok = abs(got - pred) / pred < 0.05
                                    lines.append(f"diffusion width {got:.3f} against √(2H/3) = {pred:.3f} -> {'PASS' if ok else 'FAIL'}")
                                    inner = (col0 >= 10) & (col0 <= W - 11)
                                    lines.append(f"    diagnostic, not registered (added 2026-09-08 after the run): diffusion offset sd, contacts 10 or more columns from a lattice edge ({inner.sum()}): {np.std(off[inner]):.3f}; within 10 columns of an edge ({(~inner).sum()}): {np.std(off[~inner]):.3f}")
                                if rule_c != "greedy":
                                    err = np.nanmax(np.abs(mass[ii] - 1)); ok = err < 1e-10
                                    lines.append(f"{rule_c} ({'fixed' if fixed else 'fresh'}) charge conservation: max |mass − 1| = {err:.1e} -> {'PASS' if ok else 'FAIL'}")
                        trials_primary = dict(x0=x0, x_arr=x_arr, reg=reg, x_c=x_c, rot=rot, w=w, bd=bd, record=record,
                                              col0=np.full(M, -1, int), cur=cur)
                        trials_primary["col0"][inpix] = col0
        # C1: the field reading, every atom of both pixels in play
        if r == 4.0:
            for n_s in (121, 5):
                xa = np.concatenate([np.linspace(*PIX[1], n_s), np.linspace(*PIX[2], n_s)])
                pix_a = np.r_[np.ones(n_s, int), np.full(n_s, 2)]
                A = np.interp(xa, X, rho_d)
                res = {"pure_least": [], "R3_deposit": [], "R2": []}
                for _ in range(M):
                    dphi = rng.uniform(-np.pi, np.pi, 2 * n_s); wv = np.cos(dphi / 2) ** 2; E = A * wv
                    res["pure_least"].append(pix_a[np.abs(dphi).argmin()])
                    res["R3_deposit"].append(pix_a[E.argmax()])
                    res["R2"].append(pix_a[rng.choice(2 * n_s, p=E / E.sum())])
                for k, v in res.items():
                    v = np.array(v); p1 = np.mean(v == 1)
                    field_rows.append(dict(atoms_per_pixel=n_s, rule=k, share_pixel1=p1, se=np.sqrt(p1 * (1 - p1) / M), born=0.8, count_share=0.5))

    # ------------------------------------------------------------------------- files
    with open(OUT / "summary.csv", "w", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=list(summary[0].keys())); wr.writeheader()
        for row in summary:
            wr.writerow({k: (f"{v:.4f}" if isinstance(v, float) else v) for k, v in row.items()})
    with open(OUT / "field_reading.csv", "w", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=list(field_rows[0].keys())); wr.writeheader()
        for row in field_rows:
            wr.writerow({k: (f"{v:.4f}" if isinstance(v, float) else v) for k, v in row.items()})
    with open(OUT / "current.csv", "w", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=list(current_rows[0].keys())); wr.writeheader()
        for row in current_rows:
            wr.writerow({k: (f"{v:.4f}" if isinstance(v, float) else v) for k, v in row.items()})
    tp = trials_primary
    with open(OUT / "trials.csv", "w", newline="") as f:
        wr = csv.writer(f)
        wr.writerow(["x_start", "x_arrival", "region", "x_contact", "rotation_deg", "w", "carriers_at_electrode", "record",
                     "contact_col", "greedy_end_col", "greedy_T", "dbm_centroid", "dbm_width", "dbm_T", "diff_centroid", "diff_width", "diff_T"])
        g, d, df = tp["cur"][("greedy", True)], tp["cur"][("dbm", True)], tp["cur"][("diffusion", True)]
        for i in range(M):
            wr.writerow([f"{tp['x0'][i]:.4f}", f"{tp['x_arr'][i]:.4f}", tp["reg"][i], f"{tp['x_c'][i]:.4f}", f"{np.degrees(tp['rot'][i]):.3f}",
                         f"{tp['w'][i]:.4f}", tp["bd"][i], int(tp["record"][i]), tp["col0"][i],
                         f"{g[0][i]:.0f}", f"{g[2][i]:.3f}", f"{d[0][i]:.3f}", f"{d[1][i]:.3f}", f"{d[2][i]:.3f}",
                         f"{df[0][i]:.3f}", f"{df[1][i]:.3f}", f"{df[2][i]:.3f}"])

    # ------------------------------------------------------------------------- figures
    xs, x_arr, rho_d, reg = fig1_data[4.0]
    fig, axes = plt.subplots(1, 2, figsize=(9.4, 4.2), gridspec_kw={"width_ratios": [1.15, 1]})
    ax = axes[0]
    for pix, (lo, hi) in PIX.items():
        ax.axvspan(lo, hi, color=SLOT[2], alpha=0.10, lw=0)
    sel = rng.choice(M, 90, replace=False)
    for i in sel:
        ax.plot(xs[:, i], ts, color=SLOT[0], lw=0.6, alpha=0.55)
    ax.set_xlim(-12, 12); ax.set_xlabel("x across the detector plane (Compton lengths)"); ax.set_ylabel("time (Compton times)")
    ax.set_title("Bohm–Dirac trajectories to two SPAD pixels, r = 4 (90 of 4000)")
    ax = axes[1]
    for pix, (lo, hi) in PIX.items():
        ax.axvspan(lo, hi, color=SLOT[2], alpha=0.10, lw=0)
    ax.plot(X, rho_d, color=TEXT, lw=1.4, label="ρ(x, t_d)")
    bins = np.arange(-12, 12.01, 0.4)
    ax.hist(x_arr, bins=bins, density=True, color=SLOT[0], alpha=0.45, label="arrivals")
    xr = tp["x_c"][tp["record"]]
    ax.hist(xr, bins=bins, histtype="step", color=SLOT[1], lw=1.4, label="records per quantum (least rotation, D4a)", weights=np.full(xr.size, 1.0 / (M * 0.4)))
    fr = {k: np.mean(reg == v) for k, v in (("p1", 1), ("p2", 2), ("gap", 0), ("out", -1))}
    ax.set_xlim(-12, 12); ax.set_xlabel("x (Compton lengths)"); ax.set_ylabel("density")
    ax.set_title(f"pixel 1 {fr['p1']:.3f}, pixel 2 {fr['p2']:.3f}, gap {fr['gap']:.3f}, outside {fr['out']:.3f}")
    ax.legend(loc="upper left")
    fig.tight_layout(); fig.savefig(OUT / "fig1_arrivals_pixels.png", dpi=150); plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(9.4, 3.9))
    ax = axes[0]
    for rule, c in (("least", SLOT[0]), ("R2", SLOT[1])):
        rows = [s for s in summary if s["r"] == 4.0 and s["rule"] == rule and s["deposit"] == "D4a"]
        ax.plot([s["n"] for s in rows], [s["rotation_mean_deg"] for s in rows], "o-", color=c, label=f"{'least rotation' if rule == 'least' else 'R2, sampled'}")
    ns = np.array(N_LIST); ax.plot(ns, np.degrees(np.pi / (ns + 1)), ":", color=TEXT2, label="π/(n + 1)")
    ax.axhline(53.5, color=TEXT2, lw=0.8, ls="--"); ax.text(1.1, 55, "53.5°", color=TEXT2, fontsize=8)
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("atoms in reach n"); ax.set_ylabel("rotation used (degrees)")
    ax.set_title("Rotation needed to match the second spinor"); ax.legend()
    ax = axes[1]
    for rule, c in (("least", SLOT[0]), ("R2", SLOT[1])):
        rows = [s for s in summary if s["r"] == 4.0 and s["rule"] == rule and s["deposit"] == "D4b"]
        ax.plot([s["n"] for s in rows], [max(s["miss_rate"], 1e-4) for s in rows], "o-", color=c, label=f"{'least rotation' if rule == 'least' else 'R2, sampled'}")
    ax.plot(ns, np.maximum(0.468 ** ns, 1e-4), ":", color=TEXT2, label="0.468ⁿ")
    ax.axhline(0.152, color=TEXT2, lw=0.8, ls="--"); ax.text(1.1, 0.17, "0.152", color=TEXT2, fontsize=8)
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_ylim(5e-5, 1); ax.set_xlabel("atoms in reach n"); ax.set_ylabel("no pair (deposit below the gap), D4b")
    ax.set_title("The miss under the weighted deposit, θ = 0.45"); ax.legend(loc="lower left")
    fig.tight_layout(); fig.savefig(OUT / "fig2_rotation_miss.png", dpi=150); plt.close(fig)

    fig, axes = plt.subplots(1, 3, figsize=(12.5, 3.9))
    ax = axes[0]
    # one example quantum in pixel 1, the three transfer rules on the fixed landscape
    i_ex = np.nonzero((tp["reg"] == 1) & tp["record"])[0][0]; c0 = tp["col0"][i_ex]
    ax.imshow(np.log(res_fixed[1]), aspect="auto", cmap="Greys", alpha=0.35, extent=[-0.5, W - 0.5, H + 0.5, -0.5])
    for rule_c, c in zip(RULES, SLOT):
        n = np.zeros((1, W)); n[0, c0] = 1.0; col = np.array([c0])
        path_c, path_w = [c0], [0.0]
        for k in range(H):
            if rule_c == "greedy":
                cands = np.array([col[0] - 1, col[0], col[0] + 1]); valid = (cands >= 0) & (cands < W)
                cc = np.clip(cands, 0, W - 1); rr = np.where(valid, res_fixed[1][k + 1][cc], np.inf); col = np.array([cc[rr.argmin()]])
                path_c.append(col[0]); path_w.append(0.0)
            else:
                n = transfer_step(n, res_fixed[1][k + 1], rule_c); cen = (n[0] * np.arange(W)).sum(); wid = np.sqrt((n[0] * (np.arange(W) - cen) ** 2).sum())
                path_c.append(cen); path_w.append(wid)
        path_c, path_w = np.array(path_c), np.array(path_w)
        ax.plot(path_c, np.arange(H + 1), color=c, lw=1.6, label=rule_c)
        if rule_c != "greedy":
            ax.fill_betweenx(np.arange(H + 1), path_c - path_w, path_c + path_w, color=c, alpha=0.15, lw=0)
    ax.set_ylim(H, 0); ax.set_xlim(max(0, c0 - 30), min(W - 1, c0 + 30)); ax.set_xlabel("column"); ax.set_ylabel("row (depth)")
    ax.set_title("One quantum: the charge's path down the pixel"); ax.legend(loc="lower right")
    ax = axes[1]
    g = tp["cur"][("greedy", True)][0]; ii = (tp["reg"] == 1)
    ax.hist(g[ii], bins=np.arange(-0.5, W, 1), color=SLOT[0], alpha=0.8, label="greedy, fixed landscape")
    ax.hist(tp["col0"][ii], bins=np.arange(-0.5, W, 1), histtype="step", color=TEXT2, lw=1.0, label="contact columns")
    ax.set_xlabel("electrode column, pixel 1"); ax.set_ylabel("quanta"); ax.set_title("Where the charge lands: device channels"); ax.legend()
    ax = axes[2]
    labels, wid_m, off_s = [], [], []
    for row in current_rows:
        labels.append(f"{row['rule']}\n{row['landscape']}"); wid_m.append(row["width_mean"]); off_s.append(row["offset_sd"])
    xx = np.arange(len(labels))
    ax.bar(xx - 0.18, wid_m, 0.36, color=SLOT[0], label="electrode width (cells)")
    ax.bar(xx + 0.18, off_s, 0.36, color=SLOT[1], label="offset from contact, sd (cells)")
    ax.set_xticks(xx); ax.set_xticklabels(labels, fontsize=7.5); ax.set_title("Spread and wander of the record"); ax.legend()
    fig.tight_layout(); fig.savefig(OUT / "fig3_current.png", dpi=150); plt.close(fig)

    text = "\n".join(lines) + "\n"
    (OUT / "verification.txt").write_text(text)
    ok = all("FAIL" not in l for l in lines)
    (OUT / "run_info.txt").write_text(f"spad_contact.py, run {time.strftime('%Y-%m-%d %H:%M:%S')}, seed 20260908, {time.time() - t0:.1f} s, numpy {np.__version__}\nverification: {'PASS' if ok else 'FAIL'}\n")
    print(text)
    print("summary (r = 4, D4a and D4b):")
    for s in summary:
        print(f"  r={s['r']:g} n={s['n']:>4} {s['rule']:<5} {s['deposit']}: records {s['records']:>4} eff {s['efficiency']:.3f} miss {s['miss_rate']:.4f} P_b {s['P_b']:.3f} share {s['share_pixel1']:.4f} ± {s['share_se']:.4f} (ρ {s['share_rho']:.4f}) rot {s['rotation_mean_deg']:.2f}°")
    print("field reading (C1):")
    for fr_ in field_rows:
        print(f"  {fr_['atoms_per_pixel']:>3} atoms/pixel {fr_['rule']:<12} share {fr_['share_pixel1']:.4f} ± {fr_['se']:.4f}")
    print("current:")
    for c in current_rows:
        print(f"  {c['rule']:<9} {c['landscape']}: width {c['width_mean']:.2f} offset sd {c['offset_sd']:.2f} corr {c['corr_contact_centroid']:.3f} time/row {c['time_per_row']:.3f} distinct {c['distinct_end_columns_pixel1']}")
    print(f"{time.time() - t0:.1f} s")
    return 0 if ok else 1

if __name__ == "__main__":
    sys.exit(main())
