#!/usr/bin/env python3
"""
sink_field.py: sink field and channel growth in two dimensions (PREREGISTRATION.md, read it first)
==================================================================================================

A transport model with no physics labels: a density ρ on a lattice of atoms with random phases;
the second study's deposit race seeds one sink; every atom charges from ρ × phase weight and
latches into an absorber at a threshold; charged sites source a screened field; the density flows
along the field with a conductance tensor that grows along the flow it has carried.
"""
import csv, sys, datetime, platform
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

HERE = Path(__file__).resolve().parent
OUT = HERE / "results"
OUT.mkdir(exist_ok=True)

N, D_DIFF, SIGMA, XL, XR, YC = 96, 0.5, 6.0, 32.0, 64.0, 48.0
DT, T_END, TAU_C, TAU_N, GAMMA0, KAPPA, V_CAP = 0.1, 60.0, 10.0, 1.0, 1.0, 10.0, 0.5
PARAM_SETS = {"registered": dict(mu=200.0, cth=0.01), "addendum": dict(mu=2000.0, cth=0.004)}
xs = np.arange(N); XX, YY = np.meshgrid(xs, xs, indexing="ij")
kx = 2 * np.pi * np.fft.fftfreq(N); KX, KY = np.meshgrid(kx, kx, indexing="ij")
LAP_K = -(2 - 2 * np.cos(KX)) - (2 - 2 * np.cos(KY))          # discrete Laplacian symbol

def initial_density(r):
    g = lambda x0: np.exp(-((XX - x0) ** 2 + (YY - YC) ** 2) / (2 * SIGMA ** 2))
    rho = r * g(XR) + g(XL)
    return rho / rho.sum()

def solve_field(q, ell):
    """(∇² − 1/ℓ²)Φ = −q on the periodic grid; returns Φ and the field +∇Φ, which points toward the
    charges (addendum 1: the registered formula had the sign of an electrostatic field, which points away)."""
    Phi = np.real(np.fft.ifft2(np.fft.fft2(q) / (1.0 / ell ** 2 - LAP_K)))
    Fx = (np.roll(Phi, -1, 0) - np.roll(Phi, 1, 0)) / 2
    Fy = (np.roll(Phi, -1, 1) - np.roll(Phi, 1, 1)) / 2
    return Phi, Fx, Fy

def step_transport(rho, vx, vy):
    """Forward Euler: diffusion (5-point) plus upwind drift; conservative fluxes."""
    lap = np.roll(rho, 1, 0) + np.roll(rho, -1, 0) + np.roll(rho, 1, 1) + np.roll(rho, -1, 1) - 4 * rho
    # face velocities (average to faces), upwind densities
    vxf = 0.5 * (vx + np.roll(vx, -1, 0)); vyf = 0.5 * (vy + np.roll(vy, -1, 1))
    fx = np.where(vxf > 0, rho, np.roll(rho, -1, 0)) * vxf
    fy = np.where(vyf > 0, rho, np.roll(rho, -1, 1)) * vyf
    div = (fx - np.roll(fx, 1, 0)) + (fy - np.roll(fy, 1, 1))
    return rho + DT * (D_DIFF * lap - div), fx, fy

def run(r, ell, tensor, mu, M0, s, rng, C_TH, snapshots=None):
    rho = initial_density(r)
    delta = rng.uniform(-np.pi, np.pi, (N, N)); w = np.cos(delta / 2) ** 2
    # seed by R2: one atom sampled ∝ ρ w, initial charge w × threshold
    p = (rho * w).ravel(); seed = rng.choice(N * N, p=p / p.sum()); si, sj = divmod(seed, N)
    C = np.zeros((N, N)); C[si, sj] = w[si, sj] * C_TH
    E = np.zeros((N, N)); eta = np.zeros((N, N)); latched = np.zeros((N, N), bool)
    Mxx = np.full((N, N), float(M0)); Myy = np.full((N, N), float(M0)); Mxy = np.zeros((N, N))
    t_first = [np.nan, np.nan]; mass_err = 0.0; snaps = {}
    n_steps = int(round(T_END / DT))
    for k in range(n_steps + 1):
        t = k * DT
        if snapshots is not None and any(abs(t - ts) < DT / 2 for ts in snapshots):
            snaps[round(t)] = (rho.copy(), (Mxx + Myy).copy(), latched.copy())
        if k == n_steps:
            break
        q = C + E
        if M0 > 0 or mu > 0:
            _, Fx, Fy = solve_field(q, ell)
            vx = Mxx * Fx + Mxy * Fy; vy = Mxy * Fx + Myy * Fy
            sp = np.hypot(vx, vy); fac = np.where(sp > V_CAP, V_CAP / np.maximum(sp, 1e-30), 1.0)
            vx, vy = vx * fac, vy * fac
        else:
            vx = vy = np.zeros((N, N))
        before = rho.sum()
        rho, fx, fy = step_transport(rho, vx, vy)
        # conductance growth with the drift flow that passed (cell-centred magnitude and direction)
        if mu > 0:
            jx = 0.5 * (fx + np.roll(fx, 1, 0)); jy = 0.5 * (fy + np.roll(fy, 1, 1)); jm = np.hypot(jx, jy)
            if tensor:
                with np.errstate(invalid="ignore", divide="ignore"):
                    ux = np.where(jm > 0, jx / jm, 0.0); uy = np.where(jm > 0, jy / jm, 0.0)
                Mxx += mu * DT * jm * ux * ux; Myy += mu * DT * jm * uy * uy; Mxy += mu * DT * jm * ux * uy
            else:
                Mxx += mu * DT * jm; Myy += mu * DT * jm
        # absorption at latched sites
        G = np.where(latched, GAMMA0 * (1 + KAPPA * (C + E)), 0.0)
        keep = np.exp(-G * DT); absorbed = rho * (1 - keep); rho = rho * keep; E += absorbed
        mass_err = max(mass_err, abs(rho.sum() + absorbed.sum() - before))
        # charging, noise, latch
        eta += -eta * DT / TAU_N + s * np.sqrt(DT) * rng.standard_normal((N, N))
        C += DT * (rho * w - C / TAU_C)
        new = (~latched) & (C + eta >= C_TH); latched |= new
        for side, mask in ((0, XX < N // 2), (1, XX > N // 2)):
            if np.isnan(t_first[side]) and (new & mask).any():
                t_first[side] = t
    # metrics
    coll = [E[XX < N // 2].sum(), E[XX > N // 2].sum()]
    n_lat = [int(latched[XX < N // 2].sum()), int(latched[XX > N // 2].sum())]
    both = all(n > 0 for n in n_lat)
    first = -1 if all(np.isnan(t_first)) else int(np.nanargmin(np.where(np.isnan(t_first), np.inf, t_first)))
    tr = Mxx + Myy; heated = tr > 3 * M0 if M0 > 0 else np.zeros((N, N), bool)
    if heated.any():
        ii, jj = np.where(heated); Lx = ii.max() - ii.min() + 1; Ly = jj.max() - jj.min() + 1
        index = Lx / Ly; area = int(heated.sum())
        lam1 = 0.5 * (tr + np.sqrt((Mxx - Myy) ** 2 + 4 * Mxy ** 2)); lam2 = tr - lam1
        aniso = float(np.mean(((lam1 - lam2) / np.maximum(tr, 1e-30))[heated]))
    else:
        index, area, aniso = np.nan, 0, np.nan
    return dict(first=first, both=both, n_lat=n_lat, coll=coll, t_first=t_first, index=index, area=area, aniso=aniso,
                mass_err=mass_err, seed_side=int(si > N // 2), snaps=snaps, rho=rho, tr=tr, latched=latched)

def settings(mu):
    return [
        ("base",        4.0, 16.0, True,  mu,  10.0, 1e-3),
        ("range_4",     4.0, 4.0,  True,  mu,  10.0, 1e-3),
        ("range_64",    4.0, 64.0, True,  mu,  10.0, 1e-3),
        ("no_field",    4.0, 16.0, True,  0.0, 0.0,  1e-3),
        ("scalar_M",    4.0, 16.0, False, mu,  10.0, 1e-3),
        ("no_growth",   4.0, 16.0, True,  0.0, 10.0, 1e-3),
        ("r_1",         1.0, 16.0, True,  mu,  10.0, 1e-3),
        ("r_1_range_64",1.0, 64.0, True,  mu,  10.0, 1e-3),
        ("no_noise",    4.0, 16.0, True,  mu,  10.0, 0.0),
        ("r_1_no_noise",1.0, 16.0, True,  mu,  10.0, 0.0),
    ]
N_TRIALS = 20

def verification(lines):
    # Poisson solve against a point source: (∇² − 1/ℓ²)Φ = −δ  ->  Σ(−LΦ + Φ/ℓ²) = 1 and symmetry
    q = np.zeros((N, N)); q[48, 48] = 1.0; Phi, _, _ = solve_field(q, 16.0)
    lap = np.roll(Phi, 1, 0) + np.roll(Phi, -1, 0) + np.roll(Phi, 1, 1) + np.roll(Phi, -1, 1) - 4 * Phi
    res = np.abs(lap - Phi / 16.0 ** 2 + q).max()
    lines.append(f"Poisson solve, point source at the centre: max residual of (∇² − 1/ℓ²)Φ + q = {res:.1e}; Φ(48,48) = {Phi[48,48]:.4f}, Φ(64,48) = {Phi[64,48]:.4f}, Φ(48,64) = {Phi[48,64]:.4f}   -> {'PASS' if res < 1e-10 and abs(Phi[64,48]-Phi[48,64]) < 1e-12 else 'FAIL'}")
    # diffusion only against the Gaussian solution
    rho = initial_density(1.0); rho = rho / rho.sum()
    for _ in range(int(round(T_END / DT))):
        rho, _, _ = step_transport(rho, np.zeros((N, N)), np.zeros((N, N)))
    s2 = SIGMA ** 2 + 2 * D_DIFF * T_END
    g = lambda x0: np.exp(-((XX - x0) ** 2 + (YY - YC) ** 2) / (2 * s2)); ref = g(XR) + g(XL); ref /= ref.sum()
    derr = np.abs(rho - ref).max() / ref.max()
    lines.append(f"Diffusion only for t = 60 against the Gaussian solution (σ² → σ² + 2Dt): max relative error {derr:.2e}   -> {'PASS' if derr < 0.02 else 'FAIL'}")
    return res < 1e-10 and derr < 0.02

def main():
    rng = np.random.default_rng(20260907 + 11)
    lines = ["Verification (PREREGISTRATION.md §4)", ""]
    ok = verification(lines)
    for tag, prm in PARAM_SETS.items():
        rows, summary, worst_mass = [], [], 0.0
        fig_snaps = None
        lines.append(""); lines.append(f"Parameter set '{tag}': μ = {prm['mu']:g}, C_th = {prm['cth']:g}")
        for name, r, ell, tensor, mu, M0, s in settings(prm["mu"]):
            outs = []
            for j in range(N_TRIALS):
                res = run(r, ell, tensor, mu, M0, s, rng, prm["cth"], snapshots=[0, 10, 30, 60] if (name == "base" and j == 0) else None)
                if name == "base" and j == 0:
                    fig_snaps = res
                outs.append(res); worst_mass = max(worst_mass, res["mass_err"])
                rows.append([tag, name, j, res["seed_side"], res["first"], int(res["both"]), res["n_lat"][0], res["n_lat"][1],
                             f"{res['coll'][0]:.4f}", f"{res['coll'][1]:.4f}", f"{res['t_first'][0]:.1f}", f"{res['t_first'][1]:.1f}",
                             f"{res['index']:.2f}", res["area"], f"{res['aniso']:.3f}"])
            both = np.mean([o["both"] for o in outs]); bright_first = np.mean([o["first"] == 1 for o in outs])
            none = np.mean([o["first"] == -1 for o in outs]); seed_bright = np.mean([o["seed_side"] == 1 for o in outs])
            win_coll = np.mean([o["coll"][o["first"]] if o["first"] >= 0 else 0.0 for o in outs])
            tot_coll = np.mean([sum(o["coll"]) for o in outs])
            idx = np.nanmean([o["index"] for o in outs]) if any(not np.isnan(o["index"]) for o in outs) else np.nan
            an = np.nanmean([o["aniso"] for o in outs]) if any(not np.isnan(o["aniso"]) for o in outs) else np.nan
            area = np.mean([o["area"] for o in outs]); nl = np.mean([o["n_lat"][0] for o in outs]); nr = np.mean([o["n_lat"][1] for o in outs])
            summary.append([tag, name, r, ell, int(tensor), mu, M0, s, N_TRIALS, f"{seed_bright:.2f}", f"{bright_first:.2f}", f"{both:.2f}", f"{none:.2f}",
                            f"{win_coll:.3f}", f"{tot_coll:.3f}", f"{idx:.2f}", f"{area:.0f}", f"{an:.3f}", f"{nl:.1f}", f"{nr:.1f}"])
            lines.append(f"  {name}: seed in bright lobe {seed_bright:.2f}, bright lobe first {bright_first:.2f}, double record {both:.2f}, none {none:.2f}, "
                         f"winner collects {win_coll:.3f} (total {tot_coll:.3f}), latched cells L/R {nl:.1f}/{nr:.1f}, channel index {idx:.2f}, heated area {area:.0f}, anisotropy {an:.3f}")
        lines.append(f"  mass conservation per step (ρ + absorbed), worst over the set: {worst_mass:.1e}   -> {'PASS' if worst_mass < 1e-6 else 'FAIL'}")
        ok = ok and worst_mass < 1e-6
        with open(OUT / f"summary_{tag}.csv", "w", newline="") as f:
            wr = csv.writer(f)
            wr.writerow(["set", "setting", "r", "range_ell", "tensor", "mu", "M0", "noise_s", "trials", "P_seed_bright", "P_bright_first", "P_double_record",
                         "P_no_record", "winner_collected_t60_mean", "total_collected_t60_mean", "channel_index_mean", "heated_area_mean", "anisotropy_mean",
                         "latched_left_mean", "latched_right_mean"])
            wr.writerows(summary)
        with open(OUT / f"trials_{tag}.csv", "w", newline="") as f:
            wr = csv.writer(f)
            wr.writerow(["set", "setting", "trial", "seed_side_0left_1right", "first_lobe_0left_1right_-1none", "double_record", "latched_left", "latched_right",
                         "collected_left", "collected_right", "t_first_left", "t_first_right", "channel_index", "heated_area", "anisotropy"])
            wr.writerows(rows)
        figures(fig_snaps, summary, tag, prm)
    (OUT / "verification.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    info = [f"run: {datetime.datetime.now().isoformat(timespec='seconds')}", f"python {platform.python_version()}, numpy {np.__version__}",
            f"verification: {'PASS' if ok else 'FAIL'}"]
    (OUT / "run_info.txt").write_text("\n".join(info) + "\n"); print("\n".join(info))
    return 0 if ok else 1

SURFACE, TEXT, TEXT2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e6e5e1"
SLOT = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]
SEQ = ["#fcfcfb", "#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]

def figures(res, summary, tag, prm):
    plt.rcParams.update({
        "font.family": "sans-serif", "font.size": 9, "axes.titlesize": 9.5, "axes.titleweight": "bold",
        "axes.edgecolor": GRID, "axes.labelcolor": TEXT2, "xtick.color": TEXT2, "ytick.color": TEXT2,
        "axes.titlecolor": TEXT, "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
        "savefig.facecolor": SURFACE, "axes.grid": False, "legend.frameon": False, "legend.fontsize": 7.5,
    })
    cmap = LinearSegmentedColormap.from_list("seqblue", SEQ)
    times = sorted(res["snaps"])
    fig, axes = plt.subplots(2, len(times), figsize=(3.1 * len(times), 6.2))
    vmax = max(res["snaps"][t][0].max() for t in times)
    for a, t in enumerate(times):
        rho, tr, lat = res["snaps"][t]
        ax = axes[0, a]; ax.imshow(rho.T, origin="lower", cmap=cmap, vmin=0, vmax=vmax)
        ii, jj = np.where(lat); ax.plot(ii, jj, ".", color=SLOT[1], ms=2.5)
        ax.set_title(f"density and latched sites, t = {t}"); ax.set_xticks([]); ax.set_yticks([])
        ax = axes[1, a]; ax.imshow(np.log10(np.maximum(tr, 1e-9)).T, origin="lower", cmap=cmap, vmin=1, vmax=4)
        ax.set_title(f"conductance trace (log10), t = {t}"); ax.set_xticks([]); ax.set_yticks([])
    fig.suptitle(f"Base setting, first trial, {tag} set (μ = {prm['mu']:g}, C_th = {prm['cth']:g}): r = 4, range 16, tensor conductance", color=TEXT, fontsize=10, fontweight="bold")
    fig.tight_layout(); fig.savefig(OUT / f"fig1_fields_{tag}.png", dpi=150); plt.close(fig)
    names = [s[1] for s in summary]; x = np.arange(len(names))
    fig, axes = plt.subplots(1, 3, figsize=(13.2, 4.0))
    ax = axes[0]
    ax.bar(x - 0.2, [float(s[11]) for s in summary], 0.4, color=SLOT[1], label="double record")
    ax.bar(x + 0.2, [float(s[10]) for s in summary], 0.4, color=SLOT[0], label="bright lobe first")
    ax.axhline(0.8, color=TEXT2, lw=0.8, alpha=0.6); ax.text(len(names) - 0.5, 0.82, "Born 0.80", color=TEXT2, fontsize=7.5, ha="right")
    ax.set_xticks(x); ax.set_xticklabels([n.replace("_", " ") for n in names], rotation=35, ha="right"); ax.set_ylim(0, 1.05); ax.set_title("Records per quantum"); ax.legend(loc="upper left")
    ax = axes[1]
    ax.bar(x - 0.2, [float(s[13]) for s in summary], 0.4, color=SLOT[2], label="collected by the winning lobe")
    ax.bar(x + 0.2, [float(s[14]) for s in summary], 0.4, color=SLOT[3], label="collected in total")
    ax.set_xticks(x); ax.set_xticklabels([n.replace("_", " ") for n in names], rotation=35, ha="right"); ax.set_ylim(0, 1.05); ax.set_title("Collection by t = 60"); ax.legend(loc="lower left")
    ax = axes[2]
    idx = [float(s[15]) if s[15] != "nan" else 0 for s in summary]; an = [float(s[17]) if s[17] != "nan" else 0 for s in summary]
    ax.bar(x - 0.2, idx, 0.4, color=SLOT[0], label="channel index (length / width)")
    ax.bar(x + 0.2, [10 * a for a in an], 0.4, color=SLOT[1], label="anisotropy × 10")
    ax.axhline(3, color=TEXT2, lw=0.8, alpha=0.6); ax.text(len(names) - 0.5, 3.15, "line above 3", color=TEXT2, fontsize=7.5, ha="right")
    ax.set_xticks(x); ax.set_xticklabels([n.replace("_", " ") for n in names], rotation=35, ha="right"); ax.set_title("Is the heated region a line?"); ax.legend(loc="upper right")
    fig.suptitle(f"{tag} set (μ = {prm['mu']:g}, C_th = {prm['cth']:g})", color=TEXT, fontsize=10, fontweight="bold")
    fig.tight_layout(); fig.savefig(OUT / f"fig2_outcomes_{tag}.png", dpi=150); plt.close(fig)

if __name__ == "__main__":
    sys.exit(main())
