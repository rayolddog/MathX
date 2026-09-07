#!/usr/bin/env python3
"""
bohm_dirac_flight.py: free flight by the Bohm–Dirac velocity, arrival, and the local race (Part C)
=================================================================================================

Definitions D8 and D9 of PREREGISTRATION.md. Conventions of the record's
../dk_origin/dirac_arnold_tongue.py:  ħ = c = m = 1,  H = σ_x p + σ_z m,  velocity ψ†σ_x ψ / ψ†ψ.
The free propagator is exact in momentum space, e^(−iHt) = cos(Et) − i sin(Et) H/E, so ψ(x, t)
is available at any t and the trajectories are integrated by RK4 against the exact field.

  first spinor   positive-energy component, rotates as e^(−iEt)
  second spinor  negative-energy component, rotates as e^(+iEt)
  offset δ       their relative phase; the free equation does not change it
"""
import csv, sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
OUT = HERE / "results"
OUT.mkdir(exist_ok=True)
rng = np.random.default_rng(20260907)

NX, DX = 4096, 0.02
X = (np.arange(NX) - NX // 2) * DX
P = 2 * np.pi * np.fft.fftfreq(NX, d=DX)
M_ = 1.0
E = np.sqrt(P ** 2 + M_ ** 2)

def gauss(x, x0, sigma=1.0):
    return (2 * np.pi * sigma ** 2) ** (-0.25) * np.exp(-(x - x0) ** 2 / (4 * sigma ** 2))

def to_p(psi):                      # psi: (2, NX) in x
    return np.fft.fft(psi, axis=1)

def psi_at(psi0_p, t):
    c, s = np.cos(E * t), np.sin(E * t)
    a, b = psi0_p
    u = (c - 1j * s * M_ / E) * a + (-1j * s * P / E) * b
    d = (-1j * s * P / E) * a + (c + 1j * s * M_ / E) * b
    return np.fft.ifft(np.array([u, d]), axis=1)

def rho_v(psi):
    rho = (np.abs(psi[0]) ** 2 + np.abs(psi[1]) ** 2)
    j = 2 * np.real(np.conj(psi[0]) * psi[1])
    v = np.where(rho > 1e-14, j / np.maximum(rho, 1e-300), 0.0)
    return rho, v

def normalise(psi):
    n = np.sqrt(((np.abs(psi) ** 2).sum()) * DX)
    return psi / n

def sample_positions(rho, M):
    cdf = np.cumsum(rho) * DX; cdf /= cdf[-1]
    return np.interp(rng.uniform(0, 1, M), cdf, X)

def integrate(psi0_p, x0, t_end, dt):
    """RK4 trajectories against the exact field; returns times and x(t) for every sample."""
    n = int(round(t_end / dt)); xs = np.empty((n + 1, x0.size)); xs[0] = x0; x = x0.copy()
    vel = lambda xx, tt: np.interp(xx, X, rho_v(psi_at(psi0_p, tt))[1])
    for k in range(n):
        t = k * dt
        k1 = vel(x, t); k2 = vel(x + 0.5 * dt * k1, t + 0.5 * dt)
        k3 = vel(x + 0.5 * dt * k2, t + 0.5 * dt); k4 = vel(x + dt * k3, t + dt)
        x = x + dt * (k1 + 2 * k2 + 2 * k3 + k4) / 6
        xs[k + 1] = x
    return np.arange(n + 1) * dt, xs

# ----------------------------------------------------------------------------- style
SURFACE, TEXT, TEXT2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e6e5e1"
SLOT = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]
plt.rcParams.update({
    "font.family": "sans-serif", "font.size": 9, "axes.titlesize": 10, "axes.titleweight": "bold",
    "axes.edgecolor": GRID, "axes.labelcolor": TEXT2, "xtick.color": TEXT2, "ytick.color": TEXT2,
    "axes.titlecolor": TEXT, "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
    "savefig.facecolor": SURFACE, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.8,
    "grid.linestyle": "-", "axes.spines.top": False, "axes.spines.right": False,
    "legend.frameon": False, "legend.fontsize": 8, "axes.axisbelow": True,
})

# ----------------------------------------------------------------------------- Part C: two lobes
R_INT, X1, T_D, DT, M_TRAJ = 4.0, 4.0, 4.0, 0.01, 4000
LATTICE, REACH = 0.05, 0.5

def two_lobe_flight(lines):
    g = np.sqrt(R_INT) * gauss(X, X1) + gauss(X, -X1)
    psi0 = normalise(np.array([g, np.zeros_like(g)], dtype=complex))
    psi0_p = to_p(psi0)
    rho0, _ = rho_v(psi0)
    x0 = sample_positions(rho0, M_TRAJ)
    ts, xs = integrate(psi0_p, x0, T_D, DT)
    psi_d = psi_at(psi0_p, T_D); rho_d, _ = rho_v(psi_d)
    norm_err = abs(rho_d.sum() * DX - 1)
    frac_arr = np.mean(xs[-1] > 0); frac_rho = rho_d[X > 0].sum() * DX
    se = np.sqrt(frac_rho * (1 - frac_rho) / M_TRAJ)
    # whole-distribution check: max |CDF difference|
    cdf_rho = np.cumsum(rho_d) * DX
    cdf_arr = np.searchsorted(np.sort(xs[-1]), X, side="right") / M_TRAJ
    ks = np.max(np.abs(cdf_rho - cdf_arr))
    lines.append(f"Part C: norm at t_d: |∫ρ − 1| = {norm_err:.1e};  arrivals at x > 0: {frac_arr:.4f}, ∫_(x>0) ρ(t_d) = {frac_rho:.4f}"
                 f" (s.e. {se:.4f}), |diff| = {abs(frac_arr-frac_rho):.4f}   -> {'PASS' if abs(frac_arr-frac_rho) < 3*se else 'FAIL'}")
    lines.append(f"Part C: max |CDF(arrivals) − CDF(ρ(t_d))| = {ks:.4f}   (tolerance 0.03)   -> {'PASS' if ks < 0.03 else 'FAIL'}")
    lines.append(f"Part C: r/(1+r) = {R_INT/(1+R_INT):.3f}; the lobes' spreading across x = 0 by t_d moves ∫_(x>0) ρ(t_d) to {frac_rho:.4f}")
    # local race at arrival: atoms on the lattice within REACH
    n_reach = int(2 * REACH / LATTICE) + 1
    arr = xs[-1]
    base = np.round(arr / LATTICE) * LATTICE
    offsets = (np.arange(n_reach) - n_reach // 2) * LATTICE
    atoms = base[:, None] + offsets[None, :]
    dphi = rng.uniform(-np.pi, np.pi, atoms.shape)
    w = np.cos(dphi / 2) ** 2
    cdf = np.cumsum(w, axis=1); u = rng.uniform(0, 1, M_TRAJ) * cdf[:, -1]
    k2 = (cdf < u[:, None]).sum(axis=1); k3 = w.argmax(axis=1)
    x_r2 = atoms[np.arange(M_TRAJ), k2]; x_r3 = atoms[np.arange(M_TRAJ), k3]
    w_r2 = w[np.arange(M_TRAJ), k2]; w_r3 = w.max(axis=1)
    lines.append(f"Part C: winner's energy weight, mean: R2 {w_r2.mean():.4f} (predicted 0.75), R3 {w_r3.mean():.4f} (predicted ≈ 0.99);"
                 f" blur x_record − x_arrival: R2 mean {np.mean(x_r2-arr):+.4f}, sd {np.std(x_r2-arr):.4f}; R3 mean {np.mean(x_r3-arr):+.4f}, sd {np.std(x_r3-arr):.4f}"
                 f" (uniform on ±{REACH}: sd {REACH/np.sqrt(3):.4f})")
    with open(OUT / "flight_arrivals.csv", "w", newline="") as f:
        wr = csv.writer(f)
        wr.writerow(["x_start", "x_arrival_t_d", "x_record_R2", "w_winner_R2", "x_record_R3", "w_winner_R3"])
        for row in zip(x0, arr, x_r2, w_r2, x_r3, w_r3):
            wr.writerow([f"{v:.4f}" for v in row])
    # figure
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.2), gridspec_kw={"width_ratios": [1.15, 1]})
    ax = axes[0]
    sel = rng.choice(M_TRAJ, 90, replace=False)
    for i in sel:
        ax.plot(xs[:, i], ts, color=SLOT[0], lw=0.6, alpha=0.55)
    ax.set_xlim(-12, 12); ax.set_xlabel("x across the aperture (Compton lengths)"); ax.set_ylabel("time (Compton times)")
    ax.set_title("Bohm–Dirac trajectories, two lobes r = 4 (90 of 4000)")
    ax = axes[1]
    ax.plot(X, rho_d, color=TEXT, lw=1.4, label="ρ(x, t_d)")
    bins = np.arange(-12, 12.01, 0.4)
    ax.hist(arr, bins=bins, density=True, color=SLOT[0], alpha=0.45, label="arrivals x(t_d)")
    ax.hist(x_r3, bins=bins, density=True, histtype="step", color=SLOT[1], lw=1.4, label="record position, R3")
    ax.set_xlim(-12, 12); ax.set_xlabel("x (Compton lengths)"); ax.set_ylabel("density")
    ax.set_title(f"Arrivals and records at t_d = {T_D:g}: {frac_arr:.3f} at x > 0, ∫ρ = {frac_rho:.3f}")
    ax.legend(loc="upper left")
    fig.tight_layout(); fig.savefig(OUT / "fig4_flight.png", dpi=150); plt.close(fig)

# ----------------------------------------------------------------------------- Part C: zitter
THETA_DEG, DELTAS_DEG, T_Z, DT_Z = 60.0, [0.0, 90.0], 12.0, 0.005
SIGMA_LIST = [(1.0, "registered"), (3.0, "addendum 2026-09-07")]

def fit_cos(t, v, w_fixed=None):
    """v ≈ R cos(ω t + φ) + C by least squares; ω free on a grid unless fixed. Returns R, φ, ω."""
    def resid(w):
        G = np.column_stack([np.cos(w * t), np.sin(w * t), np.ones_like(t)])
        coef, *_ = np.linalg.lstsq(G, v, rcond=None)
        r = v - G @ coef
        return coef, r @ r
    w_best = w_fixed if w_fixed is not None else min(np.arange(1.5, 3.0, 0.002), key=lambda w: resid(w)[1])
    (A, B, C), _ = resid(w_best)
    return np.hypot(A, B), np.arctan2(-B, A), w_best

def mean_energy(psi0_p):
    wgt = np.abs(psi0_p[0]) ** 2 + np.abs(psi0_p[1]) ** 2
    return (E * wgt).sum() / wgt.sum()

def pm_components(psi_p):
    """Amplitudes on the ±E eigenvectors of H(p) = σ_x p + σ_z m, for every p."""
    nrm = np.sqrt(2 * E * (E + M_))
    up = np.array([(E + M_) / nrm, P / nrm]); um = np.array([-P / nrm, (E + M_) / nrm])
    return (up.conj() * psi_p).sum(axis=0), (um.conj() * psi_p).sum(axis=0)

def offset_check(psi0_p, delta):
    """The offset arg(c₋/c₊) − 2E(p)t, weighted over the packet: its spread over p at t = 0 and its
    drift from t = 0 to t = 12 (the free equation should leave it exactly unchanged)."""
    cp0, cm0 = pm_components(psi0_p)
    wgt = np.abs(cp0) ** 2 + np.abs(cm0) ** 2; wgt /= wgt.sum()
    keep = wgt > 1e-6 * wgt.max()
    off0 = np.angle(cm0 / cp0)
    mean0 = np.angle((wgt[keep] * np.exp(1j * off0[keep])).sum())
    spread0 = np.degrees(np.sqrt((wgt[keep] * ((off0[keep] - mean0 + np.pi) % (2 * np.pi) - np.pi) ** 2).sum() / wgt[keep].sum()))
    drift = 0.0
    for t in (6.0, 12.0):
        cp, cm = pm_components(np.fft.fft(psi_at(psi0_p, t), axis=1))
        off_t = np.angle(cm / cp) - 2 * E * t
        d = (off_t - off0 + np.pi) % (2 * np.pi) - np.pi
        drift = max(drift, np.abs(d[keep]).max())
    return np.degrees(mean0), spread0, np.degrees(drift)

def zitter_flight(lines):
    rows = []
    fig, axes = plt.subplots(2, 1, figsize=(9.2, 5.6), sharex=True)
    for sigma, tag in SIGMA_LIST:
        cases = [(THETA_DEG, d, c) for d, c in zip(DELTAS_DEG, SLOT)]
        if sigma == 1.0:
            cases = [(0.0, 0.0, TEXT2)] + cases
        for th_deg, d_deg, c in cases:
            th, d = np.radians(th_deg), np.radians(d_deg)
            g = gauss(X, 0.0, sigma)
            psi0 = normalise(np.array([np.cos(th / 2) * g, np.sin(th / 2) * np.exp(1j * d) * g], dtype=complex))
            psi0_p = to_p(psi0)
            Em = mean_energy(psi0_p)
            x0 = np.array([0.0 if th_deg > 0 else 1.0])
            ts, xs = integrate(psi0_p, x0, T_Z, DT_Z)
            v = np.array([np.interp(xs[k, 0], X, rho_v(psi_at(psi0_p, t))[1]) for k, t in enumerate(ts)])
            if sigma == 1.0:
                lab = "θ = 0: no second spinor (start x = 1)" if th_deg == 0 else f"θ = 60°, δ = {d_deg:g}°"
                axes[0].plot(ts, v, color=c, lw=1.3, label=lab)
                axes[1].plot(ts, xs[:, 0], color=c, lw=1.3, label=lab)
                if th_deg > 0:
                    axes[0].plot(ts, np.sin(th) * np.cos(2 * Em * ts + d), color=c, lw=0.9, ls="--", alpha=0.7)
            half_range = (xs[:, 0].max() - xs[:, 0].min()) / 2
            if th_deg == 0:
                rows.append([tag, sigma, th_deg, d_deg, "", "", "", "", "", "", "", "", f"{np.abs(v).max():.4f}", f"{half_range:.4f}", "0"])
                lines.append(f"Part C zitter σ = {sigma:g}, θ = 0 (reference, start x = 1): max |v| = {np.abs(v).max():.4f}, no oscillation")
                continue
            v00 = v[0]; v00_pred = np.sin(th) * np.cos(d)
            mean_off, spread_off, drift_off = offset_check(psi0_p, d)
            early = ts <= 3.0; late = ts >= T_Z - 3.0
            R1, ph1, w1 = fit_cos(ts[early], v[early])
            _, _, w2 = fit_cos(ts[late], v[late])
            ok_v0 = abs(v00 - v00_pred) < 1e-6
            ok_off = drift_off < 1e-6
            ok_w = abs(w1 - 2 * Em) < 0.10 * 2 * Em
            rows.append([tag, sigma, th_deg, d_deg, f"{v00:.6f}", f"{v00_pred:.6f}", f"{mean_off:.2f}", f"{spread_off:.2f}", f"{drift_off:.1e}",
                         f"{2*Em:.4f}", f"{w1:.4f}", f"{w2:.4f}", f"{R1:.4f}", f"{half_range:.4f}", f"{np.sin(th)/2:.4f}"])
            lines.append(f"Part C zitter σ = {sigma:g} ({tag}), θ = {th_deg:g}°, δ = {d_deg:g}°: v(0, 0) = {v00:.6f} against sin θ cos δ = {v00_pred:.6f} ({'ok' if ok_v0 else 'MISS'});"
                         f" ±E offset over the packet: mean {mean_off:.2f}°, spread {spread_off:.2f}° (δ = {d_deg:g}°), drift to t = 12: {drift_off:.1e}° ({'maintained' if ok_off else 'NOT maintained'});"
                         f" trajectory beat: fitted {w1:.4f} early, {w2:.4f} late, against the registered 2 ({'MISS' if abs(w1 - 2) > 0.02 else 'ok'}) and 2⟨E⟩ = {2*Em:.4f} ({'ok' if ok_w else 'MISS'});"
                         f" amplitude {R1:.4f} (sin θ = {np.sin(th):.4f}); position half-range {half_range:.4f} (sin θ/2 = {np.sin(th)/2:.4f})"
                         f"   -> {'PASS' if ok_v0 and ok_off and ok_w else 'FAIL'}")
    axes[0].set_ylabel("Bohm–Dirac velocity v/c")
    axes[0].set_title("Zitterbewegung of the trajectory, σ = 1; the offset δ is the phase at t = 0. Dashed: sin θ · cos(2⟨E⟩t + δ)")
    axes[0].legend(loc="upper right", ncol=3); axes[0].set_ylim(-1.05, 1.05)
    axes[1].set_ylabel("position x (Compton lengths)"); axes[1].set_xlabel("time (Compton times, ħ/mc²)")
    axes[1].set_title("Position along the trajectory")
    fig.tight_layout(); fig.savefig(OUT / "fig5_zitter.png", dpi=150); plt.close(fig)
    with open(OUT / "flight_zitter.csv", "w", newline="") as f:
        wr = csv.writer(f)
        wr.writerow(["case", "sigma", "theta_deg", "delta_deg", "v_at_t0", "sin_theta_cos_delta", "pm_offset_mean_deg", "pm_offset_spread_deg",
                     "pm_offset_drift_to_t12_deg", "two_mean_E", "beat_fit_early", "beat_fit_late", "beat_amplitude_early",
                     "position_half_range", "predicted_sin_theta_over_2"])
        wr.writerows(rows)

def main():
    lines = ["Verification, Part C (PREREGISTRATION.md §2, free flight)", ""]
    two_lobe_flight(lines)
    zitter_flight(lines)
    text = "\n".join(lines) + "\n"
    with open(OUT / "verification.txt", "a") as f:
        f.write("\n" + text)
    print(text)
    ok = all("FAIL" not in l for l in lines)
    with open(OUT / "run_info.txt", "a") as f:
        f.write(f"verification (bohm_dirac_flight.py): {'PASS' if ok else 'FAIL'}\n")
    return 0 if ok else 1

if __name__ == "__main__":
    sys.exit(main())
