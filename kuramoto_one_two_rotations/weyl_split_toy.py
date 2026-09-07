#!/usr/bin/env python3
"""
weyl_split_toy.py: the two-Weyl-spinor split, as posited by JB on 2026-09-06
=============================================================================

Definition D8 of PREREGISTRATION.md, nothing more:

  upper component (ψ_L)  phase set equal to the bulk phase at t = 0 and held there (slaved,
                         instantaneous);
  lower component (ψ_R)  phase starts at Φ_bulk(0) + θ₀ and follows the DK repo's Kuramoto
                         form  dφ_R/dt = ω + K sin(Φ_bulk − φ_R).

Recorded: the L–R relative phase φ_R − φ_L over two rotations of the clock, and the rotation
at which the lower component enters the lock band |φ_R − Φ_bulk| ≤ ε.  The RK45 lock time
is compared with the closed form of kuramoto_rotations.py.

Limitation L1 (stated in the pre-registration): the closed Dirac mass coupling between the
two components, which the DK record found to be a cosine with no attractor, is switched off.
The two components are treated as independently coupled to the bulk.
"""
import csv, sys
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
from kuramoto_rotations import time_to_band, wrap, T_ROT, EPS_PRIMARY, SURFACE, TEXT, TEXT2, GRID, SLOT, OUT

OMEGA = 1.0
X_LIST = [0.3, 0.7]
THETA0_DEG = [30, 90, 150, 175]

def lower_component(theta0, K, omega_bulk, t_end):
    """φ_R(t) in the lab frame, and the bulk phase; upper component φ_L(t) = Φ_bulk(t)."""
    rhs = lambda t, y: [OMEGA + K * np.sin(omega_bulk * t - y[0])]
    t_eval = np.linspace(0, t_end, 4001)
    sol = solve_ivp(rhs, (0, t_end), [theta0], t_eval=t_eval, rtol=1e-11, atol=1e-13, max_step=0.02)
    phi_R = sol.y[0]
    phi_L = omega_bulk * sol.t                      # slaved from the first instant
    return sol.t, phi_L, phi_R

def main():
    rows = []
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.0), sharey=True)
    for ax, x in zip(axes, X_LIST):
        K = x * OMEGA
        ax.axhspan(-np.degrees(EPS_PRIMARY), np.degrees(EPS_PRIMARY), color=TEXT2, alpha=0.08, lw=0)
        ax.text(1.98, np.degrees(EPS_PRIMARY) + 3, "lock band ±ε", ha="right", color=TEXT2, fontsize=7.5)
        ax.axhline(0, color=TEXT2, lw=0.8, alpha=0.6)
        ax.text(0.02, 4, "upper component: slaved at t = 0", color=TEXT2, fontsize=7.5)
        for th_deg, c in zip(THETA0_DEG, SLOT):
            th0 = np.radians(th_deg)
            t, phi_L, phi_R = lower_component(th0, K, OMEGA, 2 * T_ROT)
            rel = np.degrees(wrap(phi_R - phi_L))
            ax.plot(t / T_ROT, rel, color=c, lw=1.6, solid_capstyle="round", label=f"θ₀ = {th_deg}°")
            ax.annotate(f"θ₀ = {th_deg}°", (t[0] / T_ROT, rel[0]), xytext=(6, 5), textcoords="offset points",
                        ha="left", va="bottom", color=TEXT, fontsize=8)
            # lock time: first entry into the band (RK45 trajectory) and closed form
            inside = np.abs(np.radians(rel)) <= EPS_PRIMARY
            t_rk = t[np.argmax(inside)] if inside.any() else np.inf
            t_cf = float(time_to_band(th0, K, 0.0, EPS_PRIMARY))
            rows.append([x, th_deg, f"{t_rk/T_ROT:.4f}" if np.isfinite(t_rk) else "> 2 (outside the window)",
                         f"{t_cf/T_ROT:.4f}" if np.isfinite(t_cf) else "inf",
                         int(t_cf <= T_ROT), int(t_cf <= 2 * T_ROT)])
        for n in (1, 2):
            ax.axvline(n, color=TEXT2, lw=0.8, alpha=0.5)
        ax.set_xlim(-0.02, 2.0); ax.set_ylim(-20, 190); ax.set_yticks([0, 45, 90, 135, 180])
        ax.set_xlabel("time (rotations of the clock)")
        ax.set_title(f"L–R relative phase, x = K/ω = {x}, Δ = 0")
    axes[0].set_ylabel("φ_R − φ_L (degrees)")
    axes[0].legend(loc="upper right")
    fig.tight_layout()
    fig.savefig(OUT / "fig4_weyl_split_toy.png", dpi=150)
    plt.close(fig)
    with open(OUT / "weyl_split_lock_times.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["x_K_over_omega", "theta0_deg", "lower_lock_rotations_rk45", "lower_lock_rotations_closed_form",
                    "within_1_rotation", "within_2_rotations"])
        w.writerows(rows)
    print("Weyl-split toy (upper component locked at t = 0 by construction); lower component lock time, ε = 0.35 rad")
    print("   x    θ₀     RK45 (rot, 2-rotation window)  closed form (rot)  ≤1  ≤2")
    for r in rows:
        print(f"  {r[0]:.1f}  {r[1]:4d}°   {r[2]:>26}   {r[3]:>9}          {r[4]}   {r[5]}")

if __name__ == "__main__":
    main()
