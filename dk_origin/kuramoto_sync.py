"""
kuramoto_sync.py — Task 1: Phase Synchronization Simulation
============================================================

Models the internal "time clock" of a particle as a Kuramoto oscillator.
At creation, two entangled particles share an initial phase φ_0.
Free evolution keeps their phases locked (same ω).
Upon detection, each particle synchronizes to its detector's bulk phase Φ_bulk.

Kuramoto single-oscillator equation:
    dφ/dt = ω + K · sin(Φ_target − φ)

Free flight: K = 0   →   dφ/dt = ω  (uniform precession)
Interaction: K > 0   →   φ → Φ_target  (synchronization)

The "synchronization equation" the theory draws on is the Kuramoto model
(Yoshiki Kuramoto, 1975). Related phase-locking phenomena appear in the
Josephson junction equations and Adler's phase-lock equation.
"""

import numpy as np
from pathlib import Path
_HERE = Path(__file__).resolve().parent
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# ─── Parameters ──────────────────────────────────────────────────────────────

OMEGA       = 2.0 * np.pi        # natural frequency (rad/s), represents E/ℏ
K_COUPLE    = 15.0               # Kuramoto coupling constant (strong = fast sync)
PHI_BULK    = 0.0                # bulk detector phase (reference = 0)

T_CREATE    = 0.0                # pair created at t=0
T_DETECT_A  = 3.0                # particle A hits detector D1
T_DETECT_B  = 3.5                # particle B hits detector D2 (slightly later)
T_END       = 6.0

DT          = 0.005              # integration step


# ─── Dynamics ────────────────────────────────────────────────────────────────

def dphidt(t, phi, omega, K, Phi_target, t_detect):
    """Kuramoto equation: free until t_detect, then couples to bulk."""
    coupling = K if t >= t_detect else 0.0
    return [omega + coupling * np.sin(Phi_target - phi[0])]


def simulate_particle(phi_0, t_detect, label):
    """Integrate phase trajectory for one particle."""
    t_span = (T_CREATE, T_END)
    t_eval = np.arange(T_CREATE, T_END, DT)
    sol = solve_ivp(
        fun=lambda t, y: dphidt(t, y, OMEGA, K_COUPLE, PHI_BULK, t_detect),
        t_span=t_span,
        y0=[phi_0],
        t_eval=t_eval,
        method='RK45',
        rtol=1e-8,
    )
    return sol.t, sol.y[0]


# ─── Simulation ──────────────────────────────────────────────────────────────

def run():
    phi_0 = np.pi * 0.7          # shared initial clock phase (arbitrary)

    t_A, phi_A = simulate_particle(phi_0, T_DETECT_A, 'A')
    t_B, phi_B = simulate_particle(phi_0, T_DETECT_B, 'B')

    # Wrap phase to [−π, π] for display
    phi_A_mod = (phi_A % (2 * np.pi)) - np.pi
    phi_B_mod = (phi_B % (2 * np.pi)) - np.pi

    # ── Phase difference (should stay ~0 throughout) ──────────────────────
    delta_phi = phi_A - phi_B

    # ── Plot ──────────────────────────────────────────────────────────────
    fig, axes = plt.subplots(3, 1, figsize=(10, 9), sharex=True)
    fig.suptitle('Kuramoto Phase Synchronization of Entangled Particles', fontsize=14)

    # Panel 1: Phase trajectories
    ax = axes[0]
    ax.plot(t_A, phi_A_mod, label='Particle A (clock phase mod 2π)', color='steelblue')
    ax.plot(t_B, phi_B_mod, label='Particle B (clock phase mod 2π)', color='tomato', linestyle='--')
    ax.axhline(PHI_BULK, color='gray', linestyle=':', linewidth=1.5, label='Φ_bulk (detector reference)')
    ax.axvline(T_DETECT_A, color='steelblue', linestyle=':', alpha=0.6, label=f'A hits D1 (t={T_DETECT_A})')
    ax.axvline(T_DETECT_B, color='tomato', linestyle=':', alpha=0.6, label=f'B hits D2 (t={T_DETECT_B})')
    ax.set_ylabel('φ (rad)')
    ax.legend(fontsize=8)
    ax.set_title('Clock phases: free precession → synchronization to bulk')

    # Panel 2: Phase difference between A and B
    ax = axes[1]
    ax.plot(t_A, delta_phi % (2 * np.pi), color='purple')
    ax.axhline(0, color='gray', linestyle=':', linewidth=1)
    ax.axvline(T_DETECT_A, color='steelblue', linestyle=':', alpha=0.6)
    ax.axvline(T_DETECT_B, color='tomato', linestyle=':', alpha=0.6)
    ax.set_ylabel('φ_A − φ_B (rad)')
    ax.set_title('Phase difference: near-zero throughout (clocks stay synchronized)')

    # Panel 3: Unit-circle view at key times (snapshots)
    ax = axes[2]
    times_of_interest = {
        'Creation':    0.0,
        'Pre-detect':  T_DETECT_A - 0.1,
        'Post-A sync': T_DETECT_A + 0.5,
        'Post-B sync': T_DETECT_B + 0.5,
    }
    theta = np.linspace(0, 2 * np.pi, 300)
    ax.plot(np.cos(theta), np.sin(theta), 'k-', linewidth=0.5, alpha=0.3)
    colors = ['navy', 'steelblue', 'darkorange', 'tomato']
    for (name, t_snap), col in zip(times_of_interest.items(), colors):
        idx_A = np.searchsorted(t_A, t_snap)
        idx_B = np.searchsorted(t_B, t_snap)
        pA = phi_A[min(idx_A, len(phi_A)-1)]
        pB = phi_B[min(idx_B, len(phi_B)-1)]
        ax.annotate('', xy=(np.cos(pA)*0.9, np.sin(pA)*0.9), xytext=(0, 0),
                    arrowprops=dict(arrowstyle='->', color=col, lw=2))
        ax.annotate('', xy=(np.cos(pB)*0.7, np.sin(pB)*0.7), xytext=(0, 0),
                    arrowprops=dict(arrowstyle='->', color=col, lw=1.5, linestyle='dashed'))
        ax.text(np.cos(pA)*1.05, np.sin(pA)*1.05, f'A:{name[:4]}', fontsize=7, color=col)
    ax.set_aspect('equal')
    ax.set_xlim(-1.3, 1.3)
    ax.set_ylim(-1.3, 1.3)
    ax.set_title('Phase vectors on unit circle (solid=A, dashed=B)')
    ax.set_xlabel('Re(e^{iφ})')
    ax.set_ylabel('Im(e^{iφ})')

    plt.tight_layout()
    plt.savefig(_HERE / 'kuramoto_sync.png', dpi=150)
    print("Saved: kuramoto_sync.png")
    plt.show()

    # ── Console report ────────────────────────────────────────────────────
    print(f"\nφ_0 (shared initial phase):         {phi_0:.4f} rad")
    print(f"φ_bulk (detector reference):        {PHI_BULK:.4f} rad")

    for label, t_arr, phi_arr, t_det in [('A', t_A, phi_A, T_DETECT_A),
                                          ('B', t_B, phi_B, T_DETECT_B)]:
        i_pre  = np.searchsorted(t_arr, t_det) - 1
        i_post = np.searchsorted(t_arr, t_det + 0.8)
        print(f"\nParticle {label}:")
        print(f"  φ just before detection:  {phi_arr[i_pre] % (2*np.pi):.4f} rad")
        print(f"  φ after synchronization:  {phi_arr[i_post] % (2*np.pi):.4f} rad  "
              f"(target Φ_bulk = {PHI_BULK:.4f})")


if __name__ == '__main__':
    run()
