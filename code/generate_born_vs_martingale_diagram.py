#!/usr/bin/env python3
"""
Generates a static publication-quality comparison diagram
(fig_v2_7_born_vs_martingale_comparison.png) oriented with the Arrow of Time
pointing UPWARD (Bottom -> Top):

  - LEFT FRAME: Standard Analytical Calculation
      (Bottom) Feynman multipath integral propagation A(x) = sum_gamma a_gamma exp(i S[gamma]/hbar)
      (Middle) Direct application of the postulated Born Rule P(x) = |A(x)|^2 / integral |A|^2 dx
      (Top)    Resultant continuous analytical Born distribution curve at the top.

  - RIGHT FRAME: Physical Field-Matter Selector & Conserving Martingale Mechanism
      (Bottom) Wave packet A_m(z_0^-) plotted as it reaches a short distance just underneath
               the detector surface.
      (Lower-Middle Box) Pre-modulus 3-way interference calculation between:
               (1) Quantum vacuum phase at the detector surface phi_m^vac,
               (2) Phase distribution of bound electrons at the detector surface phi_m^det, and
               (3) Incident wave packet A_m = |A_m| exp(i phi_m^inc).
      (Upper-Middle Surface + Adjacent Martingale Box)
               - Surface Box: Sub-threshold atomic excitations s_m(t) < theta and single-atom
                 threshold registration s_{m*}(tau_theta) >= theta.
               - Adjacent Box Next to Surface: Martingale energy-conserving exchange game
                 Delta s_m = 2 eps sqrt(s_m s_n) cos(Delta Phi_mn), sum_m Delta s_m = 0,
                 Var(Delta s_m) = 2 eps^2 s_m(1 - s_m), producing the stochastic spread of
                 individual fixation events.
      (Top)    Accumulated discrete single-particle detection clicks approximating the Born distribution.
"""

from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

CODE_DIR = Path(__file__).resolve().parent


def draw_box(ax, x, y, w, h, facecolor="#1e293b", edgecolor="#475569", lw=1.6, radius=0.12):
    patch = FancyBboxPatch(
        (x, y), w, h,
        boxstyle=f"round,pad=0.02,rounding_size={radius}",
        facecolor=facecolor,
        edgecolor=edgecolor,
        linewidth=lw,
        zorder=2,
    )
    ax.add_patch(patch)
    return patch


def draw_up_arrow(ax, x, y0, y1, color="#38bdf8", lw=2.2):
    ax.annotate(
        "",
        xy=(x, y1),
        xytext=(x, y0),
        arrowprops=dict(arrowstyle="-|>,head_length=0.65,head_width=0.45", color=color, lw=lw),
        zorder=5,
    )


def build_comparison_diagram():
    fig, ax = plt.subplots(figsize=(14.2, 9.6), dpi=160)
    fig.patch.set_facecolor("#0f172a")
    ax.set_facecolor("#0f172a")
    ax.set_xlim(0.0, 14.2)
    ax.set_ylim(0.0, 10.0)
    ax.axis("off")

    # Main Title
    ax.text(
        7.1, 9.72,
        "Comparison of Ensemble Born-Rule Calculation vs. Pre-Modulus Surface Coupling & Martingale Fixation",
        color="#f8fafc", fontsize=13.5, fontweight="bold", ha="center", va="center",
    )

    # --- FAR LEFT MARGIN: UPWARD ARROW OF TIME (t ↑) ---
    draw_up_arrow(ax, 0.42, 0.55, 9.25, color="#38bdf8", lw=2.8)
    ax.text(
        0.18, 4.9,
        "ARROW OF TIME  (t ↑)",
        color="#38bdf8", fontsize=11, fontweight="bold", rotation=90, ha="center", va="center",
    )

    # Frame Headers
    ax.text(
        3.55, 9.32,
        "LEFT: Analytical Path Integral + Born Postulate",
        color="#38bdf8", fontsize=11.5, fontweight="bold", ha="center",
    )
    ax.text(
        10.25, 9.32,
        "RIGHT: Microscopic Surface-Coupled Martingale Selector",
        color="#10b981", fontsize=11.5, fontweight="bold", ha="center",
    )

    # Center Vertical Divider
    ax.plot([6.45, 6.45], [0.35, 9.45], color="#334155", linestyle="--", linewidth=1.5)

    # Common x-grid for wave / Born profiles
    u = np.linspace(-1.0, 1.0, 240)
    beta = u * 1.15 * np.pi
    alpha = 3.5 * beta
    born_env = (np.sinc(beta / np.pi) ** 2) * (np.cos(alpha) ** 2) + 0.015 * (np.sinc(beta / np.pi) ** 2)
    born_env /= np.max(born_env)
    amp_env = np.sqrt(born_env)

    # =========================================================================
    # LEFT FRAME (x = 0.85 .. 6.15): FEYNMAN MULTIPATH INTEGRAL + BORN POSTULATE
    # =========================================================================

    # 1. Left-Bottom Box (y = 0.45 .. 3.15): Feynman Multipath Integral Calculation
    draw_box(ax, 0.85, 0.45, 5.30, 2.70, facecolor="#1e293b", edgecolor="#38bdf8", lw=1.7)
    ax.text(
        3.50, 2.88,
        "1. Feynman Multipath Propagation Integral",
        color="#f8fafc", fontsize=10.5, fontweight="bold", ha="center",
    )

    # Schematic of source -> two slits -> multiple interfering paths gamma_1, gamma_2, ...
    src_x, src_y = 3.50, 0.72
    slit_y = 1.25
    det_y = 1.92
    s1_x, s2_x = 2.65, 4.35
    ax.plot([1.15, s1_x - 0.22], [slit_y, slit_y], color="#94a3b8", lw=3.0, zorder=3)
    ax.plot([s1_x + 0.22, s2_x - 0.22], [slit_y, slit_y], color="#94a3b8", lw=3.0, zorder=3)
    ax.plot([s2_x + 0.22, 5.85], [slit_y, slit_y], color="#94a3b8", lw=3.0, zorder=3)
    ax.scatter([src_x], [src_y], color="#38bdf8", s=38, zorder=4)
    ax.text(src_x + 0.18, src_y - 0.05, "Source", color="#cbd5e1", fontsize=8.5)

    for tx in (2.05, 2.80, 3.50, 4.20, 4.95):
        ax.plot([src_x, s1_x, tx], [src_y, slit_y, det_y], color="#38bdf8", alpha=0.45, lw=1.1, zorder=3)
        ax.plot([src_x, s2_x, tx], [src_y, slit_y, det_y], color="#c084fc", alpha=0.45, lw=1.1, linestyle="--", zorder=3)

    ax.text(
        3.50, 2.46,
        r"$A(x) = \sum_{\mathrm{paths}\;\gamma} a_\gamma\,\exp\!\left(\frac{i}{\hbar}S[\gamma]\right) = |A(x)|\,e^{i\phi^{\mathrm{inc}}(x)}$",
        color="#e2e8f0", fontsize=10.2, ha="center", zorder=4,
    )
    ax.text(
        3.50, 2.08,
        "Deterministic complex wave amplitude across trajectories",
        color="#94a3b8", fontsize=8.6, ha="center", style="italic",
    )

    draw_up_arrow(ax, 3.50, 3.18, 3.72, color="#38bdf8")

    # 2. Left-Middle Box (y = 3.75 .. 5.95): Direct Application of the Born Rule
    draw_box(ax, 0.85, 3.75, 5.30, 2.20, facecolor="#1e293b", edgecolor="#38bdf8", lw=1.7)
    ax.text(
        3.50, 5.65,
        "2. Postulated Born Rule Applied Directly to A(x)",
        color="#f8fafc", fontsize=10.5, fontweight="bold", ha="center",
    )
    ax.text(
        3.50, 5.05,
        r"$P(x) = \frac{|A(x)|^2}{\int |A(x')|^2\,dx'}$",
        color="#38bdf8", fontsize=12.0, fontweight="bold", ha="center",
    )
    ax.text(
        3.50, 4.42,
        r"Modulus squared $|A(x)|^2 = A^*(x)A(x)$ is postulated to be",
        color="#cbd5e1", fontsize=9.0, ha="center",
    )
    ax.text(
        3.50, 4.10,
        "the ensemble probability density (no single-event mechanism).",
        color="#94a3b8", fontsize=8.8, ha="center", style="italic",
    )

    draw_up_arrow(ax, 3.50, 5.98, 6.52, color="#38bdf8")

    # 3. Left-Top Box (y = 6.55 .. 9.05): Resultant Analytical Born Distribution
    draw_box(ax, 0.85, 6.55, 5.30, 2.50, facecolor="#1e293b", edgecolor="#38bdf8", lw=1.7)
    ax.text(
        3.50, 8.78,
        "3. Resultant Analytical Ensemble Distribution P(x)",
        color="#f8fafc", fontsize=10.5, fontweight="bold", ha="center",
    )
    x_left_plot = 1.20 + (u + 1.0) * 0.5 * 4.60
    y_left_base = 6.85
    y_left_curve = y_left_base + born_env * 1.62
    ax.fill_between(x_left_plot, y_left_base, y_left_curve, color="#38bdf8", alpha=0.28, zorder=3)
    ax.plot(x_left_plot, y_left_curve, color="#38bdf8", lw=2.2, zorder=4)
    ax.text(
        3.50, 6.66,
        r"Continuous Ensemble Prediction $P(x) \propto \mathrm{sinc}^2(\beta)\cos^2(\alpha)$",
        color="#cbd5e1", fontsize=8.8, ha="center",
    )

    # =========================================================================
    # RIGHT FRAME (x = 6.75 .. 13.95): PRE-MODULUS SURFACE COUPLING & MARTINGALE
    # =========================================================================

    # 1. Right-Bottom Box (y = 0.45 .. 2.35): Wave Function Just Underneath Detector Surface
    draw_box(ax, 6.75, 0.45, 7.20, 1.90, facecolor="#1e293b", edgecolor="#38bdf8", lw=1.7)
    ax.text(
        10.35, 2.10,
        r"1. Wave Packet Reaching Just Underneath Detector Surface ($z = z_{\mathrm{surf}} - \delta z$)",
        color="#f8fafc", fontsize=10.2, fontweight="bold", ha="center",
    )
    x_right_wave = 7.10 + (u + 1.0) * 0.5 * 3.55
    y_rw_base = 0.68
    osc_wave = amp_env * np.cos(14.0 * u)
    ax.fill_between(x_right_wave, y_rw_base, y_rw_base + amp_env * 1.15, color="#38bdf8", alpha=0.22, zorder=3)
    ax.plot(x_right_wave, y_rw_base + amp_env * 1.15, color="#38bdf8", lw=1.8, zorder=4, label=r"$|A_m|$")
    ax.plot(x_right_wave, y_rw_base + 0.55 + 0.42 * osc_wave, color="#7dd3fc", lw=1.0, linestyle="--", zorder=4)

    ax.text(
        12.35, 1.48,
        r"$A_m(z_0^-) = |A_m|\,e^{i\phi^{\mathrm{inc}}_m}$",
        color="#e2e8f0", fontsize=10.2, ha="center",
    )
    ax.text(
        12.35, 1.05,
        r"$\sum_{m=1}^M |A_m|^2 = 1$",
        color="#38bdf8", fontsize=10.0, ha="center",
    )
    ax.text(
        12.35, 0.68,
        "(Continuous wave below surface)",
        color="#94a3b8", fontsize=8.4, ha="center", style="italic",
    )

    draw_up_arrow(ax, 10.35, 2.38, 2.76, color="#c084fc")

    # 2. Right Lower-Middle Box (y = 2.78 .. 4.68): Three-Way Pre-Modulus Interference Box
    #    (Between Propagation Calculation and the Detector Surface Box)
    draw_box(ax, 6.75, 2.78, 7.20, 1.90, facecolor="#1e293b", edgecolor="#c084fc", lw=1.8)
    ax.text(
        10.35, 4.44,
        "2. Pre-Modulus Three-Way Phase Interference at the Detector Boundary",
        color="#f8fafc", fontsize=10.2, fontweight="bold", ha="center",
    )
    ax.text(
        10.35, 4.10,
        r"Three complex phases couple at atom $m$ BEFORE any modulus is squared:",
        color="#cbd5e1", fontsize=8.9, ha="center",
    )
    ax.text(
        10.35, 3.74,
        r"(i) Fast Uniform Vacuum $e^{i\phi^{\mathrm{vac}}_m}$  +  (ii) Slow Detector Electron $e^{i\phi^{\mathrm{det}}_m}$  +  (iii) Incident Wave $|A_m|e^{i\phi^{\mathrm{inc}}_m}$",
        color="#e9d5ff", fontsize=8.0, ha="center", fontweight="bold",
    )
    ax.text(
        10.35, 3.30,
        r"$C_m = \frac{1}{2}[\cos(\phi^{\mathrm{det}}_m - \phi^{\mathrm{inc}}_m) + \cos(\phi^{\mathrm{det}}_m - \phi^{\mathrm{vac}}_m)], \quad F_m = 1 + \alpha_{\mathrm{seed}} C_m, \quad s_m(0) = S_r\frac{|A_m|^2 F_m}{\sum_{k \in r} |A_k|^2 F_k}$",
        color="#f8fafc", fontsize=8.7, ha="center",
    )
    ax.text(
        10.35, 2.93,
        r"$2\pi$-periodic in every phase; fast uniform QED vacuum $\phi^{\mathrm{vac}}_m \sim \mathcal{U}(-\pi, \pi)$ provides the white-noise dice",
        color="#94a3b8", fontsize=8.3, ha="center", style="italic",
    )

    # Upward arrow from 3-way interference box into the Surface Box
    draw_up_arrow(ax, 8.45, 4.70, 5.08, color="#10b981")

    # 3A. Right Upper-Middle LEFT Sub-Box (x = 6.75 .. 10.15, y = 5.10 .. 7.12):
    #     Detector Surface Atomic Layer & Whole-Quantum Registration
    draw_box(ax, 6.75, 5.10, 3.40, 2.02, facecolor="#1e293b", edgecolor="#10b981", lw=1.8)
    ax.text(
        8.45, 6.88,
        "3a. Detector Surface Atoms",
        color="#f8fafc", fontsize=10.0, fontweight="bold", ha="center",
    )
    ax.text(
        8.45, 6.56,
        r"Whole-quantum absorption ($\hbar\omega \geq E_{\mathrm{gap}},\,\theta = 1$)",
        color="#cbd5e1", fontsize=8.3, ha="center",
    )
    ax.text(
        8.45, 6.22,
        r"$\tau_{\mathrm{game}} = \inf\{t \geq 0 : \max_m s_m(t) = 1\}$",
        color="#10b981", fontsize=9.1, ha="center",
    )
    ax.text(
        8.45, 5.84,
        r"$\Delta E_{\mathrm{ex}} = \hbar\omega - E_{\mathrm{gap}} \to E_{\mathrm{kin}}/\mathrm{phonons}\;\mathrm{or}\;\hbar\omega'$",
        color="#f43f5e", fontsize=8.2, ha="center", fontweight="bold",
    )
    ax.text(
        8.45, 5.48,
        r"Winning site $m^*$ absorbs full quantum $\hbar\omega$;",
        color="#cbd5e1", fontsize=8.3, ha="center",
    )
    ax.text(
        8.45, 5.23,
        r"excess energy conserved locally at $m^*$.",
        color="#94a3b8", fontsize=8.2, ha="center", style="italic",
    )

    # Bidirectional coupling arrows between Surface Box (3a) and Martingale Game Box (3b)
    ax.annotate(
        "", xy=(10.52, 6.32), xytext=(10.16, 6.32),
        arrowprops=dict(arrowstyle="<->,head_length=0.45,head_width=0.35", color="#fbbf24", lw=2.0),
        zorder=5,
    )
    ax.annotate(
        "", xy=(10.52, 5.75), xytext=(10.16, 5.75),
        arrowprops=dict(arrowstyle="<->,head_length=0.45,head_width=0.35", color="#fbbf24", lw=2.0),
        zorder=5,
    )

    # 3B. Right Upper-Middle RIGHT Sub-Box (x = 10.54 .. 13.95, y = 5.10 .. 7.12):
    #     Box Next to Surface: Martingale Energy-Conserving Game
    draw_box(ax, 10.54, 5.10, 3.41, 2.02, facecolor="#1e293b", edgecolor="#f59e0b", lw=1.8)
    ax.text(
        12.24, 6.88,
        "3b. Conserving Martingale Game",
        color="#fbbf24", fontsize=10.0, fontweight="bold", ha="center",
    )
    ax.text(
        12.24, 6.48,
        r"$\Delta s_m = -\Delta s_n = 2\epsilon\sqrt{s_m s_n}\sin\Delta\Phi_{mn}$",
        color="#f8fafc", fontsize=8.9, ha="center",
    )
    ax.text(
        12.24, 6.08,
        r"$\sum_{m=1}^M \Delta s_m = 0, \quad \mathbb{E}[\Delta s_m \mid \mathcal{F}_t] = 0$",
        color="#10b981", fontsize=8.9, ha="center",
    )
    ax.text(
        12.24, 5.66,
        r"$\mathrm{Var}(\Delta s_m \mid \mathcal{F}_t) = 2\epsilon^2 s_m(1 - s_m)$",
        color="#f8fafc", fontsize=8.9, ha="center",
    )
    ax.text(
        12.24, 5.26,
        r"Fixation: $P(\mathrm{fix}\;m^*) = s_{m^*}(0)$",
        color="#fbbf24", fontsize=9.1, fontweight="bold", ha="center",
    )

    # Upward arrows from Surface + Martingale Boxes into Top Stochastic Spread Box
    draw_up_arrow(ax, 8.45, 7.14, 7.52, color="#f59e0b")
    draw_up_arrow(ax, 12.24, 7.14, 7.52, color="#f59e0b")

    # 4. Right-Top Box (y = 7.55 .. 9.05): Stochastic Spread of Individual Events Matching Born Curve
    draw_box(ax, 6.75, 7.55, 7.20, 1.50, facecolor="#1e293b", edgecolor="#f59e0b", lw=1.8)
    ax.text(
        10.35, 8.82,
        "4. Stochastic Spread of Individual Fixation Events (1 Click / Quantum → Born Distribution)",
        color="#f8fafc", fontsize=9.1, fontweight="bold", ha="center",
    )

    # Plot Born target curve + discrete piled dots from simulated single-quantum fixation events
    x_right_top = 7.15 + (u + 1.0) * 0.5 * 6.40
    y_rt_base = 7.72
    y_rt_curve = y_rt_base + born_env * 0.92
    ax.plot(x_right_top, y_rt_curve, color="#38bdf8", lw=1.9, zorder=4, label="Born Reference Curve")

    rng = np.random.default_rng(20261029)
    n_bins = 29
    u_bins = np.linspace(-0.95, 0.95, n_bins)
    b_beta = u_bins * 1.15 * np.pi
    b_alpha = 3.5 * b_beta
    p_bins = (np.sinc(b_beta / np.pi) ** 2) * (np.cos(b_alpha) ** 2) + 0.015 * (np.sinc(b_beta / np.pi) ** 2)
    p_bins /= np.sum(p_bins)
    counts = rng.multinomial(125, p_bins)
    max_c = max(1, int(np.max(counts)))
    for b_idx, c in enumerate(counts):
        bx = 7.15 + ((u_bins[b_idx] + 1.0) * 0.5) * 6.40
        for r in range(c):
            by = y_rt_base + ((r + 0.5) / max_c) * 0.92
            ax.scatter([bx], [by], color="#f59e0b", s=16, zorder=5, edgecolors="none")

    plt.tight_layout()
    out_png = CODE_DIR / "fig_v2_7_born_vs_martingale_comparison.png"
    fig.savefig(out_png, dpi=160, facecolor=fig.get_facecolor())
    plt.close(fig)
    print(f"Saved static comparison diagram to {out_png}")


if __name__ == "__main__":
    build_comparison_diagram()
