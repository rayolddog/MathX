import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.patches as patches

plt.rcParams['font.sans-serif'] = 'DejaVu Sans'
plt.rcParams['axes.edgecolor'] = '#334155'
plt.rcParams['axes.linewidth'] = 1.0

fig = plt.figure(figsize=(16.5, 12.6), dpi=220, facecolor='#F8FAFC')
gs = fig.add_gridspec(2, 2, height_ratios=[1.12, 1.0], width_ratios=[1.05, 1.0],
                      hspace=0.28, wspace=0.20, left=0.05, right=0.97, top=0.89, bottom=0.06)

fig.suptitle(
    "Translating a QED Feynman Diagram into Equations & Tracing Wave Peaks (↑) and Troughs (↓)",
    fontsize=15.0, fontweight='bold', color='#0F172A', y=0.968
)
fig.text(
    0.5, 0.935,
    "Every piece of a Feynman diagram (external legs, vertices x₁ & x₂, photon propagator) is a shorthand for multiplying and summing wave peaks and troughs in 3D + time",
    ha='center', va='center', fontsize=10.6, color='#334155', style='italic'
)

# ==============================================================================
# PANEL A (Top-Left): The QED Feynman Diagram & Equation Translation Dictionary
# ==============================================================================
ax_a = fig.add_subplot(gs[0, 0])
ax_a.set_facecolor('#0B1120')

# Coordinates for Vertex 1 (left/bottom) and Vertex 2 (right/top)
v1_x, v1_y = 3.2, 2.35
v2_x, v2_y = 6.8, 2.35

# Helper to draw a fermion line with wave ribbon + arrow
def draw_fermion_leg_with_wave(ax, x0, y0, x1, y1, col, wave_col, label_str, lbl_xy, n_cycles=3.5):
    dx, dy = x1 - x0, y1 - y0
    L = np.hypot(dx, dy)
    ux, uy = dx / L, dy / L
    nx, ny = -uy, ux

    # Main fermion line
    ax.plot([x0, x1], [y0, y1], color=col, lw=2.6, zorder=6)
    # Midpoint direction arrow
    xm, ym = x0 + 0.52 * dx, y0 + 0.52 * dy
    ax.annotate("", xy=(xm + 0.15 * ux, ym + 0.15 * uy), xytext=(xm - 0.15 * ux, ym - 0.15 * uy),
                arrowprops=dict(arrowstyle='->,head_width=0.32,head_length=0.28', color=col, lw=2.6), zorder=8)

    # Superimposed sinusoidal wave of peaks and troughs along the external leg
    s = np.linspace(0.2, L - 0.2, 300)
    w = 0.22 * np.cos(2 * np.pi * n_cycles * s / L)
    wx = x0 + s * ux + w * nx
    wy = y0 + s * uy + w * ny
    ax.plot(wx, wy, color=wave_col, lw=1.5, alpha=0.82, zorder=7)

    # Draw 2 Peak arrows (↑) and 1 Trough arrow (↓) along the leg
    for m_c in range(1, int(n_cycles)):
        s_pk = (m_c / n_cycles) * L
        if 0.4 < s_pk < L - 0.4:
            px, py = x0 + s_pk * ux, y0 + s_pk * uy
            ax.scatter([px], [py], s=26, color='#22D3EE', zorder=9)
            ax.annotate("", xy=(px + 0.30 * nx, py + 0.30 * ny), xytext=(px, py),
                        arrowprops=dict(arrowstyle='->,head_width=0.2,head_length=0.18', color='#22D3EE', lw=1.8), zorder=10)
        s_tr = ((m_c + 0.5) / n_cycles) * L
        if 0.4 < s_tr < L - 0.4:
            tx, ty = x0 + s_tr * ux, y0 + s_tr * uy
            ax.scatter([tx], [ty], s=26, color='#FB7185', zorder=9)
            ax.annotate("", xy=(tx - 0.30 * nx, ty - 0.30 * ny), xytext=(tx, ty),
                        arrowprops=dict(arrowstyle='->,head_width=0.2,head_length=0.18', color='#FB7185', lw=1.8), zorder=10)

    ax.text(lbl_xy[0], lbl_xy[1], label_str, color='white', fontsize=8.4, fontweight='bold',
            ha='center', va='center', zorder=12,
            bbox=dict(boxstyle='round,pad=0.24', facecolor='#1E293B', edgecolor=col, lw=1.2, alpha=0.94))

# Incoming Electron 1 (bottom-left -> Vertex 1)
draw_fermion_leg_with_wave(
    ax_a, 0.5, 0.85, v1_x, v1_y, col='#38BDF8', wave_col='#7DD3FC',
    label_str="1. Incoming e⁻(p₁)\nu(p₁) e^{i(k₁·r₁ − ω₁t₁)}",
    lbl_xy=(1.35, 1.15), n_cycles=3.0
)

# Outgoing Electron 3 (Vertex 1 -> top-left)
draw_fermion_leg_with_wave(
    ax_a, v1_x, v1_y, 0.5, 4.15, col='#38BDF8', wave_col='#7DD3FC',
    label_str="2. Outgoing e⁻(p₃)\nū(p₃) e^{−i(k₃·r₁ − ω₃t₁)}",
    lbl_xy=(1.45, 4.05), n_cycles=3.0
)

# Incoming Electron 2 (bottom-right -> Vertex 2)
draw_fermion_leg_with_wave(
    ax_a, 9.5, 0.85, v2_x, v2_y, col='#FBBF24', wave_col='#FDE047',
    label_str="Incoming e⁻(p₂)\nu(p₂) e^{i(k₂·r₂ − ω₂t₂)}",
    lbl_xy=(8.65, 1.15), n_cycles=3.0
)

# Outgoing Electron 4 (Vertex 2 -> top-right)
draw_fermion_leg_with_wave(
    ax_a, v2_x, v2_y, 9.5, 4.15, col='#FBBF24', wave_col='#FDE047',
    label_str="Outgoing e⁻(p₄)\nū(p₄) e^{−i(k₄·r₂ − ω₄t₂)}",
    lbl_xy=(8.55, 4.05), n_cycles=3.0
)

# Internal Wavy Photon Propagator between Vertex 1 (v1_x, v1_y) and Vertex 2 (v2_x, v2_y)
xp = np.linspace(v1_x, v2_x, 400)
yp = v1_y + 0.25 * np.sin(2 * np.pi * 4.5 * (xp - v1_x) / (v2_x - v1_x))
ax_a.plot(xp, yp, color='#4ADE80', lw=2.8, zorder=8)

# Place Peak (↑) and Trough (↓) arrows on the internal photon propagator
for m_p in range(4):
    x_pk = v1_x + (m_p + 0.25) / 4.5 * (v2_x - v1_x)
    x_tr = v1_x + (m_p + 0.75) / 4.5 * (v2_x - v1_x)
    ax_a.annotate("", xy=(x_pk, v1_y + 0.52), xytext=(x_pk, v1_y),
                  arrowprops=dict(arrowstyle='->,head_width=0.22,head_length=0.18', color='#4ADE80', lw=2.0), zorder=10)
    ax_a.annotate("", xy=(x_tr, v1_y - 0.52), xytext=(x_tr, v1_y),
                  arrowprops=dict(arrowstyle='->,head_width=0.22,head_length=0.18', color='#FB7185', lw=2.0), zorder=10)

# Draw Vertices x1 and x2
ax_a.scatter([v1_x, v2_x], [v1_y, v2_y], s=170, color='#F43F5E', edgecolors='white', lw=2.2, zorder=12)
ax_a.text(v1_x + 0.15, v1_y - 0.72, "Vertex x₁ = (r₁, t₁)\nFactor: −ieγ^μ / ħc\n(90° Shift −i = e^{−iπ/2})",
          color='#FDA4AF', fontsize=8.0, fontweight='bold', ha='center', va='top',
          bbox=dict(boxstyle='round,pad=0.20', facecolor='#4C0519', edgecolor='#F43F5E', lw=1.2))
ax_a.text(v2_x - 0.15, v2_y - 0.72, "Vertex x₂ = (r₂, t₂)\nFactor: −ieγ^ν / ħc\n∫ d⁴x₂ Beat Overlap",
          color='#FDA4AF', fontsize=8.0, fontweight='bold', ha='center', va='top',
          bbox=dict(boxstyle='round,pad=0.20', facecolor='#4C0519', edgecolor='#F43F5E', lw=1.2))

# Annotate Photon Propagator in center
ax_a.text(5.0, 3.62,
          "3. Photon Propagator D_F(x₂ − x₁) ↔ −ig_μν / (q² + iε)\nSpherical Wave Relay of Peaks (↑) & Troughs (↓) from x₁ to x₂",
          color='#A7F3D0', fontsize=8.1, fontweight='bold', ha='center', va='center',
          bbox=dict(boxstyle='round,pad=0.22', facecolor='#064E3B', edgecolor='#4ADE80', lw=1.3))

# Complete Translated Equation Banner at bottom of Panel A
ax_a.text(5.0, 0.30,
          "Translated Eq:  ℳ₁ = (−ie)² [ū(p₃)γ^μ u(p₁)] · [−ig_μν / (p₁−p₃)²] · [ū(p₄)γ^ν u(p₂)]",
          color='#F8FAFC', fontsize=8.2, fontweight='bold', ha='center', va='center',
          bbox=dict(boxstyle='round,pad=0.24', facecolor='#0F172A', edgecolor='#38BDF8', lw=1.3))

ax_a.set_xlim(0.0, 10.0)
ax_a.set_ylim(-0.05, 4.55)
ax_a.set_xticks([])
ax_a.set_yticks([])
ax_a.set_title("A. QED Feynman Diagram Annotated with Wave Peaks (↑), Troughs (↓) & Translated Equations",
               fontsize=11.0, fontweight='bold', color='#0F172A', pad=8)

# ==============================================================================
# PANEL B (Top-Right): 2D Spatial Map Around Vertex x₁ (Matched vs. Mismatched Waves)
# ==============================================================================
ax_b = fig.add_subplot(gs[0, 1])
ax_b.set_facecolor('#0B1120')

# Left half of Panel B (x in [-3, 0]): Matched wavelengths k_1 + q - k_3 = 0 -> Uniform Peak (+1) everywhere!
# Right half of Panel B (x in [0, +3]): Mismatched wavelengths Delta k != 0 -> Alternating Peak (+1) and Trough (-1) stripes!
xb = np.linspace(-3.0, 3.0, 400)
yb = np.linspace(-2.0, 2.0, 260)
XXb, YYb = np.meshgrid(xb, yb)

# Beat product wave Re[exp(i Delta_k . r_1)]:
# Left side (x < 0): Delta_k = 0 -> cos(0) = +0.85 (constructive peak plateau)
# Right side (x > 0): Delta_k_x = 2*pi / 1.0 -> cos(2*pi*x) alternating peaks (+1) and troughs (-1)
beat_field = np.where(XXb < 0.0, 0.82 * np.ones_like(XXb), np.cos(2 * np.pi * (XXb - 0.25) / 1.0))

ax_b.imshow(beat_field, extent=[-3.0, 3.0, -2.0, 2.0], origin='lower',
            cmap='RdBu_r', vmin=-1.0, vmax=1.0, alpha=0.55, aspect='auto')

ax_b.axvline(0.0, color='white', lw=2.5, ls='--')

# Draw grid of Feynman arrows on the Left Half (Matched Momentum: Δk = k₁ + q − k₃ = 0)
for gx in [-2.4, -1.6, -0.8]:
    for gy in [-1.1, -0.2, 0.7]:
        ax_b.scatter([gx], [gy], s=42, color='#4ADE80', edgecolors='#0F172A', lw=0.8, zorder=9)
        ax_b.annotate("", xy=(gx, gy + 0.42), xytext=(gx, gy),
                      arrowprops=dict(arrowstyle='->,head_width=0.25,head_length=0.22', color='#4ADE80', lw=2.3), zorder=10)

ax_b.text(-1.5, 1.58,
          "ON-SHELL / MOMENTUM CONSERVED\nk₁ + q = k₃  (Δk = 0)\nPeaks Match Peaks Everywhere (↑↑↑)\n∫ d³r₁ e^{iΔk·r₁} = V_overlap (HUGE!)",
          color='#A7F3D0', fontsize=8.5, fontweight='bold', ha='center', va='center',
          bbox=dict(boxstyle='round,pad=0.25', facecolor='#064E3B', edgecolor='#4ADE80', lw=1.3))

# Draw grid of Feynman arrows on the Right Half (Mismatched Momentum: Δk = k₁ + q − k₃ ≠ 0)
# Peaks at x = 0.25, 1.25, 2.25 (UP arrows ↑); Troughs at x = 0.75, 1.75, 2.75 (DOWN arrows ↓)
for gx_pk in [0.35, 1.35, 2.35]:
    for gy in [-1.1, -0.2, 0.7]:
        ax_b.scatter([gx_pk], [gy], s=38, color='#38BDF8', edgecolors='#0F172A', lw=0.8, zorder=9)
        ax_b.annotate("", xy=(gx_pk, gy + 0.40), xytext=(gx_pk, gy),
                      arrowprops=dict(arrowstyle='->,head_width=0.23,head_length=0.20', color='#38BDF8', lw=2.1), zorder=10)

for gx_tr in [0.85, 1.85, 2.80]:
    for gy in [-1.1, -0.2, 0.7]:
        ax_b.scatter([gx_tr], [gy], s=38, color='#FB7185', edgecolors='#0F172A', lw=0.8, zorder=9)
        ax_b.annotate("", xy=(gx_tr, gy - 0.40), xytext=(gx_tr, gy),
                      arrowprops=dict(arrowstyle='->,head_width=0.23,head_length=0.20', color='#FB7185', lw=2.1), zorder=10)

ax_b.text(1.5, 1.58,
          "MISMATCHED WAVELENGTHS\nk₁ + q ≠ k₃  (Δk ≠ 0)\nPeak Stripe (↑) Next to Trough Stripe (↓)\n∫ d³r₁ e^{iΔk·r₁} = 0 (CANCELS!)",
          color='#FECDD3', fontsize=8.5, fontweight='bold', ha='center', va='center',
          bbox=dict(boxstyle='round,pad=0.25', facecolor='#4C0519', edgecolor='#F43F5E', lw=1.3))

ax_b.set_xlim(-3.0, 3.0)
ax_b.set_ylim(-2.0, 2.0)
ax_b.set_xlabel("Vertex Position r₁ Across Overlap Region (in beat wavelengths)", fontsize=9.8, fontweight='bold')
ax_b.set_ylabel("Transverse Vertex Coordinate y₁", fontsize=9.8, fontweight='bold')
ax_b.set_title("B. Why a Feynman Vertex Conserves Momentum: Peak Alignment vs. Peak-Trough Cancellation",
               fontsize=10.6, fontweight='bold', color='#0F172A', pad=8)

# ==============================================================================
# PANEL C (Bottom-Left): Tracing the 3 Waves Across Vertex x₁ (1D Wave & Arrow Profile)
# ==============================================================================
ax_c = fig.add_subplot(gs[1, 0])
ax_c.set_facecolor('#FFFFFF')

x_v = np.linspace(0.0, 3.0, 600)
# Incoming wave (k1 = 1 cycle/unit) + Photon wave (q = 1 cycle/unit) = 2 cycles/unit
w_in_plus_photon = np.cos(2 * np.pi * 2.0 * x_v)
# Outgoing wave Matched (k3 = 2 cycles/unit) vs Mismatched (k3' = 1 cycle/unit)
w_out_matched = np.cos(2 * np.pi * 2.0 * x_v)
w_beat_matched = 0.5 * (1.0 + np.cos(4 * np.pi * 2.0 * x_v))  # positive everywhere, average = +1/2!
w_beat_mismatched = np.cos(2 * np.pi * 1.0 * x_v)             # oscillates +1 (Peak ↑) to -1 (Trough ↓), average = 0!

y_c1 = 1.65
y_c2 = 0.0
y_c3 = -1.75

ax_c.axhline(y_c1, color='#CBD5E1', lw=1.0, ls=':')
ax_c.axhline(y_c2, color='#CBD5E1', lw=1.0, ls=':')
ax_c.axhline(y_c3, color='#CBD5E1', lw=1.0, ls=':')

ax_c.plot(x_v, y_c1 + 0.50 * w_in_plus_photon, color='#0284C7', lw=2.0)
ax_c.plot(x_v, y_c2 + 0.55 * np.ones_like(x_v), color='#16A34A', lw=2.6)
ax_c.fill_between(x_v, y_c2, y_c2 + 0.55, color='#22C55E', alpha=0.22)
ax_c.plot(x_v, y_c3 + 0.55 * w_beat_mismatched, color='#DC2626', lw=2.3)
ax_c.fill_between(x_v, y_c3, y_c3 + 0.55 * np.maximum(w_beat_mismatched, 0), color='#38BDF8', alpha=0.25)
ax_c.fill_between(x_v, y_c3, y_c3 + 0.55 * np.minimum(w_beat_mismatched, 0), color='#F43F5E', alpha=0.25)

# Draw arrows along the Matched Beat Line (all UP ↑) and Mismatched Beat Line (alternating UP ↑ and DOWN ↓)
for xv_pt in [0.5, 1.0, 1.5, 2.0, 2.5]:
    # Matched line (y_c2): all arrows point UP (+0.52)
    ax_c.annotate("", xy=(xv_pt, y_c2 + 0.55), xytext=(xv_pt, y_c2),
                  arrowprops=dict(arrowstyle='->,head_width=0.22,head_length=0.18', color='#15803D', lw=2.3))
    # Mismatched line (y_c3): alternates Peak (UP) at 1.0, 2.0 and Trough (DOWN) at 0.5, 1.5, 2.5
    if np.isclose(xv_pt % 1.0, 0.0):
        ax_c.annotate("", xy=(xv_pt, y_c3 + 0.55), xytext=(xv_pt, y_c3),
                      arrowprops=dict(arrowstyle='->,head_width=0.22,head_length=0.18', color='#0284C7', lw=2.2))
    else:
        ax_c.annotate("", xy=(xv_pt, y_c3 - 0.55), xytext=(xv_pt, y_c3),
                      arrowprops=dict(arrowstyle='->,head_width=0.22,head_length=0.18', color='#DB2777', lw=2.2))

ax_c.text(0.05, y_c1 + 0.62, "1. Incoming + Photon Wave at Vertex x₁: exp[i(k₁ + q)·r₁]",
          color='#0284C7', fontsize=8.8, fontweight='bold')
ax_c.text(0.05, y_c2 + 0.66, "2. Matched Outgoing Wave (k₃ = k₁ + q): Beat Phase ΔΦ = 0° → All Peak Arrows (↑↑↑) Reinforce!",
          color='#15803D', fontsize=8.8, fontweight='bold')
ax_c.text(0.05, y_c3 + 0.66, "3. Mismatched Outgoing Wave (k₃ ≠ k₁ + q): Alternates Peak (↑) & Trough (↓) → Integrates to 0!",
          color='#DC2626', fontsize=8.8, fontweight='bold')

ax_c.set_xlim(0.0, 3.05)
ax_c.set_ylim(-2.65, 2.50)
ax_c.set_yticks([y_c3, y_c2, y_c1])
ax_c.set_yticklabels(["Mismatched\nBeat (Δk≠0)", "Matched\nBeat (Δk=0)", "Incoming\nWave (k₁+q)"], fontsize=8.8, fontweight='bold')
ax_c.set_xlabel("Vertex Position r₁ Along Interaction Region (in beat wavelengths)", fontsize=9.8, fontweight='bold')
ax_c.set_title("C. 1D Wave & Arrow Profile at Vertex x₁: Why ∫ d⁴x₁ Gives (2π)⁴ δ⁴(p₁ + q − p₃)",
               fontsize=10.6, fontweight='bold', color='#0F172A', pad=8)

# ==============================================================================
# PANEL D (Bottom-Right): Interference of Two Competing Feynman Diagrams (ℳ₁ + ℳ₂)
# ==============================================================================
ax_d = fig.add_subplot(gs[1, 1])
ax_d.set_facecolor('#FFFFFF')

# Center-of-mass scattering angle theta_cm in degrees (from 25 deg to 155 deg)
theta_deg = np.linspace(28, 152, 500)
theta_rad = np.deg2rad(theta_deg)

# In non-relativistic / low-energy identical-electron (Møller / Mott) scattering:
# Direct (t-channel) amplitude: A_1(theta) ~ 1 / sin^2(theta/2)
# Exchange (u-channel) amplitude: A_2(theta) ~ 1 / cos^2(theta/2)
A1 = 1.0 / (np.sin(theta_rad / 2.0)**2)
A2 = 1.0 / (np.cos(theta_rad / 2.0)**2)

# Normalize at theta = 90 deg so A1(90) = A2(90) = 0.5:
norm_90 = 1.0 / (np.sin(np.pi / 4.0)**2) * 2.0
A1 = A1 / norm_90
A2 = A2 / norm_90

# Classical / Distinguishable sum of probabilities: |M_1|^2 + |M_2|^2
I_incoherent = A1**2 + A2**2
# Antiparallel spins (symmetric spatial wave -> Peak + Peak constructive interference at 90 deg: (A1 + A2)^2):
I_constructive = (A1 + A2)**2
# Parallel spins (antisymmetric spatial wave / Fermi minus sign = 180 deg Trough shift -> Peak + Trough = (A1 - A2)^2 = 0 at 90 deg!):
I_destructive = (A1 - A2)**2

ax_d.plot(theta_deg, I_constructive, color='#16A34A', lw=2.6,
          label='Constructive Diagram Sum |ℳ₁ + ℳ₂|²: Peak (↑) + Peak (↑) = 4|ℳ₁|² at 90°')
ax_d.plot(theta_deg, I_incoherent, color='#64748B', lw=2.0, ls='--',
          label='No Interference (|ℳ₁|² + |ℳ₂|²): Classical Sum without Phase Arrows')
ax_d.plot(theta_deg, I_destructive, color='#DC2626', lw=2.8,
          label='Destructive Exchange Interference |ℳ₁ − ℳ₂|²: Peak (↑) + Trough (↓) = 0 at 90°!')

ax_d.fill_between(theta_deg, I_destructive, I_constructive, color='#38BDF8', alpha=0.14)

# Highlight 90 deg where Diagram 1 and Diagram 2 have equal magnitude |M1| = |M2|
ax_d.axvline(90.0, color='#475569', lw=1.5, ls=':')
ax_d.scatter([90.0], [1.0], s=75, color='#16A34A', edgecolors='white', lw=1.5, zorder=10)
ax_d.scatter([90.0], [0.0], s=75, color='#DC2626', edgecolors='white', lw=1.5, zorder=10)

ax_d.annotate(
    "Constructive Peak (θ_cm = 90°)\nDiagram 1 (↑) + Diagram 2 (↑) = 2↑\n|ℳ₁ + ℳ₂|² = 1.0 (2× Incoherent!)",
    xy=(90.0, 1.0), xytext=(90.0, 1.72),
    ha='center', va='center', color='#15803D', fontsize=8.4, fontweight='bold',
    bbox=dict(boxstyle='round,pad=0.24', facecolor='#DCFCE7', edgecolor='#16A34A', lw=1.2),
    arrowprops=dict(arrowstyle='->', color='#16A34A', lw=1.6)
)

ax_d.annotate(
    "Complete Destructive Node (θ_cm = 90°)\nDiagram 1 Peak (↑) + Diagram 2 Trough (↓) = 0\nIdentical Parallel-Spin Electrons Never Scatter at 90°!",
    xy=(90.0, 0.0), xytext=(90.0, -0.58),
    ha='center', va='center', color='#B91C1C', fontsize=8.4, fontweight='bold',
    bbox=dict(boxstyle='round,pad=0.24', facecolor='#FEE2E2', edgecolor='#DC2626', lw=1.2),
    arrowprops=dict(arrowstyle='->', color='#DC2626', lw=1.6)
)

ax_d.set_xlim(28, 152)
ax_d.set_ylim(-0.88, 2.35)
ax_d.set_xlabel("Center-of-Mass Scattering Angle θ_cm (degrees)", fontsize=9.8, fontweight='bold')
ax_d.set_ylabel("Normalized Scattering Probability |ℳ_total|²", fontsize=9.8, fontweight='bold')
ax_d.set_title("D. Adding Two Feynman Diagrams (ℳ₁ ± ℳ₂) = Two-Slit Interference of Wave Peaks & Troughs",
               fontsize=10.5, fontweight='bold', color='#0F172A', pad=8)
ax_d.legend(loc='upper center', fontsize=8.0, framealpha=0.95)

out_path = "/Users/john-bramble/.gemini/antigravity/brain/122a2ab8-7d76-4d47-b0c9-7485eef51086/fig_feynman_diagram_to_wave_peaks_troughs.png"
plt.savefig(out_path, dpi=220)
print("Saved Feynman diagram to wave peaks/troughs comparison figure to:", out_path)
