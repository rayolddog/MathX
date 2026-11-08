import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

plt.rcParams['font.sans-serif'] = 'DejaVu Sans'
plt.rcParams['axes.edgecolor'] = '#334155'
plt.rcParams['axes.linewidth'] = 1.0

fig = plt.figure(figsize=(16.5, 12.6), dpi=220, facecolor='#F8FAFC')
gs = fig.add_gridspec(2, 3, height_ratios=[1.42, 1.0], width_ratios=[1.0, 1.0, 1.05],
                      hspace=0.30, wspace=0.22, left=0.05, right=0.97, top=0.89, bottom=0.06)

# Overall Title & Subtitle
fig.suptitle(
    "Feynman's Path Arrows on Wave Peaks & Troughs: Least-Action Straight Paths vs. Curved Non-Least-Action Paths",
    fontsize=14.8, fontweight='bold', color='#0F172A', y=0.968
)
fig.text(
    0.5, 0.935,
    "Straight least-action paths (δS = 0) march coherently along Peaks (↑) or Troughs (↓); curved paths (δS ≠ 0) jump across wavefronts and cancel out",
    ha='center', va='center', fontsize=10.8, color='#334155', style='italic'
)

# ==============================================================================
# TOP PANEL (spanning all 3 columns): 2D Double-Slit Wavefront View
# ==============================================================================
ax_top = fig.add_subplot(gs[0, :])
ax_top.set_facecolor('#0B1120')

lam = 1.0
k = 2 * np.pi / lam
x_screen = 7.0
d_slit = 3.0
y_s1 = +d_slit / 2.0   # +1.5 lambda
y_s2 = -d_slit / 2.0   # -1.5 lambda

ys_grid = np.linspace(-4.2, 4.2, 8001)
r1_grid = np.sqrt(x_screen**2 + (ys_grid - y_s1)**2)
r2_grid = np.sqrt(x_screen**2 + (ys_grid - y_s2)**2)
dr_grid = r2_grid - r1_grid

# Target A: Constructive Bright Fringe (Delta r = +1.0 * lambda, y > 0)
idx_bright = np.argmin(np.abs(dr_grid - 1.0 * lam))
y_bright = ys_grid[idx_bright]

# Target B: Destructive Dark Fringe (Delta r = -0.5 * lambda, y < 0)
idx_dark = np.argmin(np.abs(dr_grid - (-0.5 * lam)))
y_dark = ys_grid[idx_dark]

r1_bright = np.sqrt(x_screen**2 + (y_bright - y_s1)**2)
r2_bright = np.sqrt(x_screen**2 + (y_bright - y_s2)**2)
r1_dark = np.sqrt(x_screen**2 + (y_dark - y_s1)**2)
r2_dark = np.sqrt(x_screen**2 + (y_dark - y_s2)**2)

xg = np.linspace(0.0, x_screen, 500)
yg = np.linspace(-4.2, 4.2, 360)
XX, YY = np.meshgrid(xg, yg)
R1 = np.sqrt(XX**2 + (YY - y_s1)**2)
R2 = np.sqrt(XX**2 + (YY - y_s2)**2)

phi0 = -k * r1_bright
wave_field = 0.5 * (np.cos(k * R1 + phi0) + np.cos(k * R2 + phi0))

ax_top.imshow(
    wave_field, extent=[0.0, x_screen, -4.2, 4.2], origin='lower',
    cmap='RdBu_r', vmin=-1.0, vmax=1.0, alpha=0.25, aspect='auto'
)

# Circular wavefront rings from Slit 1 (cyan) and Slit 2 (amber)
theta_arc = np.linspace(-np.pi/2.3, np.pi/2.3, 300)
for m in range(1, 9):
    r_peak = r1_bright - (8 - m) * lam
    if 0.3 < r_peak < 8.2:
        ax_top.plot(r_peak * np.cos(theta_arc), y_s1 + r_peak * np.sin(theta_arc),
                    color='#38BDF8', lw=0.85, alpha=0.35, ls='-')
        ax_top.plot(r_peak * np.cos(theta_arc), y_s2 + r_peak * np.sin(theta_arc),
                    color='#FBBF24', lw=0.85, alpha=0.35, ls='-')
    r_trough = r_peak - 0.5 * lam
    if 0.3 < r_trough < 8.2:
        ax_top.plot(r_trough * np.cos(theta_arc), y_s1 + r_trough * np.sin(theta_arc),
                    color='#38BDF8', lw=0.65, alpha=0.22, ls='--')
        ax_top.plot(r_trough * np.cos(theta_arc), y_s2 + r_trough * np.sin(theta_arc),
                    color='#FBBF24', lw=0.65, alpha=0.22, ls='--')

# Incident Single-Photon Wave Train on the left (x < 0)
for m_inc in range(3):
    x_crest = -0.35 - m_inc * 0.55
    ax_top.plot([x_crest, x_crest], [-2.5, 2.5], color='#38BDF8', lw=2.0, alpha=0.75 - 0.18*m_inc)
    ax_top.plot([x_crest - 0.275, x_crest - 0.275], [-2.5, 2.5], color='#F472B6', lw=1.3, ls='--', alpha=0.55 - 0.12*m_inc)

ax_top.annotate(
    "Single-Photon\nWave Train\n(E = ħω)", xy=(-0.2, 0.0), xytext=(-1.18, 0.0),
    ha='center', va='center', color='#F8FAFC', fontsize=9.5, fontweight='bold',
    arrowprops=dict(arrowstyle='->', color='#38BDF8', lw=2.0)
)

# Double-Slit Barrier at x = 0
barrier_kw = dict(color='#94A3B8', lw=6.0, solid_capstyle='butt', zorder=5)
ax_top.plot([0, 0], [1.85, 4.2], **barrier_kw)
ax_top.plot([0, 0], [-1.15, 1.15], **barrier_kw)
ax_top.plot([0, 0], [-4.2, -1.85], **barrier_kw)

ax_top.scatter([0], [y_s1], s=95, color='#38BDF8', edgecolors='white', lw=1.5, zorder=8)
ax_top.scatter([0], [y_s2], s=95, color='#FBBF24', edgecolors='white', lw=1.5, zorder=8)
ax_top.text(-0.18, y_s1 + 0.18, "Slit 1", color='#38BDF8', fontsize=10.5, fontweight='bold', ha='right', va='bottom')
ax_top.text(-0.18, y_s2 - 0.18, "Slit 2", color='#FBBF24', fontsize=10.5, fontweight='bold', ha='right', va='top')

def draw_traced_path(ax, x0, y0, x1, y1, r_ref, mode, path_color, arrow_color,
                     label_text, label_pos_frac=0.4, label_offset=(0, 0.35)):
    dx = x1 - x0
    dy = y1 - y0
    L = np.hypot(dx, dy)
    ux, uy = dx / L, dy / L
    nx, ny = -uy, ux

    ls_main = '-' if mode == 'peaks' else '--'
    ax.plot([x0, x1], [y0, y1], color=path_color, lw=2.2, ls=ls_main, alpha=0.92, zorder=6)

    s_vals = np.linspace(0.25, L - 0.12, 400)
    phase_shift = 0.0 if mode == 'peaks' else np.pi
    phase_vals = k * (s_vals - r_ref) + phase_shift
    wave_amp = 0.16
    wx = x0 + s_vals * ux + wave_amp * np.cos(phase_vals) * nx
    wy = y0 + s_vals * uy + wave_amp * np.cos(phase_vals) * ny
    ax.plot(wx, wy, color=path_color, lw=1.3, alpha=0.70, zorder=6)

    arrow_len = 0.36
    for m in range(-12, 12):
        s_pt = r_ref + m * lam
        if 0.65 <= s_pt <= L - 0.15:
            px = x0 + s_pt * ux
            py = y0 + s_pt * uy
            ax.scatter([px], [py], s=36, color=arrow_color, edgecolors='#0F172A', lw=0.8, zorder=9)
            if mode == 'peaks':
                ax.annotate(
                    "", xy=(px, py + arrow_len), xytext=(px, py),
                    arrowprops=dict(arrowstyle='->,head_width=0.26,head_length=0.22', color=arrow_color, lw=2.2),
                    zorder=10
                )
            else:
                ax.annotate(
                    "", xy=(px, py - arrow_len), xytext=(px, py),
                    arrowprops=dict(arrowstyle='->,head_width=0.26,head_length=0.22', color=arrow_color, lw=2.2),
                    zorder=10
                )

    lx = x0 + label_pos_frac * dx + label_offset[0]
    ly = y0 + label_pos_frac * dy + label_offset[1]
    ax.text(
        lx, ly, label_text, color='white', fontsize=8.4, fontweight='bold',
        ha='center', va='center', zorder=12,
        bbox=dict(boxstyle='round,pad=0.23', facecolor='#0F172A', edgecolor=arrow_color, alpha=0.93, lw=1.2)
    )

# Draw the 4 Straight Least-Action (delta S = 0) Paths:
draw_traced_path(
    ax_top, 0, y_s1, x_screen, y_bright, r1_bright,
    mode='peaks', path_color='#38BDF8', arrow_color='#22D3EE',
    label_text="Least-Action Path 1A (δS=0, Slit 1 → y_A): Peaks (↑)",
    label_pos_frac=0.38, label_offset=(-0.35, 0.42)
)

draw_traced_path(
    ax_top, 0, y_s2, x_screen, y_bright, r2_bright,
    mode='peaks', path_color='#FBBF24', arrow_color='#FDE047',
    label_text="Least-Action Path 2A (δS=0, Slit 2 → y_A, Δr=1λ): Peaks (↑)",
    label_pos_frac=0.38, label_offset=(0.15, 0.46)
)

draw_traced_path(
    ax_top, 0, y_s1, x_screen, y_dark, r1_dark,
    mode='peaks', path_color='#38BDF8', arrow_color='#22D3EE',
    label_text="Least-Action Path 1B (δS=0, Slit 1 → y_B): Peaks (↑)",
    label_pos_frac=0.55, label_offset=(0.25, -0.45)
)

draw_traced_path(
    ax_top, 0, y_s2, x_screen, y_dark, r2_dark,
    mode='troughs', path_color='#F472B6', arrow_color='#FB7185',
    label_text="Least-Action Path 2B (δS=0, Slit 2 → y_B, Δr=0.5λ): Troughs (↓)",
    label_pos_frac=0.45, label_offset=(-0.05, -0.56)
)

# ==============================================================================
# NEW: Draw 2 Curved Non-Least-Action Paths (delta S != 0) that jump across
# different wave peaks and troughs with rapidly rotating / scrambling arrows!
# ==============================================================================
def draw_curved_non_stationary_path(ax, x0, y0, x1, y1, arch_amp, freq_wiggle, r_ref,
                                    curve_color, label_text, label_xy):
    t = np.linspace(0.0, 1.0, 500)
    xc = x0 + (x1 - x0) * t
    # Smooth curved trajectory that arches/wiggles across wavefronts
    yc = y0 + (y1 - y0) * t + arch_amp * np.sin(np.pi * t) + 0.35 * arch_amp * np.sin(freq_wiggle * np.pi * t)

    # Plot curved trajectory
    ax.plot(xc, yc, color=curve_color, lw=2.3, ls='-.', alpha=0.95, zorder=11)

    # Compute accumulated arc-length action S(t) = k * int ds AND local wave field phase
    ds = np.hypot(np.diff(xc), np.diff(yc))
    s_arc = np.concatenate([[0.0], np.cumsum(ds)])

    # Sample 8 points along the curve to show arrows jumping between peaks, zero-crossings, and troughs
    sample_indices = np.linspace(55, 455, 8, dtype=int)
    arrow_len = 0.38
    for idx in sample_indices:
        px, py = xc[idx], yc[idx]
        # The phase along the curved path relative to the stationary wavefronts jumps rapidly:
        # Phase angle phi = k * (s_arc[idx] - r_ref) * 1.35
        r_from_slit = np.hypot(px - x0, py - y0)
        # Local phase of the wave field crossed by the curve + excess arc length phase
        phi_local = k * (s_arc[idx] + 0.65 * (s_arc[idx] - r_from_slit) - r_ref)
        # Arrow direction: phi_local = 0 -> UP (+y), phi_local = pi -> DOWN (-y), intermediate -> tilted!
        dx_arr = arrow_len * np.sin(phi_local)
        dy_arr = arrow_len * np.cos(phi_local)

        # Color code dot by whether it is near a peak (cyan), trough (pink), or jumping in between (violet)
        cos_val = np.cos(phi_local)
        if cos_val > 0.55:
            pt_col = '#22D3EE'   # near a peak
        elif cos_val < -0.55:
            pt_col = '#FB7185'   # near a trough
        else:
            pt_col = '#C084FC'   # jumping between peak and trough!

        ax.scatter([px], [py], s=40, color=pt_col, edgecolors='white', lw=0.8, zorder=13)
        ax.annotate(
            "", xy=(px + dx_arr, py + dy_arr), xytext=(px, py),
            arrowprops=dict(arrowstyle='->,head_width=0.26,head_length=0.22', color='#E879F9', lw=2.2),
            zorder=14
        )

    ax.text(
        label_xy[0], label_xy[1], label_text, color='#FAF5FF', fontsize=8.3, fontweight='bold',
        ha='center', va='center', zorder=15,
        bbox=dict(boxstyle='round,pad=0.24', facecolor='#3B0764', edgecolor='#E879F9', alpha=0.94, lw=1.3)
    )

# Curved Path C1 (from Slit 1 -> Bright Fringe y_A, arching high across extra wavefronts)
draw_curved_non_stationary_path(
    ax_top, 0, y_s1, x_screen, y_bright, arch_amp=+1.42, freq_wiggle=2.0, r_ref=r1_bright,
    curve_color='#C084FC',
    label_text="Curved Non-Least-Action Path C1 (δS ≠ 0)\nJumps across Peaks & Troughs (↑ ↗ → ↘ ↓) → Cancels out!",
    label_xy=(3.45, 3.68)
)

# Curved Path C2 (from Slit 2 -> Dark Node y_B, bowing downward across wavefronts)
draw_curved_non_stationary_path(
    ax_top, 0, y_s2, x_screen, y_dark, arch_amp=-1.45, freq_wiggle=3.0, r_ref=r2_dark,
    curve_color='#C084FC',
    label_text="Curved Non-Least-Action Path C2 (δS ≠ 0)\nRapid Phase Changes Across Wavefronts → Net Sum ≈ 0",
    label_xy=(3.75, -3.45)
)

# Detector Screen at x = x_screen and Intensity Envelope |psi(y)|^2
ax_top.plot([x_screen, x_screen], [-4.2, 4.2], color='#E2E8F0', lw=3.0, zorder=7)
I_profile = np.cos(0.5 * k * dr_grid)**2
x_intensity = x_screen + 1.15 * I_profile
ax_top.plot(x_intensity, ys_grid, color='#4ADE80', lw=2.5, zorder=8)
ax_top.fill_betweenx(ys_grid, x_screen, x_intensity, color='#4ADE80', alpha=0.22, zorder=5)
ax_top.text(x_screen + 0.65, 3.82, "Born Weight\n|ψ(y)|² = |ψ₁ + ψ₂|²", color='#4ADE80',
            fontsize=9.3, fontweight='bold', ha='center', va='center')

ax_top.scatter([x_screen], [y_bright], s=135, color='#4ADE80', edgecolors='white', lw=2.0, zorder=12)
ax_top.annotate(
    "BRIGHT FRINGE (y_A)\nLeast-Action Peaks: ↑ + ↑ = 2↑\nCurved Paths Cancel (∑ e^{iS/ħ} ≈ 0)\n|ψ|² = 4a² (Constructive)",
    xy=(x_screen + 1.15, y_bright), xytext=(x_screen + 1.40, y_bright + 0.15),
    color='#4ADE80', fontsize=8.7, fontweight='bold', va='center',
    bbox=dict(boxstyle='round,pad=0.3', facecolor='#064E3B', edgecolor='#4ADE80', lw=1.3),
    arrowprops=dict(arrowstyle='->', color='#4ADE80', lw=1.8)
)

ax_top.scatter([x_screen], [y_dark], s=135, color='#F87171', edgecolors='white', lw=2.0, zorder=12)
ax_top.annotate(
    "DARK FRINGE NODE (y_B)\nLeast-Action: Peak ↑ + Trough ↓ = 0\nCurved Paths Also Cancel ≈ 0\n|ψ|² = 0 (Destructive)",
    xy=(x_screen, y_dark), xytext=(x_screen + 1.40, y_dark - 0.15),
    color='#FCA5A5', fontsize=8.7, fontweight='bold', va='center',
    bbox=dict(boxstyle='round,pad=0.3', facecolor='#450A0A', edgecolor='#F87171', lw=1.3),
    arrowprops=dict(arrowstyle='->', color='#F87171', lw=1.8)
)

legend_elements = [
    Line2D([0], [0], color='#22D3EE', marker='^', markersize=8, lw=2.0, label='Slit 1 Least-Action Peak Paths (δS = 0, Φ = 0°, ↑)'),
    Line2D([0], [0], color='#FDE047', marker='^', markersize=8, lw=2.0, label='Slit 2 Least-Action Peak Path (δS = 0, Φ = 0°, ↑)'),
    Line2D([0], [0], color='#FB7185', marker='v', markersize=8, lw=2.0, ls='--', label='Slit 2 Least-Action Trough Path (δS = 0, Φ = 180°, ↓)'),
    Line2D([0], [0], color='#E879F9', marker='o', markersize=7, lw=2.2, ls='-.', label='Curved Non-Least-Action Paths (δS ≠ 0: Arrows Jump & Cancel)'),
]
ax_top.legend(handles=legend_elements, loc='lower left', facecolor='#0F172A', edgecolor='#475569',
              labelcolor='white', fontsize=8.3, framealpha=0.94)

ax_top.set_xlim(-1.85, 10.55)
ax_top.set_ylim(-4.2, 4.2)
ax_top.set_xlabel("Propagation Distance x (in wavelengths λ)", fontsize=10.5, fontweight='bold', color='#1E293B')
ax_top.set_ylabel("Transverse Screen Position y (in λ)", fontsize=10.5, fontweight='bold', color='#1E293B')
ax_top.set_title(
    "A. 2D Wavefront View: Straight Least-Action Paths (δS = 0) vs. Curved Non-Least-Action Paths (δS ≠ 0) Jumping Across Peaks & Troughs",
    fontsize=11.2, fontweight='bold', color='#0F172A', pad=9
)

# ==============================================================================
# BOTTOM-LEFT PANEL (B): Constructive Least-Action Paths (1A + 2A)
# ==============================================================================
ax_bl = fig.add_subplot(gs[1, 0])
ax_bl.set_facecolor('#FFFFFF')

s_plot = np.linspace(0, 3.0, 500)
w1_bright = np.cos(2 * np.pi * (s_plot - 3.0))
w2_bright = np.cos(2 * np.pi * (s_plot - 3.0))

y_off1 = 1.75
y_off2 = 0.0
y_off_sum = -1.95

ax_bl.axhline(y_off1, color='#CBD5E1', lw=1.0, ls=':')
ax_bl.axhline(y_off2, color='#CBD5E1', lw=1.0, ls=':')
ax_bl.axhline(y_off_sum, color='#CBD5E1', lw=1.0, ls=':')

ax_bl.plot(s_plot, y_off1 + 0.55 * w1_bright, color='#0284C7', lw=2.2)
ax_bl.plot(s_plot, y_off2 + 0.55 * w2_bright, color='#D97706', lw=2.2)
ax_bl.plot(s_plot, y_off_sum + 0.55 * (w1_bright + w2_bright), color='#16A34A', lw=2.6)

for s_pk in [0.5, 1.0, 1.5, 2.0, 2.5, 3.0]:
    is_peak = np.isclose((3.0 - s_pk) % 1.0, 0.0) or np.isclose((3.0 - s_pk) % 1.0, 1.0)
    if is_peak:
        ax_bl.annotate("", xy=(s_pk, y_off1 + 0.55), xytext=(s_pk, y_off1),
                       arrowprops=dict(arrowstyle='->,head_width=0.2,head_length=0.18', color='#0284C7', lw=2.0))
        ax_bl.annotate("", xy=(s_pk, y_off2 + 0.55), xytext=(s_pk, y_off2),
                       arrowprops=dict(arrowstyle='->,head_width=0.2,head_length=0.18', color='#D97706', lw=2.0))
        ax_bl.annotate("", xy=(s_pk, y_off_sum + 1.10), xytext=(s_pk, y_off_sum),
                       arrowprops=dict(arrowstyle='->,head_width=0.24,head_length=0.22', color='#16A34A', lw=2.4))
    else:
        ax_bl.annotate("", xy=(s_pk, y_off1 - 0.55), xytext=(s_pk, y_off1),
                       arrowprops=dict(arrowstyle='->,head_width=0.2,head_length=0.18', color='#DB2777', lw=1.8))
        ax_bl.annotate("", xy=(s_pk, y_off2 - 0.55), xytext=(s_pk, y_off2),
                       arrowprops=dict(arrowstyle='->,head_width=0.2,head_length=0.18', color='#DB2777', lw=1.8))
        ax_bl.annotate("", xy=(s_pk, y_off_sum - 1.10), xytext=(s_pk, y_off_sum),
                       arrowprops=dict(arrowstyle='->,head_width=0.24,head_length=0.22', color='#15803D', lw=2.4))

ax_bl.axvline(3.0, color='#16A34A', lw=2.0, ls='--')
ax_bl.text(3.06, 2.35, "Screen y_A\n↑ + ↑ = 2↑\n|ψ|² = 4a²",
           color='#15803D', fontsize=8.4, fontweight='bold', va='top')

ax_bl.text(0.05, y_off1 + 0.66, "Path 1A (δS=0): Peaks (↑) & Troughs (↓)", color='#0284C7', fontsize=8.4, fontweight='bold')
ax_bl.text(0.05, y_off2 + 0.66, "Path 2A (δS=0, Δr=1λ): Peaks Align!", color='#D97706', fontsize=8.4, fontweight='bold')
ax_bl.text(0.05, y_off_sum + 1.15, "Constructive Sum: 2× Amplitude", color='#15803D', fontsize=8.4, fontweight='bold')

ax_bl.set_xlim(0.0, 3.95)
ax_bl.set_ylim(-3.25, 2.65)
ax_bl.set_yticks([y_off_sum, y_off2, y_off1])
ax_bl.set_yticklabels(["Sum", "Slit 2", "Slit 1"], fontsize=8.8, fontweight='bold')
ax_bl.set_xlabel("Distance Along Path (in λ)", fontsize=9.3, fontweight='bold')
ax_bl.set_title("B. Least-Action Constructive (1A + 2A)\nPeak Arrows (↑) Reinforce Peak Arrows (↑)",
                fontsize=9.8, fontweight='bold', color='#0F172A', pad=6)

# ==============================================================================
# BOTTOM-MIDDLE PANEL (C): Destructive Least-Action Paths (1B + 2B)
# ==============================================================================
ax_bm = fig.add_subplot(gs[1, 1])
ax_bm.set_facecolor('#FFFFFF')

w1_dark = np.cos(2 * np.pi * (s_plot - 3.0))
w2_dark = np.cos(2 * np.pi * (s_plot - 3.0) + np.pi)

ax_bm.axhline(y_off1, color='#CBD5E1', lw=1.0, ls=':')
ax_bm.axhline(y_off2, color='#CBD5E1', lw=1.0, ls=':')
ax_bm.axhline(y_off_sum, color='#CBD5E1', lw=1.0, ls=':')

ax_bm.plot(s_plot, y_off1 + 0.55 * w1_dark, color='#0284C7', lw=2.2)
ax_bm.plot(s_plot, y_off2 + 0.55 * w2_dark, color='#DB2777', lw=2.2, ls='--')
ax_bm.plot(s_plot, y_off_sum + 0.0 * s_plot, color='#DC2626', lw=2.6)

for s_pk in [0.5, 1.0, 1.5, 2.0, 2.5, 3.0]:
    is_peak_1 = np.isclose((3.0 - s_pk) % 1.0, 0.0) or np.isclose((3.0 - s_pk) % 1.0, 1.0)
    if is_peak_1:
        ax_bm.annotate("", xy=(s_pk, y_off1 + 0.55), xytext=(s_pk, y_off1),
                       arrowprops=dict(arrowstyle='->,head_width=0.2,head_length=0.18', color='#0284C7', lw=2.0))
        ax_bm.annotate("", xy=(s_pk, y_off2 - 0.55), xytext=(s_pk, y_off2),
                       arrowprops=dict(arrowstyle='->,head_width=0.2,head_length=0.18', color='#DB2777', lw=2.0))
        ax_bm.annotate("", xy=(s_pk - 0.04, y_off_sum + 0.42), xytext=(s_pk - 0.04, y_off_sum),
                       arrowprops=dict(arrowstyle='->,head_width=0.16,head_length=0.14', color='#0284C7', lw=1.5, alpha=0.75))
        ax_bm.annotate("", xy=(s_pk + 0.04, y_off_sum - 0.42), xytext=(s_pk + 0.04, y_off_sum),
                       arrowprops=dict(arrowstyle='->,head_width=0.16,head_length=0.14', color='#DB2777', lw=1.5, alpha=0.75))
        ax_bm.scatter([s_pk], [y_off_sum], s=24, color='#DC2626', zorder=10)
    else:
        ax_bm.annotate("", xy=(s_pk, y_off1 - 0.55), xytext=(s_pk, y_off1),
                       arrowprops=dict(arrowstyle='->,head_width=0.2,head_length=0.18', color='#DB2777', lw=2.0))
        ax_bm.annotate("", xy=(s_pk, y_off2 + 0.55), xytext=(s_pk, y_off2),
                       arrowprops=dict(arrowstyle='->,head_width=0.2,head_length=0.18', color='#D97706', lw=2.0))
        ax_bm.annotate("", xy=(s_pk - 0.04, y_off_sum - 0.42), xytext=(s_pk - 0.04, y_off_sum),
                       arrowprops=dict(arrowstyle='->,head_width=0.16,head_length=0.14', color='#DB2777', lw=1.5, alpha=0.75))
        ax_bm.annotate("", xy=(s_pk + 0.04, y_off_sum + 0.42), xytext=(s_pk + 0.04, y_off_sum),
                       arrowprops=dict(arrowstyle='->,head_width=0.16,head_length=0.14', color='#D97706', lw=1.5, alpha=0.75))
        ax_bm.scatter([s_pk], [y_off_sum], s=24, color='#DC2626', zorder=10)

ax_bm.axvline(3.0, color='#DC2626', lw=2.0, ls='--')
ax_bm.text(3.06, 2.35, "Screen y_B\n↑ + ↓ = 0\n|ψ|² = 0",
           color='#DC2626', fontsize=8.4, fontweight='bold', va='top')

ax_bm.text(0.05, y_off1 + 0.66, "Path 1B (δS=0): Arrives on PEAK (↑)", color='#0284C7', fontsize=8.4, fontweight='bold')
ax_bm.text(0.05, y_off2 + 0.66, "Path 2B (δS=0, Δr=0.5λ): TROUGH (↓)", color='#DB2777', fontsize=8.4, fontweight='bold')
ax_bm.text(0.05, y_off_sum + 0.62, "Destructive Sum: Peak ↑ Cancels Trough ↓", color='#DC2626', fontsize=8.4, fontweight='bold')

ax_bm.set_xlim(0.0, 3.95)
ax_bm.set_ylim(-3.25, 2.65)
ax_bm.set_yticks([y_off_sum, y_off2, y_off1])
ax_bm.set_yticklabels(["Sum", "Slit 2", "Slit 1"], fontsize=8.8, fontweight='bold')
ax_bm.set_xlabel("Distance Along Path (in λ)", fontsize=9.3, fontweight='bold')
ax_bm.set_title("C. Least-Action Destructive (1B + 2B)\nPeak Arrows (↑) Cancel Trough Arrows (↓)",
                fontsize=9.8, fontweight='bold', color='#0F172A', pad=6)

# ==============================================================================
# BOTTOM-RIGHT PANEL (D): Curved Non-Least-Action Paths (C1 & C2)
# Showing how jumping across peaks & troughs scrambles arrow angles & cancels!
# ==============================================================================
ax_br = fig.add_subplot(gs[1, 2])
ax_br.set_facecolor('#FAF5FF')

# Top half of Panel D: Wave profile experienced along a curved path (chirped/irregular phase!)
s_curve = np.linspace(0, 3.0, 600)
# Because the curved path cuts obliquely & variably across wavefronts, its phase jumps non-linearly:
phi_curved = 2 * np.pi * (s_curve + 0.42 * np.sin(1.7 * np.pi * s_curve) + 0.18 * s_curve**2)
w_curved = np.cos(phi_curved)

y_top_d = 1.35
ax_br.axhline(y_top_d, color='#D8B4FE', lw=1.0, ls=':')
ax_br.plot(s_curve, y_top_d + 0.55 * w_curved, color='#9333EA', lw=2.2, ls='-.')

# Sample 7 points along the curved path showing arrows pointing at scrambled angles!
sample_s = np.linspace(0.3, 2.85, 7)
angles_list = []
for sp in sample_s:
    ph = 2 * np.pi * (sp + 0.42 * np.sin(1.7 * np.pi * sp) + 0.18 * sp**2)
    angles_list.append(ph)
    dx_a = 0.28 * np.sin(ph)
    dy_a = 0.48 * np.cos(ph)
    ax_br.scatter([sp], [y_top_d + 0.55 * np.cos(ph)], s=30, color='#9333EA', zorder=9)
    ax_br.annotate("", xy=(sp + dx_a, y_top_d + dy_a), xytext=(sp, y_top_d),
                   arrowprops=dict(arrowstyle='->,head_width=0.2,head_length=0.18', color='#7E22CE', lw=2.0))

ax_br.text(0.08, y_top_d + 0.72, "Curved Path (δS ≠ 0): Cuts Across Wavefronts → Arrows Jump!",
           color='#6B21A8', fontsize=8.4, fontweight='bold')

# Bottom half of Panel D: Head-to-tail vector addition of neighboring curved paths vs Least-Action paths
y_bot_d = -1.45
ax_br.axhline(0.0, color='#CBD5E1', lw=1.0)

# Left sub-diagram: Least-action bundle (arrows aligned -> large resultant)
x_la_start = 0.45
y_la_start = -2.45
for k_step in range(4):
    ax_br.annotate("", xy=(x_la_start, y_la_start + (k_step + 1) * 0.42),
                   xytext=(x_la_start, y_la_start + k_step * 0.42),
                   arrowprops=dict(arrowstyle='->,head_width=0.18,head_length=0.15', color='#0284C7', lw=2.2))
ax_br.text(x_la_start, -2.85, "Least-Action\nBundle (δS = 0)\nAligned Arrows!", ha='center', va='top',
           color='#0284C7', fontsize=8.2, fontweight='bold')

# Right sub-diagram: Curved non-least-action bundle (arrows spin in a closed circle -> net = 0!)
cx_poly, cy_poly = 2.35, -1.55
r_poly = 0.62
N_poly = 6
for m_p in range(N_poly):
    ang1 = 2 * np.pi * m_p / N_poly
    ang2 = 2 * np.pi * (m_p + 1) / N_poly
    p1x, p1y = cx_poly + r_poly * np.cos(ang1), cy_poly + r_poly * np.sin(ang1) * 1.25
    p2x, p2y = cx_poly + r_poly * np.cos(ang2), cy_poly + r_poly * np.sin(ang2) * 1.25
    ax_br.annotate("", xy=(p2x, p2y), xytext=(p1x, p1y),
                   arrowprops=dict(arrowstyle='->,head_width=0.2,head_length=0.18', color='#9333EA', lw=2.1))

ax_br.scatter([cx_poly + r_poly], [cy_poly], s=45, color='#DC2626', zorder=12)
ax_br.text(cx_poly, cy_poly, "Net Sum\n∑ e^{iS/ħ} ≈ 0\n(Cancels!)", ha='center', va='center',
           color='#6B21A8', fontsize=8.3, fontweight='bold')
ax_br.text(cx_poly, -2.85, "Curved Paths (δS ≠ 0)\nScrambled Arrows Wind\nin a Circle & Cancel Out", ha='center', va='top',
           color='#7E22CE', fontsize=8.2, fontweight='bold')

ax_br.set_xlim(0.0, 3.55)
ax_br.set_ylim(-3.45, 2.65)
ax_br.set_yticks([y_bot_d, y_top_d])
ax_br.set_yticklabels(["Vector Sum\nof Bundle", "Curved Path\nWave & Arrows"], fontsize=8.4, fontweight='bold')
ax_br.set_xlabel("Distance Along Curved Path / Phasor Addition", fontsize=9.0, fontweight='bold')
ax_br.set_title("D. Curved Non-Least-Action Paths (C1, C2)\nJumping Between Peaks & Troughs → Self-Cancellation",
                fontsize=9.6, fontweight='bold', color='#581C87', pad=6)

out_path = "/Users/john-bramble/.gemini/antigravity/brain/122a2ab8-7d76-4d47-b0c9-7485eef51086/fig_feynman_arrows_peaks_troughs_two_slit.png"
plt.savefig(out_path, dpi=220)
print("Saved updated diagram with curved non-least-action paths to:", out_path)
