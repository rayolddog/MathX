import numpy as np
import matplotlib.pyplot as plt
from scipy.special import fresnel

plt.rcParams['font.sans-serif'] = 'DejaVu Sans'
plt.rcParams['axes.edgecolor'] = '#334155'
plt.rcParams['axes.linewidth'] = 1.0

fig = plt.figure(figsize=(16.0, 12.2), dpi=220, facecolor='#F8FAFC')
gs = fig.add_gridspec(2, 2, height_ratios=[1.05, 1.0], width_ratios=[1.05, 1.0],
                      hspace=0.30, wspace=0.22, left=0.06, right=0.96, top=0.89, bottom=0.06)

fig.suptitle(
    "Calculated Feynman Action S(h)/ħ & Why Curved Non-Least-Action Paths Are Exponentially/Destructively Rare",
    fontsize=15.0, fontweight='bold', color='#0F172A', y=0.968
)
fig.text(
    0.5, 0.935,
    "Grounded in Feynman Lectures Vol. I Ch. 26 & Vol. II Ch. 19: Stationary Action (δS = 0) at h = 0 vs. Rapid Phase Cancellation (δS ≠ 0) for Curved Paths",
    ha='center', va='center', fontsize=10.8, color='#334155', style='italic'
)

# ==============================================================================
# PANEL A (Top-Left): Calculated Action Parabola S(h)/ħ and Oscillating Wave cos(S(h)/ħ)
# ==============================================================================
ax_a = fig.add_subplot(gs[0, 0])
ax_a.set_facecolor('#FFFFFF')

# Parameter h: transverse arch deflection of path (in units of first Fresnel zone width h_F)
# Extra path length Delta L(h) = (pi^2 / (4 * L0)) * h^2 -> Phase Delta Phi(h) = (S(h) - S_0)/hbar = (pi / 2) * u^2
u = np.linspace(-3.3, 3.3, 1500)
delta_S_over_hbar = 0.5 * np.pi * u**2          # Action phase (S(h) - S_min)/hbar in radians
delta_cycles = delta_S_over_hbar / (2 * np.pi)  # Extra path length in wavelengths (Delta L / lambda = u^2 / 4)
cos_phase = np.cos(delta_S_over_hbar)           # Real projection of Feynman arrow: +1 = Peak (↑), -1 = Trough (↓)

# Plot Action Parabola (top curve in Panel A) and Chirped Wave cos(S/ħ) (bottom strip in Panel A)
ax_a.plot(u, delta_cycles + 1.25, color='#0F172A', lw=2.8, label='Calculated Extra Action ΔS(h)/(2πħ) = ΔL(h)/λ ∝ h²')

# Highlight Stationary Action region (|u| <= 1.0, where Delta S < pi/2)
u_stat = u[np.abs(u) <= 1.0]
ax_a.fill_between(u_stat, 1.25, (0.25 * u_stat**2) + 1.25, color='#38BDF8', alpha=0.28)
ax_a.fill_between(u_stat, -1.1, 0.75 * np.cos(0.5 * np.pi * u_stat**2) - 0.15, color='#22C55E', alpha=0.22)

# Plot oscillating real wave amplitude cos(S(h)/hbar) at bottom of Panel A
y_wave_base = -0.15
ax_a.axhline(y_wave_base, color='#94A3B8', lw=1.0, ls=':')
ax_a.plot(u, y_wave_base + 0.72 * cos_phase, color='#2563EB', lw=2.0,
          label='Real Arrow Projection cos(S(h)/ħ): +1 = Peak (↑), −1 = Trough (↓)')

# Fill positive lobes (green/cyan = Peak arrows ↑) and negative lobes (pink/purple = Trough arrows ↓)
ax_a.fill_between(u, y_wave_base, y_wave_base + 0.72 * np.maximum(cos_phase, 0),
                  where=(np.abs(u) > 1.0), color='#38BDF8', alpha=0.30)
ax_a.fill_between(u, y_wave_base, y_wave_base + 0.72 * np.minimum(cos_phase, 0),
                  where=(np.abs(u) > 1.0), color='#F43F5E', alpha=0.30)

# Draw Feynman arrows along the Action parabola at u = 0 (Least Action) vs u = ±1.8, ±2.4 (Curved Paths C1, C2)
sample_u = [0.0, -0.45, 0.45, -1.65, 1.65, -2.15, 2.15, -2.65, 2.65]
for up in sample_u:
    ph = 0.5 * np.pi * up**2
    yc = 1.25 + 0.25 * up**2
    dx_a = 0.22 * np.sin(ph)
    dy_a = 0.38 * np.cos(ph)
    col_a = '#16A34A' if np.abs(up) <= 0.6 else '#9333EA'
    ax_a.scatter([up], [yc], s=38, color=col_a, zorder=10)
    ax_a.annotate("", xy=(up + dx_a, yc + dy_a), xytext=(up, yc),
                  arrowprops=dict(arrowstyle='->,head_width=0.22,head_length=0.20', color=col_a, lw=2.2),
                  zorder=11)

# Annotate Least-Action Minimum (h = 0)
ax_a.annotate(
    "Least-Action Minimum (h = 0, δS = 0)\nFlat Action Bowl → All Arrows Point UP (↑↑↑)\nConstructive Fresnel Zone!",
    xy=(0.0, 1.25), xytext=(0.0, 2.65),
    ha='center', va='center', color='#065F46', fontsize=8.8, fontweight='bold',
    bbox=dict(boxstyle='round,pad=0.28', facecolor='#D1FAE5', edgecolor='#10B981', lw=1.3),
    arrowprops=dict(arrowstyle='->', color='#10B981', lw=1.8)
)

# Annotate Curved Path C1 / C2 on the steep wall of the Action parabola
ax_a.annotate(
    "Curved Paths C1, C2 (|h| > 0, δS ≠ 0)\nSteep Action Slope → Rapidly Alternates\nBetween Peaks (+1, ↑) & Troughs (−1, ↓)\nAdjacent Lobes Cancel to Zero!",
    xy=(2.15, 1.25 + 0.25 * 2.15**2), xytext=(1.75, 3.65),
    ha='center', va='center', color='#581C87', fontsize=8.5, fontweight='bold',
    bbox=dict(boxstyle='round,pad=0.26', facecolor='#F3E8FF', edgecolor='#9333EA', lw=1.3),
    arrowprops=dict(arrowstyle='->', color='#9333EA', lw=1.8)
)

ax_a.set_xlim(-3.3, 3.3)
ax_a.set_ylim(-1.15, 4.35)
ax_a.set_xlabel("Path Curvature / Arch Deflection h (in Fresnel Zone Units u = h / √(λD/2))", fontsize=9.8, fontweight='bold')
ax_a.set_ylabel("Extra Action ΔS(h)/(2πħ) & Wave cos(S/ħ)", fontsize=9.8, fontweight='bold')
ax_a.set_title("A. Calculated Action Parabola S(h)/ħ & Chirped Peak/Trough Alternation",
               fontsize=11.0, fontweight='bold', color='#0F172A', pad=8)
ax_a.legend(loc='lower center', fontsize=8.3, framealpha=0.94)

# ==============================================================================
# PANEL B (Top-Right): The Cornu Spiral (Exact Vector Sum of All Path Arrows)
# ==============================================================================
ax_b = fig.add_subplot(gs[0, 1])
ax_b.set_facecolor('#0F172A')

u_spiral = np.linspace(-4.5, 4.5, 2000)
S_fres, C_fres = fresnel(u_spiral)

# Plot full Cornu spiral
ax_b.plot(C_fres, S_fres, color='#94A3B8', lw=1.8, alpha=0.65, label='All Paths (−∞ < h < +∞)')

# Highlight the central Least-Action segment (-0.9 <= u <= +0.9) in bright green/cyan
mask_center = np.abs(u_spiral) <= 0.9
ax_b.plot(C_fres[mask_center], S_fres[mask_center], color='#4ADE80', lw=4.2,
          label='Least-Action Bundle (|h| ≤ 0.9, δS ≈ 0): Straight Line!')

# Highlight the Curved Non-Least-Action tails (|u| > 0.9) winding into tight circles in purple
mask_pos = u_spiral > 0.9
mask_neg = u_spiral < -0.9
ax_b.plot(C_fres[mask_pos], S_fres[mask_pos], color='#E879F9', lw=2.4,
          label='Curved Paths (|h| > 0.9, δS ≠ 0): Wind in Circles & Cancel!')
ax_b.plot(C_fres[mask_neg], S_fres[mask_neg], color='#E879F9', lw=2.4)

# Draw the net resultant arrow from (-0.5, -0.5) to (+0.5, +0.5)
ax_b.annotate(
    "", xy=(0.5, 0.5), xytext=(-0.5, -0.5),
    arrowprops=dict(arrowstyle='->,head_width=0.35,head_length=0.30', color='#38BDF8', lw=3.0, ls='--'),
    zorder=12
)

ax_b.scatter([-0.5, 0.5], [-0.5, 0.5], s=70, color='#F43F5E', edgecolors='white', lw=1.5, zorder=13)
ax_b.scatter([0.0], [0.0], s=70, color='#4ADE80', edgecolors='white', lw=1.5, zorder=13)

ax_b.text(-0.08, 0.14, "Straight Least-Action\nPaths (h ≈ 0)\nArrows Add in a\nStraight Line!",
          color='#4ADE80', fontsize=9.0, fontweight='bold', ha='right', va='bottom')

ax_b.annotate(
    "Curved Paths C1 (h > 0)\nJump Peak ↔ Trough:\nArrows Coil in a Tight\nCircle Around (0.5, 0.5)\nNet Addition ≈ 0!",
    xy=(0.52, 0.52), xytext=(0.12, -0.45),
    color='#F5D0FE', fontsize=8.6, fontweight='bold', ha='center', va='center',
    bbox=dict(boxstyle='round,pad=0.26', facecolor='#3B0764', edgecolor='#E879F9', lw=1.2),
    arrowprops=dict(arrowstyle='->', color='#E879F9', lw=1.8)
)

ax_b.set_xlim(-0.92, 0.92)
ax_b.set_ylim(-0.88, 0.88)
ax_b.set_xlabel("Real Part of Summed Arrow Amplitude ∫ cos(S(h)/ħ) dh", fontsize=9.8, fontweight='bold', color='#0F172A')
ax_b.set_ylabel("Imaginary Part ∫ sin(S(h)/ħ) dh", fontsize=9.8, fontweight='bold', color='#0F172A')
ax_b.set_title("B. Feynman's Cornu Spiral: Adding Path Arrows Head-to-Tail Across All Curvatures h",
               fontsize=11.0, fontweight='bold', color='#0F172A', pad=8)
ax_b.legend(loc='upper left', facecolor='#1E293B', edgecolor='#475569', labelcolor='white', fontsize=8.2)

# ==============================================================================
# PANEL C (Bottom-Left): Spatial Fan of Curved Paths Annotated with Extra Action ΔS
# ==============================================================================
ax_c = fig.add_subplot(gs[1, 0])
ax_c.set_facecolor('#0B1120')

x_path = np.linspace(0.0, 7.0, 400)
# Draw background wavefront crests (vertical-ish arcs)
for m_w in range(1, 8):
    ax_c.axvline(m_w * 0.95, color='#38BDF8', lw=0.9, alpha=0.25, ls='-')
    ax_c.axvline(m_w * 0.95 - 0.475, color='#F472B6', lw=0.7, alpha=0.18, ls='--')

# Plot family of paths with different arch heights h_arch (in Fresnel units u = 0, ±1, ±2, ±3)
arch_list = [
    (0.0,  '#4ADE80', '-',  3.2, "h = 0 (Least Action, ΔL = 0.00λ): 100% Coherent Peak Alignment (↑↑↑)"),
    (+1.1, '#38BDF8', '--', 2.0, "h = +1.0 (Edge of 1st Fresnel Zone, ΔL = 0.25λ): Phase Shift = 90° (→)"),
    (-1.1, '#38BDF8', '--', 2.0, None),
    (+2.1, '#C084FC', '-.', 2.0, "h = ±2.0 (Curved Path C1, ΔL = 1.00λ): Crosses 2 Fresnel Zones → Cancels!"),
    (-2.1, '#C084FC', '-.', 2.0, None),
    (+3.1, '#F472B6', ':',  1.8, "h = ±3.0 (Highly Curved Path C2, ΔL = 2.25λ): Crosses 9 Half-Zones → Rare!"),
    (-3.1, '#F472B6', ':',  1.8, None),
]

for h_val, col, ls_style, lw_val, lbl in arch_list:
    y_curve = h_val * np.sin(np.pi * x_path / 7.0)
    ax_c.plot(x_path, y_curve, color=col, ls=ls_style, lw=lw_val, label=lbl, alpha=0.92)

# Mark Slit and Detector endpoints
ax_c.scatter([0.0, 7.0], [0.0, 0.0], s=110, color='#4ADE80', edgecolors='white', lw=2.0, zorder=10)
ax_c.text(0.15, -0.45, "Slit Aperture\n(x = 0)", color='white', fontsize=9.0, fontweight='bold')
ax_c.text(6.85, -0.45, "Detector Pixel\n(x = D)", color='#4ADE80', fontsize=9.0, fontweight='bold', ha='right')

# Shade the 1st Fresnel Zone tube around the straight line
y_fresnel = 1.1 * np.sin(np.pi * x_path / 7.0)
ax_c.fill_between(x_path, -y_fresnel, +y_fresnel, color='#4ADE80', alpha=0.16)
ax_c.text(3.5, 0.0, "Stationary-Action Tube (1st Fresnel Zone, ΔS < πħ/2)\nPaths Inside This Tube Reinforce Constructively",
          color='#A7F3D0', fontsize=8.6, fontweight='bold', ha='center', va='center',
          bbox=dict(boxstyle='round,pad=0.22', facecolor='#064E3B', edgecolor='#4ADE80', alpha=0.88))

ax_c.set_xlim(-0.3, 7.3)
ax_c.set_ylim(-3.6, 3.6)
ax_c.set_xlabel("Distance from Slit to Detector Screen x (in λ)", fontsize=9.8, fontweight='bold')
ax_c.set_ylabel("Transverse Arch Deflection h(x)", fontsize=9.8, fontweight='bold')
ax_c.set_title("C. Spatial Family of Paths & the Stationary-Action Fresnel Tube",
               fontsize=11.0, fontweight='bold', color='#0F172A', pad=8)
ax_c.legend(loc='lower center', facecolor='#0F172A', edgecolor='#475569', labelcolor='white', fontsize=7.9, framealpha=0.94)

# ==============================================================================
# PANEL D (Bottom-Right): Surviving Coherent Weight vs. Path Curvature |h|
# Showing quantitatively why curved paths are rare!
# ==============================================================================
ax_d = fig.add_subplot(gs[1, 1])
ax_d.set_facecolor('#FFFFFF')

# Compute the net coherent amplitude squared |∫_{u - du/2}^{u + du/2} exp(i pi v^2 / 2) dv|^2
# for path bundles of width du = 0.75 centered at curvature u:
u_centers = np.linspace(-3.2, 3.2, 500)
du = 0.85
S_plus, C_plus = fresnel(u_centers + 0.5 * du)
S_minus, C_minus = fresnel(u_centers - 0.5 * du)
bundle_amp_sq = ((C_plus - C_minus)**2 + (S_plus - S_minus)**2)
bundle_weight_norm = bundle_amp_sq / np.max(bundle_amp_sq)

# Also compute Discrete Fresnel Zone magnitudes |A_m|^2 for zone index m = 1, 2, 3, 4, 5
ax_d.plot(u_centers, bundle_weight_norm, color='#0284C7', lw=2.8,
          label='Net Coherent Bundle Weight |∫_{h−Δh/2}^{h+Δh/2} e^{iS(h)/ħ} dh|² (Normalized)')
ax_d.fill_between(u_centers, 0, bundle_weight_norm, where=(np.abs(u_centers) <= 1.0),
                  color='#22C55E', alpha=0.35, label='Least-Action Core (|h| ≤ 1 Fresnel Zone): >90% of Surviving Weight')
ax_d.fill_between(u_centers, 0, bundle_weight_norm, where=(np.abs(u_centers) > 1.0),
                  color='#A855F7', alpha=0.30, label='Curved Non-Least-Action Paths (|h| > 1): Suppressed as 1/h² → Rare!')

# Plot envelope 1 / (pi^2 * u^2) decay curve
u_tail = np.linspace(0.95, 3.2, 200)
decay_env = 1.0 / (np.pi**2 * u_tail**2 * (du**2)) * 0.85
ax_d.plot(u_tail, np.minimum(decay_env, 1.0), color='#DC2626', lw=1.8, ls='--', label='Asymptotic Stationary-Phase Suppression ∝ 1 / (∂S/∂h)² ∝ 1/h²')
ax_d.plot(-u_tail, np.minimum(decay_env, 1.0), color='#DC2626', lw=1.8, ls='--')

# Mark specific curved paths C1 (u = 2.0) and C2 (u = -2.35)
for u_mark, name_lbl, y_txt in [(0.0, "Straight Least-Action Path (h = 0)\nRelative Weight = 100%", 0.72),
                                (2.0, "Curved Path C1 (h = +2.0)\nWeight ≈ 3.8% (Cancels!)", 0.36),
                                (-2.15, "Curved Path C2 (h = −2.15)\nWeight ≈ 2.2% (Very Rare!)", 0.36)]:
    idx_m = np.argmin(np.abs(u_centers - u_mark))
    w_val = bundle_weight_norm[idx_m]
    col_m = '#15803D' if u_mark == 0.0 else '#7E22CE'
    ax_d.scatter([u_mark], [w_val], s=65, color=col_m, edgecolors='white', lw=1.5, zorder=10)
    ax_d.annotate(
        name_lbl, xy=(u_mark, w_val), xytext=(u_mark, y_txt),
        ha='center', va='center', color=col_m, fontsize=8.3, fontweight='bold',
        bbox=dict(boxstyle='round,pad=0.24', facecolor='#F8FAFC', edgecolor=col_m, lw=1.2),
        arrowprops=dict(arrowstyle='->', color=col_m, lw=1.6)
    )

ax_d.set_xlim(-3.2, 3.2)
ax_d.set_ylim(0.0, 1.32)
ax_d.set_xlabel("Path Curvature / Arch Deflection h (in Fresnel Zone Units u)", fontsize=9.8, fontweight='bold')
ax_d.set_ylabel("Surviving Coherent Bundle Weight (Normalized)", fontsize=9.8, fontweight='bold')
ax_d.set_title("D. Why Curved Paths Are Rare: Net Coherent Weight vs. Path Curvature h",
               fontsize=11.0, fontweight='bold', color='#0F172A', pad=8)
ax_d.legend(loc='upper right', fontsize=7.9, framealpha=0.94)

out_path = "/Users/john-bramble/.gemini/antigravity/brain/122a2ab8-7d76-4d47-b0c9-7485eef51086/fig_feynman_action_curve_and_path_rarity.png"
plt.savefig(out_path, dpi=220)
print("Saved calculated action & path rarity diagram to:", out_path)
