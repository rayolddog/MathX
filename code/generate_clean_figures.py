import os
os.environ["MPLCONFIGDIR"] = "/tmp/mplconfig"
import numpy as np
import matplotlib.pyplot as plt
import json

rng = np.random.default_rng(20261014)

with open("/Users/john-bramble/.gemini/antigravity/scratch/selector_verification/detailed_results.json") as f:
    det = json.load(f)

# -------------------------------------------------------------------------
# FIGURE 1: Sample Event-by-Event Microsimulation Trajectories & Born Convergence
# -------------------------------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.6), dpi=200)

# Panel A: Sample trajectories of s1(t) starting at s1(0) = 0.80
ax = axes[0]
s0 = 0.80
sigma = 0.22
for tr in range(14):
    s = s0
    traj = [s]
    for step in range(50):
        if s <= 1e-6 or s >= 1.0 - 1e-6:
            traj.append(s)
            continue
        std = sigma * np.sqrt(s * (1.0 - s))
        if std <= min(s, 1.0 - s):
            s += rng.choice([-1.0, 1.0]) * std
        elif s > 0.5:
            v = 1.0 - s
            a_high = (sigma**2) * s
            p1 = a_high / (v + a_high)
            s = 1.0 if rng.random() < p1 else s - a_high
        else:
            a_low = (sigma**2) * (1.0 - s)
            p0 = a_low / (s + a_low)
            s = 0.0 if rng.random() < p0 else s + a_low
        s = np.clip(s, 0.0, 1.0)
        traj.append(s)
    col = "#1b9e77" if traj[-1] > 0.5 else "#d95f02"
    alp = 0.75 if tr < 7 else 0.45
    ax.plot(range(len(traj)), traj, color=col, lw=1.6, alpha=alp)

ax.axhline(0.80, color="black", ls=":", lw=1.5, label="Initial Share $s_1(0) = 0.80$")
ax.plot([], [], color="#1b9e77", lw=2, label="Site 1 Registers ($s_1 \\to 1$, 80% of events)")
ax.plot([], [], color="#d95f02", lw=2, label="Site 2 Registers ($s_1 \\to 0$, 20% of events)")
ax.set_xlim(0, 42)
ax.set_ylim(-0.03, 1.03)
ax.set_xlabel("Microsimulation Exchange Step $t$", fontsize=11)
ax.set_ylabel("Site 1 Energy Share $s_1(t) = e_1(t)/\\hbar\\omega$", fontsize=11)
ax.set_title("A. Event-by-Event Trajectories ($\\Delta s \\propto \\sqrt{s_1 s_2}\\cos\\Delta\\phi$)", fontsize=11.5, fontweight="bold")
ax.grid(True, alpha=0.3)
ax.legend(loc="center right", fontsize=8.8)

# Panel B: Microsimulation Registration Frequency vs. Born Rule Across S in [0.1, 0.9]
ax2 = axes[1]
s_grid = np.array(det["s_sweep_theta_072_lam_004"]["s_grid"])
p_vac = np.array(det["s_sweep_theta_072_lam_004"]["sqrt_plus_vacuum"])
# Also run pure ungated conserving walk across s_grid
p_ungated = []
for s_init in s_grid:
    # 20,000 trials of pure conserving walk
    s = np.full(20000, s_init)
    active = np.ones(20000, dtype=bool)
    for _ in range(1000):
        idx = np.where(active)[0]
        if len(idx) == 0:
            break
        u = s[idx]
        std = 0.25 * np.sqrt(u * (1.0 - u))
        min_s = np.minimum(u, 1.0 - u)
        signs = rng.choice([-1.0, 1.0], size=len(idx))
        du_int = signs * std
        a_low = (0.25**2) * (1.0 - u)
        p0_low = a_low / np.maximum(u + a_low, 1e-15)
        du_low = np.where(rng.random(len(idx)) < p0_low, -u, a_low)
        v = 1.0 - u
        a_high = (0.25**2) * u
        p0_high = a_high / np.maximum(v + a_high, 1e-15)
        du_high = np.where(rng.random(len(idx)) < p0_high, v, -a_high)
        near_low = (std > min_s) & (u <= 0.5) & (u > 0.0)
        near_high = (std > min_s) & (u > 0.5) & (u < 1.0)
        du = np.where(near_low, du_low, np.where(near_high, du_high, du_int))
        u_new = np.clip(u + du, 0.0, 1.0)
        s[idx] = u_new
        active[idx] = (u_new > 0.005) & (u_new < 0.995)
    p_ungated.append(float(np.mean(s >= 0.5)))
p_ungated = np.array(p_ungated)

s_line = np.linspace(0, 1, 100)
ax2.plot(s_line, s_line, 'k--', lw=1.8, label="Exact Born Rule ($P_1 = S$)")
ax2.plot(s_grid, p_ungated, 'o', color="#1b9e77", markersize=7.5, label="Microsimulation: Discrete Absorber ($\\theta=1$, $p=2.001$)")
ax2.plot(s_grid, p_vac, 's', color="#2b5c8f", markersize=6.5, alpha=0.85, label="Microsimulation: Gated + Vacuum Term ($\\theta=0.72$, $p=2.031$)")
ax2.set_xlim(0, 1)
ax2.set_ylim(0, 1)
ax2.set_xlabel("Initial Wave-Energy Share $S = |A_1|^2 / (|A_1|^2 + |A_2|^2)$", fontsize=11)
ax2.set_ylabel("Observed Registration Probability $P_1$", fontsize=11)
ax2.set_title("B. Emergent Born Statistics from Microsimulation", fontsize=11.5, fontweight="bold")
ax2.grid(True, alpha=0.3)
ax2.legend(loc="upper left", fontsize=8.5)

plt.tight_layout()
plt.savefig("/Users/john-bramble/.gemini/antigravity/scratch/selector_verification/fig_clean_1_trajectories_and_born.png", dpi=200)
plt.savefig("/Users/john-bramble/.gemini/antigravity/brain/f91ca8ab-6741-4e72-b1f0-af94b132a121/fig_clean_1_trajectories_and_born.png", dpi=200)
plt.close()

# -------------------------------------------------------------------------
# FIGURE 2: Vacuum Cross-Term Compensation & 10-Site Aperture Test
# -------------------------------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.6), dpi=200)

# Panel A: Gated Deviation With vs. Without Non-Random Vacuum Cross-Term
ax = axes[0]
p_sqrt_only = np.array(det["s_sweep_theta_072_lam_004"]["sqrt_step"])
ax.axhline(0.0, color="black", ls=":", lw=1.5, label="Exact Born Comparator ($P_1 - S = 0$)")
ax.plot(s_grid, p_sqrt_only - s_grid, 's--', color="#d95f02", lw=2, markersize=6.5,
        label="Gated Walk Without Vacuum Drift ($\\epsilon_{\\mathrm{vac}} = 0$, $p=2.21$)")
ax.plot(s_grid, p_vac - s_grid, 'D-', color="#2b5c8f", lw=2.5, markersize=7,
        label="Gated Walk + Vacuum Cross-Term ($\\epsilon_{\\mathrm{vac}} = 0.018$, $p=2.03$)")
ax.set_ylim(-0.04, 0.04)
ax.set_xlabel("Initial Wave-Energy Share $S = |A_1|^2 / (|A_1|^2 + |A_2|^2)$", fontsize=11)
ax.set_ylabel("Deviation from Born Rule ($P_1 - S$)", fontsize=11)
ax.set_title("A. Vacuum Cross-Term Compensation of Gate Bias", fontsize=11.5, fontweight="bold")
ax.grid(True, alpha=0.3)
ax.legend(loc="upper left", fontsize=8.5)

# Panel B: 10-Site Detector Aperture Test (1 Bright Site A0=3 => Born 0.500; 9 Dim Sites Ai=1)
ax2 = axes[1]
categories = [
    "Discrete Absorber\n($\\theta = 1.0$,\n$\\epsilon_{\\mathrm{vac}} = 0$)",
    "Gated Continuum\n($\\theta = 0.72$,\n$\\epsilon_{\\mathrm{vac}} = 0$)",
    "Gated + Vacuum\n($\\theta = 0.72$,\n$\\epsilon_{\\mathrm{vac}} = 0.012$)",
    "Snapshot Threshold\n(Unconserved\nExtreme-Value)"
]
values = [0.5001, 0.5389, 0.5023, 0.7910]
bar_colors = ["#1b9e77", "#d95f02", "#2b5c8f", "#e7298a"]
bars = ax2.bar(categories, values, color=bar_colors, width=0.52, edgecolor="black", lw=0.8)
ax2.axhline(0.5000, color="black", ls="--", lw=1.5, label="Exact Born Weight ($P_{\\mathrm{bright}} = 0.500$)")
for b, v in zip(bars, values):
    ax2.text(b.get_x() + b.get_width()/2, v + 0.02, f"{v:.4f}", ha="center", va="bottom", fontsize=9.5, fontweight="bold")
ax2.set_ylim(0, 0.92)
ax2.set_ylabel("Bright-Site Registration Frequency $P_{\\mathrm{bright}}$", fontsize=11)
ax2.set_title("B. 10-Site Aperture Test (1 Bright + 9 Dim Sites)", fontsize=11.5, fontweight="bold")
ax2.grid(axis="y", alpha=0.3)
ax2.legend(loc="upper left", fontsize=9)

plt.tight_layout()
plt.savefig("/Users/john-bramble/.gemini/antigravity/scratch/selector_verification/fig_clean_2_vacuum_and_10site.png", dpi=200)
plt.savefig("/Users/john-bramble/.gemini/antigravity/brain/f91ca8ab-6741-4e72-b1f0-af94b132a121/fig_clean_2_vacuum_and_10site.png", dpi=200)
plt.close()
print("Clean figures generated successfully.")
