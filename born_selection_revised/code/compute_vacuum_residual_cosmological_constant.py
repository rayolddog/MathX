#!/usr/bin/env python3
"""
Quantitative Verification & Figure Generation:
Quantum Vacuum Time-Orientation Cancellation and the Baryon-Asymmetry Residual
for the Cosmological Constant (Article 0004 — Stage 02 Revision).
"""

import math
import os
import numpy as np
import matplotlib.pyplot as plt

# Fundamental physical constants (CODATA 2022 / PDG 2024 / Planck 2018)
c = 2.99792458e8            # Speed of light (m/s)
hbar = 1.054571817e-34      # Reduced Planck constant (J s)
G = 6.67430e-11             # Newton gravitational constant (m^3 kg^-1 s^-2)
eV = 1.602176634e-19        # Electron-volt in Joules

# Cosmological parameters (Planck 2018: TT,TE,EE+lowE+lensing)
H0_kmsMpc = 67.4            # Hubble constant (km/s/Mpc)
Mpc_m = 3.085677581491367e22
H0 = (H0_kmsMpc * 1e3) / Mpc_m
Omega_Lambda = 0.6847
Omega_b = 0.0493
n_gamma_m3 = 4.107e8        # CMB photon number density (m^-3) at T_CMB = 2.7255 K
eta_B = 6.12e-10            # Baryon-to-photon ratio (n_B - n_Bbar)/n_gamma

# Observed energy densities
rho_crit = (3.0 * H0**2 / (8.0 * math.pi * G)) * c**2  # J/m^3
rho_Lambda_obs = Omega_Lambda * rho_crit               # J/m^3
rho_baryon_obs = Omega_b * rho_crit                    # J/m^3
E_Lambda_obs_eV = (rho_Lambda_obs * (hbar * c)**3)**0.25 / eV

# Particle physics mass scales (PDG 2024 MSbar at 2 GeV for light quarks)
E_Planck_eV = math.sqrt(hbar * c**5 / G) / eV          # 1.2209e28 eV
v_EW_eV = 246.22e9                                     # Electroweak VEV (eV)
Lambda_QCD_eV = 217e6                                  # QCD scale (eV)
m_p_eV = 938.272088e6                                  # Proton mass (eV)
m_u_eV = 2.16e6                                        # Up quark mass (eV)
m_d_eV = 4.67e6                                        # Down quark mass (eV)

# 1. Standard QFT (+t forward light cone only, no +-t cancellation)
rho_QFT_Planck = (E_Planck_eV * eV)**4 / (hbar * c)**3
rho_QFT_EW = (v_EW_eV * eV)**4 / (hbar * c)**3
rho_QFT_QCD = (Lambda_QCD_eV * eV)**4 / (hbar * c)**3

# 2. Route A: Microscopic Quark-Vacuum Residual under +-t cancellation weighted by eta_B
m_q_iso_eV = 0.5 * (m_u_eV + m_d_eV)                   # Isospin average (3.415 MeV)
m_q_proton_eV = (2.0 * m_u_eV + m_d_eV) / 3.0          # Proton valence average (2.997 MeV)
m_q_neutron_eV = (m_u_eV + 2.0 * m_d_eV) / 3.0         # Neutron valence average (3.833 MeV)
# BBN primordial nucleon mix (~87.5% p, 12.5% n by number -> valence quark mix)
m_q_eff_eV = 0.5 * (m_u_eV + m_d_eV)

E_res_iso_eV = eta_B * m_q_iso_eV
E_res_n_eV = eta_B * m_q_neutron_eV
E_res_p_eV = eta_B * m_q_proton_eV

rho_RouteA_iso = (E_res_iso_eV * eV)**4 / (hbar * c)**3
rho_RouteA_n = (E_res_n_eV * eV)**4 / (hbar * c)**3
rho_RouteA_p = (E_res_p_eV * eV)**4 / (hbar * c)**3

# 3. Route B: Macroscopic Relic Baryon Polarization of the Isotropic 4*pi Vacuum
rho_B_relic = eta_B * n_gamma_m3 * (m_p_eV * eV)
rho_RouteB_4pi = 4.0 * math.pi * rho_B_relic
rho_RouteB_C4 = 0.5 * (16.0 / math.pi)**2 * rho_B_relic

# Monte Carlo simulation of N vacuum oscillator pairs with random time orientation s_i in {+1, -1}
# biased by CP/T baryon asymmetry eta_B
rng = np.random.default_rng(20261109)
N_modes = 200000
# Unbiased random time orientations (+t vs -t) vs biased by eta_B
signs_unbiased = rng.choice([-1.0, 1.0], size=N_modes)
cum_avg_unbiased = np.abs(np.cumsum(signs_unbiased)) / np.arange(1, N_modes + 1)

out_dir = os.path.dirname(os.path.abspath(__file__))
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13.5, 5.6), dpi=200)

# Panel A: Vacuum Energy Density Comparison Across 130 Orders of Magnitude
models = [
    "Standard QFT\n(Planck $+t$ only)",
    "Standard QFT\n(EW $+t$ only)",
    "Standard QFT\n(QCD $+t$ only)",
    "Route A: $\\pm t$ Cancel\n$(\\eta_B m_q c^2)^4/(\\hbar c)^3$",
    "Route B: $\\pm t$ Cancel\n$4\\pi \\eta_B n_\\gamma m_p c^2$",
    "Observed $\\rho_\\Lambda$\n(Planck 2018)"
]
log_rhos = [
    math.log10(rho_QFT_Planck),
    math.log10(rho_QFT_EW),
    math.log10(rho_QFT_QCD),
    math.log10(rho_RouteA_iso),
    math.log10(rho_RouteB_4pi),
    math.log10(rho_Lambda_obs)
]
colors = ["#c0392b", "#d35400", "#e67e22", "#1f77b4", "#2980b9", "#27ae60"]

bars = ax1.barh(models, log_rhos, color=colors, edgecolor="black", height=0.58, alpha=0.88)
ax1.axvline(math.log10(rho_Lambda_obs), color="#27ae60", linestyle="--", linewidth=1.8,
            label=f"Observed $\\log_{{10}}\\rho_\\Lambda = {math.log10(rho_Lambda_obs):.2f}$")
for bar, val in zip(bars, log_rhos):
    offset = 3.0 if val > 0 else -18.0
    ax1.text(val + offset, bar.get_y() + bar.get_height()/2,
             f"$10^{{{val:.1f}}}$", va="center", ha="left", fontsize=9.5, fontweight="bold")

ax1.set_xlim(-28, 130)
ax1.set_xlabel("$\\log_{10}\\!\\left(\\rho_{\\mathrm{vac}} \\; [\\mathrm{J/m^3}]\\right)$", fontsize=11)
ax1.set_title("A. Standard $+t$ Light-Cone Sum vs. $\\pm t$ Cancellation + $\\eta_B$ Residual", fontsize=11.5, fontweight="bold")
ax1.grid(axis="x", linestyle=":", alpha=0.5)
ax1.legend(loc="lower right", fontsize=9.5)

# Panel B: Zoom-In Linear Comparison Around the Observed Dark Energy Density (in units of 10^-10 J/m^3)
zoom_labels = [
    "Route A: Quark\n($\\bar{m}_q = 3.42$ MeV)",
    "Route A: Neutron\n($\\bar{m}_{q,n} = 3.83$ MeV)",
    "Route B: Isotropic\n($4\\pi\\,\\rho_B$)",
    "Route B: $\\mathbb{C}^4$ Proj.\n($\\frac{1}{2}(16/\\pi)^2\\rho_B$)",
    "Observed $\\rho_\\Lambda$\n(Planck 2018)"
]
zoom_vals = np.array([
    rho_RouteA_iso,
    rho_RouteA_n,
    rho_RouteB_4pi,
    rho_RouteB_C4,
    rho_Lambda_obs
]) * 1e10  # in units of 10^-10 J/m^3
zoom_colors = ["#1f77b4", "#3498db", "#2980b9", "#5dade2", "#27ae60"]

x_pos = np.arange(len(zoom_labels))
rects = ax2.bar(x_pos, zoom_vals, color=zoom_colors, edgecolor="black", width=0.55, alpha=0.9)
ax2.axhline(rho_Lambda_obs * 1e10, color="#27ae60", linestyle="--", linewidth=1.8,
            label=f"Planck 2018: $5.25\\times 10^{{-10}}\\;\\mathrm{{J/m^3}}$")

for rect, val in zip(rects, zoom_vals):
    ax2.text(rect.get_x() + rect.get_width()/2, val + 0.15,
             f"{val:.2f}", ha="center", va="bottom", fontsize=9.5, fontweight="bold")

ax2.set_xticks(x_pos)
ax2.set_xticklabels(zoom_labels, fontsize=9)
ax2.set_ylim(0, 7.5)
ax2.set_ylabel("Energy Density $\\rho_\\Lambda$ ($10^{-10}\\;\\mathrm{J/m^3}$)", fontsize=11)
ax2.set_title("B. Linear Comparison with Measured Expansion ($10^{-10}\\;\\mathrm{J/m^3}$)", fontsize=11.5, fontweight="bold")
ax2.grid(axis="y", linestyle=":", alpha=0.5)
ax2.legend(loc="upper right", fontsize=9.5)

plt.tight_layout()
fig_path = os.path.join(out_dir, "fig_vacuum_time_cancellation_cosmological_constant.png")
plt.savefig(fig_path, dpi=200)
print(f"Saved figure to {fig_path}")
print(f"Observed E_Lambda = {E_Lambda_obs_eV*1e3:.3f} meV, rho_Lambda = {rho_Lambda_obs:.3e} J/m^3")
print(f"Route A (isospin m_q = 3.415 MeV): E_res = {E_res_iso_eV*1e3:.3f} meV, rho = {rho_RouteA_iso:.3e} J/m^3 (ratio = {rho_RouteA_iso/rho_Lambda_obs:.3f})")
print(f"Route A (neutron m_q = 3.833 MeV): E_res = {E_res_n_eV*1e3:.3f} meV, rho = {rho_RouteA_n:.3e} J/m^3 (ratio = {rho_RouteA_n/rho_Lambda_obs:.3f})")
print(f"Route B (4*pi*rho_B): rho = {rho_RouteB_4pi:.3e} J/m^3 (ratio = {rho_RouteB_4pi/rho_Lambda_obs:.3f})")
print(f"Route B (0.5*(16/pi)^2*rho_B): rho = {rho_RouteB_C4:.3e} J/m^3 (ratio = {rho_RouteB_C4/rho_Lambda_obs:.3f})")
