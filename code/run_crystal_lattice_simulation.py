#!/usr/bin/env python3
"""
2D Condensed-Matter Crystal Lattice, Off-Center Wave-Packet Alignment (Tie-Breaking),
and Surface-Conditioned Vacuum Phase Simulation for Aperture Article 0002 (Revision 2.3 — 2026-10-29)

Tests whether the atom-by-atom pre-modulus phase mechanism holds up when the absorber
sites are placed on a physical 2D condensed-matter crystal lattice with realistic
off-center wave-packet alignment:
1. Spatial 2D Lattice Geometry:
   - Each detector region/fringe is a 2D square crystal lattice of size L_x x L_y
     (4 x 4 = 16 lattice atoms in Patch 1 and 4 x 4 = 16 lattice atoms in Patch 2,
     total M_total = 32 lattice atoms).
2. Breaking Multi-Atom Ties via Realistic Off-Center Wave-Packet Alignment:
   - When an incident particle wave packet is modeled as a uniform flat-top across the grid
     (16-fold symmetry) or perfectly aligned with the geometric center r_grid = (1.5, 1.5) of
     an even 4x4 detector lattice (4-fold symmetry of the four inner sites (1,1), (1,2), (2,1), (2,2)
     at identical distance r = 1/sqrt(2) a_lat), simulations exhibit an uncharacteristically
     high frequency of ties among competing atoms in the bright spot.
   - In a real experiment, the incident wave-packet center r_inc = r_grid + delta_r0 has a
     random sub-lattice spatial offset delta_r0 = (delta_x0, delta_y0) ~ Uniform(-0.5 a_lat, +0.5 a_lat).
     This generic alignment variance breaks the 4-fold/16-fold lattice symmetry at t = 0 so a
     single atom lies closest to the wave-packet center, sharply reducing tie frequency and
     accelerating consolidation onto the dominant atom.
3. Coherent Incident Bloch/Wavevector Phase Across Lattice Sites:
   - The incident wave has a physical spatial wavevector phase across the lattice:
     phi_inc(x, y) = phi_0 + k_x * x + k_y * y.
4. Surface-Conditioned Quantum Vacuum Phase Distribution (Casimir-Polder + Bias Field E_bias):
   - Near the condensed-matter surface under bias field E_bias, the virtual vacuum phase
     phi_vac(x, y, t) follows a surface-conditioned von Mises distribution (kappa_surf > 0)
     and is dynamically pulled toward the local patch's dominant absorbing atom.
"""

from pathlib import Path
import json
import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import chi2

CODE_DIR = Path(__file__).resolve().parent


def build_2d_lattice_neighbors(Lx, Ly):
    """
    Returns array of lattice coordinates (M, 2) and array of nearest-neighbor pairs
    for a 2D Lx x Ly crystal lattice patch.
    """
    coords = np.array([(ix, iy) for ix in range(Lx) for iy in range(Ly)], dtype=float)
    M = len(coords)
    nn_pairs = []
    for a in range(M):
        for b in range(a + 1, M):
            dist2 = np.sum((coords[a] - coords[b]) ** 2)
            if np.isclose(dist2, 1.0) or np.isclose(dist2, 2.0):
                nn_pairs.append((a, b))
    return coords, np.array(nn_pairs, dtype=int)


def fit_born_exponent(s_grid, p_obs, n_trials):
    """
    Fits P_1(S) = S^(p/2) / (S^(p/2) + (1 - S)^(p/2)) by weighted nonlinear least squares.
    """
    s_arr = np.asarray(s_grid, dtype=float)
    p_arr = np.asarray(p_obs, dtype=float)

    def model(s, p):
        a = np.power(s, 0.5 * p)
        b = np.power(1.0 - s, 0.5 * p)
        return a / (a + b)

    sigma_binom = np.sqrt(np.maximum(p_arr * (1.0 - p_arr), 1e-6) / n_trials)
    popt, pcov = curve_fit(model, s_arr, p_arr, p0=[2.0], sigma=sigma_binom, absolute_sigma=True)
    p_hat = float(popt[0])
    se_p = float(np.sqrt(pcov[0, 0]))
    ci95_half = 1.95996 * se_p

    rms_born = float(np.sqrt(np.mean((p_arr - s_arr) ** 2)))
    chi2_born = float(np.sum(n_trials * ((p_arr - s_arr) ** 2) / (s_arr * (1.0 - s_arr))))
    pval_born = float(1.0 - chi2.cdf(chi2_born, df=len(s_arr)))

    return {
        "p_hat": round(p_hat, 4),
        "se_p": round(se_p, 4),
        "ci95": [round(p_hat - ci95_half, 4), round(p_hat + ci95_half, 4)],
        "rms_vs_born": round(rms_born, 5),
        "chi2_vs_born": round(chi2_born, 2),
        "dof": int(len(s_arr)),
        "p_value_vs_born": round(pval_born, 4),
    }


def simulate_2d_crystal_lattice_aperture(
    s_patch1=0.5,
    Lx=4,
    Ly=4,
    theta=0.45,
    sigma=0.25,
    commit_rate=0.02,
    kappa_surf=2.5,
    kappa_pull=0.35,
    alpha_seed=None,
    surface_vac_eps=0.006,
    align_mode="off_center_random",
    beam_waist=1.40,
    n_trials=12000,
    max_steps=8000,
    seed=20261029,
):
    """
    Simulates two 2D crystal lattice patches (each Lx x Ly = M_patch individual lattice atoms,
    total M_tot = 2 * M_patch atoms).
    align_mode:
      - "centered_flat": uniform flat-top envelope across all 16 lattice sites (16-fold symmetry)
      - "centered_gaussian": Gaussian wave packet perfectly aligned with grid center (1.5, 1.5)
        (4-fold tie symmetry of the 4 central atoms (1,1), (1,2), (2,1), (2,2))
      - "off_center_random": realistic sub-lattice random spatial offset delta_r0 ~ U(-0.5, 0.5) a_lat
        breaking exact grid-centering symmetry.
    """
    rng = np.random.default_rng(seed)
    coords, nn_arr = build_2d_lattice_neighbors(Lx, Ly)
    M_patch = Lx * Ly
    M_tot = 2 * M_patch
    grid_center = np.array([(Lx - 1) / 2.0, (Ly - 1) / 2.0])

    # 1. Spatial wave-packet envelope across the 2D detector lattice:
    if align_mode == "centered_flat":
        env1 = np.ones((n_trials, M_patch))
        env2 = np.ones((n_trials, M_patch))
        phase_mod = 0.25
    elif align_mode == "centered_gaussian":
        d2 = np.sum((coords - grid_center[None, :]) ** 2, axis=1)
        env1 = np.tile(np.exp(-d2 / (2.0 * (beam_waist ** 2))), (n_trials, 1))
        env2 = np.tile(np.exp(-d2 / (2.0 * (beam_waist ** 2))), (n_trials, 1))
        phase_mod = 0.15
    else:
        # Realistic random sub-lattice offset delta_r0 in [-0.5, +0.5]^2 a_lat:
        off1 = rng.uniform(-0.5, 0.5, size=(n_trials, 1, 2))
        off2 = rng.uniform(-0.5, 0.5, size=(n_trials, 1, 2))
        c1 = grid_center[None, None, :] + off1
        c2 = grid_center[None, None, :] + off2
        d2_1 = np.sum((coords[None, :, :] - c1) ** 2, axis=2)
        d2_2 = np.sum((coords[None, :, :] - c2) ** 2, axis=2)
        env1 = np.exp(-d2_1 / (2.0 * (beam_waist ** 2)))
        env2 = np.exp(-d2_2 / (2.0 * (beam_waist ** 2)))
        phase_mod = 0.15
    if alpha_seed is not None:
        phase_mod = float(alpha_seed)

    # 2. Coherent incident wavevector phase across the 2D lattice: phi_inc(r) = phi_0 + k . r
    phi_0 = rng.uniform(0.0, 2.0 * np.pi, size=(n_trials, 1))
    kx, ky = 0.6, 0.4
    k_dot_r = kx * coords[:, 0] + ky * coords[:, 1]
    phi_inc = np.tile(k_dot_r, (n_trials, 2)) + phi_0

    # 3. Surface-conditioned quantum vacuum phase distribution (von Mises with concentration kappa_surf):
    phi_surf = rng.uniform(0.0, 2.0 * np.pi, size=(n_trials, 1))
    if kappa_surf > 0.0:
        phi_vac = phi_surf + rng.vonmises(0.0, kappa_surf, size=(n_trials, M_tot))
    else:
        phi_vac = rng.uniform(0.0, 2.0 * np.pi, size=(n_trials, M_tot))

    # 4. Initial atomic/electronic phases on the crystal lattice:
    phi_det = rng.uniform(0.0, 2.0 * np.pi, size=(n_trials, M_tot))

    # 5. Pre-modulus 3-way phase alignment combined with spatial wave-packet envelope:
    align_1 = env1 * (
        1.0 + phase_mod * np.cos(phi_det[:, :M_patch] - 0.5 * (phi_inc[:, :M_patch] + phi_vac[:, :M_patch]))
    )
    align_2 = env2 * (
        1.0 + phase_mod * np.cos(phi_det[:, M_patch:] - 0.5 * (phi_inc[:, M_patch:] + phi_vac[:, M_patch:]))
    )
    w1 = align_1 / np.sum(align_1, axis=1, keepdims=True)
    w2 = align_2 / np.sum(align_2, axis=1, keepdims=True)

    # Quantify initial bright-spot tie frequency (top 2 atoms within 5% or 2% relative share):
    w1_sorted = np.sort(w1, axis=1)[:, ::-1]
    rel_gap = (w1_sorted[:, 0] - w1_sorted[:, 1]) / np.maximum(w1_sorted[:, 0], 1e-30)
    tie_rate_5pct = round(float(np.mean(rel_gap < 0.05)), 4)
    tie_rate_2pct = round(float(np.mean(rel_gap < 0.02)), 4)
    mean_max_initial_atom_share = round(float(np.mean(s_patch1 * w1_sorted[:, 0])), 4)

    S = np.zeros((n_trials, M_tot), dtype=float)
    S[:, :M_patch] = s_patch1 * w1
    S[:, M_patch:] = (1.0 - s_patch1) * w2

    active = np.ones(n_trials, dtype=bool)
    winner_atom = np.full(n_trials, -1, dtype=int)
    steps = np.zeros(n_trials, dtype=int)
    n_nn = len(nn_arr)

    for step_idx in range(1, max_steps + 1):
        idx = np.where(active)[0]
        if len(idx) == 0:
            break
        cur_S = S[idx]
        cur_phi = phi_det[idx]
        cur_vac = phi_vac[idx]
        K = len(idx)
        steps[idx] += 1

        # Check individual lattice-atom bandgap threshold gate: s_m > theta
        if commit_rate > 0.0 and theta < 1.0:
            exp_mask = cur_S > theta
            if np.any(exp_mask):
                rates = np.where(exp_mask, commit_rate * cur_S, 0.0)
                r_tot = np.sum(rates, axis=1)
                committed = rng.random(K) < r_tot
                if np.any(committed):
                    c_idx = idx[committed]
                    c_rates = rates[committed]
                    c_rtot = r_tot[committed]
                    cum_probs = np.cumsum(c_rates / c_rtot[:, None], axis=1)
                    u_pick = rng.random(len(c_idx))[:, None]
                    winner_atom[c_idx] = np.argmax(u_pick < cum_probs, axis=1)
                    active[c_idx] = False
                    keep = ~committed
                    idx = idx[keep]
                    cur_S = cur_S[keep]
                    cur_phi = cur_phi[keep]
                    cur_vac = cur_vac[keep]
                    K = len(idx)
                    if K == 0:
                        break

        # Patch-local dominant lattice atom pulls the local surface-conditioned vacuum phase:
        dom1 = np.argmax(cur_S[:, :M_patch], axis=1)
        dom2 = M_patch + np.argmax(cur_S[:, M_patch:], axis=1)
        dom_phi1 = cur_phi[np.arange(K), dom1][:, None]
        dom_phi2 = cur_phi[np.arange(K), dom2][:, None]
        if kappa_surf > 0.0:
            cur_vac[:, :M_patch] = (
                cur_vac[:, :M_patch]
                + kappa_pull * np.sin(dom_phi1 - cur_vac[:, :M_patch])
                + rng.vonmises(0.0, kappa_surf, size=(K, M_patch)) * 0.25
            )
            cur_vac[:, M_patch:] = (
                cur_vac[:, M_patch:]
                + kappa_pull * np.sin(dom_phi2 - cur_vac[:, M_patch:])
                + rng.vonmises(0.0, kappa_surf, size=(K, M_patch)) * 0.25
            )
        else:
            cur_vac[:, :M_patch] = (
                cur_vac[:, :M_patch]
                + kappa_pull * np.sin(dom_phi1 - cur_vac[:, :M_patch])
                + rng.uniform(-np.pi, np.pi, size=(K, M_patch))
            )
            cur_vac[:, M_patch:] = (
                cur_vac[:, M_patch:]
                + kappa_pull * np.sin(dom_phi2 - cur_vac[:, M_patch:])
                + rng.uniform(-np.pi, np.pi, size=(K, M_patch))
            )

        # Lattice acoustic/optical phonon thermal phase modulation on the crystal sites:
        cur_phi = cur_phi + 0.20 * np.sin(cur_vac - cur_phi) + rng.uniform(-np.pi, np.pi, size=(K, M_tot))

        # Tier A (every 3rd step): Intra-lattice nearest-neighbor exciton consolidation
        if step_idx % 3 == 0:
            bond_idx = rng.integers(0, n_nn, size=K)
            a1 = nn_arr[bond_idx, 0]
            b1 = nn_arr[bond_idx, 1]
            a2 = a1 + M_patch
            b2 = b1 + M_patch
            for col_i, col_j in [(a1, b1), (a2, b2)]:
                si = cur_S[np.arange(K), col_i]
                sj = cur_S[np.arange(K), col_j]
                tot = si + sj
                u = np.where(tot > 0.0, si / np.maximum(tot, 1e-30), 0.0)
                dphi = (cur_phi[np.arange(K), col_i] - cur_vac[np.arange(K), col_i]) - (
                    cur_phi[np.arange(K), col_j] - cur_vac[np.arange(K), col_j]
                )
                signs = np.where(np.cos(dphi) >= 0.0, 1.0, -1.0)
                std = sigma * np.sqrt(u * (1.0 - u))
                min_s = np.minimum(u, 1.0 - u)
                du_int = signs * std
                rand_u = rng.random(K)
                a_low = (sigma ** 2) * (1.0 - u)
                du_low = np.where(rand_u < a_low / np.maximum(u + a_low, 1e-30), -u, a_low)
                v = 1.0 - u
                a_high = (sigma ** 2) * u
                du_high = np.where(rand_u < a_high / np.maximum(v + a_high, 1e-30), v, -a_high)
                near_low = (std > min_s) & (u <= 0.5) & (u > 0.0)
                near_high = (std > min_s) & (u > 0.5) & (u < 1.0)
                in_open = (u > 0.0) & (u < 1.0)
                du = np.where(near_low, du_low, np.where(near_high, du_high, du_int))
                u_new = np.where(in_open, u + du, u)
                u_new = np.where(u_new <= 1e-15, 0.0, np.where(u_new >= 1.0 - 1e-15, 1.0, u_new))
                cur_S[np.arange(K), col_i] = u_new * tot
                cur_S[np.arange(K), col_j] = (1.0 - u_new) * tot

        # Tier B: Shared wave / virtual vacuum cross-term exchange across the aperture
        perm = rng.permuted(np.tile(np.arange(M_tot), (K, 1)), axis=1)
        for p in range(0, M_tot - 1, 2):
            i_col = perm[:, p]
            j_col = perm[:, p + 1]
            si = cur_S[np.arange(K), i_col]
            sj = cur_S[np.arange(K), j_col]
            tot = si + sj
            u = np.where(tot > 0.0, si / np.maximum(tot, 1e-30), 0.0)
            if surface_vac_eps != 0.0:
                du_vac = surface_vac_eps * np.sqrt(u * (1.0 - u)) * (np.sqrt(1.0 - u) - np.sqrt(u))
                u = np.where((u > 0.0) & (u < 1.0), u + du_vac, u)
            dphi = (cur_phi[np.arange(K), i_col] - cur_vac[np.arange(K), i_col]) - (
                cur_phi[np.arange(K), j_col] - cur_vac[np.arange(K), j_col]
            )
            signs = np.where(np.cos(dphi) >= 0.0, 1.0, -1.0)
            std = sigma * np.sqrt(u * (1.0 - u))
            min_s = np.minimum(u, 1.0 - u)
            du_int = signs * std
            rand_u = rng.random(K)
            a_low = (sigma ** 2) * (1.0 - u)
            du_low = np.where(rand_u < a_low / np.maximum(u + a_low, 1e-30), -u, a_low)
            v = 1.0 - u
            a_high = (sigma ** 2) * u
            du_high = np.where(rand_u < a_high / np.maximum(v + a_high, 1e-30), v, -a_high)
            near_low = (std > min_s) & (u <= 0.5) & (u > 0.0)
            near_high = (std > min_s) & (u > 0.5) & (u < 1.0)
            in_open = (u > 0.0) & (u < 1.0)
            du = np.where(near_low, du_low, np.where(near_high, du_high, du_int))
            u_new = np.where(in_open, u + du, u)
            u_new = np.where(u_new <= 1e-15, 0.0, np.where(u_new >= 1.0 - 1e-15, 1.0, u_new))
            cur_S[np.arange(K), i_col] = u_new * tot
            cur_S[np.arange(K), j_col] = (1.0 - u_new) * tot

        max_s = np.max(cur_S, axis=1)
        done = max_s >= 1.0 - 1e-14
        if np.any(done):
            d_idx = idx[done]
            winner_atom[d_idx] = np.argmax(cur_S[done], axis=1)
            active[d_idx] = False

        S[idx] = cur_S
        phi_det[idx] = cur_phi
        phi_vac[idx] = cur_vac

    n_timeouts = int(np.sum(active))
    if n_timeouts > 0:
        a_idx = np.where(active)[0]
        winner_atom[a_idx] = np.argmax(S[a_idx], axis=1)

    win_patch1 = (winner_atom >= 0) & (winner_atom < M_patch)
    return {
        "p_patch1": round(float(np.mean(win_patch1)), 4),
        "s_patch1_born": round(float(s_patch1), 4),
        "bias": round(float(np.mean(win_patch1) - s_patch1), 4),
        "mean_steps": round(float(np.mean(steps)), 2),
        "n_timeouts": n_timeouts,
        "tie_rate_5pct": tie_rate_5pct,
        "tie_rate_2pct": tie_rate_2pct,
        "mean_max_initial_atom_share": mean_max_initial_atom_share,
        "align_mode": align_mode,
        "lattice_size_per_patch": f"{Lx}x{Ly}={M_patch}",
        "total_lattice_atoms": M_tot,
    }


def main():
    out = {}

    # 1. Symmetry & Tie-Frequency Comparison at S = 0.50, theta = 0.45:
    #    Compare (a) flat-top uniform grid (16-fold tie symmetry),
    #            (b) perfectly centered Gaussian wave packet at (1.5, 1.5) (4-fold inner tie symmetry),
    #            (c) realistic random off-center wave-packet alignment delta_r0 ~ U(-0.5, +0.5) a_lat.
    out["alignment_symmetry_tie_comparison_s050_theta045"] = {
        "centered_flat_16fold": simulate_2d_crystal_lattice_aperture(
            s_patch1=0.50, Lx=4, Ly=4, theta=0.45, commit_rate=0.02, kappa_surf=2.5,
            surface_vac_eps=0.006, align_mode="centered_flat", n_trials=12000, seed=4101
        ),
        "centered_gaussian_4fold": simulate_2d_crystal_lattice_aperture(
            s_patch1=0.50, Lx=4, Ly=4, theta=0.45, commit_rate=0.02, kappa_surf=2.5,
            surface_vac_eps=0.006, align_mode="centered_gaussian", n_trials=12000, seed=4102
        ),
        "off_center_random_broken_symmetry": simulate_2d_crystal_lattice_aperture(
            s_patch1=0.50, Lx=4, Ly=4, theta=0.45, commit_rate=0.02, kappa_surf=2.5,
            surface_vac_eps=0.006, align_mode="off_center_random", n_trials=12000, seed=4103
        ),
    }

    # 2. Nine-angle sweep (S = 0.1 .. 0.9) on the 32-atom 2D crystal lattice with realistic
    #    random off-center wave-packet alignment (delta_r0 ~ U(-0.5, +0.5) a_lat):
    s_grid = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
    p_lat_th072 = []
    p_lat_th045 = []
    steps_lat_072 = []
    steps_lat_045 = []
    tie5_list = []
    tie2_list = []

    n_trials_sweep = 12000
    for i, s0 in enumerate(s_grid):
        r72 = simulate_2d_crystal_lattice_aperture(
            s_patch1=s0, Lx=4, Ly=4, theta=0.72, commit_rate=0.02, kappa_surf=2.5,
            surface_vac_eps=0.005, align_mode="off_center_random", n_trials=n_trials_sweep, seed=5000 + i
        )
        r45 = simulate_2d_crystal_lattice_aperture(
            s_patch1=s0, Lx=4, Ly=4, theta=0.45, commit_rate=0.02, kappa_surf=2.5,
            surface_vac_eps=0.005, align_mode="off_center_random", n_trials=n_trials_sweep, seed=6000 + i
        )
        p_lat_th072.append(r72["p_patch1"])
        p_lat_th045.append(r45["p_patch1"])
        steps_lat_072.append(r72["mean_steps"])
        steps_lat_045.append(r45["mean_steps"])
        tie5_list.append(r45["tie_rate_5pct"])
        tie2_list.append(r45["tie_rate_2pct"])

    fit_072 = fit_born_exponent(s_grid, p_lat_th072, n_trials_sweep)
    fit_045 = fit_born_exponent(s_grid, p_lat_th045, n_trials_sweep)

    out["2d_crystal_lattice_4x4_per_patch"] = {
        "s_grid": s_grid,
        "n_trials_per_point": n_trials_sweep,
        "lattice_geometry": (
            "Two 4x4 2D crystal lattice patches (16 atoms/patch, 32 atoms total) with random "
            "sub-lattice wave-packet alignment offset delta_r0 ~ U(-0.5, +0.5) a_lat, nearest-neighbor "
            "exciton hopping, and shared vacuum bus"
        ),
        "mean_tie_rate_5pct": round(float(np.mean(tie5_list)), 4),
        "mean_tie_rate_2pct": round(float(np.mean(tie2_list)), 4),
        "theta_072_p_patch1": p_lat_th072,
        "theta_072_rms_vs_born": fit_072["rms_vs_born"],
        "theta_072_fit": fit_072,
        "theta_045_p_patch1": p_lat_th045,
        "theta_045_rms_vs_born": fit_045["rms_vs_born"],
        "theta_045_fit": fit_045,
        "mean_steps_theta_072": steps_lat_072,
        "mean_steps_theta_045": steps_lat_045,
    }

    # 3. Ablation Study at S = 0.70, theta = 0.45 (testing what each postulated mechanism changes):
    out["ablation_study_s070_theta045"] = {
        "full_baseline": simulate_2d_crystal_lattice_aperture(
            s_patch1=0.70, Lx=4, Ly=4, theta=0.45, commit_rate=0.02, kappa_surf=2.5,
            kappa_pull=0.35, alpha_seed=0.15, surface_vac_eps=0.005,
            align_mode="off_center_random", n_trials=12000, seed=7101
        ),
        "no_dominant_atom_phase_pulling_kappa_pull_0": simulate_2d_crystal_lattice_aperture(
            s_patch1=0.70, Lx=4, Ly=4, theta=0.45, commit_rate=0.02, kappa_surf=2.5,
            kappa_pull=0.0, alpha_seed=0.15, surface_vac_eps=0.005,
            align_mode="off_center_random", n_trials=12000, seed=7101
        ),
        "no_surface_vonmises_kappa_surf_0": simulate_2d_crystal_lattice_aperture(
            s_patch1=0.70, Lx=4, Ly=4, theta=0.45, commit_rate=0.02, kappa_surf=0.0,
            kappa_pull=0.35, alpha_seed=0.15, surface_vac_eps=0.005,
            align_mode="off_center_random", n_trials=12000, seed=7101
        ),
        "no_premodulus_phase_seeding_alpha_seed_0": simulate_2d_crystal_lattice_aperture(
            s_patch1=0.70, Lx=4, Ly=4, theta=0.45, commit_rate=0.02, kappa_surf=2.5,
            kappa_pull=0.35, alpha_seed=0.0, surface_vac_eps=0.005,
            align_mode="off_center_random", n_trials=12000, seed=7101
        ),
        "no_vacuum_drift_surface_vac_eps_0": simulate_2d_crystal_lattice_aperture(
            s_patch1=0.70, Lx=4, Ly=4, theta=0.45, commit_rate=0.02, kappa_surf=2.5,
            kappa_pull=0.35, alpha_seed=0.15, surface_vac_eps=0.0,
            align_mode="off_center_random", n_trials=12000, seed=7101
        ),
        "no_off_center_alignment_flat_16fold": simulate_2d_crystal_lattice_aperture(
            s_patch1=0.70, Lx=4, Ly=4, theta=0.45, commit_rate=0.02, kappa_surf=2.5,
            kappa_pull=0.35, alpha_seed=0.25, surface_vac_eps=0.005,
            align_mode="centered_flat", n_trials=12000, seed=7101
        ),
    }

    with open(CODE_DIR / "crystal_lattice_results.json", "w") as f:
        json.dump(out, f, indent=2)

    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
