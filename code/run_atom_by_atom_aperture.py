#!/usr/bin/env python3
"""
Atom-by-Atom Pre-Modulus Phase-Coupled Microsimulation of a Detector Aperture
with Realistic Off-Center Wave-Packet Alignment (Tie-Breaking)
for Aperture Article 0002 (Revision 2.3 — 2026-10-29)

Physical Architecture (John M. Bramble, MD):
1. Individual Atoms in the Aperture (Not Aperture-Integrated Ports):
   Instead of collapsing each detector region/fringe into a single macroscopic port variable S
   (the "integral of the aperture"), the detector aperture consists of M_total individual
   absorbing atoms/electrons m = 1..M_total distributed across the interference pattern.
   Because any physical fringe illuminates many individual atoms (e.g., M_1 atoms in Region 1
   and M_2 atoms in Region 2), each individual atom initially holds a small amplitude
   |psi_m| << sqrt(theta), so NO single atom starts above the bandgap threshold gate theta
   (even when theta = 0.45!).

2. Realistic Off-Center Wave-Packet Alignment (Breaking Artificial Multi-Atom Ties):
   In idealized simulations where the detector array is perfectly centered on the incident
   wave packet or modeled as a flat-top uniform profile, multiple symmetric atoms in the bright
   spot start with nearly identical shares, producing an uncharacteristically high frequency of
   ties. Introducing a realistic random sub-spacing spatial offset delta_x0 ~ U(-0.5, +0.5) a_lat
   of the incident wave packet relative to the discrete detector array breaks multi-atom tie
   symmetry at t = 0 and significantly reduces tie frequency in the bright spot.

3. Surface-Conditioned Quantum Vacuum Phase & Three-Way Pre-Modulus Phase Coupling:
   - At the condensed-matter detector surface (or under an applied detector bias field E_bias),
     the local virtual vacuum modes V_m = |V_m| exp(i * phi_vac,m) have an altered, surface-
     conditioned phase distribution around a local surface reference phase phi_surf, rather than
     free-space uniform white noise.
   - During capture, the incident particle wave A_m = |A_m| exp(i * phi_inc,m), the individual
     detector atom's internal clock phase exp(i * phi_det,m), and the surface-conditioned
     virtual vacuum phase exp(i * phi_vac,m) couple BEFORE any modulus squared is evaluated.
"""

from pathlib import Path
import json
import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import chi2

CODE_DIR = Path(__file__).resolve().parent


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


def run_atom_by_atom_aperture_sim(
    s_region1=0.5,
    n_atoms_region1=10,
    n_atoms_region2=10,
    theta=0.45,
    sigma=0.25,
    commit_rate=0.02,
    surface_vac_eps=0.0,
    align_mode="off_center_random",
    n_trials=20000,
    max_steps=8000,
    seed=20261029,
):
    """
    Simulates an aperture of M_total = n_atoms_region1 + n_atoms_region2 individual atoms
    with random sub-spacing spatial offset delta_x0 ~ U(-0.5, +0.5) of the incident wave packet
    relative to the detector array (`align_mode = "off_center_random"`).
    """
    rng = np.random.default_rng(seed)
    M1 = n_atoms_region1
    M2 = n_atoms_region2
    M_tot = M1 + M2

    # 1. Pre-modulus initial complex phases at each individual atom m = 0..M_tot-1:
    phi_inc = rng.uniform(0.0, 2.0 * np.pi, size=(n_trials, M_tot))
    phi_det = rng.uniform(0.0, 2.0 * np.pi, size=(n_trials, M_tot))
    phi_surf = rng.uniform(0.0, 2.0 * np.pi, size=(n_trials, 1))
    phi_vac = phi_surf + rng.normal(0.0, 0.5, size=(n_trials, M_tot))

    # 2. Spatial wave-packet envelope with realistic random alignment offset delta_x0 ~ U(-0.5, +0.5):
    if align_mode == "off_center_random" and M1 > 1:
        x1 = np.arange(M1, dtype=float) - (M1 - 1) / 2.0
        off1 = rng.uniform(-0.5, 0.5, size=(n_trials, 1))
        waist1 = 0.55 * M1
        env1 = np.exp(-((x1[None, :] - off1) ** 2) / (2.0 * (waist1 ** 2)))
    else:
        env1 = np.ones((n_trials, M1))

    if align_mode == "off_center_random" and M2 > 1:
        x2 = np.arange(M2, dtype=float) - (M2 - 1) / 2.0
        off2 = rng.uniform(-0.5, 0.5, size=(n_trials, 1))
        waist2 = 0.55 * M2
        env2 = np.exp(-((x2[None, :] - off2) ** 2) / (2.0 * (waist2 ** 2)))
    else:
        env2 = np.ones((n_trials, M2))

    # 3. Pre-modulus 3-way phase alignment combined with spatial wave-packet envelope:
    align_score1 = env1 * (1.0 + 0.20 * np.cos(phi_det[:, :M1] - 0.5 * (phi_inc[:, :M1] + phi_vac[:, :M1])))
    align_score2 = env2 * (1.0 + 0.20 * np.cos(phi_det[:, M1:] - 0.5 * (phi_inc[:, M1:] + phi_vac[:, M1:])))
    w1 = align_score1 / np.sum(align_score1, axis=1, keepdims=True)
    w2 = align_score2 / np.sum(align_score2, axis=1, keepdims=True)

    if M1 > 1:
        w1_sorted = np.sort(w1, axis=1)[:, ::-1]
        rel_gap = (w1_sorted[:, 0] - w1_sorted[:, 1]) / np.maximum(w1_sorted[:, 0], 1e-30)
        tie_rate_5pct = round(float(np.mean(rel_gap < 0.05)), 4)
        tie_rate_2pct = round(float(np.mean(rel_gap < 0.02)), 4)
        max_init_atom = round(float(np.mean(s_region1 * w1_sorted[:, 0])), 4)
    else:
        tie_rate_5pct = 0.0
        tie_rate_2pct = 0.0
        max_init_atom = round(float(s_region1), 4)

    S = np.zeros((n_trials, M_tot), dtype=float)
    S[:, :M1] = s_region1 * w1
    S[:, M1:] = (1.0 - s_region1) * w2

    active = np.ones(n_trials, dtype=bool)
    winner_atom = np.full(n_trials, -1, dtype=int)
    steps = np.zeros(n_trials, dtype=int)

    for _ in range(1, max_steps + 1):
        idx = np.where(active)[0]
        if len(idx) == 0:
            break
        cur_S = S[idx]
        cur_phi = phi_det[idx]
        cur_vac = phi_vac[idx]
        K = len(idx)
        steps[idx] += 1

        # Check individual-atom threshold gate: does any SINGLE atom m have s_m > theta?
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

        # Dominant atom pulls local surface vacuum phase:
        dom_atom = np.argmax(cur_S, axis=1)
        dom_phi = cur_phi[np.arange(K), dom_atom][:, None]
        cur_vac = cur_vac + 0.3 * np.sin(dom_phi - cur_vac) + rng.normal(0.0, 0.4, size=(K, M_tot))
        cur_phi = cur_phi + 0.15 * np.sin(cur_vac - cur_phi) + rng.uniform(-np.pi, np.pi, size=(K, M_tot))

        # Pairwise pre-modulus complex interference exchange between active individual atoms:
        perm = rng.permuted(np.tile(np.arange(M_tot), (K, 1)), axis=1)
        for p in range(0, M_tot - 1, 2):
            i_col = perm[:, p]
            j_col = perm[:, p + 1]
            si = cur_S[np.arange(K), i_col]
            sj = cur_S[np.arange(K), j_col]
            tot = si + sj
            act_pair = tot > 0.0
            u = np.where(act_pair, si / np.maximum(tot, 1e-30), 0.0)

            if surface_vac_eps != 0.0:
                du_vac = surface_vac_eps * np.sqrt(u * (1.0 - u)) * (np.sqrt(1.0 - u) - np.sqrt(u))
                u = np.where((u > 0.0) & (u < 1.0), u + du_vac, u)

            dphi = (cur_phi[np.arange(K), i_col] - cur_vac[np.arange(K), i_col]) - (
                cur_phi[np.arange(K), j_col] - cur_vac[np.arange(K), j_col]
            )
            signs = np.where(np.cos(dphi) >= 0.0, 1.0, -1.0)

            std = sigma * np.sqrt(u * (1.0 - u))
            min_side = np.minimum(u, 1.0 - u)
            du_int = signs * std

            rand_u = rng.random(K)
            a_low = (sigma ** 2) * (1.0 - u)
            p0_low = a_low / np.maximum(u + a_low, 1e-30)
            du_low = np.where(rand_u < p0_low, -u, a_low)

            v = 1.0 - u
            a_high = (sigma ** 2) * u
            p1_high = a_high / np.maximum(v + a_high, 1e-30)
            du_high = np.where(rand_u < p1_high, v, -a_high)

            near_low = (std > min_side) & (u <= 0.5) & (u > 0.0)
            near_high = (std > min_side) & (u > 0.5) & (u < 1.0)
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

    win_region1 = (winner_atom >= 0) & (winner_atom < M1)
    return {
        "p_region1": round(float(np.mean(win_region1)), 4),
        "s_region1_born": round(float(s_region1), 4),
        "bias": round(float(np.mean(win_region1) - s_region1), 4),
        "mean_steps": round(float(np.mean(steps)), 2),
        "n_timeouts": n_timeouts,
        "M1_atoms": M1,
        "M2_atoms": M2,
        "align_mode": align_mode,
        "tie_rate_5pct": tie_rate_5pct,
        "tie_rate_2pct": tie_rate_2pct,
        "max_initial_atom_share": max_init_atom,
    }


def main():
    results = {}

    # 1. Compare Aperture-Integrated (1 giant "atom" per region) vs. Multi-Atom Aperture
    #    with realistic off-center wave-packet alignment at theta = 0.45!
    results["theta_045_aperture_resolution_test"] = {
        "lumped_1_bright_vs_9_dim_no_vac": run_atom_by_atom_aperture_sim(
            s_region1=0.50, n_atoms_region1=1, n_atoms_region2=9, theta=0.45,
            commit_rate=0.02, surface_vac_eps=0.0, align_mode="off_center_random", seed=1001
        ),
        "resolved_5_bright_vs_15_dim_no_vac": run_atom_by_atom_aperture_sim(
            s_region1=0.50, n_atoms_region1=5, n_atoms_region2=15, theta=0.45,
            commit_rate=0.02, surface_vac_eps=0.0, align_mode="off_center_random", seed=1002
        ),
        "resolved_5_bright_vs_15_dim_with_surface_vac": run_atom_by_atom_aperture_sim(
            s_region1=0.50, n_atoms_region1=5, n_atoms_region2=15, theta=0.45,
            commit_rate=0.02, surface_vac_eps=0.008, align_mode="off_center_random", seed=1003
        ),
        "resolved_10_bright_vs_10_dim_no_vac": run_atom_by_atom_aperture_sim(
            s_region1=0.50, n_atoms_region1=10, n_atoms_region2=10, theta=0.45,
            commit_rate=0.02, surface_vac_eps=0.0, align_mode="off_center_random", seed=1004
        ),
        "resolved_10_bright_vs_10_dim_centered_flat_comparison": run_atom_by_atom_aperture_sim(
            s_region1=0.50, n_atoms_region1=10, n_atoms_region2=10, theta=0.45,
            commit_rate=0.02, surface_vac_eps=0.0, align_mode="centered_flat", seed=1004
        ),
    }

    # 2. Nine-angle sweep (S in 0.1 .. 0.9) on a 20-atom detector aperture with realistic
    #    off-center wave-packet alignment (M1 = 10 atoms in Port 1, M2 = 10 atoms in Port 2)
    s_grid = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
    p_20atom_disc = []
    p_20atom_th072_novac = []
    p_20atom_th072_vac = []
    p_20atom_th045_novac = []
    p_20atom_th045_vac = []

    n_trials_sweep = 12000
    for i, s0 in enumerate(s_grid):
        r_d = run_atom_by_atom_aperture_sim(
            s_region1=s0, n_atoms_region1=10, n_atoms_region2=10, theta=1.0, commit_rate=0.0,
            surface_vac_eps=0.0, align_mode="off_center_random", n_trials=n_trials_sweep, seed=2000 + i
        )
        r_72_0 = run_atom_by_atom_aperture_sim(
            s_region1=s0, n_atoms_region1=10, n_atoms_region2=10, theta=0.72, commit_rate=0.02,
            surface_vac_eps=0.0, align_mode="off_center_random", n_trials=n_trials_sweep, seed=2100 + i
        )
        r_72_v = run_atom_by_atom_aperture_sim(
            s_region1=s0, n_atoms_region1=10, n_atoms_region2=10, theta=0.72, commit_rate=0.02,
            surface_vac_eps=0.006, align_mode="off_center_random", n_trials=n_trials_sweep, seed=2200 + i
        )
        r_45_0 = run_atom_by_atom_aperture_sim(
            s_region1=s0, n_atoms_region1=10, n_atoms_region2=10, theta=0.45, commit_rate=0.02,
            surface_vac_eps=0.0, align_mode="off_center_random", n_trials=n_trials_sweep, seed=2300 + i
        )
        r_45_v = run_atom_by_atom_aperture_sim(
            s_region1=s0, n_atoms_region1=10, n_atoms_region2=10, theta=0.45, commit_rate=0.02,
            surface_vac_eps=0.006, align_mode="off_center_random", n_trials=n_trials_sweep, seed=2400 + i
        )
        p_20atom_disc.append(r_d["p_region1"])
        p_20atom_th072_novac.append(r_72_0["p_region1"])
        p_20atom_th072_vac.append(r_72_v["p_region1"])
        p_20atom_th045_novac.append(r_45_0["p_region1"])
        p_20atom_th045_vac.append(r_45_v["p_region1"])

    results["twenty_atom_aperture_sweep"] = {
        "s_grid": s_grid,
        "M1_atoms_per_port": 10,
        "M2_atoms_per_port": 10,
        "align_mode": "off_center_random",
        "n_trials_per_point": n_trials_sweep,
        "discrete_theta_1": {
            "p_region1": p_20atom_disc,
            "fit": fit_born_exponent(s_grid, p_20atom_disc, n_trials_sweep),
            "rms_vs_born": fit_born_exponent(s_grid, p_20atom_disc, n_trials_sweep)["rms_vs_born"],
        },
        "gated_theta_072_no_vac": {
            "p_region1": p_20atom_th072_novac,
            "fit": fit_born_exponent(s_grid, p_20atom_th072_novac, n_trials_sweep),
            "rms_vs_born": fit_born_exponent(s_grid, p_20atom_th072_novac, n_trials_sweep)["rms_vs_born"],
        },
        "gated_theta_072_vac_0006": {
            "p_region1": p_20atom_th072_vac,
            "fit": fit_born_exponent(s_grid, p_20atom_th072_vac, n_trials_sweep),
            "rms_vs_born": fit_born_exponent(s_grid, p_20atom_th072_vac, n_trials_sweep)["rms_vs_born"],
        },
        "gated_theta_045_no_vac": {
            "p_region1": p_20atom_th045_novac,
            "fit": fit_born_exponent(s_grid, p_20atom_th045_novac, n_trials_sweep),
            "rms_vs_born": fit_born_exponent(s_grid, p_20atom_th045_novac, n_trials_sweep)["rms_vs_born"],
        },
        "gated_theta_045_vac_0006": {
            "p_region1": p_20atom_th045_vac,
            "fit": fit_born_exponent(s_grid, p_20atom_th045_vac, n_trials_sweep),
            "rms_vs_born": fit_born_exponent(s_grid, p_20atom_th045_vac, n_trials_sweep)["rms_vs_born"],
        },
    }

    out_path = CODE_DIR / "atom_by_atom_aperture_results.json"
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)

    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
