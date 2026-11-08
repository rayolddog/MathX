#!/usr/bin/env python3
"""
Revision 2.9 Verification Suite (`run_rev29_verification_suite.py`)
for Aperture Article 0002 (2026-11-02)

Implements and verifies the mathematical and numerical refinements from the
Claude Fable 5.1 and ChatGPT Codex (Sol) reviews:
  1. 2pi-periodic, shift-invariant three-phase seed:
       C_m = 0.5 * [cos(phi_det - phi_inc) + cos(phi_det - phi_vac)],
       F_m = 1 + alpha_seed * C_m,  tilde_psi_m(0) = A_m * sqrt(F_m).
  2. Whole-quantum absorption (theta = 1, eps_vac = 0) 9-angle sweeps on the
     20-atom aperture and 32-atom 2D crystal lattice under the true physical
     timescale hierarchy: fast uniform QED vacuum phases (kappa_surf = 0,
     kappa_pull = 0) + slow detector electron phase drift + antisymmetric
     sin(Delta_phi_mn) sign, with strict censoring discipline (zero forced argmax).
  3. Exact analytical two-site eps_vac drift evaluation (q'(u) proportional to
     exp[-(4*eps_vac/sigma^2)*(sqrt(u) + sqrt(1-u))]), quantifying why eps_vac
     is dropped in favor of whole-quantum (theta = 1) registration.
  4. Backward-equation threshold balance residual R_r = b_r + R_r - R * S_r,
     nonlocal exact compensation b_m = R * s_m - r_m, and step-size refinement
     (Delta_t = 1, 1/2, 1/4) for the 5-vs-15 unequal-subdivision comparator at
     theta = 0.45.
"""

import json
from pathlib import Path
import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.optimize import curve_fit
from scipy.stats import binomtest, chi2

CODE_DIR = Path(__file__).resolve().parent


def periodic_phase_score(
    phi_det: np.ndarray,
    phi_inc: np.ndarray,
    phi_vac: np.ndarray,
    alpha_seed: float = 0.15,
) -> np.ndarray:
    """Strictly 2pi-periodic and shift-invariant three-phase seed factor F_m > 0."""
    c_m = 0.5 * (np.cos(phi_det - phi_inc) + np.cos(phi_det - phi_vac))
    return 1.0 + alpha_seed * c_m


def check_seed_invariances() -> dict:
    """Verifies 2pi periodicity in each phase and common-shift invariance."""
    rng = np.random.default_rng(20261102)
    pd = rng.uniform(-np.pi, np.pi, size=64)
    pi = rng.uniform(-np.pi, np.pi, size=64)
    pv = rng.uniform(-np.pi, np.pi, size=64)
    shift = 1.234567

    f0 = periodic_phase_score(pd, pi, pv, alpha_seed=0.15)
    err_inc_2pi = float(np.max(np.abs(periodic_phase_score(pd, pi + 2.0 * np.pi, pv, 0.15) - f0)))
    err_vac_2pi = float(np.max(np.abs(periodic_phase_score(pd, pi, pv - 2.0 * np.pi, 0.15) - f0)))
    err_det_2pi = float(np.max(np.abs(periodic_phase_score(pd + 2.0 * np.pi, pi, pv, 0.15) - f0)))
    err_common_shift = float(
        np.max(np.abs(periodic_phase_score(pd + shift, pi + shift, pv + shift, 0.15) - f0))
    )
    return {
        "max_error_inc_plus_2pi": err_inc_2pi,
        "max_error_vac_minus_2pi": err_vac_2pi,
        "max_error_det_plus_2pi": err_det_2pi,
        "max_error_common_phase_shift": err_common_shift,
        "min_F_m": float(np.min(f0)),
        "max_F_m": float(np.max(f0)),
    }


def apply_bounded_matching_step_2d(
    s: np.ndarray,
    ic: np.ndarray,
    jc: np.ndarray,
    signs: np.ndarray,
    sigma: float,
    rng: np.random.Generator,
) -> None:
    """
    Applies the exact two-point boundary-spliced neutral step simultaneously across
    a set of disjoint pairs (ic, jc) of shape (K, n_pairs). Zero interior clipping.
    """
    K, n_pairs = ic.shape
    rows = np.arange(K)[:, None]
    sm = s[rows, ic]
    sn = s[rows, jc]
    tot = sm + sn
    act = tot > 1e-15
    u = np.where(act, sm / np.maximum(tot, 1e-30), 0.0)
    v = 1.0 - u

    std = sigma * np.sqrt(u * v)
    du_int = signs * std

    rand_u = rng.random((K, n_pairs))
    a_low = (sigma ** 2) * v
    p0_low = a_low / np.maximum(u + a_low, 1e-30)
    du_low = np.where(rand_u < p0_low, -u, a_low)

    a_high = (sigma ** 2) * u
    p1_high = a_high / np.maximum(v + a_high, 1e-30)
    du_high = np.where(rand_u < p1_high, v, -a_high)

    du = np.where(u <= v, np.where(std > u, du_low, du_int), np.where(std > v, du_high, du_int))
    u_new = u + np.where(act, du, 0.0)
    if np.any((u_new < -1e-14) | (u_new > 1.0 + 1e-14)):
        raise RuntimeError("Pair step left [0, 1] simplex interval.")
    u_new = np.where(u_new <= 1e-15, 0.0, np.where(u_new >= 1.0 - 1e-15, 1.0, u_new))

    s[rows, ic] = np.where(act, tot * u_new, sm)
    s[rows, jc] = np.where(act, tot * (1.0 - u_new), sn)


def build_4x4_patch_bonds() -> np.ndarray:
    xy = np.array([(x, y) for x in range(4) for y in range(4)])
    bonds = [
        (i, j)
        for i in range(16)
        for j in range(i + 1, 16)
        if 0 < np.sum((xy[i] - xy[j]) ** 2) <= 2
    ]
    return np.array(bonds, dtype=int)


PATCH_BONDS_4X4 = build_4x4_patch_bonds()


def simulate_whole_quantum_trial_batch(
    M1: int,
    M2: int,
    s_port1: float,
    sigma: float = 0.25,
    use_lattice_tierA: bool = False,
    n_trials: int = 12000,
    max_steps: int = 12000,
    seed: int = 20261102,
) -> dict:
    """
    Simulates whole-quantum vertex fixation (theta = 1, eps_vac = 0) under the physical
    timescale hierarchy:
      - Periodic 3-phase seed at t = 0 with off-center spatial wave-packet envelope.
      - Slow detector electron phase drift (sigma_det = 0.02 rad/step).
      - Fast, uniform, share-blind quantum vacuum phase draws phi_vac ~ U(-pi, pi) (kappa_surf = 0).
      - Antisymmetric transfer direction sign(sin(Delta_phi_mn)) drawn freshly at each bond.
      - Strict censoring discipline: never assigns an argmax winner to an unresolved run.
    """
    rng = np.random.default_rng(seed)
    M_tot = M1 + M2

    if use_lattice_tierA and M1 == 16 and M2 == 16:
        coords = np.array([(ix, iy) for ix in range(4) for iy in range(4)], dtype=float)
        center = np.array([1.5, 1.5])
        d1 = rng.uniform(-0.5, 0.5, size=(n_trials, 1, 2))
        d2 = rng.uniform(-0.5, 0.5, size=(n_trials, 1, 2))
        dist2_1 = np.sum((coords[None, :, :] - (center[None, None, :] + d1)) ** 2, axis=2)
        dist2_2 = np.sum((coords[None, :, :] - (center[None, None, :] + d2)) ** 2, axis=2)
        env1 = np.exp(-dist2_1 / (2.0 * (1.4 ** 2)))
        env2 = np.exp(-dist2_2 / (2.0 * (1.4 ** 2)))
        phi_inc1 = 0.6 * coords[None, :, 0] + 0.4 * coords[None, :, 1]
        phi_inc2 = 0.6 * coords[None, :, 0] + 0.4 * coords[None, :, 1] + 0.4
    else:
        x1 = np.arange(M1, dtype=float)[None, :]
        x2 = np.arange(M2, dtype=float)[None, :]
        c1 = 0.5 * (M1 - 1) + rng.uniform(-0.5, 0.5, size=(n_trials, 1))
        c2 = 0.5 * (M2 - 1) + rng.uniform(-0.5, 0.5, size=(n_trials, 1))
        w1 = max(0.55 * M1, 1.0)
        w2 = max(0.55 * M2, 1.0)
        env1 = np.exp(-((x1 - c1) ** 2) / (2.0 * (w1 ** 2)))
        env2 = np.exp(-((x2 - c2) ** 2) / (2.0 * (w2 ** 2)))
        phi_inc1 = 0.45 * x1
        phi_inc2 = 0.45 * x2 + 0.3

    phi_det = rng.uniform(-np.pi, np.pi, size=(n_trials, M_tot))
    phi_vac = rng.uniform(-np.pi, np.pi, size=(n_trials, M_tot))

    f1 = env1 * periodic_phase_score(phi_det[:, :M1], phi_inc1, phi_vac[:, :M1], alpha_seed=0.15)
    f2 = env2 * periodic_phase_score(phi_det[:, M1:], phi_inc2, phi_vac[:, M1:], alpha_seed=0.15)
    f1 /= np.sum(f1, axis=1, keepdims=True)
    f2 /= np.sum(f2, axis=1, keepdims=True)

    s = np.zeros((n_trials, M_tot), dtype=float)
    s[:, :M1] = s_port1 * f1
    s[:, M1:] = (1.0 - s_port1) * f2

    active = np.ones(n_trials, dtype=bool)
    steps = np.zeros(n_trials, dtype=int)
    n_bonds = len(PATCH_BONDS_4X4)

    for step_idx in range(1, max_steps + 1):
        idx = np.where(active)[0]
        if len(idx) == 0:
            break
        K = len(idx)
        steps[idx] += 1
        cur_s = s[idx]
        cur_det = phi_det[idx] + rng.normal(0.0, 0.02, size=(K, M_tot))
        rows = np.arange(K)[:, None]

        # Optional Tier-A intra-patch lattice exchange every 3rd step (1 bond per patch, preserving S_r)
        if use_lattice_tierA and (step_idx % 3 == 0):
            b_sel = rng.integers(0, n_bonds, size=K)
            a1 = PATCH_BONDS_4X4[b_sel, 0]
            b1 = PATCH_BONDS_4X4[b_sel, 1]
            ic_b = np.column_stack([a1, a1 + M1])
            jc_b = np.column_stack([b1, b1 + M1])
            vac_b = rng.uniform(-np.pi, np.pi, size=(K, 2, 2))
            dphi_b = (cur_det[rows, ic_b] - vac_b[:, :, 0]) - (cur_det[rows, jc_b] - vac_b[:, :, 1])
            signs_b = np.where(np.sin(dphi_b) >= 0.0, 1.0, -1.0)
            apply_bounded_matching_step_2d(cur_s, ic_b, jc_b, signs_b, sigma, rng)

        # Tier-B global random-matching exchange with fresh fast uniform vacuum phases
        perm = np.argsort(rng.random((K, M_tot)), axis=1)
        ic = perm[:, 0::2]
        jc = perm[:, 1::2]
        cur_vac = rng.uniform(-np.pi, np.pi, size=(K, M_tot))
        dphi = (cur_det[rows, ic] - cur_vac[rows, ic]) - (cur_det[rows, jc] - cur_vac[rows, jc])
        signs = np.where(np.sin(dphi) >= 0.0, 1.0, -1.0)
        apply_bounded_matching_step_2d(cur_s, ic, jc, signs, sigma, rng)

        s[idx] = cur_s
        phi_det[idx] = (cur_det + np.pi) % (2.0 * np.pi) - np.pi

        fixed_now = np.max(cur_s, axis=1) >= 1.0 - 1e-12
        if np.any(fixed_now):
            active[idx[fixed_now]] = False

    n_unresolved = int(np.sum(active))
    max_cons_err = float(np.max(np.abs(np.sum(s, axis=1) - 1.0)))
    resolved_mask = ~active
    winners = np.argmax(s[resolved_mask], axis=1)
    n_port1_wins = int(np.sum(winners < M1))
    p_low = n_port1_wins / n_trials
    p_high = (n_port1_wins + n_unresolved) / n_trials
    p_point = p_low if n_unresolved == 0 else None

    return {
        "n_trials": n_trials,
        "n_unresolved": n_unresolved,
        "n_port1_wins": n_port1_wins,
        "p_port1": round(p_point, 4) if p_point is not None else None,
        "p_interval": [round(p_low, 4), round(p_high, 4)],
        "mean_steps": round(float(np.mean(steps)), 2),
        "max_conservation_error": max_cons_err,
    }


def run_nine_angle_whole_quantum_sweep(
    M1: int,
    M2: int,
    sigma: float,
    use_lattice_tierA: bool,
    n_trials: int = 12000,
    seed_base: int = 20261102,
) -> dict:
    s_grid = [round(0.1 * k, 1) for k in range(1, 10)]
    p_obs = []
    mean_steps = []
    total_unresolved = 0
    max_cons_err = 0.0

    for idx_s, s_val in enumerate(s_grid):
        res = simulate_whole_quantum_trial_batch(
            M1=M1,
            M2=M2,
            s_port1=s_val,
            sigma=sigma,
            use_lattice_tierA=use_lattice_tierA,
            n_trials=n_trials,
            max_steps=12000,
            seed=seed_base + 100 * idx_s,
        )
        total_unresolved += res["n_unresolved"]
        max_cons_err = max(max_cons_err, res["max_conservation_error"])
        p_obs.append(res["p_port1"])
        mean_steps.append(res["mean_steps"])

    if total_unresolved > 0:
        raise RuntimeError(f"Unresolved trials encountered: {total_unresolved}")

    s_arr = np.array(s_grid)
    p_arr = np.array(p_obs, dtype=float)

    def model(s, p):
        sp = np.power(s, p / 2.0)
        qp = np.power(1.0 - s, p / 2.0)
        return sp / (sp + qp)

    sigma_y = np.sqrt(s_arr * (1.0 - s_arr) / n_trials)
    popt, pcov = curve_fit(model, s_arr, p_arr, p0=[2.0], sigma=sigma_y, absolute_sigma=True)
    p_hat = float(popt[0])
    se_p = float(np.sqrt(pcov[0, 0]))
    chi2_val = float(np.sum(n_trials * ((p_arr - s_arr) ** 2) / (s_arr * (1.0 - s_arr))))
    pval = float(1.0 - chi2.cdf(chi2_val, df=len(s_arr)))
    rms = float(np.sqrt(np.mean((p_arr - s_arr) ** 2)))

    return {
        "sigma": sigma,
        "s_grid": s_grid,
        "p_observed": p_obs,
        "mean_steps": mean_steps,
        "total_unresolved": total_unresolved,
        "max_conservation_error": max_cons_err,
        "p_hat": round(p_hat, 4),
        "se_p": round(se_p, 4),
        "ci_95": [round(p_hat - 1.96 * se_p, 3), round(p_hat + 1.96 * se_p, 3)],
        "chi2_born": round(chi2_val, 2),
        "p_value": round(pval, 4),
        "rms_error_pp": round(100.0 * rms, 2),
    }


def compute_analytical_eps_vac_table(sigma: float = 0.25) -> list[dict]:
    """
    Computes the exact two-site exit probability q(u) under the pairwise drift
    Delta_u_vac = eps_vac * sqrt(u(1-u)) * (sqrt(1-u) - sqrt(u)) and variance sigma^2 u(1-u):
        q'(u) proportional to exp[-(4 * eps_vac / sigma^2) * (sqrt(u) + sqrt(1-u))].
    """
    u_grid = np.linspace(0.0, 1.0, 20001)
    s_inputs = np.array([0.1 * k for k in range(1, 10)])
    rows = []

    def model(s, p):
        sp = np.power(s, p / 2.0)
        qp = np.power(1.0 - s, p / 2.0)
        return sp / (sp + qp)

    for eps_vac in (0.0, 0.005, 0.006, 0.008, 0.012, 0.018):
        if eps_vac == 0.0:
            q_vals = s_inputs.copy()
        else:
            integrand = np.exp(-(4.0 * eps_vac / (sigma ** 2)) * (np.sqrt(u_grid) + np.sqrt(1.0 - u_grid)))
            cdf = cumulative_trapezoid(integrand, u_grid, initial=0.0)
            cdf /= cdf[-1]
            q_vals = np.interp(s_inputs, u_grid, cdf)
        popt, _ = curve_fit(model, s_inputs, q_vals, p0=[2.0])
        rows.append(
            {
                "sigma": sigma,
                "eps_vac": eps_vac,
                "q_010": round(float(q_vals[0]), 4),
                "q_030": round(float(q_vals[2]), 4),
                "q_070": round(float(q_vals[6]), 4),
                "fitted_p": round(float(popt[0]), 3),
            }
        )
    return rows


def verify_threshold_balance_and_step_refinement() -> dict:
    """
    1. Verifies Sol's backward-equation balance identity R_r = b_r + R_r - R * S_r
       on random simplex points and confirms exact cancellation under b_m = R*s_m - r_m.
    2. Runs the 5-vs-15 unequal-subdivision threshold comparator at S = 0.50, theta = 0.45
       across step sizes dt in {1.0, 0.5, 0.25} with scaled sigma_h = sigma*sqrt(dt)
       and lambda_h = lambda*dt (10,000 trials per step size).
    """
    rng = np.random.default_rng(20261103)
    M1, M2 = 5, 15
    M_tot = M1 + M2
    theta = 0.45
    lam = 0.02
    x = rng.exponential(scale=1.0, size=(500, M_tot))
    s_rand = x / np.sum(x, axis=1, keepdims=True)
    s_rand[:, 0] = 0.55
    s_rand[:, 1:] *= 0.45 / np.sum(s_rand[:, 1:], axis=1, keepdims=True)

    r_m = lam * s_rand * (s_rand > theta)
    R_tot = np.sum(r_m, axis=1, keepdims=True)
    R_1 = np.sum(r_m[:, :M1], axis=1)
    S_1 = np.sum(s_rand[:, :M1], axis=1)

    uncomp_residual = R_1 - R_tot[:, 0] * S_1
    b_comp = R_tot * s_rand - r_m
    b_1 = np.sum(b_comp[:, :M1], axis=1)
    comp_residual = b_1 + R_1 - R_tot[:, 0] * S_1

    refinement_rows = []
    for dt in (1.0, 0.5, 0.25):
        sigma_h = 0.25 * np.sqrt(dt)
        lam_h = 0.02 * dt
        n_trials = 10000
        max_steps = int(3500 / dt)

        x1 = np.arange(M1, dtype=float)[None, :]
        x2 = np.arange(M2, dtype=float)[None, :]
        c1 = 0.5 * (M1 - 1) + rng.uniform(-0.5, 0.5, size=(n_trials, 1))
        c2 = 0.5 * (M2 - 1) + rng.uniform(-0.5, 0.5, size=(n_trials, 1))
        env1 = np.exp(-((x1 - c1) ** 2) / (2.0 * (max(0.55 * M1, 1.0) ** 2)))
        env2 = np.exp(-((x2 - c2) ** 2) / (2.0 * (max(0.55 * M2, 1.0) ** 2)))
        pd = rng.uniform(-np.pi, np.pi, size=(n_trials, M_tot))
        pv = rng.uniform(-np.pi, np.pi, size=(n_trials, M_tot))
        f1 = env1 * periodic_phase_score(pd[:, :M1], 0.45 * x1, pv[:, :M1], 0.15)
        f2 = env2 * periodic_phase_score(pd[:, M1:], 0.45 * x2 + 0.3, pv[:, M1:], 0.15)
        f1 /= np.sum(f1, axis=1, keepdims=True)
        f2 /= np.sum(f2, axis=1, keepdims=True)

        s = np.zeros((n_trials, M_tot), dtype=float)
        s[:, :M1] = 0.50 * f1
        s[:, M1:] = 0.50 * f2

        active = np.ones(n_trials, dtype=bool)
        winners = np.full(n_trials, -1, dtype=int)

        for _ in range(max_steps):
            idx = np.where(active)[0]
            if len(idx) == 0:
                break
            K = len(idx)
            cur_s = s[idx]

            rates = lam_h * cur_s * (cur_s > theta)
            r_tot = np.sum(rates, axis=1)
            p_commit = 1.0 - np.exp(-r_tot)
            commit_now = (r_tot > 0.0) & (rng.random(K) < p_commit)
            if np.any(commit_now):
                c_idx = idx[commit_now]
                c_rates = rates[commit_now]
                c_rtot = r_tot[commit_now]
                cum_probs = np.cumsum(c_rates / c_rtot[:, None], axis=1)
                u_pick = rng.random(len(c_idx))[:, None]
                winners[c_idx] = np.argmax(u_pick < cum_probs, axis=1)
                active[c_idx] = False

            rem_mask = ~commit_now
            if not np.any(rem_mask):
                continue
            idx_rem = idx[rem_mask]
            K_rem = len(idx_rem)
            cur_s_rem = s[idx_rem]
            perm = np.argsort(rng.random((K_rem, M_tot)), axis=1)
            ic = perm[:, 0::2]
            jc = perm[:, 1::2]
            signs = np.where(rng.random((K_rem, M_tot // 2)) < 0.5, 1.0, -1.0)
            apply_bounded_matching_step_2d(cur_s_rem, ic, jc, signs, sigma_h, rng)
            s[idx_rem] = cur_s_rem

            v_fix = np.max(cur_s_rem, axis=1) >= 1.0 - 1e-12
            if np.any(v_fix):
                v_rows = np.where(v_fix)[0]
                winners[idx_rem[v_rows]] = np.argmax(cur_s_rem[v_rows], axis=1)
                active[idx_rem[v_rows]] = False

        n_unresolved = int(np.sum(active))
        n_port1 = int(np.sum((winners >= 0) & (winners < M1)))
        bt = binomtest(n_port1, n_trials - n_unresolved, 0.5)
        ci = bt.proportion_ci(confidence_level=0.95)
        refinement_rows.append(
            {
                "dt": dt,
                "sigma_h": round(sigma_h, 5),
                "lambda_h": round(lam_h, 5),
                "n_trials": n_trials,
                "n_unresolved": n_unresolved,
                "observed_port1_fraction": round(n_port1 / (n_trials - n_unresolved), 4),
                "ci_95": [round(ci.low, 4), round(ci.high, 4)],
                "p_value_vs_050": float(bt.pvalue),
            }
        )

    return {
        "max_uncompensated_backward_residual": round(float(np.max(np.abs(uncomp_residual))), 6),
        "max_compensated_backward_residual": float(np.max(np.abs(comp_residual))),
        "max_sum_of_compensating_drifts": float(np.max(np.abs(np.sum(b_comp, axis=1)))),
        "five_vs_fifteen_step_refinement_s050_theta045": refinement_rows,
    }


def main():
    seed_checks = check_seed_invariances()
    sweep_20atom = run_nine_angle_whole_quantum_sweep(
        M1=10, M2=10, sigma=0.25, use_lattice_tierA=False, n_trials=12000, seed_base=20261110
    )
    sweep_32lattice = run_nine_angle_whole_quantum_sweep(
        M1=16, M2=16, sigma=0.35, use_lattice_tierA=True, n_trials=12000, seed_base=20261150
    )
    eps_vac_table = compute_analytical_eps_vac_table(sigma=0.25)
    threshold_checks = verify_threshold_balance_and_step_refinement()

    results = {
        "periodic_seed_invariance_checks": seed_checks,
        "whole_quantum_theta1_20atom_sweep": sweep_20atom,
        "whole_quantum_theta1_32lattice_sweep": sweep_32lattice,
        "analytical_two_site_eps_vac_drift_table": eps_vac_table,
        "threshold_balance_and_5vs15_refinement": threshold_checks,
    }

    out_path = CODE_DIR / "rev29_verification_results.json"
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)

    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
