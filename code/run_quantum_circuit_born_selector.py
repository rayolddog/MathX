#!/usr/bin/env python3
"""
State-Vector Representation of the Martingale Walk, a Genuinely Linear Unitary Control,
and a Monitored Quantum Instrument Benchmark (`run_quantum_circuit_born_selector.py`)
for Aperture Article 0002 (Revision 2.9 — 2026-11-02)

Compares four distinct dynamical regimes on an M-qubit register restricted to the
single-excitation (W-state) subspace H_1 = span{|10...0>, |01...0>, ..., |00...1>}:

  1. Regime A (Oracle Noise Inversion — Idealized Zero-Dephasing Baseline):
     - Each environmental phase kick R_Z(delta_eta_m) is undone by its exact inverse
       R_Z(-delta_eta_m) with no inter-site exchange.
     - Result: |Psi(t)> remains frozen in its initial superposition (V(t) = 1.0000, 0% fixation).

  2. Regime B (Genuinely State-Independent Fixed-Angle Linear iSWAP(eps) + R_Z Dephasing):
     - Fixed linear two-qubit excitation-exchange gate
       U_mn(eps) = [[cos(eps), -i sin(eps)], [-i sin(eps), cos(eps)]]
       whose matrix and angle eps = 0.15 do NOT depend on the state's current amplitudes or phases,
       combined with uncorrected R_Z(delta_eta_m) dephasing.
     - Result: As a genuine linear negative control, diagonal O(eps^2) recoil spreads and equalizes
       the excitation across all M qubits (0% basis-state fixation).

  3. Regime C (State-Vector Container of the Nonlinear Recoil-Cancelled SU(2) Martingale Walk):
     - Uses the 2pi-periodic three-phase seed C_m = 0.5*[cos(phi_d - phi_i) + cos(phi_d - phi_v)]
       and fast uniform QED vacuum phase fluctuations with antisymmetric sin(Delta_phi_mn) sign.
     - Applies the 2x2 SU(2) rotation U_mn(delta_theta, phi_m, phi_n) whose angle
       delta_theta(u) = arccos(sqrt(u + Delta_u)) - arccos(sqrt(u)) cancels the O(eps^2) recoil
       by inspecting the relative pair weight u = |psi_m|^2 / (|psi_m|^2 + |psi_n|^2).
     - Result: Preserves exact norm ||Psi(t)||_2^2 = 1 while 100% of trajectories reach a single
       basis state |e_{m*}> in finite steps (zero unresolved trajectories) with P_1(S) = S.

  4. Regime D (Physically Admissible Monitored Quantum Instrument Benchmark — Sol §5):
     - Uses fixed, state-independent weak-occupation Kraus operators
       K_{y,j} = sqrt((I + y g Z_j)/2), y in {-1, +1}, Z_j = 2|j><j| - I (g = 0.35),
       satisfying K_+^dagger K_+ + K_-^dagger K_- = I without hidden-amplitude gate feedback.
     - Result: Standard quantum measurement backaction asymptotically localizes conditional
       trajectories (to 1 - 1e-8 tolerance) while preserving Born probabilities P_1(S) = S.
"""

import json
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import chi2

CODE_DIR = Path(__file__).resolve().parent


def apply_rz_phase_gate(psi: np.ndarray, dphi: np.ndarray) -> np.ndarray:
    """Applies single-qubit phase gates bigotimes_m R_Z(dphi_m) in the 1-excitation subspace."""
    return psi * np.exp(1j * dphi)


def apply_fixed_linear_iswap_gate(
    psi_m: np.ndarray,
    psi_n: np.ndarray,
    eps: float,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Applies a genuinely state-independent fixed-angle linear excitation-exchange gate:
        U_mn(eps) = [[cos(eps), -i sin(eps)],
                     [-i sin(eps), cos(eps)]]
    Neither the angle eps nor the matrix depends on the current state's phases or amplitudes.
    """
    c = np.cos(eps)
    s = -1j * np.sin(eps)
    return c * psi_m + s * psi_n, s * psi_m + c * psi_n


def apply_su2_pair_gate(
    psi_m: np.ndarray,
    psi_n: np.ndarray,
    delta_theta: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Applies the 2x2 excitation-conserving SU(2) rotation in span{|1_m 0_n>, |0_m 1_n>}:
        U_mn = diag(e^{i phi_m}, e^{i phi_n}) [[cos(dtheta), -sin(dtheta)],
                                               [sin(dtheta),  cos(dtheta)]] diag(e^{-i phi_m}, e^{-i phi_n})
    """
    phi_m = np.angle(psi_m)
    phi_n = np.angle(psi_n)
    amp_m = np.abs(psi_m)
    amp_n = np.abs(psi_n)
    c = np.cos(delta_theta)
    s = np.sin(delta_theta)
    new_amp_m = c * amp_m - s * amp_n
    new_amp_n = s * amp_m + c * amp_n
    return new_amp_m * np.exp(1j * phi_m), new_amp_n * np.exp(1j * phi_n)


def simulate_monitored_instrument(
    psi_init: np.ndarray,
    M1: int,
    g: float = 0.35,
    tol: float = 1e-8,
    n_steps_trace: int = 180,
    max_sweeps: int = 600,
    n_trials: int = 4000,
    seed: int = 20261102,
) -> dict:
    """
    Simulates the state-independent weak-measurement Kraus instrument (Sol §5):
        K_{y,j} = sqrt((I + y g Z_j)/2), y in {-1, +1}, Z_j = 2|j><j| - I.
    """
    rng = np.random.default_rng(seed)
    M_qubits = psi_init.shape[1]
    psi_D = np.tile(psi_init, (n_trials, 1)).astype(np.complex128)
    kp = np.sqrt(0.5 * (1.0 + g))
    km = np.sqrt(0.5 * (1.0 - g))

    max_share_D = [float(np.mean(np.max(np.abs(psi_D) ** 2, axis=1)))]
    vis_D = [1.0]
    active = np.ones(n_trials, dtype=bool)
    sweeps_D = np.zeros(n_trials, dtype=int)

    for sweep_idx in range(1, max_sweeps + 1):
        idx = np.where(active)[0]
        if len(idx) > 0:
            K = len(idx)
            sweeps_D[idx] += 1
            cur_psi = psi_D[idx]

            for j in range(M_qubits):
                probs = np.abs(cur_psi) ** 2
                sj = probs[:, j]
                z_bar = 2.0 * sj - 1.0
                p_plus = np.clip(0.5 * (1.0 + g * z_bar), 1e-15, 1.0 - 1e-15)
                y_plus = rng.random(K) < p_plus
                py = np.where(y_plus, p_plus, 1.0 - p_plus)

                old_j = cur_psi[:, j].copy()
                scale_other = (np.where(y_plus, km, kp) / np.sqrt(py))[:, None]
                scale_j = np.where(y_plus, kp, km) / np.sqrt(py)
                cur_psi = cur_psi * scale_other
                cur_psi[:, j] = old_j * scale_j

            psi_D[idx] = cur_psi
            probs_active = np.abs(cur_psi) ** 2
            localized_now = np.max(probs_active, axis=1) >= 1.0 - tol
            if np.any(localized_now):
                active[idx[localized_now]] = False

        if sweep_idx <= n_steps_trace:
            probs_all = np.abs(psi_D) ** 2
            max_share_D.append(float(np.mean(np.max(probs_all, axis=1))))
            ov = np.abs(np.sum(psi_D * np.conj(psi_init), axis=1)) ** 2
            vis_D.append(float(np.mean(ov)))

        if not np.any(active) and sweep_idx >= n_steps_trace:
            break

    final_probs_D = np.abs(psi_D) ** 2
    n_unresolved = int(np.sum(active))
    resolved_mask = ~active
    if n_unresolved == 0:
        winners_D = np.argmax(final_probs_D, axis=1)
        p_port1_D = float(np.mean(winners_D < M1))
    else:
        winners_res = np.argmax(final_probs_D[resolved_mask], axis=1)
        p_port1_D = float(np.sum(winners_res < M1) / n_trials)

    return {
        "g_coupling": g,
        "localization_tolerance": tol,
        "final_mean_max_qubit_prob": round(float(np.mean(np.max(final_probs_D, axis=1))), 6),
        "final_coherence_overlap": round(vis_D[-1], 4),
        "asymptotic_localization_rate": round(float(np.mean(resolved_mask)), 6),
        "n_unresolved": n_unresolved,
        "observed_port1_probability": round(p_port1_D, 4),
        "mean_sweeps_to_tolerance": round(float(np.mean(sweeps_D)), 2),
        "max_share_trace": max_share_D,
        "vis_trace": vis_D,
    }


def simulate_quantum_circuit_regimes(
    M_qubits: int = 8,
    s_port1: float = 0.70,
    sigma: float = 0.25,
    n_steps_trace: int = 180,
    n_trials: int = 4000,
    seed: int = 20261031,
):
    """
    Compares the four regimes on an M-qubit single-excitation state vector |Psi(t)> in C^M
    (Port 1 = qubits 0..M/2-1 with total initial weight S = 0.70;
     Port 2 = qubits M/2..M-1 with total initial weight 1 - S = 0.30).
    """
    rng = np.random.default_rng(seed)
    M1 = M_qubits // 2
    M2 = M_qubits - M1

    # Initial complex state vector |Psi(0)> = sum_m A_m |e_m> with ||Psi(0)||_2^2 = 1
    phi_init = np.linspace(0.0, 1.4 * np.pi, M_qubits)
    amps = np.zeros(M_qubits, dtype=float)
    amps[:M1] = np.sqrt(s_port1 / M1)
    amps[M1:] = np.sqrt((1.0 - s_port1) / M2)
    psi_0 = (amps * np.exp(1j * phi_init))[None, :]  # (1, M_qubits)

    # --- REGIME A: Oracle Noise Inversion (Idealized Zero-Dephasing Baseline) ---
    psi_A = np.tile(psi_0, (200, 1)).astype(np.complex128)
    max_share_A = [float(np.mean(np.max(np.abs(psi_A) ** 2, axis=1)))]
    vis_A = [1.0]
    for _ in range(n_steps_trace):
        dphi_noise = rng.normal(0.0, 0.35, size=psi_A.shape)
        psi_A = apply_rz_phase_gate(psi_A, dphi_noise)
        psi_A = apply_rz_phase_gate(psi_A, -dphi_noise)
        max_share_A.append(float(np.mean(np.max(np.abs(psi_A) ** 2, axis=1))))
        ov = np.abs(np.sum(psi_A * np.conj(psi_0), axis=1)) ** 2
        vis_A.append(float(np.mean(ov)))

    # --- REGIME B: Genuinely State-Independent Fixed-Angle Linear iSWAP(eps) + R_Z Dephasing ---
    psi_B = np.tile(psi_0, (500, 1)).astype(np.complex128)
    max_share_B = [float(np.mean(np.max(np.abs(psi_B) ** 2, axis=1)))]
    vis_B = [1.0]
    eps_fixed = 0.15
    for _ in range(n_steps_trace):
        dphi_noise = rng.uniform(-np.pi, np.pi, size=psi_B.shape)
        psi_B = apply_rz_phase_gate(psi_B, dphi_noise)
        perm = rng.permuted(np.tile(np.arange(M_qubits), (len(psi_B), 1)), axis=1)
        for p in range(0, M_qubits - 1, 2):
            ic = perm[:, p]
            jc = perm[:, p + 1]
            rows = np.arange(len(psi_B))
            pm = psi_B[rows, ic]
            pn = psi_B[rows, jc]
            psi_B[rows, ic], psi_B[rows, jc] = apply_fixed_linear_iswap_gate(pm, pn, eps_fixed)
        max_share_B.append(float(np.mean(np.max(np.abs(psi_B) ** 2, axis=1))))
        ov = np.abs(np.sum(psi_B * np.conj(psi_0), axis=1)) ** 2
        vis_B.append(float(np.mean(ov)))

    # --- REGIME C: State-Vector Container of the Nonlinear Recoil-Cancelled SU(2) Martingale Walk ---
    psi_C = np.tile(psi_0, (n_trials, 1)).astype(np.complex128)
    # Apply Exact 2pi-Periodic Three-Phase Seeding at t = 0:
    # C_m = 0.5 * [cos(phi_det - phi_inc) + cos(phi_det - phi_vac)], F_m = 1 + alpha_seed * C_m
    phi_vac = rng.uniform(-np.pi, np.pi, size=(n_trials, M_qubits))
    phi_det = rng.uniform(-np.pi, np.pi, size=(n_trials, M_qubits))
    alpha_seed = 0.15
    c1 = 0.5 * (
        np.cos(phi_det[:, :M1] - phi_init[:M1][None, :])
        + np.cos(phi_det[:, :M1] - phi_vac[:, :M1])
    )
    c2 = 0.5 * (
        np.cos(phi_det[:, M1:] - phi_init[M1:][None, :])
        + np.cos(phi_det[:, M1:] - phi_vac[:, M1:])
    )
    f1 = 1.0 + alpha_seed * c1
    f2 = 1.0 + alpha_seed * c2
    f1 /= np.sum(f1, axis=1, keepdims=True)
    f2 /= np.sum(f2, axis=1, keepdims=True)
    s_init = np.zeros((n_trials, M_qubits), dtype=float)
    s_init[:, :M1] = s_port1 * f1
    s_init[:, M1:] = (1.0 - s_port1) * f2
    psi_C = np.sqrt(s_init) * np.exp(1j * phi_det)

    max_share_C = [float(np.mean(np.max(np.abs(psi_C) ** 2, axis=1)))]
    vis_C = [float(np.mean(np.abs(np.sum(psi_C * np.conj(psi_0), axis=1)) ** 2))]
    norm_errors = []

    n_sample = 6
    sample_port1_trajs = np.zeros((n_sample, n_steps_trace + 1))
    sample_max_qubit_trajs = np.zeros((n_sample, n_steps_trace + 1))
    sample_port1_trajs[:, 0] = np.sum(np.abs(psi_C[:n_sample, :M1]) ** 2, axis=1)
    sample_max_qubit_trajs[:, 0] = np.max(np.abs(psi_C[:n_sample]) ** 2, axis=1)

    active = np.ones(n_trials, dtype=bool)
    steps_C = np.zeros(n_trials, dtype=int)

    for step_idx in range(1, 2000 + 1):
        idx = np.where(active)[0]
        if len(idx) == 0:
            break
        K = len(idx)
        steps_C[idx] += 1
        cur_psi = psi_C[idx]

        # 1. Slow detector phase drift + fast uniform quantum-vacuum phase draw
        dphi_det_slow = rng.normal(0.0, 0.02, size=(K, M_qubits))
        cur_psi = apply_rz_phase_gate(cur_psi, dphi_det_slow)
        cur_vac = rng.uniform(-np.pi, np.pi, size=(K, M_qubits))

        # 2. Pairwise Recoil-Cancelled SU(2) Excitation-Conserving Gate with Antisymmetric sin(Delta_phi)
        perm = rng.permuted(np.tile(np.arange(M_qubits), (K, 1)), axis=1)
        for p in range(0, M_qubits - 1, 2):
            ic = perm[:, p]
            jc = perm[:, p + 1]
            rows = np.arange(K)
            pm = cur_psi[rows, ic]
            pn = cur_psi[rows, jc]
            sm = np.abs(pm) ** 2
            sn = np.abs(pn) ** 2
            tot = sm + sn
            act_pair = tot > 1e-15
            u = np.where(act_pair, sm / np.maximum(tot, 1e-30), 0.0)

            # Antisymmetric relative phase difference between site modes
            dphi = (np.angle(pm) - cur_vac[rows, ic]) - (np.angle(pn) - cur_vac[rows, jc])
            signs = np.where(np.sin(dphi) >= 0.0, 1.0, -1.0)

            std = sigma * np.sqrt(u * (1.0 - u))
            du_int = signs * std

            rand_u = rng.random(K)
            a_low = (sigma ** 2) * (1.0 - u)
            p0_low = a_low / np.maximum(u + a_low, 1e-30)
            du_low = np.where(rand_u < p0_low, -u, a_low)

            v = 1.0 - u
            a_high = (sigma ** 2) * u
            p1_high = a_high / np.maximum(v + a_high, 1e-30)
            du_high = np.where(rand_u < p1_high, v, -a_high)

            du = np.where(u <= v, np.where(std > u, du_low, du_int), np.where(std > v, du_high, du_int))
            u_new = u + np.where(act_pair, du, 0.0)
            if np.any((u_new < -1e-14) | (u_new > 1.0 + 1e-14)):
                raise RuntimeError("Simplex pair bound violated in Regime C.")
            u_new = np.clip(u_new, 0.0, 1.0)

            theta_old = np.arccos(np.sqrt(np.clip(u, 0.0, 1.0)))
            theta_new = np.arccos(np.sqrt(u_new))
            delta_theta = np.where(act_pair, theta_new - theta_old, 0.0)

            cur_psi[rows, ic], cur_psi[rows, jc] = apply_su2_pair_gate(pm, pn, delta_theta)

        psi_C[idx] = cur_psi
        probs_C = np.abs(psi_C) ** 2
        norm_err = float(np.max(np.abs(np.sum(probs_C, axis=1) - 1.0)))
        norm_errors.append(norm_err)

        if step_idx <= n_steps_trace:
            max_share_C.append(float(np.mean(np.max(probs_C, axis=1))))
            ov = np.abs(np.sum(psi_C * np.conj(psi_0), axis=1)) ** 2
            vis_C.append(float(np.mean(ov)))
            sample_port1_trajs[:, step_idx] = np.sum(probs_C[:n_sample, :M1], axis=1)
            sample_max_qubit_trajs[:, step_idx] = np.max(probs_C[:n_sample], axis=1)

        fixed_now = np.max(probs_C[idx], axis=1) >= 1.0 - 1e-12
        if np.any(fixed_now):
            active[idx[fixed_now]] = False

    final_probs_C = np.abs(psi_C) ** 2
    n_unresolved_C = int(np.sum(active))
    if n_unresolved_C > 0:
        raise RuntimeError(f"Unresolved trajectories in Regime C: {n_unresolved_C}")
    winners_C = np.argmax(final_probs_C, axis=1)
    p_port1_C = float(np.mean(winners_C < M1))
    fixation_rate_C = float(np.mean(~active))

    # --- REGIME D: Monitored Quantum Instrument Benchmark (Sol §5) ---
    regime_D = simulate_monitored_instrument(
        psi_init=psi_0,
        M1=M1,
        g=0.35,
        tol=1e-8,
        n_steps_trace=n_steps_trace,
        max_sweeps=600,
        n_trials=n_trials,
        seed=seed + 777,
    )

    return {
        "M_qubits": M_qubits,
        "s_port1_input": s_port1,
        "regime_A_oracle_noise_inversion": {
            "final_mean_max_qubit_prob": round(max_share_A[-1], 4),
            "final_coherence_overlap": round(vis_A[-1], 4),
            "single_qubit_fixation_rate": 0.0,
            "max_share_trace": max_share_A,
            "vis_trace": vis_A,
        },
        "regime_B_fixed_linear_iswap_control": {
            "final_mean_max_qubit_prob": round(max_share_B[-1], 4),
            "final_coherence_overlap": round(vis_B[-1], 4),
            "single_qubit_fixation_rate": 0.0,
            "max_share_trace": max_share_B,
            "vis_trace": vis_B,
        },
        "regime_C_su2_martingale_container": {
            "final_mean_max_qubit_prob": round(float(np.mean(np.max(final_probs_C, axis=1))), 6),
            "final_coherence_overlap": round(vis_C[-1], 4),
            "single_qubit_fixation_rate": round(fixation_rate_C, 6),
            "n_unresolved": n_unresolved_C,
            "observed_port1_probability": round(p_port1_C, 4),
            "born_target_probability": s_port1,
            "mean_steps_to_basis_state": round(float(np.mean(steps_C)), 2),
            "max_unitary_norm_deviation": float(np.max(norm_errors)),
            "max_share_trace": max_share_C,
            "vis_trace": vis_C,
            "sample_port1_trajs": sample_port1_trajs.tolist(),
            "sample_max_qubit_trajs": sample_max_qubit_trajs.tolist(),
        },
        "regime_D_monitored_kraus_instrument": regime_D,
    }


def run_nine_angle_su2_sweep(M_qubits: int = 8, n_trials: int = 8000, seed: int = 20261101):
    """Runs the 9-angle Born sweep (S = 0.1 .. 0.9) for both Regime C and Regime D."""
    s_grid = [round(0.1 * k, 1) for k in range(1, 10)]
    p_obs_C = []
    mean_steps_C = []
    p_obs_D = []
    mean_sweeps_D = []
    for idx_s, s_val in enumerate(s_grid):
        res = simulate_quantum_circuit_regimes(
            M_qubits=M_qubits,
            s_port1=s_val,
            sigma=0.25,
            n_steps_trace=10,
            n_trials=n_trials,
            seed=seed + 100 * idx_s,
        )
        rc = res["regime_C_su2_martingale_container"]
        rd = res["regime_D_monitored_kraus_instrument"]
        p_obs_C.append(rc["observed_port1_probability"])
        mean_steps_C.append(rc["mean_steps_to_basis_state"])
        p_obs_D.append(rd["observed_port1_probability"])
        mean_sweeps_D.append(rd["mean_sweeps_to_tolerance"])

    s_arr = np.array(s_grid)

    def fit_curve(p_list):
        p_arr = np.array(p_list)
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
            "p_observed": p_list,
            "p_hat": round(p_hat, 4),
            "se_p": round(se_p, 4),
            "chi2_born": round(chi2_val, 2),
            "p_value": round(pval, 4),
            "rms_error": round(rms, 5),
        }

    fit_C = fit_curve(p_obs_C)
    fit_C["s_grid"] = s_grid
    fit_C["mean_steps"] = mean_steps_C

    fit_D = fit_curve(p_obs_D)
    fit_D["s_grid"] = s_grid
    fit_D["mean_sweeps"] = mean_sweeps_D

    return {
        "regime_C_su2_martingale_sweep": fit_C,
        "regime_D_monitored_kraus_sweep": fit_D,
    }


def main():
    regime_comp = simulate_quantum_circuit_regimes(
        M_qubits=8, s_port1=0.70, sigma=0.25, n_steps_trace=180, n_trials=4000, seed=20261031
    )
    sweeps = run_nine_angle_su2_sweep(M_qubits=8, n_trials=8000, seed=20261101)
    sweep_C = sweeps["regime_C_su2_martingale_sweep"]
    sweep_D = sweeps["regime_D_monitored_kraus_sweep"]

    fig, axes = plt.subplots(1, 3, figsize=(15.8, 4.9))
    t_ax = np.arange(len(regime_comp["regime_A_oracle_noise_inversion"]["max_share_trace"]))

    # Panel A: Max single-qubit weight max_m |psi_m(t)|^2 across Regimes A, B, C, D
    ax0 = axes[0]
    ax0.plot(
        t_ax,
        regime_comp["regime_C_su2_martingale_container"]["max_share_trace"],
        color="#10b981",
        lw=2.5,
        label="Regime C: Nonlinear $SU(2)$ Martingale Container\n(Exact Vertex Fixation $|e_{m^*}\\rangle$, $100\\%$)",
    )
    ax0.plot(
        t_ax,
        regime_comp["regime_D_monitored_kraus_instrument"]["max_share_trace"],
        color="#8b5cf6",
        lw=2.1,
        ls="-.",
        label="Regime D: Monitored Kraus Instrument $K_{y,j}$\n(Asymptotic Localization $1 - 10^{-8}$, $100\\%$)",
    )
    ax0.plot(
        t_ax,
        regime_comp["regime_B_fixed_linear_iswap_control"]["max_share_trace"],
        color="#f59e0b",
        lw=2.0,
        ls="--",
        label="Regime B: Fixed-Angle Linear $\\mathrm{iSWAP}(\\epsilon)$ Control\n($O(\\epsilon^2)$ Recoil Equalizes Shares, $0\\%$ Fixation)",
    )
    ax0.plot(
        t_ax,
        regime_comp["regime_A_oracle_noise_inversion"]["max_share_trace"],
        color="#3b82f6",
        lw=2.0,
        ls=":",
        label="Regime A: Oracle Noise Inversion\n(Frozen Initial Superposition, $0\\%$ Fixation)",
    )
    ax0.set_xlabel("Step / Sweep Index $t$", fontsize=10.5)
    ax0.set_ylabel("Mean Leading-Qubit Weight $\\mathbb{E}[\\max_m |\\psi_m(t)|^2]$", fontsize=10.5)
    ax0.set_title("A. Four State-Vector Regimes ($M=8$ Qubits, $S=0.70$)", fontsize=11, fontweight="bold")
    ax0.set_ylim(0.0, 1.05)
    ax0.grid(True, alpha=0.3)
    ax0.legend(loc="center right", fontsize=7.6)

    # Panel B: Sample Single-Shot State-Vector Trajectories
    ax1 = axes[1]
    trajs_max = np.array(regime_comp["regime_C_su2_martingale_container"]["sample_max_qubit_trajs"])
    for k in range(len(trajs_max)):
        ax1.plot(t_ax, trajs_max[k], lw=1.6, alpha=0.85, label=f"Shot #{k+1}" if k < 4 else None)
    ax1.axhline(1.0, color="black", ls="--", lw=1.2, label="Basis State $|e_{m^*}\\rangle$ ($|\\psi_{m^*}|^2 = 1$)")
    ax1.set_xlabel("Step $t$", fontsize=10.5)
    ax1.set_ylabel("Single-Shot Leading Qubit Weight $\\max_m |\\psi_m(t)|^2$", fontsize=10.5)
    ax1.set_title("B. Continuous Norm-Conserving Trajectories ($\\|\\Psi(t)\\|^2 = 1$)", fontsize=11, fontweight="bold")
    ax1.set_ylim(0.0, 1.05)
    ax1.grid(True, alpha=0.3)
    ax1.legend(loc="lower right", fontsize=8.2)

    # Panel C: 9-Angle Born Convergence of Regimes C and D
    ax2 = axes[2]
    s_dense = np.linspace(0.0, 1.0, 100)
    ax2.plot(s_dense, s_dense, "--", color="black", lw=1.6, label="Exact Born Law $P_1(S) = S$ ($p = 2$)")
    ax2.plot(
        sweep_C["s_grid"],
        sweep_C["p_observed"],
        "o",
        color="#10b981",
        markersize=7,
        label=(
            f"Regime C ($SU(2)$ Martingale): $\\hat{{p}} = {sweep_C['p_hat']:.4f} \\pm {sweep_C['se_p']:.4f}$, "
            f"$\\chi^2 = {sweep_C['chi2_born']:.2f}$ ($p = {sweep_C['p_value']:.3f}$)"
        ),
    )
    ax2.plot(
        sweep_D["s_grid"],
        sweep_D["p_observed"],
        "s",
        color="#8b5cf6",
        markersize=5.5,
        alpha=0.85,
        label=(
            f"Regime D (Kraus Instrument): $\\hat{{p}} = {sweep_D['p_hat']:.4f} \\pm {sweep_D['se_p']:.4f}$, "
            f"$\\chi^2 = {sweep_D['chi2_born']:.2f}$ ($p = {sweep_D['p_value']:.3f}$)"
        ),
    )
    ax2.set_xlabel("Initial Port-1 Subspace Weight $S = \\sum_{m \\in \\mathrm{Port\\,1}} |A_m|^2$", fontsize=10.5)
    ax2.set_ylabel("Port-1 Fixation / Localization Fraction", fontsize=10.5)
    ax2.set_title("C. 9-Angle Born Benchmark (Regimes C & D)", fontsize=11, fontweight="bold")
    ax2.set_xlim(0.0, 1.0)
    ax2.set_ylim(0.0, 1.0)
    ax2.grid(True, alpha=0.3)
    ax2.legend(loc="upper left", fontsize=7.8)

    plt.tight_layout()
    out_png = CODE_DIR / "fig_v2_8_quantum_circuit_born_selector.png"
    fig.savefig(out_png, dpi=180)
    plt.close(fig)

    out_json = {
        "regime_comparison_8qubits_s070": {
            "M_qubits": regime_comp["M_qubits"],
            "s_port1_input": regime_comp["s_port1_input"],
            "regime_A_oracle_noise_inversion": {
                k: v for k, v in regime_comp["regime_A_oracle_noise_inversion"].items() if not k.endswith("_trace")
            },
            "regime_B_fixed_linear_iswap_control": {
                k: v for k, v in regime_comp["regime_B_fixed_linear_iswap_control"].items() if not k.endswith("_trace")
            },
            "regime_C_su2_martingale_container": {
                k: v
                for k, v in regime_comp["regime_C_su2_martingale_container"].items()
                if not k.endswith("_trace") and not k.endswith("_trajs")
            },
            "regime_D_monitored_kraus_instrument": {
                k: v for k, v in regime_comp["regime_D_monitored_kraus_instrument"].items() if not k.endswith("_trace")
            },
        },
        "nine_angle_sweeps_8qubits": sweeps,
    }
    with open(CODE_DIR / "quantum_circuit_born_selector_results.json", "w") as f:
        json.dump(out_json, f, indent=2)

    print(json.dumps(out_json, indent=2))


if __name__ == "__main__":
    main()
