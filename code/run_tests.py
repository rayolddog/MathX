import numpy as np
import matplotlib.pyplot as plt
import json

rng = np.random.default_rng(20261014)

def simulate_two_site_stalling_min(s0=0.8, theta=0.45, step=0.25, commit_rate_per_step=0.0, n_trials=20000, max_steps=5000):
    """
    Replicates the paper's v0.9 discrete engine:
    delta = +- step * min(s1, s2)
    with optional threshold-gated linear hazard: rate_i = commit_rate_per_step * s_i if s_i > theta else 0.
    """
    s1 = np.full(n_trials, s0, dtype=float)
    active = np.ones(n_trials, dtype=bool)
    winner = np.full(n_trials, -1, dtype=int) # 1 if site 1 wins, 0 if site 2 wins
    total_steps = np.zeros(n_trials, dtype=int)
    exposed_steps_1 = np.zeros(n_trials, dtype=int)
    exposed_steps_2 = np.zeros(n_trials, dtype=int)

    absorb_tol = 0.995 # absorbing boundary as in Appendix D

    for t in range(1, max_steps + 1):
        idx = np.where(active)[0]
        if len(idx) == 0:
            break
        cur_s1 = s1[idx]
        cur_s2 = 1.0 - cur_s1

        # Track exposure above theta
        exp1 = cur_s1 > theta
        exp2 = cur_s2 > theta
        exposed_steps_1[idx] += exp1.astype(int)
        exposed_steps_2[idx] += exp2.astype(int)
        total_steps[idx] += 1

        # Gated commit check if commit_rate_per_step > 0
        if commit_rate_per_step > 0.0:
            p1 = np.where(exp1, commit_rate_per_step * cur_s1, 0.0)
            p2 = np.where(exp2, commit_rate_per_step * cur_s2, 0.0)
            p_tot = p1 + p2
            u_commit = rng.random(len(idx))
            committed = u_commit < p_tot
            if np.any(committed):
                c_idx = idx[committed]
                # Conditional pick between exposed sites proportional to p1, p2
                u_pick = rng.random(len(c_idx)) * p_tot[committed]
                w1 = u_pick < p1[committed]
                winner[c_idx] = np.where(w1, 1, 0)
                active[c_idx] = False
                # Filter remaining active
                keep = ~committed
                idx = idx[keep]
                cur_s1 = cur_s1[keep]
                cur_s2 = cur_s2[keep]
                if len(idx) == 0:
                    break

        # Discrete step: delta = +- step * min(s1, s2)
        signs = rng.choice([-1.0, 1.0], size=len(idx))
        delta = signs * step * np.minimum(cur_s1, cur_s2)
        new_s1 = cur_s1 + delta

        # Check absorption at 1 - tol or tol
        win1 = new_s1 >= absorb_tol
        win2 = new_s1 <= (1.0 - absorb_tol)
        done = win1 | win2
        if np.any(done):
            d_idx = idx[done]
            winner[d_idx] = np.where(win1[done], 1, 0)
            active[d_idx] = False

        s1[idx] = new_s1

    # For any still active at max_steps, assign by argmax
    if np.any(active):
        a_idx = np.where(active)[0]
        winner[a_idx] = (s1[a_idx] >= 0.5).astype(int)

    # Calculate winner's exposed fraction
    win1_mask = (winner == 1)
    win2_mask = (winner == 0)
    winner_exp_steps = np.where(win1_mask, exposed_steps_1, exposed_steps_2)
    exp_fraction = np.mean(winner_exp_steps / np.maximum(total_steps, 1))

    return {
        "p_win1": float(np.mean(winner == 1)),
        "mean_total_steps": float(np.mean(total_steps)),
        "mean_winner_exp_steps": float(np.mean(winner_exp_steps)),
        "mean_exp_fraction": float(exp_fraction),
    }


def simulate_two_site_conserving_sqrt(s0=0.8, theta=0.45, sigma=0.25, commit_rate_per_step=0.0, nonrandom_vac_eps=0.0, n_trials=20000, max_steps=5000):
    """
    True conserving Wright-Fisher / cross-term step on the simplex s1 + s2 = 1:
    ds1 = nonrandom_vac_eps * sqrt(s1*s2)*(sqrt(s2) - sqrt(s1)) + sigma * sqrt(s1 * s2) * dW
    Notice that at s1=0.5, sqrt(s1*s2) = 0.5, matching min(s1,s2)=0.5 at the midpoint so initial speed is identical!
    Also implements sticky/absorbing boundaries at 0 and 1.
    """
    s1 = np.full(n_trials, s0, dtype=float)
    active = np.ones(n_trials, dtype=bool)
    winner = np.full(n_trials, -1, dtype=int)
    total_steps = np.zeros(n_trials, dtype=int)
    exposed_steps_1 = np.zeros(n_trials, dtype=int)
    exposed_steps_2 = np.zeros(n_trials, dtype=int)

    absorb_tol = 0.995

    for t in range(1, max_steps + 1):
        idx = np.where(active)[0]
        if len(idx) == 0:
            break
        cur_s1 = s1[idx]
        cur_s2 = 1.0 - cur_s1

        exp1 = cur_s1 > theta
        exp2 = cur_s2 > theta
        exposed_steps_1[idx] += exp1.astype(int)
        exposed_steps_2[idx] += exp2.astype(int)
        total_steps[idx] += 1

        if commit_rate_per_step > 0.0:
            p1 = np.where(exp1, commit_rate_per_step * cur_s1, 0.0)
            p2 = np.where(exp2, commit_rate_per_step * cur_s2, 0.0)
            p_tot = p1 + p2
            u_commit = rng.random(len(idx))
            committed = u_commit < p_tot
            if np.any(committed):
                c_idx = idx[committed]
                u_pick = rng.random(len(c_idx)) * p_tot[committed]
                w1 = u_pick < p1[committed]
                winner[c_idx] = np.where(w1, 1, 0)
                active[c_idx] = False
                keep = ~committed
                idx = idx[keep]
                cur_s1 = cur_s1[keep]
                cur_s2 = cur_s2[keep]
                if len(idx) == 0:
                    break

        # True conserving sqrt(s1 * s2) step + optional non-random vacuum cross-term drift
        dW = rng.choice([-1.0, 1.0], size=len(idx)) # Rademacher or Gaussian increment with unit variance
        drift = nonrandom_vac_eps * np.sqrt(cur_s1 * cur_s2) * (np.sqrt(cur_s2) - np.sqrt(cur_s1))
        delta = drift + sigma * np.sqrt(cur_s1 * cur_s2) * dW
        new_s1 = np.clip(cur_s1 + delta, 0.0, 1.0)

        win1 = new_s1 >= absorb_tol
        win2 = new_s1 <= (1.0 - absorb_tol)
        done = win1 | win2
        if np.any(done):
            d_idx = idx[done]
            winner[d_idx] = np.where(win1[done], 1, 0)
            active[d_idx] = False

        s1[idx] = new_s1

    if np.any(active):
        a_idx = np.where(active)[0]
        winner[a_idx] = (s1[a_idx] >= 0.5).astype(int)

    win1_mask = (winner == 1)
    winner_exp_steps = np.where(win1_mask, exposed_steps_1, exposed_steps_2)
    exp_fraction = np.mean(winner_exp_steps / np.maximum(total_steps, 1))

    return {
        "p_win1": float(np.mean(winner == 1)),
        "mean_total_steps": float(np.mean(total_steps)),
        "mean_winner_exp_steps": float(np.mean(winner_exp_steps)),
        "mean_exp_fraction": float(exp_fraction),
    }


def simulate_ten_site_comparison(theta=0.45, step_scale=0.15, commit_rate_per_step=0.0, nonrandom_vac_eps=0.0, mode="sqrt", n_trials=10000, max_steps=4000):
    """
    10-site test configuration from the paper:
    1 bright site of amplitude 3 (weight 9/18 = 0.500),
    9 dim sites of amplitude 1 (weight 1/18 each).
    Compares:
      - mode="min": pairwise stakes-scaled exchange +- step * min(s_i, s_j)
      - mode="sqrt": conserving Wright-Fisher exchange ds_i = sum_j sigma * (delta_ij sqrt(s_i) - s_i sqrt(s_j)) dW_j
                     + nonrandom_vac_eps * (sqrt(s_i) - s_i * sum_j sqrt(s_j))
    """
    N = 10
    s_init = np.array([9.0/18.0] + [1.0/18.0]*9)
    S = np.tile(s_init, (n_trials, 1))
    active = np.ones(n_trials, dtype=bool)
    winner = np.full(n_trials, -1, dtype=int)
    total_steps = np.zeros(n_trials, dtype=int)
    exposed_steps = np.zeros((n_trials, N), dtype=int)
    absorb_tol = 0.99

    for t in range(1, max_steps + 1):
        idx = np.where(active)[0]
        if len(idx) == 0:
            break
        cur_S = S[idx] # (M, N)
        M = len(idx)

        exp_mask = cur_S > theta
        exposed_steps[idx] += exp_mask.astype(int)
        total_steps[idx] += 1

        if commit_rate_per_step > 0.0:
            rates = np.where(exp_mask, commit_rate_per_step * cur_S, 0.0)
            r_tot = np.sum(rates, axis=1)
            u_commit = rng.random(M)
            committed = u_commit < r_tot
            if np.any(committed):
                c_idx = idx[committed]
                c_rates = rates[committed]
                c_rtot = r_tot[committed]
                # categorical draw
                cum_probs = np.cumsum(c_rates / c_rtot[:, None], axis=1)
                u_pick = rng.random(len(c_idx))[:, None]
                chosen = np.argmax(u_pick < cum_probs, axis=1)
                winner[c_idx] = chosen
                active[c_idx] = False
                keep = ~committed
                idx = idx[keep]
                cur_S = cur_S[keep]
                M = len(idx)
                if M == 0:
                    break

        if mode == "min":
            # Pick random pair of alive sites per trial or all-pairs exchange
            # In multi-site discrete gambler, do N/2 random disjoint pair exchanges per step
            perm = rng.permuted(np.tile(np.arange(N), (M, 1)), axis=1)
            for p in range(0, N, 2):
                i_col = perm[:, p]
                j_col = perm[:, p+1]
                si = cur_S[np.arange(M), i_col]
                sj = cur_S[np.arange(M), j_col]
                signs = rng.choice([-1.0, 1.0], size=M)
                # optional nonrandom vacuum drift toward smaller site in pair
                drift = nonrandom_vac_eps * np.sqrt(si * sj) * (np.sqrt(sj) - np.sqrt(si))
                delta = drift + signs * (step_scale * 2.5) * np.minimum(si, sj)
                delta = np.clip(delta, -si, sj)
                cur_S[np.arange(M), i_col] += delta
                cur_S[np.arange(M), j_col] -= delta
        elif mode == "sqrt":
            # Conserving Wright-Fisher / projected amplitude-linear SDE on simplex
            dW = rng.normal(0.0, 1.0, size=(M, N))
            sqrt_S = np.sqrt(np.maximum(cur_S, 0.0))
            sum_sqrt_dW = np.sum(sqrt_S * dW, axis=1, keepdims=True)
            dS_noise = step_scale * (sqrt_S * dW - cur_S * sum_sqrt_dW)
            # Nonrandom vacuum interference cross-term: +eps * sqrt(s_i) projected onto simplex
            sum_sqrt = np.sum(sqrt_S, axis=1, keepdims=True)
            dS_vac = nonrandom_vac_eps * (sqrt_S - cur_S * sum_sqrt)
            cur_S = cur_S + dS_noise + dS_vac
            cur_S = np.maximum(cur_S, 0.0)
            row_sums = np.sum(cur_S, axis=1, keepdims=True)
            cur_S = np.where(row_sums > 0, cur_S / row_sums, s_init)

        max_s = np.max(cur_S, axis=1)
        done = max_s >= absorb_tol
        if np.any(done):
            d_idx = idx[done]
            winner[d_idx] = np.argmax(cur_S[done], axis=1)
            active[d_idx] = False

        S[idx] = cur_S

    if np.any(active):
        a_idx = np.where(active)[0]
        winner[a_idx] = np.argmax(S[a_idx], axis=1)

    win_exp = exposed_steps[np.arange(n_trials), winner]
    exp_fraction = np.mean(win_exp / np.maximum(total_steps, 1))

    return {
        "p_bright": float(np.mean(winner == 0)),
        "mean_total_steps": float(np.mean(total_steps)),
        "mean_winner_exp_steps": float(np.mean(win_exp)),
        "mean_exp_fraction": float(exp_fraction),
    }

results = {}

# 1. Pure absorbing game (commit_rate = 0): compare exposure time above theta
results["two_site_exposure"] = {}
for th in [0.45, 0.72, 0.90, 0.94]:
    res_min = simulate_two_site_stalling_min(s0=0.5, theta=th, step=0.25, commit_rate_per_step=0.0)
    res_sqrt = simulate_two_site_conserving_sqrt(s0=0.5, theta=th, sigma=0.25, commit_rate_per_step=0.0)
    results["two_site_exposure"][str(th)] = {"min_step": res_min, "sqrt_step": res_sqrt}

results["ten_site_exposure"] = {}
for th in [0.45, 0.72, 0.90, 0.94]:
    res_min10 = simulate_ten_site_comparison(theta=th, step_scale=0.15, commit_rate_per_step=0.0, mode="min")
    res_sqrt10 = simulate_ten_site_comparison(theta=th, step_scale=0.15, commit_rate_per_step=0.0, mode="sqrt")
    results["ten_site_exposure"][str(th)] = {"min_step": res_min10, "sqrt_step": res_sqrt10}

# 2. Gated commitment bias on 80/20 split (s0 = 0.80, theta = 0.45 and theta = 0.72) at matched physical commit rates lambda
results["gated_bias_80_20"] = {}
for th in [0.45, 0.72]:
    results["gated_bias_80_20"][str(th)] = {}
    for lam in [0.0, 0.01, 0.03, 0.05, 0.10, 0.20]:
        r_min = simulate_two_site_stalling_min(s0=0.8, theta=th, step=0.25, commit_rate_per_step=lam)
        r_sqrt = simulate_two_site_conserving_sqrt(s0=0.8, theta=th, sigma=0.25, commit_rate_per_step=lam)
        results["gated_bias_80_20"][str(th)][str(lam)] = {
            "min_step_p1": r_min["p_win1"],
            "min_step_bias": r_min["p_win1"] - 0.80,
            "sqrt_step_p1": r_sqrt["p_win1"],
            "sqrt_step_bias": r_sqrt["p_win1"] - 0.80,
        }

# 3. Non-random vacuum interference cross-term (+eps * sqrt(e_i)) boosting lower-amplitude sites:
# Test how nonrandom_vac_eps affects:
# (a) The ungated game (lam = 0)
# (b) The threshold-gated game (theta = 0.72, lam = 0.05) where the gate otherwise biases toward the bright site!
results["vacuum_compensation_80_20"] = {}
for eps in [0.0, 0.005, 0.010, 0.015, 0.020, 0.030, 0.050]:
    r_ungated = simulate_two_site_conserving_sqrt(s0=0.8, theta=0.72, sigma=0.25, commit_rate_per_step=0.0, nonrandom_vac_eps=eps)
    r_gated = simulate_two_site_conserving_sqrt(s0=0.8, theta=0.72, sigma=0.25, commit_rate_per_step=0.05, nonrandom_vac_eps=eps)
    results["vacuum_compensation_80_20"][str(eps)] = {
        "ungated_p1": r_ungated["p_win1"],
        "gated_p1": r_gated["p_win1"],
        "gated_bias": r_gated["p_win1"] - 0.80
    }

with open("/Users/john-bramble/.gemini/antigravity/scratch/selector_verification/results.json", "w") as f:
    json.dump(results, f, indent=2)

print(json.dumps(results, indent=2))
