import os
os.environ["MPLCONFIGDIR"] = "/tmp/mplconfig"
import numpy as np
import matplotlib.pyplot as plt
import json

rng = np.random.default_rng(20261014)

def run_10site_exact_pairwise(theta=0.45, step=0.25, commit_rate_per_step=0.0, nonrandom_vac_eps=0.0, mode="sqrt", n_trials=15000, max_steps=4000):
    """
    Exact pairwise conserving exchange on 10 sites:
    Site 0 has initial share 9/18 = 0.500; Sites 1..9 have initial share 1/18 each.
    At each step, a random pair of alive sites (s_i > 0, s_j > 0) exchanges energy:
      - mode="min":  delta = +- step * min(s_i, s_j)
      - mode="sqrt": delta = +- step * sqrt(s_i * s_j), capped symmetrically at min(s_i, s_j) so expectation is strictly 0!
                     Wait: if capped at min(s_i, s_j) with symmetric +- sign, E[delta] = 0 identically!
                     Let's use step * sqrt(s_i * s_j) where step <= 1, and if step * sqrt(s_i * s_j) > min(s_i, s_j),
                     to keep E[delta] = 0 (exact martingale) while preventing negative shares, use asymmetric absorbing jump:
                     with prob p_+ = min(s_i, s_j) / (min(s_i, s_j) + step*sqrt(s_i*s_j)) or simply take sub-steps so step*sqrt(s_i*s_j) is a true diffusion!
    Even better: in continuous Wright-Fisher diffusion between a pair (s_i, s_j) with total s_pair = s_i + s_j,
    let u = s_i / s_pair in (0, 1).
    Then du = sigma * sqrt(u(1-u)) dW.
    To make E[du] = 0 EXACLY at every discrete step without ever stepping outside [0, 1]:
    Let delta_u_max = min(u, 1-u). If we use a two-point or three-point zero-mean distribution on [-u, 1-u] with variance sigma^2 u(1-u):
       specifically, if sigma * sqrt(u(1-u)) <= min(u, 1-u), delta_u = +- sigma * sqrt(u(1-u)) with prob 1/2 each!
      if sigma * sqrt(u(1-u)) > min(u, 1-u) (near the boundary where u < sigma^2 / (1 + sigma^2)),
      the exact zero-mean boundary-hitting step jumps to 0 with prob (1-u) and to 1 with prob u (or to 2u / 0 with prob 1/2)!
      Jumping to 0 or 2u (if u < 1/2) has step size u, wait: to have variance sigma^2 u(1-u) AND hit 0 directly:
      Jump to u'=0 (delta = -u) with prob p_0, and jump to u' = u + a (delta = +a) with prob 1 - p_0,
      where p_0 / (1 - p_0) = a / u (zero mean!) and p_0 * u^2 + (1 - p_0) * a^2 = u * a = sigma^2 u(1 - u) ==> a = sigma^2 (1 - u)!
    WAIT A SECOND! Look at that exact identity:
      If u < sigma^2 / (1 + sigma^2), set:
        delta_u = -u                 with probability p_0 = sigma^2(1-u) / (u + sigma^2(1-u))
        delta_u = +sigma^2(1 - u)    with probability 1 - p_0 = u / (u + sigma^2(1-u))
      Check the mean: E[delta_u] = p_0 * (-u) + (1 - p_0) * sigma^2(1-u) = 0 IDENTICALLY!
      Check the variance: E[delta_u^2] = p_0 * u^2 + (1 - p_0) * sigma^4(1-u)^2 = u * sigma^2(1-u) = sigma^2 u(1-u) IDENTICALLY!
      And u + delta_u is ALWAYS in [0, 1], hitting 0 cleanly with ZERO boundary clipping bias!
    """
    pass

def exact_wf_step(s1, s2, sigma, vac_eps, rng):
    """
    Takes an exact conserving, zero-drift (when vac_eps=0) Wright-Fisher step between shares s1 and s2,
    with exact variance sigma^2 * s1 * s2 and exact vacuum drift vac_eps * sqrt(s1*s2)*(sqrt(s2) - sqrt(s1)),
    NEVER stepping outside [0, s1+s2] and NEVER introducing clipping drift!
    """
    tot = s1 + s2
    mask = tot > 1e-15
    u = np.where(mask, s1 / np.maximum(tot, 1e-15), 0.0)
    
    # Non-random vacuum cross-term drift on u = s1 / (s1 + s2):
    # du_vac = vac_eps * sqrt(u(1-u)) * (sqrt(1-u) - sqrt(u))
    du_vac = vac_eps * np.sqrt(u * (1.0 - u)) * (np.sqrt(1.0 - u) - np.sqrt(u))
    u = np.clip(u + du_vac, 0.0, 1.0)

    # Exact zero-mean, variance = sigma^2 * u * (1-u) step on [0, 1]:
    std = sigma * np.sqrt(u * (1.0 - u))
    min_side = np.minimum(u, 1.0 - u)

    # Case 1: interior (std <= min_side): symmetric +- std step
    signs = rng.choice([-1.0, 1.0], size=len(s1))
    du_interior = signs * std

    # Case 2: near 0 (u < std and u <= 0.5):
    # Jump to 0 (du = -u) with prob p0 = std^2 / (u^2 + std^2), or +std^2/u with prob u^2 / (u^2 + std^2)
    # Wait: std^2 = sigma^2 u (1-u), so std^2 / u = sigma^2 (1-u).
    a_low = (sigma**2) * (1.0 - u)
    p0_low = a_low / np.maximum(u + a_low, 1e-15)
    rand_u = rng.random(len(s1))
    du_low = np.where(rand_u < p0_low, -u, a_low)

    # Case 3: near 1 (1-u < std and u > 0.5):
    # By symmetry, v = 1 - u jumps to 0 (du = +(1-u)) with prob p0_high, or -sigma^2 * u
    v = 1.0 - u
    a_high = (sigma**2) * u
    p0_high = a_high / np.maximum(v + a_high, 1e-15)
    du_high = np.where(rand_u < p0_high, v, -a_high)

    near_low = (std > min_side) & (u <= 0.5) & (u > 0.0)
    near_high = (std > min_side) & (u > 0.5) & (u < 1.0)

    du = np.where(near_low, du_low, np.where(near_high, du_high, du_interior))
    du = np.where((u <= 1e-12) | (u >= 1.0 - 1e-12), 0.0, du)

    u_new = np.clip(u + du, 0.0, 1.0)
    s1_new = u_new * tot
    s2_new = (1.0 - u_new) * tot
    return s1_new, s2_new


def simulate_two_site_exact(s0=0.8, theta=0.45, sigma=0.25, commit_rate_per_step=0.0, vac_eps=0.0, mode="sqrt", n_trials=30000, max_steps=5000):
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

        if mode == "min":
            du_vac = vac_eps * np.sqrt(cur_s1 * cur_s2) * (np.sqrt(cur_s2) - np.sqrt(cur_s1))
            signs = rng.choice([-1.0, 1.0], size=len(idx))
            delta = du_vac + signs * sigma * np.minimum(cur_s1, cur_s2)
            new_s1 = np.clip(cur_s1 + delta, 0.0, 1.0)
        else:
            new_s1, _ = exact_wf_step(cur_s1, cur_s2, sigma, vac_eps, rng)

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


def simulate_ten_site_exact(theta=0.45, sigma=0.25, commit_rate_per_step=0.0, vac_eps=0.0, mode="sqrt", n_trials=15000, max_steps=4000):
    N = 10
    s_init = np.array([9.0/18.0] + [1.0/18.0]*9)
    S = np.tile(s_init, (n_trials, 1))
    active = np.ones(n_trials, dtype=bool)
    winner = np.full(n_trials, -1, dtype=int)
    total_steps = np.zeros(n_trials, dtype=int)
    exposed_steps = np.zeros((n_trials, N), dtype=int)
    absorb_tol = 0.995

    for t in range(1, max_steps + 1):
        idx = np.where(active)[0]
        if len(idx) == 0:
            break
        cur_S = S[idx]
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

        # Pairwise random matchings among the 10 sites
        perm = rng.permuted(np.tile(np.arange(N), (M, 1)), axis=1)
        for p in range(0, N, 2):
            i_col = perm[:, p]
            j_col = perm[:, p+1]
            si = cur_S[np.arange(M), i_col]
            sj = cur_S[np.arange(M), j_col]
            if mode == "min":
                signs = rng.choice([-1.0, 1.0], size=M)
                du_vac = vac_eps * np.sqrt(si * sj) * (np.sqrt(sj) - np.sqrt(si))
                delta = du_vac + signs * sigma * np.minimum(si, sj)
                delta = np.clip(delta, -si, sj)
                cur_S[np.arange(M), i_col] = si + delta
                cur_S[np.arange(M), j_col] = sj - delta
            else:
                si_new, sj_new = exact_wf_step(si, sj, sigma, vac_eps, rng)
                cur_S[np.arange(M), i_col] = si_new
                cur_S[np.arange(M), j_col] = sj_new

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

out = {}
out["two_site_exposure"] = {}
for th in [0.45, 0.72, 0.90, 0.94]:
    out["two_site_exposure"][str(th)] = {
        "min": simulate_two_site_exact(s0=0.5, theta=th, sigma=0.25, mode="min"),
        "sqrt": simulate_two_site_exact(s0=0.5, theta=th, sigma=0.25, mode="sqrt")
    }

out["ten_site_exposure"] = {}
for th in [0.45, 0.72, 0.90, 0.94]:
    out["ten_site_exposure"][str(th)] = {
        "min": simulate_ten_site_exact(theta=th, sigma=0.25, mode="min"),
        "sqrt": simulate_ten_site_exact(theta=th, sigma=0.25, mode="sqrt")
    }

# Sweep S in [0.1..0.9] for 2-site gated commitment at theta=0.72, lambda=0.04:
# Compare:
# 1. Stalling min() step (gated)
# 2. Conserving sqrt() step (gated, vac_eps=0)
# 3. Conserving sqrt() step + non-random vacuum cross-term (gated, vac_eps=0.018)
s_grid = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
curve_min = []
curve_sqrt = []
curve_sqrt_vac = []
for s0 in s_grid:
    r_m = simulate_two_site_exact(s0=s0, theta=0.72, sigma=0.25, commit_rate_per_step=0.04, vac_eps=0.0, mode="min", n_trials=20000)
    r_s = simulate_two_site_exact(s0=s0, theta=0.72, sigma=0.25, commit_rate_per_step=0.04, vac_eps=0.0, mode="sqrt", n_trials=20000)
    r_v = simulate_two_site_exact(s0=s0, theta=0.72, sigma=0.25, commit_rate_per_step=0.04, vac_eps=0.018, mode="sqrt", n_trials=20000)
    curve_min.append(r_m["p_win1"])
    curve_sqrt.append(r_s["p_win1"])
    curve_sqrt_vac.append(r_v["p_win1"])

out["s_sweep_theta_072_lam_004"] = {
    "s_grid": s_grid,
    "min_step": curve_min,
    "sqrt_step": curve_sqrt,
    "sqrt_plus_vacuum": curve_sqrt_vac
}

# Also test 10-site gated commitment at theta=0.45 and theta=0.72 (lambda=0.02):
out["ten_site_gated"] = {}
for th in [0.45, 0.72]:
    r10_m = simulate_ten_site_exact(theta=th, sigma=0.25, commit_rate_per_step=0.02, vac_eps=0.0, mode="min", n_trials=12000)
    r10_s = simulate_ten_site_exact(theta=th, sigma=0.25, commit_rate_per_step=0.02, vac_eps=0.0, mode="sqrt", n_trials=12000)
    r10_v = simulate_ten_site_exact(theta=th, sigma=0.25, commit_rate_per_step=0.02, vac_eps=0.012, mode="sqrt", n_trials=12000)
    out["ten_site_gated"][str(th)] = {
        "min_p_bright": r10_m["p_bright"],
        "sqrt_p_bright": r10_s["p_bright"],
        "sqrt_plus_vac_p_bright": r10_v["p_bright"]
    }

with open("/Users/john-bramble/.gemini/antigravity/scratch/selector_verification/detailed_results.json", "w") as f:
    json.dump(out, f, indent=2)

# Generate Figure: 2-panel comparison
fig, axes = plt.subplots(1, 2, figsize=(12, 5.2), dpi=200)

# Panel A: Exposed Steps vs Threshold theta (Stalling min vs Conserving sqrt)
thetas = [0.45, 0.72, 0.90, 0.94]
exp_min_2 = [out["two_site_exposure"][str(th)]["min"]["mean_winner_exp_steps"] for th in thetas]
exp_sqrt_2 = [out["two_site_exposure"][str(th)]["sqrt"]["mean_winner_exp_steps"] for th in thetas]
exp_min_10 = [out["ten_site_exposure"][str(th)]["min"]["mean_winner_exp_steps"] for th in thetas]
exp_sqrt_10 = [out["ten_site_exposure"][str(th)]["sqrt"]["mean_winner_exp_steps"] for th in thetas]

ax = axes[0]
ax.plot(thetas, exp_min_2, 'o--', color='#d95f02', lw=2, label='2-site: Stalling min(s1,s2) step (v0.9)')
ax.plot(thetas, exp_sqrt_2, 's-', color='#1b9e77', lw=2.5, label='2-site: Conserving sqrt(s1*s2) step')
ax.plot(thetas, exp_min_10, '^--', color='#e7298a', lw=2, label='10-site: Stalling min(si,sj) step (v0.9)')
ax.plot(thetas, exp_sqrt_10, 'D-', color='#7570b3', lw=2.5, label='10-site: Conserving sqrt(si*sj) step')
ax.set_yscale('log')
ax.set_xlabel('Gap Threshold $\\theta = E_{\\mathrm{gap}} / E_{\\mathrm{photon}}$', fontsize=11)
ax.set_ylabel('Mean Exposed Steps Above $\\theta$ (log scale)', fontsize=11)
ax.set_title('A. Elimination of Zeno Slowdown Above $\\theta$', fontsize=12, fontweight='bold')
ax.grid(True, which='both', alpha=0.3)
ax.legend(fontsize=8.5, loc='upper right')

# Panel B: Deviation from Born across S in [0.1, 0.9]
ax2 = axes[1]
s_arr = np.array(s_grid)
ax2.axhline(0.0, color='black', ls=':', lw=1.5, label='Exact Born Rule ($P_1 - S = 0$)')
ax2.plot(s_arr, np.array(curve_min) - s_arr, 'o--', color='#d95f02', lw=2, label='Stalling min() step (gated $\\theta=0.72$)')
ax2.plot(s_arr, np.array(curve_sqrt) - s_arr, 's-', color='#1b9e77', lw=2, label='Conserving $\\sqrt{s_1 s_2}$ step (gated $\\theta=0.72$)')
ax2.plot(s_arr, np.array(curve_sqrt_vac) - s_arr, 'D-', color='#2b5c8f', lw=2.5, label='Conserving $\\sqrt{s_1 s_2}$ + Non-random vacuum ($\\epsilon=0.018$)')
ax2.set_xlabel('Initial Share $S = |A_1|^2 / (|A_1|^2 + |A_2|^2)$', fontsize=11)
ax2.set_ylabel('Deviation from Born Rule ($P_1 - S$)', fontsize=11)
ax2.set_title('B. Gated Bias & Vacuum Cross-Term Restoration', fontsize=12, fontweight='bold')
ax2.grid(True, alpha=0.3)
ax2.legend(fontsize=8.5, loc='lower left')

plt.tight_layout()
plt.savefig("/Users/john-bramble/.gemini/antigravity/brain/f91ca8ab-6741-4e72-b1f0-af94b132a121/conserving_step_and_vacuum_test.png", dpi=200)
plt.savefig("/Users/john-bramble/.gemini/antigravity/scratch/selector_verification/conserving_step_and_vacuum_test.png", dpi=200)
print(json.dumps(out, indent=2))
