#!/usr/bin/env python3
"""
Generates a single-panel bottom-to-top Feynman-style animated GIF of the
Two-Slit Self-Interference Experiment across 25 lateral detectors (500 individual
atoms, ~1-2 pixels wide per atom):
  - Bottom Zone: Two-slit barrier (equal slit width a, spacing d = 3.5a) and
    intermediate/Fraunhofer cylindrical wavefronts as each photon interferes with itself.
  - Atomic Surface Zone: 500 individual detector atoms (20 atoms per detector well),
    plotting the two-slit incident wave modulus |A_m|, the 1-px atomic electron
    modulus spikes |psi_m(t)| = sqrt(s_m(t)), and the surface quantum vacuum modulus |V_m(t)|.
  - Detector Energy Well Zone: 25 side-by-side detector energy wells; once fixation occurs
    at a single atom m*, the excitation rises upward through that detector's well.
  - Top Zone: Piled click dots above each of the 25 detectors accumulating against the
    two-slit Born-rule distribution curve I(x) ~ sinc^2(beta) cos^2(alpha).
"""

from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Arc
from matplotlib.animation import FuncAnimation, PillowWriter

CODE_DIR = Path(__file__).resolve().parent


def compute_two_slit_profile(num_det=25, atoms_per_det=20, d_over_a=3.5, jitter_x=0.0):
    m_tot = num_det * atoms_per_det
    norm_x = ((np.arange(m_tot, dtype=float) + 0.5) / m_tot - 0.5) * 2.0 + jitter_x
    beta = norm_x * 1.18 * np.pi
    alpha = d_over_a * beta
    env = np.sinc(beta / np.pi)
    interf = np.cos(alpha)
    intensity = (env ** 2) * (interf ** 2) + 0.0015 * (env ** 2)
    i_born = intensity / np.sum(intensity)
    a_mod = np.sqrt(i_born)
    det_born = np.sum(i_born.reshape(num_det, atoms_per_det), axis=1)
    return m_tot, a_mod, i_born, det_born


def simulate_two_slit_photon(
    rng,
    num_det=25,
    atoms_per_det=20,
    d_over_a=3.5,
    theta=0.45,
    sigma=0.28,
    commit_rate=0.075,
    surface_vac_eps=0.0045,
    record_history=True,
):
    m_tot = num_det * atoms_per_det
    jitter_x = rng.uniform(-0.6 / m_tot, 0.6 / m_tot)
    _, a_mod, i_born, det_born = compute_two_slit_profile(num_det, atoms_per_det, d_over_a, jitter_x)

    phi_inc = rng.uniform(0.0, 2.0 * np.pi, size=m_tot)
    phi_det = rng.uniform(0.0, 2.0 * np.pi, size=m_tot)
    phi_surf = rng.uniform(0.0, 2.0 * np.pi, size=(num_det, 1))
    phi_vac = (phi_surf + rng.normal(0.0, 0.5, size=(num_det, atoms_per_det))).ravel()

    factor = 1.0 + 0.20 * np.cos(phi_det - 0.5 * (phi_inc + phi_vac))
    S = i_born * factor
    S /= np.sum(S)

    psi_hist = []
    vac_hist = []
    vac_base = 0.20
    winner_atom = -1

    def record(cur_s, cur_phi, cur_vac):
        if not record_history:
            return
        psi = np.sqrt(np.maximum(0.0, cur_s))
        align = 0.5 * (1.0 + np.cos(cur_phi - cur_vac))
        vac = vac_base * (0.45 + 0.65 * align + 0.45 * psi * align)
        psi_hist.append(psi.copy())
        vac_hist.append(vac.copy())

    record(S, phi_det, phi_vac)

    active = np.where(S > 1e-12)[0]
    for _ in range(1, 550):
        dom_atom = int(active[np.argmax(S[active])])
        max_s = float(S[dom_atom])

        if max_s >= theta:
            above = active[S[active] >= theta]
            rates = commit_rate * S[above]
            r_tot = float(np.sum(rates))
            if rng.random() < r_tot or max_s > 0.99 or len(active) <= 1:
                winner_atom = int(rng.choice(above, p=rates / r_tot)) if r_tot > 0 else dom_atom
                if record_history:
                    for k in range(1, 4):
                        alpha_c = k / 3.0
                        s_col = (1.0 - alpha_c) * S
                        s_col[winner_atom] = (1.0 - alpha_c) * S[winner_atom] + alpha_c * 1.0
                        phi_vac[winner_atom] = phi_det[winner_atom]
                        record(s_col, phi_det, phi_vac)
                break

        dom_phi = phi_det[dom_atom]
        phi_vac[active] += 0.30 * np.sin(dom_phi - phi_vac[active]) + rng.normal(0.0, 0.35, size=len(active))
        phi_det[active] += 0.15 * np.sin(phi_vac[active] - phi_det[active]) + rng.uniform(-np.pi, np.pi, size=len(active))

        perm = rng.permutation(active)
        n_act = len(perm)
        for p in range(0, n_act - 1, 2):
            i, j = int(perm[p]), int(perm[p + 1])
            tot = S[i] + S[j]
            if tot <= 1e-14:
                continue
            u = S[i] / tot
            if 0.0 < u < 1.0 and surface_vac_eps != 0.0:
                u = float(np.clip(u + surface_vac_eps * np.sqrt(u * (1.0 - u)) * (np.sqrt(1.0 - u) - np.sqrt(u)), 0.0, 1.0))
            dphi = (phi_det[i] - phi_vac[i]) - (phi_det[j] - phi_vac[j])
            sign = 1.0 if np.cos(dphi) >= 0.0 else -1.0
            std = sigma * np.sqrt(u * (1.0 - u))
            min_side = min(u, 1.0 - u)
            if std <= min_side:
                u_new = u + sign * std
            else:
                p_abs = (std ** 2) / (std ** 2 + min_side ** 2)
                if rng.random() < p_abs:
                    u_new = 0.0 if u < 0.5 else 1.0
                else:
                    step_away = (std ** 2) / max(min_side, 1e-15)
                    u_new = min(1.0, u + step_away) if u < 0.5 else max(0.0, u - step_away)
            u_new = float(np.clip(u_new, 0.0, 1.0))
            S[i] = tot * u_new
            S[j] = tot * (1.0 - u_new)

        active = active[S[active] > 1e-12]
        record(S, phi_det, phi_vac)

    if winner_atom == -1:
        winner_atom = int(np.argmax(S))

    return {
        "a_mod": a_mod,
        "det_born": det_born,
        "psi_hist": np.asarray(psi_hist) if record_history else None,
        "vac_hist": np.asarray(vac_hist) if record_history else None,
        "winner_atom": winner_atom,
        "winner_det": winner_atom // atoms_per_det,
    }


def build_two_slit_gif():
    rng = np.random.default_rng(20261042)
    num_det = 25
    atoms_per_det = 20
    m_tot = num_det * atoms_per_det
    d_over_a = 3.5
    theta = 0.45
    sqrt_theta = float(np.sqrt(theta))

    _, a_mod_ref, _, det_born = compute_two_slit_profile(num_det, atoms_per_det, d_over_a, 0.0)

    # Seed with 75 prior photon clicks so the multi-fringe two-slit pattern is visible and growing
    det_clicks = np.zeros(num_det, dtype=int)
    for _ in range(75):
        res = simulate_two_slit_photon(rng, num_det=num_det, atoms_per_det=atoms_per_det, d_over_a=d_over_a, theta=theta, record_history=False)
        det_clicks[res["winner_det"]] += 1

    frames = []
    for p_idx in range(4):
        run = simulate_two_slit_photon(rng, num_det=num_det, atoms_per_det=atoms_per_det, d_over_a=d_over_a, theta=theta, record_history=True)
        n_steps = len(run["psi_hist"])
        counts_before = det_clicks.copy()

        # 1. Wave propagation from the two slits across the Fraunhofer zone (4 frames)
        for k in range(4):
            frames.append({
                "p_num": 76 + p_idx,
                "stage": "propagate",
                "prog": (k + 1) / 4.0,
                "run": run,
                "step": 0,
                "counts": counts_before,
                "just_clicked": False,
            })

        # 2. 500-atom conserving martingale exchange at the detector entrance (10 frames)
        sub_steps = np.unique(
            np.concatenate([
                np.linspace(0, max(0, n_steps - 4), 8, dtype=int),
                np.arange(max(0, n_steps - 3), n_steps, dtype=int),
            ])
        )
        for s_i in sub_steps:
            frames.append({
                "p_num": 76 + p_idx,
                "stage": "exchange",
                "prog": 1.0,
                "run": run,
                "step": int(s_i),
                "counts": counts_before,
                "just_clicked": False,
            })

        # 3. Upward excitation transit inside the winning detector's energy well (4 frames)
        for k in range(4):
            frames.append({
                "p_num": 76 + p_idx,
                "stage": "ascend",
                "prog": (k + 1) / 4.0,
                "run": run,
                "step": n_steps - 1,
                "counts": counts_before,
                "just_clicked": False,
            })

        # 4. Click dot added to pile above winning detector (3 frames)
        det_clicks[run["winner_det"]] += 1
        counts_after = det_clicks.copy()
        for _ in range(3):
            frames.append({
                "p_num": 76 + p_idx,
                "stage": "recorded",
                "prog": 1.0,
                "run": run,
                "step": n_steps - 1,
                "counts": counts_after,
                "just_clicked": True,
            })

    x_atoms = np.arange(m_tot, dtype=float) + 0.5  # 0.5 .. 499.5
    x_dets = (np.arange(num_det, dtype=float) + 0.5) * atoms_per_det

    fig, ax = plt.subplots(figsize=(11.6, 7.5), dpi=108)
    fig.patch.set_facecolor("#0f172a")

    def update(f_idx):
        ax.clear()
        ax.set_facecolor("#1e293b")
        for spine in ax.spines.values():
            spine.set_color("#334155")
        ax.set_xlim(-22.0, m_tot + 8.0)
        ax.set_ylim(0.0, 10.6)
        ax.set_yticks([0.5, 1.5, 2.8, 4.3, 5.8, 6.8, 8.0])
        ax.set_yticklabels(["Two Slits", "Fraunhofer ↑", "Modulus 0.0", "Modulus 0.5", "Modulus 1.0", "25 Wells ↑", "Click Pile"], color="#cbd5e1", fontsize=8.5)
        ax.set_xticks([10, 110, 250, 390, 490])
        ax.set_xticklabels(["Det 1 (Atoms 1–20)", "Det 6", "Det 13 (Central Fringe)", "Det 20", "Det 25 (Atoms 481–500)"], color="#cbd5e1", fontsize=8.5)

        f = frames[f_idx]
        run = f["run"]
        step = f["step"]
        stage = f["stage"]
        prog = f["prog"]
        counts = f["counts"]
        win_atom = run["winner_atom"]
        win_det = run["winner_det"]
        psi = run["psi_hist"][step]
        vac = run["vac_hist"][step]
        a_mod = run["a_mod"]
        tot_clicks = int(np.sum(counts))

        # Left Margin: Upward Arrow of Time (t ↑)
        ax.annotate("", xy=(-11.0, 10.1), xytext=(-11.0, 0.35), arrowprops=dict(arrowstyle="->", color="#38bdf8", lw=2.5))
        ax.text(-17.5, 5.2, "ARROW OF TIME (t ↑)", color="#38bdf8", fontsize=9.5, fontweight="bold", rotation=90, va="center")

        # --- ZONE 1 (BOTTOM, y = 0.2 .. 2.6): TWO-SLIT BARRIER & FRAUNHOFER SELF-INTERFERENCE ---
        slit_y = 0.55
        prop_top = 2.65
        slit_w = 14.0
        slit_sep = slit_w * d_over_a  # 49.0
        cx = m_tot * 0.5              # 250.0
        s1_x = cx - 0.5 * slit_sep
        s2_x = cx + 0.5 * slit_sep

        # Opaque Two-Slit Screen
        ax.plot([0, s1_x - 0.5 * slit_w], [slit_y, slit_y], color="#cbd5e1", linewidth=4.0)
        ax.plot([s1_x + 0.5 * slit_w, s2_x - 0.5 * slit_w], [slit_y, slit_y], color="#cbd5e1", linewidth=4.0)
        ax.plot([s2_x + 0.5 * slit_w, m_tot], [slit_y, slit_y], color="#cbd5e1", linewidth=4.0)
        ax.text(s1_x, 0.18, "Slit 1 (width a)", color="#38bdf8", fontsize=8.5, ha="center", fontweight="bold")
        ax.text(s2_x, 0.18, f"Slit 2 (width a, spacing d={d_over_a:.1f}a)", color="#38bdf8", fontsize=8.5, ha="center", fontweight="bold")

        # Cylindrical wavefront arcs from Slit 1 and Slit 2 across the intermediate/Fraunhofer zone
        reach = (0.2 + 0.8 * prog) if stage == "propagate" else 1.0
        for r_frac in np.linspace(0.18, reach, 6):
            rx = r_frac * 240.0
            ry = r_frac * (prop_top - slit_y)
            arc1 = Arc((s1_x, slit_y), width=2.0 * rx, height=2.0 * ry, angle=0, theta1=8, theta2=172, color="#38bdf8", alpha=0.28, lw=1.2)
            arc2 = Arc((s2_x, slit_y), width=2.0 * rx, height=2.0 * ry, angle=0, theta1=8, theta2=172, color="#38bdf8", alpha=0.28, lw=1.2)
            ax.add_patch(arc1)
            ax.add_patch(arc2)

        ax.text(cx, 2.38, "↑ Single Photon Self-Interfering Across Intermediate / Fraunhofer Zone ↑", color="#94a3b8", fontsize=8.5, ha="center")

        # --- ZONE 2 (y = 2.8 .. 5.8): 500 DETECTOR ATOMS (1-px WIDE EACH), |A_m|, |psi_m(t)|, |V_m(t)| ---
        atom_base = 2.8
        atom_scale = 3.0
        y_thresh = atom_base + sqrt_theta * atom_scale
        ax.axhline(y_thresh, xmin=0.04, xmax=0.98, color="#f43f5e", linestyle="--", linewidth=1.3, label=rf"Single-Atom Threshold $\sqrt{{\theta}}={sqrt_theta:.2f}$")

        # Shaded Two-Slit Incident Wave Modulus Envelope |A_m| (scaled for visibility)
        vis_a = np.minimum(0.92, a_mod * (0.42 / np.max(a_mod)))
        ax.fill_between(x_atoms, atom_base, atom_base + vis_a * atom_scale, color="#38bdf8", alpha=0.24, label=r"Two-Slit Self-Interference Wave $|A_m|$")
        ax.plot(x_atoms, atom_base + vis_a * atom_scale, color="#38bdf8", linewidth=1.2)

        # Surface Quantum Vacuum Modulus |V_m(t)|
        ax.plot(x_atoms, atom_base + vac * atom_scale, color="#c084fc", linewidth=1.0, alpha=0.85, label=r"Surface Quantum Vacuum $|V_m(t)|$")

        # 500 Individual 1-px wide Atomic Electron Modulus Spikes |psi_m(t)| = sqrt(s_m(t))
        psi_vis = (0.25 + 0.75 * prog) if stage == "propagate" else 1.0
        vals = psi * psi_vis
        active_m = np.where(vals > 0.003)[0]
        ax.vlines(x_atoms[active_m], atom_base, atom_base + vals[active_m] * atom_scale, colors="#10b981", linewidth=1.2, label=r"500 Atomic Electrons $|\psi_m(t)| = \sqrt{s_m(t)}$")

        if stage in ("ascend", "recorded") or vals[win_atom] >= sqrt_theta:
            ax.vlines([x_atoms[win_atom]], atom_base, atom_base + vals[win_atom] * atom_scale, colors="#f59e0b", linewidth=2.6)

        # --- ZONE 3 (y = 5.95 .. 7.25): 25 DETECTOR ENERGY WELLS & UPWARD EXCITATION PULSE ---
        well_bot, well_top = 5.95, 7.25
        for d in range(num_det + 1):
            xw = d * atoms_per_det
            ax.plot([xw, xw], [well_bot, well_top], color="#94a3b8", linewidth=1.2, alpha=0.65)
            ax.plot([xw, xw], [atom_base, atom_base + atom_scale], color="#64748b", linestyle=":", linewidth=0.6, alpha=0.35)

        if stage == "ascend":
            wx_start = x_atoms[win_atom]
            wx_end = x_dets[win_det]
            wx = wx_start + prog * (wx_end - wx_start)
            wy = well_bot + prog * (well_top - well_bot)
            ax.fill_between([win_det * atoms_per_det, (win_det + 1) * atoms_per_det], well_bot, well_top, color="#f59e0b", alpha=0.18)
            ax.plot([wx_start, wx], [atom_base + atom_scale, wy], color="#fbbf24", linestyle="--", linewidth=1.8)
            ax.scatter([wx], [wy], s=120, color="#fde047", edgecolors="#f59e0b", linewidths=1.8, zorder=6)

        # --- ZONE 4 (TOP, y = 7.45 .. 10.2): PILED DETECTOR CLICKS & TWO-SLIT BORN CURVE ---
        pile_base = 7.45
        dot_step = 0.14
        for d in range(num_det):
            c = counts[d]
            if c > 0:
                ys = pile_base + np.arange(1, c + 1) * dot_step
                cols = ["#fde047" if (d == win_det and f["just_clicked"] and idx == c - 1) else "#f59e0b" for idx in range(c)]
                ax.scatter(np.full(c, x_dets[d]), ys, c=cols, s=28, zorder=5, edgecolors="none")

        # Smooth Born Two-Slit Interference Curve across the top zone
        born_curve_y = pile_base + (tot_clicks * (a_mod_ref ** 2) * atoms_per_det) * dot_step
        ax.plot(x_atoms, born_curve_y, color="#38bdf8", linewidth=2.2, label=r"Born Rule Two-Slit Fringe Line ($I(x) \propto \mathrm{sinc}^2\cdot\cos^2$)")

        ax.text(cx, 10.15, f"Piled Clicks Across 25 Detectors (Total = {tot_clicks} Photons) Approximating the Two-Slit Born Distribution",
                color="#f8fafc", fontsize=9.5, fontweight="bold", ha="center")

        ax.set_title(
            f"Two-Slit Self-Interference Across 25 Detectors (500 Atoms, 1-px Resolution)  |  Photon #{f['p_num']}  |  "
            f"Σ|ψ_m|² = {float(np.sum(psi**2)):.4f}",
            color="#f8fafc", fontsize=10.5, fontweight="bold",
        )
        ax.legend(loc="lower right", bbox_to_anchor=(0.99, 0.24), fontsize=8, facecolor="#0f172a", edgecolor="#334155", labelcolor="#e2e8f0", ncol=2)

    anim = FuncAnimation(fig, update, frames=len(frames), interval=85, blit=False)
    out_gif = CODE_DIR / "fig_v2_6_two_slit_multidetector_animation.gif"
    anim.save(out_gif, writer=PillowWriter(fps=12))
    plt.close(fig)
    print(f"Saved Two-Slit 25-detector (500-atom) animation to {out_gif}")


if __name__ == "__main__":
    build_two_slit_gif()
