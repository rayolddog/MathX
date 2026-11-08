#!/usr/bin/env python3
"""
Generates a single-panel Feynman-style bottom-to-top animated GIF
(fig_v2_5_feynman_born_accumulation.gif) across 20 atoms in Detector 1 and
20 atoms in Detector 2 (40 atoms total):
  - Arrow of Time flows from BOTTOM to TOP (t ↑).
  - Bottom Zone: Incident particle enters as an upward-propagating wave packet |A_m|.
  - Atomic Surface Zone: 20+20 detector atoms with time-unconstrained virtual quantum
    vacuum loops |V_m| and conserving atomic electron modulus |psi_m(t)| = sqrt(s_m(t))
    exchanging energy until one atom crosses the threshold sqrt(theta).
  - Detector Energy Well Zone: Once fixation occurs at atom m*, the excitation travels
    upward through the energy well bounded by the side walls of that detector.
  - Top Readout Zone: At the top of the detector, a click dot is recorded and piled
    above the fixation site, gradually approximating the Born-rule distribution line.
"""

from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse, FancyArrowPatch
from matplotlib.animation import FuncAnimation, PillowWriter

CODE_DIR = Path(__file__).resolve().parent


def simulate_single_event(
    rng,
    s_region1=0.70,
    n_atoms_per_det=20,
    theta=0.45,
    sigma=0.25,
    commit_rate=0.06,
    surface_vac_eps=0.005,
    record_history=True,
):
    M1 = n_atoms_per_det
    M2 = n_atoms_per_det
    M_tot = M1 + M2

    x1 = np.arange(M1, dtype=float) - (M1 - 1) / 2.0
    x2 = np.arange(M2, dtype=float) - (M2 - 1) / 2.0
    off1 = rng.uniform(-0.45, 0.45)
    off2 = rng.uniform(-0.45, 0.45)
    waist = 0.55 * M1
    env1 = np.exp(-((x1 - off1) ** 2) / (2.0 * waist ** 2))
    env2 = np.exp(-((x2 - off2) ** 2) / (2.0 * waist ** 2))

    A_mod = np.zeros(M_tot, dtype=float)
    A_mod[:M1] = np.sqrt(s_region1 * (env1 / np.sum(env1)))
    A_mod[M1:] = np.sqrt((1.0 - s_region1) * (env2 / np.sum(env2)))

    phi_inc = rng.uniform(0.0, 2.0 * np.pi, size=M_tot)
    phi_det = rng.uniform(0.0, 2.0 * np.pi, size=M_tot)
    phi_surf1 = rng.uniform(0.0, 2.0 * np.pi)
    phi_surf2 = rng.uniform(0.0, 2.0 * np.pi)
    phi_vac = np.zeros(M_tot, dtype=float)
    phi_vac[:M1] = phi_surf1 + rng.normal(0.0, 0.5, size=M1)
    phi_vac[M1:] = phi_surf2 + rng.normal(0.0, 0.5, size=M2)

    align1 = env1 * (1.0 + 0.20 * np.cos(phi_det[:M1] - 0.5 * (phi_inc[:M1] + phi_vac[:M1])))
    align2 = env2 * (1.0 + 0.20 * np.cos(phi_det[M1:] - 0.5 * (phi_inc[M1:] + phi_vac[M1:])))

    S = np.zeros(M_tot, dtype=float)
    S[:M1] = s_region1 * (align1 / np.sum(align1))
    S[M1:] = (1.0 - s_region1) * (align2 / np.sum(align2))

    psi_hist = []
    vac_hist = []
    vph_hist = []
    vac_base = 0.24
    winner_atom = -1

    def record(cur_s, cur_phi, cur_vac):
        if not record_history:
            return
        psi = np.sqrt(np.maximum(0.0, cur_s))
        phase_align = 0.5 * (1.0 + np.cos(cur_phi - cur_vac))
        vac = vac_base * (0.55 + 0.65 * phase_align + 0.35 * psi * phase_align)
        psi_hist.append(psi.copy())
        vac_hist.append(vac.copy())
        vph_hist.append(cur_vac.copy())

    record(S, phi_det, phi_vac)

    for _ in range(1, 550):
        dom_atom = int(np.argmax(S))
        max_s = float(S[dom_atom])

        if max_s >= theta:
            rates = np.where(S >= theta, commit_rate * S, 0.0)
            r_tot = float(np.sum(rates))
            if rng.random() < r_tot or max_s > 0.99:
                winner_atom = int(rng.choice(M_tot, p=rates / r_tot)) if r_tot > 0 else dom_atom
                if record_history:
                    for k in range(1, 4):
                        alpha = k / 3.0
                        s_col = (1.0 - alpha) * S
                        s_col[winner_atom] = (1.0 - alpha) * S[winner_atom] + alpha * 1.0
                        phi_vac[winner_atom] = phi_det[winner_atom]
                        record(s_col, phi_det, phi_vac)
                break

        dom_phi = phi_det[dom_atom]
        phi_vac = phi_vac + 0.30 * np.sin(dom_phi - phi_vac) + rng.normal(0.0, 0.35, size=M_tot)
        phi_det = phi_det + 0.15 * np.sin(phi_vac - phi_det) + rng.uniform(-np.pi, np.pi, size=M_tot)

        perm = rng.permutation(M_tot)
        for p in range(0, M_tot - 1, 2):
            i, j = int(perm[p]), int(perm[p + 1])
            tot = S[i] + S[j]
            if tot <= 1e-12:
                continue
            u = S[i] / tot
            if 0.0 < u < 1.0 and surface_vac_eps != 0.0:
                u = np.clip(u + surface_vac_eps * np.sqrt(u * (1.0 - u)) * (np.sqrt(1.0 - u) - np.sqrt(u)), 0.0, 1.0)
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

        record(S, phi_det, phi_vac)

    if winner_atom == -1:
        winner_atom = int(np.argmax(S))

    return {
        "A_mod": A_mod,
        "psi_hist": np.asarray(psi_hist) if record_history else None,
        "vac_hist": np.asarray(vac_hist) if record_history else None,
        "vph_hist": np.asarray(vph_hist) if record_history else None,
        "winner_atom": winner_atom,
    }


def build_feynman_animation():
    rng = np.random.default_rng(20261035)
    M1 = 20
    M2 = 20
    M_tot = 40
    S1 = 0.70
    theta = 0.45
    sqrt_theta = float(np.sqrt(theta))

    # Born expected distribution across the 40 atoms
    x = np.arange(M1, dtype=float) - (M1 - 1) / 2.0
    waist = 0.55 * M1
    w = np.exp(-(x ** 2) / (2.0 * waist ** 2))
    w /= np.sum(w)
    born_probs = np.concatenate([S1 * w, (1.0 - S1) * w])

    # Seed with 32 prior clicks so the piled dots already show the emerging shape, then animate 4 sequential particles
    click_counts = np.zeros(M_tot, dtype=int)
    for _ in range(32):
        res = simulate_single_event(rng, s_region1=S1, n_atoms_per_det=M1, theta=theta, record_history=False)
        click_counts[res["winner_atom"]] += 1

    # X positions with a gap between Detector 1 (1..20) and Detector 2 (23..42)
    x_pos = np.array([m + 1 if m < M1 else m + 3 for m in range(M_tot)], dtype=float)

    frames = []
    for p_idx in range(4):
        run = simulate_single_event(rng, s_region1=S1, n_atoms_per_det=M1, theta=theta, record_history=True)
        n_steps = len(run["psi_hist"])
        counts_before = click_counts.copy()

        # Stage 1: Incoming wave packet rising from bottom (4 frames)
        for k in range(4):
            frames.append({
                "p_num": 33 + p_idx,
                "stage": "incoming",
                "prog": (k + 1) / 4.0,
                "run": run,
                "step": 0,
                "counts": counts_before,
                "just_clicked": False,
            })

        # Stage 2: Pre-modulus martingale exchange at the 40 detector atoms (10 frames)
        sub_steps = np.unique(
            np.concatenate([
                np.linspace(0, max(0, n_steps - 4), 8, dtype=int),
                np.arange(max(0, n_steps - 3), n_steps, dtype=int),
            ])
        )
        for s_i in sub_steps:
            frames.append({
                "p_num": 33 + p_idx,
                "stage": "exchange",
                "prog": 1.0,
                "run": run,
                "step": int(s_i),
                "counts": counts_before,
                "just_clicked": False,
            })

        # Stage 3: Fixation excitation traveling upward through the detector energy well (5 frames)
        for k in range(5):
            frames.append({
                "p_num": 33 + p_idx,
                "stage": "ascend",
                "prog": (k + 1) / 5.0,
                "run": run,
                "step": n_steps - 1,
                "counts": counts_before,
                "just_clicked": False,
            })

        # Stage 4: Click dot recorded at top of detector energy well (3 frames)
        click_counts[run["winner_atom"]] += 1
        counts_after = click_counts.copy()
        for _ in range(3):
            frames.append({
                "p_num": 33 + p_idx,
                "stage": "recorded",
                "prog": 1.0,
                "run": run,
                "step": n_steps - 1,
                "counts": counts_after,
                "just_clicked": True,
            })

    # Single unified vertical Feynman-style canvas (NO Panel B!)
    fig, ax = plt.subplots(figsize=(11.5, 7.2), dpi=108)
    fig.patch.set_facecolor("#0f172a")

    # Vertical coordinate system y in [0.0, 10.0] (Arrow of Time bottom -> top):
    # y = 0.2 .. 1.8 : Incoming Particle Wave Packet |A_m| propagating upward
    # y = 2.1 .. 5.1 : 20+20 Detector Atoms, Electron Modulus |psi_m(t)|, Virtual Vacuum Loops |V_m|, Threshold sqrt(theta)
    # y = 5.2 .. 7.0 : Detector 1 & Detector 2 Conduction Energy Wells (side walls + upward excitation propagation)
    # y = 7.2 .. 9.6 : Piled Clicks (Dots) above each fixation site + Born Rule Distribution Line

    def update(f_idx):
        ax.clear()
        ax.set_facecolor("#1e293b")
        for spine in ax.spines.values():
            spine.set_color("#334155")
        ax.set_xlim(-1.8, 43.8)
        ax.set_ylim(0.0, 10.1)
        ax.set_yticks([2.1, 3.6, 5.1, 6.1, 7.3])
        ax.set_yticklabels(["Modulus 0.0", "Modulus 0.5", "Modulus 1.0", "Energy Well ↑", "Click Pile"], color="#cbd5e1", fontsize=8.5)
        ax.set_xticks([1, 10, 20, 23, 32, 42])
        ax.set_xticklabels(["Atom 1", "Atom 10", "Atom 20", "Atom 21", "Atom 30", "Atom 40"], color="#cbd5e1", fontsize=8.5)

        f = frames[f_idx]
        run = f["run"]
        step = f["step"]
        stage = f["stage"]
        prog = f["prog"]
        counts = f["counts"]
        winner = run["winner_atom"]
        psi = run["psi_hist"][step]
        vac = run["vac_hist"][step]
        vph = run["vph_hist"][step]
        A_mod = run["A_mod"]
        tot_clicks = int(np.sum(counts))
        det1_pct = 100.0 * np.sum(counts[:M1]) / max(1, tot_clicks)

        # Left Margin: Upward Arrow of Time (t ↑)
        ax.annotate(
            "", xy=(-0.8, 9.6), xytext=(-0.8, 0.4),
            arrowprops=dict(arrowstyle="->", color="#38bdf8", lw=2.5),
        )
        ax.text(-1.35, 5.0, "ARROW OF TIME (t ↑)", color="#38bdf8", fontsize=9.5, fontweight="bold", rotation=90, va="center")

        # --- ZONE 1 (BOTTOM, y = 0.2 .. 1.8): INCOMING WAVE PACKET |A_m| ---
        wave_base_y = 0.35 + (prog if stage == "incoming" else 1.0) * 0.85
        for det_slice in (slice(0, M1), slice(M1, M_tot)):
            xs = x_pos[det_slice]
            osc = 0.08 * np.sin(0.7 * xs - prog * 6.0)
            yw = wave_base_y + A_mod[det_slice] * 1.45 + osc
            ax.fill_between(xs, wave_base_y - 0.15, yw, color="#38bdf8", alpha=0.28)
            ax.plot(xs, yw, color="#38bdf8", linewidth=1.8)
        ax.text(21.5, 0.18, "↑ Incoming Particle Wave Packet |A_m| (Propagating Upward Into Detector 1 & 2 Apertures)",
                color="#38bdf8", fontsize=9, fontweight="bold", ha="center")

        # --- ZONE 2 (y = 2.1 .. 5.1): 40 DETECTOR ATOMS, ELECTRON MODULUS & VIRTUAL VACUUM LOOPS ---
        atom_base = 2.1
        atom_scale = 3.0  # maps modulus [0, 1] -> height [0, 3.0]
        y_thresh = atom_base + sqrt_theta * atom_scale
        ax.hlines(y_thresh, 0.5, 20.5, colors="#f43f5e", linestyles="--", linewidth=1.4, label=rf"Threshold $\sqrt{{\theta}}={sqrt_theta:.2f}$")
        ax.hlines(y_thresh, 22.5, 42.5, colors="#f43f5e", linestyles="--", linewidth=1.4)

        # Incident wave modulus backdrop bars |A_m|
        ax.bar(x_pos, A_mod * atom_scale, bottom=atom_base, width=0.68, color="#38bdf8", alpha=0.25, edgecolor="#38bdf8", linewidth=0.8, label=r"Incident Wave $|A_m|$")

        # Atomic electron modulus bars |psi_m(t)|
        psi_scale = (0.3 + 0.7 * prog) if stage == "incoming" else 1.0
        bar_cols = [
            "#f59e0b" if (m == winner and stage in ("ascend", "recorded")) else ("#fbbf24" if psi[m] * psi_scale >= sqrt_theta else "#10b981")
            for m in range(M_tot)
        ]
        ax.bar(x_pos, psi * psi_scale * atom_scale, bottom=atom_base, width=0.42, color=bar_cols, label=r"Atomic Electron $|\psi_m(t)| = \sqrt{s_m(t)}$")

        # Time-unconstrained Virtual Quantum Vacuum Loops (ellipses + surface envelope |V_m|)
        for m in range(0, M_tot, 2):
            if m + 1 < M_tot and (m + 1 != M1):
                cx = 0.5 * (x_pos[m] + x_pos[m + 1])
                cy = atom_base + vac[m] * atom_scale
                ell = Ellipse((cx, cy), width=0.95, height=0.22 + 0.18 * abs(np.sin(vph[m])),
                              edgecolor="#c084fc", facecolor="none", alpha=0.45, linewidth=1.1)
                ax.add_patch(ell)
        ax.plot(x_pos[:M1], atom_base + vac[:M1] * atom_scale, "-o", color="#c084fc", markersize=3.0, linewidth=1.5, label=r"Virtual Vacuum $|V_m|$ (Time-Unconstrained)")
        ax.plot(x_pos[M1:], atom_base + vac[M1:] * atom_scale, "-o", color="#c084fc", markersize=3.0, linewidth=1.5)

        # --- ZONE 3 (y = 5.2 .. 7.0): DETECTOR ENERGY WELLS & UPWARD EXCITATION PULSE ---
        well_bot, well_top = 5.2, 7.0
        for xL, xR, lbl in ((0.4, 20.6, "Detector 1 Energy Well ↑"), (22.4, 42.6, "Detector 2 Energy Well ↑")):
            ax.fill_between([xL, xR], well_bot, well_top, color="#0f172a", alpha=0.45)
            ax.plot([xL, xL], [well_bot, well_top], color="#94a3b8", linewidth=2.8)
            ax.plot([xR, xR], [well_bot, well_top], color="#94a3b8", linewidth=2.8)
            ax.text(0.5 * (xL + xR), 6.65, lbl, color="#94a3b8", fontsize=8.5, ha="center")

        if stage == "ascend":
            wx = x_pos[winner]
            wy = well_bot + prog * (well_top - well_bot)
            ax.plot([wx, wx], [atom_base + atom_scale, wy], color="#fbbf24", linestyle="--", linewidth=2.0)
            ax.scatter([wx], [wy], s=140, color="#fde047", edgecolors="#f59e0b", linewidths=2.0, zorder=6)
            xL, xR = (0.4, 20.6) if winner < M1 else (22.4, 42.6)
            xs_arc = np.linspace(xL + 0.2, xR - 0.2, 30)
            ys_arc = wy - 0.18 * ((xs_arc - wx) / 10.0) ** 2
            ax.plot(xs_arc, ys_arc, color="#fde047", linestyle=":", linewidth=1.5, alpha=0.75)

        # --- ZONE 4 (TOP, y = 7.2 .. 9.6): PILED DETECTION CLICKS & BORN DISTRIBUTION LINE ---
        pile_base = 7.25
        dot_step = 0.24
        for m in range(M_tot):
            c = counts[m]
            if c > 0:
                ys = pile_base + np.arange(1, c + 1) * dot_step
                cols = [
                    "#fde047" if (m == winner and f["just_clicked"] and idx == c - 1) else ("#f59e0b" if m < M1 else "#38bdf8")
                    for idx in range(c)
                ]
                ax.scatter(np.full(c, x_pos[m]), ys, c=cols, s=36, zorder=5, edgecolors="none")

        # Expected Born Rule Distribution Line across the top
        ax.plot(x_pos[:M1], pile_base + (tot_clicks * born_probs[:M1]) * dot_step, color="#38bdf8", linewidth=2.2, label="Born Rule Distribution Line")
        ax.plot(x_pos[M1:], pile_base + (tot_clicks * born_probs[M1:]) * dot_step, color="#38bdf8", linewidth=2.2)

        ax.text(10.5, 9.65, f"Detector 1 Piled Clicks: {int(np.sum(counts[:M1]))} ({det1_pct:.1f}% vs Born 70.0%)",
                color="#f8fafc", fontsize=9.5, fontweight="bold", ha="center")
        ax.text(32.5, 9.65, f"Detector 2 Piled Clicks: {int(np.sum(counts[M1:]))} ({100.0-det1_pct:.1f}% vs Born 30.0%)",
                color="#f8fafc", fontsize=9.5, fontweight="bold", ha="center")

        ax.set_title(
            f"Feynman-Oriented Bottom-to-Top Detector Selection  |  Particle #{f['p_num']}  |  "
            f"Conserved Σ|ψ_m|² = {float(np.sum(psi**2)):.4f}",
            color="#f8fafc", fontsize=11, fontweight="bold",
        )
        ax.legend(loc="lower right", bbox_to_anchor=(0.99, 0.21), fontsize=8, facecolor="#0f172a", edgecolor="#334155", labelcolor="#e2e8f0", ncol=2)

    anim = FuncAnimation(fig, update, frames=len(frames), interval=85, blit=False)
    out_gif = CODE_DIR / "fig_v2_5_feynman_born_accumulation.gif"
    anim.save(out_gif, writer=PillowWriter(fps=12))
    plt.close(fig)
    print(f"Saved Feynman bottom-to-top animation to {out_gif}")


if __name__ == "__main__":
    build_feynman_animation()
