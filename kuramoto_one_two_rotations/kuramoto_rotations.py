#!/usr/bin/env python3
"""
kuramoto_rotations.py: lock within one or two rotations of the clock
====================================================================

Definitions are fixed in PREREGISTRATION.md (read it first).

  clock       dφ/dt = ω + K sin(Φ_bulk − φ),  Φ_bulk = ω_bulk t      (EQUATIONS.md §2 form)
  error       θ = φ − Φ_bulk,   dθ/dt = Δ − K sin θ,   Δ = ω − ω_bulk (Adler form)
  rotation    one 2π turn of the clock, T = 2π/ω;  x = K/ω;  δ = Δ/ω
  locked      |θ − θ*| ≤ ε,  θ* = arcsin(Δ/K);  ε = 0.35 rad primary, 0.10 and π/4 sensitivity
  in n rot.   t_ε ≤ nT,  n ∈ {1, 2}
  likelihood  P_n = fraction of θ₀ ~ Uniform(−π, π] with t_ε ≤ nT   (floor ε/π at K = 0)

Method. The flow is one-dimensional and monotone between fixed points, so the time to the
band is a difference of the antiderivative F(θ) of 1/(Δ − K sin θ), which has a closed form.
F is checked against RK45 integration on a random sample, and the mutual two-clock case
(K₁, K₂ pulling each other) against the Adler equation with K_eff = K₁ + K₂.

Origin: the equation and the integration pattern come from
DiracKuramotoFramework/tests/kuramoto_sync.py (copy in ../dk_origin/).
Time unit 1/ω, so T = 2π and K = x, Δ = δ.
"""
import csv, sys, datetime, platform
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, BoundaryNorm

HERE = Path(__file__).resolve().parent
OUT = HERE / "results"
OUT.mkdir(exist_ok=True)

OMEGA = 1.0
T_ROT = 2 * np.pi / OMEGA
EPS_PRIMARY = 0.35
EPS_LIST = [0.10, 0.35, np.pi / 4]
EPS_LABEL = {0.10: "0.10", 0.35: "0.35", np.pi / 4: "pi/4"}
N_LIST = [1, 2]
DELTA_LIST = [0.0, 0.1, 0.25, 0.5]
X_TABLE = [0.02, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0, 1.25, 1.5]
P_TARGETS = [0.5, 0.9, 0.95, 0.99]

# ----------------------------------------------------------------------------- closed form
def wrap(a):
    """Map an angle to (−π, π]."""
    return (a + np.pi) % (2 * np.pi) - np.pi

def fixed_points(K, D):
    ts = np.arcsin(D / K)
    return ts, np.pi - ts

def F(theta, K, D):
    """Antiderivative of 1/(Δ − K sin θ) for |Δ| < K, so F(θ(t)) = F(θ₀) + t along the flow.
    Built on tan(θ/2): 2π-periodic; +∞ at the stable point, −∞ at the unstable one."""
    t = np.tan(np.asarray(theta, dtype=float) / 2)
    with np.errstate(divide="ignore", invalid="ignore"):
        if D == 0:
            return -np.log(np.abs(t)) / K
        Om = np.sqrt(K * K - D * D)
        tp = (K + Om) / D          # tan(θ_u/2)
        tm = (K - Om) / D          # tan(θ*/2)
        return np.log(np.abs((t - tp) / (t - tm))) / Om

def time_to_band(theta0, K, D, eps):
    """First time |θ − θ*| ≤ ε starting from θ₀ (vectorised). inf where no lock is possible."""
    theta0 = np.asarray(theta0, dtype=float)
    if K <= 0 or abs(D) >= K:
        if K == 0 and D == 0:                       # frozen clock: the criterion's floor only
            return np.where(np.abs(wrap(theta0)) <= eps, 0.0, np.inf)
        return np.full(theta0.shape, np.inf)
    ts, tu = fixed_points(K, D)
    if eps >= tu - ts:                              # band would swallow the unstable point
        return np.full(theta0.shape, np.nan)
    d = wrap(theta0 - ts)                           # signed error from θ* on the circle
    inband = np.abs(d) <= eps
    above = (d > eps) & (d < tu - ts)               # arc θ* → θ_u: flows down to θ* + ε
    below = ~inband & ~above                        # the rest: flows up to θ* − ε
    F0 = F(theta0, K, D)
    t = np.empty(theta0.shape)
    t[inband] = 0.0
    t[above] = F(ts + eps, K, D) - F0[above]
    t[below] = F(ts - eps, K, D) - F0[below]
    t[np.isclose(d, tu - ts, atol=1e-12)] = np.inf  # exactly on the unstable point
    return t

def lock_edges(K, D, eps, n):
    """Unwrapped edges (lower, upper) around θ* of the set of θ₀ locking within n rotations.
    The locked set is the arc [lower, upper]; P_n = (upper − lower)/2π. nan if no lock."""
    if K <= 0 or abs(D) >= K:
        return (-eps, eps) if (K == 0 and D == 0) else (np.nan, np.nan)
    ts, tu = fixed_points(K, D)
    if eps >= tu - ts:
        return np.nan, np.nan
    Tn = n * T_ROT
    if D == 0:                                       # exact inversion
        up = 2 * np.arctan(np.tan(eps / 2) * np.exp(K * Tn))
        return -up, up
    Fa, Fb = F(ts + eps, K, D), F(ts - eps, K, D)
    tiny = 1e-9
    g_up = lambda th: (Fa - F(th, K, D)) - Tn        # increasing on (θ*+ε, θ_u)
    g_lo = lambda th: (Fb - F(th, K, D)) - Tn        # decreasing on (θ_u − 2π, θ*−ε)
    upper = brentq(g_up, ts + eps + tiny, tu - tiny, xtol=1e-13, maxiter=500)
    lower = brentq(g_lo, tu - 2 * np.pi + tiny, ts - eps - tiny, xtol=1e-13, maxiter=500)
    return lower, upper

def likelihood(K, D, eps, n):
    lo, up = lock_edges(K, D, eps, n)
    if np.isnan(lo):
        return 0.0
    return (up - lo) / (2 * np.pi)

def x_needed(P, eps, n):
    """δ = 0 closed form: x with P_n(x) = P."""
    return np.log(np.tan(np.pi * P / 2) / np.tan(eps / 2)) / (2 * np.pi * n)

# ----------------------------------------------------------------------------- RK45 checks
def rk_time_to_band(theta0, K, D, eps, tmax):
    ts, _ = fixed_points(K, D)
    if abs(wrap(theta0 - ts)) <= eps:
        return 0.0
    ev = lambda t, y: abs(wrap(y[0] - ts)) - eps
    ev.terminal, ev.direction = True, -1
    sol = solve_ivp(lambda t, y: [D - K * np.sin(y[0])], (0, tmax), [theta0],
                    events=ev, rtol=1e-11, atol=1e-13, max_step=0.02)
    return sol.t_events[0][0] if len(sol.t_events[0]) else np.inf

def rk_mutual_time_to_band(theta0, K1, K2, w1, w2, eps, tmax):
    """Two clocks pulling each other; time for φ₁ − φ₂ to enter the band around θ*(K₁+K₂, ω₁−ω₂)."""
    ts, _ = fixed_points(K1 + K2, w1 - w2)
    rhs = lambda t, y: [w1 + K1 * np.sin(y[1] - y[0]), w2 + K2 * np.sin(y[0] - y[1])]
    ev = lambda t, y: abs(wrap(y[0] - y[1] - ts)) - eps
    ev.terminal, ev.direction = True, -1
    sol = solve_ivp(rhs, (0, tmax), [theta0, 0.0], events=ev, rtol=1e-11, atol=1e-13, max_step=0.02)
    return sol.t_events[0][0] if len(sol.t_events[0]) else np.inf

def run_verification():
    rng = np.random.default_rng(20260906)
    lines = ["Verification of the closed-form time to the band (PREREGISTRATION.md §3)", ""]
    worst = 0.0; n_cmp = 0
    for D in (0.0, 0.25):
        for _ in range(40):
            x = rng.uniform(max(0.05, D + 0.03), 1.5)
            th0 = rng.uniform(-np.pi, np.pi)
            t_cf = float(time_to_band(th0, x, D, EPS_PRIMARY))
            if not np.isfinite(t_cf) or t_cf > 20 * T_ROT:
                continue
            t_rk = rk_time_to_band(th0, x, D, EPS_PRIMARY, 25 * T_ROT)
            err = abs(t_rk - t_cf) / T_ROT
            worst = max(worst, err); n_cmp += 1
    lines.append(f"RK45 vs closed form: {n_cmp} samples (δ = 0 and 0.25, x in [0.05, 1.5], θ₀ uniform),")
    lines.append(f"  worst |t_RK45 − t_closed| = {worst:.2e} rotations   (tolerance 1e-6)   -> {'PASS' if worst < 1e-6 else 'FAIL'}")
    # edges against a grid measure
    worstP = 0.0
    th_grid = np.linspace(-np.pi, np.pi, 720001)[:-1]
    for D in DELTA_LIST:
        for x in (0.3, 0.7, 1.2):
            if x <= D: continue
            for n in N_LIST:
                t = time_to_band(th_grid, x, D, EPS_PRIMARY)
                P_grid = np.mean(t <= n * T_ROT)
                P_cf = likelihood(x, D, EPS_PRIMARY, n)
                worstP = max(worstP, abs(P_grid - P_cf))
    lines.append(f"Edge inversion vs grid measure of the locked set (720000 angles): worst |ΔP| = {worstP:.2e}   -> {'PASS' if worstP < 1e-5 else 'FAIL'}")
    # mutual two-clock check
    K1, K2, dw, th0 = 0.3, 0.2, 0.1, 2.5
    t_mut = rk_mutual_time_to_band(th0, K1, K2, 1.0, 1.0 - dw, EPS_PRIMARY, 25 * T_ROT)
    t_eff = float(time_to_band(th0, K1 + K2, dw, EPS_PRIMARY))
    err_m = abs(t_mut - t_eff) / T_ROT
    lines.append(f"Mutual clocks K₁ = {K1}, K₂ = {K2}, Δω = {dw}, θ₀ = {th0}: RK45 {t_mut/T_ROT:.6f} rot,"
                 f" Adler with K_eff = K₁+K₂ {t_eff/T_ROT:.6f} rot, |diff| = {err_m:.2e} rot -> {'PASS' if err_m < 1e-6 else 'FAIL'}")
    # the predicted table of §3
    lines.append("")
    lines.append("Pre-registered numbers (§3) recomputed by the script:")
    for n in N_LIST:
        lines.append(f"  n = {n}: x for P = 50/90/95/99 % at ε = 0.35: " +
                     ", ".join(f"{x_needed(P, EPS_PRIMARY, n):.3f}" for P in P_TARGETS))
    (OUT / "verification.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    return worst < 1e-6 and worstP < 1e-5 and err_m < 1e-6

# ----------------------------------------------------------------------------- tables
def write_tables():
    with open(OUT / "boundary_angles.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["eps_rad", "n_rotations", "x_K_over_omega", "K_times_T", "theta0_max_deg", "P_n", "P_floor_eps_over_pi"])
        for eps in EPS_LIST:
            for n in N_LIST:
                for x in X_TABLE:
                    lo, up = lock_edges(x, 0.0, eps, n)
                    w.writerow([f"{eps:.4f}", n, f"{x:.3f}", f"{2*np.pi*x:.3f}", f"{np.degrees(up):.2f}",
                                f"{(up-lo)/(2*np.pi):.4f}", f"{eps/np.pi:.4f}"])
    with open(OUT / "thresholds.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["eps_rad", "n_rotations", "target_P", "x_needed", "K_times_T_needed"])
        for eps in EPS_LIST:
            for n in N_LIST:
                for P in P_TARGETS:
                    x = x_needed(P, eps, n)
                    w.writerow([f"{eps:.4f}", n, f"{P:.2f}", f"{x:.4f}", f"{2*np.pi*x:.3f}"])
    with open(OUT / "detuned_edges.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["delta", "eps_rad", "n_rotations", "x_K_over_omega", "theta_star_deg", "lower_edge_deg", "upper_edge_deg", "P_n"])
        xs = np.round(np.arange(0.0, 1.5001, 0.01), 3)
        for D in DELTA_LIST:
            for n in N_LIST:
                for x in xs:
                    lo, up = lock_edges(float(x), D, EPS_PRIMARY, n)
                    if np.isnan(lo):
                        w.writerow([D, EPS_PRIMARY, n, f"{x:.2f}", "", "", "", "0.0000"]); continue
                    ts = 0.0 if x == 0 else np.degrees(np.arcsin(D / x))
                    w.writerow([D, EPS_PRIMARY, n, f"{x:.2f}", f"{ts:.2f}", f"{np.degrees(wrap(lo)):.2f}",
                                f"{np.degrees(up):.2f}", f"{(up-lo)/(2*np.pi):.4f}"])
    xs = np.round(np.arange(0.0, 1.5001, 0.01), 3)
    th = np.radians(np.arange(-180, 181, 1.0))
    with open(OUT / "lock_map_delta0.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["x_K_over_omega", "theta0_deg", "rotations_to_lock_eps0.35"])
        for x in xs:
            t = time_to_band(th, float(x), 0.0, EPS_PRIMARY) / T_ROT
            for a, r in zip(np.degrees(th), t):
                w.writerow([f"{x:.2f}", f"{a:.0f}", "inf" if not np.isfinite(r) else f"{r:.4f}"])

# ----------------------------------------------------------------------------- figures
SURFACE, TEXT, TEXT2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e6e5e1"
SLOT = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]
SEQ = ["#cde2fb", "#b7d3f6", "#9ec5f4", "#86b6ef", "#6da7ec", "#5598e7", "#3987e5",
       "#2a78d6", "#256abf", "#1c5cab", "#184f95", "#104281", "#0d366b"]
plt.rcParams.update({
    "font.family": "sans-serif", "font.size": 9, "axes.titlesize": 10, "axes.titleweight": "bold",
    "axes.edgecolor": GRID, "axes.labelcolor": TEXT2, "xtick.color": TEXT2, "ytick.color": TEXT2,
    "axes.titlecolor": TEXT, "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
    "savefig.facecolor": SURFACE, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.8,
    "grid.linestyle": "-", "axes.spines.top": False, "axes.spines.right": False,
    "legend.frameon": False, "legend.fontsize": 8, "axes.axisbelow": True,
})

def fig_lock_map():
    xs = np.linspace(0, 1.5, 301)
    th = np.radians(np.linspace(-180, 180, 721))
    R = np.empty((th.size, xs.size))
    for j, x in enumerate(xs):
        R[:, j] = time_to_band(th, float(x), 0.0, EPS_PRIMARY) / T_ROT
    Rc = np.where(np.isfinite(R), np.minimum(R, 4.0), 4.0)
    cmap = LinearSegmentedColormap.from_list("seqblue", SEQ)
    fig, ax = plt.subplots(figsize=(8.6, 4.8))
    im = ax.imshow(Rc, origin="lower", aspect="auto", extent=[0, 1.5, -180, 180], cmap=cmap, vmin=0, vmax=4)
    for xr in (1 / 64, 1 / 8):
        ax.axvline(xr, color=SURFACE, lw=0.9, alpha=0.9)
    ax.text(1 / 8 + 0.012, 8, "DK record, Stage 2: ω/K = 64 to 8", color=TEXT, fontsize=7.5, va="bottom")
    for n, c in zip(N_LIST, SLOT[1:3]):
        up = np.array([lock_edges(float(x), 0.0, EPS_PRIMARY, n)[1] for x in xs])
        for sgn in (1, -1):
            ax.plot(xs, sgn * np.degrees(up), color=c, lw=1.6, solid_capstyle="round",
                    label=f"locks within {n} rotation{'s' if n > 1 else ''}" if sgn == 1 else None)
        k = min(np.searchsorted(np.degrees(up), 120), xs.size - 1)
        sgn_lab = 1 if n == 1 else -1
        ax.annotate(f"within {n} rotation{'s' if n > 1 else ''}", (xs[k], sgn_lab * np.degrees(up[k])),
                    xytext=(7, -4), textcoords="offset points", va="center", color=TEXT, fontsize=8)
    ax.set_xlabel("x = K / ω  (pull per radian of the clock's turn;  K·T = 2πx)")
    ax.set_ylabel("starting phase difference θ₀ (degrees)")
    ax.set_yticks([-180, -90, 0, 90, 180])
    ax.set_title("Rotations to lock, Δ = 0, ε = 0.35 rad")
    ax.legend(loc="center right")
    cb = fig.colorbar(im, ax=ax, pad=0.02)
    cb.set_label("rotations to lock (4 and above shown darkest)", color=TEXT2)
    cb.ax.tick_params(colors=TEXT2); cb.outline.set_edgecolor(GRID)
    fig.tight_layout()
    fig.savefig(OUT / "fig1_lock_map.png", dpi=150)
    plt.close(fig)

def fig_likelihood():
    xs = np.linspace(0, 1.5, 301)
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.0), sharey=True)
    for ax, n in zip(axes, N_LIST):
        for eps, c in zip([0.35, 0.10, np.pi / 4], SLOT[:3]):
            P = np.array([likelihood(float(x), 0.0, eps, n) for x in xs])
            lab = f"ε = {EPS_LABEL[eps]} rad" + (" (primary)" if eps == 0.35 else "")
            ax.plot(xs, P, color=c, lw=1.6, solid_capstyle="round", label=lab)
            k = np.searchsorted(P, 0.40 if eps == np.pi / 4 else (0.62 if eps == 0.35 else 0.75))
            k = min(k, xs.size - 1)
            ax.annotate(f"ε = {EPS_LABEL[eps]}", (xs[k], P[k]), xytext=(5, -9), textcoords="offset points", color=TEXT, fontsize=8)
        for lev in (0.5, 0.95):
            ax.axhline(lev, color=TEXT2, lw=0.8, alpha=0.5)
            ax.text(1.49, lev + 0.012, f"P = {lev:.2f}", ha="right", color=TEXT2, fontsize=7.5)
        ax.set_ylim(0, 1.02); ax.set_xlim(0, 1.5)
        ax.set_xlabel("x = K / ω")
        ax.set_title(f"P{'₁' if n == 1 else '₂'}: lock within {n} rotation{'s' if n > 1 else ''}, Δ = 0")
    axes[0].set_ylabel("likelihood over a uniform starting angle")
    axes[0].legend(loc="lower right")
    fig.tight_layout()
    fig.savefig(OUT / "fig2_likelihood.png", dpi=150)
    plt.close(fig)

def fig_detuned():
    xs = np.linspace(0, 1.5, 301)
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.0), sharey=True)
    for ax, n in zip(axes, N_LIST):
        for D, c in zip(DELTA_LIST, SLOT):
            P = np.array([likelihood(float(x), D, EPS_PRIMARY, n) for x in xs])
            ax.plot(xs, P, color=c, lw=1.6, solid_capstyle="round", label=f"δ = Δ/ω = {D}")
            k = min(np.searchsorted(P, 0.55 + 0.1 * DELTA_LIST.index(D)), xs.size - 1)
            ax.annotate(f"δ = {D}", (xs[k], P[k]), xytext=(5, -9), textcoords="offset points", color=TEXT, fontsize=8)
        ax.set_ylim(0, 1.02); ax.set_xlim(0, 1.5)
        ax.set_xlabel("x = K / ω")
        ax.set_title(f"P{'₁' if n == 1 else '₂'} with detuning, ε = 0.35 rad (no lock for x ≤ δ)")
    axes[0].set_ylabel("likelihood over a uniform starting angle")
    axes[0].legend(loc="lower right")
    fig.tight_layout()
    fig.savefig(OUT / "fig3_detuned.png", dpi=150)
    plt.close(fig)

# ----------------------------------------------------------------------------- main
def main():
    ok = run_verification()
    write_tables()
    fig_lock_map(); fig_likelihood(); fig_detuned()
    print()
    print("Boundary angles and likelihood at ε = 0.35 rad, Δ = 0 (floor P = %.3f):" % (EPS_PRIMARY / np.pi))
    print("   x=K/ω   K·T   θ₀max(n=1)  P₁     θ₀max(n=2)  P₂")
    for x in X_TABLE:
        r = []
        for n in N_LIST:
            lo, up = lock_edges(x, 0.0, EPS_PRIMARY, n)
            r += [np.degrees(up), (up - lo) / (2 * np.pi)]
        print(f"  {x:5.2f}  {2*np.pi*x:5.2f}   {r[0]:6.1f}°   {r[1]:.3f}    {r[2]:6.1f}°   {r[3]:.3f}")
    info = [f"run: {datetime.datetime.now().isoformat(timespec='seconds')}",
            f"python {platform.python_version()}, numpy {np.__version__}",
            f"origin: DiracKuramotoFramework c9bbf12, tests/kuramoto_sync.py (equation and integration pattern)",
            f"verification: {'PASS' if ok else 'FAIL'}"]
    (OUT / "run_info.txt").write_text("\n".join(info) + "\n")
    print("\n".join(info))
    return 0 if ok else 1

if __name__ == "__main__":
    sys.exit(main())
