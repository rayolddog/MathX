"""
ARNOLD TONGUE: do two CONFINED Dirac particles phase-lock their Zitterbewegung
when detuned (unequal masses)?  Genuine Kuramoto entrainment shows a wedge:
locked inside, drifting outside, with threshold coupling g* proportional to the
detuning  Dw = 2(m1-m2)c^2/hbar.

Same engine as dirac_bound_sync.py but with m1 != m2.  Units hbar=c=1, reference
mass = 1.  m1 = 1+d/2, m2 = 1-d/2  ->  Dw = 2d.  Scalar harmonic well binds both.

Order parameter: Kuramoto phasor  R = |<exp(i(phi1-phi2))>|_late  from the
analytic signal of the (high-passed) ZBW velocity <alpha_i>(t).
  R ~ 1 locked,  R ~ 0 drifting.
"""
import numpy as np

N  = 128
dx = 0.125
x  = (np.arange(N) - N//2)*dx
p  = 2*np.pi*np.fft.fftfreq(N, d=dx)
dt = 0.02
T  = 70.0
nt = int(T/dt)
ts = np.arange(nt)*dt
sigma = 1.0
w     = 1.5
kappa = 0.25

X1 = x[:, None]; X2 = x[None, :]
Ur = np.exp(-(X1 - X2)**2/(2*w**2))
Vc1 = 0.5*kappa*X1**2 + 0.0*X2
Vc2 = 0.5*kappa*X2**2 + 0.0*X1
d1W = Vc1 + Vc2
d2W = Vc1 - Vc2

def prop(m):
    E = np.sqrt(p**2 + m**2)
    ctt, stt = np.cos(E*dt), np.sin(E*dt)
    return ctt - 1j*stt*m/E, -1j*stt*p/E, -1j*stt*p/E, ctt + 1j*stt*m/E

def make_psi():
    g = (2*np.pi*sigma**2)**(-0.25)*np.exp(-x**2/(4*sigma**2))
    s1 = np.array([1, 1])/np.sqrt(2); s2 = np.array([1, 1j])/np.sqrt(2)
    return (g[:, None, None, None]*g[None, :, None, None]
            * s1[None, None, :, None]*s2[None, None, None, :]).astype(complex)

def kick_H1(psi, P):
    u11, u12, u21, u22 = P
    ph = np.fft.fft(psi, axis=0)
    u = ph[:, :, 0, :].copy(); d = ph[:, :, 1, :].copy()
    ph[:, :, 0, :] = u11[:, None, None]*u + u12[:, None, None]*d
    ph[:, :, 1, :] = u21[:, None, None]*u + u22[:, None, None]*d
    return np.fft.ifft(ph, axis=0)

def kick_H2(psi, P):
    u11, u12, u21, u22 = P
    ph = np.fft.fft(psi, axis=1)
    u = ph[:, :, :, 0].copy(); d = ph[:, :, :, 1].copy()
    ph[:, :, :, 0] = u11[None, :, None]*u + u12[None, :, None]*d
    ph[:, :, :, 1] = u21[None, :, None]*u + u22[None, :, None]*d
    return np.fft.ifft(ph, axis=1)

def apply_W(psi, c, tau):
    out = np.empty_like(psi)
    for (ia, ib, d) in [((0, 0), (1, 1), d1W), ((0, 1), (1, 0), d2W)]:
        Om = np.sqrt(d**2 + c**2) + 1e-30
        co = np.cos(Om*tau); si = np.sin(Om*tau)
        fd = si*d/Om; fc = si*c/Om
        a = psi[:, :, ia[0], ia[1]]; b = psi[:, :, ib[0], ib[1]]
        out[:, :, ia[0], ia[1]] = (co - 1j*fd)*a + (-1j*fc)*b
        out[:, :, ib[0], ib[1]] = (-1j*fc)*a + (co + 1j*fd)*b
    return out

def run(g, delta):
    P1 = prop(1.0 + delta/2); P2 = prop(1.0 - delta/2)
    c = g*Ur
    psi = make_psi()
    v1 = np.empty(nt); v2 = np.empty(nt)
    for n in range(nt):
        nrm = (np.abs(psi)**2).sum()*dx*dx
        v1[n] = (2*np.real(np.conj(psi[:, :, 0, :])*psi[:, :, 1, :])).sum()*dx*dx/nrm
        v2[n] = (2*np.real(np.conj(psi[:, :, :, 0])*psi[:, :, :, 1])).sum()*dx*dx/nrm
        psi = apply_W(psi, c, dt/2)
        psi = kick_H1(psi, P1); psi = kick_H2(psi, P2)
        psi = apply_W(psi, c, dt/2)
    return v1, v2, nrm

def highpass(s, k):
    return s - np.convolve(s, np.ones(k)/k, mode="same")

def analytic(s):
    n = len(s); S = np.fft.fft(s); h = np.zeros(n)
    h[0] = 1; h[1:n//2] = 2; h[n//2] = 1
    return np.fft.ifft(S*h)

def order_param(v1, v2):
    k = int(round(3.0/dt))
    z1 = analytic(highpass(v1, k)); z2 = analytic(highpass(v2, k))
    psi_rel = np.angle(z1) - np.angle(z2)
    lo, hi = int(0.5*nt), int(0.95*nt)
    return np.abs(np.mean(np.exp(1j*psi_rel[lo:hi])))

# ---------------- sweep grid ----------------
deltas = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5]
gs     = [0.0, 0.4, 0.8, 1.2, 1.6, 2.0]
SANITY = False        # set True for a quick 2x2 validation run

if SANITY:
    deltas = [0.0, 0.4]; gs = [0.0, 1.6]

R = np.zeros((len(gs), len(deltas)))
for ig, g in enumerate(gs):
    for jd, dlt in enumerate(deltas):
        v1, v2, nrm = run(g, dlt)
        R[ig, jd] = order_param(v1, v2)
        print(f"g={g:.2f} delta={dlt:.2f} (Dw={2*dlt:.2f}): R={R[ig,jd]:.3f}  norm={nrm:.3f}", flush=True)

np.savez("/Users/john-bramble/arnold_tongue.npz", R=R, deltas=np.array(deltas), gs=np.array(gs))
print("R matrix (rows=g, cols=delta):")
print(np.round(R, 3))

if not SANITY:
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(7.4, 5.6))
    ddl = deltas[1]-deltas[0]; dgg = gs[1]-gs[0]
    ext = [deltas[0]-ddl/2, deltas[-1]+ddl/2, gs[0]-dgg/2, gs[-1]+dgg/2]
    im = ax.imshow(R, origin="lower", aspect="auto", extent=ext,
                   cmap="magma", vmin=0, vmax=1)
    cs = ax.contour(deltas, gs, R, levels=[0.5], colors="cyan", linewidths=2)
    ax.clabel(cs, fmt="lock boundary R=0.5", fontsize=8)
    # Kuramoto prediction g* proportional to detuning: g* = slope * (2 delta)
    dd = np.array(deltas)
    ax.plot(dd, 2.2*2*dd, "w--", lw=1.6, label="Kuramoto guide  g* ∝ Δω")
    ax.set_xlabel("detuning  δ   (Δω = 2δ = ZBW frequency mismatch)")
    ax.set_ylabel("Breit coupling  g   (interaction energy / mc²)")
    ax.set_title("Arnold tongue: Zitterbewegung phase-locking of two confined Dirac particles\n"
                 "R = Kuramoto order parameter (1 = locked, 0 = drifting)")
    ax.legend(loc="upper right", fontsize=8, framealpha=0.3)
    fig.colorbar(im, ax=ax, label="order parameter R")
    fig.tight_layout()
    fig.savefig("/Users/john-bramble/arnold_tongue.png", dpi=140)
    print("saved arnold_tongue.png")
