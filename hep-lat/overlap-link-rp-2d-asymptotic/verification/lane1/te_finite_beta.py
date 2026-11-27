import os  # release copy
"""T-e of projects/bounties/chiral-rp/CATCH3.md (plan te-plan.md, sealed sha 3ea0fe4fc38be925, 09:37:34, 2026-11-26).

Finite-beta normalisation at fixed physical size: beta = (L/(e l))^2. Exact Wilson-measure samples on the L x L torus
(plaquette angles von Mises with the torus constraint by rejection; links from a constant-flux field + Coulomb-gauge
Fourier solve + Haar holonomies). Overlap D_ov with the H_W kernel, m = 0.7, antiperiodic both ways (T-b conventions).
Delta(L) = <ln det D_ov | Q=0> - <ln det D_ov[flat]>_Haar.

Memory: L = 24 is a 1152 x 1152 complex matrix (21 MB); eigh and D_ov, under 200 MB.
"""
import sys, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # release copy: was an absolute path
from tc_holonomy import lndet_mom  # noqa: E402

s1 = np.array([[0, 1], [1, 0]], complex); s2m = np.array([[0, -1j], [1j, 0]], complex)
s3 = np.array([[1, 0], [0, -1]], complex); I2 = np.eye(2)
g = [s1, s2m]; m = 0.7

def build_DW_links(L, U0, U1):
    N = L * L; D = np.zeros((2 * N, 2 * N), complex)
    idx = lambda t, x: (t % L) * L + (x % L)
    for t in range(L):
        for x in range(L):
            i = idx(t, x)
            D[2*i:2*i+2, 2*i:2*i+2] += (2 - m) * I2
            for mu, (dt, dx) in enumerate([(1, 0), (0, 1)]):
                tt, xx = t + dt, x + dx
                u = (U0 if mu == 0 else U1)[t, x]
                if (mu == 0 and tt == L) or (mu == 1 and xx == L):
                    u = -u
                j = idx(tt, xx)
                D[2*i:2*i+2, 2*j:2*j+2] += -0.5 * (I2 - g[mu]) * u
                D[2*j:2*j+2, 2*i:2*i+2] += -0.5 * (I2 + g[mu]) * np.conj(u)
    return D

def lndet_ov(L, U0, U1):
    N = L * L; G5 = np.kron(np.eye(N), s3)
    H = G5 @ build_DW_links(L, U0, U1)
    w, V = np.linalg.eigh(H)
    Dov = np.eye(2 * N) + G5 @ ((V * np.sign(w)) @ V.conj().T)
    s, l = np.linalg.slogdet(Dov)
    return l, s.real, np.abs(w).min()

def plaq_angles(U0, U1):
    Up = U0 * np.roll(U1, -1, axis=0) * np.conj(np.roll(U0, -1, axis=1)) * np.conj(U1)
    return np.angle(Up)

def wrap(a):
    return (a + np.pi) % (2 * np.pi) - np.pi

def sample_plaquettes(L, beta, rng):
    V = L * L; tries = 0
    while True:
        tries += 1
        th = rng.vonmises(0.0, beta, V - 1)
        last = wrap(-th.sum())
        if rng.random() < np.exp(beta * (np.cos(last) - 1)):
            th = np.append(th, last).reshape(L, L)
            return th, int(round(th.sum() / (2 * np.pi))), tries

def links_from_plaquettes(th, Q, phi):
    L = th.shape[0]; V = L * L; c = 2 * np.pi * Q / V
    A0 = np.zeros((L, L)); A1 = np.zeros((L, L))
    A1 += c * np.arange(L)[:, None]
    A0[L - 1, :] = -c * L * np.arange(L)
    thp = th - c
    F = np.fft.fft2(thp)
    kt = 2 * np.pi * np.fft.fftfreq(L)[:, None]; kx = 2 * np.pi * np.fft.fftfreq(L)[None, :]
    ax = 1 - np.exp(1j * kx); at = np.exp(1j * kt) - 1
    Dn = np.abs(ax) ** 2 + np.abs(at) ** 2; Dn[0, 0] = 1.0
    B0 = F * np.conj(ax) / Dn; B1 = F * np.conj(at) / Dn; B0[0, 0] = 0; B1[0, 0] = 0
    A0 += np.real(np.fft.ifft2(B0)); A1 += np.real(np.fft.ifft2(B1))
    A0 += phi[0] / L; A1 += phi[1] / L
    return np.exp(1j * A0), np.exp(1j * A1)

def flat_lnavg(L, n):
    ph = 2 * np.pi * (np.arange(n) + 0.5) / n
    return np.mean([lndet_mom(L, (a, b)) for a in ph for b in ph])

if __name__ == "__main__":
    rng = np.random.default_rng(1137)
    print("T-e: finite beta at fixed e l; m = 0.7, antiperiodic torus, H_W overlap")
    # C1: flat, zero holonomy, L = 8 against the T-c momentum formula
    L = 8; U0, U1 = links_from_plaquettes(np.zeros((L, L)), 0, (0.0, 0.0))
    l, s, _ = lndet_ov(L, U0, U1); lm = lndet_mom(L, (0.0, 0.0))
    print(f"C1 {'PASS' if abs(l - lm) < 1e-9 else 'FAIL'}: dense {l:.10f} vs momentum {lm:.10f}")
    # C0: plaquette reconstruction, including Q != 0
    c0 = 0.0
    for (L, beta) in ((8, 1.0), (12, 0.5), (16, 4.0)):
        for _ in range(5):
            th, Q, _ = sample_plaquettes(L, beta, rng)
            U0, U1 = links_from_plaquettes(th, Q, rng.uniform(0, 2 * np.pi, 2))
            c0 = max(c0, np.abs(wrap(plaq_angles(U0, U1) - th)).max())
    print(f"C0 {'PASS' if c0 < 1e-10 else 'FAIL'}: max plaquette reconstruction error {c0:.1e}")
    if c0 >= 1e-10 or abs(l - lm) >= 1e-9:
        sys.exit(1)
    # one Q != 0 control: det ~ 0
    L = 8; th, Q, _ = sample_plaquettes(L, 1.0, rng)
    while Q == 0:
        th, Q, _ = sample_plaquettes(L, 1.0, rng)
    U0, U1 = links_from_plaquettes(th, Q, rng.uniform(0, 2 * np.pi, 2))
    lq, sq, gq = lndet_ov(L, U0, U1)
    print(f"Q != 0 control: L = 8, Q = {Q}, ln|det| = {lq:.3f} (flat ~ 69.4), min|eig H_W| = {gq:.3e}")

    Lgrid = [8, 12, 16, 20, 24]; res = {}
    print("\n el   L   beta    acc/tries   P(Q=0)       <lndet|Q=0> - flat = Delta  +- se    flat(n=64)    flat n64-n128   min gap")
    for el in (2.0, 4.0):
        for L in Lgrid:
            beta = (L / el) ** 2; t0 = time.time()
            fl64 = flat_lnavg(L, 64); fl128 = flat_lnavg(L, 128)
            lds = []; nq0 = 0; tries_tot = 0; gaps = []
            for _ in range(40):
                th, Q, tr = sample_plaquettes(L, beta, rng); tries_tot += tr
                if Q != 0:
                    continue
                nq0 += 1
                U0, U1 = links_from_plaquettes(th, Q, rng.uniform(0, 2 * np.pi, 2))
                l, s, gp = lndet_ov(L, U0, U1); lds.append(l); gaps.append(gp)
                assert s > 0
            lds = np.array(lds); p0 = nq0 / 40; sep = np.sqrt(p0 * (1 - p0) / 40)
            D = lds.mean() - fl64; se = lds.std(ddof=1) / np.sqrt(len(lds))
            res[(el, L)] = (D, se, p0, sep)
            print(f"{el:3.0f}  {L:3d}  {beta:6.2f}   40/{tries_tot:<6d}  {p0:.3f}+-{sep:.3f}   {D:+.4f}  +- {se:.4f}   {fl64:12.6f}   {fl64-fl128:+.1e}   {min(gaps):.2e}   ({time.time()-t0:.0f}s)", flush=True)
    print("\nslopes of Delta on V for L >= 12:")
    for el in (2.0, 4.0):
        Vs = [L * L for L in Lgrid if L >= 12]; Ds = [res[(el, L)][0] for L in Lgrid if L >= 12]
        sl = np.polyfit(Vs, Ds, 1)[0]
        d24, d16 = res[(el, 24)], res[(el, 16)]; d12 = res[(el, 12)]
        print(f"  el = {el:.0f}: slope {sl:+.2e} per site;  Delta(24) - Delta(16) = {d24[0]-d16[0]:+.4f} (3se = {3*np.hypot(d24[1], d16[1]):.4f});  Delta(24) - Delta(12) = {d24[0]-d12[0]:+.4f}")
    print(f"  ratio Delta(el=4)/Delta(el=2) at L = 24: {res[(4.0, 24)][0] / res[(2.0, 24)][0]:.3f}")
    for el in (2.0, 4.0):
        var = el ** 2 / (4 * np.pi ** 2); qs = np.arange(-10, 11)
        p = np.exp(-qs ** 2 / (2 * var)); print(f"  discrete-Gaussian P(Q=0) at el = {el:.0f}: {p[10] / p.sum():.3f}")
