"""T-c of projects/bounties/chiral-rp/CATCH3.md (plan tc-plan.md, sealed sha 1159a3d5f9c318b6, 09:18:28, 2026-11-26).

<det D_ov>_W over Haar holonomies on T-b's L x L antiperiodic torus with flat plaquettes (m = 0.7, m_f = 0).
Flat fields are gauge-equivalent to constant links U_mu = e^{i phi_mu/L}. In momentum space
k_mu = (2 pi n + pi + phi_mu)/L, D_W(k) = b + i s.gamma, and det D_ov(k) = 2 (1 + b/w), w = sqrt(s^2 + b^2).
The dense builder is T-b's (tb_flux_pair.py), with the Dirac string replaced by constant holonomy links; it is the C0 control.

Memory: dense control at L = 12 is 288 x 288 complex (1.3 MB). The momentum sums are L x L arrays per holonomy
(L = 40: 1600 doubles). f_inf is done in rows of 8192 (64 KB each). Peak well under 100 MB.
"""
import numpy as np, sys, time
s1 = np.array([[0, 1], [1, 0]], complex); s2m = np.array([[0, -1j], [1j, 0]], complex)
s3 = np.array([[1, 0], [0, -1]], complex); I2 = np.eye(2)
g = [s1, s2m]; m = 0.7

def build_DW_hol(L, phi):
    N = L * L; D = np.zeros((2 * N, 2 * N), complex)
    idx = lambda t, x: (t % L) * L + (x % L)
    for t in range(L):
        for x in range(L):
            i = idx(t, x)
            D[2*i:2*i+2, 2*i:2*i+2] += (2 - m) * I2
            for mu, (dt, dx) in enumerate([(1, 0), (0, 1)]):
                tt, xx = t + dt, x + dx
                u = np.exp(1j * phi[mu] / L)
                if (mu == 0 and tt == L) or (mu == 1 and xx == L):
                    u *= -1.0  # antiperiodic fermion boundary
                j = idx(tt, xx)
                D[2*i:2*i+2, 2*j:2*j+2] += -0.5 * (I2 - g[mu]) * u
                D[2*j:2*j+2, 2*i:2*i+2] += -0.5 * (I2 + g[mu]) * np.conj(u)
    return D

def lndet_dense(L, phi):
    N = L * L; G5 = np.kron(np.eye(N), s3)
    H = G5 @ build_DW_hol(L, phi)
    assert np.abs(H - H.conj().T).max() < 1e-12
    w, V = np.linalg.eigh(H)
    Dov = np.eye(2 * N) + G5 @ ((V * np.sign(w)) @ V.conj().T)
    s, l = np.linalg.slogdet(Dov)
    return l, s

def lnfac(k0, k1):
    """ln[2(1 + b/w)] elementwise, cancellation-safe for b < 0."""
    S2 = np.sin(k0) ** 2 + np.sin(k1) ** 2
    b = (1 - np.cos(k0)) + (1 - np.cos(k1)) - m
    w = np.sqrt(S2 + b * b)
    with np.errstate(divide='ignore'):
        one_plus = np.where(b >= 0, (w + b) / w, S2 / (w * (w - b)))
        return np.log(2 * one_plus)

def lndet_mom(L, phi):
    n = np.arange(L)
    k0 = (2 * np.pi * n + np.pi + phi[0]) / L; k1 = (2 * np.pi * n + np.pi + phi[1]) / L
    return lnfac(k0[:, None], k1[None, :]).sum()

def lse(a):
    a = np.asarray(a); mx = a[np.isfinite(a)].max()
    return mx + np.log(np.exp(a - mx).sum())

def f_inf(N, off):
    tot = 0.0; k1 = 2 * np.pi * (np.arange(N) + off[1]) / N
    for j in range(N):
        tot += lnfac(np.full(N, 2 * np.pi * (j + off[0]) / N), k1).sum()
    return tot / N ** 2

if __name__ == "__main__":
    rng = np.random.default_rng(20261126)
    print("T-c: flat plaquettes, Haar holonomies, m = 0.7, antiperiodic torus")
    print("C0: momentum formula vs dense slogdet")
    c0max = 0.0
    for L in (8, 12):
        for _ in range(3):
            phi = rng.uniform(0, 2 * np.pi, 2)
            ld, sd = lndet_dense(L, phi); lm = lndet_mom(L, phi)
            c0max = max(c0max, abs(ld - lm))
            print(f"  L={L:2d} phi=({phi[0]:.3f},{phi[1]:.3f})  dense {ld:.10f} sign {sd.real:+.0f}  mom {lm:.10f}  diff {ld-lm:.2e}")
    print(f"C0 {'PASS' if c0max < 1e-8 else 'FAIL'}: max diff {c0max:.2e}")
    if c0max >= 1e-8:
        sys.exit(1)
    t0 = time.time()
    fa = f_inf(8192, (0.5, 0.5)); fb = f_inf(8192, (0.25, 0.75)); fc = f_inf(4096, (0.25, 0.75))
    print(f"f_inf: N=8192 off(.5,.5) {fa:.10f}  off(.25,.75) {fb:.10f}  N=4096 off(.25,.75) {fc:.10f}  ({time.time()-t0:.0f}s)")
    finf = fb
    print(f"  using f_inf = {finf:.10f}; c_f = {-finf:.10f}; worst-case loss per site 2ln2 - f_inf = {2*np.log(2)-finf:.6f}")
    Ls = [int(a) for a in sys.argv[1:]] or [8, 12, 16, 24, 32, 40]
    lnavg = {}; ln0 = {}
    print("\n  L    ln det(0)        r0=ln det0 - V f      ln<det> (n=32)    ln<det> (n=64)    n-conv     R=ln<det>-ln det0   r=ln<det>-V f   MC R (4000) +- se")
    for L in Ls:
        V = L * L; ln0[L] = lndet_mom(L, (0.0, 0.0))
        vals = {}
        for nq in (32, 64):
            ph = 2 * np.pi * np.arange(nq) / nq
            lds = [lndet_mom(L, (a, b)) for a in ph for b in ph]
            vals[nq] = lse(lds) - np.log(nq * nq)
        lnavg[L] = vals[64]
        mc = np.array([lndet_mom(L, rng.uniform(0, 2 * np.pi, 2)) for _ in range(4000)]) - ln0[L]
        e = np.exp(mc); mR = np.log(e.mean()); se = e.std(ddof=1) / np.sqrt(len(e)) / e.mean()
        print(f"{L:3d}  {ln0[L]:15.8f}   {ln0[L]-V*finf:10.6f}       {vals[32]:15.8f}   {vals[64]:15.8f}   {vals[64]-vals[32]:.1e}     {vals[64]-ln0[L]:+.6f}          {vals[64]-V*finf:+.6f}     {mR:+.4f} +- {se:.4f}", flush=True)
    print("\nper-pair slopes sl(L1,L2) = (ln<det>(L2) - ln<det>(L1))/(V2 - V1), minus f_inf:")
    for a, b in zip(Ls[:-1], Ls[1:]):
        sl = (lnavg[b] - lnavg[a]) / (b * b - a * a)
        print(f"  ({a:2d},{b:2d}): slope {sl:.10f}   slope - f_inf {sl - finf:+.2e}")
    print("\nconstants: ln 2 =", np.log(2), "  ln[(1/sqrt2)/(sqrt2 eta(i))^2] =",
          np.log((1 / np.sqrt(2)) / (2 * (__import__('math').gamma(0.25) / (2 * np.pi ** 0.75)) ** 2)))
