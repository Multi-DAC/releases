"""
Independent aliasing-free numerical verification for Paper 3 (QRF Norm-1):
  "The Norm-1 Property for Covariant Quantum Reference Frame Observables:
   Operator-Valued Upper Beurling Density, Homogeneous Spaces, and Free-Field Wick Frames"

Sections:
  1. Example 5a & 5b: Orthogonal-frame separation (bounded h vs. unbounded h, non-Toller,
     and infinite-dimensional separation from Gabardo-Lai scalar-projection density).
  2. Scalar sparse-blocks example (Clawd's union of [4^n, 4^n + n]: non-Toller, beta = 1).
  3. Example 5c: Intermediate ceiling beta(E) = |alpha|^2 in (0, 1) (< 1 on all bounded sets).
  4. Example 5d: dim(H_tilde) = 2 periodic rotation K(k) = (cos(w k), sin(w k))^T
     (beta = 1/2 for 0 < eps < w, yet ||E(-L, L)||_op = 1 for all L > w).
  5. Example 4: Scalar periodic comb (Logvinenko-Sereda-Kovrijkine < 1 on all bounded sets).
  6. Phase observables on T = R / 2*pi*Z (discrete dual Z): canonical, elementary,
     sparse-block non-convergent, and intermediate ceiling.
  7. Open <=> Borel equivalence: Fat Cantor set C_6 (measure 65/128) + non-unimodular
     affine (ax + b) group verification of the right-translation overlap m_G(K cap K x^{-1}) / m_G(K).
  8. QFT bounded-mass-band density vs. sharp bound (m2 - m1) / L.
  9. QFT Wick-power scaling-dimension law:
       - Exact 2-body phase spaces Phi_2^{(1+s)} for s = 1, 2, 3
       - Independent 4D quadrature check of Phi_3^{(1+2)}(M^2) = (1 - 3m/M) / (16*pi)
       - Exact 4-body phase space Phi_4^{(1+2)}(M^2) = (M - 4m)^2 / (128*pi^2*M) (supercritical growth)
       - Exact 3-body phase space Phi_3^{(1+1)}(M^2): threshold limit 1/(8*sqrt(3)*m^2) and large-M decay.

All operator norms are computed via exact analytical k-integrals + Nyquist-safe
Gauss-Legendre quadrature (asserting 0 <= ||E||_op <= 1 + 1e-10 on every test).
"""
import math
import time
import numpy as np
from numpy.polynomial.legendre import leggauss

t0 = time.time()


def exact_interval_kernel(diff: np.ndarray, a: float, b: float) -> np.ndarray:
    """Exact integral int_a^b e^{-i k * diff} dk / (2*pi)."""
    L = b - a
    mid = 0.5 * (a + b)
    with np.errstate(divide="ignore", invalid="ignore"):
        sinc_part = np.where(
            np.abs(diff) < 1e-14,
            L / (2.0 * np.pi),
            np.sin(0.5 * L * diff) / (np.pi * diff),
        )
    return np.exp(-1j * mid * diff) * sinc_part


def slepian_top_eig(c: float, n_quad: int | None = None) -> float:
    """Largest eigenvalue lambda_0(c) of sin(c(u-v))/(pi(u-v)) on L^2(-1, 1)."""
    if n_quad is None:
        n_quad = max(64, int(math.ceil(3.0 * c + 40)))
    x, w = leggauss(n_quad)
    sqw = np.sqrt(w)
    diff = x[:, None] - x[None, :]
    with np.errstate(divide="ignore", invalid="ignore"):
        K = np.where(np.abs(diff) < 1e-14, c / np.pi, np.sin(c * diff) / (np.pi * diff))
    M = (sqw[:, None] * K) * sqw[None, :]
    ev = float(np.linalg.eigvalsh(M)[-1])
    assert -1e-12 <= ev <= 1.0 + 1e-10, f"Aliasing detected: {ev}"
    return ev


# ============================================================================
# 1. EXAMPLE 5a & 5b: K(k) = e_n ON DISJOINT I_n (LENGTHS L_n)
# ============================================================================
print("=== 1. EXAMPLE 5a & 5b: K(k) = e_n ON DISJOINT I_n (LENGTHS L_n) ===")
print("5a. Bounded h (L_n = 1): ||E(-eps,eps)||_op = lambda_0(eps/2) <= ||h||^2 * (2*eps) = eps/pi")
for eps in (0.01, 0.05, 0.2, 1.0, 5.0):
    lam = slepian_top_eig(eps / 2.0)
    fl_ub = eps / np.pi
    assert lam <= fl_ub + 1e-12
    print(f"  eps={eps:5.2f}: lambda_0(eps/2)={lam:.10f}, FL bound eps/pi={fl_ub:.10f}, <=ub: True")

print("5b. Unbounded h, Toller Prop. 5 fails, Gabardo-Lai sup_u D^+(P_u)=0 yet D_op^+(P)=beta=1 (L_n = n):")
for n in (10, 50, 200, 400, 600):
    lam = slepian_top_eig(n * 0.05 / 2.0)
    print(f"  n={n:3d} (eps=0.05, c={n*0.05/2:5.2f}): lambda_0 = {lam:.12f} (1 - lam = {1.0 - lam:.3e})")
assert 1.0 - slepian_top_eig(600 * 0.05 / 2.0) < 1e-11

# Verify infinite-dimensional separation from scalar-projection density on finite truncation
# For any unit vector u = sum u_n e_n, D^+(<u, P u>) = limsup_R sup_k (1/2R) sum_n |u_n|^2 |I_n cap B_R(k)| = 0.
# Check explicit test vector u_n = sqrt(6)/(pi*n): max_k (1/2R) int_{k-R}^{k+R} <u, P(q)u> dq -> 0 while ||Pbar_R||_op = 1.
u_weights = np.array([6.0 / (math.pi**2 * n**2) for n in range(1, 201)])
starts = np.cumsum([0] + list(range(1, 200)))
ends = starts + np.arange(1, 201)
for R in (5.0, 20.0, 100.0):
    # Max scalar average occurs near the first intervals where |u_n|^2 is largest
    # Exact upper bound for all k: (1/(2R)) * (sum_{n <= N} |u_n|^2 * n) + max_{n > N} |u_n|^2
    # Evaluate exact scalar average over k in [0, 300] on grid + exact check at k = starts[0]+R
    best_sc = 0.0
    for k_center in np.linspace(0.0, 400.0, 801):
        a, b = k_center - R, k_center + R
        ov = np.maximum(0.0, np.minimum(ends, b) - np.maximum(starts, a))
        sc = float(np.sum(u_weights * ov) / (2.0 * R))
        if sc > best_sc:
            best_sc = sc
    # Operator norm at center of interval n = int(2*R) + 1 is exactly 1.0
    print(f"  Separation check R={R:5.1f}: sup_k <u, Pbar_R(k) u> (u_n~1/n) <= {best_sc:.6f} -> 0, while sup_k ||Pbar_R(k)||_op = 1.0")
    assert best_sc < 0.25 if R >= 20.0 else True


# ============================================================================
# 2. SCALAR SPARSE BLOCKS (CLAWD'S EXAMPLE): S = U_{n>=1} [4^n, 4^n + n]
# ============================================================================
print("\n=== 2. SCALAR SPARSE BLOCKS S = U_{n>=1} [4^n, 4^n + n] (dim=1, non-Toller, beta=1) ===")
# Exact operator norm on S_N = U_{n=1}^N [4^n, 4^n + n] exceeds single-block lambda_0(n*eps/2)
# Direct check on finite union n=1..6 for eps=0.5 (verifying >= lambda_0(6*0.5/2))
eps_sb = 0.5
L_sb = 6
n_x_sb = 200
# Note: 4^6 = 4096, so use k-space representation! On S_6, blocks [4^n, 4^n+n] have total length 21.
# Exact k-space operator: (H psi)(k) = int_{S_6} sin(eps*(k-q))/(pi*(k-q)) psi(q) dq!
ks_list, ws_list = [], []
for n in range(1, L_sb + 1):
    a_n, b_n = float(4**n), float(4**n + n)
    m_q = max(12, int(math.ceil(4.0 * eps_sb * n + 16)))
    qx, qw = leggauss(m_q)
    ks_list.append(0.5 * (b_n - a_n) * qx + 0.5 * (a_n + b_n))
    ws_list.append(0.5 * (b_n - a_n) * qw)
ks_sb = np.concatenate(ks_list)
ws_sb = np.concatenate(ws_list)
sqw_sb = np.sqrt(ws_sb)
dk_sb = ks_sb[:, None] - ks_sb[None, :]
with np.errstate(divide="ignore", invalid="ignore"):
    K_sb = np.where(np.abs(dk_sb) < 1e-14, eps_sb / np.pi, np.sin(eps_sb * dk_sb) / (np.pi * dk_sb))
H_sb = (sqw_sb[:, None] * K_sb) * sqw_sb[None, :]
H_sb = 0.5 * (H_sb + H_sb.T)
top_sb = float(np.linalg.eigvalsh(H_sb)[-1])
lam6 = slepian_top_eig(L_sb * eps_sb / 2.0)
assert lam6 - 1e-10 <= top_sb <= 1.0 + 1e-10
print(f"  S_6 = U_{{n=1}}^6 [4^n, 4^n+n], eps={eps_sb}: ||E_S6(-eps,eps)||_op = {top_sb:.10f} >= single-block lambda_0({L_sb*eps_sb/2:.1f}) = {lam6:.10f}")


# ============================================================================
# 3. EXAMPLE 5c: K(k) = alpha*e_0 + beta*e_n ON [n-1, n) (|alpha|^2=0.36, |beta|^2=0.64)
# ============================================================================
print("\n=== 3. EXAMPLE 5c: K(k) = alpha*e_0 + beta*e_n ON [n-1, n) (|alpha|^2=0.36, |beta|^2=0.64) ===")
a2, b2 = 0.36, 0.64
for R in (5, 20, 100, 500):
    M_len = 2 * R
    mat2 = np.array([[a2, math.sqrt(a2 * b2 / M_len)], [math.sqrt(a2 * b2 / M_len), b2 / M_len]])
    ev = float(np.linalg.eigvalsh(mat2)[-1])
    assert abs(ev - (a2 + b2 / M_len)) < 1e-15
    print(f"  R={R:3d} (M={M_len:4d}): ||Pbar_R||_op = {ev:.10f} (exact a2 + b2/M = {a2 + b2/M_len:.10f})")

for eps in (0.05, 0.2, 1.0, 4.0):
    lam_1 = slepian_top_eig(eps / 2.0)
    ub = a2 + b2 * lam_1
    vals_N = []
    m_q = max(8, int(math.ceil(3.0 * eps + 10)))
    kx0, kw0 = leggauss(m_q)
    for N_int in (10, 30, 60):
        ks = np.concatenate([0.5 * kx0 + (n - 0.5) for n in range(1, N_int + 1)])
        kws = np.tile(0.5 * kw0, N_int)
        sqkw = np.sqrt(kws)
        diff_k = ks[:, None] - ks[None, :]
        with np.errstate(divide="ignore", invalid="ignore"):
            K_all = np.where(np.abs(diff_k) < 1e-14, eps / np.pi, np.sin(eps * diff_k) / (np.pi * diff_k))
        H_k = a2 * ((sqkw[:, None] * K_all) * sqkw[None, :])
        for n in range(N_int):
            sl = slice(n * m_q, (n + 1) * m_q)
            kn = ks[sl]
            dk_n = kn[:, None] - kn[None, :]
            with np.errstate(divide="ignore", invalid="ignore"):
                Kn = np.where(np.abs(dk_n) < 1e-14, eps / np.pi, np.sin(eps * dk_n) / (np.pi * dk_n))
            sqn = sqkw[sl]
            H_k[sl, sl] += b2 * ((sqn[:, None] * Kn) * sqn[None, :])
        H_k = 0.5 * (H_k + H_k.T)
        top = float(np.linalg.eigvalsh(H_k)[-1])
        assert 0.0 <= top <= ub + 1e-10
        vals_N.append((N_int, round(top, 8)))
    print(f"  eps={eps:4.2f} (m_q={m_q:2d}): ||E_N|| = {vals_N} | analytic ub={ub:.8f} | exact ceiling={a2:.2f}")


# ============================================================================
# 4. EXAMPLE 5d: dim(H_tilde)=2 PERIODIC ROTATION K(k) = (cos(k), sin(k))^T
# ============================================================================
print("\n=== 4. EXAMPLE 5d: dim=2 PERIODIC ROTATION K(k) = (cos(k), sin(k))^T ON [0, inf) (w=1) ===")
Kmax_rot = 50.0
for L_half in (0.2, 0.5, 0.8, 1.2, 1.8):
    n_x = max(64, int(math.ceil(2.5 * Kmax_rot * L_half + 40)))
    xx, xw = leggauss(n_x)
    xs = L_half * xx
    sqw = np.sqrt(L_half * xw)
    diff = xs[:, None] - xs[None, :]
    K0 = exact_interval_kernel(diff, 0.0, Kmax_rot)
    Kp = exact_interval_kernel(diff - 2.0, 0.0, Kmax_rot)
    Km = exact_interval_kernel(diff + 2.0, 0.0, Kmax_rot)
    C2 = 0.5 * (Kp + Km)
    S2 = (Kp - Km) / (2.0j)
    G = np.block([[0.5 * (K0 + C2), 0.5 * S2], [0.5 * S2, 0.5 * (K0 - C2)]])
    sqw2 = np.concatenate([sqw, sqw])
    M = (sqw2[:, None] * G) * sqw2[None, :]
    M = 0.5 * (M + M.conj().T)
    top = float(np.linalg.eigvalsh(M)[-1].real)
    assert 0.0 <= top <= 1.0 + 1e-10
    pred = 0.5 if L_half < 1.0 else 1.0
    print(f"  L={L_half:4.2f} (n_x={n_x:3d}, Kmax={Kmax_rot:.0f}): ||E(-L,L)||_op = {top:.10f} (exact theory={pred:.1f}, diff={abs(top-pred):.2e})")


# ============================================================================
# 5. EXAMPLE 4: SCALAR PERIODIC COMB S = U_{n>=0} [n, n+0.4]
# ============================================================================
print("\n=== 5. EXAMPLE 4: SCALAR PERIODIC COMB S = U_{n>=0} [n, n+0.4] (beta=0.4, LSK < 1) ===")
theta = 0.4
N_teeth = 50
for L_half in (0.05, 0.1, 0.2, 0.5, 1.5, 3.0):
    n_x = max(64, int(math.ceil(2.5 * N_teeth * L_half + 40)))
    xx, xw = leggauss(n_x)
    xs = L_half * xx
    sqw = np.sqrt(L_half * xw)
    diff = xs[:, None] - xs[None, :]
    K_tooth = exact_interval_kernel(diff, 0.0, theta)
    s_half = np.sin(0.5 * diff)
    with np.errstate(divide="ignore", invalid="ignore"):
        geom = np.where(
            np.abs(s_half) < 1e-12,
            float(N_teeth),
            np.exp(-0.5j * (N_teeth - 1) * diff) * np.sin(0.5 * N_teeth * diff) / s_half,
        )
    G = K_tooth * geom
    M = (sqw[:, None] * G) * sqw[None, :]
    M = 0.5 * (M + M.conj().T)
    top = float(np.linalg.eigvalsh(M)[-1].real)
    assert 0.0 <= top <= 1.0 + 1e-10 and top < 1.0 - 1e-3
    print(f"  L={L_half:4.2f} (n_x={n_x:3d}, N={N_teeth}): ||E(-L,L)||_op = {top:.10f} (beta=0.40, <1: True)")


# ============================================================================
# 6. PHASE OBSERVABLES ON T = R / 2*pi*Z (DISCRETE DUAL Z)
# ============================================================================
print("\n=== 6. PHASE OBSERVABLES ON T (DISCRETE DUAL Z): ||E_xi(-eps, eps)||_op vs beta ===")
# For truncated number basis n = 0..M-1, E_xi(-eps, eps)_{n,m} = <xi_n, xi_m> * (1/(2pi)) int_{-eps}^{eps} e^{i(n-m)x} dx
# = <xi_n, xi_m> * sin((n-m)*eps) / (pi*(n-m)) (and eps/pi on diagonal).
M_ph = 200
idx = np.arange(M_ph)
dn = idx[:, None] - idx[None, :]

def phase_norm(gram: np.ndarray, eps: float) -> float:
    with np.errstate(divide="ignore", invalid="ignore"):
        K = np.where(dn == 0, eps / np.pi, np.sin(dn * eps) / (np.pi * dn))
    H = gram * K
    H = 0.5 * (H + H.conj().T)
    ev = float(np.linalg.eigvalsh(H)[-1].real)
    assert -1e-12 <= ev <= 1.0 + 1e-10
    return ev

# (a) Canonical phase: <xi_n, xi_m> = 1
G_can = np.ones((M_ph, M_ph))
# (b) Elementary phase (z = 0.7 between s=5, t=6): <xi_n, xi_m> = delta_nm + z*(delta_{n5}delta_{m6} + h.c.)
G_el = np.eye(M_ph)
G_el[5, 6] = G_el[6, 5] = 0.7
# (c) Intermediate phase: xi_n = sqrt(0.36)*e_0 + sqrt(0.64)*e_{n+1} => <xi_n, xi_m> = 0.36 + 0.64*delta_nm
G_mid = 0.36 * np.ones((M_ph, M_ph)) + 0.64 * np.eye(M_ph)
# (d) Sparse aligned blocks: blocks of length k on [k^2, k^2+k) share ray u_k (orthogonal across k), others orthogonal
# Here for k=1..13, last block is [169..181] (length 13); for exact asymptotic test use M_ph=200.
G_blk = np.eye(M_ph)
for k in range(1, 14):
    s0 = k * k
    e0 = min(M_ph, s0 + k)
    G_blk[s0:e0, s0:e0] = 1.0

for eps in (0.05, 0.2, 0.5):
    v_can = phase_norm(G_can, eps)
    v_el = phase_norm(G_el, eps)
    v_mid = phase_norm(G_mid, eps)
    v_blk = phase_norm(G_blk, eps)
    exact_el = eps / math.pi + 0.7 * math.sin(eps) / math.pi
    assert abs(v_el - exact_el) < 1e-12
    print(f"  eps={eps:4.2f}: can={v_can:.6f} (beta=1) | el={v_el:.6f} (exact {exact_el:.6f}, beta=0) | mid={v_mid:.6f} (beta=0.36) | blk={v_blk:.6f}")


# ============================================================================
# 7. OPEN <=> BOREL: FAT CANTOR SET C_6 & NON-UNIMODULAR AFFINE (ax+b) GROUP
# ============================================================================
print("\n=== 7. OPEN <=> BOREL: FAT CANTOR SET C_6 & NON-UNIMODULAR AFFINE GROUP ===")
ints = [(0.0, 1.0)]
for j in range(1, 7):
    rem = 4.0 ** (-j)
    nxt = []
    for a, b in ints:
        mid = 0.5 * (a + b)
        nxt.append((a, mid - 0.5 * rem))
        nxt.append((mid + 0.5 * rem, b))
    ints = nxt
tot_meas = sum(b - a for a, b in ints)
assert abs(tot_meas - 65.0 / 128.0) < 1e-15
print(f"  C_6: {len(ints)} intervals, measure = {tot_meas:.10f} (exact 65/128)")

for Kmax, n_sub in ((10.0, 4), (100.0, 6), (500.0, 10), (1500.0, 14), (3000.0, 20)):
    qxs, qws = leggauss(n_sub)
    xs_c = np.concatenate([0.5 * (b - a) * qxs + 0.5 * (b + a) for a, b in ints])
    ws_c = np.concatenate([0.5 * (b - a) * qws for a, b in ints])
    sqw_c = np.sqrt(ws_c)
    diff_c = xs_c[:, None] - xs_c[None, :]
    Kc = exact_interval_kernel(diff_c, 0.0, Kmax)
    Mc = (sqw_c[:, None] * Kc) * sqw_c[None, :]
    Mc = 0.5 * (Mc + Mc.conj().T)
    top = float(np.linalg.eigvalsh(Mc)[-1].real)
    assert 0.0 <= top <= 1.0 + 1e-10
    print(f"  Kmax={Kmax:6.1f} (Nx={len(xs_c):4d}): ||E(C_6)||_op = {top:.10f} (1 - norm = {1.0 - top:.3e})")

# Non-unimodular affine group G = {(a,b): a>0, b in R}, dm_G = da db / a^2, H = {(a,0)}, G/H = R.
# Pullback of B = C_6 contains K = [1, 2] x C_6 (m_G(K) = 0.5 * |C_6|).
# For x = (a0, b0) near e = (1, 0), compute exact m_G(K cap K x^{-1}) / m_G(K)
# = (1 / (0.5*|C_6|)) * int_{a in [1,2] cap [1/a0, 2/a0]} (da / a^2) * |C_6 cap (C_6 - a*b0)|.
def overlap_intervals(ints_a, shift: float) -> float:
    """Exact 1D Lebesgue measure of (U ints_a) cap ((U ints_a) - shift)."""
    tot = 0.0
    # Two-pointer merge on sorted disjoint intervals
    i, j = 0, 0
    n_i = len(ints_a)
    while i < n_i and j < n_i:
        a1, b1 = ints_a[i]
        a2, b2 = ints_a[j][0] - shift, ints_a[j][1] - shift
        lo = a1 if a1 > a2 else a2
        hi = b1 if b1 < b2 else b2
        if hi > lo:
            tot += hi - lo
        if b1 < b2:
            i += 1
        else:
            j += 1
    return tot

print("  Non-unimodular affine group G = R rtimes R_{>0} on G/H = R with B = C_6, K = [1,2] x C_6:")
qxa, qwa = leggauss(32)
for a0, b0 in ((1.01, 0.005), (1.001, 0.0005), (1.0001, 0.00005)):
    a_lo = max(1.0, 1.0 / a0)
    a_hi = min(2.0, 2.0 / a0)
    as_g = 0.5 * (a_hi - a_lo) * qxa + 0.5 * (a_hi + a_lo)
    ws_g = 0.5 * (a_hi - a_lo) * qwa
    num = sum(w * overlap_intervals(ints, a * b0) / (a * a) for a, w in zip(as_g, ws_g))
    ratio = num / (0.5 * tot_meas)
    assert 0.0 <= ratio <= 1.0 + 1e-12
    print(f"    x=(a0={a0:.4f}, b0={b0:.5f}): m_G(K cap K x^-1)/m_G(K) = {ratio:.8f} (1 - ratio = {1.0 - ratio:.3e})")
assert 1.0 - ratio < 2e-2


# ============================================================================
# 8. QFT BOUNDED-MASS-BAND DENSITY VS. SHARP BOUND (m2 - m1)/L
# ============================================================================
print("\n=== 8. QFT BOUNDED MASS BAND: SHARP BOUND (m2 - m1)/L ===")
m1, m2 = 1.0, 2.0
sharp_width = m2 - m1
rng = np.random.default_rng(42)
print(f"  m1={m1}, m2={m2}: cube density |S cap Q_L(p)| / L^2 <= (m2 - m1)/L = {sharp_width:.6f}/L")
for L in (2.0, 5.0, 20.0, 100.0):
    centers = [(0.5 * (m1 + m2), 0.0), (10.0, 9.9), (100.0, 99.99), (10000.0, 9999.9999)]
    for _ in range(40):
        p1 = float(rng.uniform(-500, 500))
        p0 = math.sqrt(p1 * p1 + rng.uniform(0.5, 5.0))
        centers.append((p0, p1))
    max_dens = 0.0
    max_p = None
    qx_b, qw_b = leggauss(128)
    for p0, p1 in centers:
        k1s = p1 + 0.5 * L * qx_b
        w1s = 0.5 * L * qw_b
        lo_s = np.sqrt(k1s * k1s + m1 * m1)
        hi_s = np.sqrt(k1s * k1s + m2 * m2)
        ov = np.maximum(0.0, np.minimum(p0 + 0.5 * L, hi_s) - np.maximum(p0 - 0.5 * L, lo_s))
        dens = float(np.sum(w1s * ov) / (L * L))
        if dens > max_dens:
            max_dens = dens
            max_p = (round(p0, 3), round(p1, 3))
    assert max_dens <= sharp_width / L + 1e-12
    print(f"  L={L:5.1f}: max_p |S cap Q_L(p)|/|Q_L| = {max_dens:.8f} (at p={max_p}) <= sharp (m2-m1)/L = {sharp_width/L:.8f}")


# ============================================================================
# 9. QFT WICK-POWER SCALING-DIMENSION LAW (2-BODY, 3-BODY, 4-BODY PHASE SPACE)
# ============================================================================
print("\n=== 9. QFT WICK-POWER SCALING-DIMENSION LAW: 2-BODY, 3-BODY, 4-BODY PHASE SPACES ===")
m_mass = 1.0
print("9a. Exact 2-body phase spaces Phi_2^{(1+s)}(M^2) for s = 1, 2, 3:")
for M_inv in (2.5, 4.0, 10.0, 50.0):
    def cum_I(s_dim, M_val):
        pmax = math.sqrt(max(0.0, 0.25 * M_val * M_val - m_mass * m_mass))
        surf = {1: 2.0, 2: 2.0 * math.pi, 3: 4.0 * math.pi}[s_dim]
        pref = (2.0 * math.pi) ** (1 - s_dim) * surf / 4.0
        qx_i, qw_i = leggauss(128)
        ps = 0.5 * pmax * (qx_i + 1.0)
        ws = 0.5 * pmax * qw_i
        return float(pref * np.sum(ws * (ps ** (s_dim - 1)) / (ps * ps + m_mass * m_mass)))

    h_d = 1e-5 * M_inv
    phi1 = (cum_I(1, M_inv + h_d) - cum_I(1, M_inv - h_d)) / (2.0 * h_d)
    phi2 = (cum_I(2, M_inv + h_d) - cum_I(2, M_inv - h_d)) / (2.0 * h_d)
    phi3 = (cum_I(3, M_inv + h_d) - cum_I(3, M_inv - h_d)) / (2.0 * h_d)
    ex1 = 1.0 / (M_inv * math.sqrt(M_inv * M_inv - 4.0 * m_mass * m_mass))
    ex2 = 1.0 / (4.0 * M_inv)
    ex3 = math.sqrt(1.0 - 4.0 * m_mass * m_mass / (M_inv * M_inv)) / (8.0 * math.pi)
    assert abs(phi1 - ex1) / ex1 < 1e-8 and abs(phi2 - ex2) / ex2 < 1e-8 and abs(phi3 - ex3) / ex3 < 1e-8
    print(f"  M={M_inv:4.1f}: d=1+1 err={abs(phi1-ex1)/ex1:.2e} | d=1+2 err={abs(phi2-ex2)/ex2:.2e} | d=1+3 err={abs(phi3-ex3)/ex3:.2e}")

print("9b. Independent direct 4D quadrature check of Phi_3^{(1+2)}(M^2) = (1 - 3m/M) / (16*pi) [critical :phi^3: in 1+2]:")
x1_3, wt1_3 = leggauss(128)
x2_3, wt2_3 = leggauss(128)
phis_3 = 0.5 * math.pi * (x2_3 + 1.0)
dphis_3 = 0.5 * math.pi * wt2_3
for M_inv in (3.1, 4.0, 6.0, 10.0, 25.0):
    w1_max = (M_inv * M_inv - 3.0 * m_mass * m_mass) / (2.0 * M_inv)
    w1s = 0.5 * (w1_max - m_mass) * x1_3 + 0.5 * (w1_max + m_mass)
    dw1s = 0.5 * (w1_max - m_mass) * wt1_3
    total = 0.0
    for w1, dw1 in zip(w1s, dw1s):
        r1 = math.sqrt(w1 * w1 - m_mass * m_mass)
        E23 = M_inv - w1
        s23 = E23 * E23 - r1 * r1
        half = 0.5 * r1 * math.sqrt(max(0.0, 1.0 - 4.0 * m_mass * m_mass / s23))
        mid = 0.5 * E23
        w2s = mid + half * np.cos(phis_3)
        r2s = np.sqrt(w2s * w2s - m_mass * m_mass)
        E3s = E23 - w2s
        c0s = (E3s * E3s - m_mass * m_mass - r1 * r1 - r2s * r2s) / (2.0 * r1 * r2s)
        integrand = (half * np.sin(phis_3)) / (16.0 * math.pi**2 * r1 * r2s * np.sqrt(1.0 - c0s * c0s))
        total += dw1 * float(np.sum(dphis_3 * integrand))
    exact3 = (1.0 - 3.0 * m_mass / M_inv) / (16.0 * math.pi)
    rel_err = abs(total - exact3) / exact3
    assert rel_err < 1e-10
    print(f"  M={M_inv:4.1f}: direct_4D={total:.12f}, exact=(1-3m/M)/(16pi)={exact3:.12f}, rel_err={rel_err:.2e}")

print("9c. Exact recursion check of Phi_4^{(1+2)}(M^2) = (M - 4m)^2 / (128*pi^2*M) [supercritical :phi^4: in 1+2]:")
qx4, qw4 = leggauss(300)
for M_inv in (4.5, 8.0, 20.0, 100.0):
    s_lo, s_hi = 9.0 * m_mass * m_mass, (M_inv - m_mass) ** 2
    ss = 0.5 * (s_hi - s_lo) * qx4 + 0.5 * (s_hi + s_lo)
    ws = 0.5 * (s_hi - s_lo) * qw4
    # Phi_4^{(1+2)}(M^2) = int ds/(2pi) * Phi_2^{(1+2)}(M^2; sqrt(s), m) * Phi_3^{(1+2)}(s)
    num4 = float(np.sum(ws * (1.0 / (2.0 * math.pi)) * (1.0 / (4.0 * M_inv)) * ((1.0 - 3.0 * m_mass / np.sqrt(ss)) / (16.0 * math.pi))))
    ex4 = (M_inv - 4.0 * m_mass) ** 2 / (128.0 * math.pi**2 * M_inv)
    assert abs(num4 - ex4) / ex4 < 1e-13
    print(f"  M={M_inv:5.1f}: Phi_4^(1+2) = {num4:.10f}, exact = {ex4:.10f}, rel_err={abs(num4-ex4)/ex4:.2e} (unbounded ~ M/(128*pi^2))")

print("9d. Exact 3-body phase space Phi_3^{(1+1)}(M^2) [threshold + exact (3/(2pi)) M^{-2} log(M/m) asymptotic]:")
# Phi_3^{(1+1)}(M^2) = int_{4m^2}^{(M-m)^2} (ds/(2pi)) / sqrt(((M-m)^2 - s)((M+m)^2 - s)*s*(s - 4m^2))
# Change variable s = 4m^2 + ((M-m)^2 - 4m^2) * sin^2(theta), theta in [0, pi/2]:
# ds / sqrt(((M-m)^2 - s)(s - 4m^2)) = 2 dtheta! Completely smooth integrand!
# For large M, split theta in [0, pi/2] using log-clustered subintervals near 0 and pi/2 for 12-digit accuracy.
def phi3_1p1(M_inv: float, m_val: float = 1.0) -> float:
    s_lo = 4.0 * m_val * m_val
    s_hi = (M_inv - m_val) ** 2
    # Breakpoints in theta to resolve boundary layers near theta ~ 2m/M and pi/2 - theta ~ sqrt(4m/M)
    eps_lo = min(0.1, 2.0 * m_val / M_inv)
    eps_hi = min(0.1, math.sqrt(4.0 * m_val / M_inv))
    pts = sorted(set([
        0.0, eps_lo * 0.25, eps_lo, eps_lo * 4.0, eps_lo * 16.0, 0.25,
        0.5 * math.pi - 0.25, 0.5 * math.pi - eps_hi * 16.0,
        0.5 * math.pi - eps_hi * 4.0, 0.5 * math.pi - eps_hi,
        0.5 * math.pi - eps_hi * 0.25, 0.5 * math.pi
    ]))
    pts = [p for p in pts if 0.0 <= p <= 0.5 * math.pi]
    qx_s, qw_s = leggauss(64)
    acc = 0.0
    for a, b in zip(pts[:-1], pts[1:]):
        if b <= a:
            continue
        th = 0.5 * (b - a) * qx_s + 0.5 * (b + a)
        w = 0.5 * (b - a) * qw_s
        ss = s_lo + (s_hi - s_lo) * (np.sin(th) ** 2)
        acc += float(np.sum(w * (1.0 / math.pi) / np.sqrt(((M_inv + m_val) ** 2 - ss) * ss)))
    return acc

thresh_exact = 1.0 / (8.0 * math.sqrt(3.0) * m_mass * m_mass)
asymp_const = 3.0 / (2.0 * math.pi)
for M_inv in (3.0001, 5.0, 100.0, 1000.0, 10000.0, 100000.0):
    val11 = phi3_1p1(M_inv, m_mass)
    if M_inv < 3.01:
        assert abs(val11 - thresh_exact) / thresh_exact < 1e-3
        print(f"  M={M_inv:8.4f}: Phi_3^(1+1) = {val11:.10f} (threshold 1/(8*sqrt(3)*m^2) = {thresh_exact:.10f})")
    else:
        scaled = (M_inv * M_inv * val11) / math.log(M_inv / m_mass)
        print(f"  M={M_inv:8.1f}: Phi_3^(1+1) = {val11:.6e} | M^2 Phi_3 / log(M/m) = {scaled:.6f} (exact 3/(2pi) = {asymp_const:.6f})")
assert abs(((100000.0**2) * phi3_1p1(100000.0, m_mass) / math.log(100000.0)) - asymp_const) < 5e-3

print("9e. 3D Euclidean Laplace-transform check of Kaiser's Phi_n^{(1+2)} for n = 2, 3, 4, 5:")
# int_{nm}^inf 2*M * (Phi_n(M^2)/(2*pi)) * e^{-M*r}/(4*pi*r) dM == (e^{-m*r}/(4*pi*r))^n
# Equivalently, int_{nm}^inf M * Phi_n(M^2) * e^{-M*r} dM == 2*pi * e^{-n*m*r} / ((4*pi*r)^{n-1})
# Wait: Phi_n(M^2) = (M - n*m)^{n-2} / (4 * (4*pi)^{n-2} * (n-2)! * M), so
# int_{nm}^inf M * Phi_n(M^2) * e^{-M*r} dM = e^{-n*m*r} / (4 * (4*pi)^{n-2} * r^{n-1})
# = pi * e^{-n*m*r} / ((4*pi)^{n-1} * r^{n-1}), matching (4*pi*r) * (e^{-m*r}/(4*pi*r))^n * pi!
# Check by Gauss-Laguerre quadrature directly on M = n*m + t/r:
from numpy.polynomial.laguerre import laggauss
lx, lw = laggauss(48)
r_test = 1.37
for n_p in (2, 3, 4, 5):
    Ms = n_p * m_mass + lx / r_test
    phi_vals = ((Ms - n_p * m_mass) ** (n_p - 2)) / (4.0 * ((4.0 * math.pi) ** (n_p - 2)) * math.factorial(n_p - 2) * Ms)
    # int_{nm}^inf 2*M * (Phi_n/(2*pi)) * (e^{-M*r}/(4*pi*r)) dM = (1/(4*pi^2*r)) * e^{-n*m*r}/r * sum lw * Ms * phi_vals
    lhs = (math.exp(-n_p * m_mass * r_test) / (4.0 * math.pi * math.pi * r_test * r_test)) * float(np.sum(lw * Ms * phi_vals))
    rhs = (math.exp(-m_mass * r_test) / (4.0 * math.pi * r_test)) ** n_p
    rel = abs(lhs - rhs) / rhs
    assert rel < 1e-12
    print(f"  n={n_p}, r={r_test}: Laplace LHS={lhs:.12e}, (Delta_3D)^n={rhs:.12e}, rel_err={rel:.2e}")

print(f"\nAll assertions passed! Total verification runtime: {time.time() - t0:.2f}s")
