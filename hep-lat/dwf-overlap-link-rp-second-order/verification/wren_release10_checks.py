"""Standalone verification of the six Wren-computed figures in Release 10:
  "Link reflection positivity at second order: domain-wall fermions at finite L5,
   a uniform inner-moment bound, and the overlap limit"

Self-contained (numpy + scipy only; no external project imports).
Covers:
  1. s1(a), s3.2, s6: tr N~_PV scaling (0.05 S to 0.076 S for (m, p1) in {0.7, 1.0} x {0, pi})
     and finite a5 -> 0 limit at fixed S.
  2. s3.2: Weight normalisation (unshifted sum W~_j / mu_j^2 vs shifted sum W~_j at S = 4,
     L5 = 8..64: 39.7 vs 0.248 at L5 = 8, 1346 vs 0.255 at L5 = 64).
  3. s3.2 Remark (m >= 1), s6: at a5 = 1 and m > 1, rho_PV(m, 1, p1) can exceed rho_*(m, p1)
     at intermediate p1 while staying <= rho_*(m) (with m = 0.7 control).
  4. s3.3 Lemma 2 (H2): persistent spatial current d_y log det D_DW at spatial twist pi/3
     (and vanishing at twists 0 and pi), checked on open, antiperiodic, and periodic time BC.
  5. s4 Corollary 3.1(2): 120 random points bounding sum_j |A1(1,j)|^2 <= 0.354 and
     sum_j |A0(1,j)|^2 <= 0.173 (analytic bound <= 4 each, C_s <= 8).
  6. s5.1 Order of limits ("OR p1 != 0"): exact Ginsparg-Wilson pair cancellation on 100
     random points, plus theta_s = 0, Ls = 4, m = 0.7 antiperiodic-time overlap vs L = inf
     DW-PV (IR) at p1 = 0, pi/2, pi.
"""
import platform
import time
import numpy as np
import scipy
import scipy.linalg as la

# Conventions matching Section 2 & Section 4.1: gamma0 = sigma1, gamma1 = sigma2, gamma5 = sigma3
s1 = np.array([[0.0, 1.0], [1.0, 0.0]], dtype=complex)
s2 = np.array([[0.0, -1j], [1j, 0.0]], dtype=complex)
s3 = np.array([[1.0, 0.0], [0.0, -1.0]], dtype=complex)
I2 = np.eye(2, dtype=complex)
g0, g1, g5 = s1, s2, s3

M_fwd_0 = (-0.5j) * (I2 - g0)
M_bwd_0 = (+0.5j) * (I2 + g0)
M_fwd_1 = (-0.5j) * (I2 - g1)
M_bwd_1 = (+0.5j) * (I2 + g1)

# Unitary to gamma0 eigenbasis: U0^dag g0 U0 = diag(+1, -1)
U0 = (1.0 / np.sqrt(2.0)) * np.array([[1.0, 1.0], [1.0, -1.0]], dtype=complex)


def E_cut(k1, m):
    c = 2.0 - m - np.cos(k1)
    if abs(c) < 1e-14:
        return np.inf
    ch = (c * c + 1.0 + np.sin(k1) ** 2) / (2.0 * abs(c))
    return float(np.arccosh(max(1.0, ch)))


def r_star(k1, m):
    ec = E_cut(k1, m)
    return 0.0 if np.isinf(ec) else float(np.exp(-ec))


def rho_star_mp(m, p1, n_grid=2048):
    ks = np.linspace(0.0, 2.0 * np.pi, n_grid, endpoint=False)
    return float(max(r_star(k, m) * r_star(k + p1, m) for k in ks))


def shamir_poles_and_residues(k1, m, L5, a5, ap=False):
    """Exact L5 real temporal poles lam in (-1, 1) of D_X(z, k1) and PSD residue factors V_orig
    (Res_{z=lam} (g0 D_X(z, k1))^{-1} = sum_{in cluster} v v^dag in original (s, spin) basis).
    """
    S = np.zeros((L5, L5), dtype=complex)
    for s in range(L5 - 1):
        S[s, s + 1] = 1.0
    if ap:
        S[L5 - 1, 0] = -1.0

    c = 2.0 - m - np.cos(k1)
    sv = np.sin(k1)
    M0 = (1.0 + a5 * c) * np.eye(L5, dtype=complex) - 0.5 * (S + S.conj().T)
    K = 0.5 * (S - S.conj().T) - (a5 * sv) * np.eye(L5, dtype=complex)
    aI = a5 * np.eye(L5, dtype=complex)
    Z = np.zeros((L5, L5), dtype=complex)

    A_pen = np.block([[aI, Z], [-K.conj().T, M0]])
    B_pen = np.block([[M0, K], [Z, aI]])

    vals, vecs = la.eig(A_pen, B_pen)
    inside = np.where((np.isfinite(vals)) & (np.abs(vals) < 1.0 - 1e-12))[0]
    if len(inside) != L5:
        raise RuntimeError(f"Expected {L5} poles inside |z|<1, got {len(inside)}")
    if float(np.max(np.abs(vals[inside].imag))) > 1e-9:
        raise RuntimeError("Non-real pole encountered")

    lams = vals[inside].real
    V_in = vecs[:, inside]
    order = np.argsort(-np.abs(lams))
    lams = lams[order]
    V_in = V_in[:, order]

    T_to_orig = np.zeros((2 * L5, 2 * L5), dtype=complex)
    for s in range(L5):
        T_to_orig[2 * s:2 * s + 2, s] = U0[:, 0]
        T_to_orig[2 * s:2 * s + 2, L5 + s] = U0[:, 1]

    clusters = []
    used = np.zeros(L5, dtype=bool)
    for i in range(L5):
        if used[i]:
            continue
        cl = [i]
        used[i] = True
        for j in range(i + 1, L5):
            if not used[j] and abs(lams[j] - lams[i]) <= 1e-9 * max(1.0, abs(lams[i])):
                cl.append(j)
                used[j] = True
        clusters.append(cl)

    pole_lams = np.zeros(L5, dtype=float)
    V_orig = np.zeros((2 * L5, L5), dtype=complex)
    idx_out = 0
    for cl in clusters:
        lam = float(np.mean(lams[cl]))
        W = V_in[:, cl]
        w_diag = np.concatenate([np.full(L5, a5 / (lam * lam)), np.full(L5, a5)])
        Q_orth, _ = np.linalg.qr(W * np.sqrt(w_diag)[:, None])
        W_orig = T_to_orig @ (Q_orth / np.sqrt(w_diag)[:, None])
        for k in range(len(cl)):
            pole_lams[idx_out] = lam
            V_orig[:, idx_out] = W_orig[:, k]
            idx_out += 1
    return pole_lams, V_orig


def Q_operator(L5):
    R5 = np.zeros((L5, L5), dtype=complex)
    for s in range(L5):
        R5[s, L5 - 1 - s] = 1.0
    return np.kron(R5, g5 @ g0)


def exact_pair_spectral_weights(k1, k1p, m, L5, a5, ap=False, include_mu0=True):
    l1, V1 = shamir_poles_and_residues(k1, m, L5, a5, ap=ap)
    if abs((k1p - k1) % (2 * np.pi)) < 1e-13:
        l2, V2 = l1, V1
    else:
        l2, V2 = shamir_poles_and_residues(k1p, m, L5, a5, ap=ap)

    Q = Q_operator(L5)
    QV2 = Q @ V2
    W1_loc = a5 * (np.exp(1j * k1p) * M_fwd_1 + np.exp(-1j * k1) * M_bwd_1)
    g0W1 = np.kron(np.eye(L5, dtype=complex), g0 @ W1_loc)
    A1 = V1.conj().T @ (g0W1 @ QV2)
    W11 = (np.abs(A1) ** 2).ravel()
    Lams = np.outer(l1, l2).ravel()
    if not include_mu0:
        return Lams, W11, None

    P0m = np.kron(np.eye(L5, dtype=complex), 0.5 * (I2 - g0))
    P0p = np.kron(np.eye(L5, dtype=complex), 0.5 * (I2 + g0))
    Am = V1.conj().T @ (P0m @ QV2)
    Ap = V1.conj().T @ (P0p @ QV2)
    A0 = (-1j * a5) * (l1[:, None] * Am + l2[None, :] * Ap)
    a1 = A1.ravel()
    a0 = A0.ravel()
    W_mat = np.zeros((len(Lams), 2, 2), dtype=complex)
    W_mat[:, 0, 0] = np.abs(a1) ** 2
    W_mat[:, 0, 1] = a1 * np.conj(a0)
    W_mat[:, 1, 0] = a0 * np.conj(a1)
    W_mat[:, 1, 1] = np.abs(a0) ** 2
    return Lams, W11, W_mat


def check_1_tr_Npv_scaling():
    print("=== CHECK 1 (s1, s3.2, s6): tr N~_PV scaling in S = L5*a5 and finite a5 -> 0 limit ===")
    Ls_q = 32
    ratios_22 = []
    for m in (0.7, 1.0):
        for plbl, p1 in [("p1=0", 0.0), ("p1=pi", np.pi)]:
            for a5, L5s in [(0.25, [16, 32, 64]), (0.125, [32, 64])]:
                for L5 in L5s:
                    S = L5 * a5
                    tr11 = 0.0
                    tr22 = 0.0
                    for ix in range(Ls_q):
                        k1 = 2.0 * np.pi * ix / Ls_q
                        k1p = k1 + p1
                        _, W11, Wmat = exact_pair_spectral_weights(k1, k1p, m, L5, a5, ap=True, include_mu0=True)
                        tr11 += float(np.sum(W11)) / Ls_q
                        tr22 += float(np.sum(Wmat[:, 0, 0].real + Wmat[:, 1, 1].real)) / Ls_q
                    ratios_22.append(tr22 / S)
                    print(f"  m={m:.1f} {plbl:5s} a5={a5:.3f} L5={L5:2d} S={S:4.1f}: "
                          f"tr_N~_PV(11)={tr11:.5f} tr_N~_PV(2x2)={tr22:.5f} tr(2x2)/S={tr22/S:.5f}")
    print(f"  Range of tr N~_PV(2x2)/S across S in {{4,8,16}}: [{min(ratios_22):.5f}, {max(ratios_22):.5f}]")

    print("  Fixed-S a5 -> 0 convergence (m=0.7, p1=0, Ls_q=32):")
    for S in (4.0, 8.0):
        vals = []
        for a5 in (0.5, 0.25, 0.125, 0.0625):
            L5 = int(round(S / a5))
            tr22 = 0.0
            for ix in range(Ls_q):
                k1 = 2.0 * np.pi * ix / Ls_q
                _, _, Wmat = exact_pair_spectral_weights(k1, k1, 0.7, L5, a5, ap=True, include_mu0=True)
                tr22 += float(np.sum(Wmat[:, 0, 0].real + Wmat[:, 1, 1].real)) / Ls_q
            vals.append(f"a5={a5:.4f}(L5={L5:3d}):{tr22:.5f}")
        print(f"    S={S:.1f}: " + "  ".join(vals))


def check_2_weight_normalisation():
    print("\n=== CHECK 2 (s3.2): Weight normalisation sum W~_j / mu_j^2 vs sum W~_j (m=0.7, p1=0, S=4) ===")
    m = 0.7
    Ls_q = 16
    rho = rho_star_mp(m, 0.0)
    for a5 in (0.5, 0.25, 0.125, 0.0625):
        L5 = int(round(4.0 / a5))
        n_shifted = 0.0
        n_unshifted = 0.0
        min_mu = 1.0
        for ix in range(Ls_q):
            k1 = 2.0 * np.pi * ix / Ls_q
            Lams, W11, _ = exact_pair_spectral_weights(k1, k1, m, L5, a5, ap=True, include_mu0=False)
            n_shifted += float(np.sum(W11)) / Ls_q
            n_unshifted += float(np.sum(W11 / (Lams ** 2))) / Ls_q
            min_mu = min(min_mu, float(np.min(np.abs(Lams))))
        b_shifted = (rho ** 2) / (1.0 - rho ** 2) * n_shifted
        b_unshifted = (rho ** 4) / (1.0 - rho ** 2) * n_unshifted
        print(f"  L5={L5:2d} a5={a5:.4f} (S=4): min|mu|={min_mu:.3e} | "
              f"sum W~_j={n_shifted:.5f} | sum W~_j/mu_j^2={n_unshifted:.2f} | "
              f"t0=2 ratio={b_unshifted/b_shifted:.1f}x")


def check_3_m_gt_1_a5_1():
    print("\n=== CHECK 3 (s3.2 Remark m>=1, s6): rho_PV(m, 1, p1) vs rho_*(m, p1) and rho_*(m) at a5 = 1 ===")
    n_k = 360
    ks = np.linspace(-np.pi, np.pi, n_k, endpoint=False)
    k5s = np.linspace(-np.pi, np.pi, n_k, endpoint=False)
    K1, K5 = np.meshgrid(ks, k5s, indexing="ij")
    a5 = 1.0
    u5 = (1.0 - np.cos(K5)) / a5
    for m in (0.7, 1.0, 1.3):
        c = 2.0 - m - np.cos(K1)
        den = 2.0 * np.abs(u5 + c)
        with np.errstate(divide="ignore", invalid="ignore"):
            ch_pv = (c * c + 1.0 + np.sin(K1) ** 2 + 2.0 * u5 * (c + 1.0 / a5)) / den
            E_pv = np.where(den < 1e-14, np.inf, np.arccosh(np.maximum(1.0, ch_pv)))
        r_pv_k1 = np.exp(-E_pv)  # (n_k, n_k)
        r_cut_k1 = np.array([r_star(k, m) for k in ks])
        rho_m = float(np.max(r_cut_k1) ** 2)
        max_excess_p1 = -np.inf
        worst_p1 = 0.0
        worst_vals = (0.0, 0.0)
        max_rho_pv_all = 0.0
        for sh in range(n_k):
            # rho_PV(m, 1, p1) = max_{k1, k5} r_pv(k1, k5) * r_pv(k1 + p1, k5)
            rpv_p1 = float(np.max(r_pv_k1 * np.roll(r_pv_k1, -sh, axis=0)))
            rstar_p1 = float(np.max(r_cut_k1 * np.roll(r_cut_k1, -sh)))
            max_rho_pv_all = max(max_rho_pv_all, rpv_p1)
            diff = rpv_p1 - rstar_p1
            if diff > max_excess_p1:
                max_excess_p1 = diff
                worst_p1 = 2.0 * np.pi * sh / n_k
                worst_vals = (rpv_p1, rstar_p1)
        print(f"  m={m:.1f} a5=1: rho_*(m)={rho_m:.6f}, max_p1 rho_PV={max_rho_pv_all:.6f} "
              f"(diff vs rho_*(m)={max_rho_pv_all - rho_m:+.2e}) | "
              f"max_p1 [rho_PV(p1) - rho_*(p1)]={max_excess_p1:+.4e} at p1={worst_p1/np.pi:.3f}*pi "
              f"(rho_PV={worst_vals[0]:.6f}, rho_*={worst_vals[1]:.6f})")


def build_2d_wilson_and_5d(L, Ls, U0, U1, m, L5, a5, ap=False, time_bc="open"):
    """Direct 2D Wilson + 5D Shamir matrix builder (zero external imports).
    time_bc in {'open', 'ap', 'periodic'}.
    """
    T = 2 * L
    n2_sites = T * Ls
    n2 = 2 * n2_sites
    Dw = np.zeros((n2, n2), dtype=complex)
    for ti in range(T):
        for x in range(Ls):
            i = ti * Ls + x
            Dw[2 * i:2 * i + 2, 2 * i:2 * i + 2] += 2.0 * I2
            # Time hop ti -> ti + 1
            if ti + 1 < T or time_bc != "open":
                tj = (ti + 1) % T
                j = tj * Ls + x
                sgn = -1.0 if (ti == T - 1 and time_bc == "ap") else 1.0
                u = U0[ti, x] * sgn
                Dw[2 * i:2 * i + 2, 2 * j:2 * j + 2] += -0.5 * (I2 - g0) * u
                Dw[2 * j:2 * j + 2, 2 * i:2 * i + 2] += -0.5 * (I2 + g0) * np.conj(u)
            # Space hop x -> (x + 1) % Ls
            j = ti * Ls + ((x + 1) % Ls)
            u = U1[ti, x]
            Dw[2 * i:2 * i + 2, 2 * j:2 * j + 2] += -0.5 * (I2 - g1) * u
            Dw[2 * j:2 * j + 2, 2 * i:2 * i + 2] += -0.5 * (I2 + g1) * np.conj(u)

    G5 = np.kron(np.eye(n2_sites, dtype=complex), g5)
    Pm = 0.5 * (np.eye(n2, dtype=complex) - G5)
    Pp = 0.5 * (np.eye(n2, dtype=complex) + G5)
    B = np.eye(n2, dtype=complex) + a5 * (Dw - m * np.eye(n2, dtype=complex))
    D5 = np.zeros((L5 * n2, L5 * n2), dtype=complex)
    for s in range(L5):
        D5[s * n2:(s + 1) * n2, s * n2:(s + 1) * n2] = B
        if s + 1 < L5:
            D5[s * n2:(s + 1) * n2, (s + 1) * n2:(s + 2) * n2] = -Pm
            D5[(s + 1) * n2:(s + 2) * n2, s * n2:(s + 1) * n2] = -Pp
    if ap:
        D5[(L5 - 1) * n2:L5 * n2, 0:n2] = +Pm
        D5[0:n2, (L5 - 1) * n2:L5 * n2] = +Pp
    return D5


def check_4_twist_persistent_current():
    print("\n=== CHECK 4 (s3.3 Lemma 2 (H2)): persistent current d_y log det D_X at twist pi/3 ===")
    L, Ls, m, L5, a5 = 6, 4, 0.7, 6, 0.25
    eps = 1e-5
    for time_bc in ("open", "periodic", "ap"):
        print(f"  -- time_bc = {time_bc} (L={L}, Ls={Ls}, m={m}, L5={L5}, a5={a5}) --")
        for theta_s, lbl in [(0.0, "theta_s=0"), (np.pi, "theta_s=pi"), (np.pi / 3, "theta_s=pi/3")]:
            def logdet_pair(x_val, y_val, ap):
                U0 = np.ones((2 * L, Ls), dtype=complex)
                U1 = np.full((2 * L, Ls), np.exp(1j * theta_s / Ls), dtype=complex)
                U1[L - 1, :] *= np.exp(1j * x_val)
                U1[L, :] *= np.exp(1j * y_val)
                D5 = build_2d_wilson_and_5d(L, Ls, U0, U1, m, L5, a5, ap=ap, time_bc=time_bc)
                s, ld = np.linalg.slogdet(D5)
                return s.real, ld

            for ap, name in [(False, "DW"), (True, "PV")]:
                s0, _ = logdet_pair(0.0, 0.0, ap)
                _, ldp = logdet_pair(0.0, +eps, ap)
                _, ldm = logdet_pair(0.0, -eps, ap)
                d1_y = (ldp - ldm) / (2.0 * eps)
                _, lpp = logdet_pair(+eps, +eps, ap)
                _, lpm = logdet_pair(+eps, -eps, ap)
                _, lmp = logdet_pair(-eps, +eps, ap)
                _, lmm = logdet_pair(-eps, -eps, ap)
                d2_log = (lpp - lpm - lmp + lmm) / (4.0 * eps * eps)
                d2_det = d2_log + d1_y * d1_y
                print(f"    {lbl:12s} {name}: sign={s0:+.0f} d_y(log det)={d1_y:+.4e} "
                      f"d_x d_y(log det)={d2_log:+.4e} K_xy/K={d2_det:+.4e}")


def check_5_cor31_sums():
    print("\n=== CHECK 5 (s4 Cor 3.1(2)): 120 random points for sum_j |A1(1,j)|^2 and sum_j |A0(1,j)|^2 ===")
    rng = np.random.default_rng(12345)
    max_w1_norm = 0.0
    max_a1_sum = 0.0
    max_a0_sum = 0.0
    max_tr = 0.0
    max_a5R1 = 0.0
    for _ in range(120):
        m = rng.uniform(0.05, 0.95)
        a5_max = 2.0 * (1.0 - m) / (m * (2.0 - m))
        a5 = rng.uniform(0.02, min(1.0, 0.98 * a5_max))
        L5 = int(rng.choice([1, 2, 4, 8, 16]))
        k1 = rng.uniform(-np.pi, np.pi)
        p1 = rng.uniform(-np.pi, np.pi)
        k1p = k1 + p1
        W1_loc = np.exp(1j * k1p) * M_fwd_1 + np.exp(-1j * k1) * M_bwd_1
        max_w1_norm = max(max_w1_norm, float(np.linalg.norm(W1_loc, 2)))
        l1, V1 = shamir_poles_and_residues(k1, m, L5, a5, ap=False)
        l2, V2 = shamir_poles_and_residues(k1p, m, L5, a5, ap=False)
        Q = Q_operator(L5)
        QV2 = Q @ V2
        v1 = V1[:, 0]
        max_a5R1 = max(max_a5R1, float(a5 * np.vdot(v1, v1).real))
        g0W1 = np.kron(np.eye(L5, dtype=complex), g0 @ (a5 * W1_loc))
        A1_row = v1.conj() @ (g0W1 @ QV2)
        P0m = np.kron(np.eye(L5, dtype=complex), 0.5 * (I2 - g0))
        P0p = np.kron(np.eye(L5, dtype=complex), 0.5 * (I2 + g0))
        A0_row = (-1j * a5) * (l1[0] * (v1.conj() @ (P0m @ QV2)) + l2 * (v1.conj() @ (P0p @ QV2)))
        s1_v = float(np.sum(np.abs(A1_row) ** 2))
        s0_v = float(np.sum(np.abs(A0_row) ** 2))
        max_a1_sum = max(max_a1_sum, s1_v)
        max_a0_sum = max(max_a0_sum, s0_v)
        max_tr = max(max_tr, s1_v + s0_v)
    print(f"  120 random points: max ||W1/a5||_op = {max_w1_norm:.10f}, max ||a5 R1||_op = {max_a5R1:.6f}")
    print(f"  max sum_j |A1(1,j)|^2 = {max_a1_sum:.6f} (bound <= 4), "
          f"max sum_j |A0(1,j)|^2 = {max_a0_sum:.6f} (bound <= 4), max sum = {max_tr:.6f} (bound <= 8)")


def overlap_Pi_per_k1_exact(L, Ls, m, p1_idx, theta_s=0.0, mu_mass=0.0):
    nt = 2 * L
    p0s = 2.0 * np.pi * (np.arange(nt) + 0.5) / nt
    k1s = (2.0 * np.pi * np.arange(Ls) + theta_s) / Ls
    w_grid = np.zeros((nt, Ls), dtype=float)
    Pp_grid = np.zeros((nt, Ls, 2, 2), dtype=complex)
    Pm_grid = np.zeros((nt, Ls, 2, 2), dtype=complex)
    G_grid = np.zeros((nt, Ls, 2, 2), dtype=complex)
    c_ov = 0.5 * (1.0 - 0.5 * mu_mass)
    for it, p0 in enumerate(p0s):
        for ix, k1 in enumerate(k1s):
            b = (1.0 - np.cos(p0)) + (1.0 - np.cos(k1)) - m
            X = b * I2 + 1j * np.sin(p0) * g0 + 1j * np.sin(k1) * g1
            H = g5 @ X
            w = np.sqrt(b * b + np.sin(p0) ** 2 + np.sin(k1) ** 2)
            S = H / w
            w_grid[it, ix] = w
            Pp_grid[it, ix] = 0.5 * (I2 + S)
            Pm_grid[it, ix] = 0.5 * (I2 - S)
            Dov = c_ov * (I2 + g5 @ S) + (0.5 * mu_mass) * I2
            G_grid[it, ix] = np.linalg.inv(Dov)

    Pi_per_k1 = np.zeros((Ls, nt, 2, 2), dtype=complex)
    for ix in range(Ls):
        ixp = (ix + p1_idx) % Ls
        k1 = k1s[ix]
        k1p = k1s[ixp]
        for m0 in range(nt):
            itp = (np.arange(nt) + m0) % nt
            p0 = p0s
            p0p = p0s[itp]
            V_plus = np.zeros((2, nt, 2, 2), dtype=complex)
            V1_mat = g5 @ (np.exp(1j * k1) * M_fwd_1 + np.exp(-1j * k1p) * M_bwd_1)
            V_plus[0] = V1_mat[None, :, :]
            V_plus[1] = (np.exp(1j * p0)[:, None, None] * (g5 @ M_fwd_0)[None, :, :]
                         + np.exp(-1j * p0p)[:, None, None] * (g5 @ M_bwd_0)[None, :, :])
            V_minus = np.conj(np.swapaxes(V_plus, -1, -2))

            w1 = w_grid[:, ix]
            w2 = w_grid[itp, ixp]
            Pp1 = Pp_grid[:, ix]
            Pm1 = Pm_grid[:, ix]
            Pp2 = Pp_grid[itp, ixp]
            Pm2 = Pm_grid[itp, ixp]
            G1 = G_grid[:, ix]
            G2 = G_grid[itp, ixp]

            f1 = (2.0 * c_ov / (w1 + w2))[:, None, None]
            dDa = np.zeros((2, nt, 2, 2), dtype=complex)
            dDb = np.zeros((2, nt, 2, 2), dtype=complex)
            for mu in range(2):
                dDa[mu] = f1 * (g5[None] @ (Pp2 @ V_plus[mu] @ Pm1 + Pm2 @ V_plus[mu] @ Pp1))
                dDb[mu] = f1 * (g5[None] @ (Pp1 @ V_minus[mu] @ Pm2 + Pm1 @ V_minus[mu] @ Pp2))

            c_same1 = (2.0 / ((w1 + w2) ** 2))[:, None, None]
            c_flip1 = (1.0 / (w1 * (w1 + w2)))[:, None, None]
            c_flip2 = (1.0 / (w2 * (w1 + w2)))[:, None, None]

            for mu in range(2):
                for nu in range(2):
                    bub = -np.einsum("nij,njk,nkl,nli->", G2, dDa[mu], G1, dDb[nu])
                    Vb = V_minus[nu]
                    Va = V_plus[mu]
                    A1 = Vb @ Pm2 @ Va
                    B1 = Vb @ Pp2 @ Va
                    dd_1 = (-c_same1 * (Pp1 @ A1 @ Pp1)
                            + c_same1 * (Pm1 @ B1 @ Pm1)
                            - c_flip1 * (Pp1 @ B1 @ Pm1 + Pm1 @ B1 @ Pp1)
                            + c_flip1 * (Pp1 @ A1 @ Pm1 + Pm1 @ A1 @ Pp1))
                    tad1 = c_ov * np.einsum("nij,jk,nki->", G1, g5, dd_1)
                    A2 = Va @ Pm1 @ Vb
                    B2 = Va @ Pp1 @ Vb
                    dd_2 = (-c_same1 * (Pp2 @ A2 @ Pp2)
                            + c_same1 * (Pm2 @ B2 @ Pm2)
                            - c_flip2 * (Pp2 @ B2 @ Pm2 + Pm2 @ B2 @ Pp2)
                            + c_flip2 * (Pp2 @ A2 @ Pm2 + Pm2 @ A2 @ Pp2))
                    tad2 = c_ov * np.einsum("nij,jk,nki->", G2, g5, dd_2)
                    Pi_per_k1[ix, m0, mu, nu] = (bub + tad1 + tad2) / (nt * Ls)
    return Pi_per_k1


def check_6_order_of_limits_p1_ne_0():
    print("\n=== CHECK 6 (s5.1 Order of limits): GW pair cancellation & AP vs IR at theta_s=0, Ls=4, m=0.7 ===")
    rng = np.random.default_rng(77)
    max_comm = 0.0
    max_tr_id = 0.0
    max_total = 0.0
    for _ in range(100):
        m = rng.uniform(0.1, 0.9)
        q0 = rng.uniform(-np.pi, np.pi)
        p1 = rng.uniform(0.1, np.pi)
        b2 = (1.0 - np.cos(q0)) + (1.0 - np.cos(p1)) - m
        X2 = b2 * I2 + 1j * np.sin(q0) * g0 + 1j * np.sin(p1) * g1
        H2 = g5 @ X2
        w2 = np.sqrt(b2 * b2 + np.sin(q0) ** 2 + np.sin(p1) ** 2)
        S2 = H2 / w2
        Pp2 = 0.5 * (I2 + S2)
        Pm2 = 0.5 * (I2 - S2)
        D2 = 0.5 * (I2 + g5 @ S2)
        D2inv = np.linalg.inv(D2)
        S1 = -g5
        w1 = m
        Pp1 = 0.5 * (I2 + S1)
        Pm1 = 0.5 * (I2 - S1)
        Vp = [g5 @ (M_fwd_1 + np.exp(-1j * p1) * M_bwd_1),
              g5 @ (M_fwd_0 + np.exp(-1j * q0) * M_bwd_0)]
        Vm = [v.conj().T for v in Vp]
        for mu in range(2):
            for nu in range(2):
                dSa = (2.0 / (w1 + w2)) * (Pp2 @ Vp[mu] @ Pm1 + Pm2 @ Vp[mu] @ Pp1)
                dSb = (2.0 / (w1 + w2)) * (Pp1 @ Vm[nu] @ Pm2 + Pm1 @ Vm[nu] @ Pp2)
                K = dSa @ dSb
                max_comm = max(max_comm, float(np.linalg.norm(S2 @ K - K @ S2)))
                max_tr_id = max(max_tr_id, float(abs(np.trace(D2inv @ K) - np.trace(K))))
                dDa = 0.5 * g5 @ dSa
                dDb = 0.5 * g5 @ dSb
                bub = -np.trace(D2inv @ dDa @ dDb)
                c_same = 2.0 / ((w1 + w2) ** 2)
                c_flip = 1.0 / (w1 * (w1 + w2))
                A1 = Vm[nu] @ Pm2 @ Vp[mu]
                B1 = Vm[nu] @ Pp2 @ Vp[mu]
                dd1 = (-c_same * (Pp1 @ A1 @ Pp1) + c_same * (Pm1 @ B1 @ Pm1)
                       - c_flip * (Pp1 @ B1 @ Pm1 + Pm1 @ B1 @ Pp1)
                       + c_flip * (Pp1 @ A1 @ Pm1 + Pm1 @ A1 @ Pp1))
                tad = 0.5 * np.trace(g5 @ dd1)
                max_total = max(max_total, float(abs(bub + tad)))
    print(f"  GW pair identity (100 random pts): max||[S2,K]||={max_comm:.2e}, "
          f"max|Tr(D2^-1 K)-Tr(K)|={max_tr_id:.2e}, max|bub+tad|={max_total:.2e}")

    m = 0.7
    Ls = 4
    ta = tb = 2
    L = 64
    nt = 2 * L
    q0s = 2.0 * np.pi * np.arange(nt) / nt
    for p1_idx, plbl in [(0, "p1=0"), (1, "p1=pi/2"), (2, "p1=pi")]:
        Pi_ap = overlap_Pi_per_k1_exact(L, Ls, m, p1_idx, theta_s=0.0, mu_mass=0.0)
        Pi_tot = np.sum(Pi_ap, axis=0)
        H_ap = np.zeros((2, 2), dtype=complex)
        for ia in range(2):
            t_ref = (1 - ta) if ia == 0 else (-ta)
            sgn = 1.0 if ia == 0 else -1.0
            for ib in range(2):
                H_ap[ia, ib] = sgn * np.mean(np.exp(1j * q0s * ((tb - t_ref) % nt)) * Pi_tot[:, ia, ib])
        for a5 in (0.125, 0.0625):
            L5 = int(round(16.0 / a5))
            H_ir = np.zeros((2, 2), dtype=complex)
            for ix in range(Ls):
                k1 = 2.0 * np.pi * ix / Ls
                k1p = 2.0 * np.pi * ((ix + p1_idx) % Ls) / Ls
                Ld, _, Wd = exact_pair_spectral_weights(k1, k1p, m, L5, a5, ap=False, include_mu0=True)
                Lp, _, Wp = exact_pair_spectral_weights(k1, k1p, m, L5, a5, ap=True, include_mu0=True)
                H_ir += (np.einsum("n,nij->ij", Ld ** (ta + tb - 2), Wd)
                         - np.einsum("n,nij->ij", Lp ** (ta + tb - 2), Wp)) / Ls
            diff = float(np.max(np.abs(H_ir.real - H_ap.real)))
            print(f"  theta_s=0 Ls=4 m=0.7 {plbl:7s} a5={a5:.4f} (S=16): "
                  f"AP=[{H_ap[0,0].real:+.7f}, {H_ap[0,1].real:+.7f}, {H_ap[1,1].real:+.7f}] | "
                  f"IR=[{H_ir[0,0].real:+.7f}, {H_ir[0,1].real:+.7f}, {H_ir[1,1].real:+.7f}] | "
                  f"max|diff|={diff:.3e}")


if __name__ == "__main__":
    t0 = time.time()
    print(f"cmd: python -B wren_release10_checks.py")
    print(f"env: Python {platform.python_version()}, numpy {np.__version__}, scipy {scipy.__version__}")
    check_1_tr_Npv_scaling()
    check_2_weight_normalisation()
    check_3_m_gt_1_a5_1()
    check_4_twist_persistent_current()
    check_5_cor31_sums()
    check_6_order_of_limits_p1_ne_0()
    print(f"exit=0 elapsed={time.time() - t0:.2f}s")
