"""Numerical verification for Paper 6: Relativized Causality in Quantum Field Theory.

Accompanies `documentation/paper.tex` alongside Clawd's independent blind checks
`claim1_check.py` and `claim3c_check.py`.

What this script DOES and DOES NOT verify:
  - 1A: Exact generating-function remainder bound for the O(lambda^2) Luders
        expansion vs. truncated series and the tighter 5/32 * lambda^4 / (1 - lambda^2) bound.
  - 1B: Exact 4x4 2-qubit non-commuting relativized Luders channels (BOTH effect-first
        Theorem 3.1(a) and channel-first mixture Theorem 3.1(b)) vs. O(lambda^2)
        double-commutator expansion, plus a commuting zero control.
  - 1C: Step 4's alpha -> 0+ expansion of S_2^{eff}(alpha)/alpha^8 and
        S_2^{mix}(alpha)/alpha^8 in a 1+1D massive Klein-Gordon model with a
        time-translating Gaussian frame (s*m = 5.0, p_psi = 0.75), plus the
        Cauchy-Schwarz hierarchy |S_2^{eff}| <= |S_2^{mix}|.
        (Does NOT simulate finite-L cluster decomposition or spatial frame shifts.)
  - 1D: Proposition 3.3 massive no-sharpening inequality on the same 1+1D massive KG
        one-particle subspace.
  - 2A: Direct matrix-norm assertion test of Theorem 4.2(b) (n=2) and Theorem 4.2(c) (n=3)
        p_psi-weighted pairwise autocorrelation bounds on an explicit discrete causal
        qubit-chain toy model with subnormalized frame density (p_psi = 0.75).
  - 2B: Reference table (formula evaluation, not an assertion test) of Theorem 4.2(b)'s
        Gaussian autocorrelation tail bound on Theorem 3.1's truncated-double-cone geometry.
  - 3:  Appendix A (Proposition A.1) commuting control, exact closed-form 2x2 square-root
        identity (for general epsilon and balanced b = delta), and reduced-block maximal-
        decoherence family approaching column-norm squared 2 and Luders disturbance 1.
        (Full-matrix N = 300, 600 checks of every POVM hypothesis live in `claim3c_check.py`.)
"""
from __future__ import annotations

import math
import numpy as np
from scipy.linalg import sqrtm


def verify_1a_remainder_bound():
    print("=== 1A. THEOREM 3.1 (STEP 2): GENERATING-FUNCTION REMAINDER BOUND ===")
    N_ser = 120
    a = [1.0]
    for n in range(1, N_ser):
        a.append(abs(a[-1] * (0.5 - (n - 1)) / n))
    for lam in [0.05, 0.10, 0.20, 0.50, 0.80]:
        gf = 4.0 - 2.0 * (math.sqrt(1.0 - lam) + math.sqrt(1.0 + lam)) - 0.5 * lam ** 2
        s_ser = 0.0
        for n in range(N_ser):
            for m in range(N_ser):
                if (n + m) >= 4 and (n + m) % 2 == 0:
                    s_ser += a[n] * a[m] * (lam ** (n + m))
        tight_bd = (5.0 / 32.0) * (lam ** 4) / (1.0 - lam ** 2)
        diff = abs(s_ser - gf)
        assert diff < 1e-10, f"GF mismatch at lam={lam}: {diff}"
        assert tight_bd >= gf - 1e-15, f"Tight bound violated at lam={lam}"
        print(
            f"  lam={lam:.2f}: series={s_ser:.10e}, closed_GF={gf:.10e} "
            f"(diff={diff:.2e}), 5/32*lam^4/(1-lam^2)={tight_bd:.10e}"
        )


def verify_1b_2qubit_exact():
    print("\n=== 1B. THEOREM 3.1(a)-(b): 2-QUBIT EXACT RELATIVIZED LUDERS VS O(lambda^2) ===")
    I2 = np.eye(2, dtype=np.complex128)
    X = np.array([[0, 1], [1, 0]], dtype=np.complex128)
    Y = np.array([[0, -1j], [1j, 0]], dtype=np.complex128)
    Z = np.array([[1, 0], [0, -1]], dtype=np.complex128)
    rotX = lambda ang: np.cos(ang) * I2 - 1j * np.sin(ang) * X
    rotY = lambda ang: np.cos(ang) * I2 - 1j * np.sin(ang) * Y

    ys = np.array([-2.0, -1.0, 0.0, 1.0, 2.0])
    rho = np.exp(-0.5 * ys ** 2)
    rho /= rho.sum()
    p_psi = 0.80
    rho_sub = p_psi * rho

    phi_plus = np.array([1, 0, 0, 1], dtype=np.complex128) / math.sqrt(2.0)
    P_bell = np.outer(phi_plus, phi_plus.conj())

    # Use non-invariant two-qubit rotations so X_list[k] genuinely non-commutes across k
    X_list = [
        np.kron(rotY(0.25 * y), rotX(0.35 * y))
        @ P_bell
        @ np.kron(rotY(-0.25 * y), rotX(-0.35 * y))
        for y in ys
    ]
    X_psi = sum(rho_sub[k] * X_list[k] for k in range(len(ys)))

    def Lambda1(B):
        out = (1.0 - p_psi) * B
        for k in range(len(ys)):
            U1 = np.kron(rotY(0.25 * ys[k]) @ X @ rotY(-0.25 * ys[k]), I2)
            out = out + rho_sub[k] * (U1.conj().T @ B @ U1)
        return out

    B3_psi = sum(
        rho_sub[k] * np.kron(I2, rotX(0.35 * ys[k]) @ Z @ rotX(-0.35 * ys[k]))
        for k in range(len(ys))
    )
    psi0 = np.array([1, 0, 0, 0], dtype=np.complex128)

    # Zero control: Alice applies identity (or commutes with Bob & Charlie)
    B_comm = np.kron(I2, I2)
    zero_ctrl = abs(np.vdot(psi0, (Lambda1(B_comm) - B_comm) @ psi0))
    assert zero_ctrl < 1e-15
    # Non-trivial zero control: Bob uses Z x I (which commutes with I x Z), Alice uses Z x I (commutes with Bob & Charlie)
    Bob_Z = np.kron(Z, I2)
    U_comm = np.kron(Z, I2)
    sE1_c = sqrtm(0.5 * (np.eye(4) + 0.3 * Bob_Z))
    sE2_c = sqrtm(0.5 * (np.eye(4) - 0.3 * Bob_Z))
    L2_c = sE1_c @ np.kron(I2, Z) @ sE1_c + sE2_c @ np.kron(I2, Z) @ sE2_c
    sig_zero = abs(np.vdot(psi0, (U_comm.conj().T @ L2_c @ U_comm - L2_c) @ psi0))
    assert sig_zero < 1e-15
    print(f"  Zero control (commuting Alice/Bob/Charlie): |signal| = {sig_zero:.2e}")

    for lam in [0.05, 0.10, 0.20, 0.50]:
        gf_bound = 2.0 * p_psi * np.linalg.norm(B3_psi, 2) * (
            4.0 - 2.0 * (math.sqrt(1.0 - lam) + math.sqrt(1.0 + lam)) - 0.5 * lam ** 2
        )

        # (a) Effect-first relativized Luders
        E1 = 0.5 * (p_psi * np.eye(4) + lam * X_psi)
        E2 = 0.5 * (p_psi * np.eye(4) - lam * X_psi)
        E0 = (1.0 - p_psi) * np.eye(4)
        sE1, sE2, sE0 = sqrtm(E1), sqrtm(E2), sqrtm(E0)
        L2_eff = sE0 @ B3_psi @ sE0 + sE1 @ B3_psi @ sE1 + sE2 @ B3_psi @ sE2
        sig_eff = np.vdot(psi0, (Lambda1(L2_eff) - L2_eff) @ psi0).real
        comm2_eff = X_psi @ (X_psi @ B3_psi - B3_psi @ X_psi) - (X_psi @ B3_psi - B3_psi @ X_psi) @ X_psi
        L2_eff_approx = B3_psi - (lam ** 2 / (8.0 * p_psi)) * comm2_eff
        sig_eff_approx = np.vdot(psi0, (Lambda1(L2_eff_approx) - L2_eff_approx) @ psi0).real
        rem_eff = abs(sig_eff - sig_eff_approx)
        assert rem_eff <= gf_bound + 1e-14, f"Eff remainder violated: {rem_eff} > {gf_bound}"

        # (b) Channel-first (mixture) relativized Luders
        L2_mix = (1.0 - p_psi) * B3_psi
        L2_mix_approx = B3_psi.copy()
        for k in range(len(ys)):
            A1 = 0.5 * (np.eye(4) + lam * X_list[k])
            A2 = 0.5 * (np.eye(4) - lam * X_list[k])
            sA1, sA2 = sqrtm(A1), sqrtm(A2)
            L2_mix = L2_mix + rho_sub[k] * (sA1 @ B3_psi @ sA1 + sA2 @ B3_psi @ sA2)
            c1 = X_list[k] @ B3_psi - B3_psi @ X_list[k]
            c2 = X_list[k] @ c1 - c1 @ X_list[k]
            L2_mix_approx = L2_mix_approx - (lam ** 2 / 8.0) * rho_sub[k] * c2
        sig_mix = np.vdot(psi0, (Lambda1(L2_mix) - L2_mix) @ psi0).real
        sig_mix_approx = np.vdot(psi0, (Lambda1(L2_mix_approx) - L2_mix_approx) @ psi0).real
        rem_mix = abs(sig_mix - sig_mix_approx)
        assert rem_mix <= gf_bound + 1e-14, f"Mix remainder violated: {rem_mix} > {gf_bound}"

        print(
            f"  lam={lam:.2f}: eff={sig_eff:+.7f} (O(l^2)={sig_eff_approx:+.7f}, |rem|={rem_eff:.2e}) | "
            f"mix={sig_mix:+.7f} (O(l^2)={sig_mix_approx:+.7f}, |rem|={rem_mix:.2e}) <= GF={gf_bound:.2e}"
        )


def verify_1c_and_1d_massive_kg():
    print("\n=== 1C. THEOREM 3.1 (STEP 4): MASSIVE 1+1D KG, TIME-TRANSLATING FRAME (s*m=5.0) ===")
    m = 1.0
    ks = np.linspace(0.0, 2.5, 160)
    om = np.sqrt(ks ** 2 + m ** 2)
    dmu = np.exp(-ks ** 2 / 0.5)
    dmu /= np.sum(dmu * om)

    # Zero-aliasing time grid: dt = 0.5 -> Nyquist = 2*pi = 6.283 >> max(om) = 2.69
    t_pts = np.linspace(-30.0, 30.0, 121)
    dt_step = t_pts[1] - t_pts[0]
    dt_mat = t_pts[:, None] - t_pts[None, :]

    k_cos = lambda w: np.sum(w[None, None, :] * np.cos(dt_mat[:, :, None] * om[None, None, :]), axis=2)
    k_sin = lambda w: np.sum(w[None, None, :] * np.sin(dt_mat[:, :, None] * om[None, None, :]), axis=2)

    ReW_ff = k_cos(dmu)
    ReW_gg = k_cos(dmu * om ** 2)
    E_hh = 4.0 * k_sin(dmu * om ** 2)  # E(h_{y1}, h_{y2}) = 2 E(g_{y1}, g_{y2})
    ReW_gf = -k_sin(dmu * om)
    E_fg = -2.0 * k_cos(dmu * om)     # G(y - x) = E(f_x, g_y)
    E_gf = +2.0 * k_cos(dmu * om)     # -G(y - z) = E(g_y, f_z)

    s_frame = 5.0  # frame energy spread 1/s = 0.2 << m = 1.0
    p_psi = 0.75
    w_pts = p_psi * np.exp(-0.5 * (t_pts / s_frame) ** 2) / (math.sqrt(2.0 * math.pi) * s_frame) * dt_step

    M_xz = np.einsum("j,ij,jk->ik", w_pts, E_fg, E_gf)
    pred_eff = -(1.0 / (8.0 * p_psi)) * float(np.sum(np.outer(w_pts, w_pts) * (M_xz ** 2)))
    pred_mix = -(1.0 / 8.0) * float(np.einsum("i,j,k,ij,jk->", w_pts, w_pts, w_pts, E_fg ** 2, E_gf ** 2))
    assert abs(pred_eff) <= abs(pred_mix) + 1e-18, "Cauchy-Schwarz hierarchy violated!"
    print(
        f"  Analytic lim alpha->0+: eff={pred_eff:+.10e}, mix={pred_mix:+.10e} "
        f"(|eff|/|mix|={abs(pred_eff)/abs(pred_mix):.6f} <= 1)"
    )

    for alpha in [0.10, 0.05, 0.02, 0.01]:
        a2 = alpha ** 2
        tot_eff = 0.0j
        tot_mix = 0.0j
        for e1 in (+1.0, -1.0):
            for e2 in (+1.0, -1.0):
                # Conjugation W(a f_1)^* W(c) W(a f_1) = exp(+i E(a f_1, c)) W(c)
                phase_x = a2 * (e1 * E_fg[:, :, None] + e2 * E_fg[:, None, :])
                F_kick = np.tensordot(w_pts, np.exp(+1j * phase_x) - 1.0, axes=([0], [0]))
                vac_hh = np.exp(-a2 * (2.0 * ReW_gg[0, 0] + 2.0 * e1 * e2 * ReW_gg))
                F_char = np.zeros_like(F_kick)
                for e3 in (+1.0, -1.0):
                    s2 = np.sin(0.5 * a2 * e2 * e3 * E_gf)
                    s1 = np.sin(0.5 * a2 * (e1 * e2 * E_hh[:, :, None] + e1 * e3 * E_gf[:, None, :]))
                    v_jl = np.exp(-a2 * e1 * e3 * ReW_gf)
                    v_kl = np.exp(-a2 * e2 * e3 * ReW_gf)
                    v_f = np.exp(-0.5 * a2 * ReW_ff[0, 0])
                    F_char += -4.0 * v_f * np.einsum("l,kl,kl,jl,jkl->jk", w_pts, s2, v_kl, v_jl, s1)
                tot_eff += (1.0 / 8.0) * np.sum(np.outer(w_pts, w_pts) * F_kick * F_char * vac_hh)
                # Diagonal y1 = y2 = y gives the channel-first mixture amplitude (with 1/p_psi removed)
                diag_kick = np.diagonal(F_kick)
                diag_char = np.diagonal(F_char)
                diag_vac = np.diagonal(vac_hh)
                tot_mix += (1.0 / 8.0) * np.sum(w_pts * diag_kick * diag_char * diag_vac)
        S2_eff = -tot_eff / (8.0 * p_psi)
        S2_mix = -tot_mix / 8.0
        assert abs(S2_eff.imag) < 1e-15 and abs(S2_mix.imag) < 1e-15
        r_eff = (S2_eff.real / (alpha ** 8)) / pred_eff
        r_mix = (S2_mix.real / (alpha ** 8)) / pred_mix
        print(
            f"    alpha={alpha:.2f}: eff/a^8={S2_eff.real/alpha**8:+.8e} (ratio {r_eff:.6f}) | "
            f"mix/a^8={S2_mix.real/alpha**8:+.8e} (ratio {r_mix:.6f})"
        )

    # 1D: Proposition 3.3 massive no-sharpening inequality check
    print("\n=== 1D. PROPOSITION 3.3: MASSIVE NO-SHARPENING INEQUALITY CHECK ===")
    # For Gaussian time density rho_psi(t) = p_psi * N(0, s^2), rho_hat(omega) = p_psi * exp(-0.5 * s^2 * omega^2)
    rho_hat_ratio = np.exp(-0.5 * (s_frame * om) ** 2)
    r_m = float(np.max(rho_hat_ratio))
    # On one-particle vector chi = Phi(f) Omega with momentum weight dmu(k):
    # ||(rho_hat(P) - p_psi) chi|| / (p_psi ||chi||)
    lhs_ratio = math.sqrt(float(np.sum(dmu * (1.0 - rho_hat_ratio) ** 2) / np.sum(dmu)))
    rhs_bound = 1.0 - r_m
    assert lhs_ratio >= rhs_bound - 1e-15 and r_m < 1.0
    print(
        f"  s*m={s_frame:.1f}: r_psi(m) = sup_{{V_m+}} |rho_hat|/p_psi = {r_m:.3e} < 1; "
        f"||T_psi(Phi(f))Om - p Phi(f)Om|| / (p ||Phi(f)Om||) = {lhs_ratio:.10f} >= 1 - r_psi(m) = {rhs_bound:.10f}"
    )


def verify_2_discrete_psni_union_bound():
    print("\n=== 2A. THEOREM 4.2(b)-(c): DISCRETE CAUSAL-NETWORK TEST OF PAIRWISE UNION BOUND ===")
    # Exact 6-qubit (dim = 64) shift-covariant discrete causal chain:
    #   - Alice applies a non-diagonal single-site Hadamard unitary H_{s1} at site s1.
    #   - Each intermediate node Bob_j applies a strictly local two-site SWAP channel
    #     on {s_j, s_j + 1}, propagating a Pauli Z observable at most one site leftward.
    #   - Charlie measures the single-site Pauli Z_{s_last}.
    # Each node is independently channel-first relativized against a subnormalized
    # discrete shift density (p_psi = 0.75). Whenever all pairwise overlap indicators
    # vanish, Alice's Hadamard commutes with the propagated Z observable identically,
    # so the full operator norm ||Om_1(...(B)) - ...(B)|| obeys Theorem 4.2(b)-(c).
    I2 = np.eye(2)
    Z2 = np.array([[1.0, 0.0], [0.0, -1.0]])
    H2 = np.array([[1.0, 1.0], [1.0, -1.0]]) / math.sqrt(2.0)

    def site_op(op, s, n_sites):
        out = np.array([[1.0]])
        for j in range(n_sites):
            out = np.kron(out, op if j == s else I2)
        return out

    n_sites = 6
    Z_ops = [site_op(Z2, s, n_sites) for s in range(n_sites)]
    H_ops = [site_op(H2, s, n_sites) for s in range(n_sites)]
    # SWAP_{s, s+1} acting on Z_k simply permutes k: s <-> s+1
    def swap_action(s, k):
        if k == s:
            return s + 1
        if k == s + 1:
            return s
        return k

    # Test n = 2 (3 regions: Alice s1=1+x1, Bob s2={2+x2, 3+x2}, Charlie s3=4+x3)
    # Shifts x in {-1, 0, +1}, subnormalized density with p_psi = 0.75
    shifts = [-1, 0, 1]
    p_psi = 0.75
    rho = { -1: 0.15, 0: 0.45, 1: 0.15 }
    assert abs(sum(rho.values()) - p_psi) < 1e-15

    # Compute exact matrix Omega_{1,psi}(Omega_{2,psi}(B_3^psi)) - Omega_{2,psi}(B_3^psi)
    B3_psi = sum(rho[x3] * Z_ops[4 + x3] for x3 in shifts)
    # Omega_{2,psi}(B) = (1-p)B + sum_{x2} rho[x2] SWAP_{2+x2, 3+x2}(B)
    Om2_B3 = (1.0 - p_psi) * B3_psi
    for x2 in shifts:
        for x3 in shifts:
            Om2_B3 = Om2_B3 + rho[x2] * rho[x3] * Z_ops[swap_action(2 + x2, 4 + x3)]
    Om1_Om2_B3 = (1.0 - p_psi) * Om2_B3
    for x1 in shifts:
        U1 = H_ops[1 + x1]
        Om1_Om2_B3 = Om1_Om2_B3 + rho[x1] * (U1 @ Om2_B3 @ U1)
    actual_norm_n2 = float(np.linalg.norm(Om1_Om2_B3 - Om2_B3, 2))

    # Exact pairwise autocorrelation bound from Theorem 4.2(b):
    # bad_13(x1 - x3): 1 + x1 == 4 + x3  <=> x1 - x3 == 3 (0 here)
    # bad_12(x1 - x2): 1 + x1 in {2 + x2, 3 + x2} <=> x1 - x2 in {1, 2}
    # bad_23(x3 - x2): 4 + x3 in {2 + x2, 3 + x2} <=> x3 - x2 in {-2, -1}
    kap_13 = sum(rho[x1] * rho[x3] for x1 in shifts for x3 in shifts if (1 + x1) == (4 + x3))
    kap_12 = sum(rho[x1] * rho[x2] for x1 in shifts for x2 in shifts if (1 + x1) in (2 + x2, 3 + x2))
    kap_23 = sum(rho[x3] * rho[x2] for x3 in shifts for x2 in shifts if (4 + x3) in (2 + x2, 3 + x2))
    bound_n2 = 2.0 * (kap_13 + p_psi * (kap_12 + kap_23))
    assert actual_norm_n2 <= bound_n2 + 1e-14 and actual_norm_n2 > 0.0
    print(
        f"  n=2 (6-qubit exact, p_psi={p_psi}): ||Om1(Om2(B3))-Om2(B3)|| = {actual_norm_n2:.6f} "
        f"<= Thm 4.2(b) bound = {bound_n2:.6f} (ratio {actual_norm_n2/bound_n2:.4f})"
    )

    # Test n = 3 (4 regions: Alice 1+x1, Bob2 {1+x2,2+x2}, Bob3 {3+x3,4+x3}, Charlie 4+x4; shifts in {0, 1})
    shifts3 = [0, 1]
    rho3 = {0: 0.55, 1: 0.20}
    assert abs(sum(rho3.values()) - p_psi) < 1e-15
    # Apply Om3 then Om2 then Om1
    def apply_swap_psi(M_terms, base_s):
        out = [((1.0 - p_psi) * w, k) for w, k in M_terms]
        for x in shifts3:
            for w, k in M_terms:
                out.append((rho3[x] * w, swap_action(base_s + x, k)))
        return out

    terms_B4 = [(rho3[x4], 4 + x4) for x4 in shifts3]
    terms_Om3 = apply_swap_psi(terms_B4, 3)
    terms_Om23 = apply_swap_psi(terms_Om3, 1)
    Om23_B4 = sum(w * Z_ops[k] for w, k in terms_Om23)
    Om123_B4 = (1.0 - p_psi) * Om23_B4
    for x1 in shifts3:
        U1 = H_ops[1 + x1]
        Om123_B4 = Om123_B4 + rho3[x1] * (U1 @ Om23_B4 @ U1)
    actual_norm_n3 = float(np.linalg.norm(Om123_B4 - Om23_B4, 2))

    # Pairwise bad indicators for n = 3 (6 pairs: 14, 12, 13, 24, 34, 23):
    # Alice = {1+x1}, Bob2 = {1+x2,2+x2}, Bob3 = {3+x3,4+x3}, Charlie = {4+x4}
    # Note: the chain moves site 4 -> 3 -> 2 when x4=0, x3=0, x2=1, which collides with Alice at x1=1!
    k14 = sum(rho3[a] * rho3[b] for a in shifts3 for b in shifts3 if (1 + a) == (4 + b))
    k12 = sum(rho3[a] * rho3[b] for a in shifts3 for b in shifts3 if (1 + a) in (1 + b, 2 + b))
    k13 = sum(rho3[a] * rho3[b] for a in shifts3 for b in shifts3 if (1 + a) in (3 + b, 4 + b))
    k24 = sum(rho3[a] * rho3[b] for a in shifts3 for b in shifts3 if (4 + a) in (1 + b, 2 + b))
    k34 = sum(rho3[a] * rho3[b] for a in shifts3 for b in shifts3 if (4 + a) in (3 + b, 4 + b))
    k23 = sum(rho3[a] * rho3[b] for a in shifts3 for b in shifts3 if len({1 + a, 2 + a} & {3 + b, 4 + b}) > 0)
    bound_n3 = 2.0 * (k14 + p_psi * (k12 + k13 + k24 + k34) + (p_psi ** 2) * k23)
    assert actual_norm_n3 <= bound_n3 + 1e-14 and actual_norm_n3 > 1e-3
    print(
        f"  n=3 (6-qubit exact, p_psi={p_psi}): ||Om1(Om2(Om3(B4)))-Om2(Om3(B4))|| = {actual_norm_n3:.6f} "
        f"<= Thm 4.2(c) bound = {bound_n3:.6f} (ratio {actual_norm_n3/bound_n3:.4f})"
    )

    print("\n=== 2B. REFERENCE TABLE (FORMULA EVALUATION): THM 4.2(b) GAUSSIAN TAIL ON THM 3.1 GEOMETRY ===")
    s_frame = 5.0
    for L in [40.0, 80.0, 120.0, 200.0, 300.0]:
        r_perp = L / 4.0
        r_io = L / 8.0
        t_perp = (p_psi ** 2) * math.erfc(r_perp / (2.0 * s_frame))
        t_io = (p_psi ** 2) * math.erfc(r_io / (2.0 * s_frame))
        thm2_bound = 2.0 * (t_perp + 2.0 * p_psi * t_io)
        print(f"    L={L:5.1f}: Thm 4.2(b) Gaussian autocorrelation upper bound = {thm2_bound:.3e}")


def verify_3_appendix_sqrt_obstruction():
    print("\n=== 3. APPENDIX A (PROPOSITION A.1): COMMUTING CONTROL VS NON-COMMUTING OBSTRUCTION ===")
    # 3A. Commuting control
    rng = np.random.default_rng(2026)
    delta_c = 1e-4
    raw = rng.uniform(0.1, 1.0, size=(20, 50))
    Es_diag = raw / raw.sum(axis=0, keepdims=True)
    Ds_diag = np.minimum(Es_diag, delta_c * rng.uniform(0, 1, size=(20, 50)))
    Bs_diag = Es_diag - Ds_diag
    diff_sq = np.sum((np.sqrt(Es_diag) - np.sqrt(Bs_diag)) ** 2, axis=0)
    deficit = np.sum(Ds_diag, axis=0)
    assert np.all(diff_sq <= deficit + 1e-15)
    print(f"  Commuting control (M=20, d=50): max(sum(sqrt(E)-sqrt(B))^2 / sum(E-B)) = {np.max(diff_sq/deficit):.6f} <= 1.0")

    # 3B. Exact closed-form 2x2 identity for general (c, eps):
    # ||((E_i^+)^{1/2} - Btilde_i^{1/2}) e0||^2 / (b c^2) == 2 * (1 - c*(1 + sqrt(eps)) / sqrt(1 + eps*c^2 + 2*sqrt(eps)*c^2))
    print("  Exact closed-form 2x2 square-root identity (general eps and balanced b=delta):")
    for delta, c, eps in [(0.1, 0.5, 1.0 / 0.25), (0.1, 0.2, 1e-3), (1e-3, 0.05, 1e-4)]:
        s = math.sqrt(1.0 - c * c)
        b = eps * delta * c * c
        B0 = b * np.array([[c * c, c * s], [c * s, s * s]])
        E0 = B0 + delta * np.array([[0.0, 0.0], [0.0, 1.0]])
        wE, vE = np.linalg.eigh(E0)
        sE = (vE * np.sqrt(wE)) @ vE.T
        # Exact square root of rank-1 B0 = b |w><w| (avoids sqrt of O(1e-22) zero-eigenvalue roundoff)
        sB = math.sqrt(b) * np.array([[c * c, c * s], [c * s, s * s]])
        val = float(((sE - sB) @ (sE - sB))[0, 0]) / (b * c * c)
        exact = 2.0 * (1.0 - c * (1.0 + math.sqrt(eps)) / math.sqrt(1.0 + eps * c * c + 2.0 * math.sqrt(eps) * c * c))
        assert abs(val - exact) < 1e-13, f"2x2 identity mismatch: {val} vs {exact}"
        print(f"    delta={delta:.0e}, c={c:.2f}, eps={eps:.1e}: num={val:.12f}, exact={exact:.12f} (diff={abs(val-exact):.2e})")

    # 3C. Reduced-block maximal-decoherence family (see claim3c_check.py for full-matrix N=300, 600 checks)
    print("  Reduced-block maximal-decoherence family (valid POVM ||sum E_i|| < 1, sum(E_i-B_i) <= delta*I):")
    for delta, c, eps_ratio, eta in [
        (1e-2, 0.50, 1e-3, 0.02),
        (1e-2, 0.20, 1e-3, 0.02),
        (1e-3, 0.10, 1e-3, 0.01),
        (1e-4, 0.05, 1e-3, 0.01),
    ]:
        s = math.sqrt(1.0 - c * c)
        b = eps_ratio * delta * (c ** 2)
        m = int((1.0 - eta) / (b * c * c))
        a00 = m * b * c * c
        a11 = delta + b * s * s
        a01 = math.sqrt(m) * b * c * s
        SumE = np.array([[a00, a01], [a01, a11]])
        SumB = np.array([[a00, a01], [a01, b * s * s]])
        top_eig_E = float(np.linalg.eigvalsh(SumE)[-1])
        assert top_eig_E < 1.0, f"POVM normalization failed: {top_eig_E}"
        B0 = b * np.array([[c * c, c * s], [c * s, s * s]])
        E0 = B0 + delta * np.array([[0.0, 0.0], [0.0, 1.0]])
        wE, vE = np.linalg.eigh(E0)
        sE = (vE * np.sqrt(np.maximum(wE, 0.0))) @ vE.T
        sB = math.sqrt(b) * np.array([[c * c, c * s], [c * s, s * s]])
        sum_dsq = m * float(((sE - sB) @ (sE - sB))[0, 0])
        B0_tilde = np.eye(2) - SumE
        E0_plus = B0_tilde
        E0_minus = np.eye(2) - SumB
        assert float(np.linalg.eigvalsh(B0_tilde)[0]) > 0.0
        assert float(np.linalg.eigvalsh(E0_minus - B0_tilde)[0]) >= -1e-15
        w1, v1 = np.linalg.eigh(E0_plus)
        s1 = (v1 * np.sqrt(np.maximum(w1, 0.0))) @ v1.T
        w2, v2 = np.linalg.eigh(E0_minus)
        s2 = (v2 * np.sqrt(np.maximum(w2, 0.0))) @ v2.T
        rem_dsq = float(((s1 - s2) @ (s1 - s2))[0, 0])
        luders_dist = 0.5 * (sum_dsq + rem_dsq)
        print(
            f"    delta={delta:.0e}, c={c:.2f}, m={m:11d}: ||sum E_i||={top_eig_E:.5f}<1 | "
            f"<e0,+|sum D_i^2|e0,+>={sum_dsq:.6f} | Luders dist={luders_dist:.6f}"
        )


if __name__ == "__main__":
    verify_1a_remainder_bound()
    verify_1b_2qubit_exact()
    verify_1c_and_1d_massive_kg()
    verify_2_discrete_psni_union_bound()
    verify_3_appendix_sqrt_obstruction()
    print("\nALL CHECKS PASSED.")
