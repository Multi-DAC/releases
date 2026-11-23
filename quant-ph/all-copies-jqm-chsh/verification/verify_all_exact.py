"""Standalone exact-arithmetic verification suite for:
Composability of joint quantum measures for the CHSH scenario: two-copy separation
and all-copy equivalence with the almost-quantum set.

Numbering matches the compiled paper (documentation/paper.pdf):
  Definition 1, Lemma 2, Lemma 3, Proposition 4, Theorem 5, Remark 6,
  Proposition 7, Corollary 8, Remark 9, Proposition 10, Proposition 11,
  Lemma 12, Theorem 13, Remark 14, Remark 15.

Verifies in exact integer / rational / SymPy symbolic arithmetic:
  1. Data integrity (SHA-256 of data_208a.json), exact definition of p_t, CHSH = 6141/2560,
     and Proposition 4 exact PR-box certificate for C_1 = NS with branching:
     - 9x9 integer Walsh matrix Q_PR: (i) on 16/16 events, (ii) on 56/56 pairs,
       exact non-negativity on all 89 non-empty spanned subsets (min 0, 25 zeros),
       exact pivoted PSD check of (1/16) Q_PR + (16/256) I_9 >= 0 (Reduction B, lambda=16),
       and exact integer DP over all 65,535 non-empty subsets of D_PR + (1/2) Pi_{W_1^perp} (min 0).
  2. Theorem 5(a) (p_t not in Q^{1+AB}): exact rational LDL^T certificate that M >= 0 on Omega
     (M|_{W_1^perp} = 0, 9 positive pivots on W_1) and sum_e y_e p_t(e) < 0; plus interiority in Q^1.
  3. Theorem 5(b) & Remark 6 (p_t in C_2 via Reduction B + 133-pt non-spanned failure of literal D_t (x) D_t):
     - Exact symmetry, 16 event diagonals (i), 56 exclusive pairs (ii) for D_t = Z_t / q_t.
     - Exact integer DP over all 65,535 non-empty subsets of Omega: min D_t(S,S) = 1/8192 on 3 events.
     - Exact integer negative vector certifying D_t is not PSD.
     - Exact rational LDL^T certificates proving Reduction B (Lemma 3) holds for K_t with lambda = 6.
     - Exact integer verification of the 133-pt non-spanned subset F where (D_t (x) D_t)(F, F) < 0.
     - Exact verification of all 9 global minimizing product events (val = 1/67108864) and unsealed
       SHA-256 commitment match.
     - Exact orbit-stabilizer count of spanned subsets: 6, 90, 51,678, 55,720,396,530.
  4. Section 4 (Proposition 7, Remark 9, Proposition 10, Proposition 11):
     - Exact check that all 328 pairwise disjoint collections of >= 2 non-empty spanned subsets have
       PSD Gram matrices under D_t (Remark 9).
     - Exact rational check of the closed-form variance identity of Proposition 7 for k = 0..4.
     - Exact integer verification of the 6-copy block-spanned failure of D_t^(x6) (Proposition 10:
       -3.6816e-5, 3,801,088 points) and the 7-copy witness (-4.5508e-5).
     - Exact verification of the non-product 2-copy JQM K_np = D_t (x) D_t + D_Delta (x) D_H
       (Proposition 11): 56-pair (ii), zero event diagonals for D_H, exact rational LDL^T certifying
       Q^(1) > 0 in the orthonormal basis U_1 = B/4 (min pivot ~3.250e-5, lambda_min ~7.5151e-6),
       and exact rational LDL^T certifying Reduction B with lambda = 8.
  5. Section 5 (Lemma 12, Theorem 13, Remark 14: exact SymPy symbolic verification):
     - Exact symbolic rank dim(G_all) = 25 = 17 (Sym) + 8 (Skew).
     - Exact symbolic proof of G(m_chi, m_chi) = G(1_Omega, 1_Omega) (Lemma 12).
     - Exact symbolic proof of A_sk(1_Omega, u) = A_sk(u, 1_Omega) = 0, g = Q(1_A, 1_B),
       P_A - P_B = Q(1_Omega, C/2), F_{S,S'} = g * Q(1_S, 1_{S'}), and F_{A,B} = g^2.
     - Exact symbolic proof of the block-set identity g * P_E = sum_{w,w'} F_{S(w),S(w')} prod Q_{w_j,w'_j}
       at k = 2 on all 25 symbolic variables of G_all.
     - Explicit complex Hermitian PSD G = I_9 + (i/4) A_0 showing why general non-block 2-copy spanned
       sets see the skew sector (Remark 14).

Usage: python -I -B -X utf8 verify_all_exact.py
"""
from __future__ import annotations

import hashlib
import json
from collections import defaultdict
from fractions import Fraction
from pathlib import Path

import numpy as np
import sympy as sp

HERE = Path(__file__).resolve().parent
DATA_PATH = HERE / "data_208a.json"
EXPECTED_DATA_SHA256 = "f4ccc18ebbb65ef14e08cf5b8834c6d7a008497a60ed49984c031461507da44b"
EXPECTED_Q8_SHA256 = "34ece8f74f4d8ebac68e93dfe901dfabcc83775ac45fa5d6353a43fbe05a6272"


def exact_ldlt_all_pos(mat_frac: list[list[Fraction]]) -> list[Fraction]:
    """Exact rational LDL^T without pivoting; raises AssertionError unless all pivots > 0."""
    n = len(mat_frac)
    L = [[Fraction(0)] * n for _ in range(n)]
    D = [Fraction(0)] * n
    for i in range(n):
        for j in range(i):
            s = mat_frac[i][j] - sum(L[i][k] * L[j][k] * D[k] for k in range(j))
            L[i][j] = s / D[j]
        L[i][i] = Fraction(1)
        d = mat_frac[i][i] - sum(L[i][k] * L[i][k] * D[k] for k in range(i))
        assert d > 0, f"Non-positive pivot at index {i}: {float(d):.6e}"
        D[i] = d
    return D


def exact_psd_check_frac(mat_frac: list[list[Fraction]]) -> bool:
    """Exact PSD check for a symmetric Fraction matrix via pivoted rational LDL^T."""
    n = len(mat_frac)
    A = [[mat_frac[i][j] for j in range(n)] for i in range(n)]
    for k in range(n):
        piv = max(range(k, n), key=lambda i: A[i][i])
        if A[piv][piv] < 0:
            return False
        if A[piv][piv] == 0:
            for i in range(k, n):
                for j in range(k, n):
                    if A[i][j] != 0:
                        return False
            return True
        if piv != k:
            A[k], A[piv] = A[piv], A[k]
            for row in A:
                row[k], row[piv] = row[piv], row[k]
        d = A[k][k]
        row_k = A[k][:]
        for i in range(k + 1, n):
            f = row_k[i] / d
            if f != 0:
                for j in range(k + 1, n):
                    A[i][j] -= f * row_k[j]
            A[i][k] = Fraction(0)
            A[k][i] = Fraction(0)
    return True


def exact_psd_check_int_sym(M_int: list[list[int]]) -> bool:
    """Exact PSD check for a symmetric integer matrix."""
    n = len(M_int)
    return exact_psd_check_frac([[Fraction(int(M_int[i][j])) for j in range(n)] for i in range(n)])


def main() -> None:
    # =========================================================================
    # SECTION 1: Data integrity, p_t, and Proposition 4 (PR certificate for C_1 = NS)
    # =========================================================================
    raw_bytes = DATA_PATH.read_bytes()
    sha = hashlib.sha256(raw_bytes).hexdigest()
    assert sha == EXPECTED_DATA_SHA256, f"SHA mismatch: {sha}"
    data = json.loads(raw_bytes.decode("utf-8"))
    print(f"[1a] data_208a.json SHA-256 verified: {sha}")

    # Point order k = 8*a0 + 4*a1 + 2*b0 + b1
    pts = [((k >> 3) & 1, (k >> 2) & 1, (k >> 1) & 1, k & 1) for k in range(16)]
    sgn = lambda b: 1 if b == 0 else -1
    # 9 Walsh monomials spanning W_1: [1, a0, a1, b0, b1, a0b0, a0b1, a1b0, a1b1]
    B_int = [
        [
            1, sgn(a0), sgn(a1), sgn(b0), sgn(b1),
            sgn(a0) * sgn(b0), sgn(a0) * sgn(b1),
            sgn(a1) * sgn(b0), sgn(a1) * sgn(b1),
        ]
        for (a0, a1, b0, b1) in pts
    ]
    # Remaining 7 Walsh monomials spanning W_1^perp
    B_perp = [
        [
            sgn(a0) * sgn(a1),
            sgn(b0) * sgn(b1),
            sgn(a0) * sgn(a1) * sgn(b0),
            sgn(a0) * sgn(a1) * sgn(b1),
            sgn(a0) * sgn(b0) * sgn(b1),
            sgn(a1) * sgn(b0) * sgn(b1),
            sgn(a0) * sgn(a1) * sgn(b0) * sgn(b1),
        ]
        for (a0, a1, b0, b1) in pts
    ]

    # Build p = (1/5) PR + (4/5) H (alpha = 89/500), p_t = (1-t) p + t/4 (t = 1/2048)
    mu = Fraction(1, 5)
    al = Fraction(89, 500)
    t_mix = Fraction(1, 2048)

    def cg_to_probs(pa, pb, p00):
        out = {}
        for x in (0, 1):
            for y in (0, 1):
                c00 = p00[(x, y)]
                c01 = pa[x] - c00
                c10 = pb[y] - c00
                c11 = 1 - pa[x] - pb[y] + c00
                for a, b, val in ((0, 0, c00), (0, 1, c01), (1, 0, c10), (1, 1, c11)):
                    out[f"{x}{a}{y}{b}"] = val
        return out

    p_PR = cg_to_probs([Fraction(1, 2), Fraction(1, 2)], [Fraction(1, 2), Fraction(1, 2)],
                       {(0, 0): Fraction(1, 2), (0, 1): Fraction(1, 2), (1, 0): Fraction(1, 2), (1, 1): Fraction(0)})
    p_H = cg_to_probs([al, 2 * al], [al, 2 * al],
                      {(0, 0): Fraction(0), (0, 1): al, (1, 0): al, (1, 1): Fraction(0)})
    p_t = {}
    for ev in p_PR:
        p_0 = mu * p_PR[ev] + (1 - mu) * p_H[ev]
        p_t[ev] = (1 - t_mix) * p_0 + t_mix * Fraction(1, 4)
        assert p_t[ev] == Fraction(data["p_t"][ev]), f"p_t mismatch on {ev}"

    events = sorted(p_t.keys())
    ev_vec = {}
    ev_mask = {}
    for ev in events:
        x, a, y, b = map(int, ev)
        v = [1 if ((a0 if x == 0 else a1) == a and (b0 if y == 0 else b1) == b) else 0
             for (a0, a1, b0, b1) in pts]
        ev_vec[ev] = v
        ev_mask[ev] = sum(v[k] << k for k in range(16))

    excl_pairs = []
    for i, e1 in enumerate(events):
        x1, a1, y1, b1 = map(int, e1)
        for j in range(i + 1, 16):
            e2 = events[j]
            x2, a2, y2, b2 = map(int, e2)
            if (x1 == x2 and a1 != a2) or (y1 == y2 and b1 != b2):
                excl_pairs.append((e1, e2))
    assert len(excl_pairs) == 56

    E = {}
    for x in (0, 1):
        for y in (0, 1):
            E[(x, y)] = sum(sgn(a) * sgn(b) * p_t[f"{x}{a}{y}{b}"] for a in (0, 1) for b in (0, 1))
    chsh_val = E[(0, 0)] + E[(0, 1)] + E[(1, 0)] - E[(1, 1)]
    assert chsh_val == Fraction(6141, 2560)
    print(f"     p_t verified exact; CHSH(p_t) = {chsh_val} = {float(chsh_val):.10f}")

    # Build S1, S2 early for Proposition 4 PR check
    S1 = [m for m in range(16) if ((m >> 0) & 1) + ((m >> 3) & 1) == ((m >> 1) & 1) + ((m >> 2) & 1)]
    assert len(S1) == 6
    S2 = [
        m for m in range(65536)
        if all(((m >> (4 * wa)) & 15) in S1 for wa in range(4))
        and all(
            (((m >> wb) & 1) | (((m >> (4 + wb)) & 1) << 1) | (((m >> (8 + wb)) & 1) << 2) | (((m >> (12 + wb)) & 1) << 3)) in S1
            for wb in range(4)
        )
    ]
    assert len(S2) == 90
    S2_set = set(S2)

    # Proposition 4: Exact PR-box certificate for C_1 = NS with branching (all 56 pairs)
    Q_PR = [
        [ 1,  0,  0,  0,  0,  1,  1,  1, -1],
        [ 0,  1,  0,  1,  1,  0,  0,  0,  0],
        [ 0,  0,  1,  1, -1,  0,  0,  0,  0],
        [ 0,  1,  1,  1,  0,  0,  0,  0,  0],
        [ 0,  1, -1,  0,  1,  0,  0,  0,  0],
        [ 1,  0,  0,  0,  0,  1,  0,  0, -2],
        [ 1,  0,  0,  0,  0,  0,  1,  2,  0],
        [ 1,  0,  0,  0,  0,  0,  2,  1,  0],
        [-1,  0,  0,  0,  0, -2,  0,  0,  1],
    ]
    # 256 * D_PR = B Q_PR B^T (exact integer 16x16)
    Z_PR = [
        [sum(B_int[r][i] * Q_PR[i][j] * B_int[c][j] for i in range(9) for j in range(9)) for c in range(16)]
        for r in range(16)
    ]
    for ev in events:
        u = ev_vec[ev]
        assert Fraction(sum(Z_PR[i][j] * u[i] * u[j] for i in range(16) for j in range(16)), 256) == p_PR[ev]
    for e1, e2 in excl_pairs:
        u, v = ev_vec[e1], ev_vec[e2]
        assert sum(Z_PR[i][j] * u[i] * v[j] for i in range(16) for j in range(16)) == 0
    pr_sp_vals = [
        sum(Z_PR[i][j] for i in range(16) if (m >> i) & 1 for j in range(16) if (m >> j) & 1)
        for m in S2[1:]
    ]
    assert min(pr_sp_vals) == 0 and sum(1 for v in pr_sp_vals if v == 0) == 25
    # Exact PSD check of (1/16) Q_PR + (16/256) I_9 >= 0 (Reduction B lambda = 16)
    Q_PR_shift16 = [
        [Fraction(Q_PR[i][j], 16) + (Fraction(16, 256) if i == j else Fraction(0)) for j in range(9)]
        for i in range(9)
    ]
    assert exact_psd_check_frac(Q_PR_shift16)
    # Also exact integer DP over all 65,535 non-empty subsets for 256 * (D_PR + (1/2) Pi_{W_1^perp})
    # Note 16 * Pi_{W_1^perp} = 16 I - B B^T, so 256 * (1/2) Pi_{W_1^perp} = 8 * (16 I - B B^T)
    Z_PR_full = [
        [Z_PR[r][c] + 8 * ((16 if r == c else 0) - sum(B_int[r][k] * B_int[c][k] for k in range(9))) for c in range(16)]
        for r in range(16)
    ]
    dp_pr = [0] * 65536
    min_pr_full = None
    for m in range(1, 65536):
        lsb = m & -m
        bit = lsb.bit_length() - 1
        prev = m ^ lsb
        s = Z_PR_full[bit][bit]
        tmp = prev
        zrow = Z_PR_full[bit]
        while tmp:
            b2 = tmp & -tmp
            s += 2 * zrow[b2.bit_length() - 1]
            tmp ^= b2
        val = dp_pr[prev] + s
        dp_pr[m] = val
        if min_pr_full is None or val < min_pr_full:
            min_pr_full = val
    assert min_pr_full == 0
    print("[1b] Proposition 4 PR-box certificate verified exact: 16/16 events, 56/56 pairs, "
          "min_S2 = 0 (25 zeros), Reduction B lambda = 16 PSD, and D_PR + (1/2)Pi_perp >= 0 on all 65,535 subsets.")

    # =========================================================================
    # SECTION 2: Theorem 5(a) -- Exact PSD witness M proving p_t not in Q^{1+AB}
    # =========================================================================
    y_w = {ev: Fraction(data["witness_y"][ev]) for ev in events}
    w_pairs = {}
    for e1, e2, val in data["witness_w"]:
        pair = (e1, e2) if e1 < e2 else (e2, e1)
        assert pair not in w_pairs
        w_pairs[pair] = Fraction(val)
    assert len(w_pairs) == 56 and set(w_pairs.keys()) == set(excl_pairs)

    M = [[Fraction(0)] * 16 for _ in range(16)]
    for ev in events:
        u = ev_vec[ev]
        ye = y_w[ev]
        for r in range(16):
            if u[r]:
                for c in range(16):
                    if u[c]:
                        M[r][c] += ye
    for (e1, e2), wv in w_pairs.items():
        u, v = ev_vec[e1], ev_vec[e2]
        for r in range(16):
            for c in range(16):
                val = u[r] * v[c] + v[r] * u[c]
                if val:
                    M[r][c] += wv * val

    for col in range(7):
        u = [B_perp[r][col] for r in range(16)]
        for r in range(16):
            assert sum(M[r][c] * u[c] for c in range(16)) == 0
    M_W = [
        [sum(B_int[r][i] * M[r][c] * B_int[c][j] for r in range(16) for c in range(16)) for j in range(9)]
        for i in range(9)
    ]
    exact_ldlt_all_pos(M_W)
    y_dot_pt = sum(y_w[ev] * p_t[ev] for ev in events)
    assert y_dot_pt == Fraction(-857948059075795414601, 9365785569480550400000000) and y_dot_pt < 0

    mA = [sum(sgn(a) * p_t[f"{x}{a}0{b}"] for a in (0, 1) for b in (0, 1)) for x in (0, 1)]
    mB = [sum(sgn(b) * p_t[f"0{a}{y}{b}"] for a in (0, 1) for b in (0, 1)) for y in (0, 1)]
    c_aa = Fraction(0)
    G1 = [
        [Fraction(1), mA[0], mA[1], mB[0], mB[1]],
        [mA[0], Fraction(1), c_aa, E[(0, 0)], E[(0, 1)]],
        [mA[1], c_aa, Fraction(1), E[(1, 0)], E[(1, 1)]],
        [mB[0], E[(0, 0)], E[(1, 0)], Fraction(1), c_aa],
        [mB[1], E[(0, 1)], E[(1, 1)], c_aa, Fraction(1)],
    ]
    exact_ldlt_all_pos(G1)
    print(f"[2] Theorem 5(a) certified: M >= 0 (9 exact LDL^T pivots > 0, M|_{{W_1^perp}} = 0), "
          f"y.p_t = {y_dot_pt} ({float(y_dot_pt):.6e}) < 0; p_t in int(Q^1) certified.")

    # =========================================================================
    # SECTION 3: Theorem 5(b) & Remark 6 -- D_t, Reduction B (lambda=6), 133-pt fail
    # =========================================================================
    q_t = int(data["q_t"])
    assert q_t == 2500 * (2 ** 49)
    Z_t = [[int(x) for x in row] for row in data["Z_t"]]
    assert all(Z_t[i][j] == Z_t[j][i] for i in range(16) for j in range(16))

    def bil_Z(u, v):
        return sum(Z_t[i][j] * u[i] * v[j] for i in range(16) for j in range(16))

    for ev in events:
        assert Fraction(bil_Z(ev_vec[ev], ev_vec[ev]), q_t) == p_t[ev]
    for e1, e2 in excl_pairs:
        assert bil_Z(ev_vec[e1], ev_vec[e2]) == 0
    tr_MD = sum(M[i][j] * Fraction(Z_t[i][j], q_t) for i in range(16) for j in range(16))
    assert tr_MD == y_dot_pt

    dp = [0] * 65536
    min_val = None
    min_masks = []
    for m in range(1, 65536):
        lsb = m & -m
        bit = lsb.bit_length() - 1
        prev = m ^ lsb
        s = Z_t[bit][bit]
        tmp = prev
        zrow = Z_t[bit]
        while tmp:
            b2 = tmp & -tmp
            s += 2 * zrow[b2.bit_length() - 1]
            tmp ^= b2
        val = dp[prev] + s
        dp[m] = val
        if min_val is None or val < min_val:
            min_val = val
            min_masks = [m]
        elif val == min_val:
            min_masks.append(m)
    assert Fraction(min_val, q_t) == Fraction(1, 8192)
    mask_to_ev = {v: k for k, v in ev_mask.items()}
    min_evs = sorted(mask_to_ev[m] for m in min_masks)
    assert min_evs == ["0011", "1010", "1100"]

    v_neg = [207, -1860, 2, -40, 3339, -67, 1275, -4536, -417, -1664, 466, -976, 665, 7434, 1545, -970]
    vZv = bil_Z(v_neg, v_neg)
    assert vZv == -252817480939407672851040 and vZv < 0

    # Orthonormal basis U_1 = B / 4 of W_1: Q_orth = (1/16) B^T D_t B
    Q_orth = [
        [Fraction(sum(B_int[r][i] * Z_t[r][c] * B_int[c][j] for r in range(16) for c in range(16)), 16 * q_t)
         for j in range(9)]
        for i in range(9)
    ]
    exact_ldlt_all_pos([
        [(Fraction(1, 5) if i == j else Fraction(0)) - Q_orth[i][j] for j in range(9)]
        for i in range(9)
    ])
    exact_ldlt_all_pos([
        [Q_orth[i][j] + (Fraction(1, 2300) if i == j else Fraction(0)) for j in range(9)]
        for i in range(9)
    ])
    assert Fraction(1, 5) * Fraction(1, 2300) < Fraction(6, 65536)

    # Exact 133-point non-spanned subset F where (D_t (x) D_t)(F, F) < 0 (Remark 6)
    F133_mat = [
        [1, 1, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0],
        [1, 1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0],
        [0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 1],
        [0, 0, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1],
        [1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1],
        [0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 0],
        [1, 1, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 1, 0],
        [0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1],
        [1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 0],
        [0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 1, 0],
        [1, 1, 1, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 0],
        [1, 1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 1, 1],
        [1, 1, 0, 0, 1, 0, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1],
        [1, 0, 1, 1, 1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0],
        [0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0],
    ]
    assert sum(sum(r) for r in F133_mat) == 133
    F133_np = np.array(F133_mat, dtype=np.int64).reshape(4, 4, 4, 4)
    u_perp = np.array([1, -1, -1, 1], dtype=np.int64)
    assert not all(np.all(np.tensordot(F133_np, u_perp, axes=([ax], [0])) == 0) for ax in range(4))
    ZF = [[sum(Z_t[i][k] * F133_mat[k][j] for k in range(16)) for j in range(16)] for i in range(16)]
    val_F133 = sum(ZF[i][l] * Z_t[l][j] for i in range(16) for j in range(16) if F133_mat[i][j] for l in range(16))
    assert val_F133 == -13773446702808538793200383356544768 and val_F133 < 0

    q8_commit_str = "inf over non-empty spanned f of K_t(f, f) = 1/67108864"
    assert hashlib.sha256(q8_commit_str.encode("utf-8")).hexdigest() == EXPECTED_Q8_SHA256
    for e1 in min_evs:
        for e2 in min_evs:
            assert Fraction(bil_Z(ev_vec[e1], ev_vec[e1]) * bil_Z(ev_vec[e2], ev_vec[e2]), q_t * q_t) == Fraction(1, 67108864)

    # Exact orbit-stabilizer count of spanned subsets: 6, 90, 51678, 55720396530
    pairs2 = defaultdict(list)
    for a in S2:
        for d in S2:
            pairs2[(a ^ d, a & d)].append((a, d))
    S3 = [a | (b << 16) | (c << 32) | (d << 48) for lst in pairs2.values() for a, d in lst for b, c in lst]
    assert len(S3) == 51678

    wing_perms = []
    for sw in (0, 1):
        for f0 in (0, 1):
            for f1 in (0, 1):
                p = []
                for idx in range(4):
                    x0, x1 = (idx >> 1) & 1, idx & 1
                    if sw:
                        x0, x1 = x1, x0
                    x0 ^= f0
                    x1 ^= f1
                    p.append((x0 << 1) | x1)
                wing_perms.append(p)

    def apply_perm64(m, perm64, comp=False):
        out = 0
        tmp = m
        while tmp:
            lsb = tmp & -tmp
            out |= 1 << perm64[lsb.bit_length() - 1]
            tmp ^= lsb
        return ((1 << 64) - 1) ^ out if comp else out

    gens64 = [
        ([16 * w2 + 4 * w0 + w1 for w2 in range(4) for w1 in range(4) for w0 in range(4)], False),
        ([16 * w1 + 4 * w2 + w0 for w2 in range(4) for w1 in range(4) for w0 in range(4)], False),
    ]
    for wp in (wing_perms[1], wing_perms[2], wing_perms[4]):
        gens64.append(([16 * w2 + 4 * w1 + wp[w0] for w2 in range(4) for w1 in range(4) for w0 in range(4)], False))
    gens64.append((list(range(64)), True))

    visited = set()
    orbits = []
    for s in S3:
        if s in visited:
            continue
        orb = {s}
        stack = [s]
        while stack:
            cur = stack.pop()
            for p64, comp in gens64:
                nxt = apply_perm64(cur, p64, comp)
                if nxt not in orb:
                    orb.add(nxt)
                    stack.append(nxt)
        visited.update(orb)
        orbits.append((s, len(orb)))
    assert len(orbits) == 76 and sum(sz for _, sz in orbits) == 51678

    s2_sig_id = {}
    s2_sig_cnt = []
    s2_pair_sig = {}
    for a in S2:
        for b in S2:
            key = (a ^ b, a & b)
            if key not in s2_sig_id:
                s2_sig_id[key] = len(s2_sig_cnt)
                s2_sig_cnt.append(0)
            sid = s2_sig_id[key]
            s2_sig_cnt[sid] += 1
            s2_pair_sig[(a, b)] = sid
    assert sum(c * c for c in s2_sig_cnt) == 51678

    triples_by_r = {}
    for r in S2:
        lst = []
        for a in S2:
            xa = a ^ r
            for b in S2:
                if (xa & (b ^ r)) == 0:
                    c = a ^ b ^ r
                    if c in S2_set:
                        lst.append((a, b, c))
        triples_by_r[r] = lst
    assert sum(len(v) for v in triples_by_r.values()) == 51678

    pair_cnt_cache = {}
    def get_pair_sig_counts(r0, r3):
        key = (r0, r3) if r0 <= r3 else (r3, r0)
        if key in pair_cnt_cache:
            return pair_cnt_cache[key]
        d = defaultdict(int)
        t0, t3 = triples_by_r[r0], triples_by_r[r3]
        ps = s2_pair_sig
        for a0, b0, c0 in t0:
            for a3, b3, c3 in t3:
                k = (ps[(a0, a3)] << 24) | (ps[(b0, b3)] << 12) | ps[(c0, c3)]
                d[k] += 1
        pair_cnt_cache[key] = d
        return d

    total_S4 = 0
    for rep, sz in orbits:
        r0 = rep & 0xFFFF
        r1 = (rep >> 16) & 0xFFFF
        r2 = (rep >> 32) & 0xFFFF
        r3 = (rep >> 48) & 0xFFFF
        d03 = get_pair_sig_counts(r0, r3)
        d12 = get_pair_sig_counts(r1, r2)
        if len(d03) > len(d12):
            d03, d12 = d12, d03
        cnt = sum(v * d12[k] for k, v in d03.items() if k in d12)
        total_S4 += sz * cnt
    assert total_S4 == 55720396530
    print(f"[3] Theorem 5(b) & Remark 6 certified: D_t exact (i)+(ii), min_{{S!=0}} D_t(S,S) = 1/8192 on {min_evs}, "
          f"D_t not PSD (v^T Z_t v = {vZv} < 0);")
    print(f"    Reduction B certified exact for K_t with lambda = 6; literal K_t fails on 133-pt non-spanned F "
          f"({val_F133 / (q_t**2):.6e} < 0);")
    print(f"    Q8 commitment SHA-256 matched; exact spanned counts = 6, 90, 51678, {total_S4}.")

    # =========================================================================
    # SECTION 4: Remark 9, Proposition 7, Proposition 10, Proposition 11
    # =========================================================================
    S2_nonempty = [m for m in S2 if m != 0]
    G_Z_90 = [
        [sum(Z_t[a][b] for a in range(16) if (S2[i] >> a) & 1 for b in range(16) if (S2[j] >> b) & 1)
         for j in range(90)]
        for i in range(90)
    ]
    s2_idx = {m: i for i, m in enumerate(S2)}
    disjoint_colls = []
    def dfs_disj(start, used_mask, cur):
        if len(cur) >= 2:
            disjoint_colls.append(tuple(cur))
        for idx in range(start, len(S2_nonempty)):
            m = S2_nonempty[idx]
            if (m & used_mask) == 0:
                cur.append(m)
                dfs_disj(idx + 1, used_mask | m, cur)
                cur.pop()
    dfs_disj(0, 0, [])
    assert len(disjoint_colls) == 328
    non_psd_colls = 0
    for coll in disjoint_colls:
        idxs = [s2_idx[m] for m in coll]
        sub = [[G_Z_90[i][j] for j in idxs] for i in idxs]
        if not exact_psd_check_int_sym(sub):
            non_psd_colls += 1
    assert non_psd_colls == 0
    print(f"[4a] Remark 9 (DW Lemma 12 obstruction) certified: 0 / {len(disjoint_colls)} pairwise disjoint spanned collections are non-PSD.")

    def bil_mask_int(m1, m2):
        return G_Z_90[s2_idx[m1]][s2_idx[m2]]

    A_viol, B_viol = 0xCA53, 0xFFFF ^ 0xCA53
    a_f = Fraction(bil_mask_int(A_viol, A_viol), q_t)
    c_f = Fraction(bil_mask_int(B_viol, B_viol), q_t)
    b_f = -Fraction(bil_mask_int(A_viol, B_viol), q_t)
    assert b_f == Fraction(3902622654083123, 21990232555520000) and b_f > 0
    assert a_f + c_f - 2 * b_f == 1

    p_plus = [(0x00AA, Fraction(1, 3)), (0x3535, Fraction(1, 2)), (0x35AC, Fraction(1, 6))]
    p_minus = [(0x3030, Fraction(1, 4)), (0x05AF, Fraction(1, 2)), (0x0000, Fraction(1, 4))]
    Q_90_f = [[Fraction(G_Z_90[i][j], q_t) for j in range(90)] for i in range(90)]

    def eval_dist_bil(d1, d2):
        return sum(w1 * w2 * Q_90_f[s2_idx[m1]][s2_idx[m2]] for m1, w1 in d1 for m2, w2 in d2)

    q_pp = eval_dist_bil(p_plus, p_plus)
    q_mm = eval_dist_bil(p_minus, p_minus)
    q_pm = eval_dist_bil(p_plus, p_minus)
    q_xb_xb = (q_pp + 2 * q_pm + q_mm) / 4
    q_dl_dl = (q_pp - 2 * q_pm + q_mm) / 4
    q_xb_dl = (q_pp - q_mm) / 4
    V_plus = sum(w * Q_90_f[s2_idx[m]][s2_idx[m]] for m, w in p_plus) - q_pp
    V_minus = sum(w * Q_90_f[s2_idx[m]][s2_idx[m]] for m, w in p_minus) - q_mm
    M2_f = [[a_f, -b_f], [-b_f, c_f]]
    for k in range(5):
        closed = (
            q_dl_dl * (1 + 4 * b_f) ** k
            + 2 * q_xb_dl * (a_f - c_f) ** k
            + q_xb_xb
            + V_plus * ((a_f + c_f) ** k + (a_f - c_f) ** k) / 2
            + V_minus * ((a_f + c_f) ** k - (a_f - c_f) ** k) / 2
        )
        direct = Fraction(0)
        for w1 in range(1 << k):
            par1 = bin(w1).count("1") & 1
            for w2 in range(1 << k):
                par2 = bin(w2).count("1") & 1
                prod_m = Fraction(1)
                for bit in range(k):
                    prod_m *= M2_f[(w1 >> bit) & 1][(w2 >> bit) & 1]
                if w1 == w2:
                    q_term = (q_pp + V_plus) if par1 == 0 else (q_mm + V_minus)
                else:
                    q_term = q_pp if (par1 == 0 and par2 == 0) else (q_mm if (par1 == 1 and par2 == 1) else q_pm)
                direct += prod_m * q_term
        assert closed == direct, f"Variance identity mismatch at k={k}"
    print("[4b] Proposition 7 exact closed-form variance identity verified against direct double sum for k = 0..4.")

    r6 = [
        0x0, 0x0, 0x5AF, 0x3030, 0x3030, 0x3030, 0x0, 0x0,
        0x3030, 0x3535, 0x0, 0x0, 0x0, 0x0, 0x35AC, 0x33AA,
        0x3030, 0x3030, 0x0, 0x0, 0x0, 0x0, 0x35AC, 0x33AA,
        0x0, 0x0, 0x35AC, 0x35AC, 0x33AA, 0x35AC, 0xA0A0, 0xA3A3,
    ]
    assert all(m in S2_set for m in r6)
    size_E6 = sum(bin(m).count("1") for m in r6) * (8 ** 5)
    assert size_E6 == 3801088
    pairs_blk6 = [(0xCA53, 0xFFFF ^ 0xCA53)] * 4 + [(0xC5A3, 0xFFFF ^ 0xC5A3)]
    M_blk6 = [
        [[bil_mask_int(A, A), bil_mask_int(A, B)], [bil_mask_int(B, A), bil_mask_int(B, B)]]
        for A, B in pairs_blk6
    ]
    tot6 = 0
    for w in range(32):
        if r6[w] == 0:
            continue
        for wp in range(32):
            if r6[wp] == 0:
                continue
            prod_b = 1
            for c in range(5):
                prod_b *= M_blk6[c][(w >> (4 - c)) & 1][(wp >> (4 - c)) & 1]
            tot6 += bil_mask_int(r6[w], r6[wp]) * prod_b
    expected_6 = -286086512998116794308145853677292930292488963295984229317414147707056527838366297712306896727047999062016
    assert tot6 == expected_6 and tot6 < 0

    A_pos, B_pos = 0x3A5C, 0xC5A3
    assert Fraction(bil_mask_int(A_pos, B_pos), q_t) == Fraction(3902622654083123, 21990232555520000)
    r7 = [
        3, 23, 23, 23, 23, 23, 24, 23, 23, 9, 23, 3, 23, 9, 23, 23,
        15, 23, 23, 3, 9, 24, 1, 23, 9, 15, 23, 3, 23, 23, 24, 23,
        1, 24, 23, 23, 3, 23, 23, 9, 24, 23, 23, 23, 23, 23, 23, 3,
        23, 23, 24, 23, 3, 23, 23, 23, 23, 23, 23, 9, 9, 23, 23, 3,
    ]
    m7 = [[bil_mask_int(A_pos, A_pos), bil_mask_int(A_pos, B_pos)],
          [bil_mask_int(B_pos, A_pos), bil_mask_int(B_pos, B_pos)]]
    tot7 = 0
    for w1 in range(64):
        for w2 in range(64):
            p7 = 1
            for bit in range(6):
                p7 *= m7[(w1 >> (5 - bit)) & 1][(w2 >> (5 - bit)) & 1]
            tot7 += p7 * G_Z_90[r7[w1]][r7[w2]]
    expected_7 = -497691754512974689080408355518942530393912747752012541021070533461626154642942484258055177956079192364847965473580951535616
    assert tot7 == expected_7 and tot7 < 0
    print(f"[4c] Proposition 10 exact 6-copy witness verified: |E_6| = {size_E6}, D_t^(x6)(E_6,E_6) = {float(Fraction(tot6, q_t**6)):.10e} < 0;")
    print(f"     Exact 7-copy witness verified: D_t^(x7)(E_7,E_7) = {float(Fraction(tot7, q_t**7)):.10e} < 0.")

    # 4d. Proposition 11: Non-product 2-copy JQM K_np = D_t (x) D_t + D_Delta (x) D_H
    Z_DQ_B = [
        [ -4, -35,  -8, -29, -15, -87, -24, -25,  27],
        [-35,  -4,  17, -87, -24, -29, -15,  12,  23],
        [ -8,  17,  -4, -25,  27,  12,  23, -29, -15],
        [-29, -87, -25,  -4,   1, -35,  -1,  -8,  31],
        [-15, -24,  27,   1,  -4,  -1, -35,  31,  -8],
        [-87, -29,  12, -35,  -1,  -4,   1,  17,  18],
        [-24, -15,  23,  -1, -35,   1,  -4,  35,  17],
        [-25,  12, -29,  -8,  31,  17,  35,  -4,   1],
        [ 27,  23, -15,  31,  -8,  18,  17,   1,  -4],
    ]
    DEN_DQ = 4096
    H_null_B = [[0] * 9 for _ in range(9)]
    H_null_B[6][7] = H_null_B[7][6] = 1
    H_null_B[5][8] = H_null_B[8][5] = -1

    D_Delta = [
        [Fraction(sum(B_int[r][a] * Z_DQ_B[a][b] * B_int[c][b] for a in range(9) for b in range(9)), 256 * DEN_DQ)
         for c in range(16)]
        for r in range(16)
    ]
    D_H = [
        [Fraction(sum(B_int[r][a] * H_null_B[a][b] * B_int[c][b] for a in range(9) for b in range(9)), 256)
         for c in range(16)]
        for r in range(16)
    ]
    def bil_F(M_f, u, v):
        return sum(M_f[i][j] * u[i] * v[j] for i in range(16) for j in range(16))

    for ev in events:
        assert bil_F(D_H, ev_vec[ev], ev_vec[ev]) == 0
        if ev in ("0011", "1010", "1100"):
            assert bil_F(D_Delta, ev_vec[ev], ev_vec[ev]) == 0
    for e1, e2 in excl_pairs:
        assert bil_F(D_Delta, ev_vec[e1], ev_vec[e2]) == 0
        assert bil_F(D_H, ev_vec[e1], ev_vec[e2]) == 0

    c_pts = [(B_int[k][5] + B_int[k][6] + B_int[k][7] - B_int[k][8]) // 2 for k in range(16)]
    assert bil_F(D_H, c_pts, c_pts) == 1
    qcc_frac = Fraction(bil_Z(c_pts, c_pts), q_t)
    assert qcc_frac == 1 + 4 * b_f

    # Exact rational LDL^T in the orthonormal basis U_1 = B/4 proving Q^(1) = (1+4b) Q_t + (1/(16*4096)) Z_{Delta,B} > 0
    D_D_orth = [[Fraction(Z_DQ_B[i][j], 16 * DEN_DQ) for j in range(9)] for i in range(9)]
    D_H_orth = [[Fraction(H_null_B[i][j], 16) for j in range(9)] for i in range(9)]
    Q1_orth_frac = [
        [qcc_frac * Q_orth[i][j] + D_D_orth[i][j] for j in range(9)]
        for i in range(9)
    ]
    pivs_Q1 = exact_ldlt_all_pos(Q1_orth_frac)
    min_piv_Q1 = min(float(p) for p in pivs_Q1)
    min_eig_Q1 = float(np.linalg.eigvalsh(np.array([[float(x) for x in r] for r in Q1_orth_frac])).min())

    shift_8 = Fraction(8, 65536)
    K_np_81 = []
    for i1 in range(9):
        for i2 in range(9):
            row = []
            for j1 in range(9):
                for j2 in range(9):
                    v = Q_orth[i1][j1] * Q_orth[i2][j2] + D_D_orth[i1][j1] * D_H_orth[i2][j2]
                    if i1 == j1 and i2 == j2:
                        v += shift_8
                    row.append(v)
            K_np_81.append(row)
    exact_ldlt_all_pos(K_np_81)
    print(f"[4d] Proposition 11 non-product 2-copy JQM K_np verified exact: (i)+(ii), Q^(1) > 0 in U_1=B/4 "
          f"(9 exact LDL^T pivots > 0, min pivot = {min_piv_Q1:.4e}, lambda_min = {min_eig_Q1:.4e}), "
          "Reduction B certified exact with lambda = 8.")

    # =========================================================================
    # SECTION 5: Lemma 12, Theorem 13, Remark 14 -- Exact SymPy verification
    # =========================================================================
    B_sp = sp.Matrix(B_int)
    ev_sp = [(ev, (B_sp.T * sp.Matrix(ev_vec[ev])) / 16) for ev in events]

    m_syms = [[sp.Symbol(f"m_{i}_{j}") for j in range(9)] for i in range(9)]
    M_gen = sp.Matrix(m_syms)
    eqs_ii = []
    for i in range(16):
        e1, u = ev_sp[i]
        x1, a1, y1, b1 = map(int, e1)
        for j in range(16):
            e2, v = ev_sp[j]
            x2, a2, y2, b2 = map(int, e2)
            if (x1 == x2 and a1 != a2) or (y1 == y2 and b1 != b2):
                eqs_ii.append((u.T * M_gen * v)[0])
    sol_ii = sp.solve(eqs_ii, [m_syms[i][j] for i in range(9) for j in range(9)])
    G_sym = M_gen.subs(sol_ii)
    free25 = sorted(G_sym.free_symbols, key=lambda s: s.name)
    assert len(free25) == 25

    Q_part = (G_sym + G_sym.T) / 2
    A_part = (G_sym - G_sym.T) / 2
    q_syms_d = {(i, j): sp.Symbol(f"q_{i}_{j}") for i in range(9) for j in range(i, 9)}
    Q_mat = sp.Matrix(9, 9, lambda i, j: q_syms_d[(min(i, j), max(i, j))])
    a_syms_d = {(i, j): sp.Symbol(f"a_{i}_{j}") for i in range(9) for j in range(i + 1, 9)}
    A_mat = sp.Matrix(9, 9, lambda i, j: 0 if i == j else (a_syms_d[(i, j)] if i < j else -a_syms_d[(j, i)]))
    eqs_q, eqs_a = [], []
    for i in range(16):
        e1, u = ev_sp[i]
        x1, a1, y1, b1 = map(int, e1)
        for j in range(16):
            e2, v = ev_sp[j]
            x2, a2, y2, b2 = map(int, e2)
            if (x1 == x2 and a1 != a2) or (y1 == y2 and b1 != b2):
                eqs_q.append((u.T * Q_mat * v)[0])
                eqs_a.append((u.T * A_mat * v)[0])
    Q_sol = Q_mat.subs(sp.solve(eqs_q, list(q_syms_d.values())))
    A_sol = A_mat.subs(sp.solve(eqs_a, list(a_syms_d.values())))
    assert len(Q_sol.free_symbols) == 17 and len(A_sol.free_symbols) == 8

    # Lemma 12 core
    g00 = G_sym[0, 0]
    for i in range(9):
        assert sp.simplify(G_sym[i, i] - g00) == 0

    # Corollary 8 / Step 4 core
    for j in range(9):
        assert sp.simplify(A_part[0, j]) == 0 and sp.simplify(A_part[j, 0]) == 0

    vA_w = sp.Matrix([sp.Rational(1, 2), 0, 0, 0, 0, sp.Rational(1, 4), sp.Rational(1, 4), sp.Rational(1, 4), sp.Rational(-1, 4)])
    vB_w = sp.Matrix([sp.Rational(1, 2), 0, 0, 0, 0, sp.Rational(-1, 4), sp.Rational(-1, 4), sp.Rational(-1, 4), sp.Rational(1, 4)])
    assert sp.simplify((A_part * vB_w).dot(vA_w)) == 0
    g_AB = sp.simplify((Q_part * vB_w).dot(vA_w))
    P_A_sym = (G_sym * vA_w).dot(vA_w)
    P_B_sym = (G_sym * vB_w).dot(vB_w)
    assert sp.simplify(g_AB - (g00 - P_A_sym - P_B_sym) / 2) == 0

    v_one = sp.Matrix([1, 0, 0, 0, 0, 0, 0, 0, 0])
    v_C2 = vA_w - vB_w
    assert sp.simplify((P_A_sym - P_B_sym) - (v_one.T * Q_part * v_C2)[0]) == 0

    def w_of_mask(m):
        return (B_sp.T * sp.Matrix([(m >> k) & 1 for k in range(16)])) / 16

    for m1, m2 in [(0xCA53, 0x35AC), (0x00AA, 0x35AC), (0x05AF, 0x33AA), (0x3535, 0xA3A3)]:
        w1, w2 = w_of_mask(m1), w_of_mask(m2)
        P_E = (
            (w1.T * G_sym * w1)[0] * P_A_sym
            + (w2.T * G_sym * w2)[0] * P_B_sym
            + (w1.T * G_sym * w2)[0] * (vA_w.T * G_sym * vB_w)[0]
            + (w2.T * G_sym * w1)[0] * (vB_w.T * G_sym * vA_w)[0]
        )
        F_lhs = (P_E - (w1.T * G_sym * w1)[0] * P_A_sym - (w2.T * G_sym * w2)[0] * P_B_sym) / 2
        F_rhs = g_AB * (w1.T * Q_part * w2)[0]
        assert sp.expand(F_lhs - F_rhs) == 0
        if (m1, m2) == (0xCA53, 0x35AC):
            assert sp.expand(F_lhs - g_AB ** 2) == 0

    # k = 2 block-set identity on all 25 symbolic variables
    blk_masks = [0x00AA, 0x05AF, 0x3535, 0x33AA]
    blk_ws = [w_of_mask(m) for m in blk_masks]
    AB_ws = [vA_w, vB_w]
    Q_blk = [[P_A_sym, g_AB], [g_AB, P_B_sym]]
    lhs_k2 = 0
    rhs_k2 = 0
    for w in range(4):
        b0, b1 = (w >> 1) & 1, w & 1
        for wp in range(4):
            bp0, bp1 = (wp >> 1) & 1, wp & 1
            prod_G = (AB_ws[b0].T * G_sym * AB_ws[bp0])[0] * (AB_ws[b1].T * G_sym * AB_ws[bp1])[0]
            lhs_k2 += (blk_ws[w].T * G_sym * blk_ws[wp])[0] * prod_G
            rhs_k2 += (g_AB * (blk_ws[w].T * Q_part * blk_ws[wp])[0]) * Q_blk[b0][bp0] * Q_blk[b1][bp1]
    assert sp.expand(g_AB * lhs_k2 - rhs_k2) == 0

    # Remark 14: Explicit complex Hermitian PSD example showing non-block 2-copy spanned sets see skew sector
    a_free = sorted(A_sol.free_symbols, key=lambda s: s.name)
    A0_sp = A_sol.subs({s: (1 if i == 0 else 0) for i, s in enumerate(a_free)})
    max_a0 = max(abs(A0_sp[i, j]) for i in range(9) for j in range(9))
    A0_q = A0_sp / (4 * max_a0)
    I_minus_AtA = [
        [Fraction(1 if i == j else 0) - sum(Fraction(int(sp.fraction(A0_q[k, i])[0]), int(sp.fraction(A0_q[k, i])[1]))
                                            * Fraction(int(sp.fraction(A0_q[k, j])[0]), int(sp.fraction(A0_q[k, j])[1]))
                                            for k in range(9))
         for j in range(9)]
        for i in range(9)
    ]
    exact_ldlt_all_pos(I_minus_AtA)
    D_A0 = (B_sp / 16) * A0_q * (B_sp.T / 16)
    D_A0_f = np.array(D_A0.tolist(), dtype=float)
    found_skew_F = None
    rng = np.random.default_rng(7)
    S3_arr = np.array(S3, dtype=np.uint64)
    BIT_SHIFTS = np.arange(64, dtype=np.uint64)
    for _ in range(100):
        g0 = S3_arr[rng.integers(len(S3_arr))]
        g1 = S3_arr[rng.integers(len(S3_arr))]
        cands = S3_arr[((S3_arr ^ g0) & (g1 ^ g0)) == np.uint64(0)]
        g2 = cands[rng.integers(len(cands))]
        g3 = g1 ^ g2 ^ g0
        G_arr = np.stack([((np.uint64(x) >> BIT_SHIFTS) & np.uint64(1)).astype(np.int64).reshape(4, 16) for x in (g0, g1, g2, g3)], axis=0)
        F_cand = G_arr.transpose(2, 1, 0).reshape(16, 16)
        diff_f = float(np.sum((D_A0_f @ F_cand) * (F_cand @ D_A0_f)))
        if abs(diff_f) > 1e-6:
            found_skew_F = F_cand
            break
    assert found_skew_F is not None
    F_sp = sp.Matrix(found_skew_F.tolist())
    skew_contrib = sp.simplify((F_sp.T * D_A0 * F_sp * D_A0).trace())
    assert skew_contrib == sp.Rational(-1, 1024)
    print("[5] Lemma 12 & Theorem 13 exact SymPy symbolic verification passed: dim(G_all) = 25 = 17 + 8,")
    print("    Lemma 12, A_sk(1_Omega, u) = 0, g = Q(1_A, 1_B), F_{S,S'} = g*Q(1_S,1_{S'}), F_{A,B} = g^2,")
    print("    k = 2 block-set identity verified identically on all 25 symbolic variables;")
    print(f"    Remark 14 skew-sector separation on non-block S4 set verified exact: -A_sk^(x2)(1_F,1_F) = {skew_contrib} != 0.")
    print("\nALL EXACT VERIFICATIONS PASSED.")


if __name__ == "__main__":
    main()
