"""Independent check (Clawd) of Proposition 11: the non-product two-copy form K_np = D_t (x) D_t + D_Delta (x) D_H (label p:nonprod).

My own code; reads only data_208a.json (my file of 11-22) and the matrix Z_{Delta,B} as printed in the paper TeX.
EXACT (fractions): (ii) on the 56 pairs for D_Delta and D_H; D_H(e, e) = 0; D_H(c, c) = 1; (i) on the 256 product
events and (ii) on the 22,400 locally orthogonal pairs for K_np; LDL pivots of the c-slice.
FLOAT, a structured subfamily only (NOT the 55.7 billion spanned sets): K_np on all products S x S' and on all
two-block sets (S x T) + (S' x T^c) and (T x S) + (T^c x S') with S, S', T spanned.
Controls: the same checks with the sign of one entry of H flipped must break (ii) or D_H(e, e) = 0.
"""
import itertools
import json
import os
import re
import sys
from fractions import Fraction as Fr

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PTS = [((k >> 3) & 1, (k >> 2) & 1, (k >> 1) & 1, k & 1) for k in range(16)]  # k = 8 a0 + 4 a1 + 2 b0 + b1
MON = [(None, None), (0, None), (1, None), (None, 0), (None, 1), (0, 0), (0, 1), (1, 0), (1, 1)]
B = [[(1 if x is None else (-1) ** p[x]) * (1 if y is None else (-1) ** p[2 + y]) for (x, y) in MON] for p in PTS]
EV = [(x, a, y, b) for x in (0, 1) for a in (0, 1) for y in (0, 1) for b in (0, 1)]
IND = {e: [1 if (p[e[0]] == e[1] and p[2 + e[2]] == e[3]) else 0 for p in PTS] for e in EV}


def excl(e, f):
    return (e[0] == f[0] and e[1] != f[1]) or (e[2] == f[2] and e[3] != f[3])


def bil(D, u, v):
    return sum(u[i] * D[i][j] * v[j] for i in range(16) for j in range(16) if u[i] and v[j])


def from_walsh(Z, den):
    """B Z B^T / den as a 16 x 16 matrix of fractions"""
    BZ = [[sum(B[i][a] * Z[a][b] for a in range(9)) for b in range(9)] for i in range(16)]
    return [[Fr(sum(BZ[i][b] * B[j][b] for b in range(9)), den) for j in range(16)] for i in range(16)]


def read_Z(tex):
    lines = open(tex, encoding="utf-8").read().split("\n")
    start = next(i for i, s in enumerate(lines) if "Z_{\\Delta,B} \\;:=\\;" in s)
    rows = []
    for s in lines[start:start + 14]:
        nums = re.findall(r"-?\d+", s)
        if "&" in s and len(nums) == 9:
            rows.append([int(v) for v in nums])
    assert len(rows) == 9, len(rows)
    assert all(rows[i][j] == rows[j][i] for i in range(9) for j in range(9)), "Z not symmetric"
    return rows


def ldl_pivots(M):
    n = len(M)
    M = [row[:] for row in M]
    piv = []
    for k in range(n):
        piv.append(M[k][k])
        if M[k][k] == 0:
            return piv
        for i in range(k + 1, n):
            f = M[i][k] / M[k][k]
            for j in range(k, n):
                M[i][j] -= f * M[k][j]
    return piv


def checks(Dt, DD, DH, c, tag):
    pairs = [(e, f) for k, e in enumerate(EV) for f in EV[k + 1:] if excl(e, f)]
    assert len(pairs) == 56
    bad_D = sum(1 for (e, f) in pairs if bil(DD, IND[e], IND[f]) != 0)
    bad_H = sum(1 for (e, f) in pairs if bil(DH, IND[e], IND[f]) != 0)
    diag_H = sum(1 for e in EV if bil(DH, IND[e], IND[e]) != 0)
    print("%s: D_Delta non-zero on %d of 56 pairs ; D_H non-zero on %d of 56 pairs ; D_H(e,e) non-zero on %d of 16 ; "
          "D_H(c,c) = %s" % (tag, bad_D, bad_H, diag_H, bil(DH, c, c)))
    return bad_D == 0 and bad_H == 0 and diag_H == 0 and bil(DH, c, c) == 1


def main():
    d = json.load(open(os.path.join(HERE, "data_208a.json")))
    q = int(d["q_t"])
    Dt = [[Fr(int(v), q) for v in row] for row in d["Z_t"]]
    pt = {tuple(int(ch) for ch in k): Fr(v) for k, v in d["p_t"].items()}
    Z = read_Z(os.path.join(HERE, "..", "documentation", "paper.tex"))
    H = [[0] * 9 for _ in range(9)]
    H[6][7] = H[7][6] = 1
    H[5][8] = H[8][5] = -1
    DH = from_walsh(H, 256)
    DD = from_walsh(Z, 256 * 4096)
    C = [B[k][5] + B[k][6] + B[k][7] - B[k][8] for k in range(16)]
    c = [Fr(v, 2) for v in C]
    A = sum(1 << k for k in range(16) if C[k] == 2)
    print("A = {C = +2} has mask 0x%04X" % A)
    ok = checks(Dt, DD, DH, c, "as printed")
    assert ok
    H2 = [r[:] for r in H]
    H2[5][8] = H2[8][5] = 1
    assert not checks(Dt, DD, from_walsh(H2, 256), c, "CONTROL sign of one entry of H flipped")

    # (i) and (ii) of K_np on two copies, exact, from the one-copy tables
    T = {n: {(e, f): bil(M, IND[e], IND[f]) for e in EV for f in EV} for n, M in (("t", Dt), ("D", DD), ("H", DH))}
    bad_i = bad_ii = npairs = 0
    for e1 in EV:
        for e2 in EV:
            for f1 in EV:
                for f2 in EV:
                    val = T["t"][(e1, f1)] * T["t"][(e2, f2)] + T["D"][(e1, f1)] * T["H"][(e2, f2)]
                    if (e1, e2) == (f1, f2):
                        bad_i += val != pt[e1] * pt[e2]
                    elif excl(e1, f1) or excl(e2, f2):
                        npairs += 1
                        bad_ii += val != 0
    print("K_np two copies: (i) fails on %d of 256 events ; (ii) fails on %d of %d ordered locally orthogonal pairs"
          % (bad_i, bad_ii, npairs))
    assert bad_i == 0 and bad_ii == 0 and npairs == 44800

    # c-slice in Walsh coordinates: Q1 = D_t(c,c) B^T D_t B + D_H(c,c) B^T D_Delta B
    def walsh_mat(D):
        return [[sum(B[i][a] * D[i][j] * B[j][b] for i in range(16) for j in range(16)) for b in range(9)] for a in range(9)]
    dcc = bil(Dt, c, c)
    Wt, WD = walsh_mat(Dt), walsh_mat(DD)
    Q1 = [[dcc * Wt[a][b] + WD[a][b] for b in range(9)] for a in range(9)]
    piv = ldl_pivots(Q1)
    ev = np.linalg.eigvalsh(np.array([[float(v) for v in r] for r in Q1]) / 16.0)
    print("D_t(c,c) = %s = %.6f ; c-slice: %d exact LDL pivots, smallest %.3e, all positive %s ; smallest eigenvalue "
          "(orthonormal Walsh basis) %.4e" % (dcc, float(dcc), len(piv), float(min(piv)), all(p > 0 for p in piv), ev.min()))
    pivt = ldl_pivots([r[:] for r in Wt])
    print("for comparison Q_t alone: smallest eigenvalue %.4e"
          % np.linalg.eigvalsh(np.array([[float(v) for v in r] for r in Wt]) / 16.0).min())

    # structured subfamily, float
    Bf = np.array(B, dtype=float)
    P = Bf @ np.linalg.pinv(Bf)
    S = [np.zeros(16)]
    for m in range(1, 1 << 16):
        v = np.array([(m >> k) & 1 for k in range(16)], dtype=float)
        if np.allclose(P @ v, v, atol=1e-9):
            S.append(v)
    S = np.array(S)
    assert len(S) == 90
    f = lambda D: S @ np.array([[float(v) for v in r] for r in D]) @ S.T  # noqa: E731
    Gt, GD, GH = f(Dt), f(DD), f(DH)
    comp = [int(np.where((S == 1 - S[i]).all(axis=1))[0][0]) for i in range(90)]
    prod = np.outer(np.diag(Gt), np.diag(Gt)) + np.outer(np.diag(GD), np.diag(GH))
    nz = prod[1:, 1:]
    print("products S x S' (89 x 89 non-empty): smallest %.6e (1/67108864 = %.6e)" % (nz.min(), 1 / 67108864))
    worst = np.inf
    count = 0
    for ti in range(1, 90):
        tc = comp[ti]
        for (X1, X2) in ((Gt, Gt), (GD, GH)):
            pass
        # orientation 1: (S x T) + (S' x T^c) ; first factor carries S, S'
        v1 = (np.diag(Gt)[:, None] * Gt[ti, ti] + 2 * Gt * Gt[ti, tc] + np.diag(Gt)[None, :] * Gt[tc, tc]
              + np.diag(GD)[:, None] * GH[ti, ti] + 2 * GD * GH[ti, tc] + np.diag(GD)[None, :] * GH[tc, tc])
        # orientation 2: (T x S) + (T^c x S') ; second factor carries S, S'
        v2 = (Gt[ti, ti] * np.diag(Gt)[:, None] + 2 * Gt[ti, tc] * Gt + Gt[tc, tc] * np.diag(Gt)[None, :]
              + GD[ti, ti] * np.diag(GH)[:, None] + 2 * GD[ti, tc] * GH + GD[tc, tc] * np.diag(GH)[None, :])
        for v in (v1, v2):
            w = v.copy()
            if tc == 0:
                w = w[1:, :]          # T = Omega: the set is S x Omega, S non-empty
            else:
                w[0, 0] = np.inf      # both blocks empty
            worst = min(worst, w.min())
            count += np.isfinite(w).sum()
    print("two-block sets, both orientations: %d evaluations, smallest %.6e" % (count, worst))
    print("DONE")


if __name__ == "__main__":
    sys.exit(main())
