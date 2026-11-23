"""Generate binary matrix data file nonprod_data.bin for VerifyNonProd.cs
certifying K_np = D_t (x) D_t + D_Delta (x) D_H over all 55,720,396,530 spanned
4-wing subsets (Proposition 11).

Uses the exact integer matrix Z_{Delta,B} printed in Proposition 11, reads data_208a.json
from the same directory, and writes nonprod_data.bin alongside it.
Usage: python -I -B -X utf8 prep_nonprod_exact.py
"""
from __future__ import annotations

import json
import struct
from collections import defaultdict
from fractions import Fraction
from pathlib import Path

import numpy as np
import sympy as sp

HERE = Path(__file__).resolve().parent
DATA_PATH = HERE / "data_208a.json"
BIN_PATH = HERE / "nonprod_data.bin"

data = json.loads(DATA_PATH.read_text(encoding="utf-8"))
q_t = int(data["q_t"])
Z_t = [[int(x) for x in r] for r in data["Z_t"]]
D_t = np.array(Z_t, dtype=np.float64) / float(q_t)

pts = [((k >> 3) & 1, (k >> 2) & 1, (k >> 1) & 1, k & 1) for k in range(16)]
sgn = lambda b: 1 if b == 0 else -1
B = np.array(
    [
        [1, sgn(a0), sgn(a1), sgn(b0), sgn(b1),
         sgn(a0) * sgn(b0), sgn(a0) * sgn(b1), sgn(a1) * sgn(b0), sgn(a1) * sgn(b1)]
        for (a0, a1, b0, b1) in pts
    ],
    dtype=np.float64,
)
B_int = B.astype(np.int64)

DEN = 4096
Z_DQ_B = np.array(
    [
        [ -4, -35,  -8, -29, -15, -87, -24, -25,  27],
        [-35,  -4,  17, -87, -24, -29, -15,  12,  23],
        [ -8,  17,  -4, -25,  27,  12,  23, -29, -15],
        [-29, -87, -25,  -4,   1, -35,  -1,  -8,  31],
        [-15, -24,  27,   1,  -4,  -1, -35,  31,  -8],
        [-87, -29,  12, -35,  -1,  -4,   1,  17,  18],
        [-24, -15,  23,  -1, -35,   1,  -4,  35,  17],
        [-25,  12, -29,  -8,  31,  17,  35,  -4,   1],
        [ 27,  23, -15,  31,  -8,  18,  17,   1,  -4],
    ],
    dtype=np.int64,
)
H_null_B = np.zeros((9, 9), dtype=np.int64)
H_null_B[6, 7] = H_null_B[7, 6] = 1
H_null_B[5, 8] = H_null_B[8, 5] = -1

Z_Delta_pts = B_int @ Z_DQ_B @ B_int.T
Z_Hnull_pts = B_int @ H_null_B @ B_int.T
g_d = int(np.gcd.reduce(np.abs(Z_Delta_pts.ravel())))
g_h = int(np.gcd.reduce(np.abs(Z_Hnull_pts.ravel())))
Z_Delta_pts //= g_d
q_Delta = (256 * DEN) // g_d
Z_Hnull_pts //= g_h
q_Hnull = 256 // g_h
print(f"Reduced exact denominators: q_Delta={q_Delta}, q_Hnull={q_Hnull}")

D_Delta = Z_Delta_pts.astype(np.float64) / float(q_Delta)
D_Hnull = Z_Hnull_pts.astype(np.float64) / float(q_Hnull)

# Verify Sylvester's criterion on Q^(1) in orthonormal basis U_1 = B/4
c_pts = [int(round(B[k, 5] + B[k, 6] + B[k, 7] - B[k, 8])) // 2 for k in range(16)]
qcc_num = sum(c_pts[i] * Z_t[i][j] * c_pts[j] for i in range(16) for j in range(16))
Q_t_Z = [[sum(int(B_int[r, i]) * Z_t[r][c] * int(B_int[c, j]) for r in range(16) for c in range(16))
          for j in range(9)] for i in range(9)]
Q1_sp = sp.Matrix(
    9, 9,
    lambda i, j: sp.Rational(qcc_num * Q_t_Z[i][j], 16 * q_t * q_t) + sp.Rational(int(Z_DQ_B[i, j]), 16 * DEN),
)
assert all(Q1_sp[:k, :k].det() > 0 for k in range(1, 10))
ev_Q1 = np.linalg.eigvalsh(np.array(Q1_sp.tolist(), dtype=np.float64))
print(f"Q^(1) in U_1=B/4 strictly PD (Sylvester exact): lambda_min = {ev_Q1.min():.8e}")

B1 = np.array([[1, 1, 1], [1, 1, -1], [1, -1, 1], [1, -1, -1]], dtype=np.float64) / 2.0
B2 = np.kron(B1, B1)
B3 = np.kron(B2, B1)
u0 = np.array([1.0, 1.0, 1.0, 1.0])
up = np.array([1.0, 0.0, 0.0, -1.0])
um = np.array([0.0, 1.0, -1.0, 0.0])


def M_wing3(Dmat, x, y):
    M4 = np.einsum("ijkl,j,l->ik", Dmat.reshape(4, 4, 4, 4), x, y)
    return B1.T @ M4 @ B1


Q_orth = B2.T @ D_t @ B2
DQ_orth = B2.T @ D_Delta @ B2

H00 = 0.25 * (np.kron(Q_orth, M_wing3(D_t, u0, u0)) + np.kron(DQ_orth, M_wing3(D_Hnull, u0, u0)))
Hpp = 0.25 * (np.kron(Q_orth, M_wing3(D_t, up, up)) + np.kron(DQ_orth, M_wing3(D_Hnull, up, up)))
Hmm = 0.25 * (np.kron(Q_orth, M_wing3(D_t, um, um)) + np.kron(DQ_orth, M_wing3(D_Hnull, um, um)))
H0p = 0.50 * (np.kron(Q_orth, M_wing3(D_t, u0, up)) + np.kron(DQ_orth, M_wing3(D_Hnull, u0, up)))
H0m = 0.50 * (np.kron(Q_orth, M_wing3(D_t, u0, um)) + np.kron(DQ_orth, M_wing3(D_Hnull, u0, um)))

Mpm_Dt = M_wing3(D_t, up, um)
Mpm_Hn = M_wing3(D_Hnull, up, um)
sigma_Dt = Mpm_Dt[1, 2]
sigma_Hn = Mpm_Hn[1, 2]
p_vec = np.array([0.0, 1.0, 0.0])
q_vec = np.array([0.0, 0.0, 1.0])
assert np.max(np.abs(Mpm_Dt - sigma_Dt * (np.outer(p_vec, q_vec) - np.outer(q_vec, p_vec)))) < 1e-14
assert np.max(np.abs(Mpm_Hn - sigma_Hn * (np.outer(p_vec, q_vec) - np.outer(q_vec, p_vec)))) < 1e-14

Q_cross = 0.5 * (sigma_Dt * Q_orth + sigma_Hn * DQ_orth)
evC, vecC = np.linalg.eigh(Q_cross)
sqrt_abs = np.sqrt(np.abs(evC))
sgn_c = np.sign(evC)
A_bal = sqrt_abs[:, None] * vecC.T
B_bal = (sgn_c * sqrt_abs)[:, None] * vecC.T

Xmap = np.zeros((18, 27), dtype=np.float64)
Ymap = np.zeros((18, 27), dtype=np.float64)
for a in range(9):
    for b in range(3):
        col = 3 * a + b
        Xmap[:9, col] = A_bal[:, a] * p_vec[b]
        Xmap[9:, col] = A_bal[:, a] * q_vec[b]
        Ymap[:9, col] = B_bal[:, a] * q_vec[b]
        Ymap[9:, col] = -(B_bal[:, a] * p_vec[b])

A00 = H00 + Hpp + 0.5 * (H0p + H0p.T)
A33 = H00 + Hpp - 0.5 * (H0p + H0p.T)
A03 = 2.0 * (H00 - Hpp) + (H0p.T - H0p)
B11 = Hmm + 0.5 * (H0m + H0m.T)
B22 = Hmm - 0.5 * (H0m + H0m.T)
B12 = -2.0 * Hmm + (H0m.T - H0m)

S1 = [m for m in range(16) if ((m >> 0) & 1) + ((m >> 3) & 1) == ((m >> 1) & 1) + ((m >> 2) & 1)]
S2 = [
    m for m in range(65536)
    if all(((m >> (4 * wa)) & 15) in S1 for wa in range(4))
    and all(
        (((m >> wb) & 1) | (((m >> (4 + wb)) & 1) << 1) | (((m >> (8 + wb)) & 1) << 2) | (((m >> (12 + wb)) & 1) << 3)) in S1
        for wb in range(4)
    )
]
pairs2 = defaultdict(list)
for a in S2:
    for d in S2:
        pairs2[(a ^ d, a & d)].append((a, d))
S3 = [a | (b << 16) | (c << 32) | (d << 48) for lst in pairs2.values() for a, d in lst for b, c in lst]
S3_arr = np.array(S3, dtype=np.uint64)
BIT_SHIFTS = np.arange(64, dtype=np.uint64)
B3_bit = np.zeros((64, 27), dtype=np.float64)
for w3 in range(4):
    for c1 in range(16):
        B3_bit[16 * w3 + c1] = B3[4 * c1 + w3]

rng = np.random.default_rng(2026)
max_err = 0.0
for _ in range(500):
    g0 = S3_arr[rng.integers(len(S3_arr))]
    g1 = S3_arr[rng.integers(len(S3_arr))]
    cands = S3_arr[((S3_arr ^ g0) & (g1 ^ g0)) == np.uint64(0)]
    g2 = cands[rng.integers(len(cands))]
    g3 = g1 ^ g2 ^ g0
    G_arr = np.stack([((np.uint64(x) >> BIT_SHIFTS) & np.uint64(1)).astype(float).reshape(4, 16) for x in (g0, g1, g2, g3)], axis=0)
    F = G_arr.transpose(2, 1, 0).reshape(16, 16)
    v_true = float(np.sum((D_t @ F) * (F @ D_t)) + np.sum((D_Delta @ F) * (F @ D_Hnull)))

    z = [B3_bit.T @ ((np.uint64(g) >> BIT_SHIFTS) & np.uint64(1)).astype(float) for g in (g0, g1, g2, g3)]
    alpha = float(z[0] @ A00 @ z[0] + z[3] @ A33 @ z[3] + z[0] @ A03 @ z[3])
    beta = float(z[1] @ B11 @ z[1] + z[2] @ B22 @ z[2] + z[1] @ B12 @ z[2])
    cross = float((Xmap @ (z[0] - z[3])) @ (Ymap @ (z[1] - z[2])))
    max_err = max(max_err, abs(v_true - (alpha + beta + cross)))
print(f"500-sample formula check max_err = {max_err:.3e}")
assert max_err < 1e-13

with BIN_PATH.open("wb") as f:
    for r in range(16):
        for c in range(16):
            f.write(struct.pack("<q", Z_t[r][c]))
    for r in range(16):
        for c in range(16):
            f.write(struct.pack("<q", int(Z_Delta_pts[r, c])))
    for r in range(16):
        for c in range(16):
            f.write(struct.pack("<q", int(Z_Hnull_pts[r, c])))
    f.write(struct.pack("<qqq", q_t, q_Delta, q_Hnull))
    for arr in (B3_bit, A00, A33, A03, B11, B22, B12, Xmap, Ymap):
        f.write(np.ascontiguousarray(arr, dtype="<f8").tobytes())
print(f"Wrote {BIN_PATH.name} cleanly.")
