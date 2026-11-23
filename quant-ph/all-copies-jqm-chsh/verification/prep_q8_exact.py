"""Generate binary matrix data file q8_data.bin for the exact global C# verifier
VerifyQ8.cs over all 55,720,396,530 spanned 4-wing subsets (Theorem 5(b)).

Reads data_208a.json from the same directory and writes q8_data.bin alongside it.
Usage: python -I -B -X utf8 prep_q8_exact.py
"""
from __future__ import annotations

import json
from collections import defaultdict
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
DATA_PATH = HERE / "data_208a.json"
OUT_BIN = HERE / "q8_data.bin"

data = json.loads(DATA_PATH.read_text(encoding="utf-8"))
q_t = int(data["q_t"])
Z_t = [[int(x) for x in r] for r in data["Z_t"]]
D = np.array(Z_t, dtype=np.float64) / float(q_t)

# Build S1, S2, S3
S1 = [m for m in range(16) if ((m >> 0) & 1) + ((m >> 3) & 1) == ((m >> 1) & 1) + ((m >> 2) & 1)]
S2 = [
    m for m in range(65536)
    if all(((m >> (4 * wa)) & 15) in S1 for wa in range(4))
    and all(
        (((m >> wb) & 1) | (((m >> (4 + wb)) & 1) << 1) | (((m >> (8 + wb)) & 1) << 2) | (((m >> (12 + wb)) & 1) << 3)) in S1
        for wb in range(4)
    )
]

# Wing-4 decomposition: f = 0.5 * (s (x) u0 + d1 (x) up + d2 (x) um)
B1 = np.array([[1, 1, 1], [1, 1, -1], [1, -1, 1], [1, -1, -1]], dtype=np.float64) / 2.0
B2 = np.kron(B1, B1)  # 16 x 9 orthonormal on copy 1
B3 = np.kron(B2, B1)  # 64 x 27 orthonormal on (A1, B1, A2)

D4 = D.reshape(4, 4, 4, 4)
u0 = np.array([1.0, 1.0, 1.0, 1.0])
up = np.array([1.0, 0.0, 0.0, -1.0])
um = np.array([0.0, 1.0, -1.0, 0.0])


def M_wing3(x, y):
    M4 = np.einsum("ijkl,j,l->ik", D4, x, y)
    return B1.T @ M4 @ B1


Q = B2.T @ D @ B2
M00 = M_wing3(u0, u0)
Mpp = M_wing3(up, up)
Mmm = M_wing3(um, um)
M0p = M_wing3(u0, up)
M0m = M_wing3(u0, um)
Mpm = M_wing3(up, um)  # 3x3 skew-symmetric

assert np.max(np.abs(Mpm + Mpm.T)) < 1e-14
U_svd, S_svd, _ = np.linalg.svd(Mpm)
sigma = S_svd[0]
assert abs(S_svd[0] - S_svd[1]) < 1e-13 and S_svd[2] < 1e-13
p_vec = U_svd[:, 0]
q_vec = -(Mpm @ p_vec) / sigma
assert abs(np.linalg.norm(q_vec) - 1.0) < 1e-13
assert abs(p_vec @ q_vec) < 1e-13
assert np.max(np.abs(Mpm - sigma * (np.outer(p_vec, q_vec) - np.outer(q_vec, p_vec)))) < 1e-14

pairs2 = defaultdict(list)
for a in S2:
    for d in S2:
        pairs2[(a ^ d, a & d)].append((a, d))
S3 = [a | (b << 16) | (c << 32) | (d << 48) for lst in pairs2.values() for a, d in lst for b, c in lst]
S3_arr = np.array(S3, dtype=np.uint64)
BIT_SHIFTS = np.arange(64, dtype=np.uint64)

H00 = 0.25 * np.kron(Q, M00)
Hpp = 0.25 * np.kron(Q, Mpp)
Hmm = 0.25 * np.kron(Q, Mmm)
H0p = 0.50 * np.kron(Q, M0p)
H0m = 0.50 * np.kron(Q, M0m)

evQ, vecQ = np.linalg.eigh(Q)
sqrt_abs_lam = np.sqrt(np.abs(evQ))
sgn_lam = np.sign(evQ)
A_bal = sqrt_abs_lam[:, None] * vecQ.T
B_bal = (sgn_lam * sqrt_abs_lam)[:, None] * vecQ.T
assert np.max(np.abs(A_bal.T @ B_bal - Q)) < 1e-14
sc = np.sqrt(0.5 * sigma)


def XY_of(z):
    Z93 = z.reshape(9, 3)
    up_v = Z93 @ p_vec
    uq_v = Z93 @ q_vec
    X = sc * np.concatenate([A_bal @ up_v, A_bal @ uq_v])
    Y = sc * np.concatenate([B_bal @ uq_v, -(B_bal @ up_v)])
    return X, Y


rng = np.random.default_rng(999)
max_err = 0.0
for _ in range(500):
    g0 = S3_arr[rng.integers(len(S3_arr))]
    g1 = S3_arr[rng.integers(len(S3_arr))]
    cands = S3_arr[((S3_arr ^ g0) & (g1 ^ g0)) == np.uint64(0)]
    g2 = cands[rng.integers(len(cands))]
    g3 = g1 ^ g2 ^ g0
    G_arr = np.stack([((np.uint64(x) >> BIT_SHIFTS) & np.uint64(1)).astype(float).reshape(4, 16) for x in (g0, g1, g2, g3)], axis=0)
    F = G_arr.transpose(2, 1, 0).reshape(16, 16)
    v_true = float(np.sum((D @ F) * (F @ D)))

    z0 = B3.T @ (((np.uint64(g0) >> BIT_SHIFTS) & np.uint64(1)).astype(float).reshape(4, 16).T.ravel())
    z1 = B3.T @ (((np.uint64(g1) >> BIT_SHIFTS) & np.uint64(1)).astype(float).reshape(4, 16).T.ravel())
    z2 = B3.T @ (((np.uint64(g2) >> BIT_SHIFTS) & np.uint64(1)).astype(float).reshape(4, 16).T.ravel())
    z3 = B3.T @ (((np.uint64(g3) >> BIT_SHIFTS) & np.uint64(1)).astype(float).reshape(4, 16).T.ravel())
    s = z0 + z3
    d1 = z0 - z3
    d2 = z1 - z2
    c0 = float(s @ H00 @ s)
    a_val = float(s @ H0p @ d1 + d1 @ Hpp @ d1)
    b_val = float(s @ H0m @ d2 + d2 @ Hmm @ d2)
    X1, _ = XY_of(d1)
    _, Y2 = XY_of(d2)
    v_form = c0 + a_val + b_val + float(X1 @ Y2)
    max_err = max(max_err, abs(v_true - v_form))

print(f"Verified V-basis formula on 500 random spanned 4-wing subsets: max_err = {max_err:.3e}")
assert max_err < 1e-13

A00 = 0.5 * H00 + 0.5 * (H0p + H0p.T) + Hpp
A33 = 0.5 * H00 - 0.5 * (H0p + H0p.T) + Hpp
A03 = H00 - H0p + H0p.T - 2.0 * Hpp
B11 = 0.5 * H00 + 0.5 * (H0m + H0m.T) + Hmm
B22 = 0.5 * H00 - 0.5 * (H0m + H0m.T) + Hmm
B12 = H00 - H0m + H0m.T - 2.0 * Hmm

X_map = np.zeros((18, 27))
Y_map = np.zeros((18, 27))
for col in range(27):
    e = np.zeros(27)
    e[col] = 1.0
    x, y = XY_of(e)
    X_map[:, col] = x
    Y_map[:, col] = y

B3_bit = np.zeros((64, 27), dtype=np.float64)
for w3 in range(4):
    for c1 in range(16):
        B3_bit[16 * w3 + c1] = B3[4 * c1 + w3]

with OUT_BIN.open("wb") as f:
    np.array(Z_t, dtype="<i8").tofile(f)
    B3_bit.astype("<f8").tofile(f)
    for M in (A00, A33, A03, B11, B22, B12):
        M.astype("<f8").tofile(f)
    X_map.astype("<f8").tofile(f)
    Y_map.astype("<f8").tofile(f)
print(f"Wrote {OUT_BIN.name} cleanly.")
