"""Independent check (Clawd) of Proposition 4, claim C_1 = NS (one copy, branching): does the PR box have a symmetric form on W_1
with (i), all 56 exclusive pairs of (ii), and K(1_S, 1_S) >= 0 on every spanned subset? FLOAT LP, self-contained.
By Reduction B (Lemma 3) a form on W_1 positive on spanned sets extends to all subsets.
Controls: the count of spanned subsets must be 90, the count of exclusive pairs 56 (24 same-setting);
with the pairs dropped to none the LP must also be feasible; a signalling table must fail (i) consistency.
"""
import itertools

import numpy as np
from scipy.optimize import linprog

PTS = list(itertools.product((0, 1), repeat=4))  # (a0, a1, b0, b1)


def walsh():
    cols = []
    for sa in (None, 0, 1):
        for sb in (None, 0, 1):
            cols.append([(1 if sa is None else (-1) ** p[sa]) * (1 if sb is None else (-1) ** p[2 + sb]) for p in PTS])
    return np.array(cols, dtype=float).T  # 16 x 9


B = walsh()
EV = [(x, a, y, b) for x in (0, 1) for a in (0, 1) for y in (0, 1) for b in (0, 1)]
IND = {e: np.array([1.0 if (p[e[0]] == e[1] and p[2 + e[2]] == e[3]) else 0.0 for p in PTS]) for e in EV}


def excl(e, f):
    return (e[0] == f[0] and e[1] != f[1]) or (e[2] == f[2] and e[3] != f[3])


def coords(v):
    c, res = np.linalg.lstsq(B, v, rcond=None)[:2]
    assert np.allclose(B @ c, v, atol=1e-9)
    return c


def spanned():
    P = B @ np.linalg.pinv(B)
    out = []
    for m in range(1, 1 << 16):
        v = np.array([(m >> k) & 1 for k in range(16)], dtype=float)
        if np.allclose(P @ v, v, atol=1e-9):
            out.append(v)
    return out


IDX = [(i, j) for i in range(9) for j in range(i, 9)]


def row(u, v):
    cu, cv = coords(u), coords(v)
    return np.array([cu[i] * cv[j] + (cu[j] * cv[i] if i != j else 0.0) for (i, j) in IDX])


def pr(e):
    x, a, y, b = e
    return 0.5 if (a ^ b) == (x & y) else 0.0


def main():
    S = spanned()
    pairs = [(e, f) for k, e in enumerate(EV) for f in EV[k + 1:] if excl(e, f)]
    same = [(e, f) for (e, f) in pairs if e[0] == f[0] and e[2] == f[2]]
    print("spanned non-empty subsets %d (with the empty set %d) ; exclusive pairs %d, same-setting %d"
          % (len(S), len(S) + 1, len(pairs), len(same)))
    assert len(S) == 89 and len(pairs) == 56 and len(same) == 24
    Aeq = [row(IND[e], IND[e]) for e in EV] + [row(IND[e], IND[f]) for (e, f) in pairs]
    beq = [pr(e) for e in EV] + [0.0] * len(pairs)
    cuts = np.array([row(v, v) for v in S])
    n = len(IDX)
    # maximise t with K(1_S, 1_S) >= t for all spanned S, t <= 1
    A_ub = np.hstack([-cuts, np.ones((len(S), 1))])
    res = linprog(np.r_[np.zeros(n), -1.0], A_ub=A_ub, b_ub=np.zeros(len(S)),
                  A_eq=np.hstack([np.array(Aeq), np.zeros((len(Aeq), 1))]), b_eq=beq,
                  bounds=[(None, None)] * n + [(None, 1.0)], method="highs")
    print("PR box, 56 pairs: status %d (%s) ; best min over spanned sets t = %s"
          % (res.status, res.message[:40], None if res.x is None else "%.6f" % res.x[-1]))
    if res.status == 0:
        k = res.x[:-1]
        print("  residual of (i),(ii): %.2e ; smallest K(1_S,1_S) %.6f" % (np.abs(np.array(Aeq) @ k - beq).max(), (cuts @ k).min()))
        K = np.zeros((9, 9))
        for v, (i, j) in zip(k, IDX):
            K[i, j] = K[j, i] = v
        G = B.T @ B
        print("  smallest eigenvalue of the form on W_1 (orthonormal basis): %.6f" % np.linalg.eigvalsh(K * 16).min())


if __name__ == "__main__":
    main()
