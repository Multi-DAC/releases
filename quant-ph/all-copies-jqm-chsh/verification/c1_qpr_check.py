"""Independent check (Clawd) of Proposition 4: the printed 9 x 9 matrix Q_PR.

My own code, exact (integers / fractions). Reads Q_PR from the paper TeX. D_PR = B Q_PR B^T / 256.
Checks: (i) on 16 events for the PR box, (ii) on 56 exclusive pairs, D_PR(1_S, 1_S) on the 89 non-empty spanned sets,
Q_PR + I positive semidefinite (exact, by Gaussian elimination with zero-pivot handling), and
D_PR + (1/2) Pi_{W_1^perp} on all 65,535 non-empty subsets.
Controls: one free entry of Q_PR changed (2 -> 3) must break something; one fixed entry changed must break (i) or (ii).
"""
import os
import re
import sys
from fractions import Fraction as Fr

HERE = os.path.dirname(os.path.abspath(__file__))
PTS = [((k >> 3) & 1, (k >> 2) & 1, (k >> 1) & 1, k & 1) for k in range(16)]
MON = [(None, None), (0, None), (1, None), (None, 0), (None, 1), (0, 0), (0, 1), (1, 0), (1, 1)]
B = [[(1 if x is None else (-1) ** p[x]) * (1 if y is None else (-1) ** p[2 + y]) for (x, y) in MON] for p in PTS]
EV = [(x, a, y, b) for x in (0, 1) for a in (0, 1) for y in (0, 1) for b in (0, 1)]


def ind(e):
    return [1 if (p[e[0]] == e[1] and p[2 + e[2]] == e[3]) else 0 for p in PTS]


def excl(e, f):
    return (e[0] == f[0] and e[1] != f[1]) or (e[2] == f[2] and e[3] != f[3])


def coords(v):
    """B^T v, integers"""
    return [sum(B[k][a] * v[k] for k in range(16)) for a in range(9)]


def qf(Q, cu, cv):
    return sum(cu[a] * Q[a][b] * cv[b] for a in range(9) for b in range(9))


def read_Q(tex):
    lines = open(tex, encoding="utf-8").read().split("\n")
    start = next(i for i, s in enumerate(lines) if "Q_{\\PR} \\;:=\\;" in s)
    rows = []
    for s in lines[start:start + 14]:
        nums = re.findall(r"-?\d+", s)
        if "&" in s and len(nums) == 9:
            rows.append([int(v) for v in nums])
    assert len(rows) == 9, len(rows)
    assert all(rows[i][j] == rows[j][i] for i in range(9) for j in range(9)), "Q_PR not symmetric"
    return rows


def psd_exact(M):
    """True iff the symmetric rational matrix M is positive semidefinite (elimination; a zero pivot needs a zero row)"""
    n = len(M)
    M = [[Fr(v) for v in r] for r in M]
    for k in range(n):
        if M[k][k] < 0:
            return False
        if M[k][k] == 0:
            if any(M[k][j] != 0 for j in range(k, n)):
                return False
            continue
        for i in range(k + 1, n):
            f = M[i][k] / M[k][k]
            for j in range(k, n):
                M[i][j] -= f * M[k][j]
    return True


def pr(e):
    x, a, y, b = e
    return Fr(1, 2) if (a ^ b) == (x & y) else Fr(0)


def run(Q, tag):
    cE = {e: coords(ind(e)) for e in EV}
    bad_i = sum(1 for e in EV if Fr(qf(Q, cE[e], cE[e]), 256) != pr(e))
    pairs = [(e, f) for i, e in enumerate(EV) for f in EV[i + 1:] if excl(e, f)]
    assert len(pairs) == 56
    bad_ii = sum(1 for (e, f) in pairs if qf(Q, cE[e], cE[f]) != 0)
    # all subsets: spanned iff B B^T v = 16 v
    neg_sp = zeros = nsp = neg_all = 0
    min_sp = None
    min_all = None
    for m in range(1, 1 << 16):
        v = [(m >> k) & 1 for k in range(16)]
        c = coords(v)
        val256 = qf(Q, c, c)                       # 256 * D_PR(1_S, 1_S)
        cc = sum(x * x for x in c)                  # 16 * |Pi_W 1_S|^2
        size = sum(v)
        spanned = (cc == 16 * size)
        if spanned:
            nsp += 1
            neg_sp += val256 < 0
            zeros += val256 == 0
            min_sp = val256 if min_sp is None else min(min_sp, val256)
        full512 = 2 * val256 + 256 * size - 16 * cc  # 512 * (D_PR + (1/2) Pi_perp)(1_S, 1_S)
        neg_all += full512 < 0
        min_all = full512 if min_all is None else min(min_all, full512)
    QI = [[Q[i][j] + (1 if i == j else 0) for j in range(9)] for i in range(9)]
    psd = psd_exact(QI)
    psd_Q = psd_exact(Q)
    print("%s: (i) fails on %d of 16 ; (ii) fails on %d of 56 ; spanned non-empty sets %d, negative %d, zero %d, "
          "smallest %s ; Q_PR + I positive semidefinite: %s (Q_PR itself: %s) ; D_PR + (1/2) Pi_perp negative on %d of "
          "65535 non-empty subsets, smallest %s"
          % (tag, bad_i, bad_ii, nsp, neg_sp, zeros, Fr(min_sp, 256), psd, psd_Q, neg_all, Fr(min_all, 512)))
    return bad_i == 0 and bad_ii == 0 and nsp == 89 and neg_sp == 0 and psd and neg_all == 0


def main():
    Q = read_Q(os.path.join(HERE, "..", "documentation", "paper.tex"))
    for r in Q:
        print("  ", r)
    assert run(Q, "as printed")
    Q2 = [r[:] for r in Q]
    Q2[6][7] = Q2[7][6] = 3
    assert not run(Q2, "CONTROL free entry (alpha0 beta1, alpha1 beta0) 2 -> 3")
    Q3 = [r[:] for r in Q]
    Q3[1][3] = Q3[3][1] = 0
    assert not run(Q3, "CONTROL fixed entry (alpha0, beta0) 1 -> 0")
    print("DONE")


if __name__ == "__main__":
    sys.exit(main())
