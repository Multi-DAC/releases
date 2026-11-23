"""Independent check (Clawd) of the claim in Section 3 and Remark 15: p_t lies in the interior of Q^1.

Exact (fractions). Reads p_t from data_208a.json (event keys 'x a y b'). Builds the 5 x 5 level-1 NPA moment matrix
in the +-1 observables (1, A_0, A_1, B_0, B_1) with the two free entries <A_0 A_1> = <B_0 B_1> = 0 and checks that
all five elimination pivots are strictly positive. A strictly positive definite completion stays positive definite
under small changes of the eight correlators, so p_t is interior to Q^1 inside the no-signalling polytope.
Controls (same test, same free entries): the noisy PR box q PR + (1 - q) I at q = 7/10 (CHSH 2.8) must pass and at
q = 3/4 (CHSH 3) must fail.
"""
import json
import os
import sys
from fractions import Fraction as Fr

HERE = os.path.dirname(os.path.abspath(__file__))


def sgn(a):
    return 1 - 2 * a


def pivots(M):
    n = len(M)
    M = [[Fr(v) for v in r] for r in M]
    out = []
    for k in range(n):
        out.append(M[k][k])
        if M[k][k] <= 0:
            return out, False
        for i in range(k + 1, n):
            f = M[i][k] / M[k][k]
            for j in range(k, n):
                M[i][j] -= f * M[k][j]
    return out, True


def gamma(p):
    """p: dict 'xayb' -> probability. Returns the 5 x 5 matrix and the CHSH value."""
    for x in (0, 1):                      # no-signalling, exactly
        for a in (0, 1):
            assert sum(p["%d%d0%d" % (x, a, b)] for b in (0, 1)) == sum(p["%d%d1%d" % (x, a, b)] for b in (0, 1))
    for y in (0, 1):
        for b in (0, 1):
            assert sum(p["0%d%d%d" % (a, y, b)] for a in (0, 1)) == sum(p["1%d%d%d" % (a, y, b)] for a in (0, 1))
    mA = [sum(sgn(a) * p["%d%d0%d" % (x, a, b)] for a in (0, 1) for b in (0, 1)) for x in (0, 1)]
    mB = [sum(sgn(b) * p["0%d%d%d" % (a, y, b)] for a in (0, 1) for b in (0, 1)) for y in (0, 1)]
    E = {(x, y): sum(sgn(a) * sgn(b) * p["%d%d%d%d" % (x, a, y, b)] for a in (0, 1) for b in (0, 1))
         for x in (0, 1) for y in (0, 1)}
    G = [[1, mA[0], mA[1], mB[0], mB[1]],
         [mA[0], 1, 0, E[(0, 0)], E[(0, 1)]],
         [mA[1], 0, 1, E[(1, 0)], E[(1, 1)]],
         [mB[0], E[(0, 0)], E[(1, 0)], 1, 0],
         [mB[1], E[(0, 1)], E[(1, 1)], 0, 1]]
    return G, E[(0, 0)] + E[(0, 1)] + E[(1, 0)] - E[(1, 1)]


def prq(q):
    return {"%d%d%d%d" % (x, a, y, b): q * (Fr(1, 2) if (a ^ b) == (x & y) else 0) + (1 - q) * Fr(1, 4)
            for x in (0, 1) for a in (0, 1) for y in (0, 1) for b in (0, 1)}


def main():
    d = json.load(open(os.path.join(HERE, "data_208a.json")))
    pt = {k: Fr(v) for k, v in d["p_t"].items()}
    assert len(pt) == 16
    G, chsh = gamma(pt)
    piv, ok = pivots(G)
    print("p_t: CHSH = %s ; pivots %s ; all positive: %s" % (chsh, [str(v) for v in piv], ok))
    assert ok and chsh == Fr(6141, 2560)
    for q, want in ((Fr(7, 10), True), (Fr(3, 4), False)):
        G, chsh = gamma(prq(q))
        piv, ok = pivots(G)
        print("CONTROL noisy PR box q = %s: CHSH = %s ; pivots %s ; all positive: %s"
              % (q, chsh, [str(v) for v in piv], ok))
        assert ok == want
    print("DONE")


if __name__ == "__main__":
    sys.exit(main())
