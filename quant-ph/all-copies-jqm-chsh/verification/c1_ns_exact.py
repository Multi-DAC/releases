"""Independent check (Clawd) of Proposition 4: an EXACT form for the PR box at one copy with branching (C_1 = NS needs PR in C_1).
Takes the float LP vertex of c1_ns.py, snaps each coordinate to a rational, then verifies in fractions:
(i) on 16 events, (ii) on 56 pairs, K(1_S, 1_S) >= 0 on the 89 non-empty spanned subsets."""
import json
import os
from fractions import Fraction as Fr

import numpy as np
from scipy.optimize import linprog

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import c1_ns as m  # noqa: E402

S = m.spanned()
pairs = [(e, f) for k, e in enumerate(m.EV) for f in m.EV[k + 1:] if m.excl(e, f)]
Aeq = [m.row(m.IND[e], m.IND[e]) for e in m.EV] + [m.row(m.IND[e], m.IND[f]) for (e, f) in pairs]
beq = [m.pr(e) for e in m.EV] + [0.0] * len(pairs)
cuts = np.array([m.row(v, v) for v in S])
n = len(m.IDX)
# push into the interior where possible: maximise the SUM of K(1_S,1_S) over spanned S subject to all >= 0
res = linprog(-cuts.sum(axis=0), A_ub=-cuts, b_ub=np.zeros(len(S)), A_eq=np.array(Aeq), b_eq=beq,
              bounds=[(-50, 50)] * n, method="highs")
assert res.status == 0, res.message
k = [Fr(float(v)).limit_denominator(4096) for v in res.x]
print("largest denominator after snapping:", max(v.denominator for v in k))


from c1_frow import frow as frow_exact


def frow_old(u, v):
    return [Fr(float(x)).limit_denominator(64) for x in m.row(u, v)]


bad_i = sum(1 for e in m.EV if sum(a * b for a, b in zip(frow_exact(m.IND[e], m.IND[e]), k)) != Fr(m.pr(e)))
bad_ii = sum(1 for (e, f) in pairs if sum(a * b for a, b in zip(frow_exact(m.IND[e], m.IND[f]), k)) != 0)
vals = [sum(a * b for a, b in zip(frow_exact(v, v), k)) for v in S]
print("EXACT: (i) fails on %d of 16 ; (ii) fails on %d of 56 ; spanned sets with K < 0: %d of %d ; smallest %s ; zeros %d"
      % (bad_i, bad_ii, sum(1 for v in vals if v < 0), len(vals), min(vals), sum(1 for v in vals if v == 0)))
K = [[None] * 9 for _ in range(9)]
for v, (i, j) in zip(k, m.IDX):
    K[i][j] = K[j][i] = str(v)
json.dump({"about": "PR box, one copy, symmetric form in Walsh coordinates (order 1,b0,b1,a0,a0b0,a0b1,a1,a1b0,a1b1 as in c1_ns.walsh); K(u,v) = cu^T K cv",
           "K": K}, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "c1_ns_exact.json"), "w"), indent=1)
