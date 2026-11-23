from fractions import Fraction as Fr

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import c1_ns as m  # noqa: E402


def coords(u):
    return [Fr(int(round(float(sum(u[p] * m.B[p][a] for p in range(16))))), 16) for a in range(9)]


def frow(u, v):
    cu, cv = coords(u), coords(v)
    return [cu[i] * cv[j] + (cu[j] * cv[i] if i != j else 0) for (i, j) in m.IDX]
