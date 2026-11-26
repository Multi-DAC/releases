"""Prop 4 at a5 in (1, a5max(m)] (beyond the sealed plan, which capped a5 at 1; written after clock reading 14:54).

At m = 0.3, a5max = 2.745, so the grid left a5 in (1, 2.745] untested. Same quantity and t grid as prop4_large_t.py.
"""
import numpy as np

from l4_assembled import a5max
from prop4_large_t import TS, ratios

print("t grid:", TS)
for m in (0.3, 0.5):
    for a5 in (1.5, 2.0, 0.999 * a5max(m)):
        if a5 > a5max(m):
            continue
        for k1 in (0.0, np.pi / 16, np.pi / 2, np.pi):
            for L5 in (16, 32, 64):
                try:
                    r, tdw, tpv, ns = ratios(L5, k1, m, a5)
                except AssertionError as e:
                    print("m=%.1f a5=%.4f k1=%.3f L5=%3d: FAILED %s" % (m, a5, k1, L5, e), flush=True)
                    continue
                i = int(np.argmax(r))
                print("m=%.1f a5=%.4f k1=%.3f L5=%3d S=%6.1f: max %.4f at t=%5d ; top DW/r* %.6f PV/r* %.6f nsurf %d ;"
                      " r(t) = %s" % (m, a5, k1, L5, L5 * a5, r[i], TS[i], tdw, tpv, ns,
                                      " ".join("%.3g" % x for x in r)), flush=True)
print("DONE")
