"""Prop 4 at large t: ||a5 Delta(t+1, k1)||_1 / r_*(k1)^t with exact pole sums, t up to 10^4 (clock reading 14:38).

Delta = G_DW^bulk - G_PV. My operator (l4_assembled.Prop, validated against position space at 13:44).
Plan sealed first: prop4-large-t-plan.md.
"""
import sys

import numpy as np

from l4_assembled import Prop, a5max, rstar

TS = [0, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000]


def ratios(L5, k1, m, a5):
    rs = rstar(k1, m)
    DW, PV = Prop(L5, k1, m, a5, 0.0), Prop(L5, k1, m, a5, 1.0)
    bulk = [(l, R) for l, R, _, _ in DW.pairs if l <= rs * (1 + 1e-12)]
    nsurf = len(DW.pairs) - len(bulk)
    pv = [(l, R) for l, R, _, _ in PV.pairs]
    top_dw = max(l for l, _ in bulk) / rs
    top_pv = max(l for l, _ in pv) / rs
    out = []
    for t in TS:
        A = np.zeros_like(bulk[0][1])
        for l, R in bulk:
            A = A + (l / rs) ** t * R
        for l, R in pv:
            A = A - (l / rs) ** t * R
        out.append(a5 * np.linalg.svd(A, compute_uv=False).sum())
    return out, top_dw, top_pv, nsurf


if __name__ == "__main__":
    ms = [0.3, 0.7, 0.9]
    k1s = [0.0, np.pi / 16, np.pi / 8, np.pi / 2, np.pi]
    L5s = [int(x) for x in sys.argv[1].split(",")] if len(sys.argv) > 1 else [16, 32, 64, 128]
    print("t grid:", TS)
    worst = (0.0, None)
    for m in ms:
        # exact ties (a5 = a5max, k1 = 0) put every bulk pole on r_* and break the pairing; 0.999 a5max sits next to them
        for a5 in sorted({0.05, 0.1, 0.2, min(0.999 * a5max(m), 1.0)}):
            if a5 > a5max(m) + 1e-12:
                continue
            for k1 in k1s:
                for L5 in L5s:
                    try:
                        r, tdw, tpv, ns = ratios(L5, k1, m, a5)
                    except AssertionError as e:
                        print("m=%.1f a5=%.4f k1=%.3f L5=%3d: FAILED pole pairing %s" % (m, a5, k1, L5, e))
                        continue
                    i = int(np.argmax(r))
                    if r[i] > worst[0]:
                        worst = (r[i], (m, a5, k1, L5, TS[i]))
                    print("m=%.1f a5=%.4f k1=%.3f L5=%3d S=%6.1f: max %.4f at t=%5d ; top DW/r* %.6f PV/r* %.6f nsurf %d ;"
                          " r(t) = %s" % (m, a5, k1, L5, L5 * a5, r[i], TS[i], tdw, tpv, ns,
                                          " ".join("%.3g" % x for x in r)))
                    sys.stdout.flush()
    print("WORST", worst)
    print("DONE")
