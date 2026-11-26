"""Prop 4 beyond t = 10^4 for the near-tie cases still rising at 10^4 (written between clock readings 14:47 and 14:52).

Same quantity as prop4_large_t.py: ||a5 Delta(t+1, k1)||_1 / r_*(k1)^t, exact pole sums, poles normalised by r_*.
Also the exact tie a5 = a5max at L5 = 16 (the earlier small run: 0.359 flat), carried to t = 10^7.
"""
import numpy as np

from l4_assembled import Prop, a5max, rstar

TT = [0, 1000, 5000, 10000, 20000, 50000, 100000, 200000, 500000, 1000000, 3000000, 10000000]


def run(L5, k1, m, a5):
    rs = rstar(k1, m)
    DW, PV = Prop(L5, k1, m, a5, 0.0), Prop(L5, k1, m, a5, 1.0)
    bulk = [(l, R) for l, R, _, _ in DW.pairs if l <= rs * (1 + 1e-12)]
    pv = [(l, R) for l, R, _, _ in PV.pairs]
    out = []
    for t in TT:
        A = sum((l / rs) ** t * R for l, R in bulk) - sum((l / rs) ** t * R for l, R in pv)
        out.append(a5 * np.linalg.svd(A, compute_uv=False).sum())
    top = max(l for l, _ in bulk) / rs
    return out, top, len(DW.unit), max(DW.spread, PV.spread)


cases = [(0.7, 0.999 * a5max(0.7), 0.0, 32), (0.7, 0.999 * a5max(0.7), 0.0, 64), (0.7, 0.999 * a5max(0.7), 0.0, 128),
         (0.9, 0.999 * a5max(0.9), 0.0, 64), (0.9, 0.999 * a5max(0.9), 0.0, 128), (0.9, 0.2, 0.0, 128),
         (0.7, 0.999 * a5max(0.7), np.pi / 16, 128), (0.3, 1.0, 0.0, 128), (0.3, 1.0, np.pi, 128)]
print("t:", TT)
for m, a5, k1, L5 in cases:
    try:
        r, top, nu, sp = run(L5, k1, m, a5)
    except AssertionError as e:
        print("m=%.1f a5=%.4f k1=%.3f L5=%3d: FAILED %s" % (m, a5, k1, L5, e), flush=True)
        continue
    i = int(np.argmax(r))
    print("m=%.1f a5=%.4f k1=%.3f L5=%3d: max %.4f at t=%d ; 1 - top/r* %.2g ; unit groups %d ; spread %.2g ;"
          " r = %s" % (m, a5, k1, L5, r[i], TT[i], 1 - top, nu, sp, " ".join("%.3g" % x for x in r)), flush=True)
for L5 in (16, 32):
    try:
        r, top, nu, sp = run(L5, 0.0, 0.7, a5max(0.7))
        print("EXACT TIE m=0.7 a5=a5max L5=%d: 1 - top/r* %.2g ; r = %s" % (L5, 1 - top, " ".join("%.3g" % x for x in r)),
              flush=True)
    except AssertionError as e:
        print("EXACT TIE m=0.7 L5=%d: FAILED %s" % (L5, e), flush=True)
print("DONE")
