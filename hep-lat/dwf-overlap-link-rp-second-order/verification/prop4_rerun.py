"""Rerun the 8 FAILED cases of prop4_large_t.py with the joint clustering (written between clock readings 14:47 and 14:52).

Controls first: four cases that passed under the old clustering must give the same numbers under the new one.
Same line format as prop4_large_t.py, so prop4_parse.py reads it; each line adds the widest merged group in
log|lam| (DW and PV), and 10^4 times it, the relative error that merging costs at the largest t.
"""
import numpy as np

from l4_assembled import Prop, a5max, rstar
from prop4_large_t import TS


def ratios(L5, k1, m, a5):
    rs = rstar(k1, m)
    DW, PV = Prop(L5, k1, m, a5, 0.0), Prop(L5, k1, m, a5, 1.0)
    bulk = [(l, R) for l, R, _, _ in DW.pairs if l <= rs * (1 + 1e-12)]
    assert not PV.unit, PV.unit
    nsurf = len(DW.pairs) - len(bulk) + len(DW.unit)   # (14:47-14:52) a unit-circle group is the surface pole, below precision
    pv = [(l, R) for l, R, _, _ in PV.pairs]
    out = []
    for t in TS:
        A = np.zeros_like(bulk[0][1])
        for l, R in bulk:
            A = A + (l / rs) ** t * R
        for l, R in pv:
            A = A - (l / rs) ** t * R
        out.append(a5 * np.linalg.svd(A, compute_uv=False).sum())
    return (out, max(l for l, _ in bulk) / rs, max(l for l, _ in pv) / rs, nsurf, max(DW.spread, PV.spread))


def tie(m):
    return min(0.999 * a5max(m), 1.0)


controls = [(0.3, 1.0, 0.0, 32), (0.7, tie(0.7), 0.0, 16), (0.9, tie(0.9), 0.0, 64), (0.9, 0.2, 0.0, 64)]
failed = [(0.3, 1.0, 0.0, 64), (0.3, 1.0, 0.0, 128), (0.7, 0.2, 0.0, 128), (0.7, tie(0.7), 0.0, 32),
          (0.7, tie(0.7), 0.0, 64), (0.7, tie(0.7), 0.0, 128), (0.9, 0.2, 0.0, 128), (0.9, tie(0.9), 0.0, 128)]
for label, cases in (("CONTROL (passed before; must match out-prop4-large-t.txt)", controls), ("RERUN", failed)):
    print(label)
    for m, a5, k1, L5 in cases:
        r, tdw, tpv, ns, sp = ratios(L5, k1, m, a5)
        i = int(np.argmax(r))
        print("m=%.1f a5=%.4f k1=%.3f L5=%3d S=%6.1f: max %.4f at t=%5d ; top DW/r* %.6f PV/r* %.6f nsurf %d ;"
              " r(t) = %s ; spread %.2g, 1e4*spread %.2g" % (m, a5, k1, L5, L5 * a5, r[i], TS[i], tdw, tpv, ns,
                                                            " ".join("%.3g" % x for x in r), sp, 1e4 * sp), flush=True)
print("DONE")
