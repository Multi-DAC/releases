"""T3c (after clock reading 16:04): letter78_t3b.py again, but with poles only. t3b went past 24 GB at L5 = 640 and was
stopped at 16:04: Prop keeps a dense (2 L5)^2 residue for every one of its ~4 L5 poles, about 256 L5^3 bytes (67 GB at
L5 = 640). Here only the pencil's eigenvalues are computed, so memory is O(L5^2). Control: the S = 8 rows must
reproduce T3 (m = 0.3, k1 = pi/16: 0.818653 at a5 = 0.2 ... 0.818065 at a5 = 0.1)."""
import numpy as np
import scipy.linalg as sla

from l4_assembled import blocks, rstar, cluster_idx
from letter78_tests import lam0_ov


def surface_poles(L5, k1, m, a5):
    B0, Bp, Bm = blocks(L5, k1, m, a5, 0.0)
    n = 2 * L5
    Z, I = np.zeros((n, n)), np.eye(n)
    w = sla.eigvals(np.block([[Z, I], [-Bm, -B0]]), np.block([[I, Z], [Z, Bp]]))
    w = w[np.isfinite(w) & (np.abs(w) > 1e-12)]
    w = w.real[(w.real > rstar(k1, m) * (1 + 1e-12)) & (w.real < 1)]
    return [w[np.asarray(g)].mean() for g in cluster_idx(w)] if len(w) else []


print("T3c m=0.3: DW surface pole against S, poles only", flush=True)
for k1 in (np.pi / 16, np.pi / 8, np.pi / 4):
    l0 = lam0_ov(k1, 0.3)
    for a5 in (0.2, 0.1, 0.05):
        row = []
        for S in (8, 16, 32, 64):
            L5 = int(round(S / a5))
            if L5 > 640:
                continue
            s = surface_poles(L5, k1, 0.3, a5)
            row.append("S=%d: %s" % (S, "%.6f (gap %.2e)" % (s[0], l0 - s[0]) if len(s) == 1 else "n=%d" % len(s)))
        print("  k1=%.3f a5=%.2f overlap pole %.6f ; %s" % (k1, a5, l0, "; ".join(row)), flush=True)
print("DONE")
