"""Independent check of Wren's letter 74 breaks v2-1 (a5 = 2/m: E_PV = 0 at k1 = 0, k5 = pi) and v2-2 (a5 = 2,
L5 = 1: open-time-BC boundary zero mode of D_DW and an O(1) artefact at the plane). Uses my own operator code."""
import sys
import numpy as np
sys.path.insert(0, ".")
from edge_modes_dw import op, mid_block


def Ginf(L5, k1, m, a5, mf, n=1 << 15):
    """Same as edge_modes_dw_exact.Ginf (copied: that module has no main guard)."""
    D3 = op(3, L5, k1, m, a5, mf, open_t=False)
    b = 2 * L5
    A0, A1, Am1 = D3[:b, :b], D3[:b, b:2 * b], D3[:b, 2 * b:3 * b]
    p = 2 * np.pi * (np.arange(n) + 0.5) / n
    Dp = A0[None] + A1[None] * np.exp(1j * p)[:, None, None] + Am1[None] * np.exp(-1j * p)[:, None, None]
    return (np.linalg.inv(Dp) * np.exp(-1j * p)[:, None, None]).mean(axis=0)


def bulk_smin(L5, k1, m, a5, mf, p0):
    D3 = op(3, L5, k1, m, a5, mf, open_t=False)
    b = 2 * L5
    A0, A1, Am1 = D3[:b, :b], D3[:b, b:2 * b], D3[:b, 2 * b:3 * b]
    Dp = A0 + A1 * np.exp(1j * p0) + Am1 * np.exp(-1j * p0)
    return np.linalg.svd(Dp, compute_uv=False)[-1]


m = 0.7
print("v2-1: a5 = 2/m, PV, k1 = 0, p0 = 0: bulk sigma_min by L5")
for L5 in (1, 2, 3, 4, 5):
    print("   L5 = %d  %.3e   (control a5 = 1: %.3e)" % (L5, bulk_smin(L5, 0.0, m, 2 / m, 1.0, 0.0),
                                                       bulk_smin(L5, 0.0, m, 1.0, 1.0, 0.0)))
print("v2-2: a5 = 2, L5 = 1, DW, k1 = 0: open sigma_min and plane-propagator error vs exact L = infinity")
Gi = Ginf(1, 0.0, m, 2.0, 0.0)
print("   bulk sigma_min (p0 = 0) %.3e" % bulk_smin(1, 0.0, m, 2.0, 0.0, 0.0))
for nt in (12, 24, 48, 96):
    D = op(nt, 1, 0.0, m, 2.0, 0.0)
    s = np.linalg.svd(D, compute_uv=False)[-1]
    G = mid_block(np.linalg.inv(D), nt, 1)
    print("   2L = %3d  sigma_min %.3e  max|G_open - G_inf| %.3e" % (nt, s, np.abs(G - Gi).max()))
print("control: a5 = 1, L5 = 1, DW (mass 1/a5 - m = +0.3)")
Gi = Ginf(1, 0.0, m, 1.0, 0.0)
for nt in (12, 24, 48):
    D = op(nt, 1, 0.0, m, 1.0, 0.0)
    G = mid_block(np.linalg.inv(D), nt, 1)
    print("   2L = %3d  sigma_min %.3e  max|G_open - G_inf| %.3e" % (nt, np.linalg.svd(D, compute_uv=False)[-1],
                                                                   np.abs(G - Gi).max()))
