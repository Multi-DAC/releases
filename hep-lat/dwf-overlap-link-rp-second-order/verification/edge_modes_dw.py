"""(H4) check, owed item 4 of PROOF-T1-2-finiteL5.md v2: with an OPEN time boundary (2L rows, no wrap), are the free
D_DW (m_f = 0, open s) and D_PV (m_f = 1, antiperiodic s) gapped uniformly in L, and do their propagators at the
reflection plane converge as L grows?  Fixed spatial momentum k1; 2d gammas g0 = s1, g1 = s2, g5 = s3.
Setting as in the draft: each slice is the 2d Wilson operator (r = 1) at mass -m; between slices
(1/a5)(delta - P_- delta_{s+1} - P_+ delta_{s-1}); PV wall term (m_f/a5)(P_- |L5><1| + P_+ |1><L5|).
Float only; a diagnostic, not a proof."""
import sys
import numpy as np

g0 = np.array([[0, 1], [1, 0]], complex)
g1 = np.array([[0, -1j], [1j, 0]], complex)
g5 = np.array([[1, 0], [0, -1]], complex)
I2 = np.eye(2)
Pp, Pm = (I2 + g5) / 2, (I2 - g5) / 2


def op(nt, L5, k1, m, a5, mf, open_t=True):
    """Return the (2 * L5 * nt)-square operator. Index (t, s, spin) -> (t * L5 + s) * 2 + spin."""
    n = 2 * L5 * nt
    D = np.zeros((n, n), complex)

    def blk(t, s, t2, s2, M):
        i, j = (t * L5 + s) * 2, (t2 * L5 + s2) * 2
        D[i:i + 2, j:j + 2] += M

    onsite = (1 - np.cos(k1) + 1 - m) * I2 + 1j * np.sin(k1) * g1 + (1 / a5) * I2
    for t in range(nt):
        for s in range(L5):
            blk(t, s, t, s, onsite)
            if s + 1 < L5:
                blk(t, s, t, s + 1, -(1 / a5) * Pm)
                blk(t, s + 1, t, s, -(1 / a5) * Pp)
            if mf != 0:
                pass
            tn = t + 1
            if tn < nt or not open_t:
                tn %= nt
                blk(t, s, tn, s, -0.5 * (I2 - g0))
                blk(tn, s, t, s, -0.5 * (I2 + g0))
        if mf != 0:
            blk(t, L5 - 1, t, 0, (mf / a5) * Pm)
            blk(t, 0, t, L5 - 1, (mf / a5) * Pp)
    return D


def mid_block(Dinv, nt, L5):
    """Propagator from row c to row c + 1 (the two rows straddling the middle), all s and spin."""
    c = nt // 2 - 1
    i0, j0 = c * L5 * 2, (c + 1) * L5 * 2
    return Dinv[i0:i0 + 2 * L5, j0:j0 + 2 * L5]


def main():
    m = float(sys.argv[1]) if len(sys.argv) > 1 else 0.7
    a5 = float(sys.argv[2]) if len(sys.argv) > 2 else 1.0
    L5 = int(sys.argv[3]) if len(sys.argv) > 3 else 4
    print("m = %.3f  a5 = %.3f  L5 = %d" % (m, a5, L5))
    for k1 in (0.0, np.pi / 8, np.pi / 4, np.pi / 2):
        for name, mf in (("DW", 0.0), ("PV", 1.0)):
            ref = None
            rows = []
            for nt in (8, 16, 32, 64, 96):
                D = op(nt, L5, k1, m, a5, mf)
                s = np.linalg.svd(D, compute_uv=False)
                smin = s[-1]
                if smin < 1e-13:
                    rows.append((nt, smin, s[-2], np.nan))
                    continue
                G = mid_block(np.linalg.inv(D), nt, L5)
                rows.append((nt, smin, s[-2], G))
            # periodic-time reference with a long ring: translation invariant, no boundary
            Dp = op(192, L5, k1, m, a5, mf, open_t=False)
            sp = np.linalg.svd(Dp, compute_uv=False)[-1]
            Gp = mid_block(np.linalg.inv(Dp), 192, L5)
            print("k1 = %.4f  %s   periodic-192 sigma_min = %.3e" % (k1, name, sp))
            for nt, smin, s2, G in rows:
                if isinstance(G, float):
                    print("   2L = %3d  sigma_min %.3e  sigma_2 %.3e  (singular to float)" % (nt, smin, s2))
                else:
                    print("   2L = %3d  sigma_min %.3e  sigma_2 %.3e  max|G_open - G_ring| = %.3e" %
                          (nt, smin, s2, np.abs(G - Gp).max()))


if __name__ == "__main__":
    main()
