"""Orbit-averaged determinant kernel: the gauge-INVARIANT form of the pair test in mirror_kernel.py.

mirror_kernel.py fixes the crossing time-links V to 1 and tests M_ij = det D[theta A_i, 1, A_j]. That is positivity
on the whole positive-time algebra in the gauge V = 1, which is more than a gauge-invariant observable needs.
For F = sum_i c_i delta_orbit(U_+ ~ A_i) (gauge invariant), integrating the crossing links gives
    <theta(F) F>  ~  sum_ij conj(c_i) c_j  Mt_ij,
    Mt_ij = int dV  h_ij(V) det D[theta A_i, V, A_j],
    h_ij(V) = exp(beta * sum of Re(crossing plaquettes)).
Mt must be Hermitian and positive semidefinite. The V integral is done on an n^(2 Ls) grid (exact for Wilson only at beta = 0, where the integrand is a
low-degree trigonometric polynomial in V; for beta > 0 and for the overlap, run several n and compare).
Controls: Wilson r = 1 (theorem: must pass), Wilson r = 2 (fails in this test; no theorem says it must;
refuter round 11-24: r = 0.5, 0.9, 1.1 pass).
Usage: mirror_orbit.py [L] [Ls] [N] [amp] [m] [n] [seed]
"""
import sys
import itertools
import numpy as np
import mirror_test as mt
import mirror_kernel as mk


def cross_rows(lat):
    return (lat.L - 1, 2 * lat.L - 1)        # row index of the links t=0->1 and t=L->-L+1


def decompose(lat, U0, U1, r):
    """wilson(U) = W0 + sum_k v_k P_k + conj(v_k) Q_k over the 2 Ls crossing links."""
    rows = cross_rows(lat)
    Z0 = U0.copy()
    for row in rows:
        Z0[row, :] = 0
    W0 = mt.wilson(lat, Z0, U1, r)
    P = []; Q = []
    for row in rows:
        for x in range(lat.Ls):
            A = Z0.copy(); A[row, x] = 1.0
            B = Z0.copy(); B[row, x] = 1.0j
            d1 = mt.wilson(lat, A, U1, r) - W0      # P + Q
            di = mt.wilson(lat, B, U1, r) - W0      # iP - iQ
            P.append(0.5 * (d1 - 1j * di)); Q.append(0.5 * (d1 + 1j * di))
    return W0, np.array(P), np.array(Q)


def dets(lat, U0, U1, V, kind, m):
    """V: (nV, 2 Ls) complex phases. Returns det for each V."""
    r = 2.0 if kind == "w2" else 1.0
    W0, P, Q = decompose(lat, U0, U1, r)
    W = W0[None] + np.einsum('vk,kab->vab', V, P) + np.einsum('vk,kab->vab', V.conj(), Q)
    dim = W.shape[1]
    if kind == "ov":
        X = W - m * np.eye(dim)[None]
        w, v = np.linalg.eigh(np.einsum('vba,vbc->vac', X.conj(), X))
        inv = np.einsum('vab,vb,vcb->vac', v, w ** -0.5, v.conj())
        D = 0.5 * (np.eye(dim)[None] + X @ inv)
        return np.linalg.det(D), w[:, 0]
    return np.linalg.det(W + 0.3 * np.eye(dim)[None]), np.full(len(V), np.nan)


def kernel(lat, fields, kind, m, n, beta):
    Ls = lat.Ls
    ang = 2 * np.pi * np.arange(n) / n
    grid = np.array(list(itertools.product(ang, repeat=2 * Ls)))      # (nV, 2 Ls): plane 0 then plane L
    V = np.exp(1j * grid)
    N = len(fields)
    M = np.zeros((N, N), complex); gap = np.inf; bad = 0.0
    for i in range(N):
        for j in range(N):
            U0, U1 = mk.pair_links(lat, fields[i], fields[j])
            d, g = dets(lat, U0, U1, V, kind, m)
            gap = min(gap, g.min()) if kind == "ov" else gap
            wt = np.ones(len(grid))
            if beta:
                # crossing plaquette at plane 0, column x: U1(0,x) V(x+1) conj U1(1,x) conj V(x); U1(0,.) is A_i row 1
                ai, bj = fields[i][1], fields[j][1]
                v0 = grid[:, :Ls]; vL = grid[:, Ls:]
                s0 = np.cos(ai[0][None] - bj[0][None] + np.roll(v0, -1, axis=1) - v0).sum(axis=1)
                # plane L: links t=L -> -L+1; U1(L,x) is B_j last row, U1(-L+1,x) is A_i last row
                sL = np.cos(bj[-1][None] - ai[-1][None] + np.roll(vL, -1, axis=1) - vL).sum(axis=1)
                wt = np.exp(beta * (s0 + sL))
                d = d * wt
            M[i, j] = d.mean()
            if kind == "ov":
                bad = max(bad, wt[g < 0.05].sum() / wt.sum())
    return M, (gap, bad)


def main():
    a = sys.argv
    L = int(a[1]) if len(a) > 1 else 2
    Ls = int(a[2]) if len(a) > 2 else 2
    N = int(a[3]) if len(a) > 3 else 12
    amp = float(a[4]) if len(a) > 4 else 1.0
    m = float(a[5]) if len(a) > 5 else 0.7
    n = int(a[6]) if len(a) > 6 else 6
    seed = int(a[7]) if len(a) > 7 else 5
    lat = mt.Lat(L, Ls)
    rng = np.random.default_rng(seed)
    shared = len(a) > 8 and a[8] == 'shared'
    betas = [float(b) for b in a[9].split(',')] if len(a) > 9 else [0.0, 1.0, 4.0]
    fields = [(amp * rng.normal(size=(L - 1, Ls)), amp * rng.normal(size=(L, Ls))) for _ in range(N)]
    if shared:      # all fields share the spatial links on the two boundary rows, so h does not depend on (i, j)
        for f in fields[1:]:
            f[1][0] = fields[0][1][0]; f[1][-1] = fields[0][1][-1]
    print("shared boundary rows: %s" % shared)
    print("L=%d Ls=%d N=%d amp=%.2f m=%.2f grid n=%d (%d points) seed=%d" % (L, Ls, N, amp, m, n, n ** (2 * Ls), seed))
    for kind, name in (("w1", "Wilson r=1 (must pass)"), ("w2", "Wilson r=2 (must fail)"), ("ov", "overlap")):
        for beta in betas:
            M, gap = kernel(lat, fields, kind, m, n, beta)
            herm, ev = mk.spectrum(mk.normalised(M))
            print("%-24s beta=%-3g herm %.1e  min %+.3e  max %+.3e  (#ev < -1e-10: %d)%s"
                  % (name, beta, herm, ev.min(), ev.max(), int((ev < -1e-10).sum()),
                     "" if kind != "ov" else "  smallest X^dX eig on grid %.4f, worst weight fraction with eig < 0.05: %.1e" % gap))
            sys.stdout.flush()


if __name__ == "__main__":
    main()
