"""Check of Wren's letter 75, (H4) for 0 < a5 < 1/m, with my own operator (edge_modes_dw.op) and my own exact
L = infinity propagator (p0 integral), not her script.

A. Step 1 and Step 2 objects built from MY op's blocks: the bulk transfer matrix T from the per-row blocks of the
   open-time operator, checked (i) Hermitian positive definite, (ii) its plane-wave solutions solve my op, (iii) the
   determinant identity det(a5 D^(L)) = (det M0)^(2L) det(bottom block of T^(2L)).
B. Threshold: plane-propagator error |G_open - G_inf| at the two rows straddling the middle, and sigma_min of the open
   operator, against 2L, at a5 below 1/m (1.3, 1.4) and above (1.45, 1.6, 2.0), m = 0.7, L5 = 1, 2, 4, k1 = 0 and 0.3.
   Above 1/m I expect an edge mode at k1 = 0 (sigma_min collapsing exponentially in 2L); below, convergence.
Float only."""
import sys
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spl

sys.path.insert(0, ".")
from edge_modes_dw import op, g0, g5, Pp, Pm, I2  # noqa: E402


def blocks(L5, k1, m, a5, mf):
    D3 = op(3, L5, k1, m, a5, mf, open_t=True)
    b = 2 * L5
    return D3[b:2 * b, b:2 * b], D3[b:2 * b, 2 * b:3 * b], D3[b:2 * b, :b]   # B0, B+ (to t+1), B- (to t-1)


def Ginf(L5, k1, m, a5, mf, n=1 << 15):
    B0, Bp, Bm = blocks(L5, k1, m, a5, mf)
    p = 2 * np.pi * (np.arange(n) + 0.5) / n
    Dp = B0[None] + Bp[None] * np.exp(1j * p)[:, None, None] + Bm[None] * np.exp(-1j * p)[:, None, None]
    return (np.linalg.inv(Dp) * np.exp(-1j * p)[:, None, None]).mean(axis=0)   # G(t, t+1)


def sparse_op(nt, L5, k1, m, a5, mf):
    B0, Bp, Bm = blocks(L5, k1, m, a5, mf)
    b = 2 * L5
    return sp.kron(sp.eye(nt), sp.csr_matrix(B0)) + sp.kron(sp.eye(nt, k=1), sp.csr_matrix(Bp)) \
        + sp.kron(sp.eye(nt, k=-1), sp.csr_matrix(Bm)), b


def mid_G(nt, L5, k1, m, a5, mf):
    D, b = sparse_op(nt, L5, k1, m, a5, mf)
    D = D.tocsc()
    c = nt // 2 - 1
    rhs = np.zeros((nt * b, b), complex)
    rhs[(c + 1) * b:(c + 2) * b, :] = np.eye(b)
    X = spl.splu(D).solve(rhs)
    return X[c * b:(c + 1) * b, :]


def smin_open(nt, L5, k1, m, a5, mf):
    D, b = sparse_op(nt, L5, k1, m, a5, mf)
    if nt * b <= 1200:
        return np.linalg.svd(D.toarray(), compute_uv=False)[-1]
    lu = spl.splu(D.tocsc())
    luH = spl.splu(D.conj().T.tocsc())
    n = nt * b
    Aop = spl.LinearOperator((n, n), matvec=lambda x: lu.solve(luH.solve(x)), dtype=complex)  # (D D^H)^-1
    ev = spl.eigsh(Aop, k=1, which="LM", return_eigenvectors=False, tol=1e-10)
    return 1 / np.sqrt(ev[0].real)


def transfer(L5, k1, m, a5, mf):
    """T from my blocks, in the gamma0 eigenbasis: v_t = (u+(t-1), u-(t)). Generic elimination, no formula from her."""
    B0, Bp, Bm = blocks(L5, k1, m, a5, mf)
    w, U = np.linalg.eigh(g0)          # columns: gamma0 = -1, +1
    Wm, Wp = U[:, [0]], U[:, [1]]
    E = np.eye(L5)
    Pplus = np.kron(E, Wp)             # (2 L5) x L5 : u+ coordinates  (index (s, spin))
    Pminus = np.kron(E, Wm)
    C = np.hstack([Pplus, Pminus])     # orthonormal change of basis, psi = C (u+, u-)
    A0 = C.conj().T @ B0 @ C
    Ap = C.conj().T @ Bp @ C
    Am = C.conj().T @ Bm @ C
    n = L5
    # equation: Am x(t-1) + A0 x(t) + Ap x(t+1) = 0, x = (u+, u-). Check hop structure.
    hop_ok = (np.abs(Ap[:, :n]).max(), np.abs(Am[:, n:]).max())
    # rows: the u+ rows (first n) and u- rows. Ap acts on u-(t+1) only; Am on u+(t-1) only.
    # Unknowns at step t: given u+(t-1), u-(t), solve u+(t), u-(t+1):
    # A0[:, :n] u+(t) + Ap[:, n:] u-(t+1) = -Am[:, :n] u+(t-1) - A0[:, n:] u-(t)
    Lhs = np.hstack([A0[:, :n], Ap[:, n:]])
    Rhs = -np.hstack([Am[:, :n], A0[:, n:]])
    T = np.linalg.solve(Lhs, Rhs)
    return T, hop_ok, A0, C


def partA():
    m = 0.7
    print("A. transfer matrix from my blocks (gamma0 basis), checks")
    worst = dict(herm=0.0, det=0.0)
    minev = 1e9
    for mf, name in ((0.0, "DW"), (1.0, "PV")):
        for L5 in (1, 2, 3):
            for k1 in (0.0, 0.7, 2.5):
                for a5 in (0.5, 1.0, 1.3, 1.42):
                    T, hop_ok, A0, C = transfer(L5, k1, m, a5, mf)
                    assert max(hop_ok) < 1e-14, hop_ok
                    herm = np.abs(T - T.conj().T).max() / np.abs(T).max()
                    ev = np.linalg.eigvalsh((T + T.conj().T) / 2)
                    minev = min(minev, ev.min())
                    worst["herm"] = max(worst["herm"], herm)
                    # determinant identity for 2L = 6
                    nt = 6
                    D = op(nt, L5, k1, m, a5, mf)
                    lhs = np.linalg.det(a5 * D)
                    M0 = a5 * A0[:L5, :L5]          # u+u+ block of a5 * B0
                    Tn = np.linalg.matrix_power(T, nt)
                    rhs = a5 ** (nt * L5) * np.linalg.det(M0) ** nt * np.linalg.det(Tn[L5:, L5:])   # a5^(2L L5): absent from letter 75
                    rel = abs(lhs - rhs) / abs(lhs)
                    worst["det"] = max(worst["det"], rel)
                    if rel > 1e-8 or abs(lhs.imag) > 1e-8 * abs(lhs) or lhs.real <= 0:
                        print("   FLAG %s L5=%d k1=%.2f a5=%.2f det lhs %s rhs %s" % (name, L5, k1, a5, lhs, rhs))
    print("   max rel |T - T^H| = %.2e ; min eigenvalue of T over the grid = %.3e ; max rel det-identity error = %.2e"
          % (worst["herm"], minev, worst["det"]))
    # plane-wave check: eigenpairs of T give solutions of my bulk equation
    T, _, _, C = transfer(2, 0.7, m, 1.3, 0.0)
    lam, V = np.linalg.eig(T)
    B0, Bp, Bm = blocks(2, 0.7, m, 1.3, 0.0)
    res = 0.0
    for j in range(len(lam)):
        v = V[:, j]
        # v_t = lam^t v with v_t = (u+(t-1), u-(t))  =>  psi(t) = lam^t C (lam v_top, v_bot)
        phi = C @ np.concatenate([lam[j] * v[:2], v[2:]])
        r = (Bm + lam[j] * B0 + lam[j] ** 2 * Bp) @ phi
        res = max(res, np.abs(r).max() / np.abs(phi).max() / max(1, abs(lam[j]) ** 2))
    print("   plane-wave residual (DW, L5=2, k1=0.7, a5=1.3), eigenpairs of T into my op: %.2e" % res)


def partB():
    m = 0.7
    print("B. threshold 1/m = %.6f: plane-propagator error and open sigma_min vs 2L" % (1 / m))
    for L5 in (1, 2, 4):
        for k1 in (0.0, 0.3):
            for mf, name in ((0.0, "DW"), (1.0, "PV")):
                for a5 in (1.3, 1.4, 1.45, 1.6, 2.0):
                    Gi = Ginf(L5, k1, m, a5, mf)
                    T, _, _, _ = transfer(L5, k1, m, a5, mf)
                    lam = np.sort(np.abs(np.linalg.eigvals(T)))
                    lmax = lam[lam < 1].max() if (lam < 1).any() else np.nan
                    row = []
                    for nt in (16, 64, 256, 1024):
                        e = np.abs(mid_G(nt, L5, k1, m, a5, mf) - Gi).max()
                        s = smin_open(nt, L5, k1, m, a5, mf) if nt <= 256 else np.nan
                        row.append("%4d: err %.1e smin %.1e" % (nt, e, s))
                    print("   L5=%d k1=%.1f %s a5=%.2f lam_max<1 %.5f | %s" % (L5, k1, name, a5, lmax, " | ".join(row)))
                    sys.stdout.flush()


if __name__ == "__main__":
    partA()
    partB()
