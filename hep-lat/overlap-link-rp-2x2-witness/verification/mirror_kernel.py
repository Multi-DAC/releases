"""Kernel test: reflection positivity of the fermion determinant across PAIRS of gauge fields, d = 2, U(1).

For observables F = sum_i c_i delta(U_+ - A_i) (no fermion fields), <theta(F) F> is, up to positive diagonal
factors, sum_ij conj(c_i) c_j  g(A_i, A_j) det D[U_+ = A_j, U_- = theta A_i], g the crossing-plaquette factor
of the gauge action. So the matrix
    M_ij = det D[A_j, theta A_i]            (and its Schur product with g)
must be Hermitian and positive semidefinite for every finite set of fields A_i at positive times.
(Corrected 11-24: that holds only after the crossing links V are integrated with their weight, as in
mirror_orbit.py and mirror_flat.py; the V = 1 matrix computed here is not itself required to be PSD.)
The diagonal i = j is the mirror-symmetric case. Theorem for Wilson r = 1 (Osterwalder-Seiler; Luscher): used as
the rounding floor. Wilson r = 2 is the control that must fail.
With one fermion field each side, the same with M_ij (x) K(A_j, theta A_i): the block matrix test.
Usage: mirror_kernel.py [L] [Ls] [N] [amp] [m]
"""
import sys
import numpy as np
import mirror_test as mt


def pair_links(lat, A, B):
    """U_+ from B = (b0, b1); U_- = image of A = (a0, a1) under rule A. Crossing links are 1."""
    L, Ls = lat.L, lat.Ls
    a0, a1 = A; b0, b1 = B
    U0 = np.ones((2 * L, Ls), complex); U1 = np.ones((2 * L, Ls), complex)
    for t in range(1, L + 1):
        U1[t + L - 1] = np.exp(1j * b1[t - 1])
        U1[(1 - t) + L - 1] = np.exp(1j * a1[t - 1])
    for t in range(1, L):
        U0[t + L - 1] = np.exp(1j * b0[t - 1])
        U0[(-t) + L - 1] = np.exp(-1j * a0[t - 1])
    return U0, U1


def gauge_cross(lat, A, B, beta):
    """exp(beta * sum of Re(crossing plaquettes)): rows t=0->1 and the wrap t=L->-L+1 (time links there are 1)."""
    s = np.sum(np.cos(B[1][0] - A[1][0])) + np.sum(np.cos(B[1][-1] - A[1][-1]))
    return np.exp(beta * s)


def normalised(M):
    d = np.sqrt(np.abs(np.real(np.diag(M))))
    return M / np.outer(d, d)


def spectrum(M):
    herm = np.abs(M - M.conj().T).max() / np.abs(M).max()
    ev = np.linalg.eigvalsh(0.5 * (M + M.conj().T))
    return herm, ev


def main():
    L = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    Ls = int(sys.argv[2]) if len(sys.argv) > 2 else 4
    N = int(sys.argv[3]) if len(sys.argv) > 3 else 24
    amp = float(sys.argv[4]) if len(sys.argv) > 4 else 0.3
    m = float(sys.argv[5]) if len(sys.argv) > 5 else 0.7
    lat = mt.Lat(L, Ls)
    rng = np.random.default_rng(5)
    fields = [(amp * rng.normal(size=(L - 1, Ls)), amp * rng.normal(size=(L, Ls))) for _ in range(N)]
    nk = 2 * L * Ls
    print("L=%d Ls=%d N=%d amp=%.2f m=%.2f; determinant kernel %dx%d, block kernel %dx%d"
          % (L, Ls, N, amp, m, N, N, N * nk, N * nk))

    ops = {
        "overlap": lambda U0, U1: mt.overlap(lat, U0, U1, m)[0],
        "Wilson r=1 (floor)": lambda U0, U1: mt.wilson(lat, U0, U1) + 0.3 * np.eye(2 * lat.n),
        "Wilson r=2 (must fail)": lambda U0, U1: mt.wilson(lat, U0, U1, r=2.0) + 0.3 * np.eye(2 * lat.n),
    }
    for name, op in ops.items():
        Mdet = np.zeros((N, N), complex)
        blocks1 = [[None] * N for _ in range(N)]
        for i in range(N):
            for j in range(N):
                U0, U1 = pair_links(lat, fields[i], fields[j])
                D = op(U0, U1)
                sgn, logdet = np.linalg.slogdet(D)
                Mdet[i, j] = sgn * np.exp(logdet)
                K1, K2 = mt.forms(lat, D)
                blocks1[i][j] = K1
        # determinant kernel, alone and with the gauge crossing factor
        for beta in (None, 1.0, 4.0):
            M = Mdet.copy()
            if beta is not None:
                G = np.array([[gauge_cross(lat, fields[i], fields[j], beta) for j in range(N)] for i in range(N)])
                M = M * G
            herm, ev = spectrum(normalised(M))
            print("%-24s det kernel %-9s herm %.1e  min %+.3e  max %+.3e  (#ev < -1e-10: %d)"
                  % (name, "" if beta is None else "beta=%g" % beta, herm, ev.min(), ev.max(), int((ev < -1e-10).sum())))
        # block kernel: one fermion field on each side
        sign = None
        big = np.zeros((N * nk, N * nk), complex)
        for i in range(N):
            for j in range(N):
                big[i * nk:(i + 1) * nk, j * nk:(j + 1) * nk] = Mdet[i, j] / np.sqrt(abs(Mdet[i, i] * Mdet[j, j])) * blocks1[i][j]
        herm, ev = spectrum(big)
        if abs(ev.min()) > abs(ev.max()):
            ev = -ev[::-1]
        print("%-24s block kernel (det x K1)  herm %.1e  min %+.3e  max %+.3e  (#ev < -1e-10: %d)"
              % (name, herm, ev.min(), ev.max(), int((ev < -1e-10).sum())))
        sys.stdout.flush()


if __name__ == "__main__":
    main()
