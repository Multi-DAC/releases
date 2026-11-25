"""FLOAT, not certified. The 2x2 witness geometry (L = 4, time sites -3..4, Ls = 1, m = 0.7, eps = 0.2) at beta = 0,
widened to the full single-impurity family A_s: a1(t = s) = eps, s = 1..4 (the positive half).

Consultant read of PROOF-2x2.md (11-25), points 3 and 5:
  (5) after the Phi integral, K_st should depend only on the ring separation d = s + t - 1 of the two impurities
      (at 1 - s and t), so K is a Hankel matrix K_st = f(s + t - 1), with f(d) = f(8 - d), and the 2x2 witness reads
      v^T K v = 2 f(3) - 2 f(4) (K11 = f(3), K22 = f(5) = f(3), K12 = K21 = f(4)).
  (3) positive control in the same geometry and pipeline: the Wilson r = 1 determinant det X (X = D_W - m, the overlap
      kernel itself) and det(D_W + 0.1), both link-RP by Osterwalder-Seiler, must give a PSD 4 x 4 Gram.
Checks printed: Hankel residual of each K, f(d) for d = 1..7, eigenvalues of each 4 x 4 K, and the 2x2 v^T K v
against the certified enclosure [-8.264675148e-9, -8.232272695e-9].
Usage: fd_table.py [n_phi]
"""
import sys
import numpy as np
import mirror_test as mt
import mirror_kernel as mk
import certify_2x2 as c2

L, LS, EPS = c2.L, c2.LS, c2.EPS


def fam():
    out = []
    for s in range(1, L + 1):
        a0 = np.zeros((L - 1, LS)); a1 = np.zeros((L, LS))
        a1[s - 1, :] = EPS
        out.append((a0, a1))
    return out


def det_wilson(lat, U0, U1, phi, mass):
    """det(D_W + mass) with the same distributed twist as c2.X_of (mass = -0.7 gives det X)."""
    f = np.exp(1j * phi / (2 * L)); b = np.exp(-1j * phi / (2 * L))
    D = c2.wilson_tw(lat, U0, U1, f, b) + mass * np.eye(2 * lat.n)
    return np.linalg.det(D)


def gram(lat, F, fn, n):
    phis = 2 * np.pi * np.arange(n) / n
    K = np.zeros((len(F), len(F)), complex)
    for i in range(len(F)):
        for j in range(len(F)):
            U0, U1 = mk.pair_links(lat, F[i], F[j])
            K[i, j] = np.mean([fn(lat, U0, U1, ph) for ph in phis])
    return K


def report(name, K):
    print("== %s" % name)
    herm = np.max(np.abs(K - K.conj().T)); imag = np.max(np.abs(K.imag))
    f = {}
    for s in range(1, L + 1):
        for t in range(1, L + 1):
            f.setdefault(s + t - 1, []).append(K[s - 1, t - 1].real)
    hank = max(max(v) - min(v) for v in f.values())
    print("  max |K - K^dag| %.2e ; max |Im K| %.2e ; Hankel residual (spread of K_st at fixed s+t-1) %.2e" % (herm, imag, hank))
    print("  f(d), d = 1..7: " + "  ".join("%d:%.13e" % (d, np.mean(f[d])) for d in sorted(f)))
    print("  f(d) - f(8-d), d = 1..3: " + "  ".join("%.2e" % (np.mean(f[d]) - np.mean(f[8 - d])) for d in (1, 2, 3)))
    ev = np.linalg.eigvalsh(0.5 * (K + K.conj().T))
    print("  eigenvalues of the 4x4 K: " + "  ".join("%.6e" % e for e in ev))
    print("  2x2 witness v^T K v on {A_2, A_3} (impurity at t = 2, 3): %.10e" % (K[1, 1] + K[2, 2] - K[1, 2] - K[2, 1]).real)
    print("  = 2 f(3) - 2 f(4) = %.10e" % (2 * np.mean(f[3]) - 2 * np.mean(f[4])))
    return ev


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 256
    lat = mt.Lat(L, LS); F = fam()
    print("n_phi = %d (periodic trapezoid; the certified run used n = 124)" % n)
    Kov = gram(lat, F, lambda lat, U0, U1, ph: c2.det_D(lat, U0, U1, ph)[0], n)
    report("overlap det D_ov, m = 0.7", Kov)
    print("  certified 2x2 enclosure: [-8.264675148e-9, -8.232272695e-9]")
    for mass, name in ((-c2.M, "control: Wilson det X = det(D_W - 0.7)"), (0.1, "control: Wilson det(D_W + 0.1)")):
        ev = report(name, gram(lat, F, lambda lat, U0, U1, ph, ms=mass: det_wilson(lat, U0, U1, ph, ms), n))
        print("  control verdict: %s (min eigenvalue / max = %.2e)" % ("PSD" if ev[0] >= -1e-12 * ev[-1] else "NOT PSD", ev[0] / ev[-1]))
    print("DONE")


if __name__ == "__main__":
    main()
