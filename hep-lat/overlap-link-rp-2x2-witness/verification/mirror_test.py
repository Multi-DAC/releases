"""Mirror-symmetric test of reflection positivity for the overlap Dirac fermion in a U(1) background, d = 2.

Spec: Research/kikukawa-usui-reflection-positivity-free-overlap-2026-11-24.md, "A checkable necessary condition".
Lattice t in {-L+1..L} (antiperiodic), x in {0..Ls-1} (periodic). Link reflection theta(t,x) = (1-t, x).
One-particle forms (Kikukawa-Usui L250-266):
  K1[(x,a),(y,b)] = [S(y, theta x) g0]_{b a}      from F = sum f psi
  K2[(x,a),(y,b)] = [g0 S(theta x, y)]_{a b}      from F = sum g psibar
with S = D^-1, x,y at positive times. Overall signs are fixed by the theorem at U = 1.
Usage: mirror_test.py [L] [Ls] [m]
"""
import sys
import numpy as np

s1 = np.array([[0, 1], [1, 0]], complex); s2 = np.array([[0, -1j], [1j, 0]]); s3 = np.array([[1, 0], [0, -1]], complex)
g0, g1, g5 = s1, s2, s3
I2 = np.eye(2)


class Lat:
    def __init__(self, L, Ls):
        self.L, self.Ls = L, Ls
        self.T = 2 * L
        self.ts = list(range(-L + 1, L + 1))
        self.n = self.T * Ls

    def idx(self, t, x):
        return ((t + self.L - 1) % self.T) * self.Ls + (x % self.Ls)

    def theta_perm(self):
        """site -> reflected site"""
        p = np.zeros(self.n, int)
        for t in self.ts:
            for x in range(self.Ls):
                p[self.idx(t, x)] = self.idx(1 - t, x)
        return p

    def plus_sites(self):
        return [self.idx(t, x) for t in range(1, self.L + 1) for x in range(self.Ls)]


def wilson(lat, U0, U1, r=1.0):
    """D_w with U(1) links. U0[t_index, x]: link (t,x)->(t+1,x); U1: (t,x)->(t,x+1). Antiperiodic in t."""
    n = lat.n
    D = np.zeros((2 * n, 2 * n), complex)
    for t in lat.ts:
        for x in range(lat.Ls):
            i = lat.idx(t, x)
            D[2*i:2*i+2, 2*i:2*i+2] += 2 * r * I2
            # time direction
            j = lat.idx(t + 1, x); sgn = -1.0 if t == lat.L else 1.0
            u = U0[t + lat.L - 1, x] * sgn
            D[2*i:2*i+2, 2*j:2*j+2] += -0.5 * (r * I2 - g0) * u
            D[2*j:2*j+2, 2*i:2*i+2] += -0.5 * (r * I2 + g0) * np.conj(u)
            # space direction
            j = lat.idx(t, x + 1)
            u = U1[t + lat.L - 1, x]
            D[2*i:2*i+2, 2*j:2*j+2] += -0.5 * (r * I2 - g1) * u
            D[2*j:2*j+2, 2*i:2*i+2] += -0.5 * (r * I2 + g1) * np.conj(u)
    return D


def overlap(lat, U0, U1, m, r=1.0):
    X = wilson(lat, U0, U1, r) - m * np.eye(2 * lat.n)
    w, v = np.linalg.eigh(X.conj().T @ X)
    return 0.5 * (np.eye(2 * lat.n) + X @ (v @ np.diag(w ** -0.5) @ v.conj().T)), w.min()


def symmetric_links(lat, a0, a1, mode):
    """Build a reflection-symmetric configuration from phases given at positive times.
    a0[t-1, x] for the links (t,x)->(t+1,x), t = 1..L-1; a1[t-1, x] for spatial links at t = 1..L.
    Crossing links (t=0 and the wrap at t=L) are 1. mode picks the candidate rule for the image:
      'A': U1(1-t,x) = U1(t,x),        U0(-t,x) = conj U0(t,x)
      'B': U1(1-t,x) = conj U1(t,x),   U0(-t,x) = U0(t,x)
    The right rule is the one for which g0 D(theta x, theta y) g0 = D(x,y)^dagger; see sym_residual."""
    L, Ls = lat.L, lat.Ls
    U0 = np.ones((2 * L, Ls), complex); U1 = np.ones((2 * L, Ls), complex)
    for t in range(1, L + 1):
        for x in range(Ls):
            u = np.exp(1j * a1[t - 1, x])
            U1[t + L - 1, x] = u
            U1[(1 - t) + L - 1, x] = u if mode == 'A' else np.conj(u)
    for t in range(1, L):
        for x in range(Ls):
            u = np.exp(1j * a0[t - 1, x])
            U0[t + L - 1, x] = u
            U0[(-t) + L - 1, x] = np.conj(u) if mode == 'A' else u
    return U0, U1


def spin_perm(lat, D):
    p = lat.theta_perm()
    pp = np.repeat(2 * p, 2) + np.tile([0, 1], lat.n)
    return D[np.ix_(pp, pp)]


def sym_residual(lat, D):
    """max | g0 D(theta x, theta y) g0 - D(x,y)^dagger |"""
    G = np.kron(np.eye(lat.n), g0)
    return np.abs(G @ spin_perm(lat, D) @ G - D.conj().T).max()


def forms(lat, D):
    S = np.linalg.inv(D)
    p = lat.theta_perm()
    plus = np.array(lat.plus_sites())
    ip = np.repeat(2 * plus, 2) + np.tile([0, 1], len(plus))
    it = np.repeat(2 * p[plus], 2) + np.tile([0, 1], len(plus))
    Gp = np.kron(np.eye(len(plus)), g0)
    K1 = (S[np.ix_(ip, it)] @ Gp).T
    K2 = Gp @ S[np.ix_(it, ip)]
    return K1, K2


def spec(K):
    herm = np.abs(K - K.conj().T).max()
    ev = np.linalg.eigvalsh(0.5 * (K + K.conj().T))
    return herm, ev.min(), ev.max()


def report(tag, lat, D, sign):
    K1, K2 = forms(lat, D)
    out = []
    for name, K, s in (("K1", K1, sign[0]), ("K2", K2, sign[1])):
        h, lo, hi = spec(s * K)
        out.append("%s herm %.1e min %+.3e max %+.3e" % (name, h, lo, hi))
    print("%-34s symres %.1e | %s" % (tag, sym_residual(lat, D), " | ".join(out)))


def plaq_max(lat, U0, U1):
    L, Ls = lat.L, lat.Ls
    worst = 0.0
    for ti in range(2 * L):
        for x in range(Ls):
            p = U0[ti, x] * U1[(ti + 1) % (2 * L), x] * np.conj(U0[ti, (x + 1) % Ls]) * np.conj(U1[ti, x])
            worst = max(worst, abs(1 - p))
    return worst


def main():
    L = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    Ls = int(sys.argv[2]) if len(sys.argv) > 2 else 6
    m = float(sys.argv[3]) if len(sys.argv) > 3 else 1.0
    lat = Lat(L, Ls)
    rng = np.random.default_rng(7)
    print("lattice T=%d Ls=%d, m=%g, dim %d; forms are %dx%d" % (2 * L, Ls, m, 2 * lat.n, 2 * L * Ls, 2 * L * Ls))
    z0 = np.zeros((L - 1, Ls)); z1 = np.zeros((L, Ls))

    # --- control (a): U = 1. Fixes the overall signs; must be one-signed (Kikukawa-Usui).
    U0, U1 = symmetric_links(lat, z0, z1, 'A')
    D, gap = overlap(lat, U0, U1, m)
    K1, K2 = forms(lat, D)
    sign = []
    for K in (K1, K2):
        ev = np.linalg.eigvalsh(0.5 * (K + K.conj().T))
        sign.append(1.0 if abs(ev.max()) >= abs(ev.min()) else -1.0)
    print("signs fixed at U=1:", sign)
    report("(a) overlap, U=1", lat, D, sign)
    Dw = wilson(lat, U0, U1) + 0.3 * np.eye(2 * lat.n)
    Kw = forms(lat, Dw); wsign = []
    for K in Kw:
        ev = np.linalg.eigvalsh(0.5 * (K + K.conj().T)); wsign.append(1.0 if abs(ev.max()) >= abs(ev.min()) else -1.0)
    report("(a) Wilson r=1 M=0.3, U=1", lat, Dw, wsign)
    Dw2 = wilson(lat, U0, U1, r=2.0) + 0.3 * np.eye(2 * lat.n)
    report("(b) Wilson r=2 M=0.3, U=1", lat, Dw2, wsign)

    # --- which image rule is the symmetric one?
    a0 = rng.normal(size=(L - 1, Ls)); a1 = rng.normal(size=(L, Ls))
    for mode in "AB":
        U0, U1 = symmetric_links(lat, 0.3 * a0, 0.3 * a1, mode)
        print("rule %s: symres Wilson %.2e" % (mode, sym_residual(lat, wilson(lat, U0, U1))))

    # --- scan
    for mode in "AB":
        print("== rule", mode)
        for amp in (0.01, 0.03, 0.1, 0.2, 0.3, 0.5, 0.8):
            U0, U1 = symmetric_links(lat, amp * a0, amp * a1, mode)
            D, gap = overlap(lat, U0, U1, m)
            tag = "amp %.2f |1-P|max %.2f gap %.1e" % (amp, plaq_max(lat, U0, U1), gap)
            report("ov " + tag[:30], lat, D, sign)
            report("W1 " + tag[:30], lat, wilson(lat, U0, U1) + 0.3 * np.eye(2 * lat.n), wsign)
        # non-symmetric control: break the symmetry on the negative side
        U0, U1 = symmetric_links(lat, 0.3 * a0, 0.3 * a1, mode)
        U1[0, 0] *= np.exp(0.4j)
        D, gap = overlap(lat, U0, U1, m)
        report("(b) ov, symmetry broken", lat, D, sign)


if __name__ == "__main__":
    main()
