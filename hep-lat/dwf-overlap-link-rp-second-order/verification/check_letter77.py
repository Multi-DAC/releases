"""Cold checks of Wren's letter 77 with MY operator conventions (edge_modes_dw / check_prop13_kernel), now with U(1)
gauge links. Sealed first read: letter77-firstread.md (13:03).

(A) Catch 2: gamma5 = -i gamma0 gamma1, and C = sigma2 gives C gamma_mu C^-1 = -gamma_mu^T, C P_pm C^-1 = P_mp^T.
(B) Bug A: the k1 = 0 DW surface pole at m = 0.7, S = 16, a5 = 0.2 (L5 = 80): 1 - z0 against 2m(1 - a5 m)^L5 and
    2m e^{-mS}.
(C) Lemma 1, the Gram proof: open time, flat crossing, random positive-half U(1) fields.
    (C1) D_-[theta U+] in mirror row order equals Gamma0 D_+[U+]^dag Gamma0 (the one step not written out by hand);
    (C2) K_ij = det D[theta x_i, 1, x_j] is Hermitian PSD on a random family; controls: crossing link != 1, and
    no conjugation of the mirrored temporal links.
(D) P5 Step 1-2 at a random U: T5 Hermitian PD; T5 = I - a5 H_W + a5^2 R; det D_N = det D_DW / det D_PV (5d, mine).
(E) P5, the Hessian: mixed second derivatives (distinct links) of log det D_N[U] at U = 1 against those of
    log det D_S[U], D_S = (1 + g5 tanh(S H_W / 2))/2 (the a5 -> 0 limit at fixed S), and of log det D_ov[U].
    AP time and AP space (no zero mode). Expect |H_N - H_S| = O(a5) uniformly in S, and |H_S - H_ov| -> 0 in S.
"""
import sys
import numpy as np

sys.path.insert(0, ".")
from edge_modes_dw import g0, g1, g5, Pp, Pm, I2  # noqa: E402
from check_letter75_h4 import blocks  # noqa: E402
from check_props2_4 import poles  # noqa: E402

rng = np.random.default_rng(77)


def wilson_hops(add, idx, nt, Ls, U1, U0, twist, time_bc, reps):
    """Add the 2d Wilson hops (r = 1) for every replica index rep in reps (the s slices)."""
    tw = np.exp(1j * twist)
    for t in range(nt):
        for x in range(Ls):
            xn, ph = (x + 1) % Ls, (tw if x == Ls - 1 else 1.0)
            u = U1[t, x] * ph
            for s in reps:
                add(idx(t, x, s), idx(t, xn, s), -0.5 * (I2 - g1) * u)
                add(idx(t, xn, s), idx(t, x, s), -0.5 * (I2 + g1) * np.conj(u))
            if t + 1 < nt or time_bc == "ap":
                tn = (t + 1) % nt
                u = U0[t, x] * (-1.0 if t + 1 == nt else 1.0)
                for s in reps:
                    add(idx(t, x, s), idx(tn, x, s), -0.5 * (I2 - g0) * u)
                    add(idx(tn, x, s), idx(t, x, s), -0.5 * (I2 + g0) * np.conj(u))


def build5(nt, Ls, L5, m, a5, mf, U1, U0, twist=np.pi, time_bc="open"):
    n = nt * Ls * L5 * 2
    idx = lambda t, x, s: ((t * Ls + x) * L5 + s) * 2  # noqa: E731
    D = np.zeros((n, n), complex)

    def add(i, j, B):
        D[i:i + 2, j:j + 2] += B

    for t in range(nt):
        for x in range(Ls):
            for s in range(L5):
                add(idx(t, x, s), idx(t, x, s), (2 - m + 1 / a5) * I2)
                if s + 1 < L5:
                    add(idx(t, x, s), idx(t, x, s + 1), -(1 / a5) * Pm)
                    add(idx(t, x, s + 1), idx(t, x, s), -(1 / a5) * Pp)
            if mf != 0:
                add(idx(t, x, L5 - 1), idx(t, x, 0), (mf / a5) * Pm)
                add(idx(t, x, 0), idx(t, x, L5 - 1), (mf / a5) * Pp)
    wilson_hops(add, idx, nt, Ls, U1, U0, twist, time_bc, range(L5))
    return D


def build2(nt, Ls, m, U1, U0, twist=np.pi, time_bc="ap"):
    """X[U] = D_W[U] - m on the 2d lattice; index (t, x, spin)."""
    n = nt * Ls * 2
    idx = lambda t, x, s: (t * Ls + x) * 2  # noqa: E731
    X = np.zeros((n, n), complex)

    def add(i, j, B):
        X[i:i + 2, j:j + 2] += B

    for t in range(nt):
        for x in range(Ls):
            add(idx(t, x, 0), idx(t, x, 0), (2 - m) * I2)
    wilson_hops(add, idx, nt, Ls, U1, U0, twist, time_bc, [0])
    return X


def G5(n2):
    return np.kron(np.eye(n2 // 2), g5)


def chiral(n2):
    """Permutation to the (P+, P-) block basis: spin 0 (gamma5 = +1) first."""
    return np.r_[np.arange(0, n2, 2), np.arange(1, n2, 2)]


def T5_of(X, a5):
    n2 = X.shape[0]
    p = chiral(n2)
    B = (np.eye(n2) + a5 * X)[np.ix_(p, p)]
    h = n2 // 2
    Bpp, Bpm, Bmp, Bmm = B[:h, :h], B[:h, h:], B[h:, :h], B[h:, h:]
    Bi = np.linalg.inv(Bpp)
    T = np.block([[Bi, -Bi @ Bpm], [Bmp @ Bi, Bmm - Bmp @ Bi @ Bpm]])
    Xc = X[np.ix_(p, p)]
    Xpp, Xpm, Xmp = Xc[:h, :h], Xc[:h, h:], Xc[h:, :h]
    R = np.block([[Xpp @ Xpp @ Bi, Xpp @ Bi @ Xpm], [-Xmp @ Xpp @ Bi, -Xmp @ Bi @ Xpm]])
    g5c = np.diag(np.r_[np.ones(h), -np.ones(h)])
    HW = g5c @ Xc
    return T, R, HW, g5c, p


def partA():
    print("(A) gamma algebra")
    print("  |gamma5 - (-i g0 g1)| = %.1e, |gamma5 - (+i g0 g1)| = %.1e"
          % (np.abs(g5 + 1j * g0 @ g1).max(), np.abs(g5 - 1j * g0 @ g1).max()))
    C = g1
    Ci = np.linalg.inv(C)
    for name, g in (("g0", g0), ("g1", g1), ("g5", g5)):
        print("  C %s C^-1 + %s^T: %.1e" % (name, name, np.abs(C @ g @ Ci + g.T).max()))
    print("  C P+ C^-1 - P-^T: %.1e" % np.abs(C @ Pp @ Ci - Pm.T).max())


def partB():
    print("(B) Bug A: k1 = 0 DW surface pole")
    m, a5, S = 0.7, 0.2, 16.0
    for a5 in (0.2, 0.1):
        L5 = int(round(S / a5))
        lam, _ = poles(blocks(L5, 0.0, m, a5, 0.0))
        z0 = lam[lam < 1].max()
        print("  a5=%.2f L5=%d: 1 - z0 = %.4e; 2m(1 - a5 m)^L5 = %.4e (ratio %.4f, 1 - m a5/2 = %.4f);"
              " 2m e^-mS = %.4e" % (a5, L5, 1 - z0, 2 * m * (1 - a5 * m) ** L5,
                                    (1 - z0) / (2 * m * (1 - a5 * m) ** L5), 1 - m * a5 / 2, 2 * m * np.exp(-m * S)))


def partC():
    print("(C) Lemma 1 Gram proof, open time, flat crossing")
    nt, Ls, L5, m, a5, twist = 4, 2, 2, 0.7, 0.5, 0.73
    c = nt // 2 - 1          # plane between rows c and c + 1; mirror row r -> 2c + 1 - r
    half = nt // 2

    def fields(pos1, pos0, cross, conj_temporal=True):
        """pos1[r', x]: spatial links on + rows r = c + 1 + r'; pos0[r', x]: temporal r -> r + 1, r = c + 1 + r'."""
        U1 = np.ones((nt, Ls), complex)
        U0 = np.ones((nt, Ls), complex)
        for rp in range(half):
            r = c + 1 + rp
            U1[r] = pos1[rp]
            U1[2 * c + 1 - r] = pos1[rp]
        for rp in range(half - 1):
            r = c + 1 + rp
            U0[r] = pos0[rp]
            U0[2 * c - r] = np.conj(pos0[rp]) if conj_temporal else pos0[rp]
        U0[c] = cross
        return U1, U0

    def rand_conf():
        return (np.exp(1j * rng.uniform(-np.pi, np.pi, (half, Ls))),
                np.exp(1j * rng.uniform(-np.pi, np.pi, (half - 1, Ls))))

    for mf in (0.0, 1.0):
        # (C1) the mirror identity
        p1, p0 = rand_conf()
        U1, U0 = fields(p1, p0, np.ones(Ls))
        D = build5(nt, Ls, L5, m, a5, mf, U1, U0, twist, "open")
        blk = Ls * L5 * 2
        Dp = D[(c + 1) * blk:, (c + 1) * blk:]
        Dm = D[:(c + 1) * blk, :(c + 1) * blk]
        perm = np.concatenate([np.arange(r * blk, (r + 1) * blk) for r in range(c, -1, -1)])
        Dm_mirror = Dm[np.ix_(perm, perm)]
        Gam0 = np.kron(np.eye(blk // 2 * half), g0)
        e1 = np.abs(Dm_mirror - Gam0 @ Dp.conj().T @ Gam0).max()
        # (C2) the kernel on a family
        fam = [rand_conf() for _ in range(8)]

        def K(cross=np.ones(Ls), conj_temporal=True):
            M = np.zeros((8, 8), complex)
            for i, (a1, a0) in enumerate(fam):
                for j, (b1, b0) in enumerate(fam):
                    U1 = np.ones((nt, Ls), complex)
                    U0 = np.ones((nt, Ls), complex)
                    Ua1, Ua0 = fields(a1, a0, cross, conj_temporal)
                    Ub1, Ub0 = fields(b1, b0, cross, conj_temporal)
                    U1[:c + 1], U0[:c] = Ua1[:c + 1], Ua0[:c]
                    U1[c + 1:], U0[c + 1:] = Ub1[c + 1:], Ub0[c + 1:]
                    U0[c] = cross
                    M[i, j] = np.linalg.det(build5(nt, Ls, L5, m, a5, mf, U1, U0, twist, "open"))
            return M

        def report(name, M):
            herm = np.abs(M - M.conj().T).max() / np.abs(M).max()
            ev = np.linalg.eigvalsh((M + M.conj().T) / 2)
            print("  mf=%.0f %-28s Hermitian dev %.1e; eigenvalues min/max = %+.3e" % (mf, name, herm, ev.min() / ev.max()))

        print("  mf=%.0f (C1) |D_-[theta U] (mirror order) - Gamma0 D_+^dag Gamma0| = %.1e" % (mf, e1))
        report("flat crossing:", K())
        report("control, crossing e^{0.9i}:", K(cross=np.exp(0.9j) * np.ones(Ls)))
        # Added at clock reading 13:13: a crossing phase constant in x is a gauge transformation of one half, so the
        # line above is NOT a control. An x-dependent crossing changes the plaquettes on the plane; that one is.
        report("control, crossing (0.9,-0.4):", K(cross=np.exp(1j * np.array([0.9, -0.4]))))
        report("control, crossing (2.0, 0.0):", K(cross=np.exp(1j * np.array([2.0, 0.0]))))
        report("control, no conj on U0:", K(conj_temporal=False))


def partD():
    print("(D) P5 Steps 1-2 at a random U (nt = 4 AP time, Ls = 2, twist pi)")
    nt, Ls, m = 4, 2, 0.7
    U1 = np.exp(1j * rng.uniform(-np.pi, np.pi, (nt, Ls)))
    U0 = np.exp(1j * rng.uniform(-np.pi, np.pi, (nt, Ls)))
    X = build2(nt, Ls, m, U1, U0, np.pi, "ap")
    for a5, L5 in ((0.5, 4), (1.2, 3), (0.1, 6)):
        T, R, HW, g5c, p = T5_of(X, a5)
        herm = np.abs(T - T.conj().T).max()
        evmin = np.linalg.eigvalsh((T + T.conj().T) / 2).min()
        mer = np.abs(T - (np.eye(len(T)) - a5 * HW + a5 ** 2 * R)).max()
        TL = np.linalg.matrix_power(T, L5)
        h = len(T) // 2
        detN = np.linalg.det(TL[h:, h:]) / np.linalg.det(np.eye(len(T)) + TL)
        DW = build5(nt, Ls, L5, m, a5, 0.0, U1, U0, np.pi, "ap")
        PV = build5(nt, Ls, L5, m, a5, 1.0, U1, U0, np.pi, "ap")
        sw, lw = np.linalg.slogdet(DW)
        sp, lp = np.linalg.slogdet(PV)
        ratio = sw / sp * np.exp(lw - lp)
        print("  a5=%.2f L5=%d (a5 m = %.2f): T5 Hermitian %.1e, min eig %.3e; T5 - (I - a5 H_W + a5^2 R) %.1e;"
              " |det D_N / (det DW / det PV) - 1| = %.1e" % (a5, L5, a5 * m, herm, evmin, mer, abs(detN / ratio - 1)))


def hessians(Ffun, nt, Ls, links, h=1e-3):
    """Mixed second derivatives (a != b) of log det F[U] at U = 1, by central differences of F itself."""
    def F(angles):
        U1 = np.ones((nt, Ls), complex)
        U0 = np.ones((nt, Ls), complex)
        for (mu, t, x), th in zip(links, angles):
            (U1 if mu == 1 else U0)[t, x] *= np.exp(1j * th)
        return Ffun(U1, U0)

    k = len(links)
    F0 = F(np.zeros(k))
    Fi = np.linalg.inv(F0)
    e = np.eye(k)
    dF = [(F(h * e[a]) - F(-h * e[a])) / (2 * h) for a in range(k)]
    H = np.zeros((k, k), complex)
    for a in range(k):
        for b in range(a + 1, k):
            d2 = (F(h * (e[a] + e[b])) - F(h * (e[a] - e[b])) - F(h * (e[b] - e[a])) + F(-h * (e[a] + e[b]))) / (4 * h * h)
            H[a, b] = np.trace(Fi @ d2) - np.trace(Fi @ dF[a] @ Fi @ dF[b])
            H[b, a] = H[a, b]
    return H


def partE():
    print("(E) P5 Hessian: log det D_N vs D_S (fixed S, a5 -> 0) vs D_ov; nt = 8 AP time, Ls = 4, twist pi, m = 0.7")
    nt, Ls, m = 8, 4, 0.7
    links = [(1, 3, 0), (1, 4, 1), (0, 3, 0), (0, 2, 2), (1, 6, 3)]

    def X_of(U1, U0):
        return build2(nt, Ls, m, U1, U0, np.pi, "ap")

    def Dfun(fn):
        def f(U1, U0):
            X = X_of(U1, U0)
            G = G5(len(X))
            HW = G @ X
            w, V = np.linalg.eigh((HW + HW.conj().T) / 2)
            return 0.5 * (np.eye(len(X)) + G @ (V * fn(w)) @ V.conj().T)
        return f

    def DN(a5, S):
        L5 = int(round(S / a5))

        def f(U1, U0):
            X = X_of(U1, U0)
            T, R, HW, g5c, p = T5_of(X, a5)
            w, V = np.linalg.eigh((T + T.conj().T) / 2)
            Ht = -np.log(w) / a5
            Dc = 0.5 * (np.eye(len(X)) + g5c @ (V * np.tanh(L5 * a5 * Ht / 2)) @ V.conj().T)
            q = np.argsort(p)
            return Dc[np.ix_(q, q)]
        return f

    Hov = hessians(Dfun(np.sign), nt, Ls, links)
    print("  |H_ov| max %.4e" % np.abs(Hov).max())
    for S in (4.0, 8.0, 16.0):
        HS = hessians(Dfun(lambda w: np.tanh(S * w / 2)), nt, Ls, links)
        line = []
        for a5 in (0.2, 0.1, 0.05, 0.025):
            HN = hessians(DN(a5, S), nt, Ls, links)
            line.append("a5=%.3f: %.2e" % (a5, np.abs(HN - HS).max()))
        print("  S=%4.0f: |H_S - H_ov| = %.2e; |H_N - H_S|: %s" % (S, np.abs(HS - Hov).max(), ", ".join(line)))
        sys.stdout.flush()


def partE2():
    """H_N against the 5d DW - PV Hessian, small lattice, by central differences of log det directly."""
    print("(E2) Hessian of log det D_N vs log det D_DW - log det D_PV (5d, mine): nt = 4 AP time, Ls = 2, L5 = 6, a5 = 0.5")
    nt, Ls, m, a5, L5 = 4, 2, 0.7, 0.5, 6
    links = [(1, 1, 0), (1, 2, 1), (0, 1, 0), (0, 2, 1)]

    def logdet(M):
        s, l = np.linalg.slogdet(M)
        return (s, l)

    def sub(A, B):
        return (A[0] / B[0], A[1] - B[1])

    def f5(angles):
        U1 = np.ones((nt, Ls), complex)
        U0 = np.ones((nt, Ls), complex)
        for (mu, t, x), th in zip(links, angles):
            (U1 if mu == 1 else U0)[t, x] *= np.exp(1j * th)
        return sub(logdet(build5(nt, Ls, L5, m, a5, 0.0, U1, U0, np.pi, "ap")),
                   logdet(build5(nt, Ls, L5, m, a5, 1.0, U1, U0, np.pi, "ap")))

    def fN(angles):
        U1 = np.ones((nt, Ls), complex)
        U0 = np.ones((nt, Ls), complex)
        for (mu, t, x), th in zip(links, angles):
            (U1 if mu == 1 else U0)[t, x] *= np.exp(1j * th)
        X = build2(nt, Ls, m, U1, U0, np.pi, "ap")
        T, R, HW, g5c, p = T5_of(X, a5)
        TL = np.linalg.matrix_power(T, L5)
        h = len(T) // 2
        return sub(logdet(TL[h:, h:]), logdet(np.eye(len(T)) + TL))

    def rel(f, ang, c0):
        s, l = f(ang)
        return np.log(s / c0[0] * np.exp(l - c0[1]))

    k, hh = len(links), 1e-3
    e = np.eye(k)
    worst = 0.0
    for a in range(k):
        for b in range(a + 1, k):
            vals = []
            for f in (f5, fN):
                c0 = f(np.zeros(k))
                v = (rel(f, hh * (e[a] + e[b]), c0) - rel(f, hh * (e[a] - e[b]), c0)
                     - rel(f, hh * (e[b] - e[a]), c0) + rel(f, -hh * (e[a] + e[b]), c0))
                vals.append(v / (4 * hh * hh))
            worst = max(worst, abs(vals[0] - vals[1]) / max(abs(vals[0]), 1e-12))
            print("  links %d,%d: 5d %+.6e%+.6ei, D_N %+.6e%+.6ei" % (a, b, vals[0].real, vals[0].imag, vals[1].real, vals[1].imag))
    print("  max relative difference %.1e" % worst)


if __name__ == "__main__":
    partA()
    partB()
    sys.stdout.flush()
    partC()
    sys.stdout.flush()
    partD()
    sys.stdout.flush()
    partE2()
    sys.stdout.flush()
    partE()
    print("DONE")
