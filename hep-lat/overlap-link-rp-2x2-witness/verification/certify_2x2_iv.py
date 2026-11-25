"""Certified 2 x 2 witness, INTERVAL SIDE (mpmath.iv). Design: CERTIFIED-PLAN.md "Update 21:30"; float side certify_2x2.py.

L = 4, Ls = 1, r = 1, m = 0.7, eps = 0.2 (all entered as intervals from decimal strings), distributed twist Phi/8.
Part 1: lambda_min(X^dag X) >= lam' at ng grid angles by interval Cholesky (real 2N embedding, every pivot > 0
certified); sigma_min is 1/(2L)-Lipschitz, so lambda_min >= (sqrt(lam') - (2 pi/ng)/(4L))^2 on the whole circle.
Part 2: strip budget from that floor (eta, delta, Delta, c0, M_D) in interval arithmetic; trapezoid bound
Q = 2 M_D^16 / (e^{n eta} - 1) per entry.
Part 3: at Phi_k = 2 pi k/n, each pair: M0 = exactly Hermitian approximation of H^{1/2} (mp.eigh at high precision);
r >= ||M0^2 - H||_F; a = certified lower bound on lambda_min(M0) (interval Cholesky of M0 - a I);
rho = r/(a sqrt(a^2-r)(a+sqrt(a^2-r))) bounds ||H^{-1/2} - M0^{-1}||; det D0 = det(M0 + X)/(2^16 det M0);
|det D - det D0| <= 16 ||E|| (||D0|| + ||E||)^15, ||E|| <= M_X rho/2, ||D0|| <= (1 + M_X/a)/2, M_X = 3.3.
Result: interval for v^T K v = mean_k [f11 + f22 - f12 - f21](Phi_k) + [-4Q, 4Q] (real and imaginary parts).
Usage: certify_2x2_iv.py [ng] [n] [eta] [prec]
"""
import sys
import time
from mpmath import iv, mp, mpf

L, N = 4, 16


def setprec(p):
    iv.prec = p; mp.prec = p


def cexp(th):
    return iv.mpc(iv.cos(th), iv.sin(th))


def cconj(z):
    return iv.mpc(z.real, -z.imag)


def H_of(A):
    """conjugate transpose (iv's own conjugate fails on complex intervals)"""
    n = A.rows
    B = iv.matrix(n, n)
    for i in range(n):
        for j in range(n):
            B[j, i] = cconj(A[i, j])
    return B


def upper(x):
    return mpf(x.b)


def X_iv(i, j, phi):
    """X = D_w - m for pair (A_i reflected, A_j), twist e^{i phi/8} on every time link (forward), e^{-i phi/8} backward."""
    eps = iv.mpf('0.2'); m = iv.mpf('0.7')
    f = cexp(phi / 8); b = cexp(-phi / 8)
    tA = -1 if i == 0 else -2
    tB = 2 if j == 0 else 3
    one = iv.mpc(1, 0); e = cexp(eps)
    X = iv.matrix(N, N)
    for k in range(N):
        X[k, k] = iv.mpc(2 - m, 0)
    for t in range(-3, 5):
        s = (t + 3) % 8
        # time hop, U0 = 1, antiperiodic sign at t = L
        jj = (t + 4) % 8
        sg = -1 if t == 4 else 1
        uf = f * sg; ub = b * sg
        # -(1/2)(I - g0) = -(1/2)[[1,-1],[-1,1]] forward; -(1/2)(I + g0) = -(1/2)[[1,1],[1,1]] backward
        for (p, q, c) in ((0, 0, 1), (0, 1, -1), (1, 0, -1), (1, 1, 1)):
            X[2 * s + p, 2 * jj + q] += -0.5 * c * uf
        for (p, q, c) in ((0, 0, 1), (0, 1, 1), (1, 0, 1), (1, 1, 1)):
            X[2 * jj + p, 2 * s + q] += -0.5 * c * ub
        # space hop to itself (Ls = 1): -(1/2)(I - g1) u - (1/2)(I + g1) conj(u); g1 = [[0,-i],[i,0]]
        u = e if t in (tA, tB) else one
        uc = cconj(u)
        Im = iv.mpc(0, 1)
        blk = [[-0.5 * u - 0.5 * uc, -0.5 * Im * u + 0.5 * Im * uc],
               [0.5 * Im * u - 0.5 * Im * uc, -0.5 * u - 0.5 * uc]]
        for p in range(2):
            for q in range(2):
                X[2 * s + p, 2 * s + q] += blk[p][q]
    return X


def embed(A):
    """real symmetric 2n embedding [[Re, -Im], [Im, Re]] of a Hermitian interval matrix"""
    n = A.rows
    E = iv.matrix(2 * n, 2 * n)
    for i in range(n):
        for j in range(n):
            re, im = A[i, j].real, A[i, j].imag
            E[i, j] = re; E[i + n, j + n] = re
            E[i, j + n] = -im; E[i + n, j] = im
    return E


def chol_pd(A):
    """interval Cholesky; True only if every pivot is certified > 0 (then every symmetric matrix in A is PD)"""
    n = A.rows
    Lm = iv.matrix(n, n)
    for j in range(n):
        s = A[j, j]
        for k in range(j):
            s = s - Lm[j, k] * Lm[j, k]
        if not (s > 0) is True:
            return False
        d = iv.sqrt(s)
        Lm[j, j] = d
        for i in range(j + 1, n):
            t = A[i, j]
            for k in range(j):
                t = t - Lm[i, k] * Lm[j, k]
            Lm[i, j] = t / d
    return True


def frob_upper(A):
    s = iv.mpf(0)
    for i in range(A.rows):
        for j in range(A.cols):
            s += A[i, j].real ** 2 + A[i, j].imag ** 2
    return upper(iv.sqrt(s))


def mid_matrix(A):
    M = mp.matrix(A.rows, A.cols)
    for i in range(A.rows):
        for j in range(A.cols):
            M[i, j] = mp.mpc(A[i, j].real.mid, A[i, j].imag.mid)
    return M


def det_enclosure(Xi):
    """returns (det D0 interval, err radius upper bound, a, r)"""
    H = H_of(Xi) * Xi
    Hm = mid_matrix(H)
    Ev, Q = mp.eigh(Hm)
    sq = [mp.sqrt(x) for x in Ev]
    M0 = Q * mp.diag(sq) * Q.H
    M0i = iv.matrix(N, N)
    for k in range(N):
        M0i[k, k] = iv.mpc(iv.mpf(mp.re(M0[k, k])), 0)
        for l in range(k + 1, N):
            z = (M0[k, l] + mp.conj(M0[l, k])) / 2
            M0i[k, l] = iv.mpc(iv.mpf(mp.re(z)), iv.mpf(mp.im(z)))
            M0i[l, k] = iv.mpc(iv.mpf(mp.re(z)), iv.mpf(-mp.im(z)))
    r = frob_upper(M0i * M0i - H)
    a = mp.mpf(min(sq)) * (1 - mp.mpf(10) ** -6)
    A = M0i - iv.mpc(iv.mpf(a), 0) * iv.eye(N)
    if not chol_pd(embed(A)):
        raise RuntimeError("lambda_min(M0) >= a not certified")
    ai = iv.mpf(a); ri = iv.mpf(r)
    # Wren letter 70 gap 3: the residual bound needs a > 0 (a is a lower bound on lambda_min(M0))
    if not (ai > 0) is True:
        raise RuntimeError("a > 0 not certified")
    if not (ai * ai > ri) is True:
        raise RuntimeError("a^2 > r not certified")
    sr = iv.sqrt(ai * ai - ri)
    rho = ri / (ai * sr * (ai + sr))
    MX = iv.mpf('3.3')
    Enorm = MX * rho / 2
    D0n = (1 + MX / ai) / 2
    err = 16 * Enorm * (D0n + Enorm) ** 15
    num = iv.det(M0i + Xi); den = iv.det(M0i)
    # mpmath's det returns 0 when its LU step calls the matrix singular: refuse anything that is not certified nonzero
    if not isinstance(den, type(iv.mpc(0, 0))) or not (den.real > 0) is True:
        raise RuntimeError("det M0 not certified positive")
    if not isinstance(num, type(iv.mpc(0, 0))) or not ((num.real > 0) is True or (num.real < 0) is True
                                                       or (num.imag > 0) is True or (num.imag < 0) is True):
        raise RuntimeError("det(M0 + X) not certified nonzero")
    dD0 = num / (den * iv.mpf(2) ** 16)
    return dD0, upper(err), a, r


def main():
    ng = int(sys.argv[1]) if len(sys.argv) > 1 else 512
    n = int(sys.argv[2]) if len(sys.argv) > 2 else 0
    eta_s = sys.argv[3] if len(sys.argv) > 3 else "0.47"
    prec = int(sys.argv[4]) if len(sys.argv) > 4 else 200
    setprec(prec)
    pairs = [(0, 0), (0, 1), (1, 0), (1, 1)]
    t0 = time.time()
    print("certify_2x2_iv ng=%d n=%d eta=%s prec=%d" % (ng, n, eta_s, prec), flush=True)

    # Part 1
    sp = iv.mpf('0.6848'); lamp = sp * sp
    for g in range(ng):
        phi = 2 * iv.pi * g / ng
        for (i, j) in pairs:
            Xi = X_iv(i, j, phi)
            H = H_of(Xi) * Xi
            if not chol_pd(embed(H - iv.mpc(lamp, 0) * iv.eye(N))):
                print("PART1 FAIL at g=%d pair %s" % (g, (i, j)), flush=True)
                return
        if g % 64 == 0:
            print("  part1 g=%d ok (%.0f s)" % (g, time.time() - t0), flush=True)
    sfloor = sp - (2 * iv.pi / ng) / (4 * L)
    # Wren letter 70 gap 3: squaring a negative sigma floor would give a false positive lambda floor
    if not (sfloor > 0) is True:
        print("PART1 FAIL: sigma floor not certified positive", flush=True); return
    lam = sfloor * sfloor
    print("PART1 certified: lambda_min >= %s on the whole circle (sigma floor %s)" % (iv.nstr(lam, 10), iv.nstr(sfloor, 10)), flush=True)

    # Part 2
    MX = iv.mpf('3.3')
    eta = iv.mpf(eta_s)
    delta = iv.exp(eta / (2 * L)) - 1
    Delta = 2 * MX * delta + delta * delta
    c0 = lam - Delta
    if not (c0 > 0) is True:
        print("PART2 FAIL: c0 not certified positive", flush=True); return
    MD = (1 + (MX + delta) / iv.sqrt(c0)) / 2
    Mbound = MD ** 16
    if n == 0:
        n = int(mp.ceil((mp.log(2 * upper(Mbound)) - mp.log(mp.mpf('1e-11'))) / mp.mpf(eta_s))) + 1
    Q = 2 * Mbound / (iv.exp(n * eta) - 1)
    print("PART2 eta=%s delta<=%s Delta<=%s c0>=%s M_D<=%s M_D^16<=%s n=%d Q<=%s"
          % (eta_s, iv.nstr(delta, 8), iv.nstr(Delta, 8), iv.nstr(c0, 8), iv.nstr(MD, 8), iv.nstr(Mbound, 6), n,
             mp.nstr(upper(Q), 6)), flush=True)

    # Part 3
    # Wren letter 70 gap 3: the error sum is accumulated in iv.mpf (outward), not round-to-nearest mp.mpf
    S = iv.mpc(0, 0); errsum = iv.mpf(0); worst_err = mp.mpf(0)
    sgn = {(0, 0): 1, (1, 1): 1, (0, 1): -1, (1, 0): -1}
    Kacc = {p: iv.mpc(0, 0) for p in pairs}
    for k in range(n):
        phi = 2 * iv.pi * k / n
        for p in pairs:
            dD0, err, a, r = det_enclosure(X_iv(p[0], p[1], phi))
            S += sgn[p] * dD0; Kacc[p] += dD0
            errsum += iv.mpf(err); worst_err = max(worst_err, err)
        if k % 10 == 0:
            print("  part3 k=%d/%d a=%s r=%s err=%s (%.0f s)" % (k, n, mp.nstr(a, 6), mp.nstr(r, 4), mp.nstr(err, 4), time.time() - t0), flush=True)
    rad_iv = errsum / n + 4 * Q          # interval; its upper end is the certified radius
    rad = upper(rad_iv)
    vKv = S / n
    # endpoints in interval arithmetic: [lo - rad, hi + rad] computed outward
    re_lo = mpf((vKv.real - iv.mpf(rad)).a); re_hi = upper(vKv.real + iv.mpf(rad))
    im_lo = mpf((vKv.imag - iv.mpf(rad)).a); im_hi = upper(vKv.imag + iv.mpf(rad))
    for p in pairs:
        Kp = Kacc[p] / n
        print("  K%d%d center %s width(re) %s" % (p[0] + 1, p[1] + 1, mp.nstr(Kp.real.mid, 14), mp.nstr(Kp.real.delta, 3)), flush=True)
    print("RESULT v^T K v in Re [%s, %s], Im [%s, %s]; total radius %s (rounding %s, quadrature %s); worst per-point err %s"
          % (mp.nstr(re_lo, 10), mp.nstr(re_hi, 10), mp.nstr(im_lo, 4), mp.nstr(im_hi, 4),
             mp.nstr(rad, 4), mp.nstr(upper(errsum / n), 4), mp.nstr(upper(4 * Q), 4), mp.nstr(worst_err, 4)), flush=True)
    print("CERTIFIED NEGATIVE" if re_hi < 0 else "NOT CERTIFIED", flush=True)
    print("DONE %.0f s" % (time.time() - t0), flush=True)


if __name__ == "__main__":
    main()
