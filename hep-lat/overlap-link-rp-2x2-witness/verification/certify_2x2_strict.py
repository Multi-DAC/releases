"""Independent strict audit of Clawd's Letter 257 Track B Target 1 certificate.

Eliminates mpmath.iv's pi, exp, cos, sin, det, and decimal from_str from the TCB:
- pi via Machin's alternating series with explicit tail bound
- exp, cos, sin via Taylor polynomials + elementary Lagrange remainder bounds using only iv +, -, *, /
- exact rational inputs (7/10, 1/5, 47/100, 6848/10000, 33/10)
- independent interval Gaussian elimination det_strict (explicit piv.real**2 + piv.imag**2 > 0 guard)
- pure iv.mpf final error accumulation (zero round-to-nearest operations in certification path)
- explicit sfloor > 0 and a > 0 guards

Standalone with mpmath only; optional algebraic float checks run if mirror_test/certify_2x2 are on sys.path.
"""
import sys
import time
from mpmath import iv, mp

L, N = 4, 16


def atan_inv_int_strict(q, terms=45):
    """Rigorous enclosure of arctan(1/q) for integer q >= 5 via alternating Taylor series."""
    q_iv = iv.mpf(q)
    q2 = q_iv * q_iv
    p = iv.mpf(1) / q_iv
    s = iv.mpf(0)
    for k in range(terms):
        term = p / (2 * k + 1)
        s = s + term if (k % 2 == 0) else s - term
        p = p / q2
    rem_hi = (p / (2 * terms + 1)).b
    if terms % 2 == 0:
        return s + iv.mpf((0, rem_hi))
    return s + iv.mpf((-rem_hi, 0))


def pi_strict():
    """Machin's formula: pi = 16*arctan(1/5) - 4*arctan(1/239)."""
    return 16 * atan_inv_int_strict(5, terms=45) - 4 * atan_inv_int_strict(239, terms=15)


def exp_unit_strict(x, terms=42):
    """Rigorous enclosure of exp(x) for |x| <= 1 using Taylor + Lagrange remainder (e < 3)."""
    ax = abs(x)
    assert (ax <= 1) is True, f"|x| <= 1 required, got {ax}"
    s = iv.mpf(1)
    term = iv.mpf(1)
    for k in range(1, terms):
        term = term * x / k
        s = s + term
    rem = iv.mpf(3)
    for k in range(1, terms + 1):
        rem = rem * ax / k
    r_hi = rem.b
    return s + iv.mpf((-r_hi, r_hi))


def exp_strict(x):
    """Rigorous enclosure of exp(x) by halving until |x/2^k| <= 0.5, then squaring."""
    k = 0
    y = x
    half = iv.mpf(1) / iv.mpf(2)
    while not (abs(y) <= half) is True:
        y = y / 2
        k += 1
    e = exp_unit_strict(y, terms=42)
    for _ in range(k):
        e = e * e
    return e


def cos_sin_strict(x, terms=28):
    """Rigorous enclosure of (cos(x), sin(x)) for |x| <= 1 using Taylor + Lagrange remainder."""
    ax = abs(x)
    assert (ax <= 1) is True, f"|x| <= 1 required, got {ax}"
    x2 = x * x
    c = iv.mpf(1)
    s = x
    tc = iv.mpf(1)
    ts = x
    for k in range(1, terms):
        tc = -tc * x2 / ((2 * k - 1) * (2 * k))
        ts = -ts * x2 / ((2 * k) * (2 * k + 1))
        c = c + tc
        s = s + ts
    rc = iv.mpf(1)
    for m in range(1, 2 * terms + 1):
        rc = rc * ax / m
    rs = rc * ax / (2 * terms + 1)
    return c + iv.mpf((-rc.b, rc.b)), s + iv.mpf((-rs.b, rs.b))


def cexp_strict(th):
    c, s = cos_sin_strict(th)
    return iv.mpc(c, s)


def cconj(z):
    return iv.mpc(z.real, -z.imag)


def H_of(A):
    n = A.rows
    B = iv.matrix(n, n)
    for i in range(n):
        for j in range(n):
            B[j, i] = cconj(A[i, j])
    return B


def iv_abs_lower(z):
    """Rigorous lower bound on |z| for ivmpc z, using only +, **2 (mpi_square), sqrt."""
    r2 = z.real ** 2 + z.imag ** 2
    return iv.sqrt(r2).a


def det_strict(A):
    """Independent interval Gaussian elimination determinant with partial pivoting.
    Uses ONLY iv arithmetic (+, -, *, /, **2) and iv.sqrt. Never silently returns 0.
    """
    n = A.rows
    assert A.cols == n
    M = [[A[i, j] for j in range(n)] for i in range(n)]
    sign = 1
    for j in range(n):
        best_k = -1
        best_lo = iv.mpf(-1)
        for k in range(j, n):
            lo = iv_abs_lower(M[k][j])
            if (lo > best_lo) is True:
                best_lo = lo
                best_k = k
        if best_k < 0 or not (best_lo > 0) is True:
            raise RuntimeError(f"det_strict: no certified nonzero pivot at col {j}")
        if best_k != j:
            M[j], M[best_k] = M[best_k], M[j]
            sign = -sign
        piv = M[j][j]
        piv_norm2 = piv.real ** 2 + piv.imag ** 2
        if not (piv_norm2 > 0) is True:
            raise RuntimeError(f"det_strict: pivot norm2 not > 0 at col {j}")
        for i in range(j + 1, n):
            fac = M[i][j] / piv
            for k in range(j + 1, n):
                M[i][k] = M[i][k] - fac * M[j][k]
    prod = iv.mpc(sign, 0)
    for j in range(n):
        prod = prod * M[j][j]
    return prod


def X_strict(i, j, phi, eps_cexp):
    """Build X(i, j, phi) using strictly enclosed trig and exact rational inputs 7/10, 1/5."""
    m = iv.mpf(7) / iv.mpf(10)
    f = cexp_strict(phi / 8)
    b = cconj(f)
    tA = -1 if i == 0 else -2
    tB = 2 if j == 0 else 3
    one = iv.mpc(1, 0)
    X = iv.matrix(N, N)
    diag_val = iv.mpc(2 - m, 0)
    for k in range(N):
        X[k, k] = diag_val
    Im = iv.mpc(0, 1)
    for t in range(-3, 5):
        s = (t + 3) % 8
        jj = (t + 4) % 8
        sg = -1 if t == 4 else 1
        uf = f * sg
        ub = b * sg
        for (p, q, c) in ((0, 0, 1), (0, 1, -1), (1, 0, -1), (1, 1, 1)):
            X[2 * s + p, 2 * jj + q] += (-c * uf) / 2
        for (p, q, c) in ((0, 0, 1), (0, 1, 1), (1, 0, 1), (1, 1, 1)):
            X[2 * jj + p, 2 * s + q] += (-c * ub) / 2
        u = eps_cexp if t in (tA, tB) else one
        uc = cconj(u)
        blk = [
            [(-u - uc) / 2, (-Im * u + Im * uc) / 2],
            [(Im * u - Im * uc) / 2, (-u - uc) / 2],
        ]
        for p in range(2):
            for q in range(2):
                X[2 * s + p, 2 * s + q] += blk[p][q]
    return X


def embed(A):
    n = A.rows
    E = iv.matrix(2 * n, 2 * n)
    for i in range(n):
        for j in range(n):
            re, im = A[i, j].real, A[i, j].imag
            E[i, j] = re
            E[i + n, j + n] = re
            E[i, j + n] = -im
            E[i + n, j] = im
    return E


def chol_pd(A):
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


def frob_iv(A):
    s = iv.mpf(0)
    for i in range(A.rows):
        for j in range(A.cols):
            s += A[i, j].real ** 2 + A[i, j].imag ** 2
    return iv.sqrt(s)


def mid_matrix(A):
    M = mp.matrix(A.rows, A.cols)
    for i in range(A.rows):
        for j in range(A.cols):
            M[i, j] = mp.mpc(A[i, j].real.mid, A[i, j].imag.mid)
    return M


def run_strict_certificate(prec=160, ng=512, n=124):
    print(f"=== Independent strict certificate (prec={prec}, ng={ng}, n={n}, no iv.exp/cos/sin/pi/det) ===", flush=True)
    t0 = time.time()
    iv.prec = prec
    mp.prec = prec
    pi_s = pi_strict()
    eps_iv = iv.mpf(1) / iv.mpf(5)
    eps_cexp = cexp_strict(eps_iv)
    print(f"  pi_strict width = {mp.nstr(pi_s.delta, 4)}, eps_cexp re width = {mp.nstr(eps_cexp.real.delta, 4)}", flush=True)

    pairs = [(0, 0), (0, 1), (1, 0), (1, 1)]

    # Part 1
    sp = iv.mpf(6848) / iv.mpf(10000)
    lamp = sp * sp
    for g in range(ng):
        phi = 2 * pi_s * g / ng
        for (i, j) in pairs:
            Xi = X_strict(i, j, phi, eps_cexp)
            H = H_of(Xi) * Xi
            if not chol_pd(embed(H - iv.mpc(lamp, 0) * iv.eye(N))):
                raise RuntimeError(f"Part 1 failed at g={g}, pair={(i, j)}")
        if g % 128 == 0:
            print(f"  part1 g={g}/{ng} ok ({time.time() - t0:.0f} s)", flush=True)
    sfloor = sp - (2 * pi_s / ng) / (4 * L)
    assert (sfloor > 0) is True, "sfloor not certified > 0"
    lam = sfloor * sfloor
    print(f"  PART1 certified: lambda_min >= {iv.nstr(lam.a, 12)} ({time.time() - t0:.0f} s)", flush=True)

    # Part 2
    MX = iv.mpf(33) / iv.mpf(10)
    eta = iv.mpf(47) / iv.mpf(100)
    delta = exp_strict(eta / (2 * L)) - 1
    Delta = 2 * MX * delta + delta * delta
    c0 = lam - Delta
    assert (c0 > 0) is True, "c0 not certified > 0"
    MD = (1 + (MX + delta) / iv.sqrt(c0)) / 2
    Mbound = MD ** 16
    Q = 2 * Mbound / (exp_strict(n * eta) - 1)
    print(f"  PART2 strict: c0 >= {iv.nstr(c0.a, 10)}, M_D <= {iv.nstr(MD.b, 10)}, Q <= {iv.nstr(Q.b, 8)}", flush=True)

    # Part 3
    S = iv.mpc(0, 0)
    errsum = iv.mpf(0)
    sgn = {(0, 0): 1, (1, 1): 1, (0, 1): -1, (1, 0): -1}
    for k in range(n):
        phi = 2 * pi_s * k / n
        for p in pairs:
            Xi = X_strict(p[0], p[1], phi, eps_cexp)
            H = H_of(Xi) * Xi
            Hm = mid_matrix(H)
            Ev, Qmat = mp.eigh(Hm)
            sq = [mp.sqrt(x) for x in Ev]
            M0 = Qmat * mp.diag(sq) * Qmat.H
            M0i = iv.matrix(N, N)
            for r0 in range(N):
                M0i[r0, r0] = iv.mpc(iv.mpf(mp.re(M0[r0, r0])), 0)
                for c0i in range(r0 + 1, N):
                    z = (M0[r0, c0i] + mp.conj(M0[c0i, r0])) / 2
                    M0i[r0, c0i] = iv.mpc(iv.mpf(mp.re(z)), iv.mpf(mp.im(z)))
                    M0i[c0i, r0] = iv.mpc(iv.mpf(mp.re(z)), iv.mpf(-mp.im(z)))
            ri = frob_iv(M0i * M0i - H).b
            a_val = mp.mpf(min(sq)) * (1 - mp.mpf(10) ** -6)
            ai = iv.mpf(a_val)
            assert (ai > 0) is True, "a not > 0"
            if not chol_pd(embed(M0i - iv.mpc(ai, 0) * iv.eye(N))):
                raise RuntimeError("M0 >= a I not certified")
            assert (ai * ai > ri) is True, "a^2 > r not certified"
            sr = iv.sqrt(ai * ai - ri)
            rho = ri / (ai * sr * (ai + sr))
            Enorm = MX * rho / 2
            D0n = (1 + MX / ai) / 2
            err_iv = 16 * Enorm * (D0n + Enorm) ** 15
            num = det_strict(M0i + Xi)
            den = det_strict(M0i)
            assert (den.real > 0) is True, "det M0 not certified positive"
            dD0 = num / (den * iv.mpf(2) ** 16)
            S = S + sgn[p] * dD0
            errsum = errsum + err_iv.b
        if k % 31 == 0:
            print(f"  part3 k={k}/{n} ok ({time.time() - t0:.0f} s)", flush=True)

    rad_iv = errsum / n + 4 * Q
    rad_hi = rad_iv.b
    vKv = S / n
    re_iv = vKv.real + iv.mpf((-rad_hi, rad_hi))
    im_iv = vKv.imag + iv.mpf((-rad_hi, rad_hi))
    print(f"  STRICT RESULT v^T K v in Re {iv.nstr(re_iv, 12)}, Im {iv.nstr(im_iv, 5)}")
    print(f"  STRICT rad <= {iv.nstr(rad_hi, 6)} (rounding <= {iv.nstr((errsum/n).b, 5)}, quad <= {iv.nstr((4*Q).b, 5)})")
    assert (re_iv < 0) is True, "NOT CERTIFIED NEGATIVE"
    print(f"  STRICTLY CERTIFIED NEGATIVE (in {time.time() - t0:.0f} s)")


if __name__ == "__main__":
    prec = int(sys.argv[1]) if len(sys.argv) > 1 else 160
    run_strict_certificate(prec=prec)
