"""Certified 2 x 2 reflection-positivity witness, FLOAT SIDE (budget and checks; the interval side comes after).

Target (CERTIFIED-PLAN.md, "Update 20:31"): L = 4 (time sites -3..4, theta t = 1 - t), Ls = 1, U(1), overlap with Wilson
kernel r = 1 and kernel mass m = 0.7. Fields A1: a1(t=2) = eps, A2: a1(t=3) = eps, eps = 0.20, a0 = 0.
Entry K_ij = (1/2pi) Int dPhi det D[theta A_i, Phi on the crossing link t=0->1, A_j]; witness v = (1, -1):
v^T K v = K11 + K22 - K12 - K21 < 0.

Distributed twist: the crossing phase e^{i Phi} on one link is gauge-equivalent to e^{i Phi/(2L)} on each of the 2L
time links (same holonomy), so det D is unchanged. Continued to z = Phi + i y, forward time hops scale by
e^{-y/(2L)}, backward by e^{+y/(2L)}. Hops are -(S x P_-) u and -(S^dag x P_+) conj(u), P_- P_+ = 0, so
||X(z) - X(Phi)|| <= delta = e^{eta/(2L)} - 1 on |y| <= eta, and ||dX/dPhi|| <= 1/(2L).
On the strip, H(z) = Xt(z) X(z) with Xt(z) = X(conj z)^dag (holomorphic), ||H(z) - H(Phi)|| <= Delta = 2 M_X delta + delta^2.
If Delta < lambda_min, Re W(H(z)) >= c0 = lambda_min - Delta, ||H^{-1/2}|| <= c0^{-1/2}, ||D(z)|| <= M_D =
(1 + (M_X + delta) c0^{-1/2}) / 2, |det D(z)| <= M_D^dim. Periodic trapezoid with n points: error <= 2 M/(e^{n eta} - 1).
M_X: ||X(Phi)|| <= |2r - m| + 1 (time hops) + 1 (space hops, Ls = 1) = 3.3.
lambda_min: sigma_min(X(Phi)) is 1/(2L)-Lipschitz, so on a grid of spacing hg, sigma_min >= min_grid - hg/(4L).
This file checks all of that in floats and prints the budget. Nothing here is certified.
Usage: certify_2x2.py [ngrid] [target]
"""
import sys
import numpy as np
import mirror_test as mt
import mirror_kernel as mk
import mirror_orbit as mo

L, LS, M, EPS, R = 4, 1, 0.7, 0.20, 1.0
I2 = np.eye(2)


def fields():
    out = []
    for t in (2, 3):
        a0 = np.zeros((L - 1, LS)); a1 = np.zeros((L, LS))
        a1[t - 1, :] = EPS
        out.append((a0, a1))
    return out


def wilson_tw(lat, U0, U1, f, b, r=R):
    """mt.wilson with every forward time hop times f and every backward time hop times b (b = conj f on the real line)."""
    n = lat.n
    D = np.zeros((2 * n, 2 * n), complex)
    for t in lat.ts:
        for x in range(lat.Ls):
            i = lat.idx(t, x)
            D[2*i:2*i+2, 2*i:2*i+2] += 2 * r * I2
            j = lat.idx(t + 1, x); sgn = -1.0 if t == lat.L else 1.0
            u = U0[t + lat.L - 1, x] * sgn
            D[2*i:2*i+2, 2*j:2*j+2] += -0.5 * (r * I2 - mt.g0) * u * f
            D[2*j:2*j+2, 2*i:2*i+2] += -0.5 * (r * I2 + mt.g0) * np.conj(u) * b
            j = lat.idx(t, x + 1)
            u = U1[t + lat.L - 1, x]
            D[2*i:2*i+2, 2*j:2*j+2] += -0.5 * (r * I2 - mt.g1) * u
            D[2*j:2*j+2, 2*i:2*i+2] += -0.5 * (r * I2 + mt.g1) * np.conj(u)
    return D


def X_of(lat, U0, U1, z):
    """X(z) = D_w - m with the distributed twist e^{i z/(2L)} on every time link, holomorphic in z."""
    return wilson_tw(lat, U0, U1, np.exp(1j * z / (2 * L)), np.exp(-1j * z / (2 * L))) - M * np.eye(2 * lat.n)


def det_D(lat, U0, U1, z):
    """det D(z), D = (1 + X H^{-1/2})/2, H = X(conj z)^dag X(z), principal branch. Also min Re W(H) and ||D||."""
    X = X_of(lat, U0, U1, z)
    Xt = X_of(lat, U0, U1, np.conj(z)).conj().T
    H = Xt @ X
    w, V = np.linalg.eig(H)
    Hm = V @ np.diag(w ** -0.5) @ np.linalg.inv(V)
    D = 0.5 * (np.eye(len(H)) + X @ Hm)
    reW = np.linalg.eigvalsh(0.5 * (H + H.conj().T))[0]
    return np.linalg.det(D), reW, np.linalg.norm(D, 2), w


def pair(lat, F, i, j):
    return mk.pair_links(lat, F[i], F[j])


def main():
    ngrid = int(sys.argv[1]) if len(sys.argv) > 1 else 4096
    target = float(sys.argv[2]) if len(sys.argv) > 2 else 1e-11
    lat = mt.Lat(L, LS); F = fields(); dim = 2 * lat.n
    pairs = [(0, 0), (0, 1), (1, 0), (1, 1)]
    links = {p: pair(lat, F, *p) for p in pairs}

    print("== (a) distributed twist against mirror_orbit.dets (V = (e^{i Phi}, 1))")
    phis = np.array([0.0, 0.7, 2.0, 3.1, 4.4, 6.0])
    V = np.stack([np.exp(1j * phis), np.ones(len(phis))], axis=1)
    worst = 0.0
    for p in pairs:
        U0, U1 = links[p]
        d_ref, _ = mo.dets(lat, U0, U1, V, "ov", M)
        d_tw = np.array([det_D(lat, U0, U1, ph)[0] for ph in phis])
        rel = np.max(np.abs(d_tw - d_ref) / np.abs(d_ref))
        worst = max(worst, rel)
        print("  pair %s  max rel diff %.2e   |det| at Phi=0: %.6e" % (p, rel, abs(d_ref[0])))
    print("  worst %.2e" % worst)

    print("== (b) sigma_min, ||X||, Lipschitz on a grid of %d points" % ngrid)
    hg = 2 * np.pi / ngrid
    grid = hg * np.arange(ngrid)
    lam_floor = np.inf; MXmax = 0.0; lip = 0.0
    for p in pairs:
        U0, U1 = links[p]
        s = np.empty(ngrid); nx = np.empty(ngrid)
        for k, ph in enumerate(grid):
            sv = np.linalg.svd(X_of(lat, U0, U1, ph), compute_uv=False)
            s[k] = sv[-1]; nx[k] = sv[0]
        ds = np.abs(np.diff(np.append(s, s[0]))) / hg
        floor = s.min() - hg / (4 * L)
        print("  pair %s  sigma_min grid min %.6f at Phi %.4f | floor %.6f | lambda floor %.6f | max ||X|| %.4f | max |d sigma/dPhi| %.4f (bound %.4f)"
              % (p, s.min(), grid[s.argmin()], floor, floor ** 2, nx.max(), ds.max(), 1 / (2 * L)))
        lam_floor = min(lam_floor, floor ** 2); MXmax = max(MXmax, nx.max()); lip = max(lip, ds.max())
    MX = abs(2 * R - M) + 2.0
    print("  lambda_min floor %.6f ; ||X|| max seen %.4f vs bound M_X = %.2f ; Lipschitz seen %.4f vs %.4f"
          % (lam_floor, MXmax, MX, lip, 1 / (2 * L)))

    print("== (c) budget, target %.1e per entry" % target)
    best = None
    for frac in np.linspace(0.05, 0.95, 91):
        Delta = frac * lam_floor
        delta = -MX + np.sqrt(MX ** 2 + Delta)
        eta = 2 * L * np.log1p(delta)
        c0 = lam_floor - Delta
        MD = 0.5 * (1 + (MX + delta) / np.sqrt(c0))
        logM = dim * np.log(MD)
        n = int(np.ceil((np.log(2) + logM - np.log(target)) / eta)) + 1
        if best is None or n < best[0]:
            best = (n, frac, Delta, delta, eta, c0, MD, logM)
    n, frac, Delta, delta, eta, c0, MD, logM = best
    bound = 2 * np.exp(logM) / np.expm1(n * eta)
    print("  best: Delta = %.3f lambda_min = %.6f ; delta = %.6e ; eta = %.6f ; c0 = %.6f ; M_D = %.6f ; M_D^%d = %.3e ; n = %d ; bound %.2e"
          % (frac, Delta, delta, eta, c0, MD, dim, np.exp(logM), n, bound))
    frac2 = 0.5
    Delta2 = frac2 * lam_floor; delta2 = -MX + np.sqrt(MX ** 2 + Delta2); eta2 = 2 * L * np.log1p(delta2)
    c02 = lam_floor - Delta2; MD2 = 0.5 * (1 + (MX + delta2) / np.sqrt(c02))
    n2 = int(np.ceil((np.log(2) + dim * np.log(MD2) - np.log(target)) / eta2)) + 1
    print("  plan's choice Delta = lambda_min/2: eta = %.4f, M_D = %.4f, n = %d" % (eta2, MD2, n2))

    print("== (d) strip check at y = +-eta (holomorphic continuation)")
    for p in pairs:
        U0, U1 = links[p]
        worst_re = np.inf; worst_det = 0.0; worst_nd = 0.0; worst_w = np.inf
        for ph in np.linspace(0, 2 * np.pi, 64, endpoint=False):
            for y in (-eta, eta):
                d, reW, nd, w = det_D(lat, U0, U1, ph + 1j * y)
                worst_re = min(worst_re, reW); worst_det = max(worst_det, abs(d))
                worst_nd = max(worst_nd, nd); worst_w = min(worst_w, w.real.min())
        print("  pair %s  min Re W(H) %.5f (>= c0 %.5f: %s) | max ||D|| %.4f (<= M_D %.4f) | max |det| %.3e (<= %.3e) | min Re eig H %.5f"
              % (p, worst_re, c0, worst_re >= c0, worst_nd, MD, worst_det, np.exp(logM), worst_w))

    print("== (e) float trapezoid with the distributed twist, n = %d and 2n" % n)
    for nn in (n, 2 * n):
        ph = 2 * np.pi * np.arange(nn) / nn
        K = np.zeros((2, 2), complex)
        for p in pairs:
            U0, U1 = links[p]
            K[p] = np.mean([det_D(lat, U0, U1, x)[0] for x in ph])
        vKv = K[0, 0] + K[1, 1] - K[0, 1] - K[1, 0]
        print("  n=%d K11 %+.12e K22 %+.12e K12 %+.12e%+.1ei K21 %+.12e%+.1ei | v^T K v %+.6e%+.1ei"
              % (nn, K[0, 0].real, K[1, 1].real, K[0, 1].real, K[0, 1].imag, K[1, 0].real, K[1, 0].imag, vKv.real, vKv.imag))
    print("  reference (two_field_check, n = 16..64): K11 = K22 = +3.7845082952e-4, K12 = +3.7845495376e-4")


if __name__ == "__main__":
    main()
