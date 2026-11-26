"""Cold audit of Wren's PROOF-T1-2-L4-and-overlap.md, Proposition 1.3, with my own position-space operator.

Build the 2d x L5 domain-wall (mf = 0) or Pauli-Villars (mf = 1) operator on an open time box of nt rows, Ls spatial
sites (twist pi), in a U(1) background equal to 1; my conventions (edge_modes_dw.op, extended to position space).
Reflection plane between rows c = nt/2 - 1 and c + 1; row c + t' is "t'" (t' = 1 is the first + row).
Links on the + side: spatial at t' = 1..T, temporal from t' to t' + 1, t' = 1..T. Mirrors: spatial at 1 - t'; temporal
from -t' to 1 - t' (orientation reversed, so the angle enters with a minus sign: theta a0(t') = -a0(-t')).
Mixed Hessian M[e, f] = sigma_e d^2 log det D / d a_(theta e) d a_f = -sigma_e Tr(G dD_(theta e) G dD_f).

Tests: (a) M Hermitian and PSD (or NSD: report both signs); (b) Fourier in x at p1: the 2 x 2 blocks C(ta, tb) depend
only on ta + tb (Hankel); (c) the sequence C(s), s = ta + tb - 2, obeys the recurrence whose characteristic roots are
lambda_i(k1) lambda_j(k1 + p1) (lambda < 1 the eigenvalues of my transfer matrix), with a control that drops one root;
(d) weights by least squares on those roots: each 2 x 2 weight PSD."""
import sys
import numpy as np

sys.path.insert(0, ".")
from edge_modes_dw import g0, g1, Pp, Pm, I2  # noqa: E402
from check_letter75_h4 import transfer  # noqa: E402


def build(nt, Ls, L5, m, a5, mf, twist=np.pi):
    n = nt * Ls * L5 * 2
    idx = lambda t, x, s: ((t * Ls + x) * L5 + s) * 2  # noqa: E731
    D = np.zeros((n, n), complex)
    links = {}   # (mu, t, x) -> (rows, cols, blocks) of dD/da at a = 0

    def add(i, j, B):
        D[i:i + 2, j:j + 2] += B

    tw = np.exp(1j * twist)
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
            # spatial hop x -> x + 1 (twist on the wrap)
            xn, ph = (x + 1) % Ls, (tw if x == Ls - 1 else 1.0)
            fw, bw = -0.5 * (I2 - g1) * ph, -0.5 * (I2 + g1) * np.conj(ph)
            ent = []
            for s in range(L5):
                add(idx(t, x, s), idx(t, xn, s), fw)
                add(idx(t, xn, s), idx(t, x, s), bw)
                ent += [(idx(t, x, s), idx(t, xn, s), 1j * fw), (idx(t, xn, s), idx(t, x, s), -1j * bw)]
            links[(1, t, x)] = ent
            if t + 1 < nt:
                fw, bw = -0.5 * (I2 - g0), -0.5 * (I2 + g0)
                ent = []
                for s in range(L5):
                    add(idx(t, x, s), idx(t + 1, x, s), fw)
                    add(idx(t + 1, x, s), idx(t, x, s), bw)
                    ent += [(idx(t, x, s), idx(t + 1, x, s), 1j * fw), (idx(t + 1, x, s), idx(t, x, s), -1j * bw)]
                links[(0, t, x)] = ent
    return D, links


def GdD(G, ent):
    """G @ dD for a link: dD has 2x2 blocks B at (i, j)."""
    Y = np.zeros_like(G)
    for i, j, B in ent:
        Y[:, j:j + 2] += G[:, i:i + 2] @ B
    return Y


def kernel(nt, Ls, L5, m, a5, mf, T, twist=np.pi):
    D, links = build(nt, Ls, L5, m, a5, mf, twist)
    G = np.linalg.inv(D)
    c = nt // 2 - 1
    plus, minus, sig = [], [], []
    for mu in (1, 0):
        for tp in range(1, T + 1):
            for x in range(Ls):
                if mu == 1:
                    plus.append((1, c + tp, x)); minus.append((1, c + 1 - tp, x)); sig.append(1.0)
                else:
                    plus.append((0, c + tp, x)); minus.append((0, c - tp, x)); sig.append(-1.0)
    Yp = [GdD(G, links[e]) for e in plus]
    Ym = [GdD(G, links[e]) for e in minus]
    N = len(plus)
    M = np.zeros((N, N), complex)
    for a in range(N):
        for b in range(N):
            M[a, b] = -sig[a] * np.sum(Ym[a] * Yp[b].T)
    return M, plus


def blocks_C(M, plus, Ls, T, p1):
    """C[ta, tb] 2x2 (components: spatial, temporal) at momentum p1, ta, tb = 1..T."""
    C = np.zeros((T, T, 2, 2), complex)
    for a, (mu, ta, xa) in enumerate(plus):
        for b, (nu, tb, xb) in enumerate(plus):
            C[ta - min(e[1] for e in plus), tb - min(e[1] for e in plus), 1 - mu, 1 - nu] += \
                np.exp(-1j * p1 * xa) * np.exp(1j * p1 * xb) * M[a, b] / Ls
    return C


def roots(Ls, L5, m, a5, mf, p1, twist=np.pi):
    mus = []
    for n in range(Ls):
        k1 = (2 * np.pi * n + twist) / Ls
        lam = [l for l in np.linalg.eigvalsh(transfer(L5, k1, m, a5, mf)[0]) if l < 1]
        lamp = [l for l in np.linalg.eigvalsh(transfer(L5, k1 + p1, m, a5, mf)[0]) if l < 1]
        mus += [x * y for x in lam for y in lamp]
    mus = np.array(sorted(mus))
    keep = [mus[0]]
    for u in mus[1:]:
        if abs(u - keep[-1]) > 1e-9:
            keep.append(u)
    return np.array(keep)


def run(Ls, L5, m, a5, mf, nt=64, T=8, twist=np.pi):
    name = "DW" if mf == 0 else "PV"
    M, plus = kernel(nt, Ls, L5, m, a5, mf, T, twist)
    herm = np.abs(M - M.conj().T).max() / np.abs(M).max()
    ev = np.linalg.eigvalsh((M + M.conj().T) / 2)
    print("%s twist=%.2f Ls=%d L5=%d m=%.2f a5=%.2f nt=%d T=%d: |M - M^H|/|M| %.1e; eig min %.3e max %.3e (ratio %.1e)"
          % (name, twist, Ls, L5, m, a5, nt, T, herm, ev.min(), ev.max(), ev.min() / ev.max()))
    for n in range(Ls):
        p1 = 2 * np.pi * n / Ls
        C = blocks_C(M, plus, Ls, T, p1)
        seq = {}
        hank = 0.0
        for ta in range(T):
            for tb in range(T):
                s = ta + tb
                if s in seq:
                    hank = max(hank, np.abs(C[ta, tb] - seq[s]).max())
                else:
                    seq[s] = C[ta, tb]
        S = np.array([seq[s] for s in range(2 * T - 1)])     # s = 0 .. 2T - 2
        scale = np.abs(S).max()
        mus = roots(Ls, L5, m, a5, mf, p1, twist)
        def rec_res(rts):
            p = np.poly(rts)          # highest power first
            d = len(rts)
            r = [np.abs(sum(p[d - k] * S[s0 + k] for k in range(d + 1))).max() for s0 in range(len(S) - d)]
            return max(r) / scale if r else np.nan
        res_all = rec_res(mus)
        res_drop = max(rec_res(np.delete(mus, i)) for i in range(len(mus))) if len(mus) > 1 else np.nan
        res_drop_min = min(rec_res(np.delete(mus, i)) for i in range(len(mus))) if len(mus) > 1 else np.nan
        # weights by least squares
        V = np.array([[u ** s for u in mus] for s in range(len(S))])
        W, *_ = np.linalg.lstsq(V, S.reshape(len(S), 4), rcond=None)
        fit = np.abs(V @ W - S.reshape(len(S), 4)).max() / scale
        wmin = [np.linalg.eigvalsh((Wn.reshape(2, 2) + Wn.reshape(2, 2).conj().T) / 2) for Wn in W]
        # moment Hankel PSD
        h = T - 1
        H0 = np.block([[S[i + j] for j in range(h)] for i in range(h)])
        H1 = np.block([[S[i + j + 1] for j in range(h)] for i in range(h)])
        e0 = np.linalg.eigvalsh((H0 + H0.conj().T) / 2)
        e1 = np.linalg.eigvalsh((H1 + H1.conj().T) / 2)
        print("   p1=%.3f: Hankel residual %.1e (scale %.2e); %d roots in (%.4f, %.4f); recurrence residual %.1e,"
              " dropping one root: %.1e .. %.1e; LS fit %.1e; weight eigs min %.2e (max %.2e); "
              "H0 min/max %.1e, H1 min/max %.1e"
              % (p1, hank, scale, len(mus), mus.min(), mus.max(), res_all, res_drop_min, res_drop, fit,
                 min(w.min() for w in wmin), max(w.max() for w in wmin), e0.min() / e0.max(), e1.min() / e1.max()))
    sys.stdout.flush()


if __name__ == "__main__":
    for mf in (0.0, 1.0):
        run(1, 2, 0.7, 1.0, mf)
        run(2, 2, 0.7, 1.0, mf)
        run(1, 3, 0.3, 0.5, mf)
