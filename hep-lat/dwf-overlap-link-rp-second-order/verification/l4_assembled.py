"""L4 as assembled, my own computation (plan sealed in l4-assembled-plan.md at clock reading 13:44).

M_in(t) = inner (|mu| <= rho_*) moment of N~_DW - N~_PV, the full 2 x 2 (spatial, temporal) matrix, 1/Ls per site.
Time is L = infinity, through exact pole sums of D_k(z)^-1 = (B0 + z B+ + B-/z)^-1 from my blocks.
Validation: the full kernel against check_prop13_kernel.py (position space, open time) at small sizes.
"""
import sys
import numpy as np
import scipy.linalg as sla

sys.path.insert(0, ".")
from check_letter75_h4 import blocks  # noqa: E402
from edge_modes_dw import g0, g1, I2  # noqa: E402

fw, bw = -0.5 * (I2 - g1), -0.5 * (I2 + g1)       # spatial hop to x + 1, to x - 1
fw0, bw0 = -0.5 * (I2 - g0), -0.5 * (I2 + g0)     # temporal hop to t + 1, to t - 1


def a5max(m):
    return 2 * (1 - m) / (m * (2 - m))


def rstar(k, m):
    c = 2 - m - np.cos(k)
    return np.exp(-np.arccosh((c * c + 1 + np.sin(k) ** 2) / (2 * c)))


def rho_star(m, p1, n=20001):
    k = np.linspace(-np.pi, np.pi, n)
    return (rstar(k, m) * rstar(k + p1, m)).max()


def cluster_idx(vals, tol=1e-7):
    """Group numerically equal poles.

    Changed at clock reading 14:07, the launch, made just before it: the test is now on log|lam| (and sign), so an inside cluster and its partner 1/lam
    outside are grouped alike. The old mixed absolute/relative test split the outside partner of a near-degenerate
    inside pair: AssertionError (63, 64) at m = 0.9, twist 0, L5 = 64, a5 = 0.05.
    """
    sg, lg = np.sign(vals), np.log(np.abs(vals))
    order = np.lexsort((lg, sg))
    groups = []
    for i in order:
        j = groups[-1][-1] if groups else None
        # never merge across |lam| = 1 (added at clock reading 14:41, for prop4_large_t.py): at S = 64, m = 0.3 the
        # surface pole and its partner sit within 1e-10 of 1 and were one cluster. No L4 run met this case; each would
        # have stopped on the pairing assert.
        if j is not None and sg[i] == sg[j] and abs(lg[i] - lg[j]) <= tol and (lg[i] < 0) == (lg[j] < 0):
            groups[-1].append(i)
        else:
            groups.append([i])
    return groups


class Prop:
    """Exact pole data of D_k(z)^-1 for one fermion momentum k."""

    def __init__(self, L5, k, m, a5, mf):
        B0, Bp, Bm = blocks(L5, k, m, a5, mf)
        self.B = (B0, Bp, Bm)
        n = 2 * L5
        self.n = n
        G0 = np.kron(np.eye(L5), g0)
        Z, I = np.zeros((n, n)), np.eye(n)
        w, V = sla.eig(np.block([[Z, I], [-Bm, -B0]]), np.block([[I, Z], [Z, Bp]]))
        ok = np.isfinite(w) & (np.abs(w) > 1e-12)
        w, V = w[ok], V[:n, ok]
        self.imag = np.abs(w.imag).max()
        w = w.real
        # Joint clustering (between clock readings 14:47 and 14:52, for prop4_large_t.py). Cluster on (sign, |log|lam||), so a pole lam and
        # its partner 1/lam always fall in one joint group; then split each group at |lam| = 1 into its inside and
        # outside parts. The two sides are paired by construction. Before this, each side was clustered on its own, and
        # near the a5max tie at k1 = 0 the crowded top poles (gaps down to 5e-9 in log|lam| at L5 = 128) chained
        # differently on the two sides: AssertionError (63, 64), (128, 127), ... in 8 of 240 grid cases. The raw
        # eigenvalues pair to 3e-13 in log|lam| (prop4_diag.py). Where the two old clusterings agreed, the groups are the
        # same as before. self.spread records the widest merged group in log|lam|; replacing its poles by their mean
        # costs at most about t * spread relative at time t.
        sg, alg = np.sign(w), np.abs(np.log(np.abs(w)))
        self.poles = []          # (lam, Res D^-1)
        self.pairs = []
        self.spread = 0.0
        self.unit = []           # groups on |lam| = 1 to double precision; see below
        for g in cluster_idx(sg * np.exp(-alg)):
            g = np.asarray(g)
            gi, go = g[np.abs(w[g]) < 1], g[np.abs(w[g]) > 1]
            # (between clock readings 14:47 and 14:52) At m = 0.3, a5 = 1, k1 = 0, L5 = 128 (S = 128) the surface pair lam, 1/lam is
            # 1 -+ about 0.7^128 ~ 1e-20, below double precision, and both copies can land on one side: (2, 0). Such a
            # group is kept out of self.pairs and listed in self.unit; anything using the surface pole must check
            # self.unit is empty. prop4_large_t.py uses only the bulk, so it may proceed. Any other unequal split
            # still stops here.
            if len(gi) != len(go) and alg[g].max() < 1e-6:
                self.unit.append(w[g].copy())
                continue
            assert len(gi) == len(go), (len(gi), len(go))
            side = []
            for gs in (gi, go):
                lam = w[gs].mean()
                self.spread = max(self.spread, np.ptp(np.log(np.abs(w[gs]))))
                Q, _ = np.linalg.qr(V[:, gs])
                Hp = G0 @ (Bp - Bm / lam ** 2)
                Hp = (Hp + Hp.conj().T) / 2
                R = Q @ np.linalg.inv(Q.conj().T @ Hp @ Q) @ Q.conj().T
                self.poles.append((lam, R @ G0))
                side.append((lam, R @ G0))
            self.pairs.append((side[0][0], side[0][1], side[1][0], side[1][1]))
        self.inside = [(p[0], p[1]) for p in self.pairs]
        self.outside = [(p[2], p[3]) for p in self.pairs]

    def piece(self, sel):
        """Propagator tau -> G(tau) built from the selected pole pairs (indices into self.pairs)."""
        pr = [self.pairs[i] for i in sel]
        memo = {}

        def G(tau):
            if tau not in memo:
                if not pr:
                    memo[tau] = 0.0
                elif tau >= 1:
                    memo[tau] = sum(l ** (tau - 1) * R for l, R, _, _ in pr)
                else:
                    memo[tau] = -sum(a ** (tau - 1) * Ra for _, _, a, Ra in pr)
            return memo[tau]
        return G

    def full(self):
        return self.piece(range(len(self.pairs)))

    def fft(self, tau, N=4096):
        B0, Bp, Bm = self.B
        z = np.exp(2j * np.pi * (np.arange(N) + 0.5) / N)
        D = B0[None] + z[:, None, None] * Bp[None] + Bm[None] / z[:, None, None]
        return (np.linalg.inv(D) * (z ** tau)[:, None, None]).mean(axis=0)


def vertices(mu, row, k_out, k_in, L5, mirror):
    """List of (out_row, in_row, matrix) blocks of dD/da for a link; mirror temporal links run -ta -> 1 - ta."""
    E = np.eye(L5)
    if mu == 1:
        v = 1j * fw * np.exp(1j * k_in) - 1j * bw * np.exp(-1j * k_out)
        return [(row, row, np.kron(E, v))]
    return [(row, row + 1, np.kron(E, 1j * fw0)), (row + 1, row, np.kron(E, -1j * bw0))]


def C_pair(Gk, Gkp, k, p, ta, tb, L5):
    """2 x 2 kernel contribution of fermion momentum k (other leg k + p) at (ta, tb), before the 1/Ls."""
    C = np.zeros((2, 2), complex)
    for ie, mue in enumerate((1, 0)):
        rowe = 1 - ta if mue == 1 else -ta
        sig = 1.0 if mue == 1 else -1.0
        X = vertices(mue, rowe, k, k + p, L5, True)
        for jf, muf in enumerate((1, 0)):
            Y = vertices(muf, tb, k + p, k, L5, False)
            acc = 0.0
            for eo, ei, x in X:
                for fo, fi, y in Y:
                    A = Gk(fi - eo)
                    Bm = Gkp(ei - fo)
                    if np.isscalar(A) or np.isscalar(Bm):
                        continue
                    acc += np.trace(A @ x @ Bm @ y)
            C[ie, jf] = -sig * acc
    return C


def kvals(Ls, twist):
    return [(2 * np.pi * n + twist) / Ls for n in range(Ls)]


def validate():
    from check_prop13_kernel import kernel, blocks_C
    print("V. full kernel vs check_prop13_kernel (position space, open time)")
    for twist in (np.pi, 0.0):
        for mf in (0.0, 1.0):
            Ls, L5, m, a5, nt, T = 2, 2, 0.7, 0.5, 48, 3
            M, plus = kernel(nt, Ls, L5, m, a5, mf, T, twist)
            ks = kvals(Ls, twist)
            props = {round(k % (2 * np.pi), 12): Prop(L5, k, m, a5, mf) for k in ks}
            worst, scale = 0.0, 0.0
            for n in range(Ls):
                p1 = 2 * np.pi * n / Ls
                Cpos = blocks_C(M, plus, Ls, T, p1)
                for ta in range(1, T + 1):
                    for tb in range(1, T + 1):
                        Cm = sum(C_pair(props[round(k % (2 * np.pi), 12)].full(),
                                        props[round((k + p1) % (2 * np.pi), 12)].full(), k, p1, ta, tb, L5)
                                 for k in ks) / Ls
                        worst = max(worst, np.abs(Cm - Cpos[ta - 1, tb - 1]).max())
                        scale = max(scale, np.abs(Cpos[ta - 1, tb - 1]).max())
            pr = props[round(ks[0] % (2 * np.pi), 12)]
            fe = max(np.abs(pr.full()(tau) - pr.fft(tau)).max() for tau in (-3, -1, 0, 1, 2, 4)
                     if not np.isscalar(pr.full()(tau)))
            print("  twist=%.2f mf=%.0f: max |C_mine - C_pos| = %.2e (scale %.2e); pole-sum G vs FFT %.1e; imag %.0e"
                  % (twist, mf, worst, scale, fe, pr.imag))
    sys.stdout.flush()


def m_in(Ls, L5, m, a5, twist, p1s, ts, control=False):
    ks = kvals(Ls, twist)
    key = lambda k: round(k % (2 * np.pi), 12)  # noqa: E731
    DW = {key(k): Prop(L5, k, m, a5, 0.0) for k in ks}
    PV = {key(k): Prop(L5, k, m, a5, 1.0) for k in ks}
    info = dict(maxPV_over_rstar=max(max(l for l, *_ in PV[key(k)].pairs) / rstar(k, m) for k in ks),
                nsurf=0, imag=max(max(P.imag for P in DW.values()), max(P.imag for P in PV.values())))
    surf = {}
    for k in ks:
        P = DW[key(k)]
        s = [i for i, (l, *_) in enumerate(P.pairs) if l > rstar(k, m) * (1 + 1e-12)]
        assert len(s) <= 1, (k, s)
        assert not P.unit and not PV[key(k)].unit, (k, P.unit)   # (14:47-14:52) a dropped unit-circle group would hide the surface pole
        surf[key(k)] = s
        info["nsurf"] += len(s)
    out = {}
    for p1 in p1s:
        rho = rho_star(m, p1)
        Mbb = {t: np.zeros((2, 2), complex) for t in ts}
        Msurf = {t: np.zeros((2, 2), complex) for t in ts}
        Mout = {t: np.zeros((2, 2), complex) for t in ts}
        neg_surf, top_surf = 0.0, 0.0
        for k in ks:
            A, Ap = DW[key(k)], DW[key(k + p1)]
            sA, sAp = surf[key(k)], surf[key(k + p1)]
            bA = [i for i in range(len(A.pairs)) if i not in sA]
            bAp = [i for i in range(len(Ap.pairs)) if i not in sAp]
            Gb, Gbp = A.piece(bA), Ap.piece(bAp)
            Pv, Pvp = PV[key(k)].full(), PV[key(k + p1)].full()
            for t in ts:
                Mbb[t] += (C_pair(Gb, Gbp, k, p1, 1, t + 1, L5) - C_pair(Pv, Pvp, k, p1, 1, t + 1, L5)) / Ls
            # surface-involving nodes, one by one: W = C_node(1, 1), node mu = lam_i lam'_j
            nodes = [(i, j) for i in sA for j in range(len(Ap.pairs))] + [(i, j) for i in bA for j in sAp]
            for i, j in nodes:
                mu = A.pairs[i][0] * Ap.pairs[j][0]
                W = C_pair(A.piece([i]), Ap.piece([j]), k, p1, 1, 1, L5) / Ls
                ev = np.linalg.eigvalsh((W + W.conj().T) / 2)
                neg_surf = min(neg_surf, ev.min())  # against the largest weight (fixed at clock reading 14:07, the launch, made just before it)
                top_surf = max(top_surf, abs(ev).max())
                for t in ts:
                    if mu <= rho * (1 + 1e-12):
                        Msurf[t] += W * mu ** t
                    else:
                        Mout[t] += W * mu ** t
        out[p1] = dict(rho=rho, bb=Mbb, surf=Msurf, outer=Mout, neg_surf=neg_surf / max(top_surf, 1e-300))
    return out, info


def opn(M):
    return np.linalg.norm(M, 2)


def sweep(cfgs, ms, twists, p1s, ts, Ls=16, label=""):
    for m in ms:
        for twist in twists:
            for (L5, a5) in cfgs(m):
                out, info = m_in(Ls, L5, m, a5, twist, p1s, ts)
                r_tot, r_tot4, r_bb, r_s, r_ctrl = 0.0, 0.0, 0.0, 0.0, 0.0
                worst_neg = 0.0
                for p1, d in out.items():
                    rho = d["rho"]
                    for t in ts:
                        tot = d["bb"][t] + d["surf"][t]
                        r_tot = max(r_tot, opn(tot) / rho ** t)
                        r_tot4 = max(r_tot4, opn(tot) / ((t + 2) ** 4 * rho ** t))
                        r_bb = max(r_bb, opn(d["bb"][t]) / rho ** t)
                        r_s = max(r_s, opn(d["surf"][t]) / rho ** t)
                        r_ctrl = max(r_ctrl, opn(tot + d["outer"][t]) / rho ** t)
                    worst_neg = min(worst_neg, d["neg_surf"])
                print("%s m=%.1f twist=%.2f L5=%3d a5=%.4f S=%.2f: max||M_in||/rho^t = %.4f ; /((t+2)^4 rho^t) = %.3e ;"
                      " bulk-PV %.4f, surf-inner %.4f ; ctrl(i) with outer nodes %.3e ; surface node weights min eig"
                      " ratio %.1e ; max PV pole / r_* = %.4f ; n_surf %d"
                      % (label, m, twist, L5, a5, L5 * a5, r_tot, r_tot4, r_bb, r_s, r_ctrl, worst_neg,
                         info["maxPV_over_rstar"], info["nsurf"]))
                sys.stdout.flush()


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    ts = [0, 1, 2, 4, 8, 16, 32]
    p1s = [0.0, np.pi / 2, np.pi]
    if which in ("validate", "all"):
        validate()
    if which in ("main", "all"):
        print("M. sweep, Ls = 16, S = 3.2 (a5 halving)")
        sweep(lambda m: [(16, 0.2), (32, 0.1), (64, 0.05)], [0.3, 0.7, 0.9], [np.pi, 0.0], p1s, ts, label="M")
        print("M2. S = 6.4 at a5 = 0.2, 0.1")
        sweep(lambda m: [(32, 0.2), (64, 0.1)], [0.3, 0.7, 0.9], [np.pi, 0.0], p1s, ts, label="M2")
    if which in ("control", "all"):
        print("C(ii). a5 above a5max(m): m = 0.9 (a5max = %.3f), a5 = 0.5, 1.0" % a5max(0.9))
        sweep(lambda m: [(8, 0.5), (4, 1.0)], [0.9], [np.pi, 0.0], p1s, ts, label="C2")
    print("DONE")
