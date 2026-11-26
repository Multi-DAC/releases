"""Does the L4 inner moment grow with S = L5 a5?  (written after the M2 rows, clock reading 14:17)

M2 (S = 6.4) gave max ||M_in|| / rho_*^t about 1.7x the S = 3.2 values, with bulk - PV doubling. Two S values cannot tell
saturation from growth. Here, per configuration and per p1:
  - bulk - PV, the surface-involving inner piece, and the surface-involving OUTER piece (mu > rho_*), at t = 0, 1, 2;
  - the full DW - PV kernel (all nodes) at the same t, which must converge in S if the method is right;
  - max over t in {0..32} of ||M_in(t)|| / rho_*^t.
Also control (ii) made direct: above a5max, how many DW poles per k sit above r_*(k), and where the PV poles sit.
"""
import sys

import numpy as np

from l4_assembled import C_pair, Prop, a5max, kvals, opn, rho_star, rstar

TS = [0, 1, 2, 4, 8, 16, 32]


def pieces(Ls, L5, m, a5, twist, p1):
    ks = kvals(Ls, twist)
    key = lambda k: round(k % (2 * np.pi), 12)  # noqa: E731
    DW = {key(k): Prop(L5, k, m, a5, 0.0) for k in ks}
    PV = {key(k): Prop(L5, k, m, a5, 1.0) for k in ks}
    surf = {key(k): [i for i, (l, *_) in enumerate(DW[key(k)].pairs) if l > rstar(k, m) * (1 + 1e-12)] for k in ks}
    rho = rho_star(m, p1)
    z = lambda: {t: np.zeros((2, 2), complex) for t in TS}  # noqa: E731
    bmp, sin, sout, full = z(), z(), z(), z()
    for k in ks:
        A, Ap = DW[key(k)], DW[key(k + p1)]
        sA, sAp = surf[key(k)], surf[key(k + p1)]
        bA = [i for i in range(len(A.pairs)) if i not in sA]
        bAp = [i for i in range(len(Ap.pairs)) if i not in sAp]
        Gb, Gbp = A.piece(bA), Ap.piece(bAp)
        Pv, Pvp = PV[key(k)].full(), PV[key(k + p1)].full()
        Fa, Fap = A.full(), Ap.full()
        for t in TS:
            pvt = C_pair(Pv, Pvp, k, p1, 1, t + 1, L5) / Ls
            bmp[t] += C_pair(Gb, Gbp, k, p1, 1, t + 1, L5) / Ls - pvt
            if t <= 2:
                full[t] += C_pair(Fa, Fap, k, p1, 1, t + 1, L5) / Ls - pvt
        nodes = [(i, j) for i in sA for j in range(len(Ap.pairs))] + [(i, j) for i in bA for j in sAp]
        for i, j in nodes:
            mu = A.pairs[i][0] * Ap.pairs[j][0]
            W = C_pair(A.piece([i]), Ap.piece([j]), k, p1, 1, 1, L5) / Ls
            tgt = sin if mu <= rho * (1 + 1e-12) else sout
            for t in TS:
                tgt[t] += W * mu ** t
    return rho, bmp, sin, sout, full


def run(m, twist, cfgs, p1s=(0.0, np.pi / 2, np.pi), Ls=16):
    for L5, a5 in cfgs:
        worst = 0.0
        lines = []
        for p1 in p1s:
            rho, bmp, sin, sout, full = pieces(Ls, L5, m, a5, twist, p1)
            worst = max(worst, max(opn(bmp[t] + sin[t]) / rho ** t for t in TS))
            for t in (0, 1, 2):
                resid = opn(full[t] - bmp[t] - sin[t] - sout[t]) / max(opn(full[t]), 1e-300)
                lines.append("    p1=%.2f t=%d: |bulk-PV| %.4f |surf-in| %.4f |surf-out| %.4f |full DW-PV| %.4f"
                             " (split residual %.1e)" % (p1, t, opn(bmp[t]), opn(sin[t]), opn(sout[t]), opn(full[t]),
                                                          resid))
        print("m=%.1f twist=%.2f L5=%3d a5=%.2f S=%5.1f: max_t ||M_in||/rho^t = %.4f" % (m, twist, L5, a5, L5 * a5, worst))
        print("\n".join(lines))
        sys.stdout.flush()


def control_ii(m, cfgs, Ls=16):
    for twist in (np.pi, 0.0):
        for L5, a5 in cfgs:
            above, pvmax = [], 0.0
            for k in kvals(Ls, twist):
                D = Prop(L5, k, m, a5, 0.0)
                P = Prop(L5, k, m, a5, 1.0)
                above.append(sum(1 for l, *_ in D.pairs if l > rstar(k, m) * (1 + 1e-12)))
                pvmax = max(pvmax, max(l for l, *_ in P.pairs) / rstar(k, m))
            print("C(ii) m=%.1f (a5max %.3f) twist=%.2f L5=%d a5=%.2f: DW poles above r_* per k = %s ; max PV pole / r_* = %.4f"
                  % (m, a5max(m), twist, L5, a5, above, pvmax))
            sys.stdout.flush()


if __name__ == "__main__":
    control_ii(0.9, [(8, 0.5), (4, 1.0), (16, 0.25), (16, 0.21)])
    run(0.3, np.pi, [(16, 0.2), (32, 0.2), (64, 0.2), (64, 0.4), (64, 0.8)])
    run(0.9, np.pi, [(16, 0.2), (32, 0.2), (64, 0.2)])
    run(0.9, 0.0, [(16, 0.2), (32, 0.2), (64, 0.2)])
    print("DONE")
