"""Tests of Wren's letter 78 with my own operator (written after clock reading 15:41; first read sealed at 15:41,
letter78-firstread.md, sha 946785697d5573a5). Her script is not run.

T1: Table T. Is the plateau a diffusive Hankel kernel (Wren) or an atom (me)? Top singular vector's spread in s.
T2: Bug B step 4. At k1 = 0, the DW poles with |lam| >= r_mid = (1 + r_max)/2: exactly one (the surface pole), in
    (r_mid, 1), once S is large? PV poles all below r_mid?
T3: Letter 78 item 3.2. The DW surface pole at k1 != 0 tends to lam0_ov = exp(-E0), sinh E0 = |sin k1|, as a5 -> 0
    at fixed S; and lam0_ov <= exp(-c0 |k1|).
T4: Bug B (a). Re w^2 >= c1(m) and |b| < 5 on the annulus, on a grid.
"""
import numpy as np

from l4_assembled import Prop, rstar


def rmax_of(m):
    """max over k of r_*(k, m): 1 - m at k = 0 for m < 2 - sqrt 2, 1/(3 - m) at k = pi above (found 15:42)."""
    return rstar(np.linspace(-np.pi, np.pi, 20001), m).max()

np.set_printoptions(precision=4, suppress=True)


def T1(L5=128, m=0.3, a5=1.0, k1=np.pi):
    rs = rstar(k1, m)
    DW, PV = Prop(L5, k1, m, a5, 0.0), Prop(L5, k1, m, a5, 1.0)
    bulk = [(l, R) for l, R, _, _ in DW.pairs if l <= rs * (1 + 1e-12)]
    pv = [(l, R) for l, R, _, _ in PV.pairs]
    s = np.arange(L5)
    print("T1  m=%.1f a5=%.1f k1=pi L5=%d (diffusive window 1 << t << 1/(1 - top/r*) ~ %.0f)"
          % (m, a5, L5, 1 / (1 - max(l for l, _ in bulk) / rs)))
    for t in (20, 50, 100, 200, 400, 800, 1600):
        A = a5 * (sum((l / rs) ** t * R for l, R in bulk) - sum((l / rs) ** t * R for l, R in pv))
        U, sv, Vh = np.linalg.svd(A)
        out = []
        for vec in (U[:, 0], Vh[0].conj()):
            w = np.abs(vec[0::2]) ** 2 + np.abs(vec[1::2]) ** 2
            out.append((w @ s, w @ np.minimum(s, L5 - 1 - s), np.sqrt(w @ (np.minimum(s, L5 - 1 - s)) ** 2)))
        e00 = np.linalg.norm(A[0:2, 0:2], 2)
        eLL = np.linalg.norm(A[-2:, -2:], 2)
        print("  t=%5d: trace norm %.4f; sv %s; |A_00| %.5f, |A_LL| %.5f, sqrt(t)*max %.4f; top left vec: <s> %.1f,"
              " <wall dist> %.2f, rms %.2f; right: <s> %.1f, <wall dist> %.2f; <wall dist>/sqrt(t) %.3f"
              % (t, sv.sum(), sv[:4], e00, eLL, np.sqrt(t) * max(e00, eLL), out[0][0], out[0][1], out[0][2],
                 out[1][0], out[1][1], out[0][1] / np.sqrt(t)), flush=True)


def T2():
    print("T2  k1 = 0: DW and PV poles at |lam| >= r_mid = (1 + r_max)/2")
    for m in (0.3, 0.7, 0.9):
        rmid = (1 + rmax_of(m)) / 2
        for a5 in (0.1, 0.05):
            for S in (2, 4, 8, 16):
                L5 = int(round(S / a5))
                DW, PV = Prop(L5, 0.0, m, a5, 0.0), Prop(L5, 0.0, m, a5, 1.0)
                hi = [l for l, *_ in DW.pairs if l >= rmid]
                print("  m=%.1f a5=%.2f S=%4.1f L5=%3d: DW inside poles >= r_mid: %d %s; unit groups %d; max other DW "
                      "pole/r_* %.6f; max PV/r_mid %.4f" % (m, a5, S, L5, len(hi), ["%.8f" % l for l in hi],
                                                             len(DW.unit), max(l for l, *_ in DW.pairs if l < rmid)
                                                             / rstar(0.0, m), max(l for l, *_ in PV.pairs) / rmid),
                      flush=True)


def lam0_ov(k1, m):
    c = 2 - m - np.cos(k1)
    if np.sqrt(1 + np.sin(k1) ** 2) <= c:
        return None
    return np.exp(-np.arcsinh(abs(np.sin(k1))))


def T3():
    print("T3  DW surface pole vs the overlap pole exp(-E0), sinh E0 = |sin k1| (letter 78 item 3.2), S = 8")
    for m in (0.3, 0.7):
        c0 = min(0.5, np.arccosh(2 - m) / np.pi)
        for k1 in (np.pi / 16, np.pi / 8, np.pi / 4, np.pi / 2):
            l0 = lam0_ov(k1, m)
            row = []
            for a5 in (0.2, 0.1, 0.05, 0.025):
                L5 = int(round(8 / a5))
                P = Prop(L5, k1, m, a5, 0.0)
                s = [l for l, *_ in P.pairs if l > rstar(k1, m) * (1 + 1e-12)]
                row.append("a5=%.3f: %s" % (a5, "%.6f" % s[0] if len(s) == 1 else "n=%d" % len(s)))
            print("  m=%.1f k1=%.3f: overlap pole %s (r_* %.4f; exp(-c0|k1|) %.4f) ; DW surface pole %s"
                  % (m, k1, "%.6f" % l0 if l0 else "none", rstar(k1, m), np.exp(-c0 * k1), "; ".join(row)),
                  flush=True)


def T4():
    print("T4  Bug B (a): min Re(w^2)/c1 and max |b| on the annulus r_mid <= |z| <= 1/r_mid")
    for m in (0.05, 0.3, 0.7, 0.9, 0.99):
        rmax = rmax_of(m)
        ks = np.linspace(-np.pi, np.pi, 721)
        assert rstar(ks, m).max() <= rmax + 1e-12
        rmid = (1 + rmax) / 2
        c1 = (1 - m) * ((rmax + 1 / rmax) - (rmid + 1 / rmid))
        rho = np.exp(np.linspace(np.log(rmid), -np.log(rmid), 81))
        th = np.linspace(0, 2 * np.pi, 361)
        z = (rho[:, None] * np.exp(1j * th)[None, :]).ravel()
        lo, bmax = np.inf, 0.0
        for k in ks:
            c = 2 - m - np.cos(k)
            w2 = c * c + 1 + np.sin(k) ** 2 - c * (z + 1 / z)
            lo = min(lo, (w2.real / c1).min())
            bmax = max(bmax, np.abs(c - (z + 1 / z) / 2).max())
        print("  m=%.2f: r_max %.4f r_mid %.4f c1 %.4g: min Re(w^2)/c1 = %.4f ; max |b| = %.3f" % (m, rmax, rmid, c1, lo,
                                                                                               bmax), flush=True)


if __name__ == "__main__":
    T4()
    T1()
    T2()
    T3()
    print("DONE")
