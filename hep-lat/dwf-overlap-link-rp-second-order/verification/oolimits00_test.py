"""Order-of-limits example, (0,0) temporal and (1,1) spatial entries of the k1 = 0 mode (after clock reading 17:32).
Sealed plan: oolimits00-plan.md (sha256 prefix f943f2a25ada635d). My operator only (check_letter77: build2, hessians).
Ls = 1 with twist 0 is exactly the k1 = 0 mode; twist pi (k1 = pi, gapped) is the control.
Reflection (merged line 655): theta a1(t) = a1(1 - t), theta a0(t) = -a0(-t). Row r < 0 is index nt + r."""
import sys
import time

import numpy as np

sys.path.insert(0, ".")
from check_letter77 import build2, G5, hessians  # noqa: E402

m, ta, tb = 0.7, 5, 5


def Dmaker(nt, twist, mu):
    def f(U1, U0):
        X = build2(nt, 1, m, U1, U0, twist, "ap")
        G = G5(len(X))
        HW = G @ X
        w, V = np.linalg.eigh((HW + HW.conj().T) / 2)
        return 0.5 * (np.eye(len(X)) + G @ (V * np.sign(w)) @ V.conj().T) + mu * np.eye(len(X))
    return f


def Q(nt, twist, mu):
    links = [(1, nt + 1 - ta, 0), (0, nt - ta, 0), (1, tb, 0), (0, tb, 0)]   # a1 mirror, a0 mirror, a1 tb, a0 tb
    H = hessians(Dmaker(nt, twist, mu), nt, 1, links)
    # Q_{mu nu}: mu = component of the reflected (x) side, nu = of the y side; temporal reflection carries a minus
    return np.array([[H[0, 2], H[0, 3]], [-H[1, 2], -H[1, 3]]])   # rows: theta-side a1, a0; columns: a1, a0 at tb


def show(label, q):
    print("%-34s Q_11 %+.6f  Q_10 %+.2e  Q_01 %+.2e  Q_00 %+.6f  (max |Im| %.1e)" % (
        label, q[0, 0].real, abs(q[0, 1]), abs(q[1, 0]), q[1, 1].real, np.abs(q.imag).max()), flush=True)


t0 = time.time()
for twist, name in ((0.0, "twist 0 (k1 = 0)"), (np.pi, "CONTROL twist pi (k1 = pi)")):
    show("AP, exact overlap, nt=64, " + name, Q(64, twist, 0.0))
    irs = []
    for eps in (0.04, 0.02, 0.01):
        q = Q(1024, twist, eps / (2 * m))
        irs.append(q)
        show("IR, eps=%.2f, nt=1024, %s" % (eps, name), q)
    e = np.array([0.04, 0.02, 0.01])
    for (i, j), lab in (((0, 0), "Q_11"), ((1, 1), "Q_00")):
        y = np.array([q[i, j].real for q in irs])
        c = np.polyfit(e, y, 2)
        print("  IR extrapolated to eps = 0 (quadratic in eps), %s: %+.6f" % (lab, c[-1]), flush=True)
    print("  elapsed %.0f s" % (time.time() - t0), flush=True)
print("DONE")
