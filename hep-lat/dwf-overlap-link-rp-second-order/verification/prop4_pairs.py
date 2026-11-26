"""Is Prop 4's cancellation pole-by-pole, after pairing? (written after clock reading 14:56)

Group the bulk DW poles and the PV poles from the top down, matching cumulative residue rank (PV poles are doubly
degenerate in k5 -> -k5, so one rank-2 PV pole meets two simple DW poles). S1 = sum over groups of
||a5 (sum R_DW - sum R_PV)||_1 at t = 0. If S1 stays bounded as L5 grows, a pairing proof of O(1) r_*^t is open; if
S1 grows like L5 while ||a5 Delta(1)||_1 stays O(1), the cancellation is collective, across groups, and only a
contour (or transfer-matrix) argument can see it. Also S1top(n): the same for the top n groups only.
"""
import numpy as np

from l4_assembled import Prop, rstar


def tn(A):
    return np.linalg.svd(A, compute_uv=False).sum()


def run(L5, k1, m, a5):
    rs = rstar(k1, m)
    DW, PV = Prop(L5, k1, m, a5, 0.0), Prop(L5, k1, m, a5, 1.0)
    rk = lambda R: int(np.linalg.matrix_rank(R, tol=1e-9 * max(np.abs(R).max(), 1e-300)))
    dw = sorted([(l, R, rk(R)) for l, R, _, _ in DW.pairs if l <= rs * (1 + 1e-12)], key=lambda x: -x[0])
    pv = sorted([(l, R, rk(R)) for l, R, _, _ in PV.pairs], key=lambda x: -x[0])
    groups, i, j = [], 0, 0
    while i < len(dw) or j < len(pv):
        gd, gp, cd, cp = [], [], 0, 0
        while True:
            if (cd <= cp and i < len(dw)) or j >= len(pv):
                gd.append(dw[i]); cd += dw[i][2]; i += 1
            else:
                gp.append(pv[j]); cp += pv[j][2]; j += 1
            if cd == cp or (i >= len(dw) and j >= len(pv)):
                break
        groups.append((gd, gp))
    diffs = [tn(a5 * (sum(R for _, R, _ in gd) if gd else 0) - a5 * (sum(R for _, R, _ in gp) if gp else 0))
             for gd, gp in groups]
    total = tn(a5 * sum(R for _, R, _ in dw) - a5 * sum(R for _, R, _ in pv))
    mass = sum(tn(a5 * R) for _, R, _ in dw)
    return dict(ngroups=len(groups), S1=sum(diffs), top=[float("%.3g" % d) for d in diffs[:5]], total=total, mass=mass,
                pvranks=sorted(set(r for _, _, r in pv)), dwranks=sorted(set(r for _, _, r in dw)))


for m, a5, k1 in [(0.7, 0.2, np.pi / 2), (0.7, 0.2, 0.0), (0.3, 1.0, np.pi)]:
    for L5 in (16, 32, 64, 128):
        d = run(L5, k1, m, a5)
        print("m=%.1f a5=%.2f k1=%.3f L5=%3d: groups %d (DW ranks %s, PV ranks %s); ||a5 Delta(1)||_1 = %.4f;"
              " S1 = sum_groups ||diff||_1 = %.3f; DW trace mass %.1f; top-5 group diffs %s"
              % (m, a5, k1, L5, d["ngroups"], d["dwranks"], d["pvranks"], d["total"], d["S1"], d["mass"], d["top"]),
              flush=True)
print("DONE")
