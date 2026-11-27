import numpy as np
from lat import *
from audit_checks import rp_kernel, report, m, mf
rng = np.random.default_rng(99)
for L, L5, spread, n in [(6, 4, 1.0, 30), (8, 2, 0.7, 20), (8, 3, 1.0, 16)]:
    halves = [rand_half(L, rng, spread) for _ in range(n)]
    lAP = logdet(dwf(L, U_AP(L), m, mf, L5))[1]
    lAPov = logdet(overlap_HT(L, U_AP(L), m, mf))[1]
    def w_dw(U):
        s, l = logdet(dwf(L, U, m, mf, L5)); return s * np.exp(l - lAP)
    def w_ov(U):
        s, l = logdet(overlap_HT(L, U, m, mf)); return s * np.exp(l - lAPov)
    print(f" L={L} L5={L5} spread={spread} n={n}", flush=True)
    report("DWF numerator det D_DW (claim: RP)", rp_kernel(L, halves, w_dw))
    report("overlap det D_ov (control)", rp_kernel(L, halves, w_ov))
