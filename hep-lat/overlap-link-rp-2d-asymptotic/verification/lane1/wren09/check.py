import numpy as np

def lam2_T(K0, K1, m):
    s2 = np.sin(K0)**2 + np.sin(K1)**2; b = (1 - np.cos(K0)) + (1 - np.cos(K1)) - m
    nn = (2 + b)**2 + s2
    bT = (b * (2 + b) + s2) / nn; sT2 = 4 * s2 / nn**2; wT = np.sqrt(bT**2 + sT2)
    lam2 = np.where(bT < 0, sT2 / (2 * wT * (wT - bT)), (1 + bT / wT) / 2)
    return lam2, wT

def grid(n):
    k = np.linspace(-np.pi, np.pi, n, endpoint=False) + np.pi / n
    return np.meshgrid(k, k, indexing="ij")

m = 0.7
vals = {}
for n in (512, 1024, 2048, 4096):
    K0, K1 = grid(n); l2, _ = lam2_T(K0, K1, m); vals[n] = -np.mean(np.log(l2))
    print("N1 even n", n, "%.9f" % vals[n])
r = [(4 * vals[2 * n] - vals[n]) / 3 for n in (512, 1024, 2048)]
print("N1 Richardson", ["%.9f" % x for x in r])
K0, K1 = grid(4001); l2, _ = lam2_T(K0, K1, m)
i = np.unravel_index(np.argmin(K0**2 + K1**2), K0.shape)
print("N2 4001: smallest |k| =", np.hypot(K0[i], K1[i]), " lam2 there =", l2[i], " log term/n^2 =", -np.log(l2[i]) / 4001**2)
full = -np.mean(np.log(l2)); mask = np.ones_like(l2, bool); mask[i] = False
print("N2 4001 full %.9f  without that point (mean over the rest) %.9f" % (full, -np.mean(np.log(l2[mask]))))
for mm in (0.5, 0.3, 1.0):
    K0, K1 = grid(2000); _, w = lam2_T(K0 + 0.0, K1 + 0.0, mm)
    k = np.linspace(-np.pi, np.pi, 2001); A, B = np.meshgrid(k, k, indexing="ij"); _, w2 = lam2_T(A, B, mm)
    f4 = np.sqrt((mm**2 + 8 * (1 - mm)) / (mm**2 + 8 * (1 - mm) + 4 * (5 - mm)))
    print("N3 m=%.1f min|H_T| (grid with k=0) %.7f  m/(2-m) %.7f  f(4)-bound %.4f  Wren mu_T %.7f" % (mm, w2.min(), mm / (2 - mm), f4, min(mm / (2 - mm), f4)))
muT = np.sqrt(2.89 / 20.09)
print("N4 mu_T(0.7) = %.7f" % muT)
print("N5 2sqrt2/(mu_T 1.69) = %.5f ; sqrt2/0.7 = %.5f ; eta_T = mu_T*1.69/8 = %.6f" % (2 * np.sqrt(2) / (muT * 1.69), np.sqrt(2) / 0.7, muT * 1.69 / 8))
