"""Free H_T (a5 = 1) constants for extending Lemma Z: c_f(H_T), the free gap, and the small-|lambda0| constant."""
import numpy as np
m = 0.7
n = 4001; k = np.linspace(-np.pi, np.pi, n, endpoint=False) + np.pi / n
K0, K1 = np.meshgrid(k, k, indexing="ij")
s2 = np.sin(K0)**2 + np.sin(K1)**2; b = (1 - np.cos(K0)) + (1 - np.cos(K1)) - m
nn = (2 + b)**2 + s2
bT = (b * (2 + b) + s2) / nn; sT2 = 4 * s2 / nn**2; wT = np.sqrt(bT**2 + sT2)
lam2 = np.where(bT < 0, sT2 / (2 * wT * (wT - bT)), (1 + bT / wT) / 2)
print("c_f(H_T) midpoint n=4001:", -np.mean(np.log(lam2)))
w = np.sqrt(b**2 + s2); lamW = np.where(b < 0, s2 / (2 * w * (w - b)), (1 + b / w) / 2)
print("c_f(H_W) same grid      :", -np.mean(np.log(lamW)))
print("free gap min |H_T| =", wT.min(), " at k=0 predicted m/(2-m) =", m / (2 - m))
k2 = K0**2 + K1**2
neg = bT < 0
print("min over b_T<0 of |lambda0|^2/|k|^2 =", (lam2[neg] / k2[neg]).min(), "; min over b_T>=0 of |lambda0|^2 =", lam2[~neg].min())
print("max |k_mu| on b_T<0:", np.abs(K0[neg]).max())
