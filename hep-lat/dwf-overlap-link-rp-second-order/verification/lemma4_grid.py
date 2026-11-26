"""Lemma 4 (ii), numerical check: min over (k1, k5) of E_PV minus min_q E_cut(q), on a grid of m and a5.

From the derivation in Lemma 4 of the paper: with c = 2 - m - cos k1 and u5 = (1 - cos k5)/a5 in [0, 2/a5],
    cosh E_PV = (c^2 + 2 c u5 + 1 + sin^2 k1 + 2 u5/a5) / (2 (c + u5)),
and min_q E_cut(q) = min(-ln(1 - m), ln(3 - m)). Lemma 4 (ii) says the difference is >= 0 for a5 <= 1.
Grids: k1 and k5 on 721 points each over [-pi, pi]; m = 0.05, 0.10, ..., 0.95; a5 = 0.05, 0.10, ..., 1.00 and 1.05.
"""
import numpy as np

k = np.linspace(-np.pi, np.pi, 721)
K1, K5 = np.meshgrid(k, k, indexing="ij")
worst_le1 = np.inf
rows = []
for m in np.round(np.arange(0.05, 0.951, 0.05), 2):
    emin = min(-np.log(1 - m), np.log(3 - m))
    for a5 in np.round(np.arange(0.05, 1.051, 0.05), 2):
        c = 2 - m - np.cos(K1)
        u5 = (1 - np.cos(K5)) / a5
        ch = (c * c + 2 * c * u5 + 1 + np.sin(K1) ** 2 + 2 * u5 / a5) / (2 * (c + u5))
        d = float(np.min(np.arccosh(ch)) - emin)
        if a5 <= 1.0:
            worst_le1 = min(worst_le1, d)
        else:
            rows.append((m, a5, d))
print("min over m, a5 <= 1, (k1, k5) of E_PV - min_q E_cut: %.2e" % worst_le1)
for m, a5, d in rows:
    print("a5 = %.2f  m = %.2f  min diff = %+.3e" % (a5, m, d))
print("DONE")
