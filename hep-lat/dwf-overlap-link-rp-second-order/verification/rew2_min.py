"""Check of the merged text's line 797 figure (Wren: min Re det X = 9.44e-3 for m in [0.1, 0.99] on |z| = r_mid).
After clock reading 16:21. Exact minimum is at k = 0, z = r_mid real, where Re w^2 = c1(m)."""
import numpy as np
from letter78_tests import rmax_of
from l4_assembled import rstar

th = np.linspace(0, 2 * np.pi, 721)
ks = np.linspace(-np.pi, np.pi, 721)
for m in (0.1, 0.11, 0.12, 0.15, 0.2, 0.3, 0.5, 0.586, 0.7, 0.9, 0.99):
    rmid = (1 + rmax_of(m)) / 2
    z = rmid * np.exp(1j * th)
    lo = min((((2 - m - np.cos(k)) ** 2 + 1 + np.sin(k) ** 2) - (2 - m - np.cos(k)) * (z + 1 / z)).real.min() for k in ks)
    print("m=%.3f r_mid %.4f: min Re det X on |z| = r_mid = %.4e" % (m, rmid, lo))
