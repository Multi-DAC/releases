"""Section 2.6, the linear switch with c > 3 and radiation:
E(z) = mu_r/(3 sqrt(mu_m)) + 4 mu_r/(3 mu_m) <= c - 3  =>  every critical point of u from then on is a strict minimum,
hence u' < 0 (w > -1) for all later times. Prints the redshift at which E = c - 3 (the time N_e of the note) and, from
the integration, the located maximum of u and the sign of u' at every grid point below the exclusion redshift.

Parameters as in sec26_table.py. Usage: python -B sec26_exclusion.py   (needs the data files only because c43_fit.py
loads them on import; see DATA.md)
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from scipy.optimize import brentq
import c43_fit as F

ob, oc, Hbar = 0.0224, 0.1196, 55.0
hb2 = (Hbar / 100) ** 2
mu_m0 = (ob + oc) / hb2
mu_r0 = (F.OG + F.onu_a4(1e-6)) / hb2
print(f"mu_m0 {mu_m0:.5f}  mu_r0 {mu_r0:.4e}")


def E(zp1):
    mm, mr = mu_m0 * zp1 ** 3, mu_r0 * zp1 ** 4
    return mr / (3 * np.sqrt(mm)) + 4 * mr / (3 * mm)


for c in (3.001, 3.01, 3.1, 3.5, 4.0, 6.0, 10.0, 20.0):
    ze = brentq(lambda s: E(s) - (c - 3), 1e-6, 1e5)
    H, u = F.H_switch(ob, oc, Hbar, c, p=1, Ns=3.0)
    zg = np.exp(F.SG)                      # 1 + z on the grid
    du = np.diff(u) / np.diff(F.SG)        # du/d ln(1+z); u' in N has the opposite sign
    mid = np.sqrt(zg[1:] * zg[:-1])
    inside = mid <= ze
    imax = [j for j in range(1, len(u) - 1) if u[j] > u[j - 1] and u[j] > u[j + 1]]
    if imax:
        j = imax[0]
        a, b, _ = np.polyfit(F.SG[j - 1:j + 2], u[j - 1:j + 2], 2)
        zmax = float(np.exp(-b / (2 * a)))
    else:
        zmax = float('nan')
    print(f"c={c}: E = c-3 at 1+z = {ze:.2f}; maxima of u on the grid: {len(imax)}, the first at 1+z = {zmax:.2f}; "
          f"grid points with 1+z below the bound: {inside.sum()}, of which du/dln(1+z) > 0 (w > -1): {(du[inside] > 0).sum()}")
