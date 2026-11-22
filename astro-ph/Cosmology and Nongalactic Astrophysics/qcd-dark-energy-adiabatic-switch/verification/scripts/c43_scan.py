"""Profile of chi^2 in c for one switch family, the other parameters refitted at each c.

Usage: python -B c43_scan.py lin|quad|lor|litlor     (see DATA.md for the data files)
"""
import os, sys, json, numpy as np; sys.stdout.reconfigure(line_buffering=True)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import warnings; warnings.filterwarnings('ignore')
import c43_fit as F
from scipy.optimize import least_squares
def prof(p, q, numer, cs, x0=(0.0224, 0.119, 60.0)):
    out = []
    x = np.array(x0, float)
    for c in cs:
        def r(t):
            try:
                H, u = F.H_switch(t[0], t[1], t[2], c, p=p, q=q, numer=numer)
                return np.concatenate(F.resid(H, t[0], t[1]))
            except Exception:
                return np.full(1836, 1e3)
        res = least_squares(r, x, bounds=([0.015, 0.05, 30], [0.03, 0.2, 100]), x_scale=[2e-4, 1e-3, 0.5], xtol=1e-11, ftol=1e-11, gtol=1e-11, diff_step=1e-5)
        x = res.x
        out.append((c, 2 * res.cost, *res.x))
        print(f'p={p} q={q} {numer} c={c:.3f} chi2={2*res.cost:.3f} ob={x[0]:.5f} oc={x[1]:.5f} Hbar={x[2]:.3f}')
    return out
which = sys.argv[1]
if which == 'lin': prof(1, 1.0, 'H', [1, 1.5, 2, 2.25, 2.5, 2.75, 3, 3.5, 4, 6, 10, 30])
if which == 'quad': prof(2, 0.5, 'H', [1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 3, 4, 6, 10])
if which == 'lor': prof(2, 1.0, 'H', [1.5, 2.0, 2.5, 2.75, 3.0, 3.25, 3.5, 4, 5, 7, 10, 20])
if which == 'litlor': prof(2, 1.0, 'Hbar', [1.0, 1.5, 2.0, 3.0, 4, 6, 10, 20, 50])
