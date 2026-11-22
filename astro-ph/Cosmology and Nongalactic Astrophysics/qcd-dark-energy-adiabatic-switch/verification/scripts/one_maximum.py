"""Sections 2.4 and 2.8: every extremum of u = rho_DE/rho_dS over 0 <= z <= 1090, with radiation, at the fitted
quadrature, Lorentzian and general-q points. Reads the fitted parameters from c43_fit_all.json (written by c43_fit.py).

Usage: python -B one_maximum.py     (run c43_fit.py first; see DATA.md for the data files)
"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import json
import numpy as np
import c43_fit as F

A = json.load(open(os.path.join(HERE, 'c43_fit_all.json')))
for name in ('quad', 'lor', 'genq'):
    p = A[name]['params']
    kw = dict(p=2, q=p.get('q', 0.5 if name == 'quad' else 1.0))
    for zsw in (20.0, 400.0):
        H, u = F.H_switch(p['ob'], p['oc'], p['Hbar'], p['c'], zsw=zsw, **kw)
        du = np.diff(u); sg = np.sign(du)
        idx = np.where(sg[1:] * sg[:-1] < 0)[0] + 1
        ext = [(round(float(F.ZG[i]), 3), 'max' if u[i] > u[i - 1] else 'min', round(float(u[i]), 4)) for i in idx]
        print(f"{name}: c {p['c']:.3f} q {kw['q']:.3f} hand-over z={zsw:.0f}: extrema of u on the whole grid (z, kind, u): "
              f"{ext}  u at z=1090: {u[-1]:.3e}  grid steps with du = 0: {int(np.sum(du == 0))}")
