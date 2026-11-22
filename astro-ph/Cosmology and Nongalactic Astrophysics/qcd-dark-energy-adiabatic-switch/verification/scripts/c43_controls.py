"""Controls on c43_fit.py: (1) the Planck 2018 LCDM point against the compressed CMB data vector; (2) the grid distance
against direct quadrature; (3) the switch ODE against its starting point and tolerance; (4) against the hand-over
redshift between the ODE and the quasi-static solution.

Usage: python -B c43_controls.py     (see DATA.md for the data files)
"""
import os, sys, time, numpy as np; sys.stdout.reconfigure(line_buffering=True)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import c43_fit as F
print('omega_gamma', F.OG, 'Tnu eV', F.TNU_EV, 'onu today', F.onu_a4(1.0), 'mnu/93.14', F.ONUM)
# control 1: Planck 2018 LCDM point -> R, l_a near the BH data vector
H,_ = F.H_lcdm(0.02237, 0.1200, 67.36)
bao, mu, cmb, rd, rs = F.observables(H, 0.02237, 0.1200)
print('Planck point: R, l_a =', cmb[:2], 'rd', rd, 'rs*', rs)
# control 2: grid distance vs direct quad at z=0.5, 2.33
from scipy.integrate import quad
E = lambda z: 100*np.sqrt(F.om_tot(1/(1+z),0.02237,0.12) + (0.6736**2 - F.om_tot(1,0.02237,0.12)))
for z in (0.025, 0.5, 2.33, 1090):
    dq = F.CKM*quad(lambda x: 1/E(x), 0, z, epsrel=1e-12, limit=500)[0]
    from scipy.interpolate import CubicSpline
    from scipy.integrate import cumulative_simpson
    DM = F.CKM*cumulative_simpson((1+F.ZG)/H, x=F.SG, initial=0.0)
    print(z, dq, float(CubicSpline(F.SG, DM)(np.log1p(z))), 'rel', float(CubicSpline(F.SG, DM)(np.log1p(z)))/dq-1)
# control 3: switch ODE, start-point and tolerance sensitivity
for kw in (dict(c=2.0,p=1,q=1.0), dict(c=2.8,p=2,q=1.0), dict(c=1.8,p=2,q=0.5), dict(c=1.0,p=1,q=1.0,numer='Hbar')):
    t=time.time(); H0,u0 = F.H_switch(0.0224,0.119,60.0,**kw); dt=time.time()-t
    H1,_ = F.H_switch(0.0224,0.119,60.0,Ns=2.5,**kw); H2,_ = F.H_switch(0.0224,0.119,60.0,Ns=3.5,rtol=1e-12,**kw)
    print(kw, 'H0', H0[0], 'sec', round(dt,2), 'max rel Ns2.5', np.max(np.abs(H1/H0-1)), 'Ns3.5+tight', np.max(np.abs(H2/H0-1)), 'u(z*)', u0[-1])
# control 4: switch-redshift sensitivity (slaved part against ODE part)
for kw in (dict(c=2.0,p=1,q=1.0), dict(c=2.8,p=2,q=1.0)):
    Ha,_ = F.H_switch(0.0224,0.119,60.0,zsw=20,**kw); Hb,_ = F.H_switch(0.0224,0.119,60.0,zsw=60,**kw); Hc,_ = F.H_switch(0.0224,0.119,60.0,zsw=8,**kw)
    print(kw, 'zsw 20 vs 60 max rel', np.max(np.abs(Hb/Ha-1)), '20 vs 8', np.max(np.abs(Hc/Ha-1)))
