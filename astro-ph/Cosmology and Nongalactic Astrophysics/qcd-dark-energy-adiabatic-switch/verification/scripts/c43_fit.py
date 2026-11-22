"""C43 replication (Clawd, 2026-11-21). Fits LCDM, w0wa and QCD-DE switch models to
CMB (Bansal-Huterer compressed, arXiv 2502.07185 eqs cmb_vec, covariance_matrix) + DES-Dovekie SNe + DESI DR2 BAO.
Written from the public files and the papers' TeX; none of Wren's code or numbers were read.

Usage: python -B c43_fit.py [model ...]     (no argument: all nine; prints one JSON line per model and writes
c43_fit_<models>.json beside this file). The data files are not shipped: see DATA.md, and set C43_DATA.

Switch models, units Hbar = 1, N = ln a, h = H/Hbar, mu = (8 pi G/3) rho_t / Hbar^2:
  numer 'H'   : h^2 = mu + h * beta(x)      (2606.20036 L325, rho_DE prop. to Hbar H beta)
  numer 'Hbar': h^2 = mu + beta(x)          (2506.14182 eq schwinger_adia read literally)
  x = |Hdot|/(H Hbar) = -dh/dN ;  beta = (1 + (x/c)^p)^(-q):  linear p=1,q=1; quadrature p=2,q=1/2; Lorentzian p=2,q=1.
Integrated backward from near de Sitter (the stable direction).
"""
import os, sys, json, numpy as np
from scipy import constants as K
from scipy.integrate import solve_ivp, quad, cumulative_simpson
from scipy.interpolate import CubicSpline
from scipy.optimize import least_squares

HERE = os.path.dirname(os.path.abspath(__file__))
# the four public data files (see DATA.md): directory named by the environment variable C43_DATA, default ./data
D = os.environ.get('C43_DATA', os.path.join(HERE, 'data')).rstrip('/') + '/'
CKM = K.c / 1e3
TCMB, NEFF, MNU = 2.7255, 3.044, 0.06
MPC = 3.0856775814913673e22
rho100 = 3 * (100e3 / MPC) ** 2 / (8 * np.pi * K.G)
OG = (4 * K.sigma / K.c) * TCMB ** 4 / K.c ** 2 / rho100          # omega_gamma
TNU_EV = (4 / 11) ** (1 / 3) * TCMB * K.k / K.e
ONUM = MNU / 93.14
ZSTAR = 1090.0

def onu_a4(a):
    """omega_nu(a) * a^4 (Bansal-Huterer App. B, Komatsu fitting formula)."""
    y = MNU * a / TNU_EV
    f = (1 + (0.3173 * y) ** 1.83) ** (1 / 1.83)
    return 0.2271 * OG * NEFF * (2 + f) / 3

def om_tot(a, ob, oc):
    """omega of matter + photons + neutrinos at scale factor a (units of (100 km/s/Mpc)^2)."""
    return (ob + oc) / a ** 3 + (OG + onu_a4(a)) / a ** 4

# ---------- data ----------
bz, bv, bq = [], [], []
for line in open(D + 'desi_gaussian_bao_ALL_GCcomb_mean.txt'):
    if line.startswith('#') or not line.strip(): continue
    z, v, q = line.split(); bz.append(float(z)); bv.append(float(v)); bq.append(q)
bz, bv = np.array(bz), np.array(bv)
bcov = np.loadtxt(D + 'desi_gaussian_bao_ALL_GCcomb_cov.txt')
assert bcov.shape == (13, 13) and len(bz) == 13
LB = np.linalg.cholesky(np.linalg.inv(bcov))

sn = [l.split() for l in open(D + 'DES-Dovekie_HD.csv') if l.startswith('SN:')]
zhd = np.array([float(r[3]) for r in sn]); zhel = np.array([float(r[4]) for r in sn]); mu_obs = np.array([float(r[5]) for r in sn])
d = np.load(D + 'STAT+SYS.npz'); n = int(d[d.files[0]][0])
assert n == len(zhd) == 1820
icov = np.zeros((n, n)); icov[np.triu_indices(n)] = d[d.files[1]]
il = np.tril_indices(n, -1); icov[il] = icov.T[il]
LS = np.linalg.cholesky(icov)
ONE_S = LS.T @ np.ones(n); CS = ONE_S @ ONE_S

cmb_v = np.array([1.7504, 301.77, 0.022371])
cmb_c = 1e-8 * np.array([[1559.83, -1325.41, -36.45], [-1325.41, 714691.80, 269.77], [-36.45, 269.77, 2.10]])
LC = np.linalg.cholesky(np.linalg.inv(cmb_c))

# ---------- background ----------
SG = np.linspace(0, np.log(1 + ZSTAR), 8001)      # s = ln(1+z)
ZG = np.expm1(SG); AG = np.exp(-SG)

def H_lcdm(ob, oc, H0, w0=-1.0, wa=0.0):
    ode = (H0 / 100) ** 2 - om_tot(1.0, ob, oc)
    de = ode * AG ** (-3 * (1 + w0 + wa)) * np.exp(-3 * wa * (1 - AG))
    return 100 * np.sqrt(om_tot(AG, ob, oc) + de), ode * AG ** (-3 * (1 + w0 + wa)) * np.exp(-3 * wa * (1 - AG))

def dom_tot(a, ob, oc):
    """d om_tot / d ln a, analytic."""
    y = MNU * a / TNU_EV; Ay = (0.3173 * y) ** 1.83
    df = Ay * (1 + Ay) ** (1 / 1.83 - 1)
    return -3 * (ob + oc) / a ** 3 - 4 * (OG + onu_a4(a)) / a ** 4 + 0.2271 * OG * NEFF * df / 3 / a ** 4

ZSW = 20.0

def H_switch(ob, oc, Hbar, c, p=1, q=1.0, numer='H', Ns=3.0, rtol=1e-10, zsw=None):
    """Dependent variable x = -dh/dN > 0.  h is algebraic in (x, mu): no subtraction of near-equal numbers for rho_DE.
    ODE backward from N = Ns to z = zsw; above zsw the quasi-static (slaved) solution, iterated on the grid."""
    zsw = ZSW if zsw is None else zsw
    hb2 = (Hbar / 100) ** 2
    mu = lambda N: om_tot(np.exp(N), ob, oc) / hb2
    dmu = lambda N: dom_tot(np.exp(N), ob, oc) / hb2
    f = lambda x: (1 + (x / c) ** p) ** (-q)
    fp = lambda x: -q * p * (x / c) ** (p - 1) / c * (1 + (x / c) ** p) ** (-q - 1)
    if numer == 'H':
        hof = lambda x, m: 0.5 * (f(x) + np.sqrt(f(x) ** 2 + 4 * m))
        rhs = lambda N, x: (x * f(x) - 2 * hof(x, mu(N)) * x - dmu(N)) / (hof(x, mu(N)) * fp(x))
    else:
        hof = lambda x, m: np.sqrt(m + f(x))
        rhs = lambda N, x: (-2 * hof(x, mu(N)) * x - dmu(N)) / fp(x)
    m = mu(Ns); k = 1.0 / c if p == 1 else 0.0
    if numer == 'H':
        a1 = 1 / (1 + 3 * k); a2 = (-1 - 9 * q / c ** 2) if p == 2 else 0.0
    else:
        a1 = 1 / (2 + 3 * k); a2 = (-(1 / 8) - 9 * q / (8 * c ** 2)) if p == 2 else 0.0
    x0 = 3 * a1 * m + 6 * a2 * m * m
    isw = int(np.searchsorted(ZG, zsw))
    sol = solve_ivp(rhs, (Ns, -SG[isw]), [x0], method='LSODA', t_eval=np.concatenate(([Ns], -SG[:isw + 1])),
                    rtol=rtol, atol=1e-14)
    if not sol.success or len(sol.y[0]) != isw + 2: raise RuntimeError('ode: ' + str(sol.message))
    x = np.empty(len(SG)); x[:isw + 1] = sol.y[0][1:]
    if np.any(x[:isw + 1] <= 0): raise RuntimeError('x <= 0')
    mg = mu(-SG); dmg = dmu(-SG)
    # slaved part: x = -(mu' + u')/(2h), u = rho_DE/rho_dS, iterated; u' by finite differences on the grid
    hi = slice(isw + 1, len(SG))
    u = np.zeros(len(SG)); xs = np.empty(len(SG))
    h = hof(x[:isw + 1], mg[:isw + 1])
    for it in range(4):
        hh = np.sqrt(mg + u)
        du = np.gradient(u, -SG) if it else np.zeros(len(SG))
        xs = -(dmg + du) / (2 * hh)
        u = (hh if numer == 'H' else 1.0) * f(xs)
    x[hi] = xs[hi]
    h_full = hof(x, mg)
    u_full = (h_full if numer == 'H' else 1.0) * f(x)
    return Hbar * h_full, u_full

def observables(H, ob, oc):
    integ = (1 + ZG) / H
    DM = CKM * cumulative_simpson(integ, x=SG, initial=0.0)
    sDM = CubicSpline(SG, DM); sH = CubicSpline(SG, H)
    # BAO
    rd = 147.05 * (ob / 0.02236) ** -0.13 * ((ob + oc) / 0.1432) ** -0.23 * (NEFF / 3.04) ** -0.1
    out = []
    for z, qn in zip(bz, bq):
        s = np.log1p(z); dm = sDM(s); dh = CKM / sH(s)
        out.append({'DM_over_rs': dm, 'DH_over_rs': dh, 'DV_over_rs': (z * dm * dm * dh) ** (1 / 3)}[qn] / rd)
    bao = np.array(out)
    # SN
    mu_th = 5 * np.log10((1 + zhel) * sDM(np.log1p(zhd))) + 25
    # CMB: sound horizon at z*, matter + radiation only above z* (DE there is < 1e-8 of the total in every model here)
    Rb = 3 * ob / (4 * OG)
    f = lambda a: 1 / np.sqrt(3 * (1 + Rb * a)) / (100 * np.sqrt(om_tot(a, ob, oc)) * a * a)
    rs = CKM * quad(f, 0, 1 / (1 + ZSTAR), epsabs=0, epsrel=1e-10, limit=200)[0]
    dms = DM[-1]
    cmb = np.array([100 * np.sqrt(ob + oc + ONUM) * dms / CKM, np.pi * dms / rs, ob])
    return bao, mu_th, cmb, rd, rs

def resid(H, ob, oc):
    bao, mu_th, cmb, rd, rs = observables(H, ob, oc)
    rb = LB.T @ (bao - bv)
    ds = LS.T @ (mu_th - mu_obs)
    rsn = ds - ONE_S * (ONE_S @ ds) / CS            # M marginalised (projected out)
    rc = LC.T @ (cmb - cmb_v)
    return rb, rsn, rc

MODELS = {
    'lcdm':    dict(names=['ob', 'oc', 'H0'], x0=[0.0224, 0.119, 68.0], H=lambda t: H_lcdm(*t)),
    'w0wa':    dict(names=['ob', 'oc', 'H0', 'w0', 'wa'], x0=[0.0224, 0.119, 67.0, -0.8, -0.7], H=lambda t: H_lcdm(*t)),
    'lin_c1':  dict(names=['ob', 'oc', 'Hbar'], x0=[0.0224, 0.119, 62.0], H=lambda t: H_switch(*t, 1.0)),
    'lin':     dict(names=['ob', 'oc', 'Hbar', 'c'], x0=[0.0224, 0.119, 58.0, 3.0], H=lambda t: H_switch(*t)),
    'litlin_c1': dict(names=['ob', 'oc', 'Hbar'], x0=[0.0224, 0.119, 56.0], H=lambda t: H_switch(*t, 1.0, numer='Hbar')),
    'litlin':  dict(names=['ob', 'oc', 'Hbar', 'c'], x0=[0.0224, 0.119, 56.0, 3.0], H=lambda t: H_switch(*t, numer='Hbar')),
    'quad':    dict(names=['ob', 'oc', 'Hbar', 'c'], x0=[0.0225, 0.1185, 51.0, 1.65], H=lambda t: H_switch(*t, p=2, q=0.5)),
    'lor':     dict(names=['ob', 'oc', 'Hbar', 'c'], x0=[0.0224, 0.1199, 50.0, 2.8], H=lambda t: H_switch(*t, p=2, q=1.0)),
    'genq':    dict(names=['ob', 'oc', 'Hbar', 'c', 'q'], x0=[0.0224, 0.1199, 50.0, 2.8, 1.0],
                    H=lambda t: H_switch(t[0], t[1], t[2], t[3], p=2, q=t[4])),
}

def fit(name, x0=None, verbose=False):
    M = MODELS[name]
    x0 = np.array(M['x0'] if x0 is None else x0, float)
    def r(t):
        try:
            H, u = M['H'](t)
            if not np.all(np.isfinite(H)): raise RuntimeError('nan')
            return np.concatenate(resid(H, t[0], t[1]))
        except Exception as e:
            return np.full(13 + 1820 + 3, 1e3)
    lo = [0.015, 0.05, 30.0] + {'w0wa': [-3, -5], 'lin': [0.05], 'litlin': [0.05], 'quad': [0.3], 'lor': [0.3], 'genq': [0.3, 0.05]}.get(name, [])
    hi = [0.03, 0.2, 100.0] + {'w0wa': [1, 3], 'lin': [200], 'litlin': [200], 'quad': [200], 'lor': [200], 'genq': [200, 20]}.get(name, [])
    res = least_squares(r, x0, bounds=(lo, hi), x_scale=np.abs(x0) * 0.01 + 1e-6, xtol=1e-12, ftol=1e-12, gtol=1e-12, diff_step=1e-5)
    t = res.x
    H, u = M['H'](t)
    rb, rsn, rc = resid(H, t[0], t[1])
    out = dict(model=name, k=len(t), params=dict(zip(M['names'], map(float, t))),
               chi2_bao=float(rb @ rb), chi2_sn=float(rsn @ rsn), chi2_cmb=float(rc @ rc))
    out['chi2'] = out['chi2_bao'] + out['chi2_sn'] + out['chi2_cmb']
    out['H0'] = float(H[0])
    w = -1 - np.gradient(np.log(np.abs(u)), -SG) / 3
    out['w_z0'] = float(w[0])
    i1 = np.searchsorted(ZG, 3.0)
    # a crossing of w = -1 between grid points; rounding noise about w = -1 exactly (LCDM) is not one
    cr = [float(ZG[i]) for i in range(1, i1) if (w[i] + 1) * (w[i + 1] + 1) < 0 and max(abs(w[i] + 1), abs(w[i + 1] + 1)) > 1e-8]
    out['w_crossings_z_lt_3'] = cr
    out['w_at'] = {str(z): float(w[np.searchsorted(ZG, z)]) for z in (0.3, 0.5, 1.0, 2.0)}
    out['ODE_at_zstar'] = float(u[-1] * ((t[2] / 100) ** 2 if 'Hbar' in M['names'] else 1) / om_tot(AG[-1], t[0], t[1]))
    return out

if __name__ == '__main__':
    names = sys.argv[1:] or list(MODELS)
    allr = {}
    for nm in names:
        o = fit(nm); allr[nm] = o
        print(json.dumps(o), flush=True)
    label = '_'.join(sys.argv[1:])[:40] or 'all'
    json.dump(allr, open(os.path.join(HERE, 'c43_fit_' + label + '.json'), 'w'), indent=1)
