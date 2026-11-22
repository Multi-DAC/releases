"""C43 primary & secondary cosmological fits, 1-sigma c profiles, BDF robustness, and DESI BAO pulls (Wren, 2026-11-21).

Independent u-evolution implementation (BDF rtol=1e-11, atol=1e-13 + CubicSpline distances and exact ODE RHS w(z))
producing Table 1, CMB rows of Table 2 (CMB(BH)+PP+DESI, CMB(BH)+DESI, CMB(DESI)+DD+DESI), 1-sigma profile CIs on c,
Section 3.3 BDF numerical robustness checks, and DESI DR2 diagonal BAO pulls (data - model) / sigma_i.

Data directory is read from environment variable C43_DATA (default: ./data beside this script, falling back to ../data).
See DATA.md for the six public data files and SHA-256 hashes.
"""
import os
import time
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline
from scipy.optimize import minimize, brentq

HERE = os.path.dirname(os.path.abspath(__file__))
_default_data = os.path.join(HERE, 'data')
if not os.path.isdir(_default_data) and os.path.isdir(os.path.join(HERE, '..', 'data')):
    _default_data = os.path.join(HERE, '..', 'data')
DATA = os.environ.get('C43_DATA', _default_data)
C_KMS = 299792.458

# Load datasets
bao_z, bao_val, bao_type = [], [], []
with open(os.path.join(DATA, 'desi_gaussian_bao_ALL_GCcomb_mean.txt')) as f:
    for line in f:
        s = line.strip()
        if not s or s.startswith('#'):
            continue
        p = s.split()
        bao_z.append(float(p[0]))
        bao_val.append(float(p[1]))
        bao_type.append(p[2])
bao_z = np.array(bao_z)
bao_val = np.array(bao_val)
bao_cov = np.loadtxt(os.path.join(DATA, 'desi_gaussian_bao_ALL_GCcomb_cov.txt'))
bao_icov = np.linalg.inv(bao_cov)
bao_sig = np.sqrt(np.diag(bao_cov))

dd_rows = []
with open(os.path.join(DATA, 'DES-Dovekie_HD.csv')) as f:
    header = None
    for line in f:
        s = line.strip()
        if not s or s.startswith('#'):
            continue
        if s.startswith('VARNAMES:'):
            header = s.split()[1:]
            continue
        if s.startswith('SN:'):
            dd_rows.append(s.split()[1:])
dd_zHD = np.array([float(r[header.index('zHD')]) for r in dd_rows])
dd_zHEL = np.array([float(r[header.index('zHEL')]) for r in dd_rows])
dd_mu = np.array([float(r[header.index('MU')]) for r in dd_rows])
d_npz = np.load(os.path.join(DATA, 'STAT+SYS.npz'))
n_dd = int(d_npz[d_npz.files[0]][0])
dd_icov = np.zeros((n_dd, n_dd))
dd_icov[np.triu_indices(n_dd)] = d_npz[d_npz.files[1]]
i_lo = np.tril_indices(n_dd, -1)
dd_icov[i_lo] = dd_icov.T[i_lo]
dd_ones = np.sum(dd_icov, axis=0)
dd_sum = float(np.sum(dd_ones))

pp_cols = None
pp_rows = []
with open(os.path.join(DATA, 'Pantheon+SH0ES.dat')) as f:
    for line in f:
        p = line.strip().split()
        if not p:
            continue
        if pp_cols is None:
            pp_cols = p
            continue
        pp_rows.append(p)
p_zHD = np.array([float(r[pp_cols.index('zHD')]) for r in pp_rows])
p_zHEL = np.array([float(r[pp_cols.index('zHEL')]) for r in pp_rows])
p_mb = np.array([float(r[pp_cols.index('m_b_corr')]) for r in pp_rows])
with open(os.path.join(DATA, 'Pantheon+SH0ES_STAT+SYS.cov')) as f:
    n_pp = int(f.readline().strip())
    pp_cov_full = np.fromfile(f, sep=' ', count=n_pp * n_pp).reshape(n_pp, n_pp)
mask_pp = p_zHD > 0.01
pp_zHD = p_zHD[mask_pp]
pp_zHEL = p_zHEL[mask_pp]
pp_mb = p_mb[mask_pp]
pp_icov = np.linalg.inv(pp_cov_full[np.ix_(mask_pp, mask_pp)])
pp_ones = np.sum(pp_icov, axis=0)
pp_sum = float(np.sum(pp_ones))

# Compressed CMB priors: Bansal & Huterer (2025, PRD 112, 023528, Eqs. 5-6) and DESI DR2 (2503.14738, App. A Eqs. 35-36)
bh_mu = np.array([1.7504, 301.77, 0.022371])
bh_icov = np.linalg.inv(1e-8 * np.array([
    [1559.83, -1325.41, -36.45],
    [-1325.41, 714691.80, 269.77],
    [-36.45, 269.77, 2.10],
]))
desi_cmb_mu = np.array([0.01041, 0.02223, 0.14208])
desi_cmb_icov = np.linalg.inv(1e-9 * np.array([
    [0.006621, 0.12444, -1.1929],
    [0.12444, 21.344, -94.001],
    [-1.1929, -94.001, 1488.4],
]))

OMEGA_GAMMA_H2 = 2.47282e-5
F_NU_1 = (7.0 / 8.0) * (4.0 / 11.0) ** (4.0 / 3.0) * (3.044 / 3.0)
OMEGA_NU_1_REL = F_NU_1 * OMEGA_GAMMA_H2
OMEGA_RAD_MASSLESS = OMEGA_GAMMA_H2 * (1.0 + 2.0 * F_NU_1)
OMEGA_NU_MASSIVE_0 = 0.06 / 93.14
XS64, WS64 = np.polynomial.legendre.leggauss(64)


def omega_nu_m(a, rad_on=True):
    if not rad_on:
        return OMEGA_NU_MASSIVE_0 / a ** 3
    p = 1.83
    return ((OMEGA_NU_MASSIVE_0 / a ** 3) ** p + (OMEGA_NU_1_REL / a ** 4) ** p) ** (1.0 / p)


def d_omega_nu_m_dN(a, rad_on=True):
    if not rad_on:
        return -3.0 * OMEGA_NU_MASSIVE_0 / a ** 3
    p = 1.83
    t3 = (OMEGA_NU_MASSIVE_0 / a ** 3) ** p
    t4 = (OMEGA_NU_1_REL / a ** 4) ** p
    return -(3.0 * t3 + 4.0 * t4) * (t3 + t4) ** (1.0 / p - 1.0)


def omega_bg_h2(a, wbc, rad_on=True):
    return wbc / a ** 3 + (OMEGA_RAD_MASSLESS / a ** 4 if rad_on else 0.0) + omega_nu_m(a, rad_on)


def d_omega_bg_h2_dN(a, wbc, rad_on=True):
    return -3.0 * wbc / a ** 3 - (4.0 * OMEGA_RAD_MASSLESS / a ** 4 if rad_on else 0.0) + d_omega_nu_m_dN(a, rad_on)


def rd_desi_eq2(wb, wbc):
    return 147.05 * (wb / 0.02236) ** (-0.13) * (wbc / 0.1432) ** (-0.23) * (3.044 / 3.04) ** (-0.1)


def r_star_int(wb, wbc, z_star=1090.0):
    x_end = 1.0 / np.sqrt(1.0 + z_star)
    x = 0.5 * x_end * (XS64 + 1.0)
    w = 0.5 * x_end * WS64
    a = x * x
    cs = 1.0 / np.sqrt(3.0 * (1.0 + (3.0 * wb / (4.0 * OMEGA_GAMMA_H2)) * a))
    denom = np.sqrt(wbc * a + OMEGA_RAD_MASSLESS + (a ** 4) * omega_nu_m(a))
    return float((C_KMS / 100.0) * np.sum(w * 2.0 * x * cs / denom))


def dm_highz(wbc, z_split=10.0, z_star=1090.0, u_hi_h2=0.0):
    x0 = 1.0 / np.sqrt(1.0 + z_star)
    x1 = 1.0 / np.sqrt(1.0 + z_split)
    x = 0.5 * (x1 - x0) * XS64 + 0.5 * (x1 + x0)
    w = 0.5 * (x1 - x0) * WS64
    a = x * x
    return float(C_KMS * np.sum(w * 2.0 / (x ** 3 * 100.0 * np.sqrt(omega_bg_h2(a, wbc) + u_hi_h2))))


N_FUT = np.linspace(4.5, 0.0, 141)
N_PAST = np.linspace(0.0, -np.log(11.0), 201)
N_ALL = np.concatenate([N_FUT[:-1], N_PAST])
A_ALL = np.exp(N_ALL)
Z_PAST = np.exp(-N_PAST) - 1.0
A_PAST = np.exp(N_PAST)
RAD_ALL = OMEGA_RAD_MASSLESS / A_ALL ** 4 + omega_nu_m(A_ALL)
DRAD_ALL = -4.0 * OMEGA_RAD_MASSLESS / A_ALL ** 4 + d_omega_nu_m_dN(A_ALL)
INV_A3_ALL = 1.0 / A_ALL ** 3


def rhs_and_dfdu(model, u, mu, dmu, c, q=1.0):
    if u <= 1e-30:
        return 0.0, 0.0
    h = (mu + u) ** 0.5
    if model == 'schwinger':
        if u >= 1.0:
            return -dmu, 0.0
        x = c * (1.0 / u - 1.0)
        return -dmu - 2.0 * h * x, -(x / h) + 2.0 * h * c / (u * u)
    y = u / h
    if y >= 1.0 - 1e-15:
        return -dmu, 0.0
    dy_du = (2.0 * mu + u) / (2.0 * h * (mu + u))
    if model == 'linear':
        x = c * (1.0 / y - 1.0)
        dx_dy = -c / (y * y)
    elif model == 'quadrature':
        sq = max(1.0 / (y * y) - 1.0, 1e-30) ** 0.5
        x = c * sq
        dx_dy = -c / (y * y * y * sq)
    elif model == 'lorentzian':
        sq = max(1.0 / y - 1.0, 1e-30) ** 0.5
        x = c * sq
        dx_dy = -0.5 * c / (y * y * sq)
    elif model == 'gaussian':
        sq = max(-np.log(max(y, 1e-300)), 1e-30) ** 0.5
        x = c * sq
        dx_dy = -0.5 * c / (max(y, 1e-300) * sq)
    elif model == 'gen_lor':
        invq = 1.0 / q
        ym = y ** (-invq)
        sq = max(ym - 1.0, 1e-30) ** 0.5
        x = c * sq
        dx_dy = -0.5 * c * invq * (ym / y) / sq
    else:
        raise ValueError(f'Unknown model: {model}')
    return -dmu - 2.0 * h * x, -(x / h) - 2.0 * h * dx_dy * dy_du


def solve_fast(model, hbar, wbc, c, q=1.0):
    if hbar <= 0.15 or hbar >= 1.5 or wbc <= 0.05 or wbc >= 0.35 or c <= 0.02 or c >= 5000.0:
        return None
    inv_hb2 = 1.0 / (hbar * hbar)
    mu_all = (wbc * INV_A3_ALL + RAD_ALL) * inv_hb2
    dmu_all = (-3.0 * wbc * INV_A3_ALL + DRAD_ALL) * inv_hb2
    mu_s = mu_all[0]
    if model == 'schwinger':
        u = 1.0 - (3.0 / (2.0 * c + 3.0)) * mu_s
    elif model == 'linear':
        u = 1.0 + ((c - 3.0) / (c + 3.0)) * mu_s
    elif model == 'quadrature':
        u = 1.0 + mu_s - (1.0 + 9.0 / (c * c)) * mu_s ** 2
    elif model in ('lorentzian', 'gaussian'):
        u = 1.0 + mu_s - (1.0 + 18.0 / (c * c)) * mu_s ** 2
    else:
        u = 1.0 + mu_s - (1.0 + 18.0 * q / (c * c)) * mu_s ** 2
    u_arr = np.empty(len(N_ALL))
    u_arr[0] = u
    f_prev, _ = rhs_and_dfdu(model, u, mu_all[0], dmu_all[0], c, q)
    for i in range(1, len(N_ALL)):
        hdN = 0.5 * (N_ALL[i] - N_ALL[i - 1])
        if u <= 1e-25:
            u = 0.0
            u_arr[i] = 0.0
            f_prev = 0.0
            continue
        mu_i = mu_all[i]
        dmu_i = dmu_all[i]
        un = u
        for _ in range(8):
            fi, dfi = rhs_and_dfdu(model, un, mu_i, dmu_i, c, q)
            step = ((un - u) - hdN * (f_prev + fi)) / (1.0 - hdN * dfi)
            un = 0.25 * un if (un - step <= 0.0) else (un - step)
            if un <= 1e-25:
                un = 0.0
                break
            if abs(step) < 1e-11 * max(1.0, un):
                break
        if not np.isfinite(un):
            return None
        u = un
        u_arr[i] = u
        f_prev, _ = rhs_and_dfdu(model, u, mu_i, dmu_i, c, q) if u > 0 else (0.0, 0.0)
    up = u_arr[140:]
    if up[0] <= 1e-6:
        return None
    hp = np.sqrt(mu_all[140:] + up)
    return N_PAST, Z_PAST, hp, up, 100.0 * hbar * hp[0]


def solve_bdf(model, hbar, wbc, c, q=1.0, N_start=4.5, rtol=1e-11, atol=1e-13, rad_on=True, z_max=10.0):
    inv_hb2 = 1.0 / (hbar * hbar)
    mu_s = omega_bg_h2(np.exp(N_start), wbc, rad_on) * inv_hb2
    if model == 'schwinger':
        u_s = 1.0 - (3.0 / (2.0 * c + 3.0)) * mu_s
    elif model == 'linear':
        u_s = 1.0 + ((c - 3.0) / (c + 3.0)) * mu_s
    elif model == 'quadrature':
        u_s = 1.0 + mu_s - (1.0 + 9.0 / (c * c)) * mu_s ** 2
    elif model in ('lorentzian', 'gaussian'):
        u_s = 1.0 + mu_s - (1.0 + 18.0 / (c * c)) * mu_s ** 2
    else:
        u_s = 1.0 + mu_s - (1.0 + 18.0 * q / (c * c)) * mu_s ** 2

    def rhs(N, ua):
        a = np.exp(N)
        mu = omega_bg_h2(a, wbc, rad_on) * inv_hb2
        dmu = d_omega_bg_h2_dN(a, wbc, rad_on) * inv_hb2
        return [rhs_and_dfdu(model, max(ua[0], 0.0), mu, dmu, c, q)[0]]

    def jac(N, ua):
        a = np.exp(N)
        mu = omega_bg_h2(a, wbc, rad_on) * inv_hb2
        dmu = d_omega_bg_h2_dN(a, wbc, rad_on) * inv_hb2
        return [[rhs_and_dfdu(model, max(ua[0], 0.0), mu, dmu, c, q)[1]]]

    N_ev = N_PAST if z_max == 10.0 else np.linspace(0.0, -np.log(1.0 + z_max), 401)
    sol = solve_ivp(rhs, (N_start, N_ev[-1]), [u_s], method='BDF', jac=jac, t_eval=N_ev, rtol=rtol, atol=atol)
    if (not sol.success) or len(sol.y[0]) != len(N_ev):
        return solve_fast(model, hbar, wbc, c, q)
    u_ev = np.maximum(sol.y[0], 0.0)
    a_ev = np.exp(N_ev)
    z_ev = np.exp(-N_ev) - 1.0
    h_ev = np.sqrt(omega_bg_h2(a_ev, wbc, rad_on) * inv_hb2 + u_ev)
    return N_ev, z_ev, h_ev, u_ev, 100.0 * hbar * h_ev[0]


def eval_exact(model, params, use_cmb='BH', sn_set='DD', use_bdf=False):
    wb, wbc = params[0], params[1]
    if wb < 0.019 or wb > 0.026 or wbc < 0.09 or wbc > 0.20:
        return (1e12, 0, 0, 0, {})
    rd = rd_desi_eq2(wb, wbc)
    if model in ('lcdm', 'w0wa', 'logistic'):
        h = params[2]
        de_p = params[3:]
        if h < 0.50 or h > 0.85:
            return (1e12, 0, 0, 0, {})
        om_de_h2 = h * h - omega_bg_h2(1.0, wbc)
        if om_de_h2 <= 0.05:
            return (1e12, 0, 0, 0, {})
        if model == 'lcdm':
            H_g = 100.0 * np.sqrt(omega_bg_h2(A_PAST, wbc) + om_de_h2)
            u_hi_h2 = om_de_h2
            zc_val = None
        elif model == 'w0wa':
            w0, wa = de_p[0], de_p[1]
            E2 = omega_bg_h2(A_PAST, wbc) + om_de_h2 * A_PAST ** (-3.0 * (1.0 + w0 + wa)) * np.exp(-3.0 * wa * (1.0 - A_PAST))
            if np.any(E2 <= 0):
                return (1e12, 0, 0, 0, {})
            H_g = 100.0 * np.sqrt(E2)
            u_hi_h2 = om_de_h2 * (1 / 11.0) ** (-3.0 * (1.0 + w0 + wa)) * np.exp(-3.0 * wa * (10.0 / 11.0))
            a_c = 1.0 + (1.0 + w0) / wa if abs(wa) > 1e-12 else -1.0
            zc_val = (1.0 / a_c - 1.0) if (0.0 < a_c < 1.0) else None
        elif model == 'logistic':
            zq, dzq = de_p[0], de_p[1]
            if dzq < 0.03:
                return (1e12, 0, 0, 0, {})
            Aq = 1.0 + np.exp((-1.0 - zq) / dzq)
            b0 = Aq / (1.0 + np.exp(-zq / dzq))
            hb_v = om_de_h2 / (h * b0)
            hb = hb_v * (Aq / (1.0 + np.exp((Z_PAST - zq) / dzq)))
            H_g = 100.0 * 0.5 * (hb + np.sqrt(hb * hb + 4.0 * omega_bg_h2(A_PAST, wbc)))
            u_hi_h2 = 0.0
            zc_val = None
        info = {'H0': 100 * h, 'Om': (wbc + OMEGA_NU_MASSIVE_0) / (h * h), 'zc': zc_val}
    else:
        hb = params[2]
        c = params[3]
        q = params[4] if model == 'gen_lor' else 1.0
        if model == 'gen_lor' and (q < 0.3 or q > 5.0):
            return (1e12, 0, 0, 0, {})
        out = solve_bdf(model, hb, wbc, c, q) if use_bdf else solve_fast(model, hb, wbc, c, q)
        if out is None:
            return (1e12, 0, 0, 0, {})
        N_ev, z_ev, h_ev, u_ev, H0_c = out
        if H0_c < 50 or H0_c > 85:
            return (1e12, 0, 0, 0, {})
        H_g = 100.0 * hb * h_ev
        u_hi_h2 = u_ev[-1] * hb * hb
        a_ev = np.exp(N_ev)
        inv_hb2 = 1.0 / (hb * hb)
        mu_ev = omega_bg_h2(a_ev, wbc, True) * inv_hb2
        dmu_ev = d_omega_bg_h2_dN(a_ev, wbc, True) * inv_hb2
        du_ev = np.array([rhs_and_dfdu(model, max(ui, 0.0), mi, dmi, c, q)[0] for ui, mi, dmi in zip(u_ev, mu_ev, dmu_ev)])
        u_s = np.maximum(u_ev, 1e-30)
        w_ev = -1.0 - du_ev / (3.0 * u_s)
        val = u_ev > 1e-6
        wv = w_ev[val]
        nv = N_ev[val]
        duv = du_ev[val]
        ic = np.where(np.diff(np.sign(wv + 1.0)))[0]
        zc_val = None
        if len(ic):
            idx0 = ic[0]
            cs_du = CubicSpline(nv[::-1], duv[::-1])
            n_lo, n_hi = float(nv[idx0 + 1]), float(nv[idx0])
            n_root = brentq(lambda n: float(cs_du(n)), n_lo, n_hi, xtol=1e-10)
            zc_val = float(np.exp(-n_root) - 1.0)
        info = {
            'H0': H0_c,
            'Om': (wbc + OMEGA_NU_MASSIVE_0) / (H0_c / 100) ** 2,
            'Hbar': 100 * hb,
            'w0': float(w_ev[0]),
            'w05': float(np.interp(0.5, z_ev, w_ev)),
            'w1': float(np.interp(1.0, z_ev, w_ev)),
            'zc': zc_val,
        }

    cs_iH = CubicSpline(Z_PAST, C_KMS / H_g)
    cs_dm = cs_iH.antiderivative()
    dm_at = lambda z: cs_dm(z) - cs_dm(0.0)
    H_at = lambda z: C_KMS / cs_iH(z)

    pred = np.empty(len(bao_z))
    for i, (zb, tp) in enumerate(zip(bao_z, bao_type)):
        dm = float(dm_at(zb))
        dh = C_KMS / float(H_at(zb))
        if tp == 'DM_over_rs':
            pred[i] = dm / rd
        elif tp == 'DH_over_rs':
            pred[i] = dh / rd
        else:
            pred[i] = (zb * dm * dm * dh) ** (1.0 / 3.0) / rd
    db = pred - bao_val
    c_bao = float(db @ bao_icov @ db)
    info['bao_pulls_data_minus_model'] = (bao_val - pred) / bao_sig

    c_sn = 0.0
    if sn_set == 'DD':
        d = 5.0 * np.log10((1.0 + dd_zHEL) * dm_at(dd_zHD)) + 25.0 - dd_mu
        B = float(np.dot(dd_ones, d))
        c_sn = float(d @ (dd_icov @ d) - B * B / dd_sum)
    elif sn_set == 'PP':
        d = 5.0 * np.log10((1.0 + pp_zHEL) * dm_at(pp_zHD)) + 25.0 - pp_mb
        B = float(np.dot(pp_ones, d))
        c_sn = float(d @ (pp_icov @ d) - B * B / pp_sum)

    c_cmb = 0.0
    if use_cmb == 'BH':
        dm_st = float(dm_at(10.0)) + dm_highz(wbc, 10.0, 1090.0, u_hi_h2)
        rs_st = r_star_int(wb, wbc, 1090.0)
        dv = np.array([100.0 * np.sqrt(wbc + OMEGA_NU_MASSIVE_0) * dm_st / C_KMS, np.pi * dm_st / rs_st, wb]) - bh_mu
        c_cmb = float(dv @ bh_icov @ dv)
    elif use_cmb == 'DESI_CMB':
        dm_st = float(dm_at(10.0)) + dm_highz(wbc, 10.0, 1089.0, u_hi_h2)
        rs_st = r_star_int(wb, wbc, 1089.0)
        dv = np.array([rs_st / dm_st, wb, wbc]) - desi_cmb_mu
        c_cmb = float(dv @ desi_cmb_icov @ dv)
    return (c_bao + c_sn + c_cmb, c_bao, c_sn, c_cmb, info)


init_map = {
    'lcdm':        [0.02237, 0.1422, 0.682],
    'w0wa':        [0.02237, 0.1422, 0.674, -0.82, -0.66],
    'logistic':    [0.02237, 0.1422, 0.677, 0.84, 0.70],
    'schwinger_1': [0.02237, 0.1422, 0.57],
    'linear_1':    [0.02237, 0.1422, 0.61],
    'schwinger':   [0.02237, 0.1422, 0.571, 200.0],
    'linear':      [0.02237, 0.1422, 0.571, 2.47],
    'quadrature':  [0.02237, 0.1422, 0.510, 1.65],
    'lorentzian':  [0.02237, 0.1422, 0.500, 2.79],
    'gaussian':    [0.02237, 0.1422, 0.490, 3.35],
    'gen_lor':     [0.02237, 0.1422, 0.501, 2.63, 0.92],
}
dk_m = {
    'lcdm': 0, 'w0wa': 2, 'logistic': 2, 'schwinger_1': 0, 'linear_1': 0,
    'schwinger': 1, 'linear': 1, 'quadrature': 1, 'lorentzian': 1, 'gaussian': 1, 'gen_lor': 2,
}


def main():
    t0 = time.time()
    print(f'Sample counts: DESI BAO = {len(bao_z)}, DES-Dovekie = {len(dd_zHD)}, Pantheon+ (zHD > 0.01) = {len(pp_zHD)}', flush=True)
    print('=== 1. PRECISION CUBICSPLINE + BDF(1e-11) VERIFIED FIT TABLE ===', flush=True)
    dd_store = {}
    for tag, use_cmb, sn_set, N_obs in [
        ('CMB(BH)+DD+DESI',   'BH',       'DD', 1836),
        ('CMB(BH)+PP+DESI',   'BH',       'PP', 1606),
        ('CMB(BH)+DESI',      'BH',       None, 16),
        ('CMB(DESI)+DD+DESI', 'DESI_CMB', 'DD', 1836),
    ]:
        print(f'\n==================== {tag} (N={N_obs}) ====================', flush=True)
        r0 = minimize(lambda p: eval_exact('lcdm', p, use_cmb, sn_set, False)[0],
                      init_map['lcdm'], method='Nelder-Mead', options={'xatol': 1e-6, 'fatol': 1e-5})
        b0 = eval_exact('lcdm', r0.x, use_cmb, sn_set, True)
        chi0, bao0, sn0, cmb0, _ = b0
        store = {}
        for m in ['lcdm', 'w0wa', 'logistic', 'schwinger_1', 'linear_1', 'schwinger', 'linear',
                  'quadrature', 'lorentzian', 'gaussian', 'gen_lor']:
            if m in ('schwinger_1', 'linear_1'):
                bm = m.split('_')[0]
                res = minimize(lambda p3: eval_exact(bm, [p3[0], p3[1], p3[2], 1.0], use_cmb, sn_set, False)[0],
                               init_map[m], method='Nelder-Mead', options={'xatol': 1e-6, 'fatol': 1e-5, 'maxiter': 800})
                px = np.append(res.x, 1.0)
                bd_f = eval_exact(bm, px, use_cmb, sn_set, False)
                bd = eval_exact(bm, px, use_cmb, sn_set, True)
            else:
                res = minimize(lambda p: eval_exact(m, p, use_cmb, sn_set, False)[0],
                               init_map[m], method='Nelder-Mead', options={'xatol': 1e-6, 'fatol': 1e-5, 'maxiter': 1200})
                px = res.x
                bd_f = eval_exact(m, px, use_cmb, sn_set, False)
                bd = eval_exact(m, px, use_cmb, sn_set, True)
            store[m] = (px, bd)
            ctot, cbao, csn, ccmb, inf = bd
            dchi = ctot - chi0
            daic = dchi + 2.0 * dk_m[m]
            dbic = dchi + dk_m[m] * np.log(N_obs)
            diff_ode = abs(bd[0] - bd_f[0])
            bstr = f'dBAO={cbao - bao0:+6.2f}, dSN={csn - sn0:+6.2f}, dCMB={ccmb - cmb0:+5.2f}'
            pstr = f"H0={inf['H0']:.2f}, Om={inf['Om']:.4f}"
            if 'w0' in inf:
                zc_s = f"{inf['zc']:.3f}" if inf['zc'] is not None else 'none'
                pstr += (f", Hbar={inf['Hbar']:.2f}, c={px[3]:.3f}"
                         + (f", q={px[4]:.3f}" if m == 'gen_lor' else '')
                         + f", w0={inf['w0']:.3f}, w05={inf['w05']:.4f}, w1={inf['w1']:.3f}, zc={zc_s}")
            elif len(px) > 3:
                zc_s = f", zc={inf['zc']:.3f}" if inf.get('zc') is not None else ''
                pstr += f", de={np.round(px[3:], 3)}{zc_s}"
            print(f'  {m:11s} (dk={dk_m[m]}): chi2={ctot:8.2f} ({ctot:9.4f}) | dchi2={dchi:+7.2f} | '
                  f'dAIC={daic:+6.2f} | dBIC={dbic:+6.2f} | {bstr} | {pstr} [ode_chk={diff_ode:.1e}]', flush=True)
        if tag == 'CMB(BH)+DD+DESI':
            dd_store = store
        if tag in ('CMB(BH)+DD+DESI', 'CMB(BH)+PP+DESI'):
            for m in ('quadrature', 'lorentzian'):
                px_b, _ = store[m]
                c_b = px_b[3]
                chi_b = eval_exact(m, px_b, use_cmb, sn_set, False)[0]

                def prof_dchi(cv):
                    r = minimize(lambda p3: eval_exact(m, [p3[0], p3[1], p3[2], cv], use_cmb, sn_set, False)[0],
                                 px_b[:3], method='Nelder-Mead', options={'xatol': 1e-5, 'fatol': 1e-4, 'maxiter': 400})
                    return r.fun - chi_b - 1.0

                c_lo = brentq(prof_dchi, c_b * 0.65, c_b, xtol=1e-3)
                c_hi = brentq(prof_dchi, c_b, c_b * 1.50, xtol=1e-3)
                print(f'  >> 1-sigma CI ({m}): c = {c_b:.3f} (+{c_hi - c_b:.3f} / -{c_b - c_lo:.3f}) '
                      f'[range {c_lo:.3f}..{c_hi:.3f}]', flush=True)

    print('\n=== 2. NUMERICAL ROBUSTNESS AUDIT (START SHIFT, TOLERANCE, RADIATION, FULL z=1090 SPLICE) ===', flush=True)
    for m in ('lorentzian', 'quadrature'):
        px, _ = dd_store[m]
        wb, wbc, hb, cv = px
        _, z0, h0, u0, _ = solve_bdf(m, hb, wbc, cv, N_start=4.5, rtol=1e-11, atol=1e-13, rad_on=True)
        _, _, h3, _, _ = solve_bdf(m, hb, wbc, cv, N_start=3.5, rtol=1e-11, atol=1e-13, rad_on=True)
        _, _, h5, _, _ = solve_bdf(m, hb, wbc, cv, N_start=5.5, rtol=1e-11, atol=1e-13, rad_on=True)
        _, _, ht, _, _ = solve_bdf(m, hb, wbc, cv, N_start=4.5, rtol=1e-12, atol=1e-14, rad_on=True)
        _, _, hnr, _, _ = solve_bdf(m, hb, wbc, cv, N_start=4.5, rtol=1e-11, atol=1e-13, rad_on=False)
        N1k, z1k, h1k, u1k, _ = solve_bdf(m, hb, wbc, cv, N_start=4.5, rtol=1e-11, atol=1e-13, rad_on=True, z_max=1090.0)
        cs_dir = CubicSpline(N1k[::-1], np.exp(-N1k[::-1]) * C_KMS / (100.0 * hb * h1k[::-1]))
        dm_dir = float(cs_dir.integrate(-np.log(1091.0), 0.0))
        cs_sp = CubicSpline(z0, C_KMS / (100.0 * hb * h0))
        dm_sp = float(cs_sp.integrate(0.0, 10.0)) + dm_highz(wbc, 10.0, 1090.0, u0[-1] * hb * hb)
        print(f'  {m:10s} (c={cv:.3f}): N_start(3.5)={np.max(np.abs(h3 - h0) / h0):.2e}, '
              f'N_start(5.5)={np.max(np.abs(h5 - h0) / h0):.2e}, rtol(1e-12)={np.max(np.abs(ht - h0) / h0):.2e}, '
              f'rad_OFF={np.max(np.abs(hnr - h0) / h0):.2e}, dm(1090) splice_err={abs(dm_dir - dm_sp) / dm_dir:.2e}',
              flush=True)

    print('\n=== 3. DESI DR2 BIN-BY-BIN DIAGONAL PULLS (data - model)/sigma_i (CMB(BH)+DD+DESI BEST-FITS) ===', flush=True)
    print(f'{"z":>6s} {"obs":>11s} | {"LCDM":>7s} {"w0wa":>7s} {"logis":>7s} {"quad":>7s} {"loren":>7s} {"lin":>7s}', flush=True)
    for i in range(len(bao_z)):
        p = [dd_store[m][1][4]['bao_pulls_data_minus_model'][i]
             for m in ('lcdm', 'w0wa', 'logistic', 'quadrature', 'lorentzian', 'linear')]
        print(f'{bao_z[i]:6.3f} {bao_type[i]:>11s} | {p[0]:+7.2f} {p[1]:+7.2f} {p[2]:+7.2f} '
              f'{p[3]:+7.2f} {p[4]:+7.2f} {p[5]:+7.2f}', flush=True)

    print(f'\nTotal elapsed: {time.time() - t0:.1f}s', flush=True)


if __name__ == '__main__':
    main()
