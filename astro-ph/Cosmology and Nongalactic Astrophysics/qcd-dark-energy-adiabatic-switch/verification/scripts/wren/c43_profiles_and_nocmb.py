"""C43 q-profile CI, linear H-form wide c-scan, no-CMB Table 2 rows, c_min thresholds, and Sec 2.6 checks (Wren, 2026-11-21).

Imports c43_precision_verify from this folder (reads C43_DATA). Reproduces:
1. 1-sigma profile CI on q for the generalized Lorentzian on CMB(BH)+DD+DESI.
2. Linear H-form wide profile scan over c in [0.05, 1000] on CMB(BH)+DD+DESI.
3. No-CMB rows of Table 2 (DD+DESI and PP+DESI, including the external calibration against Lee et al.).
4. Matter-only past-crossing thresholds c_min(Omega_m = 0.31) for quadrature and Lorentzian.
5. Independent Wren checks of Section 2.6 (radiation exclusion bound 1+z(N_e) and 9-row table) and Lemma 2.
"""
import math
import os
import sys
import time
import numpy as np
from scipy.integrate import solve_ivp, cumulative_trapezoid
from scipy.interpolate import CubicSpline
from scipy.optimize import minimize, brentq

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import c43_precision_verify as P  # noqa: E402


def run_q_profile_and_linear_scan():
    print('=== 1. GENERALIZED LORENTZIAN q PROFILE CI & LINEAR H-FORM WIDE c-SCAN (CMB(BH)+DD+DESI) ===', flush=True)
    r_lcdm = minimize(lambda p: P.eval_exact('lcdm', p, 'BH', 'DD', False)[0],
                      P.init_map['lcdm'], method='Nelder-Mead', options={'xatol': 1e-6, 'fatol': 1e-5})
    chi_lcdm = P.eval_exact('lcdm', r_lcdm.x, 'BH', 'DD', True)[0]

    r_gl = minimize(lambda p: P.eval_exact('gen_lor', p, 'BH', 'DD', False)[0],
                    P.init_map['gen_lor'], method='Nelder-Mead', options={'xatol': 1e-6, 'fatol': 1e-5, 'maxiter': 1200})
    px_gl = r_gl.x
    chi_gl_best = P.eval_exact('gen_lor', px_gl, 'BH', 'DD', True)[0]
    print(f'Generalized Lorentzian best fit: chi2={chi_gl_best:.4f}, c={px_gl[3]:.4f}, q={px_gl[4]:.4f}', flush=True)

    q_grid = [0.45, 0.50, 0.60, 0.70, 0.80, 0.915, 1.00, 1.25, 1.50, 2.00, 3.00, 5.00]
    p_cur = px_gl[:4].copy()
    for qv in q_grid:
        rq = minimize(lambda p4: P.eval_exact('gen_lor', [p4[0], p4[1], p4[2], p4[3], qv], 'BH', 'DD', False)[0],
                      p_cur, method='Nelder-Mead', options={'xatol': 1e-5, 'fatol': 1e-4, 'maxiter': 600})
        bdq = P.eval_exact('gen_lor', [rq.x[0], rq.x[1], rq.x[2], rq.x[3], qv], 'BH', 'DD', True)
        p_cur = rq.x
        print(f'  q={qv:5.3f}: chi2={bdq[0]:.4f}, dchi2_from_best={bdq[0] - chi_gl_best:+.3f}, c={rq.x[3]:.3f}', flush=True)

    def prof_q(qv):
        c_guess = px_gl[3] * math.sqrt(qv / px_gl[4])
        rq = minimize(lambda p4: P.eval_exact('gen_lor', [p4[0], p4[1], p4[2], p4[3], qv], 'BH', 'DD', False)[0],
                      [px_gl[0], px_gl[1], px_gl[2], c_guess], method='Nelder-Mead',
                      options={'xatol': 1e-5, 'fatol': 1e-4, 'maxiter': 600})
        return P.eval_exact('gen_lor', [rq.x[0], rq.x[1], rq.x[2], rq.x[3], qv], 'BH', 'DD', True)[0] - chi_gl_best - 1.0

    q_lo = brentq(prof_q, 0.55, px_gl[4], xtol=2e-3)
    q_hi = brentq(prof_q, px_gl[4], 2.50, xtol=2e-3)
    print(f'  >> 1-sigma profile CI on q: q = {px_gl[4]:.3f} (+{q_hi - px_gl[4]:.3f} / -{px_gl[4] - q_lo:.3f}) '
          f'[range {q_lo:.3f} .. {q_hi:.3f}]', flush=True)

    r_lin = minimize(lambda p: P.eval_exact('linear', p, 'BH', 'DD', False)[0],
                     P.init_map['linear'], method='Nelder-Mead', options={'xatol': 1e-6, 'fatol': 1e-5, 'maxiter': 1200})
    p_lin = r_lin.x[:3].copy()
    print('\nLinear H-form wide profile scan (CMB(BH)+DD+DESI):', flush=True)
    for cv in [0.05, 0.1, 0.2, 0.5, 1.0, 1.5, 2.0, 2.469, 3.0, 5.0, 10.0, 30.0, 100.0, 1000.0]:
        rl = minimize(lambda p3: P.eval_exact('linear', [p3[0], p3[1], p3[2], cv], 'BH', 'DD', False)[0],
                      p_lin, method='Nelder-Mead', options={'xatol': 1e-5, 'fatol': 1e-4, 'maxiter': 600})
        bdl = P.eval_exact('linear', [rl.x[0], rl.x[1], rl.x[2], cv], 'BH', 'DD', True)
        print(f'  c={cv:7.3f}: dchi2={bdl[0] - chi_lcdm:+8.2f}, H0={bdl[4]["H0"]:.2f}, Om={bdl[4]["Om"]:.4f}', flush=True)


def run_nocmb_table2():
    print('\n=== 2. NO-CMB AUXILIARY FITS (TABLE 2: DD+DESI AND PP+DESI) ===', flush=True)
    N_FUT_F = np.linspace(4.0, 0.0, 121)
    N_PAST_F = np.linspace(0.0, -np.log(11.0), 161)
    N_ALL_F = np.concatenate([N_FUT_F[:-1], N_PAST_F])
    A_ALL_F = np.exp(N_ALL_F)
    Z_PAST_F = np.exp(-N_PAST_F) - 1.0
    A_PAST_F = np.exp(N_PAST_F)
    RAD_F = P.OMEGA_RAD_MASSLESS / A_ALL_F ** 4 + P.omega_nu_m(A_ALL_F)
    DRAD_F = -4.0 * P.OMEGA_RAD_MASSLESS / A_ALL_F ** 4 + P.d_omega_nu_m_dN(A_ALL_F)
    INV_A3_F = 1.0 / A_ALL_F ** 3

    def solve_cn(model, hbar, wbc, c, q=1.0):
        if hbar <= 0.15 or hbar >= 1.5 or wbc <= 0.05 or wbc >= 0.35 or c <= 0.02 or c >= 5000.0:
            return None
        inv_hb2 = 1.0 / (hbar * hbar)
        mu_all = (wbc * INV_A3_F + RAD_F) * inv_hb2
        dmu_all = (-3.0 * wbc * INV_A3_F + DRAD_F) * inv_hb2
        mu_s = mu_all[0]
        if model == 'linear':
            u = 1.0 + ((c - 3.0) / (c + 3.0)) * mu_s
        elif model == 'quadrature':
            u = 1.0 + mu_s - (1.0 + 9.0 / (c * c)) * mu_s ** 2
        elif model == 'lorentzian':
            u = 1.0 + mu_s - (1.0 + 18.0 / (c * c)) * mu_s ** 2
        else:
            u = 1.0 + mu_s - (1.0 + 18.0 * q / (c * c)) * mu_s ** 2
        u_arr = np.empty(len(N_ALL_F))
        u_arr[0] = u
        f_prev, _ = P.rhs_and_dfdu(model, u, mu_all[0], dmu_all[0], c, q)
        for i in range(1, len(N_ALL_F)):
            hdN = 0.5 * (N_ALL_F[i] - N_ALL_F[i - 1])
            if u <= 1e-25:
                u = 0.0
                u_arr[i] = 0.0
                f_prev = 0.0
                continue
            mu_i = mu_all[i]
            dmu_i = dmu_all[i]
            un = u
            for _ in range(8):
                fi, dfi = P.rhs_and_dfdu(model, un, mu_i, dmu_i, c, q)
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
            f_prev, _ = P.rhs_and_dfdu(model, u, mu_i, dmu_i, c, q) if u > 0 else (0.0, 0.0)
        up = u_arr[120:]
        if up[0] <= 1e-6:
            return None
        hp = np.sqrt(mu_all[120:] + up)
        return N_PAST_F, Z_PAST_F, hp, up, 100.0 * hbar * hp[0]

    def solve_bdf_f(model, hbar, wbc, c, q=1.0):
        inv_hb2 = 1.0 / (hbar * hbar)
        N_start = 4.0
        mu_s = P.omega_bg_h2(np.exp(N_start), wbc) * inv_hb2
        if model == 'linear':
            u_s = 1.0 + ((c - 3.0) / (c + 3.0)) * mu_s
        elif model == 'quadrature':
            u_s = 1.0 + mu_s - (1.0 + 9.0 / (c * c)) * mu_s ** 2
        elif model == 'lorentzian':
            u_s = 1.0 + mu_s - (1.0 + 18.0 / (c * c)) * mu_s ** 2
        else:
            u_s = 1.0 + mu_s - (1.0 + 18.0 * q / (c * c)) * mu_s ** 2
        rhs = lambda N, ua: [P.rhs_and_dfdu(model, max(ua[0], 0.0), P.omega_bg_h2(np.exp(N), wbc) * inv_hb2,
                                            P.d_omega_bg_h2_dN(np.exp(N), wbc) * inv_hb2, c, q)[0]]
        jac = lambda N, ua: [[P.rhs_and_dfdu(model, max(ua[0], 0.0), P.omega_bg_h2(np.exp(N), wbc) * inv_hb2,
                                             P.d_omega_bg_h2_dN(np.exp(N), wbc) * inv_hb2, c, q)[1]]]
        sol = solve_ivp(rhs, (N_start, N_PAST_F[-1]), [u_s], method='BDF', jac=jac, t_eval=N_PAST_F, rtol=1e-9, atol=1e-12)
        if (not sol.success) or len(sol.y[0]) != len(N_PAST_F):
            return solve_cn(model, hbar, wbc, c, q)
        u_ev = np.maximum(sol.y[0], 0.0)
        h_ev = np.sqrt(P.omega_bg_h2(A_PAST_F, wbc) * inv_hb2 + u_ev)
        return N_PAST_F, Z_PAST_F, h_ev, u_ev, 100.0 * hbar * h_ev[0]

    def eval_nocmb(model, params, sn_set, precise=False):
        if model in ('lcdm', 'w0wa', 'logistic'):
            h, om_m, rd = params[0], params[1], params[2]
            if h < 0.4 or h > 1.0 or om_m < 0.05 or om_m > 0.7 or rd < 100 or rd > 200:
                return 1e12
            wbc = om_m * h * h - P.OMEGA_NU_MASSIVE_0
            de_p = params[3:]
        else:
            hbar, wbc, rd = params[0], params[1], params[2]
            if rd < 100 or rd > 200:
                return 1e12
            de_p = params[3:]
        if model == 'lcdm':
            om_de = h * h - P.omega_bg_h2(1.0, wbc)
            if om_de <= 0.05:
                return 1e12
            H_g = 100.0 * np.sqrt(P.omega_bg_h2(A_PAST_F, wbc) + om_de)
        elif model == 'w0wa':
            w0, wa = de_p[0], de_p[1]
            if w0 < -2.5 or w0 > 0.5 or wa < -4.0 or wa > 2.0:
                return 1e12
            om_de = h * h - P.omega_bg_h2(1.0, wbc)
            if om_de <= 0.05:
                return 1e12
            E2 = P.omega_bg_h2(A_PAST_F, wbc) + om_de * A_PAST_F ** (-3.0 * (1.0 + w0 + wa)) * np.exp(-3.0 * wa * (1.0 - A_PAST_F))
            if np.any(E2 <= 0):
                return 1e12
            H_g = 100.0 * np.sqrt(E2)
        elif model == 'logistic':
            zq, dzq = de_p[0], de_p[1]
            if zq < -0.5 or zq > 3.5 or dzq < 0.05 or dzq > 3.0:
                return 1e12
            om_de = h * h - P.omega_bg_h2(1.0, wbc)
            if om_de <= 0.05:
                return 1e12
            Aq = 1.0 + np.exp((-1.0 - zq) / dzq)
            b0 = Aq / (1.0 + np.exp(-zq / dzq))
            hb_v = om_de / (h * b0)
            hb = hb_v * (Aq / (1.0 + np.exp((Z_PAST_F - zq) / dzq)))
            H_g = 100.0 * 0.5 * (hb + np.sqrt(hb * hb + 4.0 * P.omega_bg_h2(A_PAST_F, wbc)))
        else:
            cv = de_p[0]
            qv = de_p[1] if model == 'gen_lor' else 1.0
            if model == 'gen_lor' and (qv < 0.30 or qv > 5.0):
                return 1e12
            out = solve_bdf_f(model, hbar, wbc, cv, qv) if precise else solve_cn(model, hbar, wbc, cv, qv)
            if out is None or out[4] < 50 or out[4] > 85:
                return 1e12
            H_g = 100.0 * hbar * out[2]
        dm_g = np.concatenate([[0.0], cumulative_trapezoid(P.C_KMS / H_g, Z_PAST_F)])
        pred = np.empty(len(P.bao_z))
        for i, (zb, tp) in enumerate(zip(P.bao_z, P.bao_type)):
            dm = np.interp(zb, Z_PAST_F, dm_g)
            dh = P.C_KMS / np.interp(zb, Z_PAST_F, H_g)
            if tp == 'DM_over_rs':
                pred[i] = dm / rd
            elif tp == 'DH_over_rs':
                pred[i] = dh / rd
            else:
                pred[i] = (zb * dm * dm * dh) ** (1.0 / 3.0) / rd
        db = pred - P.bao_val
        c2 = float(db @ P.bao_icov @ db)
        if sn_set == 'DD':
            d = 5.0 * np.log10((1.0 + P.dd_zHEL) * np.interp(P.dd_zHD, Z_PAST_F, dm_g)) + 25.0 - P.dd_mu
            B = float(np.dot(P.dd_ones, d))
            c2 += float(d @ (P.dd_icov @ d) - B * B / P.dd_sum)
        else:
            d = 5.0 * np.log10((1.0 + P.pp_zHEL) * np.interp(P.pp_zHD, Z_PAST_F, dm_g)) + 25.0 - P.pp_mb
            B = float(np.dot(P.pp_ones, d))
            c2 += float(d @ (P.pp_icov @ d) - B * B / P.pp_sum)
        return c2

    c_grids = {
        'linear':     [0.5, 1.5, 3.0, 6.0, 15.0, 50.0, 500.0],
        'quadrature': [1.2, 1.5, 1.8, 2.2, 3.0, 5.0, 15.0],
        'lorentzian': [1.8, 2.5, 3.2, 4.2, 6.0, 12.0, 40.0],
    }
    std_x0 = {
        'lcdm':     [0.68, 0.305, 148.0],
        'w0wa':     [0.68, 0.312, 147.0, -0.85, -0.50],
        'logistic': [0.68, 0.323, 147.0, 0.80, 0.53],
        'gen_lor':  [0.50, 0.142, 147.5, 2.80, 1.00],
    }
    for tag, sn_set in [('DD+DESI (no CMB, N=1833)', 'DD'), ('PP+DESI (no CMB, N=1603)', 'PP')]:
        print(f'\n-------------------- {tag} --------------------', flush=True)
        r0 = minimize(lambda p: eval_nocmb('lcdm', p, sn_set), std_x0['lcdm'], method='Powell', options={'xtol': 1e-4, 'ftol': 1e-4, 'maxiter': 400})
        r0 = minimize(lambda p: eval_nocmb('lcdm', p, sn_set), r0.x, method='Nelder-Mead', options={'xatol': 1e-5, 'fatol': 1e-4, 'maxiter': 600})
        chi0 = eval_nocmb('lcdm', r0.x, sn_set, precise=True)
        print(f'  lcdm       : chi2={chi0:8.2f} | dchi2= +0.00 | Om={r0.x[1]:.4f}', flush=True)
        for m in ['w0wa', 'logistic', 'linear', 'quadrature', 'lorentzian', 'gen_lor']:
            if m in c_grids:
                best_chi, best_p = 1e12, None
                p3_init = [0.54, 0.142, 147.5]
                for cv in c_grids[m]:
                    rr = minimize(lambda p3: eval_nocmb(m, [p3[0], p3[1], p3[2], cv], sn_set),
                                  p3_init, method='Powell', options={'xtol': 1e-4, 'ftol': 1e-4, 'maxiter': 300})
                    if rr.fun < 1e10:
                        p3_init = rr.x
                    if rr.fun < best_chi:
                        best_chi = rr.fun
                        best_p = np.append(rr.x, cv)
                rf = minimize(lambda p: eval_nocmb(m, p, sn_set), best_p, method='Nelder-Mead',
                              options={'xatol': 1e-5, 'fatol': 1e-4, 'maxiter': 500})
                px = rf.x if rf.fun < best_chi else best_p
            else:
                r1 = minimize(lambda p: eval_nocmb(m, p, sn_set), std_x0[m], method='Powell', options={'xtol': 1e-4, 'ftol': 1e-4, 'maxiter': 400})
                r2 = minimize(lambda p: eval_nocmb(m, p, sn_set), r1.x, method='Nelder-Mead', options={'xatol': 1e-5, 'fatol': 1e-4, 'maxiter': 600})
                px = r2.x
            chi = eval_nocmb(m, px, sn_set, precise=True)
            extra = ''
            if m in ('linear', 'quadrature', 'lorentzian', 'gen_lor'):
                qv = px[4] if m == 'gen_lor' else 1.0
                _, _, _, _, H0_c = solve_bdf_f(m, px[0], px[1], px[3], qv)
                extra = f' | H0={H0_c:.2f}, Hbar={px[0]*100:.2f}, c={px[3]:.3f}' + (f', q={qv:.2f}' if m == 'gen_lor' else '')
            print(f'  {m:11s}: chi2={chi:8.2f} | dchi2={chi - chi0:+6.2f}{extra}', flush=True)


def run_cmin_and_sec26_checks():
    print('\n=== 3. MATTER-ONLY c_min(Omega_m = 0.31) THRESHOLDS ===', flush=True)

    def solve_u_mat(model, c, mu0):
        N_start = max(2.0, np.log(mu0 / 1e-6) / 3.0)
        mu_s = mu0 * np.exp(-3.0 * N_start)
        u_s = 1.0 + mu_s - (1.0 + (9.0 if model == 'quadrature' else 18.0) / (c * c)) * mu_s ** 2

        def rhs(N, ua):
            u = ua[0]
            mu = mu0 * np.exp(-3.0 * N)
            h = np.sqrt(mu + u)
            y = u / h
            if y >= 1.0:
                x = 0.0
            elif model == 'quadrature':
                x = c * np.sqrt(1.0 / (y * y) - 1.0)
            else:
                x = c * np.sqrt(1.0 / y - 1.0)
            return [3.0 * mu - 2.0 * h * x]

        sol = solve_ivp(rhs, (N_start, 0.0), [u_s], method='Radau', rtol=1e-11, atol=1e-13)
        u0 = sol.y[0, -1]
        du0 = rhs(0.0, [u0])[0]
        return mu0 / (mu0 + u0), -1.0 - du0 / (3.0 * u0)

    for mod, lo, hi in [('quadrature', 0.8, 1.5), ('lorentzian', 1.2, 2.2)]:
        cmin = brentq(lambda c: solve_u_mat(mod, c, brentq(lambda m: solve_u_mat(mod, c, m)[0] - 0.31, 0.05, 5.0))[1] + 1.0, lo, hi)
        print(f'  {mod:10s}: c_min(Om=0.31) = {cmin:.5f}', flush=True)

    print('\n=== 4. INDEPENDENT CHECK OF SEC 2.6 EXCLUSION BOUND 1+z(N_e) (WREN RADIATION DENSITY) ===', flush=True)
    hb = 0.55
    mu_m0 = (0.0224 + 0.1196) / (hb ** 2)
    mu_r0 = (P.OMEGA_RAD_MASSLESS + P.OMEGA_NU_1_REL) / (hb ** 2)
    print(f'  mu_m0 = {mu_m0:.6f}, mu_r0 = {mu_r0:.6e} (T_CMB=2.7255 K, N_eff=3.044, 3 relativistic nu)', flush=True)
    for c in [3.1, 3.5, 4.0, 6.0, 10.0, 20.0]:
        f = lambda u: (mu_r0 / (3.0 * math.sqrt(mu_m0))) * (u ** 2.5) + (4.0 * mu_r0 / (3.0 * mu_m0)) * u - (c - 3.0)
        zp1 = brentq(f, 1e-6, 1e6, xtol=1e-12)
        print(f'  c={c:4.1f}: 1+z(N_e) = {zp1:.4f}', flush=True)


if __name__ == '__main__':
    t0 = time.time()
    run_q_profile_and_linear_scan()
    run_nocmb_table2()
    run_cmin_and_sec26_checks()
    print(f'\nTotal elapsed: {time.time() - t0:.1f}s', flush=True)
