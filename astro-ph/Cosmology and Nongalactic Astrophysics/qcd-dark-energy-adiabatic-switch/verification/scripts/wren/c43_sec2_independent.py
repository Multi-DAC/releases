"""Independent Wren verification of Section 2 numerical results (standalone, needs no data files).

Reproduces Wren's independent checks referenced in the release README:
1. Lemma 2 early matter-era values u(-4) and x/h(-4) for the four switches of Section 2.7.
2. Section 2.6 radiation table (nine values of c plus c = 0.002 saddle-node check), using Wren's
   Radau u-evolution solver and exact CubicSpline + Brent root-finding for extrema.
3. Section 2.6 c > 3 radiation exclusion bound 1 + z(N_e) using Wren's radiation density
   (T_CMB = 2.7255 K, N_eff = 3.044, 3 relativistic neutrino species).
4. Section 2.4 / 2.8 fitted quadrature one-maximum check (c = 1.641, Hbar = 50.98, wbc = 0.14087)
   on 0 <= z <= 1090 with radiation.
5. Section 3.1 matter-only past-crossing thresholds c_min(Omega_m = 0.31) for quadrature and Lorentzian.
"""
import math
import time
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq

OMEGA_GAMMA_H2 = 2.47282e-5
F_NU_1 = (7.0 / 8.0) * (4.0 / 11.0) ** (4.0 / 3.0) * (3.044 / 3.0)
OMEGA_NU_1_REL = F_NU_1 * OMEGA_GAMMA_H2
OMEGA_RAD_MASSLESS = OMEGA_GAMMA_H2 * (1.0 + 2.0 * F_NU_1)
OMEGA_NU_MASSIVE_0 = 0.06 / 93.14


def omega_nu_m(a):
    p = 1.83
    return ((OMEGA_NU_MASSIVE_0 / a ** 3) ** p + (OMEGA_NU_1_REL / a ** 4) ** p) ** (1.0 / p)


def d_omega_nu_m_dN(a):
    p = 1.83
    t3 = (OMEGA_NU_MASSIVE_0 / a ** 3) ** p
    t4 = (OMEGA_NU_1_REL / a ** 4) ** p
    return -(3.0 * t3 + 4.0 * t4) * (t3 + t4) ** (1.0 / p - 1.0)


def check_lemma2():
    print('=== 1. LEMMA 2 MATTER-ERA CHECK (u(-4) vs 2 c_t / 3, x/h vs 1.5) ===', flush=True)
    for k, ct in [(0.2, 1.0), (0.2, 0.5), (0.2, 3.0), (0.4, 1.0)]:
        b, d = 1.0 + k, 1.0 / ct

        def inv_beta(y):
            # y = (1 + x) / (1 + b x + d x^2) => (d y) x^2 + (b y - 1) x + (y - 1) = 0
            if y >= 1.0:
                return 0.0
            A, B, C = d * y, b * y - 1.0, y - 1.0
            return (-B + math.sqrt(max(B * B - 4.0 * A * C, 0.0))) / (2.0 * A)

        N_start = 4.0
        mu_s = math.exp(-3.0 * N_start)
        u_s = 1.0 + mu_s * (1.0 - 3.0 * k) / (1.0 + 3.0 * k)

        def rhs(N, ua):
            u = max(ua[0], 1e-30)
            mu = math.exp(-3.0 * N)
            h = math.sqrt(mu + u)
            x = inv_beta(u / h)
            return [3.0 * mu - 2.0 * h * x]

        sol = solve_ivp(rhs, (N_start, -4.0), [u_s], method='Radau', rtol=1e-12, atol=1e-15)
        u_end = float(sol.y[0, -1])
        mu_end = math.exp(12.0)
        h_end = math.sqrt(mu_end + u_end)
        x_end = inv_beta(u_end / h_end)
        print(f'  k={k:.1f}, ct={ct:.1f}: u(-4)={u_end:.5f} (2ct/3={2 * ct / 3:.5f}), x/h(-4)={x_end / h_end:.5f}', flush=True)


def check_sec26_table_and_exclusion():
    print('\n=== 2. SECTION 2.6 RADIATION TABLE (wb=0.0224, wc=0.1196, Hbar=55) ===', flush=True)
    wbc = 0.0224 + 0.1196
    hbar = 0.55
    inv_hb2 = 1.0 / (hbar * hbar)
    mu_m0 = wbc * inv_hb2
    mu_r0 = (OMEGA_RAD_MASSLESS + OMEGA_NU_1_REL) * inv_hb2
    print(f'  mu_m0 = {mu_m0:.6f}, mu_r0 = {mu_r0:.6e} (all 3 neutrino species relativistic)', flush=True)

    def bg(a):
        mu = (wbc / a ** 3 + OMEGA_RAD_MASSLESS / a ** 4 + omega_nu_m(a)) * inv_hb2
        dmu = (-3.0 * wbc / a ** 3 - 4.0 * OMEGA_RAD_MASSLESS / a ** 4 + d_omega_nu_m_dN(a)) * inv_hb2
        return mu, dmu

    N_grid = np.linspace(0.0, -math.log(1091.0), 4001)
    for c in [0.002, 0.003, 0.01, 0.05, 0.5, 1.0, 2.469, 3.0, 4.0, 10.0]:
        N_start = 4.5 if c >= 0.5 else 6.0
        mu_m_s = (wbc + OMEGA_NU_MASSIVE_0) * math.exp(-3.0 * N_start) * inv_hb2
        mu_r_s = OMEGA_RAD_MASSLESS * math.exp(-4.0 * N_start) * inv_hb2
        u_s = 1.0 + mu_m_s * (c - 3.0) / (c + 3.0) + mu_r_s * (c - 4.0) / (c + 4.0)
        if abs(c - 3.0) < 1e-12:
            u_s += (1.0 / 12.0) * mu_m_s ** 2

        def rhs(N, ua):
            u = max(ua[0], 1e-30)
            mu, dmu = bg(math.exp(N))
            h = math.sqrt(mu + u)
            y = u / h
            x = 0.0 if y >= 1.0 else c * (1.0 / y - 1.0)
            return [-dmu - 2.0 * h * x]

        def jac(N, ua):
            u = max(ua[0], 1e-30)
            mu, _ = bg(math.exp(N))
            h = math.sqrt(mu + u)
            y = u / h
            if y >= 1.0:
                return [[0.0]]
            x = c * (1.0 / y - 1.0)
            dx_dy = -c / (y * y)
            dy_du = (2.0 * mu + u) / (2.0 * h * (mu + u))
            return [[-(x / h) - 2.0 * h * dx_dy * dy_du]]

        sol = solve_ivp(rhs, (N_start, N_grid[-1]), [u_s], method='Radau', jac=jac, t_eval=N_grid, rtol=1e-12, atol=1e-15)
        u_arr = sol.y[0]
        # Finite-difference du/dN on the smooth Radau trajectory (avoids amplifying 1e-12 stiff error by 2*h*c/y^2 ~ 1e8 at high z)
        du_arr = np.gradient(u_arr, N_grid)
        extrema = []
        for i in range(1, len(N_grid) - 1):
            if du_arr[i] * du_arr[i + 1] < 0 and max(abs(du_arr[i]), abs(du_arr[i + 1])) > 1e-9 * max(1.0, abs(u_arr[i])):
                lo, hi = max(0, i - 3), min(len(N_grid), i + 5)
                cs = CubicSpline(N_grid[lo:hi][::-1], du_arr[lo:hi][::-1])
                nr = brentq(lambda x: float(cs(x)), N_grid[i + 1], N_grid[i], xtol=1e-12)
                kind = 'MAX' if du_arr[i + 1] > 0 else 'MIN'
                a_r = math.exp(nr)
                mu_r, _ = bg(a_r)
                u_r = float(np.interp(nr, N_grid[::-1], u_arr[::-1]))
                extrema.append((kind, math.exp(-nr), u_r, u_r / (mu_r + u_r)))
        est = (3.0 * c * math.sqrt(mu_m0) / mu_r0) ** 0.4
        maxs = [e for e in extrema if e[0] == 'MAX']
        if maxs:
            z1, ur, om = maxs[0][1], maxs[0][2], maxs[0][3]
            print(f'  c={c:6.3f}: integrated 1+z_max = {z1:6.2f} | (2.9) = {est:6.2f} | Om_DE = {om:.2e} | extrema = {[(k, round(z, 2)) for k, z, _, _ in extrema]}', flush=True)
        else:
            print(f'  c={c:6.3f}: integrated 1+z_max =   none | (2.9) = {est:6.2f} | extrema = none', flush=True)

    print('\n=== 3. SECTION 2.6 EXCLUSION BOUND 1+z(N_e) FOR c > 3 (WREN RADIATION DENSITY) ===', flush=True)
    for c in [3.1, 3.5, 4.0, 6.0, 10.0, 20.0]:
        f = lambda u: (mu_r0 / (3.0 * math.sqrt(mu_m0))) * (u ** 2.5) + (4.0 * mu_r0 / (3.0 * mu_m0)) * u - (c - 3.0)
        zp1 = brentq(f, 1e-6, 1e6, xtol=1e-12)
        print(f'  c={c:4.1f}: 1+z(N_e) = {zp1:.4f}', flush=True)


def check_quadrature_one_max_and_cmin():
    print('\n=== 4. FITTED QUADRATURE ONE-MAXIMUM CHECK ON 0 <= z <= 1090 (WITH RADIATION) ===', flush=True)
    wbc = 0.022502 + 0.118365
    hbar = 0.509815
    c = 1.641046
    inv_hb2 = 1.0 / (hbar * hbar)
    N_start = 4.5
    N_grid = np.linspace(0.0, -math.log(1091.0), 4001)
    mu_s = (wbc / math.exp(3 * N_start) + OMEGA_RAD_MASSLESS / math.exp(4 * N_start) + omega_nu_m(math.exp(N_start))) * inv_hb2
    u_s = 1.0 + mu_s - (1.0 + 9.0 / (c * c)) * mu_s ** 2

    def rhs(N, ua):
        u = max(ua[0], 1e-30)
        a = math.exp(N)
        mu = (wbc / a ** 3 + OMEGA_RAD_MASSLESS / a ** 4 + omega_nu_m(a)) * inv_hb2
        dmu = (-3.0 * wbc / a ** 3 - 4.0 * OMEGA_RAD_MASSLESS / a ** 4 + d_omega_nu_m_dN(a)) * inv_hb2
        h = math.sqrt(mu + u)
        y = u / h
        x = 0.0 if y >= 1.0 else c * math.sqrt(max(1.0 / (y * y) - 1.0, 0.0))
        return [-dmu - 2.0 * h * x]

    sol = solve_ivp(rhs, (N_start, N_grid[-1]), [u_s], method='Radau', t_eval=N_grid, rtol=1e-12, atol=1e-15)
    u_arr = sol.y[0]
    du_arr = np.array([rhs(n, [u])[0] for n, u in zip(N_grid, u_arr)])
    exts = []
    for i in range(len(N_grid) - 1):
        if du_arr[i] * du_arr[i + 1] < 0 and max(abs(du_arr[i]), abs(du_arr[i + 1])) > 1e-8 * max(1.0, abs(u_arr[i])):
            lo, hi = max(0, i - 3), min(len(N_grid), i + 5)
            cs = CubicSpline(N_grid[lo:hi][::-1], du_arr[lo:hi][::-1])
            nr = brentq(lambda x: float(cs(x)), N_grid[i + 1], N_grid[i], xtol=1e-12)
            kind = 'MAX' if du_arr[i + 1] > 0 else 'MIN'
            exts.append((kind, math.exp(-nr) - 1.0, float(np.interp(nr, N_grid[::-1], u_arr[::-1]))))
    print(f'  quadrature c={c:.3f}: extrema on 0<=z<=1090 = {[(k, round(z, 3), round(u, 4)) for k, z, u in exts]}, u(z=1090)={u_arr[-1]:.4f}', flush=True)

    print('\n=== 5. MATTER-ONLY c_min(Omega_m = 0.31) THRESHOLDS ===', flush=True)

    def solve_u_mat(model, cv, mu0):
        Ns = max(2.0, math.log(mu0 / 1e-6) / 3.0)
        ms = mu0 * math.exp(-3.0 * Ns)
        us = 1.0 + ms - (1.0 + (9.0 if model == 'quadrature' else 18.0) / (cv * cv)) * ms ** 2

        def r(N, ua):
            u = ua[0]
            mu = mu0 * math.exp(-3.0 * N)
            h = math.sqrt(mu + u)
            y = u / h
            x = 0.0 if y >= 1.0 else (cv * math.sqrt(1.0 / (y * y) - 1.0) if model == 'quadrature' else cv * math.sqrt(1.0 / y - 1.0))
            return [3.0 * mu - 2.0 * h * x]

        s = solve_ivp(r, (Ns, 0.0), [us], method='Radau', rtol=1e-11, atol=1e-13)
        u0 = float(s.y[0, -1])
        du0 = r(0.0, [u0])[0]
        return mu0 / (mu0 + u0), -1.0 - du0 / (3.0 * u0)

    for mod, lo, hi in [('quadrature', 0.8, 1.5), ('lorentzian', 1.2, 2.2)]:
        cmin = brentq(lambda cv: solve_u_mat(mod, cv, brentq(lambda m: solve_u_mat(mod, cv, m)[0] - 0.31, 0.05, 5.0))[1] + 1.0, lo, hi)
        print(f'  {mod:10s}: c_min(Om=0.31) = {cmin:.5f}', flush=True)


if __name__ == '__main__':
    t0 = time.time()
    check_lemma2()
    check_sec26_table_and_exclusion()
    check_quadrature_one_max_and_cmin()
    print(f'\nTotal elapsed: {time.time() - t0:.1f}s', flush=True)
