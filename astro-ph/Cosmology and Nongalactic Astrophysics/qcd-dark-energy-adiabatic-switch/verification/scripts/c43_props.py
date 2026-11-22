"""Numerical tests of the C43 propositions (Clawd, 2026-11-21). Matter only, units Hbar = 1, mu = exp(-3N).
H-form:    h^2 = mu + h*beta(x)   (u = rho_DE/rho_dS = h*beta)
literal:   h^2 = mu + beta(x)     (u = beta)
x = -dh/dN.  Integrated backward from near de Sitter in the variable x.
"""
import sys, numpy as np, warnings
warnings.filterwarnings('ignore')
sys.stdout.reconfigure(line_buffering=True)
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

def traj(f, fp, k, literal=False, Ns=4.0, Ne=-4.0, n=16001, rtol=1e-11):
    mu = lambda N: np.exp(-3 * N)
    if literal:
        hof = lambda x, m: np.sqrt(m + f(x))
        rhs = lambda N, x: (-2 * hof(x, mu(N)) * x + 3 * mu(N)) / fp(x)
        a1 = 1 / (2 + 3 * k)
    else:
        hof = lambda x, m: 0.5 * (f(x) + np.sqrt(f(x) ** 2 + 4 * m))
        rhs = lambda N, x: (x * f(x) - 2 * hof(x, mu(N)) * x + 3 * mu(N)) / (hof(x, mu(N)) * fp(x))
        a1 = 1 / (1 + 3 * k)
    N = np.linspace(Ns, Ne, n)
    sol = solve_ivp(rhs, (Ns, Ne), [3 * a1 * mu(Ns)], method='LSODA', t_eval=N, rtol=rtol, atol=1e-16)
    assert sol.success and len(sol.y[0]) == n, sol.message
    x = sol.y[0]; m = mu(N); h = hof(x, m); b = f(x)
    u = b if literal else h * b
    return dict(N=N, x=x, mu=m, h=h, beta=b, u=u)

def extrema(N, u, lo=-3.5, hi=3.0):
    """local maxima / minima of u(N) (N descending in the arrays), by sign changes of du/dN; ignores |du| below noise."""
    du = np.gradient(u, N)
    sel = (N > lo) & (N < hi)
    idx = np.where(sel)[0]
    out = []
    for i in idx[:-1]:
        a, b = du[i], du[i + 1]          # N[i] > N[i+1]
        if a * b < 0 and max(abs(a), abs(b)) > 1e-9 * max(1.0, abs(u[i])):
            out.append(('max' if b > 0 else 'min', float(N[i]), i))   # going back in time: du/dN from - to + is a max in time
    return out

lin = lambda c: (lambda x: 1 / (1 + x / c), lambda x: -1 / c / (1 + x / c) ** 2, 1 / c)
lor = lambda c, q=1.0: (lambda x: (1 + (x / c) ** 2) ** -q, lambda x: -2 * q * x / c ** 2 * (1 + (x / c) ** 2) ** (-q - 1), 0.0)
def rat(k, ct, a=1.0):
    """beta = (1 + a x)/(1 + (a+k) x + (a/ct) x^2): beta'(0) = -k, x*beta -> ct."""
    b, d = a + k, a / ct
    f = lambda x: (1 + a * x) / (1 + b * x + d * x * x)
    fp = lambda x: (a * (1 + b * x + d * x * x) - (1 + a * x) * (b + 2 * d * x)) / (1 + b * x + d * x * x) ** 2
    return f, fp, k

print('== T1 sign rule: (u-1)/mu at N = 3 against (1-3k)/(1+3k) ==')
for k in (0.0, 0.2, 0.30, 0.32, 0.34, 0.40, 0.50, 1.0):
    f, fp, kk = rat(k, 1.0) if k > 0 else lor(2.8)
    t = traj(f, fp, kk)
    i = np.argmin(np.abs(t['N'] - 3.0))
    print(f'k={k:.2f}  measured {(t["u"][i]-1)/t["mu"][i]:+.7f}   predicted {(1-3*k)/(1+3*k):+.6f}')

print('== T1b threshold k = 1/3 (linear switch, c = 3): (u-1)/mu^2 at mu = 1e-4 against 1/12; three starting points ==')
Nstar = np.log(1e4) / 3
for Ns in (5.0, 6.0, 7.0):
    t = traj(*lin(3.0), Ns=Ns, Ne=Nstar, n=2)
    print(f'start N={Ns}: mu={t["mu"][-1]:.6e}  (u-1)/mu^2 = {(t["u"][-1]-1)/t["mu"][-1]**2:.6f}   1/12 = {1/12:.6f}')

print('== T2 early limit of u against 2 ct/3 (matter era), at N = -4 ==')
for k, ct in ((0.2, 1.0), (0.2, 0.5), (0.2, 3.0), (0.4, 1.0)):
    t = traj(*rat(k, ct))
    print(f'k={k} ct={ct}: u(-4)={t["u"][-1]:.5f}  2ct/3={2*ct/3:.5f}   x/h={t["x"][-1]/t["h"][-1]:.5f}')
t = traj(*lor(2.8)); print(f'Lorentzian c=2.8 (ct=0): u(-4)={t["u"][-1]:.3e}')

print('== T3 forced crossing: k < 1/3 and ct < 3/2 -> a maximum of u; mutants outside each hypothesis ==')
cases = [('inside', 0.2, 1.0), ('inside', 0.0, 1.0), ('inside', 0.3, 1.4), ('inside', 0.32, 1.49), ('inside', 0.1, 0.2),
         ('mutant k', 0.34, 1.0), ('mutant k', 0.5, 1.0), ('mutant ct (rational)', 0.2, 1.6), ('mutant ct (rational)', 0.2, 3.0)]
for name, k, ct in cases:
    if k == 0.0:
        d = 1.0 / ct; f = lambda x, d=d: (1 + x) / (1 + x + d * x * x)
        fp = lambda x, d=d: ((1 + x + d * x * x) - (1 + x) * (1 + 2 * d * x)) / (1 + x + d * x * x) ** 2
        t = traj(f, fp, 0.0)
    else:
        t = traj(*rat(k, ct))
    ex = extrema(t['N'], t['u'])
    print(f'{name:22s} k={k:.2f} ct={ct:.2f}: extrema {[(e[0], round(e[1],3), round(float(t["u"][e[2]]),4)) for e in ex]}  u(-4)={t["u"][-1]:.4f}')
for c in (5.0,):
    t = traj(*lin(c)); print(f'mutant ct (linear c={c}: k={1/c}, ct={c}): extrema {extrema(t["N"], t["u"])}  u(-4)={t["u"][-1]:.4f}')

print('== T4 linear form: no maximum of u for any c; beta monotone for c <= 3; at the minimum beta < c/3 ==')
for c in (0.3, 0.5, 1.0, 2.0, 2.9, 3.0, 3.1, 4.0, 6.0, 20.0):
    t = traj(*lin(c))
    ex = extrema(t['N'], t['u'])
    db = np.diff(t['beta'])            # N descending: beta should fall going back in time
    mono = bool(np.all(db[(t['N'][1:] > -3.5) & (t['N'][1:] < 3.0)] < 0))
    s = f'c={c:5.2f}: extrema {[(e[0], round(e[1],3)) for e in ex]}  beta increasing in time: {mono}'
    for e in ex:
        i = e[2]; s += f'  at {e[0]}: beta={t["beta"][i]:.4f} (c/3={c/3:.4f})'
    print(s)

print('== T5 literal form (u = beta): du/dN > 0 throughout for decreasing beta (w < -1) ==')
for name, fam in (('linear c=1', lin(1.0)), ('linear c=5', lin(5.0)), ('Lorentzian c=3', lor(3.0)), ('quadrature c=1.6', lor(1.6, 0.5)), ('rational k=.2 ct=1', rat(0.2, 1.0))):
    t = traj(*fam, literal=True)
    sel = (t['N'] > -3.5) & (t['N'] < 3.0)
    dx = np.diff(t['x'])[sel[1:]]
    print(f'{name:20s}: x decreasing in time everywhere: {bool(np.all(dx > 0))}; extrema of u: {extrema(t["N"], t["u"])}')

print('== T6 c_min: the c at which the maximum of u sits where Omega_m = mu/h^2 = 0.31 ==')
def om_at_max(fam):
    t = traj(*fam, n=64001)
    ex = [e for e in extrema(t['N'], t['u']) if e[0] == 'max']
    assert len(ex) == 1, ex
    i = ex[0][2]
    du = np.gradient(t['u'], t['N'])
    # linear interpolation of the zero of du between i and i+1
    w = du[i] / (du[i] - du[i + 1])
    om = t['mu'] / t['h'] ** 2
    return om[i] + w * (om[i + 1] - om[i])
for name, mk, lo, hi in (('quadrature', lambda c: lor(c, 0.5), 0.9, 1.3), ('Lorentzian', lambda c: lor(c, 1.0), 1.4, 1.8)):
    c0 = brentq(lambda c: om_at_max(mk(c)) - 0.31, lo, hi, xtol=1e-7)
    print(f'{name}: c_min = {c0:.5f}')
