"""PV gap check for PROOF-T1-2-finiteL5.md (Clawd, 2026-11-25).

Per fifth-dimension momentum k5, the Shamir Pauli-Villars operator in 2d (free field) is
    D(p0, k1, k5) = M + i sigma1 sin p0 + i sigma2 sin k1 + i (sin k5 / a5) sigma3,
    M = (1 - cos p0) + (1 - cos k1) - m + u5,   u5 = (1 - cos k5)/a5,
i.e. Wilson r = 1 in time and space at bulk mass -m, plus the 5th-direction hop
(1/a5)(1 - e^{ik5} P_- - e^{-ik5} P_+). det over spin = M^2 + sin^2 p0 + sin^2 k1 + sin^2 k5 / a5^2.

Checks (each with a control that must fail):
 (i)  Wren's closed form cosh E_PV = (c^2 + 1 + sin^2 k1 + 2 u5 (c + 1/a5)) / (2 |u5 + c|), c = 2 - m - cos k1,
      against a direct root of det D(p0 = iE) = 0. Control: the same with r = 2 in time (Wilson term 2(1 - cos p0)),
      which must NOT match.
 (ii) min over k5 of E_PV(k1, k5) >= E_cut(k1) := E_PV(k1, k5 = 0) for a5 <= a5_max(m), on a grid. Control: a5 large
      (a5 = 5), where the endpoint at u5 -> inf is low, must show a violation for some m.
 (iii) min over k1 of E_cut(k1) = min(-ln|1 - m|, ln(3 - m)) on a grid, m in (0, 1).
"""
import numpy as np
from scipy.optimize import brentq


def detD(E, k1, k5, m, a5, r_time=1.0):
    # p0 = iE: cos p0 = cosh E, sin p0 = i sinh E, sin^2 p0 = -sinh^2 E
    u5 = (1 - np.cos(k5)) / a5
    M = r_time * (1 - np.cosh(E)) + (1 - np.cos(k1)) - m + u5
    return M * M - np.sinh(E) ** 2 + np.sin(k1) ** 2 + (np.sin(k5) / a5) ** 2


def E_formula(k1, k5, m, a5):
    c = 2 - m - np.cos(k1)
    u5 = (1 - np.cos(k5)) / a5
    ch = (c * c + 1 + np.sin(k1) ** 2 + 2 * u5 * (c + 1 / a5)) / (2 * abs(u5 + c))
    return np.arccosh(ch)


def E_root(k1, k5, m, a5, r_time=1.0):
    # smallest positive E with detD = 0 (scan then bracket)
    Es = np.linspace(1e-9, 12, 24001)
    v = detD(Es, k1, k5, m, a5, r_time)
    idx = np.where(np.sign(v[:-1]) * np.sign(v[1:]) < 0)[0]
    if len(idx) == 0:
        return np.nan
    i = idx[0]
    return brentq(detD, Es[i], Es[i + 1], args=(k1, k5, m, a5, r_time), xtol=1e-14)


rng = np.random.default_rng(1125)
print("(i) closed form vs direct root, 400 random (k1, k5, m in (0,1), a5 in [0.05, 1])")
worst, worst_ctrl = 0.0, np.inf
for _ in range(400):
    k1, k5 = rng.uniform(-np.pi, np.pi, 2)
    m = rng.uniform(0.05, 0.95)
    a5 = rng.uniform(0.05, 1.0)
    ef, er = E_formula(k1, k5, m, a5), E_root(k1, k5, m, a5)
    worst = max(worst, abs(ef - er))
    ec = E_root(k1, k5, m, a5, r_time=2.0)
    if not np.isnan(ec):
        worst_ctrl = min(worst_ctrl, abs(ef - ec)) if abs(ef - ec) > 1e-6 else worst_ctrl
print(f"   max |formula - root| = {worst:.2e}   (r=2 control: closest non-matching gap {worst_ctrl:.2e})")
n_ctrl_match = 0
for _ in range(200):
    k1, k5 = rng.uniform(-np.pi, np.pi, 2)
    m = rng.uniform(0.05, 0.95)
    a5 = rng.uniform(0.05, 1.0)
    ec = E_root(k1, k5, m, a5, r_time=2.0)
    if not np.isnan(ec) and abs(E_formula(k1, k5, m, a5) - ec) < 1e-8:
        n_ctrl_match += 1
print(f"   r=2 control: {n_ctrl_match}/200 random points agree with the r=1 formula (must be ~0)")

print("(ii) min over k5 of E_PV >= E_cut(k1), grids k1, k5 of 721 points")
k1g = np.linspace(-np.pi, np.pi, 721)
k5g = np.linspace(-np.pi, np.pi, 721)
for m in [0.1, 0.3, 0.5, 0.7, 0.9]:
    for a5 in [1.0, 0.5, 0.1, 5.0]:
        K1, K5 = np.meshgrid(k1g, k5g, indexing="ij")
        Epv = E_formula(K1, K5, m, a5)
        Ecut = E_formula(k1g, 0.0 * k1g, m, a5)
        gap = (Epv - Ecut[:, None]).min()
        c = 2 - m - np.cos(k1g)
        cond = (c + 1 / a5 >= np.cosh(Ecut)).all()
        print(f"   m={m:.1f} a5={a5:<4} min(E_PV - E_cut) = {gap:+.3e}   endpoint condition c+1/a5 >= cosh E_cut: {cond}")

print("(iii) min over k1 of E_cut vs min(-ln(1-m), ln(3-m))")
for m in [0.1, 0.3, 0.5, 0.7, 0.7945, 0.9]:
    Ecut = E_formula(k1g, 0.0 * k1g, m, 1.0)
    pred = min(-np.log(1 - m), np.log(3 - m))
    print(f"   m={m:.4f}: grid min {Ecut.min():.10f} at k1={k1g[Ecut.argmin()]:+.4f}; predicted {pred:.10f}; rho_*(m) = exp(-2 Emin) = {np.exp(-2*pred):.6f}")
