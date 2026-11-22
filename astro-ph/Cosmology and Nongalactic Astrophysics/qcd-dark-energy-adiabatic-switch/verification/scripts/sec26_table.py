"""Section 2.6: the linear switch with radiation. For each c: every extremum of u = rho_DE/rho_dS on 0 <= z <= 1090,
the early maximum against the estimate (2.9), and the dark-energy share of the total density at the maximum.

Parameters as in the note: omega_b = 0.0224, omega_c = 0.1196, Hbar = 55 km/s/Mpc; neutrinos as in c43_fit.py in the
integration; in (2.9) mu_r0 counts photons and all three neutrino species as relativistic.

The integration (c43_fit.H_switch) solves the ODE up to a hand-over redshift and uses the quasi-static solution above
it. The default hand-over is z = 20, below six of the nine maxima in the table, so each case is run a second time with
the hand-over at z = 400, which puts every maximum inside the ODE part.

Usage: python -B sec26_table.py     (needs the data files only because c43_fit.py loads them on import; see DATA.md)
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import c43_fit as F

ob, oc, Hbar = 0.0224, 0.1196, 55.0
hb2 = (Hbar / 100) ** 2
mu_m0 = (ob + oc) / hb2
mu_r0 = (F.OG + F.onu_a4(1e-6)) / hb2        # a -> 0 value of omega_nu a^4: all species relativistic
print(f"mu_m0 {mu_m0:.5f}  mu_r0 {mu_r0:.4e}  (neutrinos relativistic)")


def extrema(u):
    """(index, kind) of every strict local extremum of u on the grid."""
    out = []
    for j in range(1, len(u) - 1):
        if u[j] > u[j - 1] and u[j] > u[j + 1]:
            out.append((j, 'max'))
        elif u[j] < u[j - 1] and u[j] < u[j + 1]:
            out.append((j, 'min'))
    return out


def refine(u, j):
    """1 + z of the extremum: parabola through three grid points in s = ln(1+z)."""
    a, b, _ = np.polyfit(F.SG[j - 1:j + 2], u[j - 1:j + 2], 2)
    return float(np.exp(-b / (2 * a)))


ratios = {}
shares = {}
for c in (10.0, 4.0, 3.1, 3.0, 2.469, 1.0, 0.5, 0.05, 0.01, 0.003, 0.002):
    pred = (3 * c * np.sqrt(mu_m0) / mu_r0) ** 0.4
    keep = {}
    for zsw in (20.0, 400.0):
        H, u = F.H_switch(ob, oc, Hbar, c, p=1, Ns=(3.0 if c >= 0.5 else 6.0), zsw=zsw)
        keep[zsw] = u
        ex = extrema(u)
        desc = [f"{kind} at 1+z = {refine(u, j):.2f} (u = {u[j]:.5f})" for j, kind in ex]
        line = f"c={c:<6} hand-over z={zsw:5.0f}: extrema of u on 0<=z<=1090: {desc if desc else 'none'}"
        mx = [j for j, kind in ex if kind == 'max']
        if mx:
            j = mx[-1]                      # the earliest maximum (largest z)
            zp1 = refine(u, j)
            h = H[j] / Hbar
            line += (f" | early maximum 1+z_max {zp1:.2f}  (2.9) gives {pred:.2f}  (2.9)/integrated {pred / zp1:.4f}"
                     f"  dark-energy share u/h^2 {u[j] / h ** 2:.2e}")
            if zsw == 400.0:
                ratios[c] = pred / zp1
                shares[c] = u[j] / h ** 2
        print(line)
    print(f"c={c:<6} largest relative difference in u between the two hand-overs, whole grid: "
          f"{np.max(np.abs(keep[400.0] / keep[20.0] - 1)):.1e}")

print("== summary, hand-over z = 400 ==")
inner = [c for c in ratios if 0.05 <= c <= 10]
print("largest |(2.9)/integrated - 1| for 0.05 <= c <= 10: %.4f (at c = %s)"
      % max((abs(ratios[c] - 1), c) for c in inner))
print("(2.9)/integrated - 1 at c = 0.003: %+.4f" % (ratios[0.003] - 1))
print("dark-energy share at the early maximum, over the nine values of c in the table of the note: %.2e to %.2e"
      % (min(shares[c] for c in shares if c != 3.1), max(shares[c] for c in shares if c != 3.1)))
