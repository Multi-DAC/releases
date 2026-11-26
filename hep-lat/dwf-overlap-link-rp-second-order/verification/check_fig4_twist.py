"""Section 3.3, (H2): the persistent current at a twist outside {0, pi}. Clawd's check, written 2026-11-26.

d_y log det D_DW, with y a uniform phase on the spatial links of the first positive-time row (row index L of 2L),
at the free field with spatial twist theta_s, L = 6 (2L = 12 rows), Ls = 4, m = 0.7, L5 = 6, a5 = 0.25.
Uses build5 from check_letter77.py (5d Shamir operator, the twist on the last spatial link). Central difference.

    python -B check_fig4_twist.py
"""
import numpy as np

from check_letter77 import build5

L, Ls, m, L5, a5 = 6, 4, 0.7, 6, 0.25
nt = 2 * L
h = 1e-5


def logdet(y, twist, time_bc, mf):
    U0 = np.ones((nt, Ls), complex)
    U1 = np.ones((nt, Ls), complex)
    U1[L, :] *= np.exp(1j * y)
    return np.linalg.slogdet(build5(nt, Ls, L5, m, a5, mf, U1, U0, twist=twist, time_bc=time_bc))[1]


for time_bc in ("open", "ap"):
    for name, twist in (("0", 0.0), ("pi", np.pi), ("pi/3", np.pi / 3)):
        for kind, mf in (("DW", 0), ("PV", 1)):
            d = (logdet(h, twist, time_bc, mf) - logdet(-h, twist, time_bc, mf)) / (2 * h)
            print("time_bc %-4s twist %-4s %s: d_y log det = %+.4e" % (time_bc, name, kind, d))
