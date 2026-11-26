"""Runs test T4 of letter78_tests.py alone: on the annulus r_mid <= |z| <= 1/r_mid, the minimum of Re(w^2)/c1(m)
and the maximum of |b|, at U = 1 (Proposition 5, (a), and the numerical check at its end).

Written for this folder on 2026-11-25. letter78_tests.py as a whole (T1-T4) needs about 8 GB (dense residues of
Prop at L5 = 320 in T2 and T3) and ran about 17 minutes, so its full output is recorded, not rerun; T4 needs
neither, and is rerun here through this wrapper. letter78_tests.py has a main guard, so importing it runs nothing.
"""
from letter78_tests import T4

T4()
print("DONE")
