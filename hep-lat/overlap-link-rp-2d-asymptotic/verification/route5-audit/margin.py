"""Attack (3): compute, from the proof's own stated constants, the L at which each piece of the
Route (v) bound first becomes < 1. Unstated constants (C_3 of Lemma 2a, the path-count C in C_geom,
C_Z of Lemma Z, C_W') are set to 1, which can only move the crossover DOWN if the true values are larger."""
import mpmath as mp
mp.mp.dps = 50
m, eps, cf = mp.mpf('0.7'), mp.mpf('0.267'), mp.mpf('0.158939')
a = 1 - m - eps/4
mu0 = mp.sqrt(m**2 - eps - eps**2/8 - eps**2/(2*a))
cPV = mu0/2
CH = 128*mp.e**2
g0 = min(mp.mpf(1), cPV**2/(2*CH))
muT = mu0/(6 - m)
K = 4*(cf + 1); Kp = 3*K
print(f"mu0={mp.nstr(mu0,6)} cPV^2={mp.nstr(cPV**2,6)} gamma0={mp.nstr(g0,6)} muT={mp.nstr(muT,6)}")
def Cgeom(C=1):
    return C*mp.nsum(lambda d: d**3*mp.e**(-2*g0*d), [1, mp.inf])
for C3, Cg, CZ in [(1, 1, 1)]:
    CPV = C3*Cgeom(Cg)
    print(f"C_geom(C=1) = {mp.nstr(Cgeom(1),6)}  C_PV = {mp.nstr(CPV,6)}")
    for el in [mp.mpf('0.45'), mp.mpf('0.2')]:
        C0 = 5*CPV*Kp*el**2/muT
        def pieceB(lnL):
            L = mp.e**lnL
            r0 = mp.ceil(7*lnL/g0)
            S = 2*L*(4*r0 + 10)
            tail = muT*L**mp.mpf(1.25)/(30*el**2*CPV*lnL)
            return mp.log(4) + C0*lnL + mp.log(2)*S - tail + CZ*L*lnL**mp.mpf(1.5)
        # find crossover by bisection on log10 L
        lo, hi = mp.mpf(1), mp.mpf(2000)
        for _ in range(200):
            mid = (lo+hi)/2
            if pieceB(mid*mp.log(10)) > 0: lo = mid
            else: hi = mid
        lnL = hi*mp.log(10)
        r0 = mp.ceil(7*lnL/g0)
        print(f" e*l={el}: C_0 = {mp.nstr(C0,4)};  Piece B bound < 1 first at L ~ 10^{mp.nstr(hi,5)}"
              f"  (there r_0 = {mp.nstr(r0,4)}, L5 = {mp.nstr(mp.ceil(5*lnL/muT),4)})")
        # domain condition L > 4 r0 + 12
        lo, hi = mp.mpf(1), mp.mpf(30)
        for _ in range(200):
            mid = (lo+hi)/2
            lnL = mid*mp.log(10)
            if mp.e**lnL <= 4*mp.ceil(7*lnL/g0) + 12: lo = mid
            else: hi = mid
        print(f"   strip-definability L > 4 r_0 + 12 first at L ~ 10^{mp.nstr(hi,4)}")
        # where the proof's own margin claim L^{5/4}/lnL >> L (lnL)^{3/2} holds with ratio 2, pure functions, unit constants
        lo, hi = mp.mpf(1), mp.mpf(200)
        for _ in range(200):
            mid = (lo+hi)/2
            lnL = mid*mp.log(10); L = mp.e**lnL
            if L**mp.mpf(1.25)/lnL < 2*L*lnL**mp.mpf(1.5): lo = mid
            else: hi = mid
        print(f"   bare margin L^(5/4)/lnL > 2 L (lnL)^(3/2) (unit constants) first at L ~ 10^{mp.nstr(hi,4)}")
