import mpmath as mp
mp.mp.dps=40
m=mp.mpf('0.7'); cfT=mp.mpf('0.158939'); cfW=mp.mpf('0.3131624')
def mu0sq(e,m=m):
    a=1-m-e/4; return m**2-e-e**2/8-e**2/(2*a)
e2d=mp.findroot(lambda e: mu0sq(e), 0.29)
eN=m*(2-m)/(2+mp.sqrt(2))
print("eps2D*",e2d," epsNeu*",eN)
for nm,cf in [("H_T",cfT),("H_W",cfW)]:
    print(nm,"Neu ceiling",eN/mp.sqrt(2*cf)," 2D ceiling",e2d/mp.sqrt(2*cf), " eps=0.267:",mp.mpf('0.267')/mp.sqrt(2*cf))
# counterexample: mu0^2>0 but a<0
for mm,e in [(0.9,0.41),(0.95,0.25)]:
    mm=mp.mpf(mm); e=mp.mpf(e); print("m",mm,"eps",e,"a=",1-mm-e/4,"mu0^2=",mu0sq(e,mm))
# rev3 constants
eps=mp.mpf('0.267')
mu0=mp.sqrt(mu0sq(eps)); cPV=mu0/2
CH=16*(mp.e**0.5+mp.e); g0=min(mp.mpf(0.5), cPV**2/(2*CH))
muT=mu0/(6-m)
C3=64*mp.e**g0/cPV**2+1824*mp.e**(2*g0)/cPV**4
q=mp.e**(-2*g0); Cg=4*q*(1+q)/(1-q)**4; CPV=C3*Cg
K=4*(cfT+1); Kp=3*K
print("mu0",mu0,"g0",g0,"muT",muT,"C3",mp.nstr(C3,5),"Cgeom",mp.nstr(Cg,5),"CPV",mp.nstr(CPV,5))
for el in [mp.mpf('0.2'),mp.mpf('0.45'),mp.mpf('0.47'),mp.mpf('0.4727')]:
    C0=6*CPV*Kp*el**2/muT
    cstar=muT/(144*el**2*CPV)
    for CZ in [1,100]:
        def B(lnL):
            L=mp.e**lnL
            return mp.log(4)+C0*lnL-2*cstar*L**mp.mpf(1.25)/lnL+CZ*L*lnL**mp.mpf(1.5)
        lo,hi=mp.mpf(1),mp.mpf(3000)
        for _ in range(200):
            mid=(lo+hi)/2
            if B(mid*mp.log(10))>0: lo=mid
            else: hi=mid
        # theorem-valid: Piece B <= 1/2 exp(-c* ...)
        def T(lnL):
            L=mp.e**lnL
            return mp.log(8)+C0*lnL-cstar*L**mp.mpf(1.25)/lnL+CZ*L*lnL**mp.mpf(1.5)
        lo2,hi2=mp.mpf(1),mp.mpf(3000)
        for _ in range(200):
            mid=(lo2+hi2)/2
            if T(mid*mp.log(10))>0: lo2=mid
            else: hi2=mid
        print(f"el={el} CZ={CZ} C0={mp.nstr(C0,4)} c*={mp.nstr(cstar,4)} PieceB<1 at 10^{mp.nstr(hi,5)}; PieceB<=exp(-c*..)/2 at 10^{mp.nstr(hi2,5)}")
print("--- raw Lemma 3 exponent (4c*), with ln2|S_cut^+|, C_Z=1, C_W'=1")
for el in [mp.mpf('0.2'),mp.mpf('0.47')]:
    C0=6*CPV*Kp*el**2/muT; cstar=muT/(144*el**2*CPV)
    def R(lnL):
        L=mp.e**lnL; r0=mp.ceil(7*lnL/g0); S=2*L*(4*r0+10)
        return mp.log(4)+C0*lnL-4*cstar*L**mp.mpf(1.25)/lnL+mp.log(2)*S+L*lnL**mp.mpf(1.5)
    lo,hi=mp.mpf(1),mp.mpf(3000)
    for _ in range(200):
        mid=(lo+hi)/2
        if R(mid*mp.log(10))>0: lo=mid
        else: hi=mid
    print(f"el={el}: raw Piece B < 1 at 10^{mp.nstr(hi,5)}")
