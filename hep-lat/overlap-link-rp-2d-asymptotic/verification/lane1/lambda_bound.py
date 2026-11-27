"""Lambda(L) of LANE1-Z-LEMMA (3.6), and 1 - tanh(a_Q/2). Deterministic."""
import numpy as np, mpmath as mp
m=0.7; eps0=0.5; ce=np.cos(eps0); mul=1.0
def Gex(L):
    k=2*np.pi*np.arange(L)/L; kh=(4*np.sin(k[:,None]/2)**2+4*np.sin(k[None,:]/2)**2).ravel()[1:]
    return (1/kh).sum()/(L*L)
def Gapprox(L): return np.log(L)/(2*np.pi)+0.048765
def Lam(L, el, Gf):
    V=L*L; mf=mul/L; C2=1/(1-mf**2); c1=np.pi*np.sqrt((2+m*m)*C2)
    hE=np.sqrt(6*el**2*Gf(L)/ce)/m
    K=16*(2+m*m)*C2; d=(hE**2/(K*V))**0.25
    tau=4*C2+4*C2*(2+m*m)+(9*np.pi/4)*C2*(2+m*m)*np.log(np.sqrt(2)*np.pi*(1+1/L)/(c1*d))
    lnm=np.log(1/mf)
    Lam=2*np.sqrt(K)*hE*L*lnm+4*lnm+2*np.log(2)*hE*np.sqrt(V*tau+d**-2)
    ok = (d<=0.25) and (d>=np.sqrt(2)*np.pi/(c1*L))
    return Lam/V, d, hE, ok
for el in (0.34,0.47):
    aQ=2*mp.pi**2*ce/el**2; print(f"el {el}: a_Q = {float(aQ):.3f}, 1 - tanh(a_Q/2) = {mp.nstr(1-mp.tanh(aQ/2),3)}")
print(" L      Lambda/V(el=.34)  Lambda/V(el=.47)  delta   h_E   valid")
for L in (64,128,256,1024,4096,10**4,10**5,10**6):
    Gf = Gex if L<=4096 else Gapprox
    a=Lam(L,0.34,Gf); b=Lam(L,0.47,Gf)
    print(f"{L:8d}  {a[0]:.4f}           {b[0]:.4f}           {b[1]:.4f}  {b[2]:.3f}  {b[3]}")
