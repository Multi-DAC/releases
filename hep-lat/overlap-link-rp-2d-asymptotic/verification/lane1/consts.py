"""Deterministic constants for LANE1 Step 1 (no sampling). 2026-11-26."""
import numpy as np
from scipy.special import ive, erf
from scipy.integrate import quad

def khat2(L):
    k = 2*np.pi*np.arange(L)/L
    return (4*np.sin(k[:,None]/2)**2 + 4*np.sin(k[None,:]/2)**2)

def G(L):
    kh = khat2(L).ravel()[1:]
    return (1/kh).sum()/(L*L)

def S4(L):
    kh = khat2(L).ravel()[1:]
    return (1/kh**2).sum()

print("L   G(L)      bound (1+ln(L/2))/2   G - ln(L)/(2pi)   relsd(||a||^2) Gaussian")
for L in (8,12,16,24,32,48,64,128,256):
    g = G(L); b = (1+np.log(L/2))/2
    rel = np.sqrt(2*S4(L))/(L*L*g)   # sd/mean for ||a||^2 = sum |th_k|^2/khat^2, th iid N(0,1/beta)
    print(f"{L:4d} {g:.6f}  {b:.6f}  {'OK' if g<=b else 'FAIL'}   {g-np.log(L)/(2*np.pi):+.6f}     {rel:.4f}")

# e^{-x} I0(x) >= erf(pi sqrt(x/2))/sqrt(2 pi x)
xs = np.concatenate([np.linspace(0.05,5,200), np.logspace(0.7,5,200)])
lhs = ive(0,xs); rhs = erf(np.pi*np.sqrt(xs/2))/np.sqrt(2*np.pi*xs)
print("I0 lower bound holds on grid:", bool(np.all(lhs>=rhs)), " min ratio", (lhs/rhs).min())

def vm_var_beta(beta):
    f = lambda t: t*t*np.exp(beta*(np.cos(t)-1)); n = lambda t: np.exp(beta*(np.cos(t)-1))
    pts = [0, min(np.pi, 10/np.sqrt(beta))]
    a = quad(f,0,np.pi,points=pts[1:],limit=200)[0]; b = quad(n,0,np.pi,points=pts[1:],limit=200)[0]
    return beta*a/b
for beta in (4,16,64,256,1000,1e4):
    print(f"beta {beta:8.0f}: beta*E_vM[th^2] = {vm_var_beta(beta):.6f}")

m=0.7; eps0=0.5; ce=np.cos(eps0)
def etaC(L,el):
    V=L*L; beta=V/el**2
    return V*(1-eps0/np.pi)*2*np.pi*beta*np.exp(-beta*(1-ce))/erf(np.pi*np.sqrt(beta/2))**2
def etainf(L,el,eta=m/4):
    V=L*L; return 4*V*np.exp(-V*ce*eta**2/(el**2*(1+np.log(L/2))/2))
def tq(el):
    a=2*np.pi**2*ce/el**2; return np.tanh(a/2)
for el in (0.34,0.47,1.0,2.0,4.0):
    print(f"el {el}: tanh(a/2) = {tq(el):.6f};  etaC L=8,12,16,32: " + ", ".join(f"{etaC(L,el):.2e}" for L in (8,12,16,32))
          + ";  eta_inf L=16,32,64,128: " + ", ".join(f"{etainf(L,el):.2e}" for L in (16,32,64,128)))
