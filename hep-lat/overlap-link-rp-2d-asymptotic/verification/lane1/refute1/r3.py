import os  # release copy
import numpy as np, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # release copy
from r2 import ops
from r1 import m, cf
def lam2(K0,K1):
    b=(1-np.cos(K0))+(1-np.cos(K1))-m; s2=np.sin(K0)**2+np.sin(K1)**2; w=np.sqrt(s2+b*b)
    return (1+b/w)/2, b
# dense vs momentum at L=4
L=4; phi=np.array([0.7,2.1])
_,_,D = ops(L,np.full((L,L),phi[0]/L),np.full((L,L),phi[1]/L),0.0)
n=np.arange(L); K0,K1=np.meshgrid((2*np.pi*n+np.pi+phi[0])/L,(2*np.pi*n+np.pi+phi[1])/L,indexing='ij')
l2,_=lam2(K0,K1)
print("L=4 dense lndet D_0", np.linalg.slogdet(D)[1], " momentum", np.log(l2).sum())
cfv = cf(4000)
print("print c_f", cfv)
# free counts
c2f = lambda mf: 1/(1-mf**2)
worst_n=0; worst_T=0
rng=np.random.default_rng(3)
for L in (16,32,64,128):
    mf=1.0/L; C2=c2f(mf); c1=np.pi*np.sqrt((2+m*m)*C2)
    for trial in range(200):
        phi = rng.uniform(0,2*np.pi,2) if trial>0 else np.array([-np.pi,-np.pi])  # trial 0: exact zero mode
        n=np.arange(L)
        k0=np.angle(np.exp(1j*(2*np.pi*n+np.pi+phi[0])/L)); k1=np.angle(np.exp(1j*(2*np.pi*n+np.pi+phi[1])/L))
        K0,K1=np.meshgrid(k0,k1,indexing='ij'); l2,b=lam2(K0,K1)
        sig2=(1-mf**2)*l2+mf**2; sig=np.sqrt(sig2).ravel()
        V=L*L
        dmin=np.sqrt(2)*np.pi/(c1*L)
        for delta in np.geomspace(max(dmin,1e-3),0.25,12):
            n0=2*np.sum(sig<2*delta); b32=16*(2+m*m)*C2*delta**2*V+4
            T0=2*np.sum(sig[sig>=2*delta]**-2.0)
            tau=4*C2+4*C2*(2+m*m)+(9*np.pi/4)*C2*(2+m*m)*np.log(np.sqrt(2)*np.pi*(1+1/L)/(c1*delta))
            b33=V*tau+delta**-2
            worst_n=max(worst_n,n0/b32); worst_T=max(worst_T,T0/b33)
print("max n0/bound(3.2):",worst_n," max T0/bound(3.3):",worst_T)
# lower bound |lambda0|^2 >= min(1/2, |k|^2/(pi^2(2+m^2))) on a fine grid
k=np.linspace(-np.pi,np.pi,1601); K0,K1=np.meshgrid(k,k,indexing='ij'); l2,b=lam2(K0,K1)
kk=K0**2+K1**2
neg=b<0
print("b<0 case: min l2*pi^2(2+m^2)/|k|^2 =", np.min(l2[neg&(kk>0)]*np.pi**2*(2+m*m)/kk[neg&(kk>0)]), " max|k_mu| on b<0:", max(np.abs(K0[neg]).max(),np.abs(K1[neg]).max()), " b>=0 min l2:", l2[~neg].min())
s2=np.sin(K0)**2+np.sin(K1)**2
print("free gap min (s2+b2)/m^2:", np.min(s2+b*b)/m**2)
