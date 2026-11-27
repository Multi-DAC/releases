import numpy as np
m=0.7; el=0.47; eps0=0.5; mul=1.0; cfT=0.15893896
muT=np.sqrt(2.89/20.09); kT=(1/muT)*2*np.sqrt(2)/(2-m)**2
for L in (64,256,1024,4096,10**4,10**5,10**6):
    V=L*L; mf=mul/L; C2=1/(1-mf**2); c1=2*np.pi*np.sqrt(C2)
    n=np.arange(1,L) if L<=1024 else None
    if L<=1024:
        k=2*np.pi*np.arange(L)/L; K1,K2=np.meshgrid(k,k); kh=4*np.sin(K1/2)**2+4*np.sin(K2/2)**2
        GL=np.sum(1/kh[kh>0])/V
    else: GL=np.log(L)/(2*np.pi)+0.04876
    hE=kT*np.sqrt(3*el**2*GL/np.cos(eps0))
    delta=(hE**2/(64*C2*V))**0.25
    tau=4*C2+4*c1**2/np.pi**2+(9/(4*np.pi))*c1**2*np.log(np.sqrt(2)*np.pi*(1+1/L)/(c1*delta))
    Lam=2*np.sqrt(64*C2)*hE*L*np.log(1/mf)+4*np.log(1/mf)+2*np.log(2)*hE*np.sqrt(V*tau+delta**-2)
    print(L,'G=%.4f hE=%.3f delta=%.4f ok=%s Lambda/V=%.3f  (c_f=%.3f)'%(GL,hE,delta,delta<=0.25 and delta>=np.sqrt(2)*np.pi/(c1*L),Lam/V,cfT))
