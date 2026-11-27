import numpy as np
m=0.7; cf=0.15893895971582203; mul=1.0
def lam2(k1,k2):
    U=(1-np.cos(k1))+(1-np.cos(k2)); b=U-m; s2=np.sin(k1)**2+np.sin(k2)**2
    n=(2+b)**2+s2; bT=(b*(2+b)+s2)/n; sT2=4*s2/n**2; wT=np.sqrt(bT**2+sT2)
    return (1+bT/wT)/2
print('U_AP: ln det D_{m_f} + c_f V  (m_f=mu l/L)')
for L in (8,16,32,64,128,256):
    mf=mul/L; n=np.arange(L)
    K1,K2=np.meshgrid(2*np.pi*n/L,(2*np.pi*n+np.pi)/L)
    l2=lam2(K1,K2); sig2=(1-mf**2)*l2+mf**2
    print(L, np.sum(np.log(sig2))+cf*L*L, ' massless:',np.sum(np.log(l2))+cf*L*L)
# counts vs bounds, worst over holonomy phi
print('n0 and T0 vs bounds (H_T constants), worst phi')
for L in (32,64,128):
    V=L*L; mf=mul/L; C2=1/(1-mf**2); c1=2*np.pi*np.sqrt(C2)
    for delta in (0.25,0.1,0.05, np.sqrt(2)*np.pi/(c1*L)*1.01):
        if delta<np.sqrt(2)*np.pi/(c1*L): continue
        worstn=0; worstT=0
        for p1 in np.linspace(0,2*np.pi,13)[:-1]:
            for p2 in np.linspace(0,2*np.pi,13)[:-1]:
                n=np.arange(L); K1,K2=np.meshgrid((2*np.pi*n+p1)/L,(2*np.pi*n+np.pi+p2)/L)
                K1=(K1+np.pi)%(2*np.pi)-np.pi; K2=(K2+np.pi)%(2*np.pi)-np.pi
                sig2=(1-mf**2)*lam2(K1,K2)+mf**2
                sig=np.sqrt(sig2)
                n0=2*np.sum(sig<2*delta); T0=2*np.sum(1/sig2[sig>=2*delta])
                worstn=max(worstn,n0); worstT=max(worstT,T0)
        nb=64*C2*delta**2*V+4
        tau=4*C2+4*c1**2/np.pi**2+(9/(4*np.pi))*c1**2*np.log(np.sqrt(2)*np.pi*(1+1/L)/(c1*delta))
        Tb=V*tau+delta**-2
        print(L,'delta=%.4f'%delta,'n0 %d <= %.1f'%(worstn,nb), worstn<=nb,' T0 %.1f <= %.1f'%(worstT,Tb), worstT<=Tb)
