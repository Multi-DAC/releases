import numpy as np
m=0.7
def lam2_T(k1,k2,m=m):
    U=(1-np.cos(k1))+(1-np.cos(k2)); b=U-m; s2=np.sin(k1)**2+np.sin(k2)**2
    n=(2+b)**2+s2; bT=(b*(2+b)+s2)/n; sT2=4*s2/n**2; wT=np.sqrt(bT**2+sT2)
    return (1+bT/wT)/2, wT, bT, sT2
def lam2_W(k1,k2,m=m):
    U=(1-np.cos(k1))+(1-np.cos(k2)); b=U-m; s2=np.sin(k1)**2+np.sin(k2)**2
    w=np.sqrt(s2+b**2); return (1+b/w)/2
def cf(fn,N):
    k=(np.arange(N)+0.5)*2*np.pi/N-np.pi
    K1,K2=np.meshgrid(k,k)
    l=fn(K1,K2); l=l[0] if isinstance(l,tuple) else l
    return -np.mean(np.log(l))
for name,fn in [('T',lam2_T),('W',lam2_W)]:
    vals=[cf(fn,N) for N in (512,1024,2048,4096)]
    rich=[(4*vals[i+1]-vals[i])/3 for i in range(3)]
    print(name,vals,rich)
# free gap min |H_T| over fine grid, several m
for mm in (0.1,0.3,0.5,0.586,0.7,0.9,1.0):
    k=np.linspace(-np.pi,np.pi,2001); K1,K2=np.meshgrid(k,k)
    _,wT,_,_=lam2_T(K1,K2,mm)
    f0=4*(1-mm)/mm**2; f4=4*(5-mm)/(mm**2+8*(1-mm))
    muT=min(mm/(2-mm), np.sqrt((mm**2+8*(1-mm))/(mm**2+8*(1-mm)+4*(5-mm))))
    print('m',mm,'min|H_T|',wT.min(),'muT',muT, 'ok', wT.min()>=muT-1e-12)
# Adams threshold
c1=2+np.sqrt(2)
for mm in (0.7,):
    epsstar=mm*(2-mm)/c1
    print('eps* =',epsstar,'Adams bound at eps=0.267:',np.sqrt(1-c1*0.267)-abs(1-mm))
    cfT=0.15893896; cfW=0.3131624
    print('el0_T(0.267)=',0.267/np.sqrt(2*cfT),' el0_T(eps*)=',epsstar/np.sqrt(2*cfT))
    print('el0_W(0.267)=',0.267/np.sqrt(2*cfW),' el0_W(eps*)=',epsstar/np.sqrt(2*cfW))
muT=np.sqrt(2.89/20.09); print('muT',muT,'etaT',muT*1.3**2/8,'hT',(1/muT)*2*np.sqrt(2)/1.3**2)
