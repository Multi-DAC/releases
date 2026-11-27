import numpy as np
m=0.7
k=np.linspace(-np.pi,np.pi,1601); K0,K1=np.meshgrid(k,k,indexing='ij')
b=(1-np.cos(K0))+(1-np.cos(K1))-m; s2=np.sin(K0)**2+np.sin(K1)**2; w=np.sqrt(s2+b*b)
l2=(1+b/w)/2; l2s=s2/(2*w*(w-b))
kk=K0**2+K1**2; neg=(b<0)&(kk>0)
r=l2[neg]*np.pi**2*(2+m*m)/kk[neg]; rs=l2s[neg]*np.pi**2*(2+m*m)/kk[neg]
i=np.argmin(r); print("naive min",r[i],"at",K0[neg][i],K1[neg][i],"l2",l2[neg][i],"stable",l2s[neg][i]," stable-formula min ratio",rs.min())
