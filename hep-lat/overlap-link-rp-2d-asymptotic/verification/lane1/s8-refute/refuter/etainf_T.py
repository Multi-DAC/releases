import numpy as np
m=0.7; eps0=0.5; muT=np.sqrt(2.89/20.09); etaT=muT*(2-m)**2/8
for el in (0.34,0.47,0.4735):
    row=[]
    for L in (16,24,32,64):
        V=L*L; Gb=(1+np.log(L/2))/2
        row.append('L=%d: %.3g'%(L,4*V*np.exp(-V*np.cos(eps0)*etaT**2/(el**2*Gb))))
    print('el=',el,' eta_inf,T (proved G bound):',', '.join(row))
