import os, sys
os.environ["OPENBLAS_NUM_THREADS"]="1"
import numpy as np
from scipy.optimize import minimize
s1=np.array([[0,1],[1,0]],complex); s2=np.array([[0,-1j],[1j,0]],complex); s3=np.diag([1,-1]).astype(complex)
g=[s1,s2]; g5=s3; I2=np.eye(2)
Pp=(I2+g5)/2; Pm=(I2-g5)/2
def DW(A0,A1,M,L):
    V=L*L; D=np.zeros((2*V,2*V),complex)
    idx=lambda t,x: (t%L)*L+(x%L)
    for t in range(L):
        for x in range(L):
            i=idx(t,x)
            D[2*i:2*i+2,2*i:2*i+2]+=(2+M)*I2
            u=np.exp(1j*A0[t,x])*(-1 if t==L-1 else 1); j=idx(t+1,x)
            D[2*i:2*i+2,2*j:2*j+2]+=-0.5*(I2-g[0])*u
            D[2*j:2*j+2,2*i:2*i+2]+=-0.5*(I2+g[0])*np.conj(u)
            u=np.exp(1j*A1[t,x]); j=idx(t,x+1)
            D[2*i:2*i+2,2*j:2*j+2]+=-0.5*(I2-g[1])*u
            D[2*j:2*j+2,2*i:2*i+2]+=-0.5*(I2+g[1])*np.conj(u)
    return D
def DWF(A0,A1,m,mf,L5,L):
    V=L*L; D4=DW(A0,A1,-m,L); n=2*V
    D=np.zeros((n*L5,n*L5),complex)
    PM=np.kron(np.eye(V),Pm); PP=np.kron(np.eye(V),Pp)
    for s in range(L5):
        D[s*n:(s+1)*n,s*n:(s+1)*n]=D4+np.eye(n)
        if s+1<L5: D[s*n:(s+1)*n,(s+1)*n:(s+2)*n]+=-PM
        else: D[s*n:(s+1)*n,0:n]+=mf*PM
        if s-1>=0: D[s*n:(s+1)*n,(s-1)*n:s*n]+=-PP
        else: D[s*n:(s+1)*n,(L5-1)*n:L5*n]+=mf*PP
    return D
def lndet(A0,A1,m,mf,L5,L):
    s,l=np.linalg.slogdet(DWF(A0,A1,m,mf,L5,L)); return s,l
L=int(sys.argv[1]) if len(sys.argv)>1 else 4
m=0.7
rng=np.random.default_rng(1)
for mf,L5 in [(0.1,1),(0.1,2),(1.0,2),(0.1,3)]:
    sAP,lAP=lndet(np.zeros((L,L)),np.full((L,L),np.pi/L),m,mf,L5,L)
    # flat holonomy grid
    best=-1e9; arg=None
    for p0 in np.linspace(0,2*np.pi,9,endpoint=False):
        for p1 in np.linspace(0,2*np.pi,9,endpoint=False):
            s,l=lndet(np.full((L,L),p0/L),np.full((L,L),p1/L),m,mf,L5,L)
            if l>best: best=l; arg=(p0,p1,s)
    # uniform flux Q
    fl=[]
    for Q in range(0,L*L//2+1, max(1,L*L//8)):
        f=2*np.pi*Q/(L*L)
        A1=np.array([[f*t for x in range(L)] for t in range(L)]); A0=np.zeros((L,L)); A0[L-1,:]=-f*L*np.arange(L)
        s,l=lndet(A0,A1,m,mf,L5,L); fl.append((Q,round(l-lAP,4),round(s.real,3)))
    # adversarial ascent
    def negf(v):
        s,l=lndet(v[:L*L].reshape(L,L),v[L*L:].reshape(L,L),m,mf,L5,L); return -l
    bestasc=-1e9
    for trial in range(3):
        v0=rng.uniform(-np.pi,np.pi,2*L*L)
        r=minimize(negf,v0,method='L-BFGS-B',options={'maxiter':200})
        bestasc=max(bestasc,-r.fun)
    print(f"L={L} mf={mf} L5={L5}: lnAP={lAP:.6f} sAP={sAP.real:+.0f} | flat max-lnAP={best-lAP:+.3e} at {arg[0]:.2f},{arg[1]:.2f} | flux Q:(dl,sign) {fl} | ascent max-lnAP={bestasc-lAP:+.3e}",flush=True)
