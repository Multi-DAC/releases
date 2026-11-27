import numpy as np, scipy.optimize as so
np.random.seed(1)
m=0.7; L=4; V=L*L; N=2*V
s1=np.array([[0,1],[1,0]],complex); s2=np.array([[0,-1j],[1j,0]]); g5=np.diag([1,-1]).astype(complex)
G=[s1,s2]; I2=np.eye(2)
def idx(x,y): return (x%L)*L+(y%L)
def DW(A, M):  # A shape (2,L,L) link angles; D_W + M ; antiperiodic in dir 1 (time)
    D=np.zeros((N,N),complex)
    for x in range(L):
        for y in range(L):
            i=idx(x,y)
            D[2*i:2*i+2,2*i:2*i+2]+= (2+M)*I2
            for mu in range(2):
                xp,yp=(x+1,y) if mu==0 else (x,y+1)
                j=idx(xp,yp); U=np.exp(1j*A[mu,x,y])
                if mu==1 and y==L-1: U=-U   # antiperiodic twist
                Pm=(I2-G[mu])/2; Pp=(I2+G[mu])/2
                D[2*i:2*i+2,2*j:2*j+2]+= -Pm*U
                D[2*j:2*j+2,2*i:2*i+2]+= -Pp*np.conj(U)
    return D
G5=np.kron(np.eye(V),g5)
def HT(A):
    D=DW(A,-m); return G5@D@np.linalg.inv(np.eye(N)+np.eye(N)+D)
# 1. ||(2+D_W(-m))^{-1}|| (2-m) over arbitrary fields: maximise
def negnorm(x):
    A=x.reshape(2,L,L); D=DW(A,2-m); return -np.linalg.svd(D,compute_uv=False).min()**-1*(2-m)
best=0
for s in range(4):
    r=so.minimize(negnorm,np.random.uniform(-np.pi,np.pi,2*V),method='L-BFGS-B',options={'maxiter':60})
    best=max(best,-r.fun)
print('max (2-m)||(2+D_W)^-1|| over U (should be <=1):',best)
# 2. min gap of H_T1 with ||a||_inf<=etaT around flat holonomy; compare muT/2
muT=np.sqrt(2.89/20.09); etaT=muT*(2-m)**2/8
def gap(x):
    a=etaT*np.tanh(x[:2*V]).reshape(2,L,L); phi=x[2*V:]
    A=a+np.array(phi)[:,None,None]/L
    H=HT(A); H=(H+H.conj().T)/2
    return np.abs(np.linalg.eigvalsh(H)).min()
worst=9
for s in range(4):
    x0=np.concatenate([np.random.randn(2*V)*2,np.random.uniform(-np.pi,np.pi,2)])
    r=so.minimize(gap,x0,method='L-BFGS-B',options={'maxiter':80})
    worst=min(worst,r.fun)
print('min gap H_T1 on sup cap (bound muT/2=%.4f):'%(muT/2),worst)
# 3. ratio ||H_T1-H_T0||_op / ||a||_inf and HS / ||a||_2
def ratios(x):
    a=x[:2*V].reshape(2,L,L); phi=x[2*V:]
    A0=np.zeros((2,L,L))+np.array(phi)[:,None,None]/L
    dH=HT(A0+a)-HT(A0)
    op=np.linalg.norm(dH,2)/np.abs(a).max(); hs=np.linalg.norm(dH)/np.linalg.norm(a)
    return op,hs
mo=0;mh=0
for s in range(4):
    for f,which in [(lambda x:-ratios(x)[0],0),(lambda x:-ratios(x)[1],1)]:
        x0=np.concatenate([np.random.randn(2*V)*0.05,np.random.uniform(-np.pi,np.pi,2)])
        r=so.minimize(f,x0,method='L-BFGS-B',options={'maxiter':40})
        if which==0: mo=max(mo,-r.fun)
        else: mh=max(mh,-r.fun)
print('max op ratio (bound %.4f):'%(4/(2-m)**2),mo,' max HS ratio (bound %.4f):'%(2*np.sqrt(2)/(2-m)**2),mh)
