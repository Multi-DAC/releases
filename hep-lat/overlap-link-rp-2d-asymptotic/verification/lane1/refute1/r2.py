import os  # release copy
import numpy as np, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # release copy
from r1 import DW, coulomb, GL, cf, m, s3
def ops(L, A0, A1, mf):
    N=L*L; G5=np.kron(np.eye(N),s3)
    H = G5@DW(L,np.exp(1j*A0),np.exp(1j*A1))
    w,Vv = np.linalg.eigh(H); S = (Vv*np.sign(w))@Vv.conj().T
    D0 = 0.5*(np.eye(2*N)+G5@S); D = (1-mf)*D0+mf*np.eye(2*N)
    return H, w, D
rng = np.random.default_rng(11)
L=8; V=L*L; mf=0.1
worst = {}
for trial in range(6):
    sc = [0.05,0.1,0.2,0.3,0.4,0.6][trial]
    th = rng.uniform(-1,1,size=(L,L))*sc; th -= th.mean()
    a0,a1 = coulomb(th); phi = rng.uniform(0,2*np.pi,2)
    H1,w1,D1 = ops(L,a0+phi[0]/L,a1+phi[1]/L,mf)
    H0,w0,D0 = ops(L,0*a0+phi[0]/L,0*a1+phi[1]/L,mf)
    ainf = max(np.abs(a0).max(),np.abs(a1).max()); a2 = np.sqrt((a0**2+a1**2).sum())
    dH = H1-H0
    op = np.linalg.norm(dH,2); hs = np.linalg.norm(dH)
    Delta2 = (np.abs(np.exp(1j*a0)-1)**2).sum()+(np.abs(np.exp(1j*a1)-1)**2).sum()
    gap0 = np.abs(w0).min(); gap1=np.abs(w1).min()
    dD = np.linalg.norm(D1-D0)
    sv1 = np.linalg.svd(D1,compute_uv=False); sv0=np.linalg.svd(D0,compute_uv=False)
    mir = np.sum((sv1-sv0)**2)
    l1 = np.linalg.slogdet(D1); l0=np.linalg.slogdet(D0)
    print(f"sc {sc}: ||a||inf {ainf:.3f} op {op:.4f} <= 2ainf {2*ainf:.4f} | HS^2 {hs**2:.5f} vs 2sum|D|^2 {2*Delta2:.5f} | gap0 {gap0:.4f} gap1 {gap1:.4f} | ||dD||HS {dD:.4f} <= sqrt2|a|/m {np.sqrt(2)*a2/m:.4f} (valid only if ainf<=m/4={m/4}) | Mirsky {mir:.5f}<= {dD**2:.5f} | sv range [{sv1.min():.4f},{sv1.max():.4f}] | det signs {l1[0]:.3f},{l0[0]:.3f} dlndet {l1[1]-l0[1]:+.4f}")
