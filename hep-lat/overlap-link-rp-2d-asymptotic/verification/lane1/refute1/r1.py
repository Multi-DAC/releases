import numpy as np, sys
m = 0.7
s1 = np.array([[0,1],[1,0]],complex); s2 = np.array([[0,-1j],[1j,0]],complex); s3 = np.diag([1.,-1.]).astype(complex)
I2 = np.eye(2); g = [s1, s2]

def DW(L, U0, U1):
    N = L*L; D = np.zeros((2*N,2*N),complex); idx = lambda t,x: (t%L)*L+(x%L)
    for t in range(L):
        for x in range(L):
            i = idx(t,x); D[2*i:2*i+2,2*i:2*i+2] += (2-m)*I2
            for mu,(dt,dx) in enumerate([(1,0),(0,1)]):
                tt,xx = t+dt,x+dx; u = (U0 if mu==0 else U1)[t,x]
                if (mu==0 and tt==L) or (mu==1 and xx==L): u = -u
                j = idx(tt,xx)
                D[2*i:2*i+2,2*j:2*j+2] += -0.5*(I2-g[mu])*u
                D[2*j:2*j+2,2*i:2*i+2] += -0.5*(I2+g[mu])*np.conj(u)
    return D

def cf(n=4000):
    k = 2*np.pi*(np.arange(n)+0.5)/n - np.pi
    K0,K1 = np.meshgrid(k,k,indexing='ij')
    b = (1-np.cos(K0))+(1-np.cos(K1))-m; s2_ = np.sin(K0)**2+np.sin(K1)**2; w = np.sqrt(s2_+b*b)
    return -np.mean(np.log((1+b/w)/2))

def GL(L):
    k = 2*np.pi*np.arange(L)/L; K0,K1 = np.meshgrid(k,k,indexing='ij')
    kh = 4*np.sin(K0/2)**2+4*np.sin(K1/2)**2; kh[0,0]=np.inf
    return np.sum(1/kh)/(L*L)

print("c_f(0.7) midpoint n=4000:", cf(4000), " n=8000:", cf(8000))
for L in (16,64):
    print(f"G({L}) = {GL(L):.4f}, ln L/2pi+0.04876 = {np.log(L)/(2*np.pi)+0.04876:.4f}, proved bound {(1+np.log(L/2))/2:.4f}")

# Coulomb potential from theta (sum 0), plaquette convention of te_finite_beta
def coulomb(th):
    L = th.shape[0]
    F = np.fft.fft2(th)
    kt = 2*np.pi*np.fft.fftfreq(L)[:,None]; kx = 2*np.pi*np.fft.fftfreq(L)[None,:]
    ax = 1-np.exp(1j*kx); at = np.exp(1j*kt)-1
    Dn = np.abs(ax)**2+np.abs(at)**2; Dn[0,0]=1
    B0 = F*np.conj(ax)/Dn; B1 = F*np.conj(at)/Dn; B0[0,0]=0; B1[0,0]=0
    return np.real(np.fft.ifft2(B0)), np.real(np.fft.ifft2(B1))

def plaq(U0,U1):
    return np.angle(U0*np.roll(U1,-1,axis=0)*np.conj(np.roll(U0,-1,axis=1))*np.conj(U1))

# Lip^2 of a_l: a is linear in theta; take gradient of a_(0,0,mu) wrt theta via unit vectors projected on H
L = 8; V = L*L
M0 = np.zeros((V,V)); M1 = np.zeros((V,V))
for p in range(V):
    e = np.zeros(V); e[p]=1; e -= 1/V
    a0,a1 = coulomb(e.reshape(L,L)); M0[:,p]=a0.ravel(); M1[:,p]=a1.ravel()
# M maps theta (on H, via projection) to a; rows = links
lip0 = np.linalg.norm(M0[0]); lip1 = np.linalg.norm(M1[0])
print(f"L=8 Lip^2 a_l: dir0 {lip0**2:.6f} dir1 {lip1**2:.6f}  G/2 = {GL(8)/2:.6f}")
# check ||a||^2 = <theta,G theta>
rng = np.random.default_rng(5)
th = rng.normal(size=(L,L))*0.1; th -= th.mean()
a0,a1 = coulomb(th); U0=np.exp(1j*a0); U1=np.exp(1j*a1)
print("plaq recon err", np.abs(plaq(U0,U1)-th).max())
Ft = np.fft.fft2(th)/L; k = 2*np.pi*np.arange(L)/L; K0,K1=np.meshgrid(k,k,indexing='ij'); kh=4*np.sin(K0/2)**2+4*np.sin(K1/2)**2; kh[0,0]=np.inf
print("||a||^2", (a0**2+a1**2).sum(), " <th,G th>", np.sum(np.abs(Ft)**2/kh))
