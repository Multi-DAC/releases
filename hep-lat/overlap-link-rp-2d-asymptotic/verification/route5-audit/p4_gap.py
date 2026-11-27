"""P4 (replaces the identity check): the consequence of Lemma 2.0 itself.
On admissible fields (max |1-U_p| < eps = 0.267), min eig H_W(-0.7)^2 must be >= mu0^2 = 0.061272,
on the full torus and for the Dirichlet compression P_B H^2 P_B with the field outside B left wild."""
import numpy as np
from lat import *
rng = np.random.default_rng(7)
m, eps, mu0sq = 0.7, 0.267, 0.061272
worst_full, worst_dir = np.inf, np.inf
for trial in range(60):
    L = [6, 8, 10][trial % 3]
    A = rng.normal(size=(2, L, L))
    # scale amplitude so that max |1-U_p| sits just under eps
    for s in np.linspace(3, 0.01, 300):
        U = np.exp(1j * s * A)
        if np.max(np.abs(1 - plaq(L, U))) < eps * 0.999:
            break
    # add a random charge sector by a uniform flux q (small) and random holonomy
    q = trial % 2
    t_, x_ = np.meshgrid(np.arange(L), np.arange(L), indexing='ij')
    U[1] *= np.exp(-2j*np.pi*q*t_/(L*L)); U[0, L-1] *= np.exp(2j*np.pi*q*np.arange(L)/L)
    if np.max(np.abs(1 - plaq(L, U))) >= eps:
        continue
    H = np.kron(np.eye(L*L), G5) @ wilson(L, U, -m)
    w = np.linalg.eigvalsh(H.conj().T @ H)
    worst_full = min(worst_full, w[0])
    # Dirichlet ball radius 1 around a site, field outside B^{+2} made wild
    Uw = U.copy()
    c = L // 2
    for mu in range(2):
        for t in range(L):
            for x in range(L):
                if max(abs(t - c), abs(x - c)) > 3:
                    Uw[mu, t, x] = np.exp(2j*np.pi*rng.random())
    H2 = (lambda Hh: Hh.conj().T @ Hh)(np.kron(np.eye(L*L), G5) @ wilson(L, Uw, -m))
    idx = [ (t*L + x)*2 + s for t in range(c-1, c+2) for x in range(c-1, c+2) for s in range(2)]
    wd = np.linalg.eigvalsh(H2[np.ix_(idx, idx)])
    worst_dir = min(worst_dir, wd[0])
print(f"min eig H_W^2 on admissible full tori: {worst_full:.4f}  (claimed floor mu0^2 = {mu0sq})")
print(f"min eig P_B H_W^2 P_B, B = 3x3 ball, field wild beyond B^+2: {worst_dir:.4f}")
