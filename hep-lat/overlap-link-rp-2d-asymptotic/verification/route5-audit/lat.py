"""Independent 2D U(1) lattice operators for the Route (v) audit (Clawd, 2026-11-26).
Written from the conventions stated in LANE2-ROUTE5-AUX-RP-PROOF.md section 1, without reading
Wren's verify script. Sites (t,x), index i = t*L + x; spin index last.
gamma0 = sigma1, gamma1 = sigma2, gamma5 = sigma3. AP in time (seam on t=L-1 -> 0), periodic in space.
"""
import numpy as np

s0 = np.eye(2, dtype=complex)
s1 = np.array([[0, 1], [1, 0]], dtype=complex)
s2 = np.array([[0, -1j], [1j, 0]], dtype=complex)
s3 = np.array([[1, 0], [0, -1]], dtype=complex)
G = [s1, s2]
G5 = s3


def hop(L, U, mu):
    """T_mu as a V x V matrix: (T_mu psi)(n) = U_mu(n) psi(n+mu), AP seam in time."""
    V = L * L
    T = np.zeros((V, V), dtype=complex)
    for t in range(L):
        for x in range(L):
            i = t * L + x
            if mu == 0:
                j = ((t + 1) % L) * L + x
                ph = U[0, t, x] * (-1.0 if t == L - 1 else 1.0)
            else:
                j = t * L + (x + 1) % L
                ph = U[1, t, x]
            T[i, j] = ph
    return T


def wilson(L, U, M):
    V = L * L
    D = (2.0 + M) * np.kron(np.eye(V), s0)
    for mu in range(2):
        T = hop(L, U, mu)
        D -= 0.5 * (np.kron(T, s0 - G[mu]) + np.kron(T.conj().T, s0 + G[mu]))
    return D


def dwf(L, U, m, mf, L5):
    """Shamir DWF: (D psi)_s = (D_W(-m)+1) psi_s - P_- psi_{s+1} - P_+ psi_{s-1},
    psi_{L5+1} = -mf psi_1, psi_0 = -mf psi_{L5}."""
    V = L * L
    n = 2 * V
    Dw = wilson(L, U, -m)
    Pp = np.kron(np.eye(V), (s0 + G5) / 2)
    Pm = np.kron(np.eye(V), (s0 - G5) / 2)
    D = np.zeros((L5 * n, L5 * n), dtype=complex)
    I = np.eye(n)
    for s in range(L5):
        D[s*n:(s+1)*n, s*n:(s+1)*n] = Dw + I
        # -P_- psi_{s+1}
        if s + 1 < L5:
            D[s*n:(s+1)*n, (s+1)*n:(s+2)*n] += -Pm
        else:
            D[s*n:(s+1)*n, 0:n] += mf * Pm
        # -P_+ psi_{s-1}
        if s - 1 >= 0:
            D[s*n:(s+1)*n, (s-1)*n:s*n] += -Pp
        else:
            D[s*n:(s+1)*n, (L5-1)*n:L5*n] += mf * Pp
    return D


def herm_func(H, f):
    w, Q = np.linalg.eigh(H)
    return (Q * f(w)) @ Q.conj().T


def HT(L, U, m):
    V = L * L
    Dw = wilson(L, U, -m)
    H = np.kron(np.eye(V), G5) @ Dw @ np.linalg.inv(2 * np.eye(2 * V) + Dw)
    return 0.5 * (H + H.conj().T)


def overlap_HT(L, U, m, mf):
    V = L * L
    H = HT(L, U, m)
    g5 = np.kron(np.eye(V), G5)
    S = herm_func(H, np.sign)
    return (1 - mf) * (np.eye(2 * V) + g5 @ S) / 2 + mf * np.eye(2 * V)


def DN(L, U, m, mf, L5):
    V = L * L
    H = HT(L, U, m)
    g5 = np.kron(np.eye(V), G5)
    S = herm_func(H, lambda w: np.tanh(L5 * np.arctanh(w)))  # tanh(L5 * Htilde/2), Htilde = 2 artanh H
    return (1 - mf) * (np.eye(2 * V) + g5 @ S) / 2 + mf * np.eye(2 * V)


def logdet(A):
    s, l = np.linalg.slogdet(A)
    return s, l


def U_AP(L):
    U = np.ones((2, L, L), dtype=complex)
    U[1] = np.exp(1j * np.pi / L)
    return U


def plaq(L, U):
    """theta_p for p at (t,x): U0(t,x) U1(t+1,x) U0(t,x+1)^* U1(t,x)^*  (gauge plaquette, no fermion seam)."""
    P = np.zeros((L, L), dtype=complex)
    for t in range(L):
        for x in range(L):
            P[t, x] = U[0, t, x] * U[1, (t+1) % L, x] * np.conj(U[0, t, (x+1) % L]) * np.conj(U[1, t, x])
    return P


# ---- link reflection (proof section 1) ----
def full_from_halves(L, up, um_img):
    """Build U with U_+ = up and U_- = theta(um_img); cut links = 1.
    up, um_img: dict with 'U1' shape (L/2, L) for t=0..L/2-1 and 'U0' shape (L/2-1, L) for t=0..L/2-2.
    theta: (thetaU)_1(t,x) = U_1(L-1-t,x); (thetaU)_0(t,x) = U_0((L-2-t) mod L, x)^*.
    So U_- with theta(U_-) = um_img means U_1(L-1-t,x) = um_img.U1[t], U_0(L-2-t,x) = conj(um_img.U0[t])."""
    U = np.ones((2, L, L), dtype=complex)
    h = L // 2
    for t in range(h):
        U[1, t] = up['U1'][t]
        U[1, L - 1 - t] = um_img['U1'][t]
    for t in range(h - 1):
        U[0, t] = up['U0'][t]
        U[0, L - 2 - t] = np.conj(um_img['U0'][t])
    # cut links U0(h-1,x), U0(L-1,x) stay 1
    return U


def rand_half(L, rng, spread):
    h = L // 2
    return {'U1': np.exp(1j * spread * rng.uniform(-np.pi, np.pi, (h, L))),
            'U0': np.exp(1j * spread * rng.uniform(-np.pi, np.pi, (h - 1, L)))}
