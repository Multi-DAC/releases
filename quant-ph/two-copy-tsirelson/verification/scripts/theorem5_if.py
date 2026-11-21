# Numerical check of the "if" half of Theorem 5 as constructed in Section 5 of the paper.
# For q <= 1/sqrt(2): rho_1 = t |Phi><Phi| + (1 - t) 1/4 with t = sqrt(2) q, rho = rho_1 (x) rho_1 on
# A1 B1 A2 B2; party j measures box 0 and box 1 with the Tsirelson observables (A: Z, X; B: (Z+X)/sqrt2,
# (Z-X)/sqrt2), projector Pi^j_{b,v} = (1 + (-1)^v O^j_b)/2 for outcome value v (0 <-> +1).
# For a point gamma, V_gamma = (x)_j Pi^j_{0,gamma_j0} Pi^j_{1,gamma_j1}, and
# D(1_gamma, 1_gamma') = Re tr(rho V_gamma^dagger V_gamma'). Checks, on all of Omega:
#   the behaviour is PR_q (x) PR_q; (i) D(1_e, 1_e) = P(e) for all 256 events; (ii) D(1_e, 1_e') = 0 for all
#   22,400 locally orthogonal pairs; (iii) D is positive semidefinite (so D(1_S, 1_S) >= 0 for every S);
#   and the four-set sum of Theorem 2 equals 2 - 4 q^2.
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from base import NP, NE, EVENTS, EV_IND, ORTH_PAIRS, pr_q, val
from character_sets import S
from fractions import Fraction

I2 = np.eye(2); Z = np.diag([1.0, -1.0]); X = np.array([[0.0, 1.0], [1.0, 0.0]])
OBS = {"A": [Z, X], "B": [(Z + X) / np.sqrt(2), (Z - X) / np.sqrt(2)]}
KIND = ["A", "B", "A", "B"]


def proj(j, b, v):
    return (I2 + (-1) ** v * OBS[KIND[j]][b]) / 2


def kron(*ms):
    out = np.eye(1)
    for m_ in ms:
        out = np.kron(out, m_)
    return out


phi = np.zeros(4); phi[0] = phi[3] = 1 / np.sqrt(2)
PHI = np.outer(phi, phi)

for q in [1 / np.sqrt(2), 0.6, 0.3]:
    t = np.sqrt(2) * q
    rho1 = t * PHI + (1 - t) * np.eye(4) / 4
    rho = np.kron(rho1, rho1)
    w, U = np.linalg.eigh(rho)
    assert w.min() > -1e-12
    rhalf = U @ np.diag(np.sqrt(np.clip(w, 0, None))) @ U.T
    Xv = np.zeros((NP, 256), dtype=complex)
    for g in range(NP):
        Vg = kron(*[proj(j, 0, val(g, j, 0)) @ proj(j, 1, val(g, j, 1)) for j in range(4)])
        Xv[g] = (Vg @ rhalf).reshape(-1)
    D = np.real(Xv.conj() @ Xv.T)  # D[g, g'] = Re tr(rho V_g^dagger V_g')
    E = EV_IND.astype(float)
    De = E @ D @ E.T
    qf = Fraction(q).limit_denominator(10 ** 12)
    P = np.array([float(p) for p in pr_q(qf)])
    P_model = np.array([np.real(np.trace(rho @ kron(*[proj(j, b, o) for j, (b, o) in enumerate(zip(*EVENTS[i]))])))
                        for i in range(NE)])
    i_err = np.abs(np.diag(De) - P).max()
    ii_err = max(abs(De[i, k]) for i, k in ORTH_PAIRS)
    lam = np.linalg.eigvalsh(D).min()
    four = sum(s.astype(float) @ D @ s.astype(float) for s in S)
    print("q = %.6f: behaviour = PR_q (x) PR_q, max |diff| %.1e ; (i) max |D(1_e,1_e) - P(e)| %.1e ; "
          "(ii) max |D(1_e,1_e')| over 22,400 pairs %.1e ; (iii) least eigenvalue of D %.1e ; "
          "four-set sum %.12f, 2 - 4q^2 = %.12f"
          % (q, np.abs(P_model - P).max(), i_err, ii_err, lam, four, 2 - 4 * q * q))
print("DONE")
