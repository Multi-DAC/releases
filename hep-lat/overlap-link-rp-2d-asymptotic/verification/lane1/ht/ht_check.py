import os  # release copy
"""Sealed check of draft-section8-HT.md (plan ht-plan.md). H_T = g5 D_W(-m) (2 + D_W(-m))^{-1}, a5 = 1."""
import sys, time
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))  # release copy
pass  # release copy: lane1_check sits in the same folder as te_finite_beta
from te_finite_beta import build_DW_links, sample_plaquettes, s3  # noqa: E402
from lane1_check import coulomb  # noqa: E402

m = 0.7; eps0 = 0.5; mul = 1.0; MUT = 0.3793; KOP = 4 / (2 - m) ** 2; KHS = 2 * np.sqrt(2) / (2 - m) ** 2
rng = np.random.default_rng(20261127)
c = 4 / (25 * np.pi ** 2)

def opHT(L, A0, A1):
    N = L * L; G5 = np.kron(np.eye(N), s3); I = np.eye(2 * N)
    DW = build_DW_links(L, np.exp(1j * A0), np.exp(1j * A1))
    HT = G5 @ DW @ np.linalg.inv(2 * I + DW)
    HT = 0.5 * (HT + HT.conj().T)
    w, V = np.linalg.eigh(HT)
    mf = mul / L
    D = (1 - mf) * 0.5 * (I + G5 @ ((V * np.sign(w)) @ V.conj().T)) + mf * I
    sv = np.sort(np.linalg.svd(D, compute_uv=False))[::-1]
    return HT, np.abs(w).min(), D, sv

def bounds(L, d):
    V = L * L; mf = mul / L; C2 = 1 / (1 - mf ** 2); c1 = np.sqrt(C2 / c)
    n0b = 16 * c1 ** 2 * d * d * V / np.pi ** 2 + 4
    valid = d >= np.sqrt(2) * np.pi / (c1 * L)
    tau = 4 * C2 + 4 * c1 ** 2 / np.pi ** 2 + (9 / (4 * np.pi)) * c1 ** 2 * np.log(np.sqrt(2) * np.pi * (1 + 1 / L) / (c1 * d))
    return n0b, (V * tau + d ** -2), valid

if __name__ == "__main__":
    print("H_T check; plan ht-plan.md")
    for el in (0.47, 2.0):
        for L in (8, 12, 16):
            V = L * L; beta = V / el ** 2; mf = mul / L; t0 = time.time()
            viol = []; diffs = []; gaps0 = []; herm = 0.0; n = 0
            while n < 10:
                th, Q, _ = sample_plaquettes(L, beta, rng)
                if Q != 0 or np.abs(th).max() > eps0:
                    continue
                n += 1
                a0, a1, _ = coulomb(th); phi = rng.uniform(0, 2 * np.pi, 2)
                a2 = (a0 ** 2 + a1 ** 2).sum(); ainf = max(np.abs(a0).max(), np.abs(a1).max())
                H1, g1, D1, s1 = opHT(L, a0 + phi[0] / L, a1 + phi[1] / L)
                H0, g0, D0, s0 = opHT(L, np.full((L, L), phi[0] / L), np.full((L, L), phi[1] / L))
                gaps0.append(g0)
                if g0 < MUT:
                    viol.append(("Q1", g0))
                dH = H1 - H0
                if np.linalg.norm(dH, 2) > KOP * ainf + 1e-10:
                    viol.append(("Q2op", np.linalg.norm(dH, 2), KOP * ainf))
                if np.linalg.norm(dH) > KHS * np.sqrt(a2) + 1e-10:
                    viol.append(("Q2hs", np.linalg.norm(dH), KHS * np.sqrt(a2)))
                if g1 < g0 - KOP * ainf - 1e-10:
                    viol.append(("Q3", g1, g0 - KOP * ainf))
                hD = np.linalg.norm(D1 - D0)
                mu0 = MUT - KOP * ainf
                if mu0 > 0 and hD > (1 - mf) * KHS * np.sqrt(a2) / (2 * mu0) + 1e-10:
                    viol.append(("Q4", hD))
                if ((s1 - s0) ** 2).sum() > hD ** 2 + 1e-10:
                    viol.append(("Q6mirsky",))
                lhs = np.log(s1).sum() - np.log(s0).sum(); diffs.append(lhs)
                for d in (0.05, 0.1, 0.2):
                    n0 = int((s0 < 2 * d).sum()); T0 = (s0[s0 >= 2 * d] ** -2.0).sum()
                    rhs = -(hD ** 2 / d ** 2 + n0) * np.log(1 / mf) - 2 * np.log(2) * hD * np.sqrt(T0)
                    if lhs < rhs - 1e-9:
                        viol.append(("Q6", d))
                    n0b, T0b, valid = bounds(L, d)
                    if n0b < n0:
                        viol.append(("Q5n0", d, n0, n0b))
                    if valid and T0b < T0:
                        viol.append(("Q5T0", d, T0, T0b))
            diffs = np.array(diffs)
            print(f"el {el:4.2f} L {L:3d}: violations {len(viol)} {viol[:3]}  min free gap {min(gaps0):.4f}  "
                  f"paired diff {diffs.mean():+.4f} +- {diffs.std(ddof=1) / np.sqrt(len(diffs)):.4f}  ({time.time() - t0:.0f}s)", flush=True)
