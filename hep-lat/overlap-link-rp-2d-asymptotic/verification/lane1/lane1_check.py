import os  # release copy
"""LANE 1 numeric controls (sealed plan lane1-plan.md, sha e76f0c4fbe0f830a, 12:25:33 2026-11-26).

Checks the inequalities and normalisations of projects/bounties/chiral-rp/LANE1-Z-LEMMA.md on exact Wilson samples.
Conventions of te_finite_beta.py (m = 0.7, antiperiodic torus, r = 1). eps0 = 0.5, eta = m/4, mu l = 1 (m_f = 1/L).
"""
import sys, time
import numpy as np
from scipy.integrate import quad
from scipy.special import erf
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # release copy: was an absolute path
from te_finite_beta import build_DW_links, plaq_angles, wrap, sample_plaquettes, s3  # noqa: E402
from tc_holonomy import lndet_mom  # noqa: E402

m = 0.7; eps0 = 0.5; ce = np.cos(eps0); eta = m / 4; mul = 1.0; CF = 0.3131624
rng = np.random.default_rng(20261126)

def khat2(L):
    k = 2 * np.pi * np.fft.fftfreq(L)
    return 4 * np.sin(k[:, None] / 2) ** 2 + 4 * np.sin(k[None, :] / 2) ** 2

def Gex(L):
    kh = khat2(L).ravel()[1:]
    return (1 / kh).sum() / (L * L)

def coulomb(th):
    """Q = 0 Coulomb potential (a0, a1) with da = th exactly; T-e's Fourier solve with c = 0."""
    L = th.shape[0]
    F = np.fft.fft2(th)
    kt = 2 * np.pi * np.fft.fftfreq(L)[:, None]; kx = 2 * np.pi * np.fft.fftfreq(L)[None, :]
    ax = 1 - np.exp(1j * kx); at = np.exp(1j * kt) - 1
    Dn = np.abs(ax) ** 2 + np.abs(at) ** 2; Dn[0, 0] = 1.0
    B0 = F * np.conj(ax) / Dn; B1 = F * np.conj(at) / Dn; B0[0, 0] = 0; B1[0, 0] = 0
    return np.real(np.fft.ifft2(B0)), np.real(np.fft.ifft2(B1)), F

def vm_var(beta):
    w = lambda t: np.exp(beta * (np.cos(t) - 1))
    br = min(np.pi, 12 / np.sqrt(beta))
    a = quad(lambda t: t * t * w(t), 0, np.pi, points=[br], limit=400)[0]
    b = quad(w, 0, np.pi, points=[br], limit=400)[0]
    return a / b

def bound14(L, el):
    V = L * L; beta = V / el ** 2
    etaC = V * (1 - eps0 / np.pi) * 2 * np.pi * beta * np.exp(-beta * (1 - ce)) / erf(np.pi * np.sqrt(beta / 2)) ** 2
    etainf = 4 * V * np.exp(-V * ce * eta ** 2 / (el ** 2 * Gex(L)))
    aQ = 2 * np.pi ** 2 * ce / el ** 2
    return max(0.0, (1 - etaC)) * np.tanh(aQ / 2) * (2 / 3 - etainf), np.tanh(aQ / 2)

def dense(L, a0, a1, phi):
    mf = mul / L; N = L * L; G5 = np.kron(np.eye(N), s3)
    out = {}
    for tag, (A0, A1) in (("1", (a0 + phi[0] / L, a1 + phi[1] / L)),
                          ("0", (np.full((L, L), phi[0] / L), np.full((L, L), phi[1] / L)))):
        H = G5 @ build_DW_links(L, np.exp(1j * A0), np.exp(1j * A1))
        w, Vv = np.linalg.eigh(H)
        D = (1 - mf) * 0.5 * (np.eye(2 * N) + G5 @ ((Vv * np.sign(w)) @ Vv.conj().T)) + mf * np.eye(2 * N)
        sv = np.linalg.svd(D, compute_uv=False)
        sgn, ld = np.linalg.slogdet(D)
        out[tag] = dict(H=H, gap=np.abs(w).min(), D=D, sv=np.sort(sv)[::-1], ld=ld, sgn=sgn)
    return out

def n0T0_bounds(L, d):
    V = L * L; mf = mul / L; C2 = 1 / (1 - mf ** 2); c1 = np.pi * np.sqrt((2 + m * m) * C2)
    n0b = 16 * (2 + m * m) * C2 * d * d * V + 4
    valid = d >= np.sqrt(2) * np.pi / (c1 * L)
    tau = 4 * C2 + 4 * C2 * (2 + m * m) + (9 * np.pi / 4) * C2 * (2 + m * m) * np.log(np.sqrt(2) * np.pi * (1 + 1 / L) / (c1 * d))
    return n0b, (V * tau + d ** -2) if valid else np.nan, valid

if __name__ == "__main__":
    t00 = time.time()
    print("LANE 1 controls; plan sha e76f0c4fbe0f830a")
    # C3: Jensen identity
    for L in (8, 16):
        V = L * L
        vals = np.array([(lndet_mom(L, rng.uniform(0, 2 * np.pi, 2)) - 2 * V * np.log(2)) / V for _ in range(4000)])
        mu_, se_ = vals.mean(), vals.std(ddof=1) / np.sqrt(len(vals))
        print(f"C3 L={L}: mean (1/V) ln det D_0 = {mu_:.6f} +- {se_:.6f}; -c_f = {-CF:.6f}; "
              f"{'PASS' if abs(mu_ + CF) <= 3 * se_ else 'FAIL'} ({(mu_ + CF) / se_:+.2f} se)")
    # C0, C1 on a few samples
    c0 = 0.0; c1p = 0.0; c1d = 0.0
    for (L, el) in ((8, 2.0), (16, 1.0), (24, 0.47)):
        beta = (L / el) ** 2
        for _ in range(5):
            th, Q, _ = sample_plaquettes(L, beta, rng)
            if Q != 0:
                continue
            a0, a1, F = coulomb(th)
            phi = rng.uniform(0, 2 * np.pi, 2)
            U0 = np.exp(1j * (a0 + phi[0] / L)); U1 = np.exp(1j * (a1 + phi[1] / L))
            c0 = max(c0, np.abs(wrap(plaq_angles(U0, U1) - th)).max())
            kh = khat2(L); kh[0, 0] = np.inf
            pars = (np.abs(F) ** 2 / kh).sum() / (L * L)
            c1p = max(c1p, abs((a0 ** 2 + a1 ** 2).sum() - pars))
            div = (a0 - np.roll(a0, 1, axis=0)) + (a1 - np.roll(a1, 1, axis=1))
            c1d = max(c1d, np.abs(div).max())
    print(f"C0 {'PASS' if c0 < 1e-10 else 'FAIL'}: plaquette reconstruction {c0:.1e}")
    print(f"C1 {'PASS' if (c1p < 1e-9 and c1d < 1e-10) else 'FAIL'}: Parseval {c1p:.1e}, max |d*a| {c1d:.1e}")
    if c0 >= 1e-10 or c1p >= 1e-9 or c1d >= 1e-10:
        sys.exit(1)

    # P3, P4: MC cells
    print("\nP3/P4 MC cells (N = 400 draws each)")
    print(" el    L   beta     P(C)    P(Q=0|C)  tanh   P(E)    bound(1.4)  n_nu  r=<|a|^2>/(Evm V G)  +-se    <|a|^2>/bnd(K3)  frac|a|inf>eta")
    for el in (0.47, 1.0, 2.0):
        for L in (8, 16, 24, 32):
            V = L * L; beta = V / el ** 2; t0 = time.time()
            Evm = vm_var(beta); G = Gex(L); Ab = el ** 2 * G / ce; Astar = 3 * Ab
            nC = 0; nCQ = 0; nE = 0; a2s = []; ninf = 0
            for _ in range(400):
                th, Q, _ = sample_plaquettes(L, beta, rng)
                if np.abs(th).max() > eps0:
                    continue
                nC += 1
                if Q != 0:
                    continue
                nCQ += 1
                a0, a1, _ = coulomb(th)
                a2 = (a0 ** 2 + a1 ** 2).sum(); ainf = max(np.abs(a0).max(), np.abs(a1).max())
                a2s.append(a2); ninf += ainf > eta
                nE += (a2 <= Astar) and (ainf <= eta)
            pC = nC / 400; pQ = nCQ / nC if nC else np.nan; pE = nE / 400
            b14, tq = bound14(L, el)
            sePE = np.sqrt(max(pE * (1 - pE), 1e-12) / 400); sePQ = np.sqrt(max(pQ * (1 - pQ), 1e-12) / max(nC, 1))
            if len(a2s) >= 2:
                a2s = np.array(a2s); r = a2s.mean() / (Evm * V * G); ser = a2s.std(ddof=1) / np.sqrt(len(a2s)) / (Evm * V * G)
                k3 = a2s.mean() / Ab; sek3 = a2s.std(ddof=1) / np.sqrt(len(a2s)) / Ab
            else:
                r = ser = k3 = sek3 = np.nan
            scored = len(a2s) >= 100
            p3 = ("HIT" if abs(r - 1) <= 3 * ser else "MISS") if scored else "n/s"
            k3f = ("FIRED" if (k3 - 1) > 3 * sek3 else "ok") if scored else "n/s"
            p4a = ("ok" if pE >= b14 - 3 * sePE else "FIRED") if b14 > 0 else "vac"
            p4b = ("ok" if pQ >= tq - 3 * sePQ else "FIRED") if nC else "n/s"
            print(f"{el:4.2f} {L:4d} {beta:8.1f}  {pC:.3f}   {pQ:.4f}   {tq:.4f} {pE:.3f}   {b14:8.4f}   {nCQ:4d}  {r:.4f} +- {ser:.4f} [{p3}]  "
                  f"{k3:.4f} [{k3f}]   {ninf / max(nCQ, 1):.3f}   P4a {p4a} P4b {p4b}  ({time.time() - t0:.0f}s)", flush=True)

    # P1, P2, P5: dense cells
    print("\nP1/P2/P5 dense cells (20 samples in C and Q = 0, one phi each)")
    deltas = (0.05, 0.1, 0.2)
    for el in (0.47, 2.0):
        for L in (8, 12, 16, 24):
            V = L * L; beta = V / el ** 2; mf = mul / L; t0 = time.time()
            viol = []; ratios = []; diffs = []; slack01 = []; nskip = 0; c2bad = 0; n = 0
            while n < 20:
                th, Q, _ = sample_plaquettes(L, beta, rng)
                if Q != 0 or np.abs(th).max() > eps0:
                    continue
                n += 1
                a0, a1, _ = coulomb(th); phi = rng.uniform(0, 2 * np.pi, 2)
                a2 = (a0 ** 2 + a1 ** 2).sum(); ainf = max(np.abs(a0).max(), np.abs(a1).max())
                o = dense(L, a0, a1, phi); o1, o0 = o["1"], o["0"]
                for oo in (o1, o0):
                    if oo["sv"].min() < mf - 1e-10 or oo["sv"].max() > 1 + 1e-10 or abs(np.angle(oo["sgn"])) > 1e-8:
                        c2bad += 1
                hD = np.linalg.norm(o1["D"] - o0["D"]); hH2 = np.linalg.norm(o1["H"] - o0["H"]) ** 2
                lin = 2 * (np.abs(np.exp(1j * a0) - 1) ** 2 + np.abs(np.exp(1j * a1) - 1) ** 2).sum()
                if abs(hH2 - lin) > 1e-9 * max(1, lin):
                    viol.append(("P1b", hH2, lin))
                mu0 = m - 2 * ainf
                if mu0 > 0:
                    if o1["gap"] < mu0 - 1e-10:
                        viol.append(("P1a", o1["gap"], mu0))
                    bnd = (1 - mf) * np.sqrt(2) * np.sqrt(a2) / (2 * mu0)
                    ratios.append(hD / bnd)
                    if hD > bnd + 1e-10:
                        viol.append(("P1c", hD, bnd))
                else:
                    nskip += 1
                s1v, s0v = o1["sv"], o0["sv"]
                if ((s1v - s0v) ** 2).sum() > hD ** 2 + 1e-10:
                    viol.append(("P2a", ((s1v - s0v) ** 2).sum(), hD ** 2))
                lhs = np.log(s1v).sum() - np.log(s0v).sum(); diffs.append(lhs)
                for d in deltas:
                    n0 = int((s0v < 2 * d).sum()); T0 = (s0v[s0v >= 2 * d] ** -2.0).sum()
                    rhs = -(hD ** 2 / d ** 2 + n0) * np.log(1 / mf) - 2 * np.log(2) * hD * np.sqrt(T0)
                    if lhs < rhs - 1e-9:
                        viol.append(("P2b", d, lhs, rhs))
                    if d == 0.1:
                        slack01.append(lhs - rhs)
                    n0b, T0b, valid = n0T0_bounds(L, d)
                    if n0b < n0:
                        viol.append(("P2c-n0", d, n0, n0b))
                    if valid and T0b < T0:
                        viol.append(("P2c-T0", d, T0, T0b))
            diffs = np.array(diffs); md, sd = diffs.mean(), diffs.std(ddof=1) / np.sqrt(len(diffs))
            rr = (min(ratios), max(ratios)) if ratios else (np.nan, np.nan)
            print(f"el {el:4.2f} L {L:3d}: violations {len(viol)} {viol[:3]}  C2 bad {c2bad}  mu0<=0 skipped {nskip}  "
                  f"P1c ratio [{rr[0]:.3f}, {rr[1]:.3f}]  slack(delta=0.1) min {min(slack01):.2f}  "
                  f"paired diff {md:+.4f} +- {sd:.4f}  ({time.time() - t0:.0f}s)", flush=True)
    print(f"\ntotal {time.time() - t00:.0f}s")
