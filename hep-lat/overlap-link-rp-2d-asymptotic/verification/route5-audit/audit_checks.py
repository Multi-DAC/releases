"""Route (v) audit: P0 sanity (Schur identity, my conventions), P1 diamagnetic, P2 DWF link-RP kernel,
P3 control (overlap / DW-over-PV kernel must be able to fail), P4 local gap on admissible fields."""
import sys
import numpy as np
from lat import *

m, mf = 0.7, 0.2
rng = np.random.default_rng(20261126)


def haar(L, spread=1.0):
    return np.exp(1j * spread * rng.uniform(-np.pi, np.pi, (2, L, L)))


print("== P0: det D_DW / det D_PV = det D_N (validates my DWF/H_T conventions) ==")
for L, L5 in [(4, 2), (4, 3), (6, 3)]:
    for k in range(3):
        U = haar(L)
        s1_, l1 = logdet(dwf(L, U, m, mf, L5))
        s2_, l2 = logdet(dwf(L, U, m, 1.0, L5))
        s3_, l3 = logdet(DN(L, U, m, mf, L5))
        print(f"  L={L} L5={L5}: sign DW={s1_.real:+.0f} PV={s2_.real:+.0f}  "
              f"ln(DW/PV)-ln DN = {l1 - l2 - l3:+.2e}")

print("== P1: diamagnetic ln det D_DW[U] - ln det D_DW[U_AP] <= 0 ==")
for L, L5 in [(6, 2), (6, 3), (6, 4), (8, 3)]:
    lAP = logdet(dwf(L, U_AP(L), m, mf, L5))[1]
    worst = -np.inf
    nonpos = 0
    for k in range(40):
        spread = [1.0, 0.3, 0.05][k % 3]
        U = haar(L, spread)
        s, l = logdet(dwf(L, U, m, mf, L5))
        if s.real <= 0 or abs(s.imag) > 1e-8:
            nonpos += 1
        worst = max(worst, l - lAP)
    # flat fields with every holonomy pair: should all be <= AP, equality only at AP
    flat_worst = -np.inf
    for a in np.linspace(0, 2*np.pi, 9)[:-1]:
        for b in np.linspace(0, 2*np.pi, 9)[:-1]:
            U = np.ones((2, L, L), dtype=complex)
            U[0] *= np.exp(1j * a / L)
            U[1] *= np.exp(1j * b / L)
            l = logdet(dwf(L, U, m, mf, L5))[1]
            flat_worst = max(flat_worst, l - lAP) if abs(b - np.pi) > 1e-9 or abs(a) > 1e-9 else flat_worst
    print(f"  L={L} L5={L5}: nonpositive dets={nonpos}/40, max over random = {worst:+.4f}, "
          f"max over 63 non-AP flat fields = {flat_worst:+.4f}")


def rp_kernel(L, halves, weight):
    n = len(halves)
    K = np.zeros((n, n), dtype=complex)
    for i in range(n):
        for j in range(n):
            K[i, j] = weight(full_from_halves(L, halves[j], halves[i]))
    return K


def report(name, K):
    herm = np.max(np.abs(K - K.conj().T)) / np.max(np.abs(K))
    w = np.linalg.eigvalsh(0.5 * (K + K.conj().T))
    print(f"  {name}: herm_err={herm:.1e}, min/max eig = {w[0]/w[-1]:+.3e}, "
          f"#neg(< -1e-10 rel) = {np.sum(w < -1e-10 * w[-1])}/{len(w)}")


print("== P2/P3: link-RP kernel K_ij = W(U_+=u_j, theta U_- = u_i, U_cut=1), beta=0 ==")
for L, L5, spread, n in [(4, 2, 1.0, 30), (4, 3, 0.5, 30), (6, 2, 0.5, 30), (6, 3, 1.0, 24)]:
    halves = [rand_half(L, rng, spread) for _ in range(n)]
    # normalise by AP to keep magnitudes sane
    lAP = logdet(dwf(L, U_AP(L), m, mf, L5))[1]
    lAPpv = logdet(dwf(L, U_AP(L), m, 1.0, L5))[1]
    lAPov = logdet(overlap_HT(L, U_AP(L), m, mf))[1]
    lAPw = logdet(wilson(L, U_AP(L), 0.5))[1]

    def w_dw(U):
        s, l = logdet(dwf(L, U, m, mf, L5)); return s * np.exp(l - lAP)

    def w_ratio(U):
        s, l = logdet(dwf(L, U, m, mf, L5)); s2_, l2 = logdet(dwf(L, U, m, 1.0, L5))
        return (s / s2_) * np.exp(l - l2 - lAP + lAPpv)

    def w_ov(U):
        s, l = logdet(overlap_HT(L, U, m, mf)); return s * np.exp(l - lAPov)

    def w_wil(U):  # Wilson r=1 positive mass: known RP, a second positive control
        s, l = logdet(wilson(L, U, 0.5)); return s * np.exp(l - lAPw)

    print(f" L={L} L5={L5} spread={spread} n={n}")
    report("Wilson r=1, M=+0.5 (known RP)", rp_kernel(L, halves, w_wil))
    report("DWF numerator det D_DW (claim: RP)", rp_kernel(L, halves, w_dw))
    report("DW/PV ratio = det D_N (control)", rp_kernel(L, halves, w_ratio))
    report("overlap det D_ov (control)", rp_kernel(L, halves, w_ov))
    sys.stdout.flush()
