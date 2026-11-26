"""Consultant's item 2 (after clock reading 16:41): is the residue of D_N(z, 0)^-1 at z0 ALONE 2m P0^+?
A fixed circle about z = 1 encloses z0 and 1/z0 both. Proposed repair: D_N(z, 0) commutes with gamma0 = sigma1, so
D_N = d+ P0^+ + d- P0^-, and the inner zero z0 is in the + channel.
D_N at U = 1 from its closed form (merged text, Prop 5 Step 1):
D_N = (P+ + P- T5^L5)(I + T5^L5)^-1, with B = I + a5 X, X = b I + s0 sigma1 + sin k1 * (chiral off-diagonal part),
b = c(k1) - (z + 1/z)/2, s0 = (z - 1/z)/2, gamma5 = sigma3, gamma0 = sigma1."""
import numpy as np

s1 = np.array([[0, 1], [1, 0]], complex)
Pp, Pm = np.diag([1, 0]).astype(complex), np.diag([0, 1]).astype(complex)   # chiral P+ and P- (gamma5 = sigma3)
P0p = (np.eye(2) + s1) / 2


def DN(z, k1, m, a5, L5):
    c = 2 - m - np.cos(k1)
    b, s0 = c - (z + 1 / z) / 2, (z - 1 / z) / 2
    Xpp = Xmm = b
    Xpm, Xmp = s0 + np.sin(k1), s0 - np.sin(k1)
    Bpp, Bmm, Bpm, Bmp = 1 + a5 * Xpp, 1 + a5 * Xmm, a5 * Xpm, a5 * Xmp
    T = np.array([[1 / Bpp, -Bpm / Bpp], [Bmp / Bpp, Bmm - Bmp * Bpm / Bpp]], complex)
    TL = np.linalg.matrix_power(T, L5)
    return (Pp + Pm @ TL) @ np.linalg.inv(np.eye(2) + TL), T


def detDN(z, *a):
    return np.linalg.det(DN(z, *a)[0])


for m, a5, S in ((0.7, 0.2, 8), (0.7, 0.1, 8), (0.7, 0.05, 8), (0.7, 0.1, 16), (0.3, 0.1, 16), (0.3, 0.05, 32)):
    L5 = int(round(S / a5))
    args = (0.0, m, a5, L5)
    # sigma1 symmetry of T5 and D_N at a test point
    zt = 0.93 * np.exp(0.2j)
    D, T = DN(zt, *args)
    symT = np.abs(s1 @ T @ s1 - np.linalg.inv(T)).max()
    symD = np.abs(s1 @ D @ s1 - D).max()
    # find z0 in (r_mid, 1) on the real axis: zero of det D_N, by bisection on the + channel d+ = <+|D_N|+> in sigma1 basis
    def dpm(z):
        Dz = DN(z, *args)[0]
        U = np.array([[1, 1], [1, -1]]) / np.sqrt(2)
        Ds = U.conj().T @ Dz @ U
        return Ds[0, 0].real, Ds[1, 1].real, abs(Ds[0, 1]) + abs(Ds[1, 0])
    xs = 1 - np.logspace(-1, -14, 4000)
    vals = np.array([dpm(x)[:2] for x in xs])
    roots = []
    for ch in (0, 1):
        sg = np.sign(vals[:, ch])
        idx = np.where(sg[:-1] != sg[1:])[0]
        for i in idx:
            lo, hi = xs[i], xs[i + 1]
            for _ in range(200):
                mid = (lo + hi) / 2
                if np.sign(dpm(mid)[ch]) == np.sign(dpm(lo)[ch]):
                    lo = mid
                else:
                    hi = mid
            roots.append(("+" if ch == 0 else "-", (lo + hi) / 2))
    z0 = [r for ch, r in roots if ch == "+"]
    z0m = [r for ch, r in roots if ch == "-"]
    line = "m=%.1f a5=%.2f S=%d L5=%d: |s1 T s1 - T^-1| %.1e, |s1 D s1 - D| %.1e; zeros in (0.9, 1): + %s, - %s" % (
        m, a5, S, L5, symT, symD, ["%.10f" % r for r in z0], ["%.10f" % r for r in z0m])
    if len(z0) == 1:
        z0 = z0[0]
        rad = (1 - z0) / 4
        th = np.linspace(0, 2 * np.pi, 4001)[:-1]
        R = sum(np.linalg.inv(DN(z0 + rad * np.exp(1j * t), *args)[0]) * rad * np.exp(1j * t) for t in th) / len(th)
        Rz = sum(np.linalg.inv(DN(z0 + rad * np.exp(1j * t), *args)[0]) * rad * np.exp(1j * t) / (z0 + rad * np.exp(1j * t))
                 for t in th) / len(th)
        line += "; residue at z0 alone (dz): %s ; /(2m P0^+) dev %.2e ; (dz/z) dev %.2e" % (
            np.round(R, 5).tolist(), np.abs(R - 2 * m * P0p).max() / (2 * m), np.abs(Rz - 2 * m * P0p).max() / (2 * m))
    print(line, flush=True)
print("DONE")
