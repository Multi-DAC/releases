"""W16 desk build, before a parts list: (1) the field and its flatness at the final geometry, (2) the force and the
torque on every cube from all the others, (3) the axial force on a cube as it is pushed into its pocket last.

Geometry (W16-desk-envelope-1121.md, section 4): 10 mm cubes, N = 12 per ring, R = 30.3 mm, rings at z = +-12.72 mm,
Br = 1.45 T (N52, upper end). Two holders: S3 (plain Halbach, in-plane angle 2 theta) and S4 (largest centre field).

Method for the forces. Rigid uniform magnetisation, relative permeability 1 (as in w16_prism.py). The cube is then
two charge sheets +-M on its pole faces, M = Br/mu0, and in the field B of the OTHER blocks
    F   = (Br/mu0) [ int_top B dA - int_bottom B dA ],
    tau = (Br/mu0) [ int_top (r - r0) x B dA - int_bottom (r - r0) x B dA ]      (about the cube centre r0).
B is the exact block kernel; the face integrals are Gauss-Legendre, n x n points, with n varied.

Controls, printed first:
  C1  two cubes, both magnetised along z, general offset: against the closed form of Akoun and Yonnet (1984) as I
      typed it from recall. The two are independent, so agreement is the evidence for both; a disagreement would
      leave C2 and C3 standing alone.
  C2  two cubes far apart on a common axis: against the point-dipole force 3 mu0 m1 m2 / (2 pi d^4).
  C3  Newton's third law on the whole assembly: the 24 forces sum to zero, and so does the total torque
      sum(r0 x F + tau). Each force is an integral over a different cube, so this is not automatic.

Usage: python -I -B -X utf8 w16_forces.py
"""
import contextlib
import io
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import w16_prism as wp

with contextlib.redirect_stdout(io.StringIO()):          # the module prints its own tables on import
    import w16_desk_cmp as dc

MU0 = 4e-7 * np.pi
BR = dc.BR
S, N, R, Z = 0.010, 12, 0.0303, 0.01272


def gl_face(n, s):
    """Gauss-Legendre nodes and weights on a square face of edge s, in the face's own coordinates."""
    x, w = np.polynomial.legendre.leggauss(n)
    x = 0.5 * s * x; w = 0.5 * s * w
    X, Y = np.meshgrid(x, x, indexing="ij")
    W = np.outer(w, w)
    return X.ravel(), Y.ravel(), W.ravel()


def force_torque(i, blocks, n=24, shift=None):
    """Force [N] and torque about its centre [N m] on block i from all the others. `shift` moves block i first."""
    r0, dims, Q, Br = blocks[i]
    r0 = np.asarray(r0, float) + (0.0 if shift is None else np.asarray(shift, float))
    others = [b for j, b in enumerate(blocks) if j != i]
    X, Y, W = gl_face(n, dims[0])
    F = np.zeros(3); T = np.zeros(3)
    for sign in (1.0, -1.0):
        loc = np.stack([X, Y, np.full_like(X, sign * 0.5 * dims[2])], axis=1)
        arm = loc @ Q.T
        B = wp.field_many(r0[None, :] + arm, others)
        F += sign * (W[:, None] * B).sum(0)
        T += sign * (W[:, None] * np.cross(arm, B)).sum(0)
    return F * Br / MU0, T * Br / MU0


def mean_ext_field(i, blocks, n=8):
    """Volume average of the others' field over block i (Gauss-Legendre n^3)."""
    r0, dims, Q, _ = blocks[i]
    x, w = np.polynomial.legendre.leggauss(n)
    x = 0.5 * dims[0] * x; w = 0.5 * w
    X, Y, ZZ = np.meshgrid(x, x, x, indexing="ij")
    W = np.einsum("i,j,k->ijk", w, w, w).ravel()
    pts = np.asarray(r0)[None, :] + np.stack([X.ravel(), Y.ravel(), ZZ.ravel()], axis=1) @ Q.T
    others = [b for j, b in enumerate(blocks) if j != i]
    return (W[:, None] * wp.field_many(pts, others)).sum(0)


def akoun_yonnet(a, A, off, J1, J2):
    """Force on the second of two blocks, both polarised along z. Half-dimensions a (first), A (second);
    off = centre of the second minus centre of the first. Typed from recall of Akoun and Yonnet, IEEE Trans. Magn.
    20 (1984) 1962; see C1."""
    out = np.zeros(3)
    for i in (0, 1):
        for j in (0, 1):
            for k in (0, 1):
                for l in (0, 1):
                    for p in (0, 1):
                        for q in (0, 1):
                            u = off[0] + (-1) ** j * A[0] - (-1) ** i * a[0]
                            v = off[1] + (-1) ** l * A[1] - (-1) ** k * a[1]
                            w = off[2] + (-1) ** q * A[2] - (-1) ** p * a[2]
                            r = np.sqrt(u * u + v * v + w * w)
                            at = np.arctan(u * v / (r * w))
                            px = 0.5 * (v * v - w * w) * np.log(r - u) + u * v * np.log(r - v) + v * w * at + 0.5 * r * u
                            py = 0.5 * (u * u - w * w) * np.log(r - v) + u * v * np.log(r - u) + u * w * at + 0.5 * r * v
                            pz = -u * w * np.log(r - u) - v * w * np.log(r - v) + u * v * at - r * w
                            out += (-1) ** (i + j + k + l + p + q) * np.array([px, py, pz])
    return out * J1 * J2 / (4 * np.pi * MU0)


I3 = np.eye(3)
print("== C1 two 10 mm cubes, both along z, second at offset (3, 2, 14) mm: face quadrature against the closed form ==")
pair = [(np.zeros(3), (S, S, S), I3, BR), (np.array([0.003, 0.002, 0.014]), (S, S, S), I3, BR)]
ay = akoun_yonnet((S / 2,) * 3, (S / 2,) * 3, pair[1][0], BR, BR)
for n in (8, 16, 32):
    f, _ = force_torque(1, pair, n=n)
    print("  n = %2d: F = (%+.6f, %+.6f, %+.6f) N   closed form (%+.6f, %+.6f, %+.6f) N   largest difference %.1e N"
          % (n, *f, *ay, np.abs(f - ay).max()))
print("  same pair, second cube 1 mm above the first (offset (0, 0, 11) mm), the hardest case for the quadrature:")
pair = [(np.zeros(3), (S, S, S), I3, BR), (np.array([0.0, 0.0, 0.011]), (S, S, S), I3, BR)]
ay = akoun_yonnet((S / 2,) * 3, (S / 2,) * 3, pair[1][0] + np.array([1e-9, 1e-9, 0.0]), BR, BR)
for n in (8, 16, 32, 64):
    f, _ = force_torque(1, pair, n=n)
    print("  n = %2d: Fz = %+.5f N   closed form %+.5f N" % (n, f[2], ay[2]))

print("== C2 two cubes on a common axis, 200 mm apart, both along the axis: against the point dipole ==")
d = 0.200
pair = [(np.zeros(3), (S, S, S), I3, BR), (np.array([0.0, 0.0, d]), (S, S, S), I3, BR)]
f, _ = force_torque(1, pair, n=16)
m = BR * S ** 3 / MU0
print("  Fz = %+.6e N   dipole -3 mu0 m^2/(2 pi d^4) = %+.6e N   ratio %.6f" % (f[2], -3 * MU0 * m * m / (2 * np.pi * d ** 4), f[2] / (-3 * MU0 * m * m / (2 * np.pi * d ** 4))))


def flatness(name, blocks):
    b0 = wp.field_many(np.zeros((1, 3)), blocks)[0]
    print("  %s: B0 = %.2f mT (By %.4f, Bz %.4f mT)" % (name, b0[0] * 1e3, b0[1] * 1e3, b0[2] * 1e3))
    for travel in (1e-3, 2e-3, 5e-3):
        t = np.linspace(-travel, travel, 41)
        dev = []
        for axis in range(3):
            P = np.zeros((41, 3)); P[:, axis] = t
            dev.append(np.abs(wp.field_many(P, blocks)[:, 0] - b0[0]).max() / b0[0] * 100)
        g = np.linspace(-travel, travel, 9)
        G = np.stack(np.meshgrid(g, g, g, indexing="ij"), -1).reshape(-1, 3)
        cube = np.abs(wp.field_many(G, blocks)[:, 0] - b0[0]).max() / b0[0] * 100
        print("    |dBx|/B0 within +-%d mm: along x %.3f%%, y %.3f%%, z %.3f%%; anywhere in the cube of that half-edge %.3f%%"
              % (round(travel * 1e3), *dev, cube))
    return b0[0]


designs = {"S3 plain Halbach": dc.build(S, N, R, (Z, -Z), "halbach"),
           "S4 largest centre field": dc.build(S, N, R, (Z, -Z), "max")}

print("== 1. field and flatness, 10 mm cubes, N = 12, R = 30.3 mm, rings at z = +-12.72 mm, Br = %.2f T ==" % BR)
b = {k: flatness(k, v) for k, v in designs.items()}
print("  ratio S4/S3 = %.4f" % (b["S4 largest centre field"] / b["S3 plain Halbach"]))
print("  ratio if the sensor sits off-centre by the same vector in both holders (worst over a 9^3 grid):")
for travel in (1e-3, 2e-3):
    g = np.linspace(-travel, travel, 9)
    G = np.stack(np.meshgrid(g, g, g, indexing="ij"), -1).reshape(-1, 3)
    r = wp.field_many(G, designs["S4 largest centre field"])[:, 0] / wp.field_many(G, designs["S3 plain Halbach"])[:, 0]
    print("    within +-%d mm: ratio from %.4f to %.4f" % (round(travel * 1e3), r.min(), r.max()))

for name, blocks in designs.items():
    print("== 2. %s: force and torque on each cube from the other 23 ==" % name)
    print("  convergence on cube 0 (upper ring, theta = 0) and cube 1 (theta = 30 deg):")
    for n in (12, 24, 48):
        f0, t0 = force_torque(0, blocks, n=n); f1, t1 = force_torque(1, blocks, n=n)
        print("    n = %2d: |F0| %.5f N |tau0| %.5f N mm   |F1| %.5f N |tau1| %.5f N mm"
              % (n, np.linalg.norm(f0), np.linalg.norm(t0) * 1e3, np.linalg.norm(f1), np.linalg.norm(t1) * 1e3))
    FT = [force_torque(i, blocks, n=24) for i in range(len(blocks))]
    Fs = np.array([f for f, _ in FT]); Ts = np.array([t for _, t in FT])
    r0s = np.array([b_[0] for b_ in blocks])
    print("  C3 Newton: |sum F| = %.2e N (largest single |F| %.3f N); |sum (r0 x F + tau)| = %.2e N m (largest single |tau| %.4f N m)"
          % (np.linalg.norm(Fs.sum(0)), np.linalg.norm(Fs, axis=1).max(),
             np.linalg.norm((np.cross(r0s, Fs) + Ts).sum(0)), np.linalg.norm(Ts, axis=1).max()))
    print("  upper ring (z = +12.72 mm). F_r outward, F_t along increasing theta, F_z up; tau about the cube centre;")
    print("  m.B is the cosine between the cube's magnetisation and the others' field averaged over its volume.")
    print("   theta   F_r [N]   F_t [N]   F_z [N]    |F| [N]   tau_z [N mm]  |tau| [N mm]  <B_ext> [mT]  m.B cos")
    for i in range(N):
        th = 2 * np.pi * i / N
        er = np.array([np.cos(th), np.sin(th), 0.0]); et = np.array([-np.sin(th), np.cos(th), 0.0])
        f, t = Fs[i], Ts[i]
        be = mean_ext_field(i, blocks)
        e3 = blocks[i][2][:, 2]
        print("   %5.0f  %+8.3f  %+8.3f  %+8.3f  %8.3f   %+10.2f   %10.2f   %9.1f   %+.3f"
              % (np.degrees(th), f @ er, f @ et, f[2], np.linalg.norm(f), t[2] * 1e3, np.linalg.norm(t) * 1e3,
                 np.linalg.norm(be) * 1e3, be @ e3 / np.linalg.norm(be)))
    print("  largest |F| over all 24: %.3f N; largest |tau|: %.2f N mm; equivalent force couple on a 10 mm cube: %.2f N"
          % (np.linalg.norm(Fs, axis=1).max(), np.linalg.norm(Ts, axis=1).max() * 1e3,
             np.linalg.norm(Ts, axis=1).max() / S))
    print("  lower ring mirrors the upper: largest | |F_lower| - |F_upper| | = %.1e N"
          % np.abs(np.linalg.norm(Fs[N:], axis=1) - np.linalg.norm(Fs[:N], axis=1)).max())
    print("  net axial force on the upper ring from the lower ring and its own cubes: sum F_z = %+.3f N"
          " (negative: the two rings attract)" % Fs[:N, 2].sum())

    print("== 3. %s: the last cube in. Cube lifted by d along +z out of its upper-ring seat, other 23 seated ==" % name)
    print("  F_z > 0 pushes the cube out of the pocket, F_z < 0 pulls it in. d from 0 to 40 mm in 0.25 mm steps.")
    print("   theta   F_z at seat   max push out (at d)     max pull in (at d)     work to seat [mJ]   max |tau| on the way [N mm]")
    ds = np.arange(0.0, 0.04001, 0.00025)
    for i in range(N // 2 + 1):                      # theta = 0 ... 180 deg; the rest follow by symmetry
        fz = []; tq = []
        for d_ in ds:
            f, t = force_torque(i, blocks, n=16, shift=(0.0, 0.0, d_))
            fz.append(f[2]); tq.append(np.linalg.norm(t))
        fz = np.array(fz); tq = np.array(tq)
        work = -np.trapezoid(fz[::-1], ds[::-1]) * 1e3      # work done BY the hand moving the cube from d = 40 mm to 0
        print("   %5.0f   %+8.3f N   %+8.3f N (%5.2f mm)   %+8.3f N (%5.2f mm)   %+9.2f          %8.2f"
              % (np.degrees(2 * np.pi * i / N), fz[0], fz.max(), ds[fz.argmax()] * 1e3, fz.min(), ds[fz.argmin()] * 1e3,
                 work, tq.max() * 1e3))
