"""My kernel on three of Wren's letter-11 numbers (S3 = Halbach pair, S4 = largest-centre-field pair; 10 mm cubes,
N = 12, R = 29.3 mm, z = +-12.3 mm):
  1. smallest gap between neighbouring cubes (Wren: S3 2.820 mm, S4 1.135 mm);
  2. S4/S3 when S4 alone is printed with R larger by 0.2 mm (Wren: -1.67%), and when both are larger by 0.3 mm in R
     and in z (Wren: -3.73% "isotropic", my reading of it; printed for comparison only);
  3. S4/S3 measured 2 mm and 3 mm off centre along z (Wren: -0.439%, -0.978%).
Gap: surface points of cube A (a 61 x 61 grid on each face) against the exact point-to-box distance to cube B, both
ways, for every pair of cubes whose centres are closer than 25 mm.
Usage: python -I -B w16_wren_checks.py
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import w16_prism as wp
from w16_desk_cmp import build  # noqa: E402  (prints its own table on import)

S, N, R, Z = 0.010, 12, 0.0293, 0.0123


def face_points(s, m=61):
    g = np.linspace(-s / 2, s / 2, m)
    A, B = np.meshgrid(g, g, indexing="ij")
    a, b = A.ravel(), B.ravel()
    out = []
    for ax in range(3):
        for sign in (-1, 1):
            p = np.zeros((len(a), 3))
            p[:, ax] = sign * s / 2
            p[:, (ax + 1) % 3] = a
            p[:, (ax + 2) % 3] = b
            out.append(p)
    return np.vstack(out)


def gap(blocks):
    loc = face_points(S)
    best = (1e9, None)
    for i, (ri, _, Qi, _) in enumerate(blocks):
        for j, (rj, _, Qj, _) in enumerate(blocks):
            if i == j or np.linalg.norm(ri - rj) > 0.025:
                continue
            world = loc @ Qi.T + ri
            q = (world - rj) @ Qj
            d = np.linalg.norm(np.maximum(np.abs(q) - S / 2, 0.0), axis=1).min()
            if d < best[0]:
                best = (d, (i, j))
    return best


def bx(blocks, z=0.0):
    return wp.field_many(np.array([[0.0, 0.0, z]]), blocks)[0, 0]


print()
print("--- checks of Wren's letter 11 ---")
s3 = build(S, N, R, (Z, -Z), "halbach")
s4 = build(S, N, R, (Z, -Z), "max")
for name, b in (("S3", s3), ("S4", s4)):
    d, pair = gap(b)
    print("%s smallest gap between cubes: %.3f mm (cubes %s)" % (name, d * 1e3, pair))
r0 = bx(s4) / bx(s3)
print("S4/S3 nominal: %.6f" % r0)
r1 = bx(build(S, N, R + 0.0002, (Z, -Z), "max")) / bx(s3)
print("S4 with R + 0.2 mm, S3 nominal: %.6f (%+.2f%%)" % (r1, (r1 / r0 - 1) * 100))
r2 = bx(build(S, N, R + 0.0003, (Z + 0.0003, -Z - 0.0003), "max")) / bx(s3)
print("S4 with R + 0.3 mm and z + 0.3 mm, S3 nominal: %.6f (%+.2f%%)" % (r2, (r2 / r0 - 1) * 100))
r3 = bx(build(S, N, R + 0.0002, (Z, -Z), "max")) / bx(build(S, N, R + 0.0002, (Z, -Z), "halbach"))
print("both with R + 0.2 mm: %.6f (%+.2f%%)" % (r3, (r3 / r0 - 1) * 100))
for dz in (0.002, 0.003):
    r = bx(s4, dz) / bx(s3, dz)
    print("both sensors %.0f mm off centre in z: %.6f (%+.3f%%)" % (dz * 1e3, r, (r / r0 - 1) * 100))
    r = bx(s4, dz) / bx(s3)
    print("S4 sensor %.0f mm off in z, S3 centred: %.6f (%+.3f%%)" % (dz * 1e3, r, (r / r0 - 1) * 100))
