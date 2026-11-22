"""Same cubes, same pockets' centres, two orientation rules: the plain Halbach rule (in-plane angle 2 theta) against
the largest-centre-field rule (each cube magnetised along the field a unit x-dipole at the centre would make at the
cube: direction 3 (n.x) n - x, n the unit vector from the centre to the cube). The second rule is the known optimum
for the centre field with no homogeneity condition (reciprocity; Insinga, Bjork, Smith, Bahl 2016); it is evaluated
here with the exact block kernel, the cube turned with its magnetisation.
Designs: 10 mm cubes, N = 12, R = 29.3 mm, one ring; and the pair at z = +-12.3 mm. Also 10 mm, N = 8, R = 23.3 mm.
Usage: python -I -B w16_desk_cmp.py
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import w16_prism as wp

BR, TRAVEL = 1.45, 5e-3


def frame_for(m):
    e3 = m / np.linalg.norm(m)
    ref = np.array([0.0, 0.0, 1.0]) if abs(e3[2]) < 0.9 else np.array([1.0, 0.0, 0.0])
    e1 = np.cross(ref, e3); e1 /= np.linalg.norm(e1)
    return np.stack([e1, np.cross(e3, e1), e3], axis=1)


def build(s, n, R, zs, rule):
    blocks = []
    for z0 in zs:
        for i in range(n):
            th = 2 * np.pi * i / n
            r0 = np.array([R * np.cos(th), R * np.sin(th), z0])
            if rule == "halbach":
                m = np.array([np.cos(2 * th), np.sin(2 * th), 0.0])
            else:
                nn = r0 / np.linalg.norm(r0)
                m = 3 * nn[0] * nn - np.array([1.0, 0.0, 0.0])
            blocks.append((r0, (s, s, s), frame_for(m), BR))
    return blocks


def report(name, blocks):
    b0 = wp.field_many(np.zeros((1, 3)), blocks)[0]
    t = np.linspace(-TRAVEL, TRAVEL, 21)
    dev = []
    for axis in range(3):
        P = np.zeros((21, 3)); P[:, axis] = t
        dev.append(np.abs(wp.field_many(P, blocks)[:, 0] - b0[0]).max() / b0[0] * 100)
    tilt = max(abs(np.degrees(np.arcsin(Q[2, 2]))) for _, _, Q, _ in blocks)
    print("  %-22s B0 = %7.2f mT (By %.3f, Bz %.3f); |dBx|/B0 over +-5 mm: x %.2f%%, y %.2f%%, z %.2f%%; largest tilt out of plane %.1f deg"
          % (name, b0[0] * 1e3, b0[1] * 1e3, b0[2] * 1e3, *dev, tilt))
    return b0[0]


for label, s, n, R, zs in (("10 mm, N = 12, R = 29.3 mm, one ring", 0.010, 12, 0.0293, (0.0,)),
                           ("10 mm, N = 12, R = 29.3 mm, pair at +-12.3 mm", 0.010, 12, 0.0293, (0.0123, -0.0123)),
                           ("10 mm, N = 8, R = 23.3 mm, one ring", 0.010, 8, 0.0233, (0.0,)),
                           ("10 mm, N = 8, R = 23.3 mm, pair at +-9.9 mm", 0.010, 8, 0.0233, (0.0099, -0.0099))):
    print(label)
    a = report("Halbach 2 theta", build(s, n, R, zs, "halbach"))
    b = report("largest centre field", build(s, n, R, zs, "max"))
    print("  ratio %.4f" % (b / a))
th = np.linspace(0, 2 * np.pi, 200001)[:-1]
print("dipole ring in its own plane, continuum: mean of 0.5 sqrt(10 + 6 cos 2 theta) over 1.5 = %.4f" % (np.mean(0.5 * np.sqrt(10 + 6 * np.cos(2 * th))) / 1.5))
