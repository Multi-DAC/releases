"""How the pocket clearance moves the ratio. A cube in a pocket with clearance is pushed to one side by the force
on it and turned by the torque on it, and the two holders push and turn their cubes differently. Each seated cube
is displaced by `delta` along the in-pocket part of its net force (all three components; the lid and floor are
treated like the walls), and separately turned by `alpha` about its torque axis. The centre field of each holder
and the ratio S4/S3 are recomputed. Also: every cube moved outward or inward radially by delta (a holder printed
slightly large or small), one holder only and both; and the geometry scaled to 3/8 inch cubes.
Geometry and force code as in w16_forces.py.
Usage: python -I -B -X utf8 w16_seating.py
"""
import contextlib
import io
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import w16_prism as wp

with contextlib.redirect_stdout(io.StringIO()):
    import w16_forces as wf
    import w16_desk_cmp as dc


def rot(axis, ang):
    a = axis / np.linalg.norm(axis)
    K = np.array([[0, -a[2], a[1]], [a[2], 0, -a[0]], [-a[1], a[0], 0]])
    return np.eye(3) + np.sin(ang) * K + (1 - np.cos(ang)) * K @ K


def b0(blocks):
    return wp.field_many(np.zeros((1, 3)), blocks)[0, 0]


FT = {k: [wf.force_torque(i, v, n=16) for i in range(len(v))] for k, v in wf.designs.items()}
nom = {k: b0(v) for k, v in wf.designs.items()}
k3, k4 = "S3 plain Halbach", "S4 largest centre field"
print("nominal: S3 %.3f mT, S4 %.3f mT, ratio %.4f" % (nom[k3] * 1e3, nom[k4] * 1e3, nom[k4] / nom[k3]))


def pushed(k, delta):
    return [(np.asarray(r0) + delta * f / np.linalg.norm(f), dims, Q, Br)
            for (r0, dims, Q, Br), (f, _) in zip(wf.designs[k], FT[k])]


def turned(k, alpha):
    return [(r0, dims, rot(t, alpha) @ Q, Br) for (r0, dims, Q, Br), (_, t) in zip(wf.designs[k], FT[k])]


def radial(k, delta):
    out = []
    for r0, dims, Q, Br in wf.designs[k]:
        r0 = np.asarray(r0, float).copy()
        r0[:2] *= 1 + delta / np.hypot(r0[0], r0[1])
        out.append((r0, dims, Q, Br))
    return out


print("== every cube pushed by delta along the force on it ==")
for delta in (0.05e-3, 0.1e-3, 0.2e-3):
    a, b = b0(pushed(k3, delta)), b0(pushed(k4, delta))
    print("  delta %.2f mm: S3 %+.3f%%, S4 %+.3f%%, ratio %.4f (%+.3f%%)"
          % (delta * 1e3, (a / nom[k3] - 1) * 100, (b / nom[k4] - 1) * 100, b / a, (b / a / (nom[k4] / nom[k3]) - 1) * 100))
print("== every cube turned by alpha about the torque on it ==")
for deg in (0.5, 1.0, 2.0):
    al = np.radians(deg)
    a, b = b0(turned(k3, al)), b0(turned(k4, al))
    print("  alpha %.1f deg: S3 %+.3f%%, S4 %+.3f%%, ratio %.4f (%+.3f%%)"
          % (deg, (a / nom[k3] - 1) * 100, (b / nom[k4] - 1) * 100, b / a, (b / a / (nom[k4] / nom[k3]) - 1) * 100))
print("== every cube moved radially by delta (the holder printed large or small) ==")
for delta in (0.1e-3, 0.2e-3, -0.1e-3, -0.2e-3):
    a, b = b0(radial(k3, delta)), b0(radial(k4, delta))
    print("  delta %+.1f mm: S3 alone %+.3f%%, S4 alone %+.3f%%; ratio if S4 only %+.3f%%, if S3 only %+.3f%%, if both %+.3f%%"
          % (delta * 1e3, (a / nom[k3] - 1) * 100, (b / nom[k4] - 1) * 100, (b / nom[k4] - 1) * 100,
             (nom[k3] / a - 1) * 100, (b / a / (nom[k4] / nom[k3]) - 1) * 100))
print("== random errors, 400 draws: each cube's Br, position and direction drawn independently ==")
rng = np.random.default_rng(20261121)


def draw(k, sbr, spos, sang):
    out = []
    for r0, dims, Q, Br in wf.designs[k]:
        ax = rng.normal(size=3)
        out.append((np.asarray(r0) + rng.normal(scale=spos, size=3), dims, rot(ax, rng.normal(scale=sang)) @ Q,
                    Br * (1 + rng.normal(scale=sbr))))
    return out


for sbr, spos, sdeg in ((0.01, 0.0, 0.0), (0.02, 0.0, 0.0), (0.0, 0.1e-3, 0.0), (0.0, 0.0, 1.0), (0.02, 0.1e-3, 1.0)):
    r = np.array([b0(draw(k4, sbr, spos, np.radians(sdeg))) / b0(draw(k3, sbr, spos, np.radians(sdeg))) for _ in range(400)])
    print("  sd Br %.0f%%, sd position %.2f mm per axis, sd direction %.1f deg: ratio mean %.4f, sd %.3f%%  (48 different cubes)"
          % (sbr * 100, spos * 1e3, sdeg, r.mean(), r.std() / r.mean() * 100))
print("  same 24 cubes moved from S3 to S4, pocket for pocket, Br spread only:")
for sbr in (0.01, 0.02):
    rr = []
    for _ in range(400):
        f = 1 + rng.normal(scale=sbr, size=24)
        A = [(r0, d, Q, Br * fi) for (r0, d, Q, Br), fi in zip(wf.designs[k3], f)]
        B = [(r0, d, Q, Br * fi) for (r0, d, Q, Br), fi in zip(wf.designs[k4], f)]
        rr.append(b0(B) / b0(A))
    rr = np.array(rr)
    print("    sd Br %.0f%%: ratio mean %.4f, sd %.3f%%" % (sbr * 100, rr.mean(), rr.std() / rr.mean() * 100))
print("== the same design with 3/8 inch cubes: every length times 0.9525 ==")
sc = 0.9525
for k, rule in ((k3, "halbach"), (k4, "max")):
    bl = dc.build(wf.S * sc, wf.N, wf.R * sc, (wf.Z * sc, -wf.Z * sc), rule)
    f = max(np.linalg.norm(wf.force_torque(i, bl, n=16)[0]) for i in range(wf.N))
    print("  %s: edge %.3f mm, R %.3f mm, rings at +-%.3f mm: B0 %.3f mT at Br 1.45 T, %.3f mT at Br 1.30 T; largest force on a cube %.3f N at Br 1.45 T"
          % (k, wf.S * sc * 1e3, wf.R * sc * 1e3, wf.Z * sc * 1e3, b0(bl) * 1e3, b0(bl) * 1e3 * 1.30 / 1.45, f))
