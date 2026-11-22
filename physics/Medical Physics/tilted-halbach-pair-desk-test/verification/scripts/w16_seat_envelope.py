"""Where the sensor may sit: the field near the centre of both holders at the build scale (3/8 inch cubes, 1.30 T).

Prints the centre fields and ratio, the parity of each field component under x, y, z reflections, the curvature of
B_x along each axis, the ratio's shift per mm^2 and the offset at which it reaches 0.10 %, the on-axis z scan, and
a parabola fit to seven shim heights in S4 with the true centre displaced.
My own run of the numbers in Wren's letter 39 (S2, S3), from my kernel.

Usage: python -I -B -X utf8 w16_seat_envelope.py
"""
import contextlib, io, os, sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
with contextlib.redirect_stdout(io.StringIO()):
    import w16_desk_cmp as dc
import w16_prism as wp

K = 0.9525
SCALE = 1.30 / dc.BR                       # the kernel is linear in Br
BL = {tag: dc.build(0.010 * K, 12, 0.0303 * K, (0.01272 * K, -0.01272 * K), rule)
      for tag, rule in (("S3", "halbach"), ("S4", "max"))}


def B(tag, pts_mm):
    return wp.field_many(np.atleast_2d(np.asarray(pts_mm, float)) * 1e-3, BL[tag]) * SCALE * 1e3     # mT


b0 = {t: B(t, [0, 0, 0])[0] for t in BL}
print("centre field, mT: S3 (%.4f, %.1e, %.1e)  S4 (%.4f, %.1e, %.1e)  ratio %.6f"
      % (*b0["S3"], *b0["S4"], b0["S4"][0] / b0["S3"][0]))

print("-- parity at (1.3, 1.7, 1.1) mm: value, then value under x -> -x, y -> -y, z -> -z, each over the first")
p = np.array([1.3, 1.7, 1.1])
for t in BL:
    v = B(t, p)[0]
    for c, name in enumerate(("Bx", "By", "Bz")):
        signs = []
        for ax in range(3):
            q = p.copy(); q[ax] = -q[ax]
            signs.append(B(t, q)[0][c] / v[c])
        print("   %s %s = %+.5f mT; ratios %+.9f %+.9f %+.9f" % (t, name, v[c], *signs))

print("-- curvature of B_x: (1 / 2 B0) d2Bx/di2 in %/mm^2, from second differences at h = 0.25 mm")
h = 0.25
cur = {}
for t in BL:
    cur[t] = []
    for ax in range(3):
        e = np.zeros(3); e[ax] = h
        d2 = (B(t, e)[0][0] + B(t, -e)[0][0] - 2 * b0[t][0]) / h ** 2
        cur[t].append(100 * d2 / (2 * b0[t][0]))
    print("   %s: x %+.4f  y %+.4f  z %+.4f   sum %+.1e" % (t, *cur[t], sum(cur[t])))
print("   ratio shift, %%/mm^2: x %+.4f  y %+.4f  z %+.4f" % tuple(cur["S4"][i] - cur["S3"][i] for i in range(3)))


def ratio_err(pt):
    return 100 * ((B("S4", pt)[0][0] / B("S3", pt)[0][0]) / (b0["S4"][0] / b0["S3"][0]) - 1)


print("-- offset along one axis at which the ratio (same offset in both holders) is 0.10 % off, exact kernel")
for ax, name in enumerate("xyz"):
    lo, hi = 0.0, 5.0
    for _ in range(50):
        mid = (lo + hi) / 2
        e = np.zeros(3); e[ax] = mid
        lo, hi = (mid, hi) if abs(ratio_err(e)) < 0.10 else (lo, mid)
    print("   |%s| = %.2f mm" % (name, lo))

print("-- ratio error over a box, extreme values on a grid (same offset in both holders)")
for label, (ex, ey, ez) in (("|x|<=1.2, |y|<=0.6, |z|<=0.5", (1.2, 0.6, 0.5)), ("r_xy<=0.75, |z|<=0.75 (square 0.75)", (0.75, 0.75, 0.75)),
                            ("+-2 mm on each axis alone", (None, None, None))):
    if ex is None:
        vals = [ratio_err(np.eye(3)[a] * s) for a in range(3) for s in (-2.0, 2.0)]
    else:
        g = [(x, y, z) for x in np.linspace(-ex, ex, 5) for y in np.linspace(-ey, ey, 5) for z in np.linspace(-ez, ez, 5)]
        if label.startswith("r_xy"):
            g = [q for q in g if q[0] ** 2 + q[1] ** 2 <= 0.75 ** 2 + 1e-12]
        vals = [ratio_err(np.array(q)) for q in g]
    print("   %-38s %+.3f %% to %+.3f %%" % (label, min(vals), max(vals)))

print("-- on-axis z scan, B_x(0,0,z) - B_x(0), mT")
zs = [0.5, 1.0, 2.0, 3.0, 4.0]
for t in BL:
    print("   %s: " % t + "  ".join("z=%.1f %+.3f" % (z, B(t, [0, 0, z])[0][0] - b0[t][0]) for z in zs))

print("-- S4: parabola through seven shim heights 1 mm apart, the true centre displaced by z_true")
rng = np.random.default_rng(1121)
for z_true in (0.0, 0.4, -0.7, 1.3):
    hts = np.arange(-3.0, 3.01, 1.0)
    vals = np.array([B("S4", [0, 0, zz - z_true])[0][0] for zz in hts])
    for noise in (0.0, 0.010):
        est = []
        for _ in range(1 if noise == 0 else 400):
            c = np.polyfit(hts, vals + rng.normal(0, noise, len(hts)), 2)
            est.append(-c[1] / (2 * c[0]))
        est = np.array(est)
        print("   z_true %+.1f mm, noise %.3f mT per height: fitted centre %+.3f mm%s"
              % (z_true, noise, est.mean(), "" if noise == 0 else " (sd %.3f over 400 draws)" % est.std(ddof=1)))

print("-- ambient field of 0.05 mT along x: ratio error, same heading %+.4f %%, one holder turned round %+.4f %%"
      % (100 * ((b0["S4"][0] + 0.05) / (b0["S3"][0] + 0.05) / (b0["S4"][0] / b0["S3"][0]) - 1),
         100 * ((b0["S4"][0] + 0.05) / (b0["S3"][0] - 0.05) / (b0["S4"][0] / b0["S3"][0]) - 1)))
