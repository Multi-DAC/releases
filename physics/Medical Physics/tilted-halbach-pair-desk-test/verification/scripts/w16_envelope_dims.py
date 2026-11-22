"""Holder envelope for the two desk holders: where the cubes' corners reach. For each design: the smallest and
largest distance of any cube corner from the bore axis, the axial extent, the closest approach of the two rings,
and the smallest gap between cubes (w16_wren_checks.gap). Geometry as in w16_forces.py.
Usage: python -I -B -X utf8 w16_envelope_dims.py
"""
import contextlib
import io
import itertools
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np

with contextlib.redirect_stdout(io.StringIO()):
    import w16_desk_cmp as dc
    import w16_wren_checks as wc

S, N, R, Z = 0.010, 12, 0.0303, 0.01272
corn = np.array(list(itertools.product((-0.5, 0.5), repeat=3))) * S

for name, rule in (("S3 plain Halbach", "halbach"), ("S4 largest centre field", "max")):
    blocks = dc.build(S, N, R, (Z, -Z), rule)
    pts = np.array([np.asarray(r0) + corn @ Q.T for r0, _, Q, _ in blocks])       # (24, 8, 3)
    rho = np.hypot(pts[..., 0], pts[..., 1])
    up = pts[:N]
    print(name)
    print("  corners: nearest to the axis %.3f mm, farthest %.3f mm" % (rho.min() * 1e3, rho.max() * 1e3))
    print("  upper ring spans z from %.3f to %.3f mm; clear space between the rings %.3f mm"
          % (up[..., 2].min() * 1e3, up[..., 2].max() * 1e3, 2 * up[..., 2].min() * 1e3))
    print("  smallest gap between cubes: %s" % (wc.gap(blocks),))
    for i in range(N // 4 + 1):
        r0, _, Q, _ = blocks[i]
        e3 = Q[:, 2]
        print("    theta %3d: magnetisation (%+.4f, %+.4f, %+.4f), out of plane %+.2f deg, in-plane angle %+.2f deg;"
              " cube spans z %.2f to %.2f mm, radius %.2f to %.2f mm"
              % (30 * i, *e3, np.degrees(np.arcsin(e3[2])), np.degrees(np.arctan2(e3[1], e3[0])),
                 pts[i, :, 2].min() * 1e3, pts[i, :, 2].max() * 1e3, rho[i].min() * 1e3, rho[i].max() * 1e3))
