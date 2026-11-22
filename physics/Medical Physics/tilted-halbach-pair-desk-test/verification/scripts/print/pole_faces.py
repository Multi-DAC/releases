"""Which faces of each cube show from above when a ring's cradle (outer half) lies split face up.

The cradle of the upper ring holds z >= z0; turned split face up, the viewer looks along +z, so a face shows when
its outward normal has n_z < 0. For each slot: every showing face, its slope from level (acos|n_z|), whether it is
a pole face (+e3 north, -e3 south), and whether it lies on the head or tail side of the printed arrow (the arrow is
the in-plane part of the magnetisation). Then the one-line rule is tested on every slot:
  the steepest showing face is a pole face; it is NORTH if on the arrow-head side, SOUTH if on the tail side.

Usage: python -I -B -X utf8 pole_faces.py
"""
import math
import numpy as np
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import holder as H

NAMES = ["+e1", "-e1", "+e2", "-e2", "+e3 NORTH", "-e3 SOUTH"]
for tag, rule in (("S3", "halbach"), ("S4", "max")):
    print("== %s cradle, split face up ==" % tag)
    bad = 0
    for k, (r0, Q) in enumerate(H.poses(rule)):
        th = (math.degrees(math.atan2(r0[1], r0[0])) + 360) % 360
        m = Q[:, 2]
        mxy = m[:2] / np.linalg.norm(m[:2])
        faces = []
        for j in range(3):
            for sgn in (+1, -1):
                n = sgn * Q[:, j]
                if n[2] < -1e-9:
                    side = float(n[:2] @ mxy)
                    faces.append((math.degrees(math.acos(min(1.0, abs(n[2])))), NAMES[2 * j + (0 if sgn > 0 else 1)],
                                  "head" if side > 1e-9 else ("tail" if side < -1e-9 else "level")))
        faces.sort(reverse=True)
        tilt = math.degrees(math.asin(abs(m[2])))
        row = "; ".join("%s at %.2f deg, %s side" % (nm, sl, sd) for sl, nm, sd in faces)
        if len(faces) == 1:
            # level cube: both pole faces are vertical; north is toward the arrow head by construction
            verdict = "level cube, pole faces vertical"
        else:
            sl, nm, sd = faces[0]
            ok = ("NORTH" in nm and sd == "head") or ("SOUTH" in nm and sd == "tail")
            gap = faces[0][0] - faces[1][0]
            verdict = "rule %s, steep minus flat %.2f deg" % ("HOLDS" if ok else "FAILS", gap)
            bad += 0 if ok else 1
        print("  slot %2d theta %5.1f order %2d  m tilt %5.2f deg | %s | %s"
              % (k, th, H.ORDER.index(k) + 1, tilt, row, verdict))
    print("  slots where the rule fails: %d" % bad)
