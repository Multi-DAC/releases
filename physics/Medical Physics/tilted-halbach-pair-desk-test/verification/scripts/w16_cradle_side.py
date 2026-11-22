"""Which half of a single-ring clamshell should be the cradle, and in which order should the cubes go in?

One ring alone (the upper ring's 12 cubes). side = +1: the cradle is the half nearer the assembly mid-plane and
the cubes come in from outside (+z). side = -1: the cradle is the outer half and the cubes come in from the
mid-plane side (-z); the ring is then loaded upside down. "Lift" is the axial force pulling a cube out of the
open cradle (side * F_z > 0). Orders: W = Wren's (0, 180, 30, 150, 210, 330, 60, 120, 240, 300, 90, 270),
R = the reverse grouping (90, 270, 60, 120, 240, 300, 30, 150, 210, 330, 0, 180).
Printed per stage: the largest force pushing the travelling cube back out on its way in, its lift once seated,
and the largest lift on any seated cube. 10 mm cubes at Br 1.45 T; multiply by 0.729 for 3/8 inch at 1.30 T.
Usage: python -I -B -X utf8 w16_cradle_side.py
"""
import contextlib
import io
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np

with contextlib.redirect_stdout(io.StringIO()):
    import w16_forces as wf

N = wf.N
ds = np.arange(0.0, 0.04001, 0.001)
ORD = {"W": [0, 6, 1, 5, 7, 11, 2, 4, 8, 10, 3, 9], "R": [3, 9, 2, 4, 8, 10, 1, 5, 7, 11, 0, 6]}
F = 0.9525 ** 2 * (1.30 / 1.45) ** 2

for name, blocks in wf.designs.items():
    ring = blocks[:N]
    for side in (1.0, -1.0):
        for oname, order in ORD.items():
            rows = []; w_push = 0.0; w_lift = 0.0
            for k, idx in enumerate(order):
                if k == 0:
                    continue
                present = [ring[j] for j in order[:k]] + [ring[idx]]
                fo = np.array([side * wf.force_torque(k, present, n=16, shift=(0.0, 0.0, side * d))[0][2] for d in ds])
                lifts = [side * wf.force_torque(j, present, n=16)[0][2] for j in range(k + 1)]
                jl = int(np.argmax(lifts))
                rows.append((k + 1, 30 * idx, max(fo.max(), 0.0), fo[0], lifts[jl], 30 * (order[:k] + [idx])[jl]))
                w_push = max(w_push, fo.max()); w_lift = max(w_lift, lifts[jl])
            print("== %s, cradle %s, order %s: worst push on the way in %.2f N, worst lift on a seated cube %.2f N"
                  " (3/8 inch, 1.30 T: %.2f N, %.2f N) =="
                  % (name, "inner half (cubes enter from outside)" if side > 0 else "outer half (cubes enter from the mid-plane side)",
                     oname, w_push, w_lift, w_push * F, w_lift * F))
            print("   " + "  ".join("%d:%d push %.2f seat %+.2f lift %+.2f(%d)" % r for r in rows))
