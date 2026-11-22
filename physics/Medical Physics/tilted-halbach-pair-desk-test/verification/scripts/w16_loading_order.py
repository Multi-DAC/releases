"""Forces while the holder is being filled: every stage of a loading order, not only the last cube.

Order A: the upper ring, pocket by pocket round the ring (theta = 0, 30, ... 330 deg), then the lower ring the same
way. Order B: opposite pockets in turn (0, 180, 30, 210, ...), upper ring then lower ring.
At stage k the cubes already loaded are seated and the new cube travels along the bore axis from 40 mm outside its
seat (away from the mid-plane) to its seat, orientation kept. Printed for each stage: the largest axial force
pushing it back out, the largest pulling it in, the largest sideways force, the largest torque, and the axial
force once seated. Geometry and method as in w16_forces.py (face quadrature n = 16, steps of 0.5 mm).
Usage: python -I -B -X utf8 w16_loading_order.py
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
ds = np.arange(0.0, 0.04001, 0.0005)
orders = {"A, round the ring": list(range(N)) + list(range(N, 2 * N)),
          "B, opposite pockets in turn": [j for i in range(N // 2) for j in (i, i + N // 2)]
                                         + [N + j for i in range(N // 2) for j in (i, i + N // 2)]}

for name, blocks in wf.designs.items():
    for oname, order in orders.items():
        print("== %s, order %s ==" % (name, oname))
        print("  stage  ring   theta   push out max (at d)     pull in max (at d)     sideways max   torque max    seated F_out")
        worst = np.zeros(4)
        for k, idx in enumerate(order):
            if k == 0:
                print("  %5d  %s  %5d   first cube: no force" % (1, "upper" if idx < N else "lower", 30 * (idx % N)))
                continue
            present = [blocks[j] for j in order[:k]] + [blocks[idx]]
            me = len(present) - 1
            out = 1.0 if idx < N else -1.0                   # direction away from the mid-plane
            fo = []; fs = []; tq = []
            for d in ds:
                f, t = wf.force_torque(me, present, n=16, shift=(0.0, 0.0, out * d))
                fo.append(out * f[2]); fs.append(np.hypot(f[0], f[1])); tq.append(np.linalg.norm(t))
            fo = np.array(fo); fs = np.array(fs); tq = np.array(tq)
            push = max(fo.max(), 0.0); pull = max(-fo.min(), 0.0)
            worst = np.maximum(worst, [push, pull, fs.max(), tq.max()])
            print("  %5d  %s  %5d   %6.2f N (%5.1f mm)     %6.2f N (%5.1f mm)     %6.2f N     %6.1f N mm   %+6.2f N"
                  % (k + 1, "upper" if idx < N else "lower", 30 * (idx % N), push, ds[fo.argmax()] * 1e3,
                     pull, ds[fo.argmin()] * 1e3, fs.max(), tq.max() * 1e3, fo[0]))
        print("  worst over the 24 stages: push out %.2f N, pull in %.2f N, sideways %.2f N, torque %.1f N mm"
              % (worst[0], worst[1], worst[2], worst[3] * 1e3))
