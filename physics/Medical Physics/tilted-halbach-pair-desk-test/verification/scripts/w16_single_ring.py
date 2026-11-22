"""Loading one ring on its own, in Wren's order, and then mating two loaded rings.

Part 1. One ring alone (the upper ring's 12 cubes, no lower ring). Order: theta = 0, 180, then 30, 150, 210, 330,
then 60, 120, 240, 300, then 90, 270 (anti-aligned cubes first, aligned ones last; Wren's letter 35, F7 item 4).
The new cube travels straight along the bore axis from 40 mm outside its seat to its seat, orientation kept,
the earlier cubes seated. Printed per stage, for the travelling cube: largest axial force pushing it back out,
largest pulling it in, largest sideways force, largest torque. And for the cubes ALREADY seated, once the new
one is in: the largest axial force lifting one out of an open half-nest and the largest torque on one. A cube
weighs about 0.07 N, so anything above that needs a keeper until the lid is on.

Part 2. Both rings loaded and capped. The upper ring is moved along the axis to a gap g above its final place;
printed: the net axial force on it from the lower ring.

Computed for 10 mm cubes at Br = 1.45 T (as w16_forces.py). For 3/8 inch cubes at the same Br multiply forces by
0.9525^2 and torques by 0.9525^3, distances by 0.9525; for another Br multiply both by (Br/1.45)^2.
Usage: python -I -B -X utf8 w16_single_ring.py
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
K = 0.9525
SC = {"10 mm, Br 1.45 T": (1.0, 1.0),
      "3/8 inch, Br 1.45 T": (K ** 2, K ** 3),
      "3/8 inch, Br 1.30 T": (K ** 2 * (1.30 / 1.45) ** 2, K ** 3 * (1.30 / 1.45) ** 2)}
ds = np.arange(0.0, 0.04001, 0.0005)
order = [0, 6, 1, 5, 7, 11, 2, 4, 8, 10, 3, 9]
assert sorted(order) == list(range(N))

for name, blocks in wf.designs.items():
    ring = blocks[:N]
    print("== %s: one ring alone, Wren's order ==" % name)
    print("  stage theta | travelling cube: push out max (at d)  pull in max (at d)   sideways   torque   seated F_out"
          " | seated cubes after: max lift   (theta)   max torque")
    worst = np.zeros(6)
    for k, idx in enumerate(order):
        if k == 0:
            print("  %5d %5d | first cube: no force" % (1, 30 * idx))
            continue
        present = [ring[j] for j in order[:k]] + [ring[idx]]
        me = len(present) - 1
        fo = []; fs = []; tq = []
        for d in ds:
            f, t = wf.force_torque(me, present, n=16, shift=(0.0, 0.0, d))
            fo.append(f[2]); fs.append(np.hypot(f[0], f[1])); tq.append(np.linalg.norm(t))
        fo = np.array(fo); fs = np.array(fs); tq = np.array(tq)
        push = max(fo.max(), 0.0); pull = max(-fo.min(), 0.0)
        lift = []; ttq = []
        for j in range(len(present)):
            f, t = wf.force_torque(j, present, n=16)
            lift.append(f[2]); ttq.append(np.linalg.norm(t))
        jl = int(np.argmax(lift))
        worst = np.maximum(worst, [push, pull, fs.max(), tq.max(), max(lift), max(ttq)])
        print("  %5d %5d | %6.2f N (%5.1f mm)   %6.2f N (%5.1f mm)   %6.2f N  %6.1f N mm  %+6.2f N"
              " | %+6.2f N  (%3d)   %6.1f N mm"
              % (k + 1, 30 * idx, push, ds[fo.argmax()] * 1e3, pull, ds[fo.argmin()] * 1e3, fs.max(), tq.max() * 1e3,
                 fo[0], lift[jl], 30 * (order[:k] + [idx])[jl], max(ttq) * 1e3))
    print("  full single ring, axial force on each cube (theta 0..330): "
          + " ".join("%+.2f" % wf.force_torque(j, ring, n=16)[0][2] for j in range(N)))
    for label, (sf, st) in SC.items():
        print("  worst over the 12 stages, %s: push out %.2f N, pull in %.2f N, sideways %.2f N, torque %.1f N mm;"
              " seated cube lift %.2f N, torque %.1f N mm"
              % (label, worst[0] * sf, worst[1] * sf, worst[2] * sf, worst[3] * st * 1e3, worst[4] * sf,
                 worst[5] * st * 1e3))

    print("== %s: mating two loaded rings. Upper ring a gap g above its final place ==" % name)
    lower = blocks[N:]
    gs = np.arange(0.0, 0.0801, 0.004)
    Fz = []
    for g in gs:
        tot = 0.0
        for j in range(N):
            f, _ = wf.force_torque(N, lower + [ring[j]], n=16, shift=(0.0, 0.0, g))
            tot += f[2]
        Fz.append(tot)
    Fz = np.array(Fz)
    print("   g [mm]: " + " ".join("%6.0f" % (g * 1e3) for g in gs))
    print("   F_z [N]: " + " ".join("%+6.2f" % f for f in Fz) + "   (positive: the rings push apart)")
    work = np.trapezoid(Fz, gs)
    for label, (sf, st) in SC.items():
        print("  %s: largest push apart %.2f N at g = %.0f mm, largest pull together %.2f N; force when mated %+.2f N;"
              " work to bring together from 80 mm %.1f mJ"
              % (label, Fz.max() * sf, gs[Fz.argmax()] * 1e3 * (1.0 if label.startswith("10") else K),
                 max(-Fz.min(), 0.0) * sf, Fz[0] * sf, work * sf * (1.0 if label.startswith("10") else K) * 1e3))
