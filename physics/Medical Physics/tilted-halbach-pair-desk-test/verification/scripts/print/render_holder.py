"""Picture of one holder ring from its two STL files: both halves seen from above as printed (split face up),
and the closed ring in section through the slots at theta = 0 and 180.

Usage: python render_holder.py <cradle.stl> <lid.stl> <png> <name>
"""
import sys
import numpy as np
import trimesh
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LightSource

cradle_stl, lid_stl, png, name = sys.argv[1:5]
STEP = 0.1


def height_map(m):
    (x0, y0, z0), (x1, y1, z1) = m.bounds
    xs = np.arange(x0 + STEP / 2, x1, STEP); ys = np.arange(y0 + STEP / 2, y1, STEP)
    X, Y = np.meshgrid(xs, ys)
    Hm = np.full(X.shape, -1.0)
    for t in m.triangles[m.face_normals[:, 2] > 1e-9]:
        (ax, ay, az), (bx, by, bz), (cx, cy, cz) = t
        i0 = max(int((min(ax, bx, cx) - x0) / STEP) - 1, 0); i1 = min(int((max(ax, bx, cx) - x0) / STEP) + 2, len(xs))
        j0 = max(int((min(ay, by, cy) - y0) / STEP) - 1, 0); j1 = min(int((max(ay, by, cy) - y0) / STEP) + 2, len(ys))
        if i1 <= i0 or j1 <= j0:
            continue
        px = X[j0:j1, i0:i1]; py = Y[j0:j1, i0:i1]
        d = (by - cy) * (ax - cx) + (cx - bx) * (ay - cy)
        if abs(d) < 1e-14:
            continue
        l1 = ((by - cy) * (px - cx) + (cx - bx) * (py - cy)) / d
        l2 = ((cy - ay) * (px - cx) + (ax - cx) * (py - cy)) / d
        l3 = 1 - l1 - l2
        inside = (l1 >= -1e-9) & (l2 >= -1e-9) & (l3 >= -1e-9)
        z = l1 * az + l2 * bz + l3 * cz
        sub = Hm[j0:j1, i0:i1]
        sub[inside] = np.maximum(sub[inside], z[inside])
    through = Hm < 0
    Hm[through] = 0.0
    return Hm, through, (x0, x1, y0, y1), z1


cr = trimesh.load(cradle_stl, process=True)
li = trimesh.load(lid_stl, process=True)
fig = plt.figure(figsize=(15, 10.5), dpi=120)
gs = fig.add_gridspec(2, 2, height_ratios=[2.6, 1])
ls = LightSource(azdeg=315, altdeg=50)
for k, (m, label) in enumerate(((cr, "%s cradle (outer half), split face up. Arrow: the cube's north pole points along it. "
                                     "Number: loading order." % name),
                                (li, "%s lid (inner half), split face up." % name))):
    Hm, through, ext, ztop = height_map(m)
    rgb = ls.shade(Hm, cmap=plt.cm.bone, vert_exag=1.2, dx=STEP, dy=STEP, blend_mode="soft", vmin=-2, vmax=ztop + 2)
    rgb[through] = (1, 1, 1, 1)
    a = fig.add_subplot(gs[0, k])
    a.imshow(rgb, origin="lower", extent=ext)
    a.axhline(0.0, color="tab:red", lw=1, ls="--")
    a.set_title(label, fontsize=9); a.set_xlabel("x (mm)"); a.set_ylabel("y (mm)")
    print("%s: height map %d x %d cells, top %.3f mm" % (label[:12], Hm.shape[1], Hm.shape[0], Hm.max()))

a2 = fig.add_subplot(gs[1, :])
PIN = float(sys.argv[5]) if len(sys.argv) > 5 else 2.5      # the cradle's pins stand this far above its split face
ztop = cr.bounds[1][2] - PIN + li.bounds[1][2]
for m, col, flip in ((li, "0.6", False), (cr, "0.82", True)):
    sec = m.section(plane_origin=[0, 0, 0], plane_normal=[0, 1, 0])
    for ent in sec.discrete:
        z = ztop - ent[:, 2] if flip else ent[:, 2]
        a2.fill(ent[:, 0], z, facecolor=col, edgecolor="k", lw=0.8)
a2.set_aspect("equal"); a2.set_xlim(-50, 50); a2.set_ylim(-1.5, ztop + 1.5)
a2.set_title("%s ring closed, section along the red line (x across, z up; the holder's mid-plane is z = 0). "
             "Dark: lid. Light: cradle. The white squares are the magnet pockets at theta = 180 and 0." % name, fontsize=9)
a2.set_xlabel("x (mm)"); a2.set_ylabel("z (mm)")
fig.tight_layout()
fig.savefig(png)
