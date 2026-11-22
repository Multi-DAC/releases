"""Picture of the gauge strip from its STL: a shaded height map seen from above, and one cross-section.

Usage: python render_gauge.py <stl> <png> <x_of_section_mm>
"""
import sys
import numpy as np
import trimesh
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LightSource

stl, png, xsec = sys.argv[1], sys.argv[2], float(sys.argv[3])
m = trimesh.load(stl, process=True)
(x0, y0, z0), (x1, y1, z1) = m.bounds
step = 0.1
xs = np.arange(x0 + step / 2, x1, step)
ys = np.arange(y0 + step / 2, y1, step)
X, Y = np.meshgrid(xs, ys)

# height map by rasterising each upward-facing triangle (the part has no overhangs seen from above)
Hm = np.full(X.shape, -1.0)
tri = m.triangles
nz = m.face_normals[:, 2]
for t in tri[nz > 1e-9]:
    (ax, ay, az), (bx, by, bz), (cx, cy, cz) = t
    i0 = max(int((min(ax, bx, cx) - x0) / step) - 1, 0); i1 = min(int((max(ax, bx, cx) - x0) / step) + 2, len(xs))
    j0 = max(int((min(ay, by, cy) - y0) / step) - 1, 0); j1 = min(int((max(ay, by, cy) - y0) / step) + 2, len(ys))
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
through = Hm < 0            # nothing below: a through hole
Hm[through] = 0.0

fig, (a1, a2) = plt.subplots(2, 1, figsize=(15, 8.6), dpi=120, gridspec_kw={"height_ratios": [1.45, 1]})
ls = LightSource(azdeg=315, altdeg=50)
rgb = ls.shade(Hm, cmap=plt.cm.bone, vert_exag=1.2, dx=step, dy=step, blend_mode="soft", vmin=-2, vmax=z1 + 2)
rgb[through] = (1, 1, 1, 1)
a1.imshow(rgb, origin="lower", extent=(x0, x1, y0, y1))
a1.axvline(xsec, color="tab:red", lw=1, ls="--")
a1.set_title("W16 pocket gauge, seen from above (mm). Top row: flat pockets. Bottom row: nests tilted 34.63 deg. "
             "Red line: the section below.", fontsize=10)
a1.set_xlabel("x"); a1.set_ylabel("y")

sec = m.section(plane_origin=[xsec, 0, 0], plane_normal=[1, 0, 0])
for ent in sec.discrete:
    a2.fill(ent[:, 1], ent[:, 2], facecolor="0.75", edgecolor="k", lw=0.8)
a2.set_aspect("equal")
a2.set_xlim(y0 - 2, y1 + 2); a2.set_ylim(z0 - 1.5, z1 + 6)
a2.set_title("Section at x = %.2f mm (y across, z up): tilted nest on the left, flat pocket on the right. "
             "A cube drops straight in; nothing overhangs." % xsec, fontsize=10)
a2.set_xlabel("y (mm)"); a2.set_ylabel("z (mm)")
fig.tight_layout()
fig.savefig(png)
print("height map %d x %d cells, top %.3f mm, through-hole cells %d" % (Hm.shape[1], Hm.shape[0], Hm.max(), through.sum()))
