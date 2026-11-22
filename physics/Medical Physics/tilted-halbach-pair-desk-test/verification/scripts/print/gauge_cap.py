"""Cap strip for row B (the tilted nests) of the W16 pocket gauge strip.

A cube turned about the strip's long axis is carried onto itself by a half turn about that axis, and the half turn
swaps "swept up" for "swept down". So the lid half of each tilted nest is the nest itself turned over about the line
(y = Y_TILT, z = H): the cap is the row-B slice of the gauge strip, y from 0 to 2 * Y_TILT, printed the same way up
(pockets up, nothing overhangs) and turned over onto the gauge for the test. Column i of the cap closes column i of
the gauge. One corner is cut off at the small-size end; when the cap is on, both cut corners are at the same end.

Usage: python -I -B -X utf8 gauge_cap.py <out_dir>
"""
import os, sys
import numpy as np
from manifold3d import Manifold, CrossSection

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gauge_strip as G

WC = 2 * G.Y_TILT                  # cap width: the slice the half turn maps onto itself


def build():
    body = G.box(0, G.WX, 0, WC, 0, G.H)
    c = 5.0
    wedge = CrossSection([np.array([[-1, WC - c], [c, WC + 1], [-1, WC + 1]], float)]).extrude(G.H + 2).translate([0, 0, -1])
    cut = wedge
    for i, s in enumerate(G.SIZES):
        cut += G.tilted_nest(G.col_x(i), G.Y_TILT, s)
    return body - cut


def turned_over(m):
    """The cap as it sits on the gauge: half turn about the line y = Y_TILT, z = H."""
    return m.translate([0, -G.Y_TILT, -G.H]).rotate([180, 0, 0]).translate([0, G.Y_TILT, G.H])


def main(out_dir):
    import trimesh
    os.makedirs(out_dir, exist_ok=True)
    cap = build()
    gauge = G.build()
    on = turned_over(cap)
    me = cap.to_mesh()
    V = np.asarray(me.vert_properties)[:, :3].astype(float); F = np.asarray(me.tri_verts).astype(int)
    stl = os.path.join(out_dir, "w16-gauge-cap.stl")
    trimesh.Trimesh(V, F, process=False).export(stl)
    back = trimesh.load(stl, process=True)
    print("manifold status:", cap.status(), " genus:", cap.genus(), " triangles:", len(F))
    print("volume: %.2f cm^3   bounding box: %.3f x %.3f x %.3f mm" % ((cap.volume() / 1000,) + tuple(V.max(0) - V.min(0))))
    print("STL read back: watertight %s, winding consistent %s, volume %.2f cm^3" % (
        back.is_watertight, back.is_winding_consistent, back.volume / 1000))
    n = back.face_normals; cz = back.triangles_center[:, 2]
    print("downward-facing area off the bed: %.3f mm^2" % float(back.area_faces[(n[:, 2] < -0.5) & (cz > 0.05)].sum()))
    lo, hi = on.bounding_box()[2], on.bounding_box()[5]
    print("cap on the gauge: occupies z %.3f to %.3f (gauge top face at %.3f); overlap with the gauge %.2e mm^3"
          % (lo, hi, G.H, (on ^ gauge).volume()))
    print("size   cube of that size outside the closed nest   in plastic   next size up in plastic   walls (y) cap")
    both = gauge + on
    for i, s in enumerate(G.SIZES):
        cx = G.col_x(i)
        cube = Manifold.cube([s, s, s], center=True).rotate([G.TILT, 0, 0]).translate([cx, G.Y_TILT, G.H])
        big = Manifold.cube([s + 0.05] * 3, center=True).rotate([G.TILT, 0, 0]).translate([cx, G.Y_TILT, G.H])
        slab = G.box(cx - G.PITCH / 2, cx + G.PITCH / 2, -5, 50, -5, 30)
        closed = slab - (both ^ slab)                      # empty space in this column
        outside = (cube - closed).volume()
        vy = np.asarray((cube).to_mesh().vert_properties)[:, 1]
        print("%.2f   %.2e mm^3                              %.2e     %.3f mm^3                %.3f / %.3f mm"
              % (s, outside, (cube ^ both).volume(), (big ^ both).volume(), vy.min(), WC - vy.max()))


if __name__ == "__main__":
    main(sys.argv[1])
