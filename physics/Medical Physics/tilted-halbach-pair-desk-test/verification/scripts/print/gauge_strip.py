"""W16 pocket gauge strip for a 3/8 inch (9.525 mm) cube magnet, printed flat with the pockets facing up.

Row A: ten flat half-depth square pockets, 9.45 to 9.90 mm in 0.05 mm steps (the S3 half-nest).
Row B: ten tilted half-nests of the same sizes, the cube turned 34.63 degrees about the strip's long axis and
       swept straight up, so a cube drops in vertically and nothing overhangs (the S4 half-nest).
Each pocket has a 3 mm push-out hole through the floor. The size is debossed above each column. One corner of
the strip is cut off at the small-size end.

Usage: python gauge_strip.py <out_dir>
"""
import sys, os, math, json
import numpy as np
from manifold3d import Manifold, CrossSection, FillRule, set_circular_segments
from matplotlib.textpath import TextPath
from matplotlib.font_manager import FontProperties

set_circular_segments(48)

NOM = 9.525                       # nominal cube edge, mm
SIZES = [round(9.45 + 0.05 * i, 2) for i in range(10)]
TILT = 34.63                      # degrees, the largest tilt in S4
PITCH = 14.5
H = 9.0                           # strip height; the top face stands in for the clamshell split plane
MARGIN_X = 3.5
Y_TEXT, Y_FLAT, Y_TILT = 36.5, 26.0, 9.5      # row centre lines
WY = 43.0
HOLE_R = 1.5
RELIEF_R = 0.5                    # corner relief in the flat pockets
TEXT_H, TEXT_DEPTH = 4.2, 0.6
WX = 2 * MARGIN_X + PITCH * len(SIZES)


def col_x(i):
    return MARGIN_X + PITCH * (i + 0.5)


def box(x0, x1, y0, y1, z0, z1):
    return Manifold.cube([x1 - x0, y1 - y0, z1 - z0]).translate([x0, y0, z0])


def flat_pocket(cx, cy, s):
    d = NOM / 2
    p = box(cx - s / 2, cx + s / 2, cy - s / 2, cy + s / 2, H - d, H + 1)
    for sx in (-1, 1):
        for sy in (-1, 1):
            p += Manifold.cylinder(d + 1, RELIEF_R).translate([cx + sx * s / 2, cy + sy * s / 2, H - d])
    p += Manifold.cylinder(H + 2, HOLE_R).translate([cx, cy, -1])
    return p


def tilted_nest(cx, cy, s):
    t = math.radians(TILT)
    a = s / 2
    cube = Manifold.cube([s, s, s], center=True).rotate([TILT, 0, 0])
    swept = Manifold.batch_hull([cube, cube.translate([0, 0, 12])]).translate([cx, cy, H])
    y_low = cy - a * (math.cos(t) - math.sin(t))        # the lowest edge of the tilted cube
    return swept + Manifold.cylinder(H + 2, HOLE_R).translate([cx, y_low, -1])


def text3d(s, cx, cy, height, depth):
    tp = TextPath((0, 0), s, size=height, prop=FontProperties(family="DejaVu Sans", weight="bold"))
    polys = [np.asarray(p[:-1], float) for p in tp.to_polygons() if len(p) > 3]
    cs = CrossSection(polys, FillRule.EvenOdd)
    x0, y0, x1, y1 = cs.bounds()
    cs = cs.translate([cx - (x0 + x1) / 2, cy - (y0 + y1) / 2])
    return cs.extrude(depth + 1).translate([0, 0, H - depth])


def build():
    body = box(0, WX, 0, WY, 0, H)
    # cut-off corner at the small-size end, on the label side
    c = 6.0
    wedge = CrossSection([np.array([[-1, WY - c], [c, WY + 1], [-1, WY + 1]], float)]).extrude(H + 2).translate([0, 0, -1])
    cut = wedge
    for i, s in enumerate(SIZES):
        cut += flat_pocket(col_x(i), Y_FLAT, s)
        cut += tilted_nest(col_x(i), Y_TILT, s)
        cut += text3d("%.2f" % s, col_x(i), Y_TEXT, TEXT_H, TEXT_DEPTH)
    return body - cut


def gap(m, x0, x1, y0, y1, z0, z1, along):
    """Length of the empty stretch along a thin probe box."""
    probe = box(x0, x1, y0, y1, z0, z1)
    area = {"x": (y1 - y0) * (z1 - z0), "y": (x1 - x0) * (z1 - z0)}[along]
    return (probe.volume() - (probe ^ m).volume()) / area


def main(out_dir):
    os.makedirs(out_dir, exist_ok=True)
    m = build()
    mesh = m.to_mesh()
    V = np.asarray(mesh.vert_properties)[:, :3].astype(float)
    F = np.asarray(mesh.tri_verts).astype(int)
    import trimesh
    tm = trimesh.Trimesh(V, F, process=False)
    stl = os.path.join(out_dir, "w16-gauge-strip.stl")
    tm.export(stl)
    back = trimesh.load(stl, process=True)
    print("manifold status:", m.status(), " genus:", m.genus(), " triangles:", len(F))
    print("volume: %.2f cm^3   bounding box: %.3f x %.3f x %.3f mm" % ((m.volume() / 1000,) + tuple(V.max(0) - V.min(0))))
    print("STL read back: watertight %s, winding consistent %s, volume %.2f cm^3" % (
        back.is_watertight, back.is_winding_consistent, back.volume / 1000))
    t = math.radians(TILT)
    print("size   flat x-width  flat y-width  flat depth   nest x-width  nest top y-width (expected)  floor under nest")
    e = 0.1
    rows = []
    for i, s in enumerate(SIZES):
        cx = col_x(i)
        off = 2.6      # probe away from the push-out hole and the corner reliefs
        fx = gap(m, cx - 7, cx + 7, Y_FLAT + off - e, Y_FLAT + off + e, H - 2 - e, H - 2 + e, "x")
        fy = gap(m, cx + off - e, cx + off + e, Y_FLAT - 7, Y_FLAT + 7, H - 2 - e, H - 2 + e, "y")
        # depth: a thin vertical probe from the top down, off the hole
        pr = box(cx + off - e, cx + off + e, Y_FLAT + off - e, Y_FLAT + off + e, 0, H)
        fd = (pr.volume() - (pr ^ m).volume()) / (4 * e * e)
        nx = gap(m, cx - 7, cx + 7, Y_TILT - e, Y_TILT + e, H - 0.5 - e, H - 0.5 + e, "x")
        ny = gap(m, cx + off - e, cx + off + e, Y_TILT - 9, Y_TILT + 9, H - 2 * e, H, "y")
        a = s / 2
        ny_exp = a * (math.cos(t) + math.sin(t) + 1 / math.cos(t)) - e * math.tan(t)   # mean over the 0.2 mm probe
        floor = H - a * (math.sin(t) + math.cos(t))
        rows.append((s, fx, fy, fd, nx, ny, ny_exp, floor))
        print("%.2f   %9.4f     %9.4f     %8.4f     %9.4f     %9.4f (%.4f)          %.3f" % rows[-1])
    worst = max(max(abs(r[1] - r[0]), abs(r[2] - r[0]), abs(r[4] - r[0])) for r in rows)
    print("largest |width - nominal size| over all probes: %.2e mm" % worst)
    print("flat pocket depth expected %.4f mm (half the nominal edge)" % (NOM / 2))
    # the nest top width is the mean over the top 0.2 mm; the leaning wall makes it e*tan(tilt) under the top-face value
    json.dump({"sizes": SIZES, "pitch": PITCH, "bbox": list(map(float, V.max(0) - V.min(0)))},
              open(os.path.join(out_dir, "w16-gauge-strip.json"), "w"))
    render(V, F, os.path.join(out_dir, "w16-gauge-strip.png"))


def render(V, F, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    tri = V[F]
    n = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    n /= np.linalg.norm(n, axis=1)[:, None] + 1e-30
    light = np.array([0.35, -0.5, 0.8]); light /= np.linalg.norm(light)
    shade = 0.35 + 0.65 * np.clip(n @ light, 0, 1)
    cols = np.stack([0.55 * shade, 0.68 * shade, 0.85 * shade, np.ones_like(shade)], axis=1)
    fig = plt.figure(figsize=(14, 9), dpi=130)
    for k, (el, az) in enumerate(((62, -80), (28, -60))):
        ax = fig.add_subplot(2, 1, k + 1, projection="3d")
        ax.add_collection3d(Poly3DCollection(tri, facecolors=cols, edgecolors="none"))
        ax.set_xlim(0, V[:, 0].max()); ax.set_ylim(0, V[:, 1].max()); ax.set_zlim(0, V[:, 2].max())
        ax.set_box_aspect(V.max(0) - V.min(0))
        ax.view_init(elev=el, azim=az)
        ax.set_axis_off()
    fig.subplots_adjust(0, 0, 1, 1, 0, 0)
    fig.savefig(path)


if __name__ == "__main__":
    main(sys.argv[1])
