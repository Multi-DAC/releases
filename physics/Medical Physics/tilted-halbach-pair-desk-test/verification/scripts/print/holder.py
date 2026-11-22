"""W16 magnet holders as print files: each ring a two-half clamshell, for 3/8 inch cubes.

One holder = two rings. A ring = an OUTER half (the cradle, cubes are loaded into it) and an INNER half (the lid,
it faces the other ring across the mid-plane). The lower ring's parts are the upper ring's parts turned over
(rotation by 180 degrees about the x axis maps one ring onto the other, arrows included), so each holder needs two
STL files, each printed twice.

Geometry (mm), 0.9525 x the 10 mm design: cube centres at radius R = 28.861, rings at z = +-12.116, cube
orientation from w16_desk_cmp.build(): e3 along the magnetisation, e1 = z x e3, e2 = e3 x e1.
Split at the ring's own mid-plane z0. In the outer half the pocket is the cube swept straight toward the split
plane (-z) and cut at it; in the inner half the cube swept +z and cut. So a cube drops straight into the cradle,
the lid comes straight down, and both halves print split face up with no overhang (Wren's letter 36).

Usage: python -I -B -X utf8 holder.py <out_dir> [pocket_size_mm]
"""
import contextlib, io, math, os, sys
import numpy as np
from manifold3d import Manifold, CrossSection, FillRule, set_circular_segments
from matplotlib.textpath import TextPath
from matplotlib.font_manager import FontProperties

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
with contextlib.redirect_stdout(io.StringIO()):
    import w16_desk_cmp as dc

set_circular_segments(96)
K = 0.9525
NOM = 10.0 * K                     # 9.525 nominal cube edge
N = 12
R = 30.3 * K
Z0 = 12.72 * K
R_IN, R_OUT = 19.5, 47.0           # bore radius (39 mm bore: 1.5 mm of wall at the S4 cube corners), flange radius
CAP = 2.0
Z_TOP = 19.676 * K + CAP           # outer face of the outer half (S4's highest cube corner plus the cap)
SCREW_R, SCREW_D = 42.0, 3.4
SCREW_ANG = (45, 135, 225, 315)    # STACK screws through both rings; symmetric under theta -> -theta
RING_ANG = (15, 165, 195, 345)     # RING screws closing one ring (lid to cradle); also symmetric
HEX_AF, HEX_DEPTH = 5.8, 3.0       # captive M3 nut pockets on the outer face
CB_D, CB_DEPTH = 6.5, 4.4          # head counterbore for the ring screws, on the lid's mid-plane face
PIN_ANG, PIN_D, PIN_H = (105, 255), 3.0, 2.5      # cradle-to-lid pins on the cradle's split face
ROD_D = 3.4                        # hole through the lid's floor under every pocket, for a 3 mm rod (Wren, letter 38):
                                   # rods hold the cubes down while the lid goes on, and push cubes out afterwards.
                                   # No ring-to-ring dowels (v1 had two): a dowel longer than its holes would hold the
                                   # rings apart, and the axial gap is the error that matters (letter 38, Q5).
ARROW_R, ARROW_DEPTH = 41.5, 0.6
ORDER = [0, 6, 1, 5, 7, 11, 2, 4, 8, 10, 3, 9]    # Wren's loading order, by slot index (theta = 30 deg x index)
NUM_H, NUM_OFF = 3.0, 6.2          # loading-order numbers: height, tangential offset from the arrow


def poses(rule):
    blocks = dc.build(0.010 * K, N, 0.0303 * K, (0.01272 * K,), rule)
    return [(np.asarray(r0) * 1e3, np.asarray(Q)) for r0, _, Q, _ in blocks]


def cube_at(r0, Q, s):
    A = np.zeros((3, 4)); A[:, :3] = Q; A[:, 3] = r0
    return Manifold.cube([s, s, s], center=True).transform(A)


def swept(c, dz):
    return Manifold.batch_hull([c, c.translate([0, 0, dz])])


def annulus(z_lo, z_hi):
    cs = CrossSection.circle(R_OUT) - CrossSection.circle(R_IN)
    return cs.extrude(z_hi - z_lo).translate([0, 0, z_lo])


def at_angle(m, r, deg):
    a = math.radians(deg)
    return m.translate([r * math.cos(a), r * math.sin(a), 0])


def arrow(r0, Q, z_face, into_plus_z):
    m = Q[:, 2]
    u = np.array([m[0], m[1]]); u /= np.linalg.norm(u)
    v = np.array([-u[1], u[0]])
    th = math.atan2(r0[1], r0[0])
    c = ARROW_R * np.array([math.cos(th), math.sin(th)])
    pts = [(-2.5, -0.6), (0.5, -0.6), (0.5, -1.5), (2.5, 0.0), (0.5, 1.5), (0.5, 0.6), (-2.5, 0.6)]
    poly = np.array([c + a * u + b * v for a, b in pts])
    ex = CrossSection([poly]).extrude(ARROW_DEPTH + 1.0)
    return ex.translate([0, 0, z_face - 1.0 if into_plus_z else z_face - ARROW_DEPTH])


def number(label, th):
    """Loading-order number on the cradle's split face, readable when that face is up (the model seen from -z)."""
    tp = TextPath((0, 0), label, size=NUM_H, prop=FontProperties(family="DejaVu Sans", weight="bold"))
    polys = [np.asarray(p[:-1], float) for p in tp.to_polygons() if len(p) > 3]
    lo = np.min(np.vstack(polys), axis=0); hi = np.max(np.vstack(polys), axis=0)
    mid = (lo + hi) / 2
    # printed frame (x', y') = model (x, -y). There the slot sits at angle -th; the text goes NUM_OFF further round.
    a = -th
    er = np.array([math.cos(a), math.sin(a)]); et = np.array([-math.sin(a), math.cos(a)])
    c = ARROW_R * er + NUM_OFF * et
    out = []
    for p in polys:
        q = p - mid
        w = c + np.outer(q[:, 0], et) + np.outer(q[:, 1], -er)      # baseline tangential, tops toward the axis
        w[:, 1] *= -1.0                                              # back to model coordinates
        out.append(w)
    return CrossSection(out, FillRule.EvenOdd).extrude(ARROW_DEPTH + 1.0).translate([0, 0, Z0 - 1.0])


def build(rule, s):
    P = poses(rule)
    outer = annulus(Z0, Z_TOP)
    inner = annulus(0.0, Z0)
    cav_o = Manifold(); cav_i = Manifold()
    for r0, Q in P:
        c = cube_at(r0, Q, s)
        cav_o += swept(c, -30.0)
        cav_i += swept(c, +30.0)
    hex_pocket = Manifold.cylinder(HEX_DEPTH + 1, HEX_AF / math.sqrt(3), circular_segments=6).translate([0, 0, Z_TOP - HEX_DEPTH])
    hard_o = Manifold(); hard_i = Manifold()          # every cut that is not a magnet pocket, per half
    for a in SCREW_ANG + RING_ANG:
        hole = at_angle(Manifold.cylinder(60, SCREW_D / 2).translate([0, 0, -5]), SCREW_R, a)
        hard_o += hole + at_angle(hex_pocket.rotate([0, 0, a]), SCREW_R, a)
        hard_i += hole
    for a in RING_ANG:
        hard_i += at_angle(Manifold.cylinder(CB_DEPTH + 1, CB_D / 2).translate([0, 0, -1]), SCREW_R, a)
    for a in PIN_ANG:      # holes in the lid for the cradle's pins
        hard_i += at_angle(Manifold.cylinder(PIN_H + 0.6, PIN_D / 2 + 0.15).translate([0, 0, Z0 - PIN_H - 0.6]), SCREW_R, a)
    rods = [Manifold.cylinder(Z0 + 1.0, ROD_D / 2).translate([r0[0], r0[1], -1.0]) for r0, Q in P]
    marks = Manifold()
    for k, (r0, Q) in enumerate(P):   # arrows on the cradle's split face: the cube's north face goes toward the head
        marks += arrow(r0, Q, Z0, into_plus_z=True)
        o = ORDER.index(k) + 1
        marks += number(("%d." if o in (6, 9) else "%d") % o, math.atan2(r0[1], r0[0]))
    notch = Manifold.cylinder(60, 1.0).translate([R_OUT, 0, -5])      # rim mark at theta = 0, both halves
    outer = outer - cav_o - hard_o - marks - notch
    inner = inner - cav_i - hard_i - notch
    for rod in rods:
        inner -= rod
    for a in PIN_ANG:      # pins on the cradle (outer half) split face, pointing toward the lid
        outer += at_angle(Manifold.cylinder(PIN_H, PIN_D / 2).translate([0, 0, Z0 - PIN_H]), SCREW_R, a)
    return outer, inner, P, cav_o, cav_i, hard_o + marks, hard_i, rods


def as_printed(m, split_up_is_minus_z):
    if split_up_is_minus_z:
        m = m.rotate([180, 0, 0])
    lo = m.bounding_box()[2]
    return m.translate([0, 0, -lo])


def export(m, path):
    import trimesh
    me = m.to_mesh()
    V = np.asarray(me.vert_properties)[:, :3].astype(float); F = np.asarray(me.tri_verts).astype(int)
    trimesh.Trimesh(V, F, process=False).export(path)
    back = trimesh.load(path, process=True)
    return back.is_watertight and back.is_winding_consistent, V.max(0) - V.min(0), back


def overhang_area(mesh, split_face_z):
    """Area of faces pointing downward (normal z < -0.5, steeper than 30 deg from vertical) that are not on the bed."""
    n = mesh.face_normals; c = mesh.triangles_center
    bad = (n[:, 2] < -0.5) & (c[:, 2] > 0.05)
    return float(mesh.area_faces[bad].sum())


def main(out_dir, s):
    os.makedirs(out_dir, exist_ok=True)
    print("pocket size %.3f mm for a nominal %.3f mm cube; R = %.3f, z0 = %.3f, bore %.1f, flange diameter %.1f,"
          " ring height %.3f, holder height %.3f" % (s, NOM, R, Z0, 2 * R_IN, 2 * R_OUT, Z_TOP, 2 * Z_TOP))
    for tag, rule in (("S3", "halbach"), ("S4", "max")):
        outer, inner, P, cav_o, cav_i, hard_o, hard_i, rods = build(rule, s)
        print("== %s ==" % tag)
        # 0. the lower ring is the upper ring turned over about the x axis: cubes as solids, and magnetisations
        lower = dc.build(0.010 * K, N, 0.0303 * K, (-0.01272 * K,), rule)
        up = Manifold(); lo = Manifold(); dm = 0.0
        for k, (r0, Q) in enumerate(P):
            up += cube_at(r0, Q, NOM)
            rl, _, Ql, _ = lower[(N - k) % N]
            lo += cube_at(np.asarray(rl) * 1e3, np.asarray(Ql), NOM)
            dm = max(dm, float(np.abs(np.asarray(Ql)[:, 2] - np.array([1, -1, -1]) * Q[:, 2]).max()))
        turned = up.rotate([180, 0, 0])
        print("  lower ring against the upper ring turned over: solids differ by %.2e mm^3, magnetisation by %.1e"
              % ((turned - lo).volume() + (lo - turned).volume(), dm))
        # 1. each nominal cube lies inside the two cavities together, and the cubes do not touch the bodies
        left = 0.0; clash = 0.0
        cubes = Manifold()
        for r0, Q in P:
            c = cube_at(r0, Q, NOM)
            cubes += c
            left += (c - (cav_o + cav_i)).volume()
        clash = (cubes ^ outer).volume() + (cubes ^ inner).volume()
        print("  nominal cubes outside the pockets: %.2e mm^3; cubes overlapping plastic: %.2e mm^3" % (left, clash))
        # 2. lid comes straight off, cradle releases cubes straight out
        w = max(((cubes.translate([0, 0, -t]) ^ outer).volume() for t in (0.5, 1, 2, 4, 8, 16)))
        w2 = max(((cubes.translate([0, 0, t]) ^ inner).volume() for t in (0.5, 1, 2, 4, 8, 16)))
        print("  cubes moved out of the cradle along its axis: largest overlap %.2e mm^3; lid lifted off: %.2e mm^3" % (w, w2))
        # 3. walls
        walls = []
        for half, body in (("outer", outer), ("inner", inner)):
            cs = [swept(cube_at(r0, Q, s), -30.0 if half == "outer" else 30.0) for r0, Q in P]
            zlo, zhi = (Z0, Z_TOP) if half == "outer" else (0.0, Z0)
            slab = Manifold.cube([200, 200, zhi - zlo]).translate([-100, -100, zlo])
            cs = [c ^ slab for c in cs]
            gaps = [cs[i].min_gap(cs[(i + 1) % N], 20.0) for i in range(N)]
            rad = [np.linalg.norm(np.asarray(c.to_mesh().vert_properties)[:, :2], axis=1) for c in cs]
            walls.append((half, min(gaps), min(r.min() for r in rad) - R_IN, R_OUT - max(r.max() for r in rad)))
            print("  %s half: thinnest web between pockets %.3f mm, thinnest wall to the bore %.3f mm, to the rim %.3f mm"
                  % walls[-1])
            hard = (hard_o if half == "outer" else hard_i) ^ slab
            allc = Manifold()
            for c in cs:
                allc += c
            hr = np.linalg.norm(np.asarray(hard.to_mesh().vert_properties)[:, :2], axis=1)
            print("  %s half: screw holes, nut pockets, pin holes and marks: nearest to a magnet pocket %.3f mm,"
                  " nearest to the rim %.3f mm" % (half, allc.min_gap(hard, 20.0), R_OUT - hr.max()))
        # 3b. rod holes: each lies wholly in the vertical shadow of its own pocket, so it removes floor plastic only
        stray = 0.0; plug = []
        for (r0, Q), rod in zip(P, rods):
            c = cube_at(r0, Q, s)
            shadow = Manifold.batch_hull([c.translate([0, 0, -40.0]), c.translate([0, 0, 40.0])])
            stray += (rod - shadow).volume()
            body = rod ^ Manifold.cube([200, 200, Z0]).translate([-100, -100, 0.0])
            plug.append((body - swept(c, 30.0)).volume() / (math.pi * ROD_D ** 2 / 4))
        print("  rod holes (%.1f mm, one under each pocket of the lid): volume outside the pocket's shadow %.2e mm^3;"
              " mean floor drilled through %.3f to %.3f mm" % (ROD_D, stray, min(plug), max(plug)))
        cap_o = Z_TOP - max(np.asarray(cube_at(r0, Q, s).to_mesh().vert_properties)[:, 2].max() for r0, Q in P)
        floor_i = min(np.asarray(cube_at(r0, Q, s).to_mesh().vert_properties)[:, 2].min() for r0, Q in P)
        print("  plastic over the highest cube corner (outer face) %.3f mm; under the lowest corner (mid-plane side) %.3f mm"
              % (cap_o, floor_i))
        # 4. export as printed: split face up
        for half, body, flip in (("outer-cradle", outer, True), ("inner-lid", inner, False)):
            m = as_printed(body, flip)
            path = os.path.join(out_dir, "w16-%s-%s.stl" % (tag, half))
            ok, bbox, mesh = export(m, path)
            print("  %s: status %s, genus %d, volume %.1f cm^3, box %.1f x %.1f x %.2f mm, STL watertight %s,"
                  " downward-facing area off the bed %.2f mm^2"
                  % (os.path.basename(path), m.status(), m.genus(), m.volume() / 1000, bbox[0], bbox[1], bbox[2], ok,
                     overhang_area(mesh, None)))


if __name__ == "__main__":
    main(sys.argv[1], float(sys.argv[2]) if len(sys.argv) > 2 else 9.60)
