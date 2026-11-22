"""W16 sensor seat: a plug for the holder's 39 mm bore that carries the Adafruit 6490 board, plus height shims.

The plug hangs in the bore by a flange resting on the holder's outer face (or on shims laid on that face). The
board stands on edge in a well, long axis along the bore, face normal along the holder's x (the field), cable up.
The well is set off along x so the sensor element, not the board, is on the bore axis. Unshimmed, the element is
DEEP below the mid-plane; shims of 0.6, 1, 2 and 3 mm raise it through the centre, so a height scan brackets the
peak of the tilted holder's field (Wren, letter 39). Two tabs half a turn apart end in the same half-round notch
as the holder's rim mark at theta = 0: a length of 1.75 mm filament laid in the rim groove registers either tab.

Coordinates here: the plug's contact face (underside of the flange) is z = 0, the body hangs toward -z, +x is
tab A. In the holder the contact face is at z = Z_TOP + (shim stack).

NOT SOURCED: board thickness BT = 1.6 and element height ELEM_H = 0.7 above the board (Wren's nominals). They set
only the x offset, the forgiving axis. Component height 3.0 = product page's 4.6 overall less 1.6.

Usage: python -I -B -X utf8 sensor_seat.py <out_dir>
"""
import itertools, math, os, sys
import numpy as np
from manifold3d import Manifold, CrossSection

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import holder as H

Z_TOP = H.Z_TOP                    # holder's outer face above the mid-plane
BORE_R = H.R_IN
GAP = 0.25                         # radial clearance of the two guiding lands in the bore
RELIEF = 0.40                      # the neck and the waist are this much smaller in radius again
FL_R, FL_T = 30.0, 3.0             # flange
TAB_W, TAB_R, NOTCH_R = 8.0, H.R_OUT, 1.0
A_HOLE = (38.0, 2.5)               # tab A only: a hole at this radius, of this diameter. Component side faces tab A.
BL, BW, BT, CORNER = 25.40, 17.78, 1.60, 2.54     # board: length (along the bore), width, thickness, corner radius
ELEM_H = 0.70                      # element above the board's component face
COMP_H, COMP_HALF_W = 3.0, 5.0     # component envelope used in the fit check
SLOT_T, SLOT_W = BT + 0.20, BW + 0.30
LIP = 1.0                          # the well's lips overlap each long edge of the board by LIP - 0.15
ROOM_X = (-1.0, 3.2)               # component room in front of the board (x range), between the lips
DEEP = 3.0                         # element below the mid-plane with no shim
RELIEF_Z, FLOOR = 1.0, 2.0         # component room is this much deeper than the board's seat; plastic under it
ROD_D = H.ROD_D                    # push-out hole through the floor, same rod as the holder's
NECK, WAIST = 8.0, (19.5, 28.5)    # depths below the contact face: neck 0..NECK, waist between these
CONE = RELIEF                      # 45 degree transitions
TIP = 1.0                          # lead-in chamfer at the tip
SHIMS = (0.6, 1.0, 2.0, 3.0)
SHIM_RI, SHIM_RO = 19.9, FL_R
FIL_R = 0.875                      # 1.75 mm filament
HARD_R = H.HEX_AF / math.sqrt(3)   # nut pocket's corner radius on the holder's outer face

X_MID = -(BT / 2 + ELEM_H)         # board's mid-plane, so that the element is at x = 0


def box(x0, x1, y0, y1, z0, z1):
    return Manifold.cube([x1 - x0, y1 - y0, z1 - z0]).translate([x0, y0, z0])


def body_length(deep):
    return Z_TOP + deep + BL / 2 + RELIEF_Z + FLOOR


def build(plug_r=BORE_R - GAP, fl_r=FL_R, deep=DEEP, well=True, notch_b_deg=180.0):
    L = body_length(deep)
    rn = plug_r - RELIEF
    # (depth from, depth to, radius at from, radius at to)
    prof = [(0.0, NECK, rn, rn), (NECK, NECK + CONE, rn, plug_r), (NECK + CONE, WAIST[0] - CONE, plug_r, plug_r),
            (WAIST[0] - CONE, WAIST[0], plug_r, rn), (WAIST[0], WAIST[1], rn, rn), (WAIST[1], WAIST[1] + CONE, rn, plug_r),
            (WAIST[1] + CONE, L - TIP, plug_r, plug_r), (L - TIP, L, plug_r, plug_r - TIP)]
    m = Manifold()
    for d0, d1, r0, r1 in prof:           # cylinder(height, radius_low, radius_high): low is the deeper end
        m += Manifold.cylinder(d1 - d0, r1, r0).translate([0, 0, -d1])
    m += Manifold.cylinder(FL_T, fl_r)
    m += box(-TAB_R, TAB_R, -TAB_W / 2, TAB_W / 2, 0, FL_T)
    for deg in (0.0, notch_b_deg):
        a = math.radians(deg)
        m -= Manifold.cylinder(FL_T + 2, NOTCH_R).translate([TAB_R * math.cos(a), TAB_R * math.sin(a), -1])
    if well:
        seat_z = -(Z_TOP + deep + BL / 2)
        m -= box(X_MID - SLOT_T / 2, X_MID + SLOT_T / 2, -SLOT_W / 2, SLOT_W / 2, seat_z, FL_T + 1)
        m -= box(ROOM_X[0], ROOM_X[1], -SLOT_W / 2 + LIP, SLOT_W / 2 - LIP, seat_z - RELIEF_Z, FL_T + 1)
        m -= Manifold.cylinder(L + 1 + seat_z, ROD_D / 2).translate([X_MID, 0, -L - 1])
        m -= Manifold.cylinder(FL_T + 2, A_HOLE[1] / 2).translate([A_HOLE[0], 0, -1])
    return m


def shim(t):
    cs = CrossSection.circle(SHIM_RO) - CrossSection.circle(SHIM_RI)
    m = cs.extrude(t)
    k = {0.6: 0, 1.0: 1, 2.0: 2, 3.0: 3}[t]
    for i in range(k):                    # number of rim notches = thickness in mm; the 0.6 mm shim has none
        a = math.radians(90 + 14 * (i - (k - 1) / 2))
        m -= Manifold.cylinder(t + 2, 1.5).translate([SHIM_RO * math.cos(a), SHIM_RO * math.sin(a), -1])
    return m


def board(deep=DEEP, dt=0.0, dw=0.0):
    """The board and its component envelope, in plug coordinates, resting on the well's seat."""
    w, t = BW + dw, BT + dt
    r = CORNER
    pts = []
    for cx, cy, a0 in ((BL - r, w / 2 - r, 0), (r, w / 2 - r, 90), (r, -w / 2 + r, 180), (BL - r, -w / 2 + r, 270)):
        for k in range(13):
            a = math.radians(a0 + 90 * k / 12)
            pts.append((cx + r * math.cos(a), cy + r * math.sin(a)))
    plate = CrossSection([np.array(pts)]).extrude(t)                 # x' = along the board, y' = across, z' = thickness
    comp = box(0.3, BL - 0.3, -COMP_HALF_W, COMP_HALF_W, BT, BT + COMP_H)   # on the component face (z' = BT)
    seat_z = -(Z_TOP + deep + BL / 2)
    out = []
    for part in (plate, comp):
        # a quarter turn about y: (x', y', z') -> (z', y', -x'), so the component face looks toward +x (tab A)
        p = part.rotate([0, 90, 0])
        out.append(p.translate([X_MID - BT / 2 - (dt / 2 if part is plate else 0), 0, seat_z + BL]))
    return out[0], out[1]


def verts(m):
    return np.asarray(m.to_mesh().vert_properties)[:, :3].astype(float)


def export(m, path):
    import trimesh
    me = m.to_mesh()
    V = np.asarray(me.vert_properties)[:, :3].astype(float); F = np.asarray(me.tri_verts).astype(int)
    trimesh.Trimesh(V, F, process=False).export(path)
    back = trimesh.load(path, process=True)
    return back, V.max(0) - V.min(0)


def as_printed(m):
    """Flange on the bed, body up: turn over about x, bed at z = 0."""
    p = m.rotate([180, 0, 0])
    return p.translate([0, 0, -p.bounding_box()[2]])


_HOLDERS = {}


def holder_solid(rule):
    if rule not in _HOLDERS:
        outer, inner = H.build(rule, 9.60)[:2]
        ring = outer + inner
        _HOLDERS[rule] = ring + ring.rotate([180, 0, 0])
    return _HOLDERS[rule]


def hardware():
    m = Manifold()
    for a in H.SCREW_ANG + H.RING_ANG:
        m += H.at_angle(Manifold.cylinder(40, HARD_R).translate([0, 0, Z_TOP - 20]), H.SCREW_R, a)
    return m


def in_holder(m, t, turn=0.0):
    return m.rotate([0, 0, turn]).translate([0, 0, Z_TOP + t])


def element_from_mesh(m, deep_for_len):
    """Seat height and slot centre read off the built solid, not the constants."""
    L = body_length(deep_for_len)
    col = box(X_MID - 0.1, X_MID + 0.1, 3.0, 3.2, -L - 1, FL_T + 1) ^ m           # plastic under the board's seat
    seat_z = verts(col)[:, 2].max()
    probe = box(-10, 10, 8.5, 8.6, -30.0, -29.9)                                  # in a lip: only the slot is open
    gap = verts(probe - m)
    x_mid = (gap[:, 0].min() + gap[:, 0].max()) / 2
    return seat_z, x_mid, gap[:, 0].max() - gap[:, 0].min()


def band_radius(m, z_lo, z_hi):
    v = verts(m ^ box(-60, 60, -60, 60, z_lo, z_hi))
    return float(np.linalg.norm(v[:, :2], axis=1).max())


def main(out_dir):
    os.makedirs(out_dir, exist_ok=True)
    fails = []

    def check(name, ok, text):
        print("%s  %s: %s" % ("PASS" if ok else "FAIL", name, text))
        if not ok:
            fails.append(name)

    seat = build()
    L = body_length(DEEP)
    print("plug radius %.2f in a %.1f bore (radial gap %.2f), neck and waist %.2f; body %.3f long; flange radius %.1f,"
          " %.1f thick; tabs %.0f wide to R = %.0f" % (BORE_R - GAP, 2 * BORE_R, GAP, BORE_R - GAP - RELIEF, L, FL_R, FL_T,
                                                    TAB_W, TAB_R))
    print("well: board slot %.2f x %.2f, centred at x = %.2f; component room x %.1f to %.1f between lips of %.1f;"
          " NOT SOURCED: board thickness %.1f, element height %.1f" % (SLOT_T, SLOT_W, X_MID, ROOM_X[0], ROOM_X[1], LIP,
                                                                      BT, ELEM_H))

    # A. the solid and its print file
    pr = as_printed(seat)
    back, bbox = export(pr, os.path.join(out_dir, "w16-sensor-seat.stl"))
    n = back.face_normals; cz = back.triangles_center[:, 2]
    flat = float(back.area_faces[(n[:, 2] < -0.9) & (cz > 0.05)].sum())
    slope = float(back.area_faces[(n[:, 2] < -0.5) & (n[:, 2] >= -0.9) & (cz > 0.05)].sum())
    check("A1 seat solid", seat.status().name == "NoError" and seat.genus() == 2 and back.is_watertight
          and back.is_winding_consistent,
          "status %s, genus %d (well through the floor hole, and tab A's hole), volume %.1f cm^3, as printed %.1f x %.1f"
          " x %.2f mm, STL watertight %s" % (seat.status().name, seat.genus(), seat.volume() / 1000, bbox[0], bbox[1],
                                             bbox[2], back.is_watertight))
    print("      as printed (flange on the bed): flat downward area off the bed %.1f mm^2 (the well's end, bridged"
          " across %.1f and %.1f mm); 45 degree downward area %.1f mm^2 (the two cones)"
          % (flat, SLOT_T, ROOM_X[1] - ROOM_X[0], slope))

    # B. in the holder
    hard = hardware()
    for tag, rule in (("S3", "halbach"), ("S4", "max")):
        hs = holder_solid(rule)
        for t in (0.0, sum(SHIMS)):
            s = in_holder(seat, t)
            below = s ^ box(-60, 60, -60, 60, -60, Z_TOP + t - 0.01)
            ov = (s ^ hs).volume()
            gap = below.min_gap(hs, 5.0)
            tip = verts(s)[:, 2].min() + Z_TOP
            check("B1 %s shim %.1f" % (tag, t), ov < 1e-9 and abs(gap - GAP) < 1e-3 and tip > 1.0,
                  "overlap with the holder %.2e mm^3; plug body to the bore %.3f mm; tip %.3f mm above the desk"
                  % (ov, gap, tip))
    top = in_holder(seat, 0.0) ^ box(-60, 60, -60, 60, Z_TOP - 0.01, 60)
    g = top.min_gap(hard, 30.0)
    check("B2 hardware", g > 2.0, "flange and tabs to the nearest nut pocket (corner radius %.2f at R = %.0f): %.3f mm"
          % (HARD_R, H.SCREW_R, g))
    fat = in_holder(build(plug_r=BORE_R + 0.05), 0.0)
    check("B3 mutant: plug 0.05 mm over the bore radius", (fat ^ holder_solid("max")).volume() > 1.0,
          "overlap %.1f mm^3 (the check can fail)" % (fat ^ holder_solid("max")).volume())
    wide = in_holder(build(fl_r=40.0), 0.0) ^ box(-60, 60, -60, 60, Z_TOP - 0.01, 60)
    check("B4 mutant: flange radius 40", (wide ^ hard).volume() > 1.0,
          "overlap with the nut pockets' columns %.1f mm^3" % (wide ^ hard).volume())

    # C. the board in the well
    plate, comp = board()
    ov = (plate ^ seat).volume() + (comp ^ seat).volume()
    check("C1 board fits", ov < 1e-9, "board %.2f x %.2f x %.2f with a %.1f mm component block over |y| <= %.1f:"
          " overlap with the plug %.2e mm^3" % (BL, BW, BT, COMP_H, COMP_HALF_W, ov))
    def ovl(dx=0.0, dy=0.0, dz=0.0):
        return (plate.translate([dx, dy, dz]) ^ seat).volume()
    check("C2 clearances", ovl(dx=0.09) < 1e-9 and ovl(dx=-0.09) < 1e-9 and ovl(dx=0.11) > 1e-3 and ovl(dx=-0.11) > 1e-3
          and ovl(dy=0.14) < 1e-9 and ovl(dy=0.16) > 1e-3 and ovl(dz=-0.02) > 1e-3,
          "moved 0.09 along x: %.1e, 0.11: %.3f mm^3; 0.14 along y: %.1e, 0.16: %.3f; 0.02 down: %.3f (it rests on the seat)"
          % (ovl(dx=0.09), ovl(dx=0.11), ovl(dy=0.14), ovl(dy=0.16), ovl(dz=-0.02)))
    thick, _ = board(dt=0.3); wider, _ = board(dw=0.4)
    check("C3 mutants: board 0.3 thicker, 0.4 wider", (thick ^ seat).volume() > 1e-3 and (wider ^ seat).volume() > 1e-3,
          "overlap %.2f and %.2f mm^3" % ((thick ^ seat).volume(), (wider ^ seat).volume()))
    bt = verts(plate)[:, 2].max()
    print("      board's top edge %.2f mm below the flange's top face; lips overlap each long edge by %.2f mm"
          % (FL_T - bt, LIP - (SLOT_W - BW) / 2))

    # D. where the element sits, read off the solid
    seat_z, x_mid, slot_t = element_from_mesh(seat, DEEP)
    ex = x_mid + BT / 2 + ELEM_H
    ez0 = Z_TOP + seat_z + BL / 2
    check("D1 element on the axis", abs(ex) < 1e-6 and abs(slot_t - SLOT_T) < 1e-6,
          "slot %.3f wide, centred at x = %.4f; element at x = %.4f, y = 0 by symmetry (for the nominal 1.6 and 0.7)"
          % (slot_t, x_mid, ex))
    stacks = {}
    for k in range(len(SHIMS) + 1):
        for c in itertools.combinations(SHIMS, k):
            stacks.setdefault(round(sum(c), 3), c)
    hs = sorted(stacks)
    zs = [ez0 + h for h in hs]
    print("      shim stack (mm) -> element height above the mid-plane (mm):")
    print("      " + "  ".join("%.1f -> %+.1f" % (h, z) for h, z in zip(hs, zs)))
    step = max(b - a for a, b in zip(hs, hs[1:]))
    check("D2 scan brackets the centre", abs(ez0 + DEEP) < 1e-6 and min(zs) <= -2.999 and max(zs) >= 2.999
          and all(any(abs(h - w) < 1e-9 for h in hs) for w in range(7)),
          "unshimmed %+.3f; highest %+.3f; whole millimetres 0 to 6 all present; largest step %.1f mm, so the chosen"
          " stack is within %.2f mm of any fitted centre" % (ez0, max(zs), step, step / 2))
    mz, _, _ = element_from_mesh(build(deep=0.0), 0.0)
    mut0 = Z_TOP + mz + BL / 2
    check("D3 mutant: seat with no depth", not (mut0 <= -2.999),
          "unshimmed element at %+.3f: shims could only raise it, and D2's lowest-height test fails on it" % mut0)

    # E. the half turn
    ext = build(well=False)
    d = (ext - ext.rotate([0, 0, 180])).volume() + (ext.rotate([0, 0, 180]) - ext).volume()
    check("E1 half turn maps the outside onto itself", d < 1e-6, "body, flange, tabs and notches: difference %.2e mm^3" % d)
    bad = build(well=False, notch_b_deg=181.0)
    d2 = (bad - bad.rotate([0, 0, 180])).volume() + (bad.rotate([0, 0, 180]) - bad).volume()
    check("E2 mutant: tab B's notch 1 degree off", d2 > 0.1, "difference %.2f mm^3" % d2)

    # F. registration on the rim groove with a length of filament
    hs4 = holder_solid("max")
    fil = Manifold.cylinder(60, FIL_R).translate([TAB_R, 0, -25])
    v = [(fil ^ hs4).volume()] + [(fil ^ in_holder(seat, t, turn)).volume() for t in (0.0, 3.0, 6.6) for turn in (0.0, 180.0)]
    check("F1 filament in the groove", max(v) < 1e-9, "1.75 mm filament in the rim groove: overlap with the holder and"
          " with the plug at shims 0, 3, 6.6, either tab: largest %.2e mm^3" % max(v))
    slack = 2 * (NOTCH_R - FIL_R)
    f_off = fil.translate([0, NOTCH_R - FIL_R - 0.005, 0])
    ok_in = (f_off ^ in_holder(seat, 0.0, math.degrees(slack / TAB_R) - 0.02)).volume()
    out1 = (fil ^ in_holder(seat, 0.0, 1.0)).volume()
    check("F2 registration", ok_in < 1e-9 and out1 > 1e-3,
          "slack %.2f mm at R = %.0f, +-%.2f degrees: free at %.2f degrees with the filament at the groove's side"
          " (%.1e mm^3), blocked at 1 degree (%.3f mm^3)" % (slack, TAB_R, math.degrees(slack / TAB_R),
                                                           math.degrees(slack / TAB_R) - 0.02, ok_in, out1))

    # G. shims
    for t in SHIMS:
        sh = shim(t)
        back, bb = export(sh, os.path.join(out_dir, "w16-seat-shim-%.1fmm.stl" % t))
        check("G1 shim %.1f" % t, sh.status().name == "NoError" and back.is_watertight and abs(bb[2] - t) < 1e-6,
              "thickness %.2f, bore %.1f (%.2f clear of the lands), outer diameter %.0f, rim notches %d, volume %.2f cm^3,"
              " watertight %s" % (bb[2], 2 * SHIM_RI, SHIM_RI - (BORE_R - GAP), 2 * SHIM_RO,
                                  {0.6: 0, 1.0: 1, 2.0: 2, 3.0: 3}[t], sh.volume() / 1000, back.is_watertight))
    stack = Manifold(); z = Z_TOP
    for t in SHIMS:
        stack += shim(t).translate([0, 0, z]); z += t
    s = in_holder(seat, sum(SHIMS))
    check("G2 full stack", (stack ^ s).volume() < 1e-9 and (stack ^ hs4).volume() < 1e-9 and (stack ^ hard).volume() < 1e-9,
          "all four shims under the flange: overlap with the plug %.1e, the holder %.1e, the nut pockets' columns %.1e mm^3"
          % ((stack ^ s).volume(), (stack ^ hs4).volume(), (stack ^ hard).volume()))

    # H. the reliefs sit where a printed bore is tightest (first layers on the bed: outer faces and mid-plane faces)
    ok = True; rows = []
    for t in (0.0, 3.0, sum(SHIMS)):
        s = in_holder(seat, t) ^ box(-60, 60, -60, 60, -60, Z_TOP + t - 0.01)
        r_top = band_radius(s, Z_TOP - 0.6, Z_TOP + 0.0)
        r_mid = band_radius(s, -0.6, 0.6)
        up = band_radius(s, 6.0, 10.0); lo = band_radius(s, -12.0, -9.0)
        ok &= abs(r_top - (BORE_R - GAP - RELIEF)) < 1e-6 and abs(r_mid - (BORE_R - GAP - RELIEF)) < 1e-6
        ok &= abs(up - (BORE_R - GAP)) < 1e-6 and abs(lo - (BORE_R - GAP)) < 1e-6
        rows.append("shim %.1f: radius %.2f at the outer face, %.2f at the mid-plane; %.2f in the upper ring, %.2f in the lower"
                    % (t, r_top, r_mid, up, lo))
    check("H1 reliefs and lands", ok, "; ".join(rows))

    print("%d checks failed%s" % (len(fails), (": " + ", ".join(fails)) if fails else ""))
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1]))
