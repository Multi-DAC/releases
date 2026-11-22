"""Read an Eagle .brd (XML): board outline, holes, and where each part sits, with the pads of one named part.

Usage: python -I -B -X utf8 brd_read.py <file.brd> [part-name-substring]
"""
import math, sys
import xml.etree.ElementTree as ET

root = ET.parse(sys.argv[1]).getroot()
want = sys.argv[2].upper() if len(sys.argv) > 2 else "TMAG"
board = root.find("drawing/board")

xs, ys = [], []
print("-- layer 20 (Dimension) in <plain>:")
for el in board.find("plain"):
    if el.get("layer") == "20":
        a = dict(el.attrib)
        print("  ", el.tag, {k: a[k] for k in a if k != "layer"})
        for k in ("x1", "x2", "x"):
            if k in a: xs.append(float(a[k]))
        for k in ("y1", "y2", "y"):
            if k in a: ys.append(float(a[k]))
if xs:
    print("   outline box: x %.3f to %.3f, y %.3f to %.3f  -> %.3f x %.3f mm" % (min(xs), max(xs), min(ys), max(ys),
                                                                                 max(xs) - min(xs), max(ys) - min(ys)))
print("-- holes in <plain>:")
for el in board.find("plain").iter("hole"):
    print("  ", el.attrib)

packs = {}
for lib in board.iter("library"):
    for p in lib.iter("package"):
        packs[(lib.get("name"), p.get("name"))] = p

print("-- elements:")
for e in board.find("elements"):
    rot = e.get("rot", "R0")
    print("   %-8s %-28s %-22s x %8s y %8s %s" % (e.get("name"), e.get("value", "")[:28], e.get("package")[:22],
                                                 e.get("x"), e.get("y"), rot))
    p = packs.get((e.get("library"), e.get("package")))
    if p is None:
        continue
    for h in p.iter("hole"):
        print("        package hole", h.attrib)
    if want in (e.get("value", "") + e.get("name") + e.get("package")).upper():
        mirror = rot.startswith("M")
        ang = math.radians(float(rot.lstrip("MSR") or 0))
        ex, ey = float(e.get("x")), float(e.get("y"))
        print("     pads of %s (board coordinates):" % e.get("name"))
        for s in p.iter("smd"):
            x, y = float(s.get("x")), float(s.get("y"))
            if mirror: x = -x
            bx = ex + x * math.cos(ang) - y * math.sin(ang); by = ey + x * math.sin(ang) + y * math.cos(ang)
            print("        pad %-3s package (%.3f, %.3f)  board (%.3f, %.3f)  size %s x %s"
                  % (s.get("name"), float(s.get("x")), float(s.get("y")), bx, by, s.get("dx"), s.get("dy")))
        for w in p:
            if w.tag in ("wire", "circle") and w.get("layer") in ("21", "51"):
                print("        outline %s layer %s %s" % (w.tag, w.get("layer"), {k: w.get(k) for k in w.attrib if k not in ("layer", "width")}))
