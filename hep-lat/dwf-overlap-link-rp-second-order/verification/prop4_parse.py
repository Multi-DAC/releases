"""Summarise out-prop4-large-t*.txt (clock reading 14:42): worst cases, growth in t, saturation in L5."""
import re
import sys
from collections import defaultdict

TS = [0, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000]
rows = []
for path in sys.argv[1:]:
    for line in open(path, encoding="utf8"):
        if not line.startswith("m=") or "FAILED" in line:
            continue
        g = re.match(r"m=([\d.]+) a5=([\d.]+) k1=([\d.]+) L5=\s*(\d+) S=\s*([\d.]+): max ([\d.e+-]+) at t=\s*(\d+) ;"
                     r" top DW/r\* ([\d.]+) PV/r\* ([\d.]+) nsurf (\d+) ; r\(t\) = (.*)$", line.strip())
        if not g:
            continue
        m, a5, k1, L5, S, mx, ta, tdw, tpv, ns, rs = g.groups()
        r = [float(x) for x in rs.split(";")[0].split()]   # prop4_rerun.py adds "; spread ..." (edit between clock readings 14:47 and 14:52)
        rows.append(dict(m=float(m), a5=float(a5), k1=float(k1), L5=int(L5), S=float(S), max=float(mx), targ=int(ta),
                         tdw=float(tdw), tpv=float(tpv), ns=int(ns), r=r))
print("rows:", len(rows))
w = max(rows, key=lambda d: d["max"])
print("worst sup_t: %.4f at m=%.1f a5=%.4f k1=%.3f L5=%d (t=%d)" % (w["max"], w["m"], w["a5"], w["k1"], w["L5"], w["targ"]))
grow = sorted(rows, key=lambda d: -d["max"] / max(d["r"][0], 1e-300))
print("largest sup_t / r(0):")
for d in grow[:12]:
    print("  %.3f  m=%.1f a5=%.4f k1=%.3f L5=%3d S=%6.1f argmax t=%5d top DW/r* %.6f PV/r* %.6f r(10^4)=%.3g"
          % (d["max"] / d["r"][0], d["m"], d["a5"], d["k1"], d["L5"], d["S"], d["targ"], d["tdw"], d["tpv"], d["r"][-1]))
print("cases with argmax t > 0:", sum(1 for d in rows if d["targ"] > 0), "of", len(rows))
by = defaultdict(list)
for d in rows:
    by[(d["m"], d["a5"], d["k1"])].append(d)
print("sup_t by L5 (saturation), for k1 = 0 and k1 = pi/16:")
for key in sorted(by):
    if key[2] > 0.2:
        continue
    ds = sorted(by[key], key=lambda d: d["L5"])
    print("  m=%.1f a5=%.4f k1=%.3f: %s" % (key + ("; ".join("L5=%d %.4f@t=%d" % (d["L5"], d["max"], d["targ"]) for d in ds),)))
