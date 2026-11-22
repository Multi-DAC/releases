"""W16 host side, v2: collect sensor lines from the board over USB serial, and work out the ratio S4/S3.

One SEATING is the sensor seat resting in a holder at mark A, or half a turn round at mark B.
One REPEAT is an A seating and a B seating of the same holder, the seat taken right out in between.
The board sends lines on two ranges in turn (+-266 mT and +-133 mT); a seating collects n lines of each.

  python w16_read.py ports
  python w16_read.py collect --port COM5 --out run.csv --holder S3 --rep 1 --seating A [--n 16] [--shim 3.0]
  python w16_read.py analyse run.csv [--predict 1.1562]
  python w16_read.py zfit scan.csv            (S4, seating A only, one --shim value per repeat number)

Arithmetic (analyse), separately for each range. Per seating, the mean reading vector over its samples. Per repeat,
  h = (A - B) / 2      the field across the bore in the sensor's own axes, sensor offset removed
  s = (A + B) / 2      the sensor's offset (plus any field along the bore, which a half turn does not reverse)
  m = |h|              the size of the field across the bore
The half turn reverses both components across the bore whatever way the sensor sits. A half turn that misses by
delta scales m by cos(delta / 2). The half turn does NOT remove the room's own field or a position error.
Per holder, the mean of m over its repeats; the ratio is S4 over S3. The +-266 mT range is the stated result; the
+-133 mT range is a check on it (its linear range is a typical figure with no stated minimum).
The three sensor axes have their own gains (5 % typical), so the seat must hold one sensor axis along the field:
"share" below is the part of m carried by the largest axis, and a warning is printed under 0.996 (5 degrees).
"""
import argparse, csv, math, os, sys, time

FIELDS = ["holder", "rep", "seating", "shim_mm", "host_time", "ms", "range_mT", "x_uT", "y_uT", "z_uT", "temp_C",
          "conv_status", "dev_status"]
RANGES = (266, 133)
PRIMARY = 266


def parse_line(line):
    """A data line 'D,ms,range,x,y,z,temp,conv,dev' -> (ms, range, x, y, z, temp, conv, dev), anything else -> None."""
    p = line.strip().split(",")
    if len(p) != 9 or p[0] != "D":
        return None
    try:
        d = (int(p[1]), int(p[2]), float(p[3]), float(p[4]), float(p[5]), float(p[6]), int(p[7]), int(p[8]))
    except ValueError:
        return None
    return d if d[1] in RANGES else None


def take(lines, n, skip=2):
    """n data lines of EACH range from an iterable of text lines. The first `skip` data lines after the start and
    after every change of range are dropped (a stale buffer; a conversion begun on the old range)."""
    out = {r: [] for r in RANGES}
    last = None; drop = skip
    for line in lines:
        d = parse_line(line)
        if d is None:
            continue
        if d[1] != last:
            if last is not None:
                drop = skip
            last = d[1]
        if drop > 0:
            drop -= 1
            continue
        if len(out[d[1]]) < n:
            out[d[1]].append(d)
        if all(len(v) == n for v in out.values()):
            return out[RANGES[0]] + out[RANGES[1]]
    raise RuntimeError("the stream ended with %s of %d data lines per range" % ({r: len(v) for r, v in out.items()}, n))


def append_rows(path, holder, rep, seating, shim, samples):
    new = not os.path.exists(path)
    with open(path, "a", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        if new:
            w.writerow(FIELDS)
        for d in samples:
            w.writerow([holder, rep, seating, shim, "%.3f" % time.time()] + list(d))


def load(path):
    """{(holder, rep, seating): [(range, x, y, z, temp, conv, dev, shim), ...]}"""
    rows = {}
    with open(path, newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            key = (r["holder"], int(r["rep"]), r["seating"])
            rows.setdefault(key, []).append((int(r["range_mT"]), float(r["x_uT"]), float(r["y_uT"]), float(r["z_uT"]),
                                             float(r["temp_C"]), int(r["conv_status"]), int(r["dev_status"]),
                                             float(r["shim_mm"])))
    return rows


def mean_sd(vals):
    n = len(vals)
    mu = sum(vals) / n
    sd = math.sqrt(sum((v - mu) ** 2 for v in vals) / (n - 1)) if n > 1 else float("nan")
    return mu, sd


def seating_stats(samples, rng):
    """samples: (range, x, y, z, temp, conv, dev, ...) -> statistics of those on range `rng`, or None if there are none."""
    sel = [s for s in samples if s[0] == rng]
    if not sel:
        return None
    cols = list(zip(*sel))
    return {"n": len(sel), "mean": [mean_sd(c)[0] for c in cols[1:4]], "sd": [mean_sd(c)[1] for c in cols[1:4]],
            "temp": mean_sd(cols[4])[0], "peak": max(abs(v) for c in cols[1:4] for v in c),
            "status": sorted(set(zip(cols[5], cols[6])))}


def repeat_result(a, b):
    """a, b: seating_stats of the A and B seatings of one repeat."""
    h = [(x - y) / 2 for x, y in zip(a["mean"], b["mean"])]
    s = [(x + y) / 2 for x, y in zip(a["mean"], b["mean"])]
    m = math.sqrt(sum(v * v for v in h))
    var_h = [(sa ** 2 / a["n"] + sb ** 2 / b["n"]) / 4 for sa, sb in zip(a["sd"], b["sd"])]
    se_m = math.sqrt(sum((hi / m) ** 2 * vi for hi, vi in zip(h, var_h))) if m > 0 and a["n"] > 1 and b["n"] > 1 else float("nan")
    return {"h": h, "s": s, "m": m, "se_m_samples": se_m, "temp": (a["temp"] + b["temp"]) / 2,
            "share": max(abs(v) for v in h) / m if m > 0 else float("nan")}


def analyse(rows, rng=PRIMARY, predict=None):
    """rows: {(holder, rep, seating): [samples]} -> per-holder results, the ratio and warnings, for one range."""
    res = {"range": rng, "holders": {}, "warnings": []}
    for holder in sorted({k[0] for k in rows}):
        reps = []
        for rep in sorted({k[1] for k in rows if k[0] == holder}):
            a = seating_stats(rows.get((holder, rep, "A"), []), rng); b = seating_stats(rows.get((holder, rep, "B"), []), rng)
            if a is None or b is None:
                res["warnings"].append("%s repeat %d has no A and B pair on this range: left out" % (holder, rep))
                continue
            for tag, st in (("A", a), ("B", b)):
                if st["peak"] > 0.9 * rng * 1000.0:
                    res["warnings"].append("%s repeat %d seating %s: a reading above 90%% of full scale" % (holder, rep, tag))
            r = repeat_result(a, b); r["rep"] = rep; r["a"] = a; r["b"] = b
            if r["share"] < 0.996:
                res["warnings"].append("%s repeat %d: the largest sensor axis carries %.4f of the field; the seat is "
                                       "more than 5 degrees off a sensor axis" % (holder, rep, r["share"]))
            reps.append(r)
        if not reps:
            continue
        mu, sd = mean_sd([r["m"] for r in reps])
        if len(reps) > 1:
            se = sd / math.sqrt(len(reps)); basis = "spread of the repeats"
        else:
            se = reps[0]["se_m_samples"]; basis = "sample scatter only (one repeat)"
        res["holders"][holder] = {"reps": reps, "m": mu, "sd_reps": sd, "se": se, "basis": basis,
                                  "temp": sum(r["temp"] for r in reps) / len(reps)}
    H = res["holders"]
    if "S3" in H and "S4" in H:
        ratio = H["S4"]["m"] / H["S3"]["m"]
        rel = math.sqrt((H["S3"]["se"] / H["S3"]["m"]) ** 2 + (H["S4"]["se"] / H["S4"]["m"]) ** 2)
        res["ratio"] = ratio; res["ratio_se"] = ratio * rel
        dT = H["S4"]["temp"] - H["S3"]["temp"]
        res["temp_diff"] = dT
        if abs(dT) > 1.0:
            res["warnings"].append("sensor temperature differs between the holders by %+.2f C" % dT)
        if predict:
            res["predict"] = predict; res["diff_pct"] = 100 * (ratio / predict - 1)
    return res


def zfit(rows, rng=PRIMARY):
    """S4, seating A, one shim height per repeat: parabola through (shim, reading along the mean direction)."""
    pts = []
    for (holder, rep, seating), samples in sorted(rows.items()):
        st = seating_stats(samples, rng)
        if holder == "S4" and seating == "A" and st is not None:
            pts.append((samples[0][7], st["mean"]))
    if len(pts) < 3:
        raise RuntimeError("a fit needs at least three shim heights; found %d" % len(pts))
    mean_dir = [sum(p[1][i] for p in pts) / len(pts) for i in range(3)]
    nrm = math.sqrt(sum(v * v for v in mean_dir)); u = [v / nrm for v in mean_dir]
    xs = [p[0] for p in pts]; ys = [sum(p[1][i] * u[i] for i in range(3)) for p in pts]
    # least squares for y = c0 + c1 x + c2 x^2 by the normal equations
    S = [sum(x ** k for x in xs) for k in range(5)]
    T = [sum(y * x ** k for x, y in zip(xs, ys)) for k in range(3)]
    A = [[S[0], S[1], S[2], T[0]], [S[1], S[2], S[3], T[1]], [S[2], S[3], S[4], T[2]]]
    for i in range(3):
        piv = max(range(i, 3), key=lambda r: abs(A[r][i])); A[i], A[piv] = A[piv], A[i]
        for r in range(3):
            if r != i:
                f = A[r][i] / A[i][i]
                A[r] = [a - f * b for a, b in zip(A[r], A[i])]
    c0, c1, c2 = (A[i][3] / A[i][i] for i in range(3))
    vertex = -c1 / (2 * c2)
    resid = max(abs(y - (c0 + c1 * x + c2 * x * x)) for x, y in zip(xs, ys))
    return {"points": list(zip(xs, ys)), "vertex_shim": vertex, "peak": c0 + c1 * vertex + c2 * vertex ** 2,
            "c2": c2, "is_maximum": c2 < 0, "inside": min(xs) < vertex < max(xs), "worst_residual": resid}


def report(res):
    print("---- range +-%d mT%s ----" % (res["range"], " (the stated result)" if res["range"] == PRIMARY else " (check)"))
    for holder, H in res["holders"].items():
        print("== %s ==" % holder)
        for r in H["reps"]:
            print("  repeat %d: field across the bore %.3f mT (from samples: +- %.4f); h = (%.1f, %.1f, %.1f) uT, share %.5f;"
                  " half-sum s = (%.1f, %.1f, %.1f) uT; temperature %.2f C; n = %d, %d; status values %s %s"
                  % (r["rep"], r["m"] / 1000, r["se_m_samples"] / 1000, r["h"][0], r["h"][1], r["h"][2], r["share"],
                     r["s"][0], r["s"][1], r["s"][2], r["temp"], r["a"]["n"], r["b"]["n"], r["a"]["status"], r["b"]["status"]))
        print("  mean %.3f mT, standard error %.4f mT (%s), %d repeat(s), temperature %.2f C"
              % (H["m"] / 1000, H["se"] / 1000, H["basis"], len(H["reps"]), H["temp"]))
    if "ratio" in res:
        print("ratio S4 / S3 = %.5f +- %.5f   (sensor temperature S4 minus S3: %+.2f C)" % (res["ratio"], res["ratio_se"], res["temp_diff"]))
        if "predict" in res:
            print("predicted %.4f; measured minus predicted %+.3f %%" % (res["predict"], res["diff_pct"]))
    for w in res["warnings"]:
        print("WARNING:", w)


def main(argv):
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("ports")
    c = sub.add_parser("collect")
    c.add_argument("--port", required=True); c.add_argument("--out", required=True)
    c.add_argument("--holder", required=True, choices=["S3", "S4", "EMPTY"]); c.add_argument("--rep", type=int, required=True)
    c.add_argument("--seating", required=True, choices=["A", "B"]); c.add_argument("--n", type=int, default=16)
    c.add_argument("--shim", type=float, default=3.0)
    a = sub.add_parser("analyse")
    a.add_argument("csv"); a.add_argument("--predict", type=float, default=None)
    z = sub.add_parser("zfit")
    z.add_argument("csv")
    args = ap.parse_args(argv)
    if args.cmd == "ports":
        from serial.tools import list_ports
        for p in list_ports.comports():
            print(p.device, "|", p.description, "|", p.hwid)
    elif args.cmd == "collect":
        import serial
        with serial.Serial(args.port, 115200, timeout=5) as ser:
            ser.reset_input_buffer()

            def lines():
                while True:
                    raw = ser.readline()
                    if not raw:
                        return
                    yield raw.decode("utf-8", "replace")
            samples = take(lines(), args.n)
        append_rows(args.out, args.holder, args.rep, args.seating, args.shim, samples)
        for rng in RANGES:
            st = seating_stats([s[1:] for s in samples], rng)
            print("%s repeat %d seating %s shim %.1f, +-%d mT: %d samples, mean (%.1f, %.1f, %.1f) uT, sd (%.1f, %.1f, %.1f), %.2f C"
                  % (args.holder, args.rep, args.seating, args.shim, rng, st["n"], *st["mean"], *st["sd"], st["temp"]))
    elif args.cmd == "analyse":
        rows = load(args.csv)
        out = {}
        for rng in RANGES:
            out[rng] = analyse(rows, rng, args.predict)
            report(out[rng])
        if all("ratio" in out[r] for r in RANGES):
            print("ratio on +-133 over ratio on +-266, minus 1: %+.4f %%" % (100 * (out[133]["ratio"] / out[266]["ratio"] - 1)))
    else:
        f = zfit(load(args.csv))
        for x, y in f["points"]:
            print("shim %.2f mm: %.1f uT" % (x, y))
        print("fitted extremum at shim %.2f mm (%s, %s the heights tried), value %.1f uT, worst residual %.1f uT"
              % (f["vertex_shim"], "a maximum" if f["is_maximum"] else "A MINIMUM: not the S4 peak",
                 "inside" if f["inside"] else "OUTSIDE", f["peak"], f["worst_residual"]))


if __name__ == "__main__":
    main(sys.argv[1:])
