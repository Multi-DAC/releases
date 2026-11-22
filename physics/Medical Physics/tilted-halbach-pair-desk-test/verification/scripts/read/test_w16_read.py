"""Self-test of w16_read.py (v2) on made-up readings with a known field, offset and ratio. No hardware.

Usage: python -I -B -X utf8 test_w16_read.py
"""
import math, os, random, sys, tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import w16_read as W

B3, B4 = 89270.0, 103212.0          # uT, the centre fields at 3/8 inch and 1.30 T (out-seat-envelope.txt)
TRUE = B4 / B3
EYE = [[1.0, 0, 0], [0, 1.0, 0], [0, 0, 1.0]]


def rand_rot(rng):
    """A random rotation matrix (rows = sensor axes in the holder's frame)."""
    while True:
        a = [rng.gauss(0, 1) for _ in range(3)]; b = [rng.gauss(0, 1) for _ in range(3)]
        na = math.sqrt(sum(v * v for v in a)); a = [v / na for v in a]
        d = sum(x * y for x, y in zip(a, b)); b = [y - d * x for x, y in zip(a, b)]
        nb = math.sqrt(sum(v * v for v in b))
        if nb > 1e-3:
            b = [v / nb for v in b]
            c = [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
            return [a, b, c]


def edge_frame(phi_deg):
    """Board on edge, its normal (sensor z) at phi from the holder's x; sensor y along the bore."""
    p = math.radians(phi_deg)
    return [[-math.sin(p), math.cos(p), 0.0], [0.0, 0.0, 1.0], [math.cos(p), math.sin(p), 0.0]]


def reading(Q, field, offset, turn_deg, gain=(1.0, 1.0, 1.0)):
    """What the sensor shows: the seat turned by turn_deg about the bore (z) axis, then the sensor's own axes Q."""
    t = math.radians(turn_deg)
    f = [field[0] * math.cos(t) + field[1] * math.sin(t), -field[0] * math.sin(t) + field[1] * math.cos(t), field[2]]
    return [gain[i] * sum(Q[i][j] * f[j] for j in range(3)) + offset[i] for i in range(3)]


def make_rows(rng, Q, offset, noise=0.0, n=16, reps=3, delta=0.0, bz=(0.0, 0.0), reseat_deg=0.0, gain=(1.0, 1.0, 1.0),
              notch=(0.0, 0.0), gain133=None, shim=3.0):
    """Both ranges for every seating. gain133: per-axis gains on the +-133 range (default: the same as +-266)."""
    rows = {}
    for holder, B, z, nt in (("S3", B3, bz[0], notch[0]), ("S4", B4, bz[1], notch[1])):
        for rep in range(1, reps + 1):
            for seating, turn in (("A", nt), ("B", nt + 180.0 + delta)):
                turn += rng.gauss(0, reseat_deg) if reseat_deg else 0.0
                out = []
                for r, g in ((266, gain), (133, gain133 or gain)):
                    v = reading(Q, [B, 0.0, z], offset, turn, g)
                    out += [(r,) + tuple(c + (rng.gauss(0, noise) if noise else 0.0) for c in v) + (22.0, 1, 0, shim)
                            for _ in range(n)]
                rows[(holder, rep, seating)] = out
    return rows


def ratio_half_sum(rows):
    """MUTANT: the half-sum where the half-difference belongs."""
    m = {}
    for holder in ("S3", "S4"):
        a = W.seating_stats(rows[(holder, 1, "A")], 266); b = W.seating_stats(rows[(holder, 1, "B")], 266)
        m[holder] = math.sqrt(sum(((x + y) / 2) ** 2 for x, y in zip(a["mean"], b["mean"])))
    return m["S4"] / m["S3"]


def ratio_one_seating(rows):
    """MUTANT: the size of a single seating's reading, offset left in."""
    m = {h: math.sqrt(sum(v * v for v in W.seating_stats(rows[(h, 1, "A")], 266)["mean"])) for h in ("S3", "S4")}
    return m["S4"] / m["S3"]


def main():
    rng = random.Random(1121)
    fails = 0

    def check(name, ok, detail):
        nonlocal fails
        fails += 0 if ok else 1
        print("%s  %s: %s" % ("PASS" if ok else "FAIL", name, detail))

    # 1. exact recovery with no noise and equal gains: any sensor orientation, any offset, a field along the bore too
    worst = 0.0; worst_m = 0.0
    for _ in range(200):
        Q = rand_rot(rng); off = [rng.uniform(-1000, 1000) for _ in range(3)]
        rows = make_rows(rng, Q, off, reps=1, bz=(rng.uniform(-2000, 2000), rng.uniform(-2000, 2000)))
        for r in W.RANGES:
            res = W.analyse(rows, r)
            worst = max(worst, abs(res["ratio"] / TRUE - 1)); worst_m = max(worst_m, abs(res["holders"]["S3"]["m"] / B3 - 1))
    check("exact, 200 random orientations and offsets, both ranges", worst < 1e-12 and worst_m < 1e-12,
          "largest relative error of the ratio %.1e, of the S3 field %.1e" % (worst, worst_m))

    # 2. the two mutants must miss on the same data
    Q = rand_rot(rng); off = [300.0, -300.0, 300.0]
    rows = make_rows(rng, Q, off, reps=1)
    e_ok = abs(W.analyse(rows)["ratio"] / TRUE - 1)
    e_sum = abs(ratio_half_sum(rows) / TRUE - 1); e_one = abs(ratio_one_seating(rows) / TRUE - 1)
    check("mutant half-sum misses", e_sum > 0.01 and e_ok < 1e-12, "error %.3f against %.1e" % (e_sum, e_ok))
    check("mutant single seating misses", e_one > 1e-5 and e_ok < 1e-12,
          "error %.2e with a 300 uT offset on each axis" % e_one)

    # 3. a half turn that misses by delta scales the field by cos(delta/2) and leaves the ratio alone
    for delta in (2.0, 5.0):
        res = W.analyse(make_rows(rng, Q, off, reps=1, delta=delta))
        got = res["holders"]["S3"]["m"] / B3; want = math.cos(math.radians(delta) / 2)
        check("half turn off by %.0f deg" % delta, abs(got - want) < 1e-12 and abs(res["ratio"] / TRUE - 1) < 1e-12,
              "field scaled by %.6f (cos(delta/2) = %.6f), ratio error %.1e" % (got, want, abs(res["ratio"] / TRUE - 1)))

    # 4. noise 24 uT rms per sample, 16 samples, 5 repeats: null study
    errs = []; ses = []
    for _ in range(400):
        res = W.analyse(make_rows(rng, rand_rot(rng), [rng.uniform(-1000, 1000) for _ in range(3)], noise=24.0, reps=5))
        errs.append(res["ratio"] / TRUE - 1); ses.append(res["ratio_se"] / res["ratio"])
    mu, sd = W.mean_sd(errs)
    check("noise only, 400 trials", abs(mu) < 4 * sd / math.sqrt(len(errs)),
          "mean error %.2e, sd %.2e (relative); the script's own standard error averages %.2e" % (mu, sd, sum(ses) / len(ses)))

    # 5. re-seating scatter of 1 deg (sd) in the turn, no noise, 5 repeats
    errs = []
    for _ in range(400):
        errs.append(W.analyse(make_rows(rng, Q, off, reseat_deg=1.0, reps=5))["ratio"] / TRUE - 1)
    mu, sd = W.mean_sd(errs)
    check("re-seating 1 deg sd, 400 trials", abs(mu) < 4 * sd / math.sqrt(len(errs)) and sd < 1e-3,
          "mean error %.2e, sd %.2e (relative)" % (mu, sd))

    # 6. unequal axis gains (sensor z 5 % high), S4's seat 2 deg round from S3's: the error against the board's angle
    #    phi to the field. This is what the "share" warning is for; printed as a finding, with two bounds checked.
    found = {}
    for phi in (0.0, 5.0, 15.0, 45.0, 90.0):
        res = W.analyse(make_rows(rng, edge_frame(phi), off, reps=1, gain=(1.0, 1.0, 1.05), notch=(0.0, 2.0)))
        found[phi] = (100 * (res["ratio"] / TRUE - 1), min(r["share"] for H in res["holders"].values() for r in H["reps"]),
                      len([w for w in res["warnings"] if "largest sensor axis" in w]))
    print("      gain mismatch 5 %, notch 2 deg: " + "; ".join("phi %g: %+.4f %% (share %.4f, %d warning)" % ((p,) + found[p]) for p in found))
    check("gain mismatch is harmless on a sensor axis and not at 45 deg", abs(found[0.0][0]) < 0.01 and abs(found[45.0][0]) > 0.1,
          "%.4f %% at 0, %.4f %% at 45" % (found[0.0][0], found[45.0][0]))
    check("share warning fires at 15 and 45 deg, not at 0", found[0.0][2] == 0 and found[15.0][2] > 0 and found[45.0][2] > 0,
          "warnings %d / %d / %d" % (found[0.0][2], found[15.0][2], found[45.0][2]))

    # 7. the two ranges are analysed apart: a 0.2 % gain error on the narrow range shows in its ratio line only
    rows = make_rows(rng, edge_frame(0.0), off, reps=1, gain133=(1.0, 1.0, 1.002))
    rows[("S3", 1, "A")] = [s for s in rows[("S3", 1, "A")]]
    r266 = W.analyse(rows, 266)["holders"]["S3"]["m"]; r133 = W.analyse(rows, 133)["holders"]["S3"]["m"]
    check("ranges kept apart", abs(r266 / B3 - 1) < 1e-12 and abs(r133 / B3 - 1.002) < 1e-9, "S3 field x %.6f on 266, x %.6f on 133" % (r266 / B3, r133 / B3))

    # 8. the line parser and the collector
    good = "D,1234,266,-89270.5,12.0,300.1,22.31,1,0\r\n"
    check("parser takes a data line", W.parse_line(good) == (1234, 266, -89270.5, 12.0, 300.1, 22.31, 1, 0), repr(W.parse_line(good)))
    bad = ["# w16 code.py v2 id=0x02", "D,12,266,1.0,2.0", "D,x,266,1,2,3,4,5,6", "", "Traceback (most recent call last):",
           "E,1,266,2,3,4,5,6,7", "D,1,200,2,3,4,5,6,7", "D,1234,-89270.5,12.0,300.1,22.31,1,0"]
    check("parser refuses eight other lines", all(W.parse_line(b) is None for b in bad), "all None")
    stream = [bad[0]]
    k = 0
    for block in range(6):
        for _ in range(8):
            stream.append("D,%d,%d,%d.0,0.0,0.0,22.0,1,0\n" % (k, 266 if block % 2 == 0 else 133, k)); k += 1
    got = W.take(iter(stream), 9, skip=2)
    ids = {r: [g[0] for g in got if g[1] == r] for r in W.RANGES}
    check("collector drops 2 after each switch, takes 9 per range",
          ids[266] == [2, 3, 4, 5, 6, 7, 18, 19, 20] and ids[133] == [10, 11, 12, 13, 14, 15, 26, 27, 28], str(ids))
    try:
        W.take(iter(stream), 50); ended = False
    except RuntimeError:
        ended = True
    check("collector reports a stream that ends early", ended, "RuntimeError raised")

    # 9. through the file: write, load, analyse
    rows = make_rows(rng, edge_frame(1.0), off, noise=24.0, reps=5)
    d = tempfile.mkdtemp(prefix="w16read-")
    path = os.path.join(d, "run.csv")
    for (holder, rep, seating), samples in rows.items():
        W.append_rows(path, holder, rep, seating, 3.0, [(0,) + s[:7] for s in samples])
    r1 = W.analyse(rows, 266, 1.1562); r2 = W.analyse(W.load(path), 266, 1.1562)
    check("file round trip", abs(r1["ratio"] - r2["ratio"]) < 1e-9, "ratio %.6f from memory, %.6f from the file" % (r1["ratio"], r2["ratio"]))
    W.report(r2)

    # 10. warnings fire
    rows = make_rows(rng, EYE, off, reps=1)
    del rows[("S4", 1, "B")]
    w1 = W.analyse(rows)["warnings"]
    rows = make_rows(rng, EYE, [0, 0, 0], reps=1)
    rows[("S3", 1, "A")] = [(266, 250000.0, 0.0, 0.0, 22.0, 1, 0, 3.0)] * 5
    rows[("S4", 1, "A")] = [s[:4] + (24.0,) + s[5:] for s in rows[("S4", 1, "A")]]
    rows[("S4", 1, "B")] = [s[:4] + (24.0,) + s[5:] for s in rows[("S4", 1, "B")]]
    w2 = W.analyse(rows)["warnings"]
    check("warnings: missing pair, full scale, temperature", len(w1) == 1 and len(w2) == 2, "%s | %s" % (w1, w2))

    # 11. the height scan: S4's on-axis values from out-seat-envelope.txt (drop from the peak at |z| = 0, 1, 2, 3 mm),
    #     the true centre at shim 3.4 mm, seven shims 0 to 6 mm, seating A only, offset left in
    drop = {0: 0.0, 1: -107.0, 2: -430.0, 3: -978.0}

    def b_s4(z):         # even interpolation through the four values, uT
        a = abs(z); i = min(int(a), 2); f = a - i
        return B4 + drop[i] + f * (drop[i + 1] - drop[i])
    rows = {}
    for kk, shim in enumerate([0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0]):
        v = reading(edge_frame(0.0), [b_s4(shim - 3.4), 0.0, 0.0], off, 0.0)
        rows[("S4", kk + 1, "A")] = [(266,) + tuple(v) + (22.0, 1, 0, shim)] * 4
    f = W.zfit(rows)
    check("height scan finds the centre", abs(f["vertex_shim"] - 3.4) < 0.15 and f["is_maximum"] and f["inside"],
          "fitted shim %.3f mm for a true 3.400; maximum %s; inside %s" % (f["vertex_shim"], f["is_maximum"], f["inside"]))
    rows_min = {k: [(s[0], -s[1], -s[2], -s[3]) + s[4:] for s in v] for k, v in rows.items()}
    f2 = W.zfit({k: v for k, v in list(rows.items())[:3]})
    check("height scan says when the centre is outside the heights tried", not f2["inside"], "shims 0 to 2 only: vertex %.2f, inside %s" % (f2["vertex_shim"], f2["inside"]))

    print("failures: %d" % fails)
    return fails


if __name__ == "__main__":
    sys.exit(1 if main() else 0)
