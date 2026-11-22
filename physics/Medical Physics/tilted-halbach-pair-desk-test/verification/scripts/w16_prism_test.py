"""Tests of w16_prism.field, each with a mutant that must fail it.

T1  brute force: the block cut into cubic cells of 0.5 mm, each a point dipole, at near points (outside the block).
    (Non-cubic cells left a 5e-3 midpoint error on the bar and the slab; the reference was fixed, not the tolerance.)
T2  on-axis closed form for a block magnetised along its thickness D, pole face L x W, at distance z from the face:
    B = (Br/pi) [atan(L W/(2 z sqrt(4 z^2 + L^2 + W^2))) - atan(L W/(2 (z+D) sqrt(4 (z+D)^2 + L^2 + W^2)))].
T3  far field against the point dipole.
T4  a rotated block: field computed through Q equals the field of the unrotated block at rotated points, and the
    brute force in the laboratory frame.
T5  div B = 0 and curl B = 0 by central differences at near points.
Mutants: M1 flips the sign of the atan term; M2 swaps the roles of the two log terms; M3 drops the stable branch
sign (uses +ln for hx). Each must fail T1.
Usage: python -I -B w16_prism_test.py   (run from anywhere; imports w16_prism from this file's directory)
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import w16_prism as wp

rng = np.random.default_rng(7)
I3 = np.eye(3)


def brute(points, r0, dims, Q, Br, cell=0.0005):
    """Midpoint sum over cubic cells of edge `cell` (the block edges here are whole multiples of it)."""
    ns = [int(round(d / cell)) for d in dims]
    assert all(abs(n * cell - d) < 1e-12 for n, d in zip(ns, dims))
    g = [(np.arange(n) + 0.5) / n - 0.5 for n in ns]
    X, Y, Z = np.meshgrid(g[0] * dims[0], g[1] * dims[1], g[2] * dims[2], indexing="ij")
    loc = np.stack([X.ravel(), Y.ravel(), Z.ravel()], axis=1)
    cells = loc @ Q.T + np.asarray(r0)[None, :]
    m = Q[:, 2] * Br * cell ** 3
    out = np.zeros((len(points), 3))
    for c in cells:
        out += wp.dipole_field(points, c, m)
    return out


def near_points(dims, k, lo=0.75, hi=1.6):
    """Points outside the block, at 0.75 to 1.6 half-diagonals from the centre and at least 0.2 edge from the surface."""
    half = 0.5 * np.asarray(dims)
    pts = []
    while len(pts) < k:
        p = rng.normal(size=3)
        p *= rng.uniform(lo, hi) * np.linalg.norm(half) / np.linalg.norm(p)
        if np.any(np.abs(p) > half + 0.2 * min(dims)):
            pts.append(p)
    return np.array(pts)


def relerr(a, b):
    return np.abs(a - b).max() / np.abs(b).max()


results = []
Br = 1.45
for name, dims in (("cube 12", (0.012, 0.012, 0.012)), ("bar 12x12x50, M along 12", (0.012, 0.050, 0.012)),
                   ("slab 20x10x5, M along 5", (0.020, 0.010, 0.005))):
    P = near_points(dims, 40)
    e1 = relerr(wp.field(P, np.zeros(3), dims, I3, Br), brute(P, np.zeros(3), dims, I3, Br))
    results.append(("T1 brute force, %s" % name, e1, 2e-3))

for name, (L, W, D) in (("cube 12", (0.012, 0.012, 0.012)), ("slab 20x10x5", (0.020, 0.010, 0.005))):
    z = np.array([0.001, 0.003, 0.010, 0.030])
    P = np.stack([np.zeros_like(z), np.zeros_like(z), D / 2 + z], axis=1)
    f = lambda t: np.arctan(L * W / (2 * t * np.sqrt(4 * t * t + L * L + W * W)))
    want = (Br / np.pi) * (f(z) - f(z + D))
    got = wp.field(P, np.zeros(3), (L, W, D), I3, Br)
    results.append(("T2 on-axis formula, %s" % name, np.abs(got[:, 2] - want).max() / np.abs(want).max(), 1e-12))
    results.append(("T2 transverse components vanish on axis, %s" % name, np.abs(got[:, :2]).max() / np.abs(want).max(), 1e-12))

dims = (0.012, 0.012, 0.012)
for mult, tol in ((10, 2e-4), (40, 1e-6)):
    v = rng.normal(size=(30, 3)); P = mult * 0.012 * v / np.linalg.norm(v, axis=1)[:, None]
    m = np.array([0, 0, Br * 0.012 ** 3])
    results.append(("T3 far field at %d edges vs dipole" % mult, relerr(wp.field(P, np.zeros(3), dims, I3, Br), wp.dipole_field(P, np.zeros(3), m)), tol))

ang = 0.7
Q = wp.frame_inplane(ang)
r0 = np.array([0.03, -0.01, 0.004])
dimsb = (0.012, 0.050, 0.012)
Ploc = near_points(dimsb, 30)
P = Ploc @ Q.T + r0
results.append(("T4 rotated block vs unrotated at rotated points", relerr(wp.field(P, r0, dimsb, Q, Br), wp.field(Ploc, np.zeros(3), dimsb, I3, Br) @ Q.T), 1e-12))
results.append(("T4 rotated block vs brute force in the lab frame", relerr(wp.field(P, r0, dimsb, Q, Br), brute(P, r0, dimsb, Q, Br)), 2e-3))
results.append(("T4 magnetisation axis is (cos, sin, 0)", np.abs(Q[:, 2] - np.array([np.cos(ang), np.sin(ang), 0])).max(), 1e-15))
results.append(("T4 Q is a proper rotation", abs(np.linalg.det(Q) - 1) + np.abs(Q.T @ Q - I3).max(), 1e-14))

h = 1e-5
Pn = near_points(dims, 20, lo=1.0, hi=1.5)
J = np.zeros((len(Pn), 3, 3))
for k in range(3):
    d = np.zeros(3); d[k] = h
    J[:, :, k] = (wp.field(Pn + d, np.zeros(3), dims, I3, Br) - wp.field(Pn - d, np.zeros(3), dims, I3, Br)) / (2 * h)
scale = np.abs(J).max()
div = np.abs(J[:, 0, 0] + J[:, 1, 1] + J[:, 2, 2]).max() / scale
curl = max(np.abs(J[:, 2, 1] - J[:, 1, 2]).max(), np.abs(J[:, 0, 2] - J[:, 2, 0]).max(), np.abs(J[:, 1, 0] - J[:, 0, 1]).max()) / scale
results.append(("T5 div B = 0 (central differences)", div, 1e-5))
results.append(("T5 curl B = 0 (central differences)", curl, 1e-5))

ok = True
for name, err, tol in results:
    flag = "PASS" if err <= tol else "FAIL"
    ok &= err <= tol
    print("%-58s %.2e  (tol %.0e)  %s" % (name, err, tol, flag))

# mutants: each must FAIL T1 on the cube
P = near_points(dims, 40)
ref = brute(P, np.zeros(3), dims, I3, Br)
orig_sheet = wp._sheet


def mutant(kind):
    def sheet(x, y, w, a, b):
        hx, hy, hz = orig_sheet(x, y, w, a, b)
        if kind == "M1 atan sign flipped":
            return hx, hy, -hz
        if kind == "M2 log terms swapped":
            return hy, hx, hz
        return -hx, hy, hz
    return sheet


print("mutants (each must fail T1 at tolerance 2e-3):")
for kind in ("M1 atan sign flipped", "M2 log terms swapped", "M3 hx sign flipped"):
    wp._sheet = mutant(kind)
    e = relerr(wp.field(P, np.zeros(3), dims, I3, Br), ref)
    caught = e > 2e-3
    ok &= caught
    print("  %-30s T1 error %.2e  %s" % (kind, e, "CAUGHT" if caught else "NOT CAUGHT"))
wp._sheet = orig_sheet
print("ALL PASS" if ok else "SOMETHING FAILED")
