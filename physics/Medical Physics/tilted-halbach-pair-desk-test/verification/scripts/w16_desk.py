"""Desk-scale envelope: plain Halbach rings of catalogue cube magnets, exact block kernel.

Geometry. N cubes of edge s, centres on a circle of radius R in the plane z = z0, at azimuth theta_i = 2 pi i/N.
Cube i is magnetised at in-plane angle 2 theta_i (dipolar Halbach, field along +x at the centre) and is turned
with its magnetisation (a square pocket at angle 2 theta_i). Holder rules used to fix R:
  clear bore:   R - s/sqrt(2) - WALL >= bore/2          (a cube turned 45 degrees reaches s/sqrt(2) inward)
  neighbours:   2 R sin(pi/N) >= sqrt(2) s + RIB        (circumscribed circles of adjacent pockets, plus a rib)
R is the smallest radius meeting both, rounded up to 0.1 mm.
Two-ring stack: the same ring at z = +h and z = -h, h >= s/2 + FLOOR/2, h chosen to flatten Bx along z.

Output per design: R, B0 from blocks and from point dipoles, and the largest |Bx - B0|/B0 over a probe travel
of +-TRAVEL along x, along y and along z (through the centre), in per cent. Br = 1.45 T (N52 nominal); the field
is linear in Br, so N42 at 1.30 T is these numbers times 1.30/1.45.
Usage: python -I -B w16_desk.py [out.json]
"""
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import w16_prism as wp

WALL, RIB, FLOOR = 1.2e-3, 1.0e-3, 1.0e-3
BR = 1.45
TRAVEL = 5e-3
CUBES = (5e-3, 6.35e-3, 10e-3, 12.7e-3)
BORES = (20e-3, 30e-3, 40e-3)
NS = (8, 12, 16)


def radius(s, n, bore):
    r = max(bore / 2 + WALL + s / np.sqrt(2), (np.sqrt(2) * s + RIB) / (2 * np.sin(np.pi / n)))
    return np.ceil(r * 1e4 - 1e-9) / 1e4


def ring(s, n, R, z0=0.0, angle_of=lambda th: 2 * th):
    out = []
    for i in range(n):
        th = 2 * np.pi * i / n
        out.append((np.array([R * np.cos(th), R * np.sin(th), z0]), (s, s, s), wp.frame_inplane(angle_of(th)), BR))
    return out


def dipole_b0(blocks):
    b = np.zeros((1, 3))
    for r0, dims, Q, Br in blocks:
        b += wp.dipole_field(np.zeros((1, 3)), r0, Br * dims[0] * dims[1] * dims[2] * Q[:, 2])
    return b[0]


def lines(blocks, travel=TRAVEL, k=21):
    t = np.linspace(-travel, travel, k)
    res = {}
    b0 = wp.field_many(np.zeros((1, 3)), blocks)[0]
    for name, axis in (("x", 0), ("y", 1), ("z", 2)):
        P = np.zeros((k, 3)); P[:, axis] = t
        B = wp.field_many(P, blocks)
        res[name] = float(np.abs(B[:, 0] - b0[0]).max() / abs(b0[0]) * 100)
        res[name + "_perp"] = float(np.abs(B[:, 1:]).max() / abs(b0[0]) * 100)
    return b0, res


def best_pair(s, n, R):
    hmin = s / 2 + FLOOR / 2
    best = None
    for h in np.arange(hmin, hmin + 0.6 * R + 1e-9, 0.1e-3):
        blocks = ring(s, n, R, h) + ring(s, n, R, -h)
        t = np.linspace(-TRAVEL, TRAVEL, 11)
        P = np.zeros((11, 3)); P[:, 2] = t
        bx = wp.field_many(P, blocks)[:, 0]
        dev = np.abs(bx - bx[5]).max() / abs(bx[5])
        if best is None or dev < best[0]:
            best = (dev, h)
    return best[1]


rows = []
print("Br = %.2f T; wall %.1f mm, rib %.1f mm, floor %.1f mm; travel +-%.0f mm" % (BR, WALL * 1e3, RIB * 1e3, FLOOR * 1e3, TRAVEL * 1e3))
print("%5s %3s %5s | %6s %7s %7s %6s | %6s %6s %6s | pair: %5s %7s %6s %6s %6s"
      % ("cube", "N", "bore", "R mm", "B0 mT", "dip mT", "dip %", "dx %", "dy %", "dz %", "h mm", "B0 mT", "dx %", "dy %", "dz %"))
for s in CUBES:
    for bore in BORES:
        for n in NS:
            R = radius(s, n, bore)
            blocks = ring(s, n, R)
            b0, res = lines(blocks)
            bd = dipole_b0(blocks)
            h = best_pair(s, n, R)
            pair = ring(s, n, R, h) + ring(s, n, R, -h)
            b0p, resp = lines(pair)
            row = dict(cube_mm=s * 1e3, n=n, bore_mm=bore * 1e3, R_mm=R * 1e3, clear_bore_mm=2 * (R - s / np.sqrt(2) - WALL) * 1e3,
                       B0_mT=b0[0] * 1e3, B0_dipole_mT=bd[0] * 1e3, single=res, pair_h_mm=h * 1e3, pair_B0_mT=b0p[0] * 1e3, pair=resp,
                       b0_offaxis_mT=[b0[1] * 1e3, b0[2] * 1e3])
            rows.append(row)
            print("%5.2f %3d %5.0f | %6.1f %7.2f %7.2f %+6.1f | %6.2f %6.2f %6.2f | pair: %5.1f %7.2f %6.2f %6.2f %6.2f"
                  % (s * 1e3, n, bore * 1e3, R * 1e3, b0[0] * 1e3, bd[0] * 1e3, (bd[0] / b0[0] - 1) * 100,
                     res["x"], res["y"], res["z"], h * 1e3, b0p[0] * 1e3, resp["x"], resp["y"], resp["z"]))

if len(sys.argv) > 1:
    json.dump(rows, open(sys.argv[1], "w", encoding="utf-8"), indent=1)
    print("wrote", sys.argv[1])
