"""Section through the W16 tilted holder (S4) with the sensor seat in place, cut by the plane y = 0.

Every outline is a slice of the solids the print files are made from; nothing is drawn by hand except the labels.
Usage: python -I -B -X utf8 render_seat.py <out.png> [shim_mm]
"""
import os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.path import Path
from matplotlib.patches import PathPatch
from manifold3d import Manifold

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import holder as H
import sensor_seat as S


def section(m):
    """Polygons of m cut by y = 0, as (x, z)."""
    return [np.asarray(p) for p in m.rotate([-90, 0, 0]).slice(0.0).to_polygons()]     # (x, y, z) -> (x, z, -y)


def draw(ax, polys, face, edge="k", lw=0.6, z=1):
    if not polys:
        return
    vs, cs = [], []
    for p in polys:
        vs += list(p) + [p[0]]
        cs += [Path.MOVETO] + [Path.LINETO] * (len(p) - 1) + [Path.CLOSEPOLY]
    ax.add_patch(PathPatch(Path(vs, cs), facecolor=face, edgecolor=edge, lw=lw, zorder=z))


def main(out, t):
    outer, inner, P = H.build("max", 9.60)[:3]
    ring = outer + inner
    hold = ring + ring.rotate([180, 0, 0])
    cubes = Manifold()
    for r0, Q in P:
        cubes += H.cube_at(r0, Q, H.NOM)
    cubes += cubes.rotate([180, 0, 0])
    seat = S.in_holder(S.build(), t)
    plate, comp = S.board()
    plate = S.in_holder(plate, t); comp = S.in_holder(comp, t)
    shims = Manifold(); z = S.Z_TOP; used = []
    left = t
    for s in sorted(S.SHIMS, reverse=True):
        if left >= s - 1e-9:
            shims += S.shim(s).translate([0, 0, z]); z += s; left -= s; used.append(s)

    fig, ax = plt.subplots(figsize=(11.5, 5.2), dpi=170)
    draw(ax, section(hold), "#c9ced6")
    draw(ax, section(cubes), "#8a8f98")
    draw(ax, section(shims), "#f2c14e")
    draw(ax, section(seat), "#7fb7be")
    draw(ax, section(comp), "#e8e8e8", lw=0.4, z=3)
    draw(ax, section(plate), "#2e7d32", z=4)
    ez = S.Z_TOP + t - (S.Z_TOP + S.DEEP + S.BL / 2) + S.BL / 2
    ax.plot([0], [ez], marker="o", ms=5, color="#c62828", zorder=6)
    ax.axhline(0, color="k", lw=0.5, ls=(0, (6, 4)), zorder=5)
    ax.axvline(0, color="k", lw=0.5, ls=(0, (1, 3)), zorder=5)
    ax.plot([-70, 70], [-S.Z_TOP, -S.Z_TOP], color="k", lw=1.2)

    def note(text, xy, xt, ha="left"):
        ax.annotate(text, xy=xy, xytext=xt, fontsize=8.5, ha=ha, va="center",
                    arrowprops=dict(arrowstyle="-", lw=0.6, color="k"), zorder=7)

    note("sensor element\n(on the bore axis)", (0, ez), (9, -27), "left")
    note("board on edge,\nface toward the field (x)", (S.X_MID, ez + 9), (-97, 18), "left")
    note("shims: %s mm" % " + ".join("%.1f" % u for u in used) if used else "no shim", (26, S.Z_TOP + t / 2), (52, 30), "left")
    note("plug: flange rests on the shims", (-24, S.Z_TOP + t + 1.5), (-97, 30), "left")
    note("neck and waist: clear of\nthe bore's tight first layers", (-18.9, 0.0), (-97, 6), "left")
    note("magnet cube", (27.5, 12.1), (56, 14), "left")
    note("mid-plane", (40, 0), (56, 3), "left")
    note("3 mm push-out hole", (S.X_MID, -S.Z_TOP + 3.5 + t), (-9, -27), "right")
    ax.text(68, -S.Z_TOP - 1.6, "desk", fontsize=8.5, ha="right", va="top")
    ax.set_aspect("equal"); ax.set_xlim(-100, 85); ax.set_ylim(-30, 36)
    ax.set_xlabel("x (mm): the field direction"); ax.set_ylabel("z (mm): the bore axis")
    ax.set_title("W16 tilted holder with the sensor seat, cut through the middle (element at %+.1f mm, shim stack %.1f mm)"
                 % (ez, t), fontsize=10)
    fig.tight_layout(); fig.savefig(out)
    print("wrote %s; element at z = %+.3f with a %.1f mm stack (%s)" % (out, ez, t, used))


if __name__ == "__main__":
    main(sys.argv[1], float(sys.argv[2]) if len(sys.argv) > 2 else 3.0)
