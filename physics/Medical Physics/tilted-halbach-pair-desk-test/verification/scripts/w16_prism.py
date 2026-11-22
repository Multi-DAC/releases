"""Exact field of uniformly magnetised rectangular blocks (surface-charge closed form), outside the blocks.

A block is given by its centre r0 [m], its full edge lengths (L1, L2, L3) [m] along its own axes e1, e2, e3, and
its remanence Br [T]; the magnetisation is along e3. Q is the 3x3 matrix whose columns are e1, e2, e3 in the
laboratory frame. The relative permeability is taken as 1 (rigid magnetisation), as in the dipole scripts.

field(points, r0, dims, Q, Br) -> B [T] at each point, shape (n, 3).

Derivation: charge sheets +M and -M on the faces local z = +c and -c. For a sheet at local z0 with u = x - x',
v = y - y', w = z - z0, R = sqrt(u^2 + v^2 + w^2):
  Hx = (sigma/4 pi) [[ -ln(v + R) ]],  Hy = (sigma/4 pi) [[ -ln(u + R) ]],  Hz = (sigma/4 pi) [[ atan(u v/(w R)) ]]
with [[ ]] the double difference over u in [x - a, x + a], v in [y - b, y + b]. Outside the block B = mu0 H.
"""
import numpy as np

TINY = 1e-14


def _log_plus(v, R, rho2):
    """ln(v + R), stable for v < 0, where rho2 = R^2 - v^2 >= 0."""
    rho2 = np.maximum(rho2, TINY ** 2)
    return np.where(v >= 0, np.log(np.maximum(v + R, TINY)), np.log(rho2 / np.maximum(R - v, TINY)))


def _sheet(x, y, w, a, b):
    """Field per unit (sigma/4 pi) of a rectangular sheet |x'| <= a, |y'| <= b at height w below the field point."""
    w = np.where(np.abs(w) < TINY, TINY, w)
    hx = np.zeros_like(x); hy = np.zeros_like(x); hz = np.zeros_like(x)
    for su in (1, -1):
        for sv in (1, -1):
            u = x + su * a
            v = y + sv * b
            R = np.sqrt(u * u + v * v + w * w)
            s = su * sv
            hx += s * -_log_plus(v, R, u * u + w * w)
            hy += s * -_log_plus(u, R, v * v + w * w)
            hz += s * np.arctan(u * v / (w * R))
    return hx, hy, hz


def field_local(p, dims, Br):
    """B [T] at local points p (n, 3) of a block centred at the origin, magnetised along local z."""
    a, b, c = (0.5 * d for d in dims)
    x, y, z = p[:, 0], p[:, 1], p[:, 2]
    tx, ty, tz = _sheet(x, y, z - c, a, b)
    bx, by, bz = _sheet(x, y, z + c, a, b)
    return (Br / (4 * np.pi)) * np.stack([tx - bx, ty - by, tz - bz], axis=1)


def field(points, r0, dims, Q, Br):
    p = (np.atleast_2d(points) - np.asarray(r0)[None, :]) @ Q
    return field_local(p, dims, Br) @ Q.T


def field_many(points, blocks):
    """Sum over blocks; blocks is an iterable of (r0, dims, Q, Br)."""
    out = np.zeros((len(np.atleast_2d(points)), 3))
    for r0, dims, Q, Br in blocks:
        out += field(points, r0, dims, Q, Br)
    return out


def frame_inplane(angle):
    """Block whose magnetisation points at `angle` in the x-y plane; one edge stays along the laboratory z axis."""
    c, s = np.cos(angle), np.sin(angle)
    e3 = np.array([c, s, 0.0])
    e1 = np.array([-s, c, 0.0])
    e2 = np.cross(e3, e1)
    return np.stack([e1, e2, e3], axis=1)


def dipole_field(points, r0, moment_Tm3):
    """Point dipole with 'moment' Br*V [T m^3] (vector): B = (3 (m.n) n - m)/(4 pi r^3)."""
    d = np.atleast_2d(points) - np.asarray(r0)[None, :]
    r = np.linalg.norm(d, axis=1)
    n = d / r[:, None]
    m = np.asarray(moment_Tm3)[None, :]
    return (3 * (m * n).sum(1)[:, None] * n - m) / (4 * np.pi * r[:, None] ** 3)
