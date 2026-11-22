"""The tilted rule of this build against the printed formulas of Rehberg and Blumler (arXiv:2502.18262).

Their "focused" ring of point dipoles: a dipole at angle phi on a ring of radius R is turned to the in-plane angle
alpha_f (their Eq. for alpha^f(f)) and tilted out of the plane by theta_f (their Eq. for theta^f(f)), so as to
maximise B_x on the axis at the height f R from the ring's plane. Our rule: m parallel to 3 (n.x) n - x, with n the
unit vector between the target point and the magnet. Each ring of our build is z0/R = 12.116/28.861 from the centre.

Printed here:
  1. alpha and |theta| by both routes at our twelve angles, and the largest difference;
  2. the centre field of a two-ring stack of POINT DIPOLES, tilted rule over plain Halbach rule (alpha = 2 phi), at
     our spacing and at their f = 0.576, beside the same ratio for our cubes (1.156179, from w16_seat_envelope.py);
  3. a control: a wrong formula (the sign of f^2 in alpha_f flipped) must disagree.

Usage: python -I -B -X utf8 w16_rehberg_cmp.py
"""
import numpy as np

R, Z0 = 28.861, 12.116
N = 12


def ours(phi, f):
    """Unit magnetisation of a dipole at (cos phi, sin phi, -f), target at the origin, wanted field along x."""
    n = np.array([np.cos(phi), np.sin(phi), -f]) / np.sqrt(1 + f * f)
    u = 3 * n[0] * n - np.array([1.0, 0.0, 0.0])
    return u / np.linalg.norm(u)


def theirs(phi, f, flip=1.0):
    """(alpha, theta) from the two printed formulas; alpha by atan2 as their footnote prescribes."""
    alpha = np.arctan2(3 * np.sin(2 * phi), 6 * np.cos(phi) ** 2 - 2 - flip * 2 * f * f)
    dalpha = alpha - 2 * phi
    theta = np.arctan(-6 * f * np.cos(phi) / (3 * np.cos(dalpha) + (1 - 2 * f * f) * np.cos(alpha)))
    return alpha, theta


def bx_centre(f, rule):
    """B_x at the origin from two rings of N unit point dipoles at z = -f and +f (ring radius 1), arbitrary units."""
    total = 0.0
    for k in range(N):
        phi = 2 * np.pi * k / N
        for s in (-1.0, 1.0):
            r = np.array([np.cos(phi), np.sin(phi), s * f])
            if rule == "tilted":
                n = r / np.linalg.norm(r)
                m = 3 * n[0] * n - np.array([1.0, 0.0, 0.0]); m /= np.linalg.norm(m)
            else:
                m = np.array([np.cos(2 * phi), np.sin(2 * phi), 0.0])
            d = -r; dn = np.linalg.norm(d)
            total += (3 * np.dot(m, d) * d[0] / dn ** 5 - m[0] / dn ** 3)
    return total


def wrap(a):
    return (a + np.pi) % (2 * np.pi) - np.pi


def main():
    f = Z0 / R
    print("f = z0/R = %.5f" % f)
    print("1. our rule against the printed formulas, twelve angles")
    print("   phi    alpha ours   alpha theirs   |tilt| ours   |tilt| theirs")
    worst = 0.0; worst_bad = 0.0
    for k in range(N):
        phi = 2 * np.pi * k / N
        u = ours(phi, f)
        a_o, t_o = np.arctan2(u[1], u[0]), np.arcsin(u[2])
        a_t, t_t = theirs(phi, f)
        a_b, t_b = theirs(phi, f, flip=-1.0)
        worst = max(worst, abs(wrap(a_o - a_t)), abs(abs(t_o) - abs(t_t)))
        worst_bad = max(worst_bad, abs(wrap(a_o - a_b)), abs(abs(t_o) - abs(t_b)))
        print("   %5.1f  %10.4f   %10.4f     %8.4f      %8.4f" % (np.degrees(phi), np.degrees(a_o), np.degrees(a_t),
                                                                 abs(np.degrees(t_o)), abs(np.degrees(t_t))))
    print("   largest difference, either angle: %.2e deg" % np.degrees(worst))
    print("2. two rings of point dipoles, centre field, tilted rule over plain Halbach rule")
    for label, ff in (("our spacing, f = %.4f" % f, f), ("their f = 0.576", 0.576), ("f = 0 (one plane)", 0.0)):
        print("   %-26s ratio %.6f" % (label, bx_centre(ff, "tilted") / bx_centre(ff, "halbach")))
    print("   our cubes at our spacing (w16_seat_envelope.py): 1.156179")
    print("3. control: alpha_f with the sign of f^2 flipped differs from our rule by up to %.3f deg" % np.degrees(worst_bad))
    ok = np.degrees(worst) < 1e-9 and np.degrees(worst_bad) > 1.0
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
