"""Cold audit of Wren's Prop 5 at U = 1, with MY operator: D_eff = M2 = I - V^dag D_PV^-1 U (her a5-scaled D_PV is
a5 times mine), against the limit form 1/2 (I + tau X/w), tau = tanh(S w/2), X = b + i sin p0 gamma0 + i sin k1 gamma1,
w = sqrt(b^2 + sin^2 p0 + sin^2 k1), b = 2 - m - cos k1 - cos p0. Basis-free comparisons: det, and sorted singular
values. Claim: error O(a5) at fixed S, uniformly on [-pi, pi]^2. Also: Sylvester (det M2 = det D_DW / det D_PV) and
min over the torus of w (claim >= min(m, 2 - m))."""
import sys
import numpy as np

sys.path.insert(0, ".")
from check_letter75_h4 import blocks  # noqa: E402
from edge_modes_dw import g0, g1  # noqa: E402


def M2(L5, k1, m, a5, z):
    B0, Bp, Bm = blocks(L5, k1, m, a5, 1.0)
    D = B0 + z * Bp + Bm / z
    n = 2 * L5
    U = np.zeros((n, 2), complex); V = np.zeros((n, 2), complex)
    U[2 * (L5 - 1) + 1, 0] = 1; U[0, 1] = 1          # |L5-1> e-, |0> e+   (spin 0 = e+, 1 = e-)
    V[1, 0] = 1; V[2 * (L5 - 1), 1] = 1              # |0> e-, |L5-1> e+
    M = np.eye(2) - V.conj().T @ np.linalg.solve(D, U) / a5
    Bd = blocks(L5, k1, m, a5, 0.0)
    Dd = Bd[0] + z * Bd[1] + Bd[2] / z
    syl = np.linalg.slogdet(Dd)
    syp = np.linalg.slogdet(D)
    ratio = syl[0] / syp[0] * np.exp(syl[1] - syp[1])
    return M, abs(np.linalg.det(M) - ratio) / abs(ratio)


def limit_form(p0, k1, m, S):
    b = 2 - m - np.cos(k1) - np.cos(p0)
    w = np.sqrt(b * b + np.sin(p0) ** 2 + np.sin(k1) ** 2)
    X = b * np.eye(2) + 1j * np.sin(p0) * g0 + 1j * np.sin(k1) * g1
    tau = np.tanh(S * w / 2)
    return 0.5 * (np.eye(2) + tau * X / w), w


def main():
    grid = np.linspace(-np.pi, np.pi, 9)
    wmin = {m: min(limit_form(p, k, m, 1.0)[1] for p in np.linspace(-np.pi, np.pi, 721)
                   for k in np.linspace(-np.pi, np.pi, 181)) for m in (0.3, 0.7)}
    print("min over torus of w: m=0.3 %.6f, m=0.7 %.6f (claim >= min(m, 2 - m))" % (wmin[0.3], wmin[0.7]))
    for m in (0.3, 0.7):
        for S in (4.0, 8.0):
            line = []
            for a5 in (0.2, 0.1, 0.05, 0.025):
                L5 = int(round(S / a5))
                ed, es, sy = 0.0, 0.0, 0.0
                for p0 in grid:
                    for k1 in grid:
                        M, s = M2(L5, k1, m, a5, np.exp(1j * p0))
                        F, w = limit_form(p0, k1, m, S)
                        ed = max(ed, abs(np.linalg.det(M) - np.linalg.det(F)))
                        es = max(es, np.abs(np.linalg.svd(M, compute_uv=False)
                                            - np.linalg.svd(F, compute_uv=False)).max())
                        sy = max(sy, s)
                line.append("a5=%.3f (L5=%d): |det diff| %.2e, sv diff %.2e, Sylvester %.0e" % (a5, L5, ed, es, sy))
            print("m=%.1f S=%.0f: %s" % (m, S, "; ".join(line)))
            sys.stdout.flush()
    print("DONE")


if __name__ == "__main__":
    main()
