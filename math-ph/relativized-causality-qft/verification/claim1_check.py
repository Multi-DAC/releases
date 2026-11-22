"""Blind check of Wren's Claim 1 (letter 42): a control of MY hand formula, in a two-mode toy.

Usage: python -I -B -X utf8 claim1_check.py

Hand result being checked. A_1 = 1/2 + eps*cos(lam*phi(g)), A_2 = 1 - A_1, relativised effects
E_i = sum_x p_x alpha_x(A_i), E_0 = (1-p) 1, Lueders channel Lam(X) = sum_i E_i^1/2 X E_i^1/2,
X = sum_z p_z cos(phi(h_z)), kick U = exp(i phi(f)) with sigma(f, h_z) = 0. Then
   signal := <U* Lam(X) U> - <Lam(X)>
           = -(eps^2 lam^4 / (2p)) <cos phi(h)> sum_z p_z M_z^2  + O(eps^2 lam^5) + O(eps^4),
   M_z = sum_x p_x sigma(h_z, g_x) sigma(f, g_x).
The same with the MIXTURE relativisation Om(X) = (1-p) X + sum_x p_x Lam_x(X), Lam_x the Lueders channel of
{alpha_x(A_i)}:   signal = -(eps^2 lam^4 / 2) <cos phi(h)> sum_x sum_z p_x p_z (sigma(h_z, g_x) sigma(f, g_x))^2.

Toy: mode 1 (Alice's side), mode 3 (Charlie's side), [X, P] = i, vacuum of both. Translations = a common
phase rotation by theta of both modes (what a time translation does to a massive field's large-scale modes).
   phi(g_theta) = X_1(theta) + X_3(theta),  phi(h_theta) = eta0 * P_3(theta),  phi(f) = tau * P_1.
   sigma(f, g_theta) = -tau cos(theta),  sigma(h_z, g_x) = -eta0 cos(theta_x - theta_z),  sigma(f, h) = 0.
Frame density: a finite list of (theta_j, p_j).

Controls:
   ZERO: single point theta = pi/2. Then sigma(f, g) = 0, U commutes with every Kraus operator: signal must be 0.
   ABSOLUTE: single point theta = 0, p = 1: the ordinary Sorkin signal, formula with M = eta0 * tau.
"""
import numpy as np

N = 36  # Fock truncation per mode


def ladder(n):
    a = np.zeros((n, n))
    for k in range(1, n):
        a[k - 1, k] = np.sqrt(k)
    return a


A = ladder(N)
X1 = (A + A.T) / np.sqrt(2)
P1 = (A - A.T) / (1j * np.sqrt(2))
I1 = np.eye(N)


def two(op1, op3):
    return np.kron(op1, op3)


def herm_fun(H, fun):
    w, v = np.linalg.eigh((H + H.conj().T) / 2)
    return (v * fun(w)) @ v.conj().T


def quad(theta):
    """X(theta), P(theta) on one mode."""
    return X1 * np.cos(theta) + P1 * np.sin(theta), -X1 * np.sin(theta) + P1 * np.cos(theta)


def run(thetas, ps, eps, lam, eta0, tau):
    p = float(sum(ps))
    dim = N * N
    Id = np.eye(dim)
    vac = np.zeros(dim, dtype=complex)
    vac[0] = 1.0
    U = herm_fun(two(tau * P1, I1), lambda w: np.exp(1j * w))
    kicked = U @ vac
    ks, Xs = [], []
    for th in thetas:
        Xq, Pq = quad(th)
        ks.append(herm_fun(two(Xq, I1) + two(I1, Xq), lambda w: np.cos(lam * w)))
        Xs.append(herm_fun(two(I1, eta0 * Pq), np.cos))
    Xrel = sum(pj * Xj for pj, Xj in zip(ps, Xs))
    # Lueders channel of the relativised effects
    K = sum(pj * kj for pj, kj in zip(ps, ks))
    E1 = (p / 2) * Id + eps * K
    E2 = (p / 2) * Id - eps * K
    lo = min(np.linalg.eigvalsh(E1)[0], np.linalg.eigvalsh(E2)[0])
    S1 = herm_fun(E1, lambda w: np.sqrt(np.clip(w, 0, None)))
    S2 = herm_fun(E2, lambda w: np.sqrt(np.clip(w, 0, None)))
    LamX = (1 - p) * Xrel + S1 @ Xrel @ S1 + S2 @ Xrel @ S2
    sig = (kicked.conj() @ LamX @ kicked - vac.conj() @ LamX @ vac).real
    # mixture of translated Lueders channels
    OmX = (1 - p) * Xrel
    for pj, kj in zip(ps, ks):
        a1 = herm_fun(0.5 * Id + eps * kj, np.sqrt)
        a2 = herm_fun(0.5 * Id - eps * kj, np.sqrt)
        OmX = OmX + pj * (a1 @ Xrel @ a1 + a2 @ Xrel @ a2)
    sig_mix = (kicked.conj() @ OmX @ kicked - vac.conj() @ OmX @ vac).real
    # hand formulas
    cosh_ = np.exp(-eta0 ** 2 / 4)
    M = [sum(px * (-eta0 * np.cos(tx - tz)) * (-tau * np.cos(tx)) for tx, px in zip(thetas, ps)) for tz in thetas]
    pred = -(eps ** 2 * lam ** 4 / (2 * p)) * cosh_ * sum(pz * Mz ** 2 for pz, Mz in zip(ps, M))
    pred_mix = -(eps ** 2 * lam ** 4 / 2) * cosh_ * sum(
        px * pz * (eta0 * np.cos(tx - tz) * tau * np.cos(tx)) ** 2
        for tx, px in zip(thetas, ps) for tz, pz in zip(thetas, ps))
    # the kick must not move X itself (sigma(f, h) = 0)
    leak = abs(kicked.conj() @ Xrel @ kicked - vac.conj() @ Xrel @ vac)
    return dict(p=p, min_eig_E=lo, sig=sig, pred=pred, sig_mix=sig_mix, pred_mix=pred_mix, leak=leak)


def show(tag, thetas, ps, eps, lam, eta0=1.0, tau=1.0):
    r = run(thetas, ps, eps, lam, eta0, tau)
    ratio = r["sig"] / r["pred"] if r["pred"] != 0 else float("nan")
    ratio_m = r["sig_mix"] / r["pred_mix"] if r["pred_mix"] != 0 else float("nan")
    print("%-34s eps=%.3f lam=%.3f p=%.2f  Lueders-of-averaged: %+.6e (hand %+.6e, ratio %.4f)   "
          "mixture: %+.6e (hand %+.6e, ratio %.4f)   [min eig E %.3f, <X> shift %.1e]" % (
              tag, eps, lam, r["p"], r["sig"], r["pred"], ratio, r["sig_mix"], r["pred_mix"], ratio_m,
              r["min_eig_E"], r["leak"]))


def main():
    print("N per mode = %d" % N)
    print("== controls ==")
    show("ZERO (theta = pi/2 only)", [np.pi / 2], [1.0], 0.2, 0.3)
    show("ABSOLUTE (theta = 0 only, p = 1)", [0.0], [1.0], 0.2, 0.3)
    show("ABSOLUTE, smaller eps, lam", [0.0], [1.0], 0.05, 0.1)
    print("== relativised, three-point density, p = 0.7 ==")
    th = [-0.9, 0.2, 1.3]
    ps = [0.2, 0.3, 0.2]
    for eps, lam in [(0.2, 0.3), (0.1, 0.2), (0.05, 0.1), (0.02, 0.05)]:
        show("three points", th, ps, eps, lam)
    print("== relativised, phase spread over a full turn (8 points, p = 0.8): a wide frame ==")
    th8 = [2 * np.pi * j / 8 for j in range(8)]
    ps8 = [0.1] * 8
    for eps, lam in [(0.2, 0.3), (0.05, 0.1)]:
        show("eight points, uniform", th8, ps8, eps, lam)
    print("== not small: eps = 0.45, lam = 1 ==")
    show("three points, large", th, ps, 0.45, 1.0)
    show("eight points, large", th8, ps8, 0.45, 1.0)


if __name__ == "__main__":
    main()
