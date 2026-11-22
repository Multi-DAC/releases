"""Blind check of Wren's Claim 3(c) (letter 42): my own construction, written before seeing theirs.

Usage: python -I -B -X utf8 claim3c_check.py

Claim, non-commuting part: for every delta in (0, 1/2) and eta in (0, 1) there are a finite-dimensional
H_1 (x) C^2, a POVM {E_i}, and 0 <= B_i <= E_i with B_i = Btilde_i (x) 1_2, sum (E_i - B_i) <= delta, yet
   || sum (E_i^1/2 - B_i^1/2)^2 || >= 2 - eta      and      || Lambda(1 (x) Z) - 1 (x) Z || >= 1 - eta,
Lambda(X) = sum E_i^1/2 X E_i^1/2.

My construction. H_1 = span{e_1..e_N, w}. For i = 1..N:
   Btilde_i = beta |k_i><k_i|,  k_i = a e_i + b w  (a^2 + b^2 = 1, b small),
   E_i = Btilde_i (x) 1 + c |e_i><e_i| (x) |+><+|.
Outcome 0: Btilde_0 = 1 - sum Btilde_i - c P_N,  E_0 = 1 - sum_{i>=1} E_i = Btilde_0 (x) 1 + c P_N (x) |-><-|.
Then sum_i (E_i - B_i) = c P_N (x) 1, so delta = c. beta is the largest value keeping Btilde_0 >= 0.

Two routes:
  brute(N, c, b): every operator as a full matrix; norms by eigenvalues; every hypothesis checked numerically.
  reduced(N, c, b): the same quantities' LOWER bounds from the 4x4 block of one outcome and a 2x2 block of
      outcome 0, valid at any N:   column error >= <w,+| sum_{i>=1} T_i^2 |w,+>,
      disturbance >= |<w,+| (Lambda(Z) - Z) |w,->|.   brute() prints the same two numbers for comparison.
Control: b = 0 makes every B_i commute with its E_i; then the commuting benchmark of Claim 3(c) must hold,
      column error <= delta, and the script must report that.
"""
import numpy as np

Z2 = np.array([[1.0, 0.0], [0.0, -1.0]])
PLUS = np.array([1.0, 1.0]) / np.sqrt(2)
MINUS = np.array([1.0, -1.0]) / np.sqrt(2)


def psd_sqrt(M):
    w, v = np.linalg.eigh((M + M.T) / 2)
    return (v * np.sqrt(np.clip(w, 0, None))) @ v.T


def beta_max(N, c, b, slack=1e-9):
    """Largest beta with beta*K + c*P_N <= 1 - slack, K = sum |k_i><k_i|; decided in span{s-hat, w}."""
    a = np.sqrt(1 - b * b)

    def top(beta):
        M = np.array([[beta * a * a + c, beta * a * b * np.sqrt(N)], [beta * a * b * np.sqrt(N), beta * N * b * b]])
        return np.linalg.eigvalsh(M)[-1]

    lo, hi = 0.0, 1.0
    for _ in range(200):
        mid = (lo + hi) / 2
        if top(mid) < 1 - slack:
            lo = mid
        else:
            hi = mid
    return lo


def block(c, b, beta):
    """E_i, B_i for one outcome on span{e_i, w} (x) C^2; index = 2*h + q, h in (e_i, w), q in (up, down)."""
    a = np.sqrt(1 - b * b)
    k = np.array([a, b])
    Bt = beta * np.outer(k, k)
    B = np.kron(Bt, np.eye(2))
    D = c * np.kron(np.diag([1.0, 0.0]), np.outer(PLUS, PLUS))
    return B + D, B


def reduced(N, c, b):
    beta = beta_max(N, c, b)
    a = np.sqrt(1 - b * b)
    E, B = block(c, b, beta)
    S, R = psd_sqrt(E), psd_sqrt(B)
    T = S - R
    wp = np.kron(np.array([0.0, 1.0]), PLUS)
    wm = np.kron(np.array([0.0, 1.0]), MINUS)
    Zb = np.kron(np.eye(2), Z2)
    col = N * (wp @ T @ T @ wp)
    lam = N * (wp @ S @ Zb @ S @ wm)
    # outcome 0: E_0^1/2 = Btilde_0^1/2 (x) |+><+| + (1 - Btilde_tot)^1/2 (x) |-><-|; in span{s-hat, w}:
    M = np.array([[beta * a * a, beta * a * b * np.sqrt(N)], [beta * a * b * np.sqrt(N), beta * N * b * b]])
    r_plus = psd_sqrt(np.eye(2) - M - np.diag([c, 0.0]))      # on the |+> sector
    r_minus = psd_sqrt(np.eye(2) - M)                          # on the |-> sector
    # <w,+| E0^1/2 Z E0^1/2 |w,-> = sum_h r_plus[w,h] r_minus[h,w]   (Z maps |-> to |+>)
    lam += r_plus[1, :] @ r_minus[:, 1]
    return dict(beta=beta, c_over_beta=c / beta, g=beta * N * b * b, col_lb=col, dist_lb=abs(lam - 1.0))


def brute(N, c, b):
    beta = beta_max(N, c, b)
    a = np.sqrt(1 - b * b)
    n1 = N + 1
    dim = 2 * n1
    I = np.eye(dim)
    Zb = np.kron(np.eye(n1), Z2)
    sumE = np.zeros((dim, dim))
    sumD = np.zeros((dim, dim))
    col = np.zeros((dim, dim))
    lam = np.zeros((dim, dim))
    Bt_tot = np.zeros((n1, n1))
    worst_BE, worst_comm, tr_T2, tr_D = 0.0, 0.0, 0.0, 0.0
    for i in range(N):
        k = np.zeros(n1)
        k[i], k[N] = a, b
        Bt = beta * np.outer(k, k)
        e = np.zeros(n1)
        e[i] = 1.0
        B = np.kron(Bt, np.eye(2))
        E = B + c * np.kron(np.outer(e, e), np.outer(PLUS, PLUS))
        S, R = psd_sqrt(E), psd_sqrt(B)
        T = S - R
        sumE += E
        sumD += E - B
        col += T @ T
        lam += S @ Zb @ S
        Bt_tot += Bt
        if i < 5:
            worst_BE = min(worst_BE, np.linalg.eigvalsh(E - B)[0])
            worst_comm = max(worst_comm, abs(B @ Zb - Zb @ B).max())
        tr_T2 += np.trace(T @ T)
        tr_D += np.trace(E - B)
    PN = np.diag([1.0] * N + [0.0])
    Bt0 = np.eye(n1) - Bt_tot - c * PN
    B0 = np.kron(Bt0, np.eye(2))
    E0 = I - sumE
    S0, R0 = psd_sqrt(E0), psd_sqrt(B0)
    T0 = S0 - R0
    sumD += E0 - B0
    col += T0 @ T0
    lam += S0 @ Zb @ S0
    tr_T2 += np.trace(T0 @ T0)
    tr_D += np.trace(E0 - B0)
    wp = np.kron(np.eye(n1)[N], PLUS)
    wm = np.kron(np.eye(n1)[N], MINUS)
    return dict(
        beta=beta, c_over_beta=c / beta, dim=dim,
        min_eig_E0=np.linalg.eigvalsh(E0)[0], min_eig_Bt0=np.linalg.eigvalsh(Bt0)[0],
        min_eig_E0_minus_B0=np.linalg.eigvalsh(E0 - B0)[0], min_eig_Ei_minus_Bi=worst_BE,
        comm_B_Z=max(worst_comm, abs(B0 @ Zb - Zb @ B0).max()),
        sqrt_resid=abs(S0 @ S0 - E0).max(),
        deficit_norm=np.linalg.eigvalsh(sumD)[-1],
        col_norm=np.linalg.eigvalsh(col)[-1],
        dist_norm=abs(np.linalg.eigvalsh(lam - Zb)).max(),
        col_at_w=wp @ (col - T0 @ T0) @ wp, dist_at_w=abs(wp @ (lam - Zb) @ wm),
        trace_T2=tr_T2, trace_deficit=tr_D)


def main():
    print("== brute force, every hypothesis checked ==")
    for (N, c, b) in [(300, 0.3, 0.3), (600, 0.1, 0.3), (600, 0.3, 0.0)]:
        r = brute(N, c, b)
        tag = "CONTROL b=0 (commuting): column error must be <= deficit" if b == 0 else ""
        print("N=%d c=%.2f b=%.2f dim=%d  %s" % (N, c, b, r["dim"], tag))
        print("   beta=%.3e  c/beta=%.1f" % (r["beta"], r["c_over_beta"]))
        print("   min eig: E_0 %.2e, Btilde_0 %.2e, E_0-B_0 %.2e, E_i-B_i (first five) %.2e" % (
            r["min_eig_E0"], r["min_eig_Bt0"], r["min_eig_E0_minus_B0"], r["min_eig_Ei_minus_Bi"]))
        print("   max |[B_i, 1(x)Z]| %.2e   sqrt residual of E_0 %.2e" % (r["comm_B_Z"], r["sqrt_resid"]))
        print("   || sum (E_i - B_i) ||             = %.6f   (delta)" % r["deficit_norm"])
        print("   || sum (E_i^1/2 - B_i^1/2)^2 ||   = %.6f   (lower bound at |w,+>, outcome 0 dropped: %.6f)" % (
            r["col_norm"], r["col_at_w"]))
        print("   || Lambda(1(x)Z) - 1(x)Z ||       = %.6f   (lower bound |<w,+|.|w,->|: %.6f)" % (
            r["dist_norm"], r["dist_at_w"]))
        print("   traces: sum tr T_i^2 = %.4f <= sum tr (E_i - B_i) = %.4f  (Powers-Stormer)" % (
            r["trace_T2"], r["trace_deficit"]))
        if b > 0:
            q = reduced(N, c, b)
            print("   reduced route, same parameters: column >= %.6f, disturbance >= %.6f" % (q["col_lb"], q["dist_lb"]))
    print()
    print("== reduced route (lower bounds), toward the claimed 2 - eta and 1 - eta ==")
    print("  delta      b        N           c/beta     column >=    2-2b       disturbance >=   1-b")
    for c in (0.4, 0.1, 0.01):
        for b in (0.1, 0.025, 0.005):
            for K in (1e2, 1e4):
                N = int(np.ceil(K / (c * b * b)))
                q = reduced(N, c, b)
                print("  %-9.3g  %-7.3g  %-10.3g  %-9.3g  %.6f     %.4f     %.6f         %.4f" % (
                    c, b, N, q["c_over_beta"], q["col_lb"], 2 - 2 * b, q["dist_lb"], 1 - b))


if __name__ == "__main__":
    main()
