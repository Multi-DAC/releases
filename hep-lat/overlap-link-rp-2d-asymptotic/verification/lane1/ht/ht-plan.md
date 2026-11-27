# H_T extension (draft section 8): sealed numeric check (clock reading 12:31:52, 2026-11-26; before the code)

Object: draft-section8-HT.md. Kernel H_T = gamma5 D_W(-m)(2 + D_W(-m))^{-1}, m = 0.7, a5 = 1. Same sampler,
conventions and m_f = 1/L as lane1_check.py. Dense cells: e l in {0.47, 2}, L in {8, 12, 16}, 10 samples each in C and
Q = 0, one Haar phi each. Seed 20261127.

Killers (deterministic; any violation means the draft is wrong):
- Q1: min |eig H_T0| >= 0.3793 at every sampled phi (the analytic free gap).
- Q2: ||H_T1 - H_T0||_op <= (4/(2 - m)^2) ||a||_inf and ||H_T1 - H_T0||_HS <= (2 sqrt2/(2 - m)^2) ||a||_2 (+1e-10).
- Q3: min |eig H_T1| >= min |eig H_T0| - (4/(2 - m)^2) ||a||_inf - 1e-10.
- Q4: wherever mu0 := 0.3793 - (4/(2 - m)^2) ||a||_inf > 0,
  ||D_T1 - D_T0||_HS <= (1 - m_f)(1/(2 mu0))(2 sqrt2/(2 - m)^2) ||a||_2.
- Q5: with c = 4/(25 pi^2), n_0(2 delta) <= 16 c1^2 delta^2 V/pi^2 + 4, and T_0(2 delta) <= V tau(delta) + delta^-2 at each
  valid delta in {0.05, 0.1, 0.2}, on the exact free spectrum at the sampled phi.
- Q6: Mirsky and (3.1) with the measured h, as in P2.

Guess (0.6), not a killer: at e l = 2, the mean paired difference ln det D_T1 - ln det D_T0 is negative, at between 0.5
and 2 times the H_W value in the same cell from out-lane1-2026-11-26-1226.log (-0.042 at L = 12, -0.069 at L = 16).
