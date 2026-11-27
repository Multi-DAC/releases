## T-c (sealed before writing the code; clock reading in the sha line below): <det D_ov>_W over Haar holonomies
- Setup: T-b's L x L torus, antiperiodic in both directions, m = 0.7, m_f = 0, all plaquettes flat. Flat fields are
  gauge-equivalent to constant links U_mu = e^{i phi_mu / L}, with holonomies e^{i phi_mu} Haar-uniform (they carry
  no Wilson action). In momentum space, k_mu = (2 pi n + pi + phi_mu)/L, D_W(k) = i s.gamma + b, with
  s_mu = sin k_mu, b = sum(1 - cos k_mu) - m, and per momentum det D_ov(k) = 2 (1 + b/w), w = sqrt(s^2 + b^2).
  So det D_ov(phi) = prod_k 2(1 + b/w) <= 4^V, zero at phi = (pi, pi) (k = 0 allowed).
- Quantities: <det>_W = (2 pi)^-2 int dphi0 dphi1 det(phi) (periodic trapezoid; the integrand is analytic, so
  convergence is checked n vs 2n), plus a Haar Monte Carlo with 4000 samples as a literal check.
  f_inf = int d^2k/(2 pi)^2 ln 2(1 + b/w), from an independent fine midpoint grid.
  c_f := -f_inf. Per-pair slopes sl(L1, L2) = (ln<det>(L2) - ln<det>(L1))/(V2 - V1). Remainders
  r0(L) = ln det(0) - V f_inf and R(L) = ln<det> - ln det(0). L in {8, 12, 16, 24, 32, 40}.
- C0 (control, must pass before anything is read): the momentum formula against the dense slogdet of T-b's builder
  (constant holonomy links) at L = 8 and 12, 3 random phi each: |diff| < 1e-8.
- P1 (0.85): c_f does NOT drift. Every per-pair slope with L1 >= 16 equals f_inf within 1e-4 per site.
  Killer K1: any such slope off f_inf by > 1e-3, or |r(40) - r(24)| > 0.1 for r = ln<det> - V f_inf.
- P2 (0.65): r0(L) -> ln|theta_3(i)/eta(i)|^2 = ln 2 = 0.6931 (the massless Dirac torus partition function at
  tau = i, antiperiodic-antiperiodic). |r0(40) - 0.6931| < 0.03.
- P3 (0.60): R(L) -> ln[ <|theta[a,b](0|i)|^2>_{a,b} / theta_3(i)^2 ] = ln[(1/sqrt 2)/(sqrt 2 eta(i))^2] = -0.5124
  (eta(i) = Gamma(1/4)/(2 pi^{3/4})). |R(40) + 0.5124| < 0.03.
- P4 (0.55, a guess): f_inf in [1.10, 1.35], so the worst-case loss per site 2 ln 2 - f_inf is in [0.04, 0.29].
- Reading: if P1 holds, (b') at flat plaquettes has c_f = -f_inf, L-independent, and the volume loss of
  sup(det)/<det> is V (2 ln 2 - f_inf) + O(1). Flat plaquettes are the beta -> infinity limit only.
