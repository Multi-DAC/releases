# LANE 1 numeric controls: sealed plan (clock reading 12:25:08, 2026-11-26; written before the code)

Object: projects/bounties/chiral-rp/LANE1-Z-LEMMA.md, Steps 1 to 4. The checks test the inequalities and the
normalisations on exact Wilson samples, not the cited theorems (Prekopa, Caffarelli, Mirsky).
Conventions: those of te_finite_beta.py (m = 0.7, antiperiodic torus, r = 1). eps0 = 0.5, eta = m/4 = 0.175,
mu l = 1 (m_f = 1/L). Sampler: T-e's exact rejection sampler, with links from T-e's Coulomb solve plus Haar phi.
Seeded rng 20261126.

## Controls (must pass before anything is read)
- C0: the plaquette reconstruction from (a, phi) matches the sampled theta to 1e-10 (Q = 0 samples).
- C1: Parseval. ||a||_2^2 computed from the link angles (sum of A_l^2 after removing phi/L) equals
  (1/V) sum_{k != 0} |F_k|^2/khat^2 to 1e-9, and d*a = 0 to 1e-10.
- C2: on every dense sample, D_{m_f} = (1 - m_f)(1 + gamma5 sign H)/2 + m_f has all singular values in
  [m_f - 1e-10, 1 + 1e-10], and det is real and positive (|Im ln det| < 1e-8 mod 2 pi).
- C3: the Jensen identity. The mean over 4000 uniform phi of (1/V) ln det D_0[phi/L], from T-c's momentum formula
  minus 2 ln 2 per site, equals -c_f = -0.3131624 within 3 se, at L = 8 and 16.

## Predictions

P1, deterministic: every dense sample with mu0 := m - 2||a||_inf > 0. Killer if any sample violates it.
- (a) min |eig H1| >= mu0 - 1e-10.
- (b) ||H1 - H0||_HS^2 = 2 sum_l |e^{i a_l} - 1|^2 to 1e-9.
- (c) ||D1 - D0||_HS <= (1 - m_f) sqrt2 ||a||_2/(2 mu0) (+1e-10), with D1, D0 at the same phi.
- Tightness guess (0.6): the ratio in (c) lies in [0.1, 0.7] on every sample.

P2, deterministic: every dense sample. Killer if any sample violates it.
- (a) Mirsky: sum_i (sigma_i^1 - sigma_i^0)^2 <= ||D1 - D0||_F^2.
- (b) (3.1) holds with the measured h and the exact free n_0(2 delta), T_0(2 delta), for delta in {0.05, 0.1, 0.2}.
- (c) The bounds (3.2) and (3.3) are >= the exact n_0 and T_0 at each valid delta (delta >= sqrt2 pi/(c1 L)).
- Size guess (0.6), not a killer: the slack of (3.1) at delta = 0.1 is above 10 on every sample, so the bound is loose.

P3, MC, band from the sample se: the Coulomb norm under nu (samples in C and Q = 0).
- r := mean ||a||_2^2 / (E_vM[theta^2] V G(L)), with E_vM[theta^2] by quadrature at the cell's beta.
- Prediction (0.75): |r - 1| <= 3 se(r) in every cell.
- se, sized before sealing from the Gaussian formula: relsd per sample is 0.25, 0.18, 0.16, 0.15 at L = 8, 16, 24, 32.
  With N = 400 per cell, 3 se is about 0.022 to 0.038.
- Killer K3 (the lemma's bound): mean ||a||_2^2 > (e l)^2 G(L)/cos eps0 + 3 se in any cell.
- Cells: e l in {0.47, 1, 2}, L in {8, 16, 24, 32}. Cells with fewer than 100 samples in C and Q = 0 are reported,
  not scored.

P4, MC: the joint event, in the same cells, as fractions of ALL samples.
- (a) Measured P(E) >= the proved bound (1.4) - 3 binomial se wherever (1.4) is non-vacuous. Killer if violated.
- (b) Measured P(Q = 0 | C) >= tanh(a_Q/2) - 3 se. Killer if violated. The sharpest cell is e l = 2: bound 0.974,
  discrete-Gaussian truth about 0.985, se about 0.006 at N = 400.
- (c) Guess (0.7): measured P(E) >= 0.95 at e l = 0.47 for every L >= 16.
- (d) Guess (0.6): at e l = 0.47 no sample has ||a||_inf > eta, and P(C) = 1 in every cell.

P5, dense, a guess: the paired difference ln det D1 - ln det D0 (same phi), averaged over the dense samples.
- (0.6): |mean| <= 0.5 in every dense cell. At e l = 2 it is negative by more than 2 se in at least two of the
  L values.
- Dense cells: e l in {0.47, 2}, L in {8, 12, 16, 24}, 20 samples each (C and Q = 0), one uniform phi each.

## Reading
- If P1, P2 and K3 hold and C0 to C3 pass, then Steps 1 to 4 have no normalisation or factor error visible at these
  sizes. That is a control, not a proof check.
- A P1 or P2 violation is a proof error, and the lemma is withdrawn until it is repaired.
