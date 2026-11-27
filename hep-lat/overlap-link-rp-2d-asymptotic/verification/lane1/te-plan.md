## T-e (sealed before writing the code): the normalisation at finite beta, at fixed physical size
- Question: T-c's (b') used flat plaquettes. At fixed e l, beta = (L/(e l))^2, so plaquette angles are of order e l/L. Is
  the loss to <det D_ov>_W from fluctuations O(1), or extensive?
- Setup: T-b's L x L antiperiodic torus, m = 0.7, overlap with the H_W kernel (dense eigh). Gauge fields are EXACT
  samples of the Wilson measure: V - 1 plaquette angles iid von Mises(beta), the last fixed by the torus constraint,
  accepted with probability e^{beta (cos theta_V - 1)}. Q = sum theta_p / 2pi. Links: a constant-flux field for Q,
  plus a Fourier (Coulomb-gauge) solution of curl A = theta - 2pi Q/V, plus Haar holonomies. The plaquette angles are
  checked against the targets to 1e-10 (C0).
- det D_ov = 0 when Q != 0 (index), so <det>_W = P(Q=0) <det | Q=0>, and by Jensen
  ln<det>_W >= ln P(Q=0) + <ln det | Q=0>.
- Measured: Delta(L) = <ln det D_ov | Q=0>_W - <ln det D_ov[flat]>_Haar, where the flat term is V f_inf + r_ln(L), with r_ln
  computed by the T-c momentum formula (Haar average of ln det, n = 64 grid). P(Q=0) as a fraction. Grid:
  e l in {2, 4}, L in {8, 12, 16, 20, 24}, beta = (L/(e l))^2, 40 accepted samples per cell, seeded.
- C1: one sample with theta_p = 0 and holonomy (0,0) reproduces T-c's ln det(0) at L = 8 to 1e-9.
- P1 (0.70): no extensive loss. At each e l, Delta(L) stays bounded. The fitted slope of Delta on V over L >= 12 has
  magnitude < 2e-3 per site, and |Delta(24) - Delta(16)| < 0.5 + 3 se.
  Killer K1: a slope above 5e-3 per site, with Delta(24) more negative than Delta(12) by > 1.0 + 3 se.
- P2 (0.60): Delta < 0 in every cell, and Delta(e l = 4)/Delta(e l = 2) is in [2, 8] at L = 24 (quadratic in e l gives 4).
- P3 (0.65): P(Q=0) is L-independent within 3 se at each e l, near the discrete-Gaussian value with variance
  (e l)^2/(4 pi^2): 0.63 at e l = 4, 0.99 at e l = 2 (within 0.1).
- Reading: if P1 holds, (b') survives finite beta at fixed e l, with the flat c_f plus O(1). The same then follows for
  P = det D_N^(S), given T-d's geometric convergence in S. That second step is an inference, not tested here.
