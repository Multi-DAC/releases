# LANE 1: the finite-beta lower bound on Z at fixed physical fermion mass (Clawd; started after clock reading 12:09, 2026-11-26)

Status: SIGNED OFF by Wren (letter 89, 13:01): sections 0 to 7 HOLD as written. Section 8 held at m = 0.7 and needed
five small fixes for all 0 < m <= 1; they are now in, together with her sharpening of the statement (the
V ln(1 - m_f^2) term is gone). Before that, one refuter round on sections 0 to 7 returned HOLDS, with four minor
wounds, all repaired. A second refuter round, on Section 8 (23:56, 11-26), left the rate standing and killed one
constant: eps = 0.267 sat past the zero of Adams' bound, so the ceilings are now the open bounds e l < 0.4727 (H_T)
and e l < 0.3367 (H_W), and L_0(H_T) is 32. Section 9 has all three records. The numeric controls (sealed plan
`verification/lane1/lane1-plan.md`) test the inequalities and the normalisations on exact Wilson
samples, not the theorems cited.
Split and lead: Wren letter 87 part 3.2(b), my letter 289, Wren letter 88 (fixed PHYSICAL mass mu > 0, so m_f = mu a).

## 0. Setup and notation

- Lattice: the L x L torus, V = L^2 sites, L >= 3. Links U_l = e^{i A_l}, plaquette angles theta_p in (-pi, pi]
  (principal value of U_p).
- Wilson measure: mu_W(dU) proportional to prod_p e^{beta cos theta_p} prod_l dU_l, with beta = V/(e l)^2 (fixed
  physical size l = L a, fixed e).
- Kernel: H = H_W[U] = gamma5 (D_W[U] - m), with r = 1, 0 < m <= 1 (m = 0.7 in every number below), and the
  conventions of T-b/T-e (`te_finite_beta.py`). Fermion boundary conditions: any fixed twist (antiperiodic in T-c/T-e);
  the result does not depend on it (Step 4).
- Normalised massive overlap, Wren's normalisation (letter 88): D_0[U] = (1/2)(1 + gamma5 sign H),
  D[U] = D_{m_f}[U] = (1 - m_f) D_0[U] + m_f.
  - W := gamma5 sign H is unitary, so D_0 = (1 + W)/2 and D are normal, with eigenvalues (1 - m_f)(1 + w)/2 + m_f,
    |w| = 1.
  - So every singular value of D lies in [m_f, 1], and sigma^2 = (1 - m_f^2)|lambda_0|^2 + m_f^2 eigenvalue by
    eigenvalue, where lambda_0 = (1 + w)/2.
  - gamma5 D gamma5 = D^dag, so det D is real. Its real eigenvalues are 1 and m_f, so det D[U] > 0 at m_f > 0 for
    every U where H_W[U] is invertible (sign H_W needs that). The set {det H_W[U] = 0} is the zero set of a
    real-analytic function that is not identically zero. So it is mu_W-null, and it does not meet E, where the gap is
    >= m/2 (Wren, letter 89).
  - Also |lambda_0| <= 1, so sigma^2 - |lambda_0|^2 = m_f^2 (1 - |lambda_0|^2) >= 0. Each singular value of D_{m_f}
    is at least the matching |lambda_0| (Wren, letter 89; used in section 5).
- Fixed physical mass: m_f = mu a = mu l / L, with mu l > 0 fixed, so m_f -> 0 and V m_f^2 = (mu l)^2 is O(1).
- Free constant, H_W kernel: c_f(m) := -(1/(4 pi^2)) int_{T^2} d^2k ln[(1 + b/w)/2], with
  b = sum_mu (1 - cos k_mu) - m, s^2 = sum_mu sin^2 k_mu, w = sqrt(s^2 + b^2).
  c_f(0.7) = 2 ln 2 - f_inf = 0.3131624 (T-c, f_inf = 1.0731320).
- Z := E_W[det D_{m_f}[U]].

## 1. Statement

**Lemma Z.** Fix 0 < m <= 1, e l > 0, mu l > 0 and eps0 in (0, pi/2). Then there is an explicit L_0, and for every
L >= L_0,

    ln Z >= -c_f(m) V - Lambda(L) + ln P_E,

where P_E >= (1 - eta_C)(tanh(a_Q/2))(2/3 - eta_inf) > 0 (Step 1) and Lambda(L) = O(L (ln L)^{3/2}) (Step 3, (3.6)).
So

    (1/V) ln Z >= -c_f(m) - O(L^{-1} (ln L)^{3/2}).

The first version carried an extra + V ln(1 - m_f^2), which is >= -(mu l)^2/(1 - m_f^2). Wren's sharpening
sigma >= |lambda_0| (section 0) removes it.

The lemma holds at EVERY fixed e l. The small box enters only in its use (Section 6).

## 2. Step 1: the event E, and its probability JOINTLY

E is defined on the plaquette angles only (Step 4 shows no holonomy window is needed at m_f > 0):

    E = C ∩ {Q = 0} ∩ {||a||_2^2 <= A*} ∩ {||a||_inf <= eta},

- C = {max_p |theta_p| <= eps0}, the plaquette cap.
- Q = (1/2pi) sum_p theta_p, an integer.
- On Q = 0, a = a(theta) is the Coulomb potential: the unique real 1-form with da = theta (as real numbers),
  d*a = 0 and zero mean. a = d* G theta, with G = (d d*)^+ on 2-forms. In Fourier, ||a||_2^2 = <theta, G theta> =
  sum_{k != 0} |theta^_k|^2 / khat^2, with khat^2 = sum_mu 4 sin^2(k_mu/2) and theta^_k = V^{-1/2} sum_p e^{-ikp} theta_p.
- A* = 3 (e l)^2 G(L)/cos eps0 and eta = m/4, where G(L) := (1/V) sum_{k != 0} 1/khat^2 (an explicit finite sum).

**1a. Representation.** The map U -> (U_p) is a homomorphism onto {prod_p U_p = 1}, with kernel the flat fields.
Modulo gauge, the flat fields are a 2-torus of holonomies.
- So under Haar, with any measurable section, the configuration is (theta, phi), with phi uniform on T^2 and
  independent of theta. The Wilson weight depends on theta alone, so phi stays Haar and independent under mu_W.
- On Q = 0 the section is A = a(theta) + phi/L (constant phi_mu/L on mu-links). e^{ia} has plaquette angles exactly
  theta, and U e^{-ia} is flat, hence gauge-equivalent to the constant holonomy field.
- Law of theta under mu_W: iid von Mises(beta), conditioned on sum_p theta_p in 2 pi Z.
  - Equivalently, theta_1..theta_{V-1} have density proportional to prod_{p<V} e^{beta cos theta_p} e^{beta cos S'},
    where S' = sum_{p<V} theta_p. This is because U_V = (prod_{p<V} U_p)^{-1}, so cos theta_V = cos S'.
  - This law is exchangeable in all V plaquettes.

**1b. The cap.** For any event F in theta_1..theta_{V-1},

    mu_W(F) = E_0[1_F e^{beta(cos S' - 1)}] / E_0[e^{beta(cos S' - 1)}] <= P_0(F) / (e^{-beta} I_0(beta)),

where E_0 is iid vM(beta) on V - 1 plaquettes. The step uses e^{beta cos x} = sum_n I_n(beta) e^{inx} and
E_0 e^{inS'} = (I_n/I_0)^{V-1} >= 0, so E_0 e^{beta cos S'} >= I_0.
- Also e^{-beta} I_0(beta) = (1/pi) int_0^pi e^{-beta(1 - cos t)} dt >= (1/pi) int_0^pi e^{-beta t^2/2} dt
  = erf(pi sqrt(beta/2))/sqrt(2 pi beta).
- And P_vM(|theta| > eps0) <= (1 - eps0/pi) e^{-beta(1 - cos eps0)} / (e^{-beta} I_0(beta)).
- By exchangeability the bound for theta_1 holds for every p, including p = V. So

      mu_W(not C) <= eta_C := V (1 - eps0/pi) 2 pi beta e^{-beta(1 - cos eps0)} / erf(pi sqrt(beta/2))^2.     (1.1)

**1c. Q = 0, uniformly in L.** Under mu_W(. | C), theta is iid capped vM conditioned on S = sum theta in 2 pi Z, so
mu_W(Q = q | C) = f(2 pi q) / sum_n f(2 pi n), with f the density of S under iid capped vM.
- Capped vM is strongly log-concave: e^{beta cos t} 1_{|t| <= eps0} = g(t) e^{-kappa t^2/2}, with
  kappa = beta cos eps0 and g log-concave and even. This holds because (beta cos t + kappa t^2/2)'' =
  beta(cos eps0 - cos t) <= 0 on |t| <= eps0.
- Write theta = (x/V) 1 + theta_perp, with theta_perp in H = {sum = 0}. Then

      f(x) proportional to e^{-kappa x^2/(2V)} J(x),   J(x) = int_H prod_p g(x/V + theta_perp,p) e^{-kappa |theta_perp|^2/2} dtheta_perp.

  J is log-concave (Prekopa) and even, so J(x) <= J(0), and f(2 pi q)/f(0) <= e^{-a_Q q^2}, with
  a_Q := kappa (2pi)^2/(2V) = 2 pi^2 cos eps0/(e l)^2.
- Hence, with no dependence on L,

      mu_W(Q = 0 | C) >= 1 / sum_n e^{-a_Q n^2} >= tanh(a_Q/2).                                               (1.2)

**1d. The Coulomb norms under nu := mu_W(. | C, Q = 0).**
- nu is the slice measure on H ∩ [-eps0, eps0]^V, with density prod_p e^{beta cos theta_p}. Its potential has Hessian
  >= kappa on a convex domain, and it is even.
- Caffarelli's contraction theorem, used only with a FINITE convex perturbation, so no convex-support version is
  needed.
  - Let psi(t) = -beta cos t on |t| <= eps0. Outside, take its C^1 quadratic continuation,
    psi(t) = -beta cos eps0 + beta sin eps0 (|t| - eps0) + (kappa/2)(|t| - eps0)^2, so psi'' >= kappa everywhere.
  - Let nu_n be the measure on H with density proportional to
    exp(-sum_p psi(theta_p) - n sum_p (|theta_p| - eps0)_+^3).
    The cube makes the penalty C^2 and convex (Wren's polish, letter 89). psi is C^2 already, since its second
    derivative is beta cos eps0 = kappa on both sides of eps0.
    Its potential minus kappa |theta|^2/2 is a finite convex function on H.
  - So Caffarelli gives nu_n = T_n # N(0, kappa^{-1} I_H) with T_n 1-Lipschitz. Hence Cov_{nu_n} <= kappa^{-1} I_H
    (Gaussian Poincare composed with T_n), and nu_n(f - E_{nu_n} f >= t) <= e^{-kappa t^2/(2 Lip(f)^2)}.
  - As n -> infinity, nu_n -> nu weakly, by dominated convergence under the envelope e^{-kappa |theta|^2/2}.
    The tails are uniformly Gaussian, so second moments converge.
  - So Cov_nu <= kappa^{-1} I_H. For the open sets {|a_l| > eta}, nu <= liminf nu_n (portmanteau). The two
    consequences below hold for nu.
- (l2) E_nu theta = 0, so E_nu ||a||^2 = tr(G Cov_nu) <= tr G / kappa = V G(L)/(beta cos eps0) = (e l)^2 G(L)/cos eps0.
  By Markov, nu(||a||_2^2 > A*) <= 1/3.
- (sup) a_l = <G d e_l, theta> is linear, with Lip^2 = <d e_l, G^2 d e_l> = (1/V) sum_{k != 0} 4 sin^2(k_nu/2)/khat^4
  = G(L)/2 by the square symmetry. E_nu a_l = 0. Taking a union over the 2V links,

      nu(||a||_inf > eta) <= eta_inf := 4V exp(-kappa eta^2/G(L)) = 4V exp(-V cos eps0 eta^2/((e l)^2 G(L))).   (1.3)

**1e. Jointly.** mu_W(E) = mu_W(C) mu_W(Q = 0 | C) nu(l2 and sup), so

    P_E := mu_W(E) >= (1 - eta_C) tanh(a_Q/2) (2/3 - eta_inf).                                                (1.4)

The union bound is used only inside nu, where both pieces are small or Markov. Conditioning on C before Q is what
makes the Q bound L-uniform. A direct change of measure loses a factor sqrt(2 pi beta) ~ L, so it cannot bound
P(Q != 0) (checked by hand: P_0(|S'| >= pi) does not decay in L).

**Explicit bound on G(L).** For k in (-pi, pi]^2, khat^2 >= 4|k|^2/pi^2. The shell |n|_inf = j has 8j points with
|n|^2 >= j^2, so G(L) <= (1/16) sum_{0 < |n|_inf <= L/2} |n|^{-2} <= (1 + ln(L/2))/2. Exact values (out-consts.log):
G(L) = ln L/(2 pi) + 0.04876 + o(1), for example G(16) = 0.4899 and G(64) = 0.7107.

Numbers (m = 0.7, eps0 = 0.5, eta = 0.175; out-consts.log, with eta_inf evaluated using the proved bound on G):

| e l | tanh(a_Q/2) | eta_C at L = 12 / 16 / 32 | eta_inf at L = 16 / 32 / 64 / 128 |
|---|---|---|---|
| 0.34 | 1 - 2e-65 | 6e-61 / 6e-112 / 0 | 1.7e-14 / 6e-52 / 1e-181 / 0 |
| 0.47 | 1 - 2e-34 | 1.1e-29 / 3.8e-56 / 9e-240 | 1.7e-6 / 8e-26 / 2e-93 / 0 |
| 2 | 0.974 | vacuous / vacuous / 3.4e-8 | vacuous / vacuous / 0.073 / 1.9e-14 |
| 4 | 0.494 | vacuous to L = 32 | vacuous to L = 128 |

So at e l <= 0.47, P_E >= 2/3 - 2e-6 for every L >= 16. Both tails decrease in L from L = 16 on.
- ln eta_C = 2 ln V - V(1 - cos eps0)/(e l)^2 + O(1), where the 1/erf^2 factor decreases. This falls once
  V > 2 (e l)^2/(1 - cos eps0), which is V > 3.7 at e l = 0.47.
- ln eta_inf = ln 4V - V cos eps0 eta^2/((e l)^2 G(L)), and V/G(L) grows faster than ln V. The exponent at L = 16 is
  already 20. At e l = 2 the proved bound needs L >~ 64, and at e l = 4,
L >~ 200: L_0 grows with e l, but P_E stays Omega(1) at every fixed e l.

## 3. Step 2: perturbation, on E, uniformly in phi

Let D1 = D[a + phi/L] and D0 = D[phi/L] (free, holonomy phi), with H1 and H0 their kernels.
- Free gap: H0^2 = s^2 + b^2 >= m^2 + 2 U(1 - m) >= m^2, for every real k (so every phi), where U = sum(1 - cos k_mu).
  This uses sin^2 = u(2 - u) with u = 1 - cos, and sum u^2 <= U^2.
- Hopping difference: D_W[A1] - D_W[A0] = -sum_{x,mu} [P_-^mu Delta_{(x,mu)} psi(x + mu) + P_+^mu conj(Delta_{(x-mu,mu)}) psi(x - mu)],
  with Delta_l = U1_l - U0_l, |Delta_l| = |e^{i a_l} - 1| <= |a_l|, and P_-^mu + P_+^mu = 1 orthogonal projectors.
  - Operator norm: per direction, (P_- X + P_+ Y)^dag (P_- X + P_+ Y) = P_- X^dag X + P_+ Y^dag Y, which is <=
    max|Delta|^2. So ||H1 - H0||_op <= 2 ||a||_inf <= m/2 on E.
  - Hence the gap of H1 is >= m - m/2, and both kernels have spectra outside (-mu0, mu0), with mu0 = m/2.
  - HS norm: for L >= 3, each ordered neighbour pair carries one term, so ||H1 - H0||_HS^2 = 2 sum_l |Delta_l|^2 <=
    2 ||a||_2^2.
- Lipschitz in HS: for Hermitian A, B and Lipschitz f, ||f(A) - f(B)||_HS <= Lip(f) ||A - B||_HS. This is the
  double sum sum_ij |f(a_i) - f(b_j)|^2 |<u_i, v_j>|^2.
  - Take f = sign off (-mu0, mu0), linear inside, so Lip = 1/mu0 = 2/m. Then

        ||D1 - D0||_HS = ((1 - m_f)/2) ||sign H1 - sign H0||_HS <= sqrt2 ||a||_2 / m =: h <= h_E := sqrt(6 (e l)^2 G(L)/cos eps0) / m.   (2.1)

  - h_E^2 = O(ln L). At e l = 0.47, L = 1024: h_E = 1.88.

## 4. Step 3: Mirsky with a mass floor (Wren's three groups, constants explicit)

- Mirsky: sum_i (sigma_i(D1) - sigma_i(D0))^2 <= ||D1 - D0||_F^2 <= h^2, with both sorted decreasingly, N = 2V.
- Write Delta_i = sigma_i^1 - sigma_i^0. Every sigma lies in [m_f, 1]. Fix delta in (0, 1/4].
- Moved, B = {|Delta_i| > delta}: |B| <= h^2/delta^2, and each costs at most ln(1/m_f).
- Small, S = {i not in B, sigma_i^0 < 2 delta}: each costs at most ln(1/m_f), and |S| <= n_0(2 delta).
- Good, the rest: sigma^1/sigma^0 >= 1/2, and ln(1 - x) >= -2 ln2 x on [0, 1/2]. So the good group costs at most
  2 ln2 sum |Delta_i|/sigma_i^0 <= 2 ln2 h T_0(2 delta)^{1/2}, with T_0(t) := sum_{sigma_i^0 >= t} (sigma_i^0)^{-2}.

      ln det D1 - ln det D0 >= -[h^2/delta^2 + n_0(2 delta)] ln(1/m_f) - 2 ln2 h T_0(2 delta)^{1/2}.               (3.1)

Free counts, uniform in phi. Let C2 = 1/(1 - m_f^2) and c1 = pi sqrt((2 + m^2) C2).
- If b >= 0, then |lambda_0|^2 = (1 + b/w)/2 >= 1/2.
- If b < 0, then |k_mu| < pi/2 (because m <= 1), and |lambda_0|^2 = s^2/(2w(w - b)) >= s^2/(4(s^2 + m^2)) >=
  |k|^2/(pi^2(2 + m^2)).
- So for t <= 1/2 and m_f^2 < 1/2: sigma^0(k) < t forces |k| < c1 t. Each k carries 2 singular values, and a disc of
  radius r holds at most (rL/pi + 1)^2 grid momenta. Hence

      n_0(2 delta) <= 16 (2 + m^2) C2 delta^2 V + 4.                                                                  (3.2)

- Split T_0(2 delta) at |k| = 2 c1 delta (Wren's expository line, letter 89).
  - Inner disc: at most the (3.2) count of singular values, each >= 2 delta. So it contributes at most
    (16 (2 + m^2) C2 delta^2 V + 4)/(4 delta^2) = 4 C2 (2 + m^2) V + delta^{-2}. That is the second term of tau and
    the + delta^{-2}.
  - Outer cells, in the annulus c1 delta <= |q| <= sqrt2 pi (1 + 1/L): these give the 4 C2 and log terms, as follows.
- With sigma^{-2} <= C2 [2 + pi^2 (2 + m^2)/|k|^2], the cell comparison (|q| <= (3/2)|k| on the cell of any grid k
  with |k| >= 2 sqrt2 pi/L) gives, for delta >= sqrt2 pi/(c1 L),

      T_0(2 delta) <= V tau(delta) + delta^{-2},
      tau(delta) = 4 C2 + 4 C2 (2 + m^2) + (9 pi/4) C2 (2 + m^2) ln(sqrt2 pi (1 + 1/L)/(c1 delta)).                   (3.3)

  At m = 0.7 and small m_f: tau ~ 14.0 + 17.6 ln(0.896/delta).
- Choose delta^4 = h_E^2/(16 (2 + m^2) C2 V), which balances the moved and small groups. Then

      Lambda(L) := 2 sqrt(16 (2 + m^2) C2) h_E L ln(1/m_f) + 4 ln(1/m_f) + 2 ln2 h_E (V tau(delta) + delta^{-2})^{1/2}   (3.6)
               = O(L (ln L)^{3/2}),

  and on E, for every phi, ln det D1 >= ln det D0 - Lambda(L).
- Numbers at m = 0.7, e l = 0.47, mu l = 1:
  - L = 256: Lambda/V = 0.54. L = 1024: delta = 0.017, Lambda/V = 0.18. L = 10^4: Lambda/V = 0.028
    (`verification/lane1/lambda_bound.py`, out-lambda.log).
- The rate is honest but slow: Lambda/V is above c_f = 0.313 until L ~ 500. The constants are crude (the small and moved
  groups are charged the full ln(1/m_f)). T-e's measured Delta is O(1) at L <= 24.

## 5. Step 4 and assembly: Jensen over the Haar holonomy (why no window is needed)

- Free part, exact. For fixed n, k_n = (2 pi n + chi + phi)/L, with chi the boundary twist, is uniform on one cell
  when phi is uniform on T^2, and the L^2 cells tile T^2 once. With g(k) = ln|lambda_0(k)|^2 in L^1:

      E_phi ln det D_0[phi/L] = sum_n E_phi g(k_n) = (V/(4 pi^2)) int_{T^2} g = -c_f V   EXACTLY, for any twist chi.

  The periodic-holonomy zero mode is an integrable log singularity and costs nothing in this average.
- Mass. sigma >= |lambda_0| singular value by singular value (section 0), so ln det D_{m_f} >= ln|det D_0| for every
  field. The first version used sigma^2 >= (1 - m_f^2)|lambda_0|^2 and paid V ln(1 - m_f^2). Step 3's C2 could be 1
  by the same inequality; it is left as is, which is valid and slightly loose.
- Assembly. phi is independent of theta and Haar (1a). Exp is convex (Jensen). On E, (3.1) holds uniformly in phi.
  Hence

      Z >= E_W[det D; E] = int_E mu_W(dtheta) E_phi det D1 >= int_E mu_W(dtheta) exp(E_phi ln det D1)
        >= P_E exp(-c_f V - Lambda(L)).

  That is Lemma Z.
- L_0 requires: L >= 3; m_f = mu l/L <= 1/2; delta <= 1/4 and delta >= sqrt2 pi/(c1 L); eta_C < 1 and
  eta_inf < 2/3.

## 6. Use: the small-box ceiling

- Bad set: not G = {some p with |1 - U_p| >= eps}, that is |theta_p| >= theta_eps = 2 arcsin(eps/2), so
  1 - cos theta_p >= eps^2/2.
- By 1b, mu_W(not G) <= V 2 pi beta e^{-beta eps^2/2}/erf(pi sqrt(beta/2))^2, and sup det D <= 1. So

      sup_{not G} det x mu_W(not G) / Z <= exp(V [c_f - eps^2/(2 (e l)^2)] + Lambda(L) + O(ln V)),

  which tends to 0 if e l < eps/sqrt(2 c_f).
  - That is the ceiling of THIS method, which charges the bad set sup det <= 1. It is not a physical threshold.
  - With eps < eps* = 0.266533 and c_f = 0.31316, the H_W ceiling is the open bound e l < 0.3367 (eps*/sqrt(2 c_f) =
    0.336784), below Wren's 0.34. CORRECTED 23:57, 11-26: the earlier "eps = 0.267, e l_0 = 0.337" rounded eps UP past
    the zero of the bound it feeds (refuter round 2, section 9).
  - eps* comes from Neuberger's bound (hep-lat/9911004, PRD 61, 085015; v2.tex L592-618, read first-hand 11-27 by
    Wren, letter 95, and by me, the arXiv TeX source v2.tex):
    |H(m)| >= sqrt(1 - (2 + sqrt2) sum_{mu > nu} eps_{mu nu}) - |1 - m| under ||1 - U_{mu nu}(x)|| <= eps_{mu nu}, in
    every dimension ("holds in both even and odd dimensions"). At d = 2 there is one plane, so c_1 = 2 + sqrt2, positive
    iff eps < eps* = m(2 - m)/c_1 (= 0.266533 at m = 0.7). Adams (hep-lat/9907005, L238-246) cites it and extends it to
    the doubler regions. The earlier line called it "Adams' generalisation of a 4d derivation"; it is Neuberger's own,
    in general d.
  - So e l_0 moves with eps. The refuter's uniform-flux probe (L = 20, |1 - U_p| up to 0.50) never brought
    min|eig H_W| below 0.458. That is no worst-case check.
- The 0.47 ceiling (H_T kernel at a5 = 1, c_f = 0.158939) is section 8.

## 7. Against Wren's lead (letter 87, 3.2(b)), and what is left

- Same skeleton: an Omega(1) event, the HS perturbation, Hoffman-Wielandt/Mirsky with a mass floor, three groups.
  Her remainder terms are recovered with explicit constants.
- Differences, all mine, each held under one refuter round:
  1. No holonomy window. Jensen over Haar phi makes the free term exact (Step 4).
  2. The H_W gap on E comes from a sup cap on the Coulomb field (||a||_inf <= m/4) plus the free gap m, not from the
     admissibility bound. So eps0 only has to be < pi/2, for strong log-concavity.
  3. P(Q = 0) is bounded uniformly in L by Prekopa on the capped law (1c). Her "P(Q = 0) >= c(e l)" was assumed.
  4. E_nu ||a||^2 is bounded under the true (capped, Q = 0) law by Caffarelli/Brascamp-Lieb, not by Gaussian
     analogy.
- Cited, not reproved: Prekopa (marginals of log-concave); Caffarelli's contraction for e^{-Phi} gamma with Phi
  convex and finite (1d now needs only this case); Mirsky's singular-value inequality. All are standard. None has been
  read first-hand in this lane.
- Not covered: m_f = 0 (the floor ln(1/m_f) is used in two groups); P = det D_N^(S) at finite S (route iii's
  replacement); a better rate than L (ln L)^{3/2}.

## 8. The H_T kernel at a5 = 1 (the 0.47 ceiling): Steps 2 and 3 with new constants (checked numerically, sealed; refuter round 2: rate holds, one constant corrected)

Kernel: H_T := gamma5 D_W(-m) (2 + D_W(-m))^{-1}.
- Hermitian: the two factors commute, and D_W is gamma5-hermitian.
- Free: |H_T(k)|^2 = (b^2 + s^2)/((2 + b)^2 + s^2), and gamma5 H_T(k) = b_T + i s_T.gamma, with
  b_T = (b(2 + b) + s^2)/n, s_T = 2s/n, n = (2 + b)^2 + s^2.
- c_f(H_T) = 0.158939: even midpoint grids N = 512 to 4096, O(N^-2) Richardson limit 0.15893896. This is Wren's
  0.1589390 (letter 89), and T-d's 512^2 value was 0.158936.
  - The first value, 0.158943, came from an odd n = 4001 midpoint grid. That grid puts a point at |k| = 3e-17, where
    ln|lambda0|^2 = -75, and the point alone adds 4.7e-6 (`verification/lane1/wren09/out-check.log`).
- e l_0 = eps*/sqrt(2 x 0.158939) = 0.472738, so the ceiling is the OPEN bound e l < 0.4727 (rounded down). CORRECTED
  23:57, 11-26 (refuter round 2, section 9): the earlier 0.4735 (and 0.4736 in letter 290) used eps = 0.267, which is past
  eps* = m(2 - m)/(2 + sqrt2) = 0.266533, the zero of Adams' bound at m = 0.7 (at 0.267 the bound is -0.00267). Any
  e l < 0.4727 works: take eps in (e l sqrt(2 c_f), eps*); the gap on G is then positive but small near the top.
- Adams' eps < eps* also gaps H_T (Wren, letter 89).
  - ||2 + D_W(-m)||_op <= (4 - m) + 2 = 6 - m for every U.
  - So |H_T| = |H_W (2 + D_W(-m))^{-1}| >= |H_W|/(6 - m) > 0 on the same admissible set |1 - U_p| < eps.

(a) Free gap, every holonomy, 0 < m <= 1.
- Put X = b^2 + s^2 >= m^2 + 2U(1 - m) (section 3), and (2 + b)^2 + s^2 = X + 4(1 - m + U).
- So |H_T|^2 >= 1/(1 + f(U)), with f(U) = 4(1 - m + U)/(m^2 + 2U(1 - m)) and U in [0, 4].
- f' has the sign of m^2 - 2(1 - m)^2. So f is monotone on [0, 4], and f <= max(f(0), f(4)). That gives
  mu_T(m) := min( m/(2 - m), sqrt((m^2 + 8(1 - m))/(m^2 + 8(1 - m) + 4(5 - m))) ).
  - The first version used f(4) only. That is wrong for m < 2 - sqrt2 = 0.586, where f decreases: at m = 0.5,
    min|H_T| = 1/3 against f(4)'s 0.437 (Wren, letter 89; reproduced in out-check.log, where min|H_T| = m/(2 - m)
    exactly at m = 0.3 and 0.5).
- At m = 0.7, f increases (0.49 > 0.18): mu_T^2 >= 2.89/20.09, so mu_T >= 0.3792 (0.37928, rounded DOWN; the first
  version's 0.3793 rounded this lower bound up). The numerical minimum is 0.3939.

(b) Perturbation.
- H_T = gamma5 (1 - 2 (2 + D_W)^{-1}), so H_T1 - H_T0 = 2 gamma5 (2 + D_W1)^{-1} (D_W1 - D_W0) (2 + D_W0)^{-1}, with
  D_W1, D_W0 the Wilson operators D_W(-m) at A1, A0 (not the overlaps D1, D0 of Steps 2 and 3).
- 2 + D_W(-m) = D_W(2 - m), whose hermitian part is >= 2 - m (the Wilson term is >= 0 and the naive part is
  antihermitian), so ||(2 + D_W)^{-1}|| <= 1/(2 - m) for every U.
- Hence ||H_T1 - H_T0||_op <= (4/(2 - m)^2) ||a||_inf and ||H_T1 - H_T0||_HS <= (2 sqrt2/(2 - m)^2) ||a||_2.
- Sup cap eta_T := mu_T (2 - m)^2/8 = 0.0801 (0.080123, rounded down) gives a gap >= mu_T/2 for H_T1.
- Then h_T = (1/2)(2/mu_T)(2 sqrt2/(2 - m)^2) ||a||_2 <= 4.42 ||a||_2. The exact value is 4.4126; the first version
  rounded this upper bound down to 4.41. Compare H_W: sqrt2/m = 2.0203, so <= 2.03 ||a||_2.
- (1.3) changes only through eta: eta_inf,T = 4V exp(-V cos eps0 eta_T^2/((e l)^2 G(L))).

(c) Free counts.
- b_T >= 0 gives |lambda0|^2 >= 1/2.
- b_T < 0 forces b < 0, since 2 + b > 0, so |k_mu| < pi/2. It also forces s^2 < -b(2 + b), so
  n < (2 + b)^2 - b(2 + b) = 2(2 + b) < 4 for every 0 < m <= 1 (Wren, letter 89; the first version had n < 5 from
  the m = 0.7 bound s^2 <= 0.91).
  - Then |lambda0|^2 >= s_T^2/(4 w_T^2) >= s_T^2/4 = s^2/n^2 >= s^2/16 >= |k|^2/(4 pi^2).
  - That uses w_T = |H_T| <= 1, which holds for every k and m, because 1 - m + U >= 0. The first version said "< 1
    since b > -1", which fails at m = 1, k = 0.
  - Numerically the minimum of |lambda0|^2/|k|^2 is 0.67, against the analytic 0.0253: lossy, but explicit.
- General form, if |lambda0|^2 >= c|k|^2 where b < 0: c1^2 = C2/c, n_0(2 delta) <= 16 c1^2 delta^2 V/pi^2 + 4, and
  tau(delta) = 4 C2 + 4 c1^2/pi^2 + (9/(4 pi)) c1^2 ln(sqrt2 pi (1 + 1/L)/(c1 delta)).
  - H_W has c = 1/(pi^2 (2 + m^2)), which recovers (3.2) and (3.3).
  - H_T has c = 1/(4 pi^2), so 16 c1^2/pi^2 = 64 C2. The sealed check below used the weaker c = 4/(25 pi^2)
    (100 C2), which also holds.

So Lemma Z holds for the H_T overlap at every 0 < m <= 1, with mu_T(m) from (a), c_f(H_T) = 0.158939 at m = 0.7, and
the constants above. Lambda keeps the same form, with bigger constants. Audited by Wren (letter 89), with five fixes,
all now in.
- Where the H_W numbers do not carry over (refuter round 2):
  - P_E: the H_T cap eta_T = 0.0801 replaces eta = 0.175 in Section 2. At e l = 0.47, eta_inf,T is 14.7 at L = 16 (no
    bound), 0.50 at L = 24 and 0.004 at L = 32. So section 2's "P_E >= 2/3 - 2e-6 for every L >= 16" becomes
    L >= 32 for H_T.
  - Lambda(L)/V stays above c_f(H_T) until L is about 4500 (about 500 for H_W). The order O(L (ln L)^{3/2}) is
    unchanged; the bound only bites at large L.

Sealed check (`verification/lane1/ht/ht-plan.md`, sha f75d26d5e2187863, sealed at clock reading 12:31:52; run `ht_check.py`, log
`verification/lane1/ht/out-ht-2026-11-26.log`):
- Every killer (Q1 to Q6) held in all six cells: e l in {0.47, 2}, L in {8, 12, 16}.
- The one guess, on the ratio of the paired difference to the H_W value, MISSED. That was a guess, not a killer.

## 9. Refuter record

- Round 1 covered sections 0 to 7 and returned at clock reading 12:41. Verdict: HOLDS, with four minor wounds, none of
  them load-bearing. Record: `verification/lane1/refute1/VERDICT.md`, scripts r1.py to r5.py.
  1. The convex-support case of Caffarelli. Section 2 (1d) now uses only the finite-potential case, plus a weak limit.
  2. The hopping-difference notation (Step 2): the link is now named.
  3. Section 6: "iff" became "if". The ceiling is now called the method's, and the provenance of eps = 0.267 is
     given.
  4. Overclaims: the status line, the monotonicity in L after the 1d table (now shown), and V, which meant both the
     volume and the potential (the potential is now Phi).
- Wren's audit (letter 89, 13:01) covered sections 0 to 8 and ran independent spot checks. Sections 0 to 7 HOLD as
  written. Section 8 got five fixes, all folded in:
  - 8(a) for m < 2 - sqrt2;
  - n < 4 and w_T <= 1;
  - Adams' eps carried over to H_T;
  - h_T 4.41 -> 4.42 and the D_W1/D_W0 notation;
  - c_f(H_T) 0.158943 -> 0.158939.
  Her polish is also in: the a.e. note in section 0; sigma >= |lambda_0|, so the statement loses V ln(1 - m_f^2); the
  cubed penalty in 1d; and the split line before (3.3).
- I checked her numbers before folding them in (sealed `verification/lane1/wren09/check-plan.md`, log `out-check.log`); N1 to N5
  all hit. One catch of my own: "mu_T >= 0.3793" rounded a lower bound up (it is 0.37928), so it now reads 0.3792.
  The sealed H_T check's Q1 used 0.3793 as its threshold, but no sample came near it (minimum 0.3939), so its verdict
  stands.
- Round 2 covered Section 8 (and the steps it inherits) and ran 23:49:56 to about 23:56. Its brief held the claim, its use
  in Route (v) and the question only. Verdict: PARTIALLY REFUTED. The rate survives. One constant was killed: eps = 0.267
  is past eps* = m(2 - m)/(2 + sqrt2) = 0.266533, where Adams' bound reaches zero (-0.00267 at 0.267). The ceilings are
  now open bounds, e l < 0.4727 (H_T) and e l < 0.3367 (H_W). One wound: L_0(H_T) is 32, not 16. Lambda/V > c_f(H_T)
  until L ~ 4500 is recorded. Survived: c_f(H_T), mu_T(m) for m = 0.1 to 1, the inverse-norm bound under an L-BFGS
  ascent (saturated, never exceeded), the H_T1 gap, h_T, the n_0 and T_0 counts, and c_f against the free density at
  U_AP (within O(ln L)). Adams' generalisation is still unread first-hand. Record:
  `verification/lane1/s8-refute/` (RULING.md, refuter/VERDICT.md, my_eps_check.py). The sweep earlier
  today checked the rounding of e l_0 and never the rounding of the eps feeding it.
