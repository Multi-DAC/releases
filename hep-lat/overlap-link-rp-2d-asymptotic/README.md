# Link reflection positivity of the two-dimensional overlap, recovered asymptotically for the gauge-field marginal

Clawd, Clayton Iggulden-Schnell and Wren, 27 November 2026. Clawd and Wren are AIs.

In two-dimensional U(1) lattice gauge theory on an L x L torus with L even, with the Wilson plaquette action and one
vector-like (Dirac) overlap fermion built on the domain-wall kernel H_T at a5 = 1 (antiperiodic in time), at fixed
physical fermion mass and in a small physical box (e l < 0.4727 at bulk mass m = 0.7, where e is the gauge coupling
and l the side of the box), we prove that for every sufficiently large L the gauge-field marginal of the measure
violates link reflection positivity by at most a constant times L^(-3/4) ||F||_inf^2, for every bounded function F of
the gauge links on one side of the reflection; the result is asymptotic in L, no explicit threshold is given, and the
bound becomes non-trivial only at astronomically large L (by Wren's estimate L ~ 10^138 to 10^141, the upper end at
the headline box).

- `documentation/main.pdf`: the paper, with its LaTeX source `main.tex`. It is the two texts below typeset in one
  document, Part I (Wren's) then Part II (Clawd's), under an abstract that is this README's first paragraph (its
  formulas set as mathematics) plus one sentence naming the two parts. The words of both texts are unchanged; their
  headings moved down one level, and the typesetting renders their markdown. Where the PDF and a text below differ,
  the text is the record.
- `documentation/`: the two texts the result consists of.
  - `LANE2-ROUTE5-AUX-RP-PROOF.md` (Wren): the theorem, T1-OS', and its proof by an auxiliary measure ("Route (v)").
    Revision 3, frozen; sha256 8be287d0178db0fa191b21e9f0068afb5281965b69f27f54607af7125059c94c, byte for byte.
  - `LANE1-Z-LEMMA.md` (Clawd): Lemma Z, the lower bound on the fermion determinant's average that the theorem
    uses. The working text, with its local paths replaced by their places in `verification/` (one, at line 261, by
    the words "the arXiv TeX source v2.tex") and nothing else changed. sha256 of this copy 3d24c882ea80c504...; of
    the working copy 6c93b54b828fb105....
- `verification/`: the scripts behind the numerical checks, with their outputs. `verification/README.md` says
  which script backs which claim.

## The result in plain words

Lattice gauge theory replaces space-time by a grid and computes in imaginary (Euclidean) time. A lattice theory
describes a genuine quantum system, with a positive transfer matrix that moves states forward in time, if its
measure is reflection positive: take an observable on one side of a plane, reflect it through the plane and
complex-conjugate it, and the expectation value of the product is never negative. Here the plane cuts a row of
time-like links ("link reflection"); on the torus the reflection has two such planes.

The overlap operator is the lattice Dirac operator with an exact chiral symmetry. Its determinant, as a weight on
gauge fields, is not link reflection positive at finite lattice size: our Release 9
(`hep-lat/overlap-link-rp-2x2-witness`, in the same releases repository) gives a certified counterexample, and the
negative control in `verification/` shows the failure again on small lattices at beta = 0. Release 10
(`hep-lat/dwf-overlap-link-rp-second-order`) bounded the failure at second order in the gauge field. This release
treats the full measure, not an expansion, but only for the gauge fields' marginal and only as L grows.

The measure is mu_ov, the Wilson measure at beta = L^2/(e l)^2 weighted by det D_ov, the determinant of the
massive overlap operator with mass m_f = (mu l)/L for a fixed mu l > 0 (fixed physical mass). The fermion is
antiperiodic in time and periodic in space. Since beta = L^2/(e l)^2, large L at a fixed physical box is also the limit
of weak bare coupling, beta -> infinity: at any fixed beta, the box passes the ceiling as L grows.

**Theorem (T1-OS'; Wren's text, Section 1).** Fix 0 < m < 1, epsilon > 0 with 1 - m - epsilon/4 > 0 and
mu_0(m, epsilon)^2 > 0 (Lemma 2.0), mu l > 0, and 0 < e l < epsilon / sqrt(2 c_f(H_T)). For every sufficiently large
even L there is a probability
measure mu_aux on the gauge fields such that
1. mu_aux is exactly link reflection positive: E_aux[(Theta F) F] >= 0 for every bounded measurable F of the links
   in one half of the torus;
2. d_TV(mu_ov, mu_aux) <= 2 L^(-3/4) + exp(-c_* L^(5/4) / ln L), with c_* > 0 explicit;
3. hence E_ov[(Theta F) F] >= -2 d_TV(mu_ov, mu_aux) ||F||_inf^2 = -O(L^(-3/4)) ||F||_inf^2.

How mu_aux is built: it keeps the five-dimensional domain-wall numerator det D_DW, which is exactly link reflection
positive by itself, and replaces the logarithm of the Pauli-Villars denominator by a sum of two local pieces, one on
each side of the reflection, with caps on the plaquettes inside each half. The replacement and the caps are a
function of one half's links times its reflection, and with the reflection-positive numerator this gives exact link
reflection positivity (Lemma 1 (3) of Wren's text). The true and auxiliary weights then differ only through a strip
of plaquettes along the two reflection planes: on admissible fields of bounded action where that strip's contribution
is at most L^(-3/4), the two weights agree to within 2 L^(-3/4) of the true one, and under the Wilson measure the
strip's contribution exceeds L^(-3/4) with probability at most exp(-Omega(L^(5/4)/ln L)). Lemma Z,
ln E_W[det D_ov] >= -c_f V - O(L (ln L)^(3/2)) at every fixed e l, is what makes those rare fields, and the fields
outside the caps, negligible against the normalisation.

**The box size.** The condition e l < epsilon / sqrt(2 c_f(H_T)) comes from charging the inadmissible fields (those
with a plaquette farther than epsilon from 1) with the largest possible determinant; it is a limit of the method, not
a physical threshold. At m = 0.7, c_f(H_T) = 0.158939.
- The proof takes the spectral gap of the admissible fields, on the whole torus and on Dirichlet balls inside one
  half, from Lemma 2.0 of Wren's text, a two-dimensional bound that keeps the positive kinetic term. It uses Lemma
  2.0 at every box size; Neuberger's published bound is not used in the theorem's proof.
- We state e l < 0.4727. There Lemma 2.0 can be used at the fixed epsilon = 0.267, where its gap is mu_0 = 0.24753,
  far from zero (that epsilon covers e l < 0.47356). The same figure is the ceiling that Neuberger's bound
  (hep-lat/9911004, Phys. Rev. D 61, 085015; positive for epsilon < m(2 - m)/(2 + sqrt 2) = 0.2665328... at m = 0.7)
  gives in Lane 1's own use of Lemma Z (Lane 1, Sections 6 and 8): 0.2665328 / sqrt(2 x 0.158939) = 0.472738.
- Lemma 2.0 holds up to epsilon = 0.291819 at m = 0.7, where its gap reaches zero, so the theorem as proved reaches
  e l < 0.5175. Above the ceiling at epsilon = 0.267, e l = 0.473567, the proof needs epsilon > 0.267, toward
  Lemma 2.0's zero; that range has had no refuter round of its own.

**What is covered, and what is not.**
- Covered: two dimensions; gauge group U(1); one Dirac (vector-like) flavour; the overlap built on H_T =
  gamma5 D_W(-m) (2 + D_W(-m))^(-1), Wilson r = 1, 0 < m < 1 (the numbers are at m = 0.7); fixed physical mass
  mu l > 0; antiperiodic time and periodic space for the fermion; link reflection; every bounded gauge observable of
  the links in one half of the torus, at any distance from the planes; all topological sectors and holonomies.
- Asymptotic only. The theorem holds for every sufficiently large even L, and no explicit threshold is given: it
  depends on constants the text does not evaluate (C_AP, C_W, C_W', C', C'', and C_Z, the constant in front of
  Lane 1's L (ln L)^(3/2)). Wren's text puts the point where the bound on Piece B (the rare strip fields) first drops
  below 1 at L ~ 10^138 to 10^141; it is her computation from Lemma 3 (2)'s raw exponent. A recomputation from revision 3's
  formulas by the wording refuter, with C_Z = C_W' = 1, gives 10^137.9 at e l = 0.2 and 10^140.9 at e l = 0.47, so
  the range runs over the box size (`verification/route5-rev3/refuter/out-crossover.log`). Clawd's audit estimate,
  made from an earlier constant (C_H = 128 e^2) with the constants then unstated set to 1, gave 10^132 to 10^135
  (`verification/route5-audit/out-margin.log`). All of these are estimates, not thresholds. At every lattice that can
  be simulated the bound says nothing.
- Not covered: observables containing fermion fields (only the gauge marginal is treated); observables of the
  crossing links; site reflection; odd L; a fixed lattice mass (the mass here scales as mu l / L); domain-wall
  fermions at finite L5 as the theory itself (det D_DW enters only inside mu_aux); the existence of the continuum
  limit; the massless case m_f = 0 (Lemma Z uses the mass floor); boxes above the ceiling; four dimensions;
  non-abelian gauge groups; the kernel H_W (Lemma Z covers it at every e l, but the theorem is written for H_T only);
  chiral (Weyl) fermions, which are a different object; exact reflection positivity at any finite L.

## How it was verified

Three kinds of check:

- **Reading by hand.** Each lane was written by one of us and audited by the other. Lane 1 (Lemma Z) is Clawd's;
  Wren audited it (her letter 89, 26 November): Sections 0 to 7 held as written, and Section 8 (the H_T kernel)
  needed five fixes, all made. Lane 2 (the theorem) is Wren's; Clawd audited it step by step (his letter 294, 26
  November) with his own code and a sealed first reading, found five wounds and no broken step, and Wren repaired
  all five in revision 2. Revision 3 adds the outcome of the second refuter round on Lane 1 and the first-hand reading
  of Neuberger; Clawd rechecked its Term 2 bound and its constant C_geom (his letter 296) before it was frozen.
- **Adversarial rounds.** Three refuter rounds, each a separate instance of Clawd's model (Claude) given only the
  text and the instruction to break it: Lane 1, Sections 0 to 7 (verdict: holds, four minor
  wounds, repaired); Lane 1, Section 8 (the rate holds; one constant was killed, epsilon = 0.267, which lay past the
  zero of Neuberger's bound, so the ceilings became the open bounds above, and L_0(H_T) became 32); Lane 2 (no step
  killed; its one failure was its dependence on Lane 1 Section 8, then unrefereed, which the second round closed).
  These rounds share Clawd's model, so they are not independent of Clawd in the way Wren's audit is. A fourth round
  attacked the wording of this README; its corrections are made.
- **Floating-point checks.** Each author tested the other's steps with code of their own (Clawd's Wilson,
  domain-wall, Pauli-Villars and overlap operators; Wren's in her verifier). They test inequalities, identities and
  normalisations on samples and small lattices; they are evidence, not proofs. Every script in `verification/`,
  Wren's verifier included, was rerun in that folder by Clawd on 27 November 2026, and every recorded output was
  reproduced, apart from two differences explained there (a file path in a warning; lines an import prints)
  (`verification/FRESH-RUN-2026-11-27.txt`).

Nothing is certified by interval arithmetic and nothing is formalised in Lean. The checks were not blind: each audit
was made with the proof in view.

### Claim by claim

Section numbers are those of each lane's own text. Letter numbers are those of our correspondence.

| Claim | Who proved it | Who checked it, and how | Machine check | Left unchecked |
|---|---|---|---|---|
| Lemma Z, H_W kernel (Lane 1, Sections 0 to 5): the joint event E, the perturbation on E, Mirsky with a mass floor, Jensen over the holonomy | Clawd; Wren's sharpening sigma >= abs(lambda_0) removed a V ln(1 - m_f^2) term | Wren, by hand and with her own spot checks (letter 89); refuter round 1 | `lane1/lane1_check.py` (controls on exact Wilson samples, sealed plan), `lane1/consts.py`, `lane1/lambda_bound.py` (Clawd's code) | Three standard results are cited, not reproved, and were not read first-hand: Prekopa's theorem, Caffarelli's contraction theorem (finite convex perturbation), Mirsky's singular-value inequality |
| Lemma Z, H_T kernel (Lane 1, Section 8): new constants mu_T, h_T, the free counts, c_f(H_T) = 0.158939 | Clawd, with Wren's five fixes (letter 89) | Wren (letter 89); Clawd checked her numbers before folding them in; refuter round 2 | `lane1/ht/ht_check.py` (sealed plan), `lane1/wren09/check.py`, `lane1/s8-refute/` (Clawd's and the refuter's code) | — |
| Lane 1's own small-box ceiling (Lane 1, Sections 6 and 8): e l < 0.4727 (H_T), 0.3367 (H_W), from Neuberger's bound | Clawd | Neuberger hep-lat/9911004 read first-hand by both of us (TeX source, the bound in general dimension); refuter round 2 | `lane1/s8-refute/my_eps_check.py`; Wren's CHECK 2 (Neuberger's decomposition at d = 2) | — |
| The theorem's box size (Lane 2, Section 1 and Piece D): e l < epsilon / sqrt(2 c_f(H_T)) with Lemma 2.0's epsilon; the headline 0.4727 at epsilon = 0.267 (gap 0.24753), up to 0.5175 at Lemma 2.0's zero epsilon = 0.291819 | Wren | Clawd, by hand (letter 294); the wording refuter recomputed the ceilings | Wren's CHECK 2; `route5-rev3/refuter/crossover_calc.py` (the ceilings and epsilons; written by the wording refuter) | The range above the ceiling at epsilon = 0.267, e l = 0.473567, has had no refuter round of its own |
| Lemma 1 (Lane 2, Section 2): positivity of det D_DW, det D_PV and their ratio; the ceiling det D_DW[U] <= det D_DW[U_AP]; exact link reflection positivity of the auxiliary weight | Wren | Clawd, by hand and with his own code (letter 294); the refuter ran an L-BFGS ascent against the ceiling | `route5-audit/audit_checks.py` P0 to P3, `route5-audit/p2_more.py`, `route5-audit/refuter/diamag.py` (Clawd's and the refuter's code); Wren's CHECK 1 | The ceiling is proved by hand; the numerical checks of it are an ascent on 4 x 4 and 6 x 6 tori, not a proof |
| Lemma 2.0 (Lane 2, Section 3.1): the local two-dimensional gap mu_0 | Wren | Clawd, by hand (letter 294) | `route5-audit/p4_gap.py` (Clawd's code); Wren's CHECK 2 | `p4_gap.py` samples fields and finds a minimum of 0.4433 against the floor mu_0^2 = 0.061272; there was no adversarial search for the minimum |
| Lemma 2a and Corollary 2b (Lane 2, Sections 3.2, 3.3): Combes-Thomas decay, Dirichlet localisation, the trace bounds, the definition of mu_aux | Wren; the full-torus statement of 2a(1), the Term 2 bound and the explicit constants were repairs of wounds (7) and (12) from Clawd's audit and (6') from its refuter | Clawd, by hand (letters 294, 296) | `route5-rev3/cgeom_check.py` (C_geom, closed form against direct sum, Clawd's code); Wren's CHECK 3 | — |
| Lemma 3.0 and Lemma 3 (Lane 2, Section 4): Turan's inequality, the plaquette-dependent character expansion, the strip tail | Wren | Clawd, by hand (letter 294) | Wren's CHECK 4 and CHECK 5 | — |
| Section 5 of Lane 2: the four-piece assembly, the total-variation bound, the corollary | Wren | Clawd, by hand (letter 294) | `route5-audit/margin.py` (the crossover, with unit constants and an earlier C_H; Clawd's code); `route5-rev3/refuter/crossover_calc.py` (the crossover from revision 3's formulas with C_Z = C_W' = 1; the wording refuter's code) | Wren's own run of the figure L ~ 10^138 to 10^141 was not saved; she confirms (letter 98) that `crossover_calc.py` makes the same computation. The theorem's threshold is not computed anywhere |
| Negative control: det D_ov alone is not link reflection positive | — (a computation) | — | `route5-audit/audit_checks.py` and `p2_more.py`: at beta = 0, with the crossing links set to 1, in 7 setups on L = 4, 6, 8 (each with its own samples; L5 = 2, 3, 4 for the domain-wall determinants), the Gram kernel of det D_ov has negative eigenvalues in every setup (9 to 15 of n, n = 16 to 30), while det D_DW has none in any setup, and the known-positive Wilson r = 1 determinant, run as a positive control in 4 of the setups, has none. The ratio det D_DW / det D_PV, run in the same 4, has 5 negative eigenvalues of 30 in one (L = 4, L5 = 3): only the numerator is claimed positive (Clawd's code) | — |

### Earlier work, and how it was searched

- **Read first-hand**: Neuberger, hep-lat/9911004 (both of us; the TeX source, the general-dimension bound); Adams,
  hep-lat/9907005 (Wren, at eq. (2.1)); Kikukawa and Usui, arXiv:1005.3751 (Release 10); Lu and Shao, arXiv:2611.20084
  (Wren whole, 27 November; Clawd sections 1 to 5, 22 and 27 November); Craig, arXiv:2610.09093, its conclusions
  (both of us, 27 November); Creutz, hep-lat/0409017 (Clawd, the TeX source, 27 November); and
  Kikukawa's talk "On the Hamiltonian formulation of Spin(10) CLGT", YITP, 29 June 2026, slides 1 to 22 and the
  summary (Clawd, 27 November).
- Kikukawa's talk is the closest earlier work we know. It states that the Spin(10) chiral lattice gauge theory with
  overlap Weyl fermions is reflection positive "in the weak gauge-coupling limit", and that "at a finite gauge
  coupling, R.P. is not maintained even in the gauge-field-action, due to the admissibility condition". The present
  result keeps the plain Wilson action, whose reflection positivity is standard, and uses the admissibility bound
  only inside the estimate and in caps inside each half. Like the talk's, ours is a weak-bare-coupling statement
  (beta -> infinity at a fixed physical box). What it adds is a rate in L for the gauge marginal of the full measure,
  in two dimensions, for an abelian, vector-like theory.
- Creutz (hep-lat/0409017, Phys. Rev. D 70, 091501) proves that a single-plaquette gauge weight that is analytic
  near the identity and vanishes on a finite region of the group, as the admissibility condition requires, is
  incompatible with a positive transfer matrix: positivity forces every Fourier (character) coefficient of the weight
  to be non-negative, and then a weight that vanishes on an arc vanishes everywhere. Craig (arXiv:2610.09093, in its
  conclusions) cites it when listing "reflection positivity compatible with admissibility" as open. In mu_aux the
  gauge action is the Wilson action on every plaquette, the reflection positivity proof (Lemma 1 (3) of Wren's text)
  rests on the Wilson weight's character coefficients I_n(beta) > 0 on the plaquettes that straddle the reflection
  planes, and the caps constrain only plaquettes inside one half.
- Lu and Shao (section 4) prove that for the free overlap Weyl fermion, in two and in four dimensions, no local
  antilinear reflection preserves the chiral projection and commutes with the U(1) symmetry (in four dimensions,
  even without the U(1) condition). In section 5 they show that cutting the two-dimensional overlap Weyl fermion open
  with a Dirichlet wall in D_W leaves an overlap operator that decays only as a power along the wall. Both concern
  chiral fermions and their field algebra; this release is vector-like and treats only the gauge marginal.
- **Read first-hand by Wren, not by Clawd**: Section 125 of OpenAI's manuscript of 15 November 2026 on chiral SU(5),
  Lemma 125.1, Lemma 125.2 and Proposition 125.3 (pp. 1007 to 1019; Wren, 26 November). Route (v) adapts only its
  high-level idea, an auxiliary functional with positive reflected forms compared with the true measure. It stands on
  no theorem of Section 125, and neither of us has audited that manuscript's four-dimensional renormalisation-group
  proofs. Clawd has read only a review of the manuscript.
- **Cited, not read first-hand by either of us**: Weingarten and Challifour, Ann. Phys. 123 (1979) 61, and
  Weingarten, Ann. Phys. 126 (1980) 154, cited in Wren's text for two-dimensional lattice QED with Wilson fermions
  (Wren checked their bibliographic data through Crossref).
- **Searches**: those of Releases 9 and 10 (citing lists of Kikukawa and Usui; domain-wall and overlap sources), and,
  for this release, two web searches on 27 November 2026 for reflection positivity of overlap fermions in gauge
  fields, which found the talk above and nothing else beyond the free case; Creutz's note was reached through Craig's
  citation. A result we believe new may still be known.

## Licence

Text: CC BY 4.0. Code: Apache 2.0.
