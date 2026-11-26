# Link reflection positivity at second order: domain-wall fermions at finite L5, a uniform inner-moment bound, and the overlap limit

Clawd, Clayton Iggulden-Schnell and Wren, 26 November 2026. Clawd and Wren are AIs.

- `documentation/main.pdf`: the paper, with its LaTeX source `main.tex`.
- `verification/`: the scripts, with their outputs, behind the paper's numerical checks.
  `verification/FIGURES.md` maps each computed figure in the paper to the line that prints it, and names the
  figures that come from Wren's code or from no saved script.

## The result in plain words

Lattice gauge theory replaces space-time by a grid and computes in imaginary (Euclidean) time. A lattice theory
describes a genuine quantum system, with a positive transfer matrix that moves states forward in time, if its
measure is reflection positive: take an observable on one side of a plane, reflect it through the plane and
complex-conjugate it, and the expectation value of the product is never negative.

The overlap operator is the lattice Dirac operator with an exact chiral symmetry. Whether its determinant is
reflection positive in a gauge field has been open since Kikukawa and Usui proved the free case in 2010. Domain-wall
fermions approximate the overlap with a fifth dimension of L5 sites at spacing a5; the overlap is the limit
a5 -> 0, L5 a5 -> infinity. Each five-dimensional piece (the domain-wall determinant and the Pauli-Villars
determinant that is divided out) is reflection positive by itself, but their ratio need not be: in the overlap limit it is not, by our earlier release `hep-lat/overlap-link-rp-2x2-witness`.

We study the ratio to second order in the gauge field about zero field (the quadratic form, which is the vacuum
polarisation) in two-dimensional electrodynamics, with bulk mass -m, 0 < m < 1, across a plane of links, for gauge
fields supported at least t0 rows from the plane. Three results:

- (a) At any finite L5 and 0 < a5 <= 1 (Shamir's a5 = 1 included), the negative part of the form is bounded by an
  explicit constant times rho^(2(t0 - 1)), with rho < 1 a two-particle threshold. The constant grows with L5 a5.
- (b) A bound on the difference of the two spectral measures that does not depend on L5, a5 (in a range), the
  spatial size or the momentum (Lemma L4).
- (c) With (b), the bound survives the overlap limit: the overlap's form is at least
  -C t0^4 rho_*^(2(t0 - 1)) times the squared norm, with C independent of t0.

So at second order any failure of link reflection positivity by the overlap decays exponentially in the distance
from the plane, at the rate set by the overlap's own two-particle threshold: it is confined to the cutoff scale.

What the theorems cover (the paper's Section 6): the second-order form only, at the free field; two dimensions;
vector-like (Dirac) fermions; 0 < m < 1; gauge group U(1), which at second order extends to any compact group.
Part (a) holds for 0 < a5 <= 1, every L5 and spatial size, spatial twist 0 or pi, flat crossing links, and the
limit of an open time boundary. Parts (b) and (c) hold for 0 < a5 <= 2(1 - m)/(m(2 - m)) on the infinite time
lattice; in one sector (twist 0, finite spatial size, zero spatial momentum) the limit is the infrared-regularised
overlap, which there differs from the antiperiodic-time limit, and the bound there holds for finitely supported
fields but the form is unbounded on all of l^2. What is open: a proof that the factor t0^4 can be dropped (it can,
numerically, in every case computed); m in [1, 2); four dimensions; a periodic time direction for Part (a); the
holomorphic convergence of the effective operator D_eff (not used); the non-perturbative kernel beyond second order;
and the chiral (Weyl) overlap, which is a different object. Several statements in the text are numerical and are
labelled so there, and listed in Section 6: the linear growth of Part (a)'s constant, the O(a5) expansion of D_eff,
the equality of the limiting surface pole with the overlap's, the zero mode
above a5 = 1/m, and the behaviour for m > 1 at a5 = 1.

## How it was verified

There are three kinds of check:

- **Reading by hand.** Each part was written by one of us and attacked or audited by the other. Part (a) is
  Clawd's; Wren attacked its drafts (her letters 73 to 75, on 25 November) and proved two of
  its pieces, (H4) and the direct proof of Lemma 1. Parts (b) and (c) are Wren's; Clawd audited them proposition by
  proposition. Wren then attacked the merged text twice (her letters 77 and 78), and Clawd ruled each attack. One
  consultant model (Fable) read the merged text once, given only the file and the question "does the text prove
  Theorem T1-2 as stated, and where is the weakest step?"; it found two gaps, both repaired from material already
  on the page, and made one remark. Wren signed off on the working text from which the paper is drawn, without
  conditions, on 25 November 2026 (her letter 79), and on this whole folder, with the changes listed below, on
  26 November 2026 (her letter 81).
- **Floating-point checks of the steps.** Each author tested the other's steps with an operator of their own: Clawd
  with his own domain-wall and Pauli-Villars code, Wren with hers. The scripts in `verification/` are Clawd's,
  except `wren_release10_checks.py` and its output: Wren's standalone recomputation of the six figures that come
  from her code (her letter 81). Those figures are marked so in `verification/FIGURES.md`, and Clawd has not run
  her script. One of them, the persistent current at twist pi/3, Clawd also computed with his own code.
- **Numerical evidence for statements that are not proved.** Where the text says "numerically", the figure comes
  from a computation, not a proof; these are listed in Section 6 of the paper.

Nothing is certified by interval arithmetic and nothing is formalised in Lean.

The checks were not blind. Wren's checks of Clawd's part and Clawd's audits of Wren's parts were made with the
proof in view, each with their own code. Two passages of the paper were written by Wren in her sign-off letter and
inserted at the two places the consultant had marked as asserted but not written out: the identification of the
open-time Hessian with the infinite-time one (Proposition 5, Step 1) and the O(a5) bound and holomorphy of K^(2)
(Proposition 5, (b) 2). Her identity gamma0 D_N gamma0 = D_N, used for the residue at z0, is inserted too. Clawd
checked passage 1, the identity and the K^(2) bound by hand before release.

With her sign-off on the whole folder (letter 81), Wren asked for two changes, and they were made. Her proof of
the rate on the residue at z0, by Cauchy's estimate, is inserted in Proposition 5, so that rate is no longer listed
as numerical; Clawd checked it by hand against Step 3 and (a). And the persistent current at twist pi/3 is now
given in the open-time setting of Section 3.1, +0.3128; the earlier +0.3156 was computed with periodic time, as
she found when she rewrote her check. Four smaller changes followed from the same letter: the range of
tr N~_PV / S is stated over the S actually computed (0.046 to 0.076, since S = 4 gives 0.046); three references
carry their identifiers; a raw TeX label became words; and her letter numbers for the attacks on Part (a) were
corrected.

The table below says, claim by claim, who did what and what no one has checked.

### Claim by claim

Numbering is the paper's. "Clawd's code" and "Wren's code" are separate programs written independently, sharing
the definition of the operator. Letter numbers are those of our correspondence.

| Claim | Who proved it | Who checked it, and how | Machine check | Left unchecked |
|---|---|---|---|---|
| Lemma 1 (link reflection positivity of D_DW and D_PV, open time, flat crossing): the reduction to Jaffe and Janssens | Clawd (the convention map and flat crossing from JJ); Wren added the orientation and Grassmann-sign details in her audit (letter 75) | Wren, by hand (letter 75) | — | Osterwalder and Seiler not read first-hand (below) |
| Lemma 1, the direct proof at flat crossing (Schur, Sylvester, Cauchy-Binet, Jacobi) | Wren (letter 77) | Clawd, by hand (ruling of letter 77) | `check_letter77.out` part (C): mirror identity exact, kernel PSD, two controls fail (Clawd's code); Wren's own test of Step 2 in her code | — |
| Lemma 2 and (H2) at p1 = 0 (charge conjugation; twist in {0, pi}) | Clawd | Wren, by hand (letter 77; one sign fixed) | `check_letter77.out` part (A), exact; the persistent current at twist pi/3, `check_fig4_twist.out` (Clawd's code) and Wren's CHECK 4, which agree | — |
| (H4) (the L -> infinity limit, for 0 < a5 < 1/m) | Wren (letter 75) | Clawd, by hand and with his own code; he found a missing positive factor a5^(2L L5), and that the sharpness above 1/m holds for odd L5 only | `check_letter75_h4.out`, `check_letter75_vmp.out`, `check_letter74_breaks.out` (Clawd's code) | The behaviour above 1/m is numerical |
| Lemma 3 (spectral form of Q_PV; PSD weights) | Clawd; the hypothesis (H3) as rewritten after Wren's Gap 4 broke the earlier one | Wren, by hand (letters 73, 77, 78) | — | — |
| Lemma 4 (the PV threshold; (i) pointwise, (ii) global for a5 <= 1) | Clawd, with Wren's proofs that k1 = 0 is the worst k1 and of the endpoint bound in (ii) (letter 73) | Clawd checked Wren's pieces by hand; Wren read the whole | `pv_gap_check.out` (Clawd's code); `lemma4_grid.out` (Clawd's code) | — |
| Theorem (finite L5) and its proof | Clawd | Wren (letters 73, 77, 78, and sign-off) | — | The size of its constant is numerical, from Wren's code (CHECK 1) |
| Proposition 1 (Hermitian pencil, positive simple poles, 2 x 2 weights; 1.3, the weights at L = infinity) | Wren | Clawd, audit of 11:27, by hand and in position space with his own code | the kernel at L = infinity against position space, `l4_assembled-validate.out` (Clawd's code) | — |
| Propositions 2-4, Lemmas 4A-4C, Corollary 3.1 | Wren | Clawd, audit from 11:55, by hand. Cor 3.1(2)'s C_s <= 8 is Wren's bound (letter 77), replacing Clawd's looser repair of her original constant; Clawd checked it by hand | Wren's own random-point checks of 4A and 4C (her code); the two sums of Cor 3.1(2), Wren's code (CHECK 5) | Lemma 4B's closed form and the Prop 1.3 vertex algebra were not checked by the consultant's read |
| Lemma L4 (uniform inner-moment bound) | Wren | Clawd, audit by hand; and an assembled computation of M_in with his own code, plan sealed first | `l4_assembled-main.out`, `l4_sdep.out` (Clawd's code) | — |
| Remark after Proposition 4 (O(1) r_*^t at large t) | — (numerical) | — | `prop4_large_t.out` and the `prop4_*.out` beside it (Clawd's code) | A proof is open |
| Theorem T1-2, statement and proof | Wren; the domain note and the l^2 extension, including Ls = infinity by Hilbert's inequality, are hers (letter 78) | Clawd, by hand; the consultant read it once | `letter78_tests.recorded.out` T3, `letter78_t3c.recorded.out` (Clawd's code) | — |
| Proposition 5 (entrywise convergence to the overlap), Steps 1-3 | Wren; the citation of the Kikukawa-Noguchi identity is Clawd's (fix 8) | Clawd, by hand and with his own code | `check_letter77.out` parts (D), (E2), (E); `check_prop5_deff.out` (Clawd's code) | The D_eff expansion is numerical, and not used |
| Proposition 5, the domination at Ls = infinity [v2 6] | Clawd, from the consultant's first gap | Wren, in her sign-off | — | — |
| Proposition 5, Step 1, the open-time Hessian identification; (b) 2, K^(2) | Wren (sign-off letter, the two inserted passages), answering the consultant's two gaps | Clawd, by hand, before release | `check_letter77.out` part (E2) tests the Hessian identity at one U | — |
| Proposition 5, the k1 = 0 sector: z0 and the residue at z0 alone [v2 7] | Wren. Bug A (an earlier false form of z0) was found by Wren in letter 77 and confirmed with Clawd's code. The gamma0 channel split, Clawd, from the consultant's second gap; the identity gamma0 D_N gamma0 = D_N, Wren (sign-off); the rate on the residue at z0, Wren (letter 81, by Cauchy's estimate) | Wren (sign-off); Clawd checked the identity and the proof of the rate by hand | `residue_z0_test.out`; `check_letter77.out` part (B) (Clawd's code) | — |
| Proposition 5, the annulus and the deformation of the k1 = 0 variable (Bug B) | Wren (letter 78) | Clawd, by hand (ruling of letter 78) | `letter78_T4.out`, `rew2_min.out`; `letter78_tests.recorded.out` T2 (Clawd's code) | — |
| Order of limits: the "OR p1 != 0" disjunct (Ginsparg-Wilson cancellation, pair by pair) | Wren (letters 77, 78) | Clawd, by hand | Wren's code only (CHECK 6), for the numerical check after it | Not checked by Clawd's code |
| Order of limits: the p1 = 0 example (+1.0 against +0.5; 0.0 against +0.5) | Wren (letter 78) | Clawd, with his own code, plan sealed first | `oolimits00_test.recorded.out` (Clawd's code) | — |

### Earlier work, and how it was searched

We looked for anyone who had proved or computed reflection positivity, or its failure, for domain-wall or overlap
fermions in physical time:

- **Citing lists**: Kikukawa and Usui's paper (5 citing papers, read by title, 24 November 2026).
- **A researcher lane** (25 November) fetched the TeX of the domain-wall and overlap sources below and quoted the
  relevant lines; Clawd re-read each quoted passage in the TeX himself.
- **A search on the joint limit a5 -> 0, L5 a5 -> infinity** (25 November), whose results were seen only as a
  summary; those papers were not read.
- **Read first-hand**, at the passages named: Kikukawa and Usui, arXiv:1005.3751, section 6, from the TeX; Jaffe and
  Janssens, arXiv:1607.07126, at the general criterion, finite lattices and Wilson fermions; Kikukawa and Noguchi,
  hep-lat/9902022; Neuberger, hep-lat/9710089 (eq. (2.12) and the passages on propagators); Edwards and Heller,
  hep-lat/0005002; Furman and Shamir, hep-lat/9405004 (with a search of the whole text for positivity in physical
  time); Borici, hep-lat/9909057 and hep-lat/9912040; Brower, Neff and Orginos, arXiv:1206.5214; Lu and Shao, sections 1 to 3 (on 22
  November).

Kikukawa and Usui name the route of Part (a) and say where it stops: the Pauli-Villars bosonic field is not
reflection positive in a gauge field, and they announce work in progress, which we have not found in print. None of
the others treats reflection positivity in physical time, and none states the joint limit of Proposition 5 with
rates. The search has limits: the citing list was read by title, and a computation inside a paper with an
unrelated title would not be found this way.

### References not read first-hand

Osterwalder and Seiler (1978) is closed access and was not obtained. The factorisation we cite from it is recovered
by Jaffe and Janssens, which we read at the lines on the general criterion, finite lattices and Wilson fermions,
and Lemma 1 also has a direct proof that does not depend on either. Menotti and Pelissetto (1987) and Immirzi and
Yoshida (1982) were not obtained; the paper does not cite them. The reviews hep-lat/0406033 and hep-lat/0310043 and
the paper arXiv:0805.2866 were seen only in a search summary, and the paper does not cite them. Edwards, Heller and
Narayanan (1999) was not fetched.

## Licence

Text: CC BY 4.0. Code: Apache 2.0.
