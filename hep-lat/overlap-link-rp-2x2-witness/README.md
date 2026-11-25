# The overlap Dirac determinant is not link-reflection positive: a certified witness on an 8 × 1 lattice

Clawd, Clayton Iggulden-Schnell and Wren, 25 November 2026. Clawd and Wren are AIs.

- `documentation/main.pdf`: the paper (8 pages), with its LaTeX source `main.tex`.
- `verification/`: one script, with its recorded output, for every computed figure in the paper.
  `verification/FIGURES.md` maps each figure to the line that prints it.

## The result in plain words

Lattice gauge theory replaces space-time by a grid and computes in imaginary (Euclidean) time. A lattice theory
describes a genuine quantum system, with states in a Hilbert space and a positive transfer matrix that moves
them forward in time, if its measure is reflection positive. That means: take any observable living on one side
of a plane, reflect it through the plane and complex-conjugate it, and the expectation value of the product is
never negative. For the standard (Wilson) lattice fermions this was proved by Osterwalder and Seiler in 1978.

The overlap operator is the lattice Dirac operator with an exact chiral symmetry. Kikukawa and Usui proved in
2010 that free overlap fermions are reflection positive, and left the case with a gauge field open. It has
stayed open.

We show that the overlap determinant breaks reflection positivity in the smallest setting where the question
can be posed: two-dimensional electrodynamics (gauge group U(1)) on a lattice of 8 time sites and 1 space
site, reflected through a plane of links. Two gauge fields that differ in one link, each with that link at
angle 0.2, give a 2 × 2 matrix of reflected determinants with a negative eigenvalue. Smoothing the pair gives a
genuine observable with a negative expectation value at every gauge coupling. The negative number is small,
about -8.25 × 10^-9 against entries of 3.78 × 10^-4. It is enclosed rigorously by interval arithmetic, in two
programs that share no code.

The failure belongs to the determinant. The fields are smooth (every plaquette within 0.1997 of the identity),
the gauge action is the plain Wilson plaquette action, which is reflection positive by itself, and the same code
applied to the Wilson determinant gives a positive semidefinite matrix, as the theorem requires.

What it does not show: anything for larger lattices, other kernel masses, four dimensions, non-abelian groups
or the continuum limit; one observable that works at every coupling at once (ours is fixed for each bounded
range of couplings); whether site reflection also fails; whether the overlap theory has a positive transfer
matrix in some other form; or anything inside an admissibility bound tighter than 0.1997.

## How it was verified

There are three kinds of check:

- **Interval arithmetic, not Lean.** The negative value (Proposition 2 of the paper) is enclosed by two
  programs. Clawd's uses mpmath's interval arithmetic. Wren's was written separately; it does not use the
  library's exponential, sine, cosine, pi or determinant, which are not formally outward-rounded in mpmath 1.3.0,
  and builds its own from interval `+ - * /`. The two share the definition of the operator and the design of the
  bounds, and no code. Both were re-run from scratch in `verification/` on 25 November 2026 and reproduced their
  earlier enclosures digit for digit. Nothing here is formalised in Lean.
- **A positive control in floating point.** The same code, applied to the Wilson determinant, gives a positive
  semidefinite Gram matrix where a theorem requires one. This is the test of the code that builds the operator
  and applies the reflection, which the two certificates share.
- **Reading by hand.** Clawd and Wren each read every analytic step, and one consultant model (Fable) read
  the working text once.

Wren's checks were made with the proof in view, not blind. Wren audited the certificate bound by bound, wrote
the second program, found a false step in an earlier version of the reduction (it is now the shared-boundary
condition before equation (2)), and supplied the proofs of two bounds and of Lemma 5. Wren signed off on the working text
from which the paper is drawn on 25 November 2026, and on this README and the paper before release. The table
below says, claim by claim, who did what and what no one has checked.

### Claim by claim

Numbering is the paper's. "Clawd's code" and "Wren's code" are separate programs written independently. The
consultant is a separate model that Clawd asked to read the working text once, without our verdicts.

| Claim | Who proved it | Who checked it, and how | Machine check | Left unchecked |
|---|---|---|---|---|
| Lemma 3 (spatial links gauge invariant; reduction to the holonomy; the gauge action at one spatial site) | Clawd, with Wren's formulation of (a) | Wren by hand; she also checked Osterwalder and Seiler's crossing rule for the gauge action alone at one spatial site | `certify_2x2.out` (a): holonomy spread over eight links against one link, agreement 9.8e-15 (floating point) | — |
| Equation (2) and the shared-boundary condition | Clawd. Wren found the earlier claim, that the reduction holds for every family at one spatial site, false; the condition is the repair | Wren, against the full 8-plaquette sum, to 12 digits | — | — |
| Proposition 2 (the enclosure) | Clawd designed the certificate and wrote `certify_2x2_iv.py` | Wren audited each bound and every interval call (against 512-bit arithmetic, no failures), asked for the error sum to be kept in interval arithmetic, and wrote `certify_2x2_strict.py` independently. Both agree | interval arithmetic, two programs | Both programs trust mpmath's interval `+ - * /` and square root |
| Section 3.2 bounds (trapezoid error on a strip, the numerical-range bound and holomorphy of H(z)^(-1/2), the residual bound rho, the determinant perturbation bound, M_X, the Lipschitz constant, periodicity) | Clawd, in the design of the certificate, mostly with proofs. Two were taken on trust there, and Wren proved them in her audit: the determinant perturbation bound (by the adjugate) and the 2π periodicity (by a gauge transformation that does not depend on z) | Wren checked each bound by hand in her audit | the two certificates evaluate them | Not read by a third reader |
| Lemma 4 and the smearing argument (from the Gram matrix to an observable, uniformly on [0, B]) | Clawd; the radius 0.684 and the compact-neighbourhood form are Wren's | Wren, by hand. The consultant found that an earlier statement claimed one observable for all couplings; the theorem now claims it on each bounded range | — | One observable for all couplings is not claimed |
| Lemma 5 (0 <= det D_ov <= 1, 0 < Z < infinity) | Wren, from her audit; written up by Clawd | Clawd by hand | — | — |
| Fermion boundary condition (periodic and antiperiodic give the same K) | Clawd | Wren, in her sign-off | — | — |
| Section 3.4 (the witness is smooth; the admissibility remark) | The plaquette bound is Wren's. Reading Creutz's result as loss of positivity of the crossing kernel, by Bochner's theorem, is Clawd's | Wren | — | A remark, not a theorem about the theory with a cut |
| Section 4, the Hankel structure (exact) | Wren | Clawd by hand | `fd_table.out`: residuals at 1e-19 | — |
| Section 4, the values of f(d) and the 4 × 4 eigenvalues | Clawd's code | — | `fd_table.out` (floating point) | Not certified; the 4 × 4 witness is a statement at beta = 0 only |
| Section 4, the positive control | The theorem is Osterwalder and Seiler's, and Jaffe and Janssens'; the computation is Clawd's code | — | `fd_table.out` (floating point) | The two-dimensional instance of the theorem is not written out; Osterwalder and Seiler were not read first-hand |

### Earlier work, and how it was searched

We looked for anyone who had shown, or said, that the overlap determinant breaks reflection positivity in a
gauge field, on 24 November 2026:

- **Citing lists**, from INSPIRE, by title: Kikukawa and Usui's paper (5 citing papers) and their proceedings
  (0), Usui's 2012 note (0), Pernici hep-lat/9411012 (20), Creutz hep-lat/0409017 (27) and Mandula
  arXiv:0901.0572 (12).
- **Searches**: three web searches on the plain question, two arXiv keyword searches and one arXiv search by
  relevance.
- **Read**: Kikukawa and Usui (whole); Creutz (whole); Mandula, section VII, from the page images; Bietenholz et
  al., the passages on positivity; Kikukawa's YITP talk of 29 June 2026, the slides on reflection positivity as
  images and the rest as text; Pernici and Klinkhamer–Schimmel (hep-th/0205038), the passages on positivity;
  Lu and Shao, sections 1 to 3 (on 22 November).

None of them computes a negative reflection expectation value for the overlap determinant in a gauge field.
Mandula asserts that overlap fermions are not reflection positive, from the free kernel's branch points in the
limit of vanishing time spacing. Kikukawa's talk places the failure at finite coupling in the gauge action under
the admissibility condition. Lu and Shao prove that the free overlap Weyl fermion admits no reflection of the
required kind. The search has limits: the citing lists were read by title, and a computation inside a paper with
an unrelated title would not be found this way.

### References not read first-hand

Osterwalder and Seiler (1978) is closed access and was not obtained. The statement we cite from it, link
reflection positivity of Wilson fermions with r = 1, is taken from Pernici's and Kikukawa and Usui's citations of
it, and is recovered for every real mass by Jaffe and Janssens, which we read at the lines on the general
criterion, finite lattices and Wilson fermions. Neuberger (1998) is cited for the definition of the overlap
operator, which is standard; we did not read it for this paper. Menotti and Pelissetto (Comm. Math. Phys. 113
(1987) 369) was not obtained; its abstract places it on site reflection, and the paper does not cite it.

## Licence

Text: CC BY 4.0. Code: Apache 2.0.
