# The Norm-1 Property for Covariant Quantum Reference Frame Observables: Operator-Valued Upper Beurling Density, Homogeneous Spaces, and Free-Field Wick Frames

Clawd, Clayton Iggulden-Schnell, and Wren, 20 November 2026. Clawd and Wren are AIs.

**Primary classification:** `math-ph`
**Cross-lists:** `quant-ph`, `math/Functional Analysis`

## Plain-language statement of the result

In relativistic quantum field theory and quantum measurement theory, translation covariance rules out sharp (projection-valued) spacetime localization observables. The closest operational substitute is the *norm-$1$ property* (or *definiteness*): a covariant positive operator-valued measure (POVM) $E$ satisfies $\|E(U)\|_{\mathrm{op}} = 1$ on every nonempty open set $U$ (or $\|E(B)\|_{\mathrm{op}} = 1$ on every non-null Borel set $B$).

This paper resolves two open questions from Fewster and Loveridge (`arXiv:2611.21768v1`, Rem. 2.10 & Sec. 5) and one open question from Heinonen, Lahti, Pellonpää, Pulmannová, and Ylinen (`arXiv:quant-ph/0212091`, Sec. 5.1):

1. **Exact operator-valued upper Beurling density criterion (`Theorem 2.1`, `Remark 2.2`).** For any translation-covariant POVM $E(\Delta) = W^* M_{\chi_\Delta} W$ on $\mathbb{R}^d$ with Toller–Fewster–Loveridge kernel $K(k) : \mathcal{H}(k) \to \widetilde{\mathcal{H}}$ on separable Hilbert spaces and $P(k) := K(k)K(k)^* \mathbf{1}_{V_{\mathrm{ac}}^+}(k)$, the small-ball localization ceiling $\beta(E) := \lim_{\varepsilon \to 0^+} \|E(B_\varepsilon(0))\|_{\mathrm{op}}$ equals the operator-valued upper Beurling density
   $$\beta(E) \;=\; \lim_{R\to\infty}\sup_{k\in\mathbb{R}^d}\|\overline{P}_R(k)\|_{\mathrm{op}} \;=\; \inf_{R>0}\sup_{k\in\mathbb{R}^d}\|\overline{P}_R(k)\|_{\mathrm{op}} \;=:\; \mathcal{D}_{\mathrm{op}}^+(P),$$
   with explicit finite-scale two-sided bounds. This extends the scalar and finite-matrix Beurling-density identity of Gabardo and Lai (`arXiv:1605.00165`) to infinite-dimensional separable fibres without compactness—a genuine separation, as an explicit kernel achieves $\sup_{\|u\|=1}\mathcal{D}^+(\langle u, Pu\rangle) = 0$ yet $\mathcal{D}_{\mathrm{op}}^+(P) = \beta(E) = 1$—and completely characterizes the intermediate regime left open by Fewster and Loveridge.
2. **Open-to-Borel equivalence on arbitrary homogeneous spaces $G/H$ and phase observables (`Theorem 3.1`, `Remark 3.2`, `Corollary 3.4`).** For any transitive covariant POVM $(U,E)$ on a locally compact second-countable homogeneous space $\Omega = G/H$ (arbitrary closed subgroup $H \le G$, not necessarily normal or unimodular), a direct proof via $L^1$-continuity of right Haar translation and Tonelli's theorem gives $\|E(B)\|_{\mathrm{op}} \ge \inf_{\emptyset \neq U \text{ open}} \|E(U)\|_{\mathrm{op}}$ for every Borel set $B$ with $m_G(\pi^{-1}(B)) > 0$. Combined with the null-set characterization (valid unconditionally in the forward direction, and in both directions whenever $U$ is strongly continuous or $\mathcal{H}$ is separable—sharp by two explicit counterexamples), open-set and full Borel definiteness coincide on every $G/H$ whenever $U$ is strongly continuous or $\mathcal{H}$ is separable, resolving `0212091` Sec. 5.1. Repeating the density theorem on $\mathbb{T}^d$ (dual $\mathbb{Z}^d$) yields an exact norm-$1$ criterion for all phase observables on $\mathbb{T}$ (including the $N=2$ necessary condition $\sup_n |\langle \xi_n, \xi_{n+1}\rangle| = 1$).
3. **Obstructions and complete regime separation (`Theorems 4.1, 4.2`, `Proposition 4.3`).** We sharpen Fewster–Loveridge's Hilbert–Schmidt obstruction (`Cor. 2.8`) to an operator-norm tail bound $\beta(E) \le \Lambda_\infty(K)^2$, derive a four-line Logvinenko–Sereda–Kovrijkine obstruction giving a sharp $0$–$1$ Borel dichotomy for scalar isometric kernels, and construct five explicit normalized examples separating every regime—including a $\dim\widetilde{\mathcal{H}} = 2$ kernel with $\beta(E) = 1/2$ yet $\|E((-L,L))\|_{\mathrm{op}} = 1$ for all $L > \omega$, proving the scalar $0$–$1$ dichotomy fails at $\dim\widetilde{\mathcal{H}} \ge 2$.
4. **Relativistic QFT bounded-mass-band no-go and Wick-power scaling-dimension law (`Proposition 5.1`, `Lemma 5.2`, `Theorem 5.3`, `Remark 5.4`).** Any kernel supported in a bounded invariant-mass band $m_1^2 \le k^2 \le m_2^2$ obeys the cube-density bound $|S \cap Q_L(p)| / |Q_L(p)| \le (m_2 - m_1)/L$ (best possible uniformly in $m_1, m_2$), forcing $\beta(E) = 0$ and $\|E(\Delta)\|_{\mathrm{op}} < 1$ on every bounded Borel set. Combining Kaiser's exact $2+1$D $n$-body phase-space formula (`arXiv:2510.08704v2`) with a general rescaled monotonicity lemma for all $d \ge 3, n \ge 2$, we prove that for a free massive real scalar field $\phi$ in $d = 1+s \ge 3$ dimensions, there exists $c > 0$ making the pointwise Wick-power kernel $W_{n,c}\psi(x) = c\langle\Omega|:\!\phi(x)^n\!:|\psi\rangle$ a globally bounded definite Poincaré-covariant POVM if and only if $n(d-2) = d$, i.e. $(d,n) \in \{(1+3,2), (1+2,3)\}$ (both satisfying Toller's Prop. 5, subnormalized via the vacuum projector, and collapsing to $\beta=0$ under any $L^1$/Schwartz smearing); $(1+2,2)$ realizes the subcritical intermediate regime ($\beta=0$, unbounded $h$, $\Lambda_\infty = 1$ at $c=2\sqrt{m}$ and $\Lambda_\infty = c/(2\sqrt{m}) < 1$ for $0 < c < 2\sqrt{m}$), and all remaining pairs are UV-unbounded.

## Authorship and checking record

- **Authors:** Clawd, Clayton Iggulden-Schnell, and Wren (alphabetical order; see repository root `README.md` for authorship disclosure).
- **Origin and drafting:** Originated, proved, and drafted by **Wren** (with Proposition 4.3(b'), the sharpened $(m_2-m_1)/L$ constant in Proposition 5.1(i), the 3-line Tonelli repair and two sharpness counterexamples in Theorem 3.1/Remark 3.2, and the $N=2$ phase condition in Corollary 3.4 contributed by **Clawd** during the blind audits).
- **Checking mode:** **Blind-checked** by **Clawd** in two proof-free statement-only packages (`11-data`/`11v` for base theorems `R1`–`R7`; `13-data`/`13v` for extensions `U1`–`U3`). In both rounds Clawd proved/tested every claim from primary sources and independent code and fixed his ruling file with a SHA-256 cryptographic hash (`504a10bb...cda3` on round 2) before opening Wren's proofs. Every claim held; the one proof fault caught (Weil's continuous $\rho$-function sentence in Theorem 3.1(ii)) and the prior-art catch (Kaiser `arXiv:2510.08704v2` for Lemma 5.2(ii)) were repaired in Draft 1. Clawd then performed a claim-by-claim proof-in-view read of `paper.tex` and independent copy rebuilds of `lean/` and `verify_qrf_norm1.py` (Letter 166; five faults and nine wording rows applied before release).

## Claim-by-claim verification map

| Claim | Proof | Who proved / checked (mode) | Lean 4 (`verification/lean/`) | Script (`verification/verify_qrf_norm1.py`) | What was NOT checked |
|---|---|---|---|---|---|
| **Theorem 2.1** (`thm:beurling`) | Analytic | Wren / Clawd (**blind**, round 1) | Not formalized | §§1–5 numerical checks on $\overline{P}_R$ and $\beta(E)$ | Infinite-dim Fourier/Bochner proof not in Lean |
| **Remark 2.2** (`rem:gl-sep`, $\infty$-dim separation from Gabardo–Lai) | Analytic | Wren / Clawd (**blind**, round 2) | Not formalized | §1 explicit $u_n \sim 1/n$ check (`0.0109` vs `1.0`) | Not in Lean |
| **Remark 2.3** (Smeared position control) | Analytic + cited | Wren / Clawd (**blind**, round 1) | N/A | N/A | `Carmeli2004` cited as background |
| **Theorem 3.1 & Remark 3.2** (General $G/H$ open-to-Borel & null-set sharpness) | Analytic | Wren (3-line Tonelli & counterexamples from Clawd) / Clawd (**blind**, round 2) | Final algebraic step `tonelli_floor_step` (`AlgebraicBounds.lean`) | §7 fat Cantor $C_6$ (`1 - 1.78e-8`) + affine $ax+b$ overlap (`0.9912`) | Haar-measure / Tonelli core not in Lean; full citation list of `0212091` not exhaustively audited |
| **Corollary 3.4** ($\mathbb{T}^d/\mathbb{Z}^d$ & phase observables) | Analytic | Wren / Clawd (**blind**, round 2; independent Fejér proof + 6 periodic numerical checks) | `phase_N2_char_poly`, `phase_N2_necessary` (`AlgebraicBounds.lean`) | §6 canonical, elementary (`<1e-12`), intermediate, sparse-block checks | Harmonic proof on $\mathbb{T}^d$ not in Lean |
| **Theorem 4.1** (Operator-norm tail bound) | Analytic | Wren / Clawd (**blind**, round 1) | `cauchy_schwarz_2d` (`AlgebraicBounds.lean`) | §1 bounded-$h$ control | Plancherel/Bochner steps not in Lean |
| **Theorem 4.2** (LSK obstruction & scalar dichotomy) | Analytic (conditional on Kovrijkine 2001 Thm. 1) | Wren / Clawd (**blind**, round 1) | `lsk_factor_identity` (`AlgebraicBounds.lean`) | §5 periodic comb ($\le 0.4000 < 1$) | Kovrijkine (2001) Thm. 1 cited without re-proving |
| **Proposition 4.3** (Five separating examples) | Analytic | Wren (with (b') from Clawd) / Clawd (**blind**, round 1) | `ex5c_char_poly` (`AlgebraicBounds.lean`) | §§1–4 aliasing-free checks (Ex. 5d exact to `<7e-15`) | Slepian/operator bounds not in Lean |
| **Proposition 5.1** (Mass-band bound $(m_2-m_1)/L$, best possible uniformly in $m_1, m_2$) | Analytic | Wren (sharpened constant from Clawd) / Clawd (**blind**, round 1) | **Core inequality formalized:** `vertical_slice_bound` & `vertical_slice_eq_iff` (`MassBand.lean`, arbitrary `OrderedCommRing`) | §8 max cube densities strictly $\le 1/L$ | Fubini integration step not in Lean |
| **Lemma 5.2** (Phase space: 2-body, Kaiser $2+1$D, general monotonicity $d \ge 3$) | Analytic + cited (Part (ii) = **Kaiser 2025**) | Wren / Clawd (**blind**, round 2; independent Dalitz & Laplace checks) | `kaiserDenom_closed_form`, `kaiserDenom_values`, `kallen_diff_identity`, `kallen_monotone_in_mu_sq` (`PhaseSpaceCoeff.lean`), `alpha_integrability_bound/eq_iff` (`ScalingLaw.lean`) | §9a–9e: 2-body, 4D quadrature for $\Phi_3^{(1+2)}$ (`<6.3e-12`), $\Phi_4^{(1+2)}$ (`<4e-15`), 3D Euclidean Laplace check (`<1.2e-13`) | MCT integration step not in Lean |
| **Theorem 5.3 & Remark 5.4** (Wick-power scaling-dimension law + $1+1$D remark) | Analytic | Wren / Clawd (**blind**, rounds 1 & 2; $3/(2\pi)$ constant checked by both) | **Diophantine classification formalized:** `wick_exponent_identity`, `wick_critical_iff`, `wick_subcritical_iff`, `wick_supercritical_iff`, `dim_1p1_always_subcritical` (`ScalingLaw.lean`) | §9d $\Phi_3^{(1+1)}$ threshold + $M^2\Phi_3/\log(M/m) \to 3/(2\pi)$ | Fock-space / Riemann–Lebesgue steps not in Lean |

## Lean-to-paper matching judgements

- **Scope:** The zero-dependency Lean 4 package (`leanprover/lean4:v4.33.1`, 5 modules, 19 theorems, ~2s build, zero `sorry`, zero warnings, zero `native_decide`, zero custom axioms) formalizes the **exact algebraic, Diophantine, and finite-dimensional matrix core** listed above (most notably the complete Diophantine classification of Theorem 5.3, the integrability-exponent bound and Källén monotonicity of Lemma 5.2(iii), Kaiser's denominator recursion, and the complete algebraic proof and equality characterization of Proposition 5.1(i)'s vertical-slice inequality over an abstract `OrderedCommRing`).
- **Not in Lean:** All infinite-dimensional Hilbert-space, Fourier-analytic, Haar-measure, and Dominated/Monotone Convergence arguments are proved analytically in `paper.tex` and blind cross-audited, not formalized in Mathlib. `OrderedCommRing` is witnessed non-vacuously on `Int`; that $\mathbb{R}$ is an instance is by hand.
- **Final review and independent rebuilds:** Clawd read `paper.tex` claim by claim (Letter 166: five faults and nine wording rows, applied), rebuilt `lean/` from a copy (green, standard axioms only; three docstring labels in `AlgebraicBounds.lean` were subsequently corrected from `Corollary 3.2` to `Corollary 3.4` with code and proofs unchanged and rebuilt cleanly by Wren), and re-ran the script from a copy (output identical but for runtime).

## Literature and prior-art record

- **Read directly at cited lines/equations in the original source:**
  - N. Kaiser, `arXiv:2510.08704v2` (`Sec. 2, Eqs. (2),(5)–(10)` for the exact $2+1$D $n$-body phase space — prior art for Lemma 5.2(ii); `Sec. 5, Eqs. (39),(41)–(43)` for $1+1$D 2- and 3-body phase spaces).
  - J.-P. Gabardo & C.-K. Lai, `arXiv:1605.00165` (`Thm. 1`, `Cor. 3`, `Lem. 9` / compactness at L739, `Thms. 14–15`, `Cor. 16`).
  - C. J. Fewster & L. Loveridge, `arXiv:2611.21768v1` (`Thm. 2.4`, `Thm. 2.7`, `Cor. 2.8`, `Rem. 2.10`, `Sec. 5`).
  - R. Beneduci & F. E. Schroeck Jr., `arXiv:1211.3847` (`Thm. 2.10`).
  - M. Toller, `arXiv:quant-ph/9805030` (`Prop. 5`, L697).
  - T. Heinonen *et al.*, `arXiv:quant-ph/0212091` (`Ex. 2`, `Prop. 8`, `Sec. 5.1` L552, `Sec. 5.2`, `Prop. 10`).
- **Background references NOT read directly at line level:**
  - Busch, Lahti, Pellonpää, Ylinen (2016), *Quantum Measurement* (monograph not held on disk).
  - Kovrijkine (2001) Thm. 1 (established harmonic-analysis theorem cited without re-proving; Theorem 4.2 is conditional on it).
  - Wightman (1962), Giannitrapani (1998), Carmeli–Heinonen–Toigo (`quant-ph/0403108`).
- **Search scope:** Targeted searches were conducted for `0212091` L552, operator-valued Beurling density for covariant POVMs, and $2+1$D $n$-body phase space; the complete citing list of `quant-ph/0212091` was not exhaustively audited.

## Files and reproduction

- `documentation/paper.pdf` and `documentation/paper.tex`: the paper (compiled with `tectonic -X compile paper.tex`).
- `verification/verify_qrf_norm1.py` and `verification/verify_qrf_norm1.out`: aliasing-free numerical verification (`python verify_qrf_norm1.py`, ~1s, requires `numpy`).
- `verification/lean/` and `verification/lean/verify_lean.out`: zero-dependency Lean 4 package (`lake build` inside `verification/lean/`).

## Licence

Text: CC BY 4.0. Code: Apache 2.0.
