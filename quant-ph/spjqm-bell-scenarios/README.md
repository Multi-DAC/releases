# Strongly Positive Joint Quantum Measures in Bell Scenarios: Exact Characterization and Obstructions to Persistence-of-Zero Patching

Clawd, Clayton Iggulden-Schnell, and Wren, 20 November 2026. Clawd and Wren are AIs.

**Primary classification:** `quant-ph`
**Cross-lists:** `math-ph`

## Plain-language statement of the result

In histories-based quantum mechanics (quantum measure theory), a *strongly positive joint quantum measure* ($\mathrm{SPJQM}$) represents an experimental probability table by a positive semidefinite decoherence functional on all fine-grained histories. Dowker, Henson, and Wallden (`arXiv:1311.6287`) proved the bracket
$$\tilde{Q} = \mathrm{SPJQM}_b \;\subseteq\; \mathrm{SPJQM} \;\subseteq\; Q^1,$$
where $\mathrm{SPJQM}_b$ adds adaptive branching measurements, $\tilde{Q}$ is the almost-quantum set, and $Q^1$ is level~$1$ of the Navascués–Pironio–Acín (NPA) hierarchy, leaving open where $\mathrm{SPJQM}$ sits and whether branching strictly strengthens the condition.

We settle this completely for Bell scenarios and trace the consequences for Dowker and Sorkin's (`arXiv:2305.16828`) Persistence-of-Zero ($\mathrm{PoZ}$) and Lack-of-Novelty ($\mathrm{LoN}$) framework:

1. **Bipartite equality (`Theorem 3.2`).** In every finite bipartite Bell scenario (any numbers of settings and outcomes per party), $\mathrm{SPJQM} = Q^1$. Combined with known separations $\tilde{Q} \subsetneq Q^1$, this shows $\mathrm{SPJQM}_b \subsetneq \mathrm{SPJQM}$ and immediately transfers bipartite $Q^1 \setminus \tilde{Q}$ constructions (such as Gazi *et al.* `arXiv:2610.41448v1`) into $\mathrm{SPJQM}$.
2. **Multipartite separation and exact $n$-party characterization (`Theorems 4.1, 4.3`, `Lemma 4.2`).** For $n \ge 3$ parties, $\mathrm{SPJQM}_n \subsetneq Q^1_n$: the three-party box $a \oplus b \oplus c = xyz$ lies in $Q^1_3 \setminus \mathrm{SPJQM}_3$. For every $n \ge 2$, $\mathrm{SPJQM}_n$ is characterized exactly by Hilbert-space vectors indexed by all non-maximal ($\le n-1$-party) outcomes via an $n$-party vector-valued marginal lifting lemma.
3. **Diagnosis and refutation of Dowker–Sorkin's converse theorem (`Theorems 5.1, 5.2`).** Dowker and Sorkin's Quantum Patching Theorem yields exactly $\mathrm{SPJQM}$ (without branching). Their printed converse construction obeys $\mathrm{PoZ}$ if and only if the nonzero atom vectors of the joint functional are linearly independent (failing already on a rank-$4$ quantum functional $J_1$, while mixing with the uniform behaviour yields $\mathrm{relint}(Q^1) \subseteq \mathcal{P}_{\mathrm{PoZ+LoN}} \subseteq Q^1$ and $\sup_{\mathcal{P}_{\mathrm{PoZ+LoN}}} = \max_{Q^1}$ on every linear Bell functional). Yet the converse theorem itself is false as stated: an explicit rational rank-$4$ behaviour $\mathrm{Hardy}_{1/8} \in \mathrm{SPJQM} = Q^1$ admits no $\mathrm{PoZ}+\mathrm{LoN}$ model on any finite past event algebra (or any infinite algebra satisfying algebraic $\mathrm{LoN}$), proving that $\mathcal{P}_{\mathrm{PoZ+LoN}}$ is not closed.
4. **Elementary quantitative $Q^{1+AB}$ inequalities and minimal three-party chain (`Propositions 6.1, 6.2`, `Corollaries 6.3, 6.5`).** $\mathrm{Hardy}_{1/8}$ is the unique $Q^1$-maximizer on its four-zero face $F_4$. An elementary quantitative $Q^{1+AB}$ bound derived directly from DHW additivity and cross-context orthogonality gives a computer-free witness $\mathrm{Hardy}_{1/8} \in \mathrm{SPJQM} \setminus \mathrm{SPJQM}_b$ in the $(2,2,2)$ scenario, proves that the interior segment $P_\varepsilon = (1-\varepsilon)\mathrm{Hardy}_{1/8} + \varepsilon U$ lies in $\mathcal{P}_{\mathrm{PoZ+LoN}} \setminus Q^{1+AB}$ for all $0 < \varepsilon < 1/19$ (and $0 < \varepsilon < (11+8\sqrt{2})^{-1}$ via the symmetric square-root bound), and yields the minimal $(3,2,2)$ strict chain $\tilde{Q}_3 = \mathrm{SPJQM}_{b,3} \subsetneq \mathrm{SPJQM}_3 \subsetneq Q^1_3$.

## Authorship and checking record

- **Authors:** Clawd, Clayton Iggulden-Schnell, and Wren (alphabetical order; see repository root `README.md` for authorship disclosure).
- **Origin and drafting:** Co-originated by **Clawd** and **Wren**; paper, Python verification script, and Lean 4 package drafted by **Wren**.
  - **Clawd** originated Theorems 3.2, 4.1, 4.3, Lemma 4.2, Theorem 5.1, and the symmetric bound of Proposition 6.2.
  - **Wren** originated Theorem 5.2 (`Hardy-1/8`), Proposition 6.1, the Pythagoras sharpening of Proposition 6.2, and Corollaries 6.3 and 6.5.
- **Checking modes:**
  - *Checked with written proof in view (primary sources + independent zero-shared-code scripts):* Theorems 3.2, 4.1, 4.3, Lemma 4.2, and Theorem 5.1 (checked by Wren); Theorem 5.2 finite-past case (checked by Clawd).
  - *Checked blind (statements only, SHA-256 hash-sealed ruling recorded before opening proofs):* Proposition 6.1, Corollary 6.3 (up to $(11+8\sqrt{2})^{-1}$), and Corollary 6.5 (checked by Clawd; the 1-line Pythagoras sharpening to $1/19$ was checked by hand by Clawd).
  - *Lean 4 rebuild:* Clawd independently copied the Lean sources to a clean folder and verified `lake build` (`leanprover/lean4:v4.33.1`, 9 jobs, exit 0), zero `sorry`/`axiom`/`native_decide`/`#eval`, and `#print axioms`.

## Claim-by-claim verification map

| Claim | Proof | Who proved / checked (mode) | Lean 4 (`verification/lean/`) | Script (`verification/verify_theorem101.py`) | What was NOT checked |
|---|---|---|---|---|---|
| **Lemma 3.1** (Gram isometry) | Analytic | Clawd / Wren (proof-in-view); read by Clawd | Not formalized | Exercised in §5 Gazi check | No Lean proof |
| **Theorem 3.2** ($\mathrm{SPJQM} = Q^1$ bipartite) | Analytic | Clawd / Wren (proof-in-view); Steps 1–5 read line-by-line by Clawd | Marginal-recovery step (`lem_lifting_step_a`, not read by Clawd); Gram dilation analytic | §1 exact `Fraction` S4 Claim A; §5 Gazi $p^{\mathrm{Vec}}_5$ (`1.55e-15`, drop-$K_{xy}$ control `1.84e-2`) | Gram dilation not in Lean |
| **Remark 3.3** (Gazi *et al.*; Das Bhowmik–Banik) | Analytic + cited | Clawd (first-hand read) / Wren (§5 check) | N/A | §5 rank-11 $p^{\mathrm{Vec}}_5$ check | Citations read by Clawd only |
| **Theorem 4.1** ($a\oplus b\oplus c = xyz \in Q^1_3 \setminus \mathrm{SPJQM}_3$) | Analytic | Clawd / Wren (proof-in-view); Clawd read & rebuilt Lean | **Formalized** (`ThreePartyCounter`: marginals, correlators, `thm_3party_no_spjqm`; read & rebuilt by Clawd) | §4 exact `Fraction` check | Real/complex inner-product instance of `OrderedCommRing`/`RealInnerProd` by hand |
| **Lemma 4.2** ($n$-party lifting) | Analytic | Clawd / Wren (proof-in-view); Steps (a)–(b) read by Clawd | **Step (a) formalized** (`BipartiteLift.lem_lifting_step_a`, arbitrary $m \ge 1$ & `dims`, not read by Clawd); Step (b) induction analytic | §1 exact `Fraction` `ind == unw == marginals` up to 4 parties (`|Xi|=96`) | Unwound closed form checked by Clawd for $n=2$ only (and up to $n=4$ in Wren's script) |
| **Theorem 4.3** ($n$-party $\mathrm{SPJQM}_n$) | Analytic | Clawd / Wren (proof-in-view); both directions read by Clawd | Marginal kernel (`lem_lifting_step_a`); dilation analytic | §1 + non-quantum 3- and 4-party tests | Dilation not in Lean |
| **Theorem 5.1** (DS23 printed construction, `J1`/`J2`/`J3`, $\mathrm{relint}$ sandwich) | Analytic + exact finite certificates | Clawd / Wren (proof-in-view); exact arithmetic by both | Partial (`j3_delta_context_zero`, `j3_delta_branching_sum` `-2/125`, read & rebuilt by Clawd); (a)–(c) analytic | §3 exact `Fraction` ranks & PoZ violations on `J1, J2, J3` | Sign formula for $\Delta$ not compared entry-by-entry with Clawd's outer product; (a)–(c) not in Lean |
| **Theorem 5.2** (`Hardy-1/8` refutes DS23 converse; $\mathcal{P}_{\mathrm{PoZ+LoN}}$ not closed) | Analytic + exact finite certificate | Wren / Clawd (finite past **HOLDS**; algebraic-$\mathrm{LoN}$ infinite extension added in Draft 2, unchecked by Clawd) | **Partial:** `hardy18_spjqm_exact` (64/64 same-context entries via `decide`); `thm_ds23_no_model` (4-line operator deduction from additivity + 2 commutators). Read & rebuilt by Clawd | §2 exact `Fraction` 16x16 rank-4 Gram matrix | PSD of `Hardy-1/8` by Gram construction (by hand); reduction from PoZ+LoN to operator hypotheses by hand; algebraic-$\mathrm{LoN}$ infinite case unchecked by Clawd; topological-closure unbounded-operator case not claimed |
| **Proposition 6.1** ($Q^1$ iff on $F_4$, unique max $1/8$) | Analytic | Wren / Clawd (**blind**) | **AM–GM bound & uniqueness formalized** (`prop_face_max_bound`, `prop_face_max_unique`, not read by Clawd); $Q^1$ iff analytic | §6a 400 random 5x5 $Q^1$ checks on $F_4$ | $Q^1$ iff reduction not in Lean |
| **Proposition 6.2** (Quantitative $Q^{1+AB}$ bounds) | Analytic | Symmetric bound proved independently by Clawd; Pythagoras step by Wren, checked by hand by Clawd | **Both bounds formalized** (`prop_qab_vector_identity`, `prop_qab_pythagoras_bound`, `prop_qab_sqrt_bound`; read & rebuilt by Clawd) | §6b exact affine check; §6c 80 boundary $9\times 9$ SDP points | Real/complex Hilbert-space instance of Lean classes by hand |
| **Corollary 6.3** (Direct witness & $P_\varepsilon \in \mathcal{P}_{\mathrm{PoZ+LoN}} \setminus Q^{1+AB}$) | Analytic (thresholds); numerical SDP (exit bracket) | Wren / Clawd (**blind** up to $(11+8\sqrt{2})^{-1}$, by hand for $1/19$); numerical SDP bracket **Wren only** | **Formalized** (`cor_hardy_face_zero`, `cor_peps_violates_pyth`; read & rebuilt by Clawd) | §2 exact `Fraction` at $\varepsilon=1/23, 1/20$; §6d numerical $9\times 9$ SDP exit in `(0.0655, 0.0656)` | Numerical SDP exit `(0.0655, 0.0656)` computed by Wren only, unchecked by Clawd |
| **Corollary 6.5** (Minimal $(3,2,2)$ strict chain) | Analytic + exact finite certificate | Wren / Clawd (**blind**) | **Left 64-atom check (`hardy18_3party_spjqm_exact`, 512/512 via `decide`, read & rebuilt by Clawd) + right counterexample** | §2 exact `Fraction` 64-atom rank-4 Gram matrix (`512/512`) | Composition closure & restriction to $AB$ by hand |

## Lean-to-paper matching judgements

Recorded jointly by Wren (author of the Lean package) and Clawd (independent rebuilder):
- **Toolchain & axioms:** `leanprover/lean4:v4.33.1`, zero external dependencies, `lake build` succeeds in ~2s (9 jobs), zero `sorry`, zero `native_decide`, zero custom axioms (`#print axioms` reports only `propext`, `Quot.sound`, `Classical.choice`).
- **Non-vacuousness & real/complex instantiation:** The algebraic typeclasses (`AddCommGrp`, `OrderedCommRing`, `RealInnerProd`) are witnessed non-vacuously on `Int` (`Instances.lean`). That $\mathbb{R}$ and real inner-product spaces satisfy these axioms—and that complex Hilbert spaces satisfy them via $\operatorname{Re}\langle\cdot,\cdot\rangle$ (preserving norms and orthogonality)—is verified by hand. Square roots enter as algebraic hypotheses `s ≥ 0 ∧ s * s = p`.
- **Script independence:** `verify_theorem101.py` was written and run by Wren; Clawd read its output (`verify_theorem101.out`) but did not read or execute the script.

## Literature and prior-art record

- **Read directly at cited lines by both Wren and Clawd:** Dowker–Henson–Wallden (`arXiv:1311.6287`, `SPJDF0.22.tex` numbering verified), Dowker–Sorkin (`arXiv:2305.16828`), Dowker–Wilkes (`arXiv:2011.06120`), Henson–Sainz (`arXiv:1501.06052`).
- **Read directly by Clawd only:** Gazi *et al.* (`arXiv:2610.41448v1`), Das Bhowmik–Banik (2026 Research Square preprint).
- **Read directly by Wren only (unread by Clawd):** Rai *et al.* (`arXiv:1812.06057`), Araújo (`arXiv:2302.03529`).
- **Background references not line-audited in this study:** Sorkin (1994), Dowker–Johnston–Sorkin (2010), NPA (2007, 2008), Navascués–Wunderlich (2010), Navascués *et al.* (2015, `arXiv:1403.4621`).

## Files and reproduction

- `documentation/paper.pdf` and `documentation/paper.tex`: the paper (compiled with `tectonic -X compile paper.tex`).
- `verification/verify_theorem101.py` and `verification/verify_theorem101.out`: exact rational (`fractions.Fraction`) and numerical SDP checks (`python verify_theorem101.py`, requires `numpy` and `scipy`).
- `verification/lean/` and `verification/lean/verify_lean.out`: zero-dependency Lean 4 package (`lake build` inside `verification/lean/`).

## Licence

Text: CC BY 4.0. Code: Apache 2.0.
