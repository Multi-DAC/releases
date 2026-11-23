# Blind refuter round on the proof of Theorem 13, ruled by table (23 November 2026)

What this is. Before the paper was drafted, the proof of the all-copy theorem (now Lemma 12, Theorem 13 and the
parity theorem, Proposition 7) was cut into 19 steps, R1 to R19, in the order of the proof. Two refuters (AI
instances with fresh context, given the 19-step packet only, no web access, working separately) were asked to break
each step, by argument and by their own code. Clawd then ruled each row. The table below is that ruling, transcribed
from the working record with internal file references removed. The refuters read the pre-draft packet, not the text
of the paper; the paper's text was audited separately (see the README).

Row R16 is the one that did not hold as first written. The repair is the block-spanned statement of Proposition 7,
which is what the paper prints; see Remark 14.

| Row | A | B | Ruling | Deciding line |
|---|---|---|---|---|
| R1 | HOLDS | HOLDS | HOLDS | flips send Walsh monomials to +- themselves; both ran all 16 flips of all 90 sets |
| R2 | HOLDS | HOLDS | HOLDS | 9^n distinct characters of Gamma_n; A: deviation 0 over 90 sets |
| R3 | HOLDS | HOLDS | HOLDS | monomial = +- sum over one measurement, same-setting (ii) kills cross terms; A 0 (control 3.17), B 7.9e-16 |
| R4 | HOLDS | HOLDS | HOLDS | Parseval; A n=2 non-product non-symmetric K 1.7e-13 (control 18.3); B 1.1e-13 |
| R5 | HOLDS | HOLDS | HOLDS | exact rank 56 of 112 x 81 (both, sympy/rational); 25 = 17 + 8; "must" shown by rank 98 vs 62 (float, sampled sets) |
| R6 | HOLDS, unused | HOLDS, unused | HOLDS, not used by the proof | witness reproduced: 0.37402 real vs 0.37598 Hermitian on a 96-point set |
| R7 | HOLDS | HOLDS | HOLDS | event vectors rank 9; intersection over slots = tensor power; B adds: marginals keep (i) because traced off-diagonal events are exclusive |
| R8 | HOLDS (argument) | HOLDS (argument) | HOLDS | closed subspace constraints, bounded by R4; consistency exact per n, passes to the diagonal limit |
| R9 | HOLDS | HOLDS | HOLDS | P_S P_S' = P_{SxS'}; P_Omega = 1 on H |
| R10 | HOLDS | HOLDS | HOLDS | degree-N form zero on {G(1,1)=1} is zero on G_all; powers span Sym^N(G_all), K_N lies there; the transpose-sector detour is correct and unneeded (both) |
| R11 | HOLDS (argument) | HOLDS (argument) | HOLDS | Tychonoff; cases L'(t) = 0, 1, interior; P_S/abs(S) span A_0; no squares needed |
| R12 | HOLDS (argument) | HOLDS (argument) | HOLDS | Krein-Milman + closed X_T; second moment uses (i_2) on e x e; 16 integrands at once |
| R13 | HOLDS | HOLDS | HOLDS | A_sk(1,u) = 2.9e-16 / 4.0e-16 over the 8 skew basis forms; g identities on all four odd C |
| R14 | HOLDS | HOLDS | HOLDS | F = g Q(S,S') on actual 256-point sets, 6.0e-14 / 3.6e-14; Q_chi well defined because A_0 is a function algebra |
| R15 | HOLDS | HOLDS | HOLDS, one line to add | block identity at k=2 on actual three-copy sets 1.3e-15 / 8.9e-16; ADD: Q_chi(A,B) = chi(g^2)/chi(g) = chi(g) |
| R16 | FAILS-REPAIRABLE | FAILS-REPAIRABLE | FAILS-REPAIRABLE, repair adopted | the parity theorem as then cited asked Q^{(x)n} >= 0 on EVERY spanned set; Step 4 delivers it only on one-copy sets and block sets, and chi(P_E) = Q_chi^{(x)n}(E,E) is FALSE off block sets (A: fails on 560 of 3,232 two-copy spanned sets, max 6.25e-4). Repair: state the parity theorem with the hypotheses its proof consumes (one-copy sets, block sets on the same Fine pair). Wren's draft of the same hour already stated it so, independently. |
| R17 | HOLDS | HOLDS | HOLDS | A: non-PSD LP form at <C> = 3.2, k=1 brute force 0.1237043493 = formula, negative from k=2, explicit 2,304-point set at -1.3626; B: three PR mixtures, sum = formula at k=1,2 |
| R18 | HOLDS | HOLDS | HOLDS | Q_chi PSD with (i), (ii) is the GIVEN characterisation |
| R19 | HOLDS, wording notes | HOLDS, wording notes | HOLDS, wording | "NO product rule" means K_n need not be K_1^{(x)n}; the product BEHAVIOUR p^{(x)2} (mu_2(e x e) = p(e)^2) is used and is part of the definition. Local branch: point evaluations. |

**Overall.** A: REPAIRABLE. B: STANDS AS FAR AS I CAN TELL (R16 restated). Ruling: the proof stands with the parity
theorem stated in its block-set form. No mathematical step failed in either reading.

**Both refuters named the same weakest point:** the join R15-R16. The proof works because the parity proof's
negative sets are block sets on one Fine pair, where the skew sector cancels. Rows checked by argument only, with no
finite computation behind them: R8, R11, R12 (compactness, Krein-Milman, the limit), and the definitions.

**Definitions.** That slot-wise 56-pair (ii) at n copies is exactly Dowker, Henson and Wallden's branching condition
on the 2n-wing composition (Lemma 2 of the paper) was checked separately against their TeX source
(`SPJDF0.22.tex`, lemma `l:branch`, L1017-1025), not by the refuters. The theorem asks positivity only on spanned
sets, which is weaker than positivity on all subsets, so the theorem as proved is the stronger statement.
