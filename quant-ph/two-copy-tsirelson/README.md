# A joint quantum measure on two copies of the CHSH scenario gives Tsirelson's bound, and the converse of Henson's theorem fails

Clawd, Clayton Iggulden-Schnell and Wren, 20 November 2026. Clawd and Wren are AIs.

- `documentation/main.pdf`: the paper (8 pages), with its LaTeX source `main.tex`.
- `verification/lean/`: the machine-checked statements (Lean 4 with Mathlib).
- `verification/scripts/`: one script, with its recorded output, for every computed figure in the paper.
  `verification/scripts/FIGURES.md` maps each figure to the line that prints it.

## The result in plain words

Quantum theory allows correlations between distant measurements that no classical model allows, but not
every correlation that would still forbid signalling faster than light. The CHSH experiment measures this
with one number: at most 2 classically, at most 2√2 (Tsirelson's bound) in quantum theory, and up to 4 for
hypothetical "PR boxes". A long-standing question is which simple principle explains the limit 2√2.

Sorkin's quantum measure theory offers a candidate. A quantum measure assigns a non-negative weight to
every set of outcomes. It may break the ordinary rule that the weights of two disjoint sets add, but it
obeys the corresponding rule for three sets: there is no "third-order interference". Dowker, Henson and
Wallden turned this into a test of a correlation, called JQM. A correlation passes if one quantum measure
on the space of all outcomes of all measurements reproduces every observed probability. Henson proved that
JQM implies an older principle, Consistent Exclusivity: the probabilities of mutually exclusive events add
up to at most 1. Through Consistent Exclusivity, two copies of a box are excluded when the box's CHSH value
is above 2.883. Henson asked whether the converse holds, that is, whether JQM is the same as Consistent
Exclusivity.

It is not. We take two copies of the CHSH experiment, four parties, with the measurements Dowker, Henson
and Wallden include in which one party's choice of box depends on outcomes others have already obtained.
On that scenario we write down four explicit sets of 128 outcomes each. We prove an identity: for any
correlation and any candidate quantum measure, the four sets' weights add up to 2 − ⟨C₁C₂⟩/4, where C₁ and
C₂ are the CHSH expressions of the two copies. A quantum measure's weights are not negative, so
⟨C₁C₂⟩ ≤ 8. For two copies of the same box this is exactly Tsirelson's bound 2√2, where Consistent
Exclusivity gives only 2.883. For the standard family of noisy PR boxes the condition is also sufficient:
two copies pass JQM exactly when the box is quantum. Boxes with CHSH value between 2√2 and 2.883 obey
Consistent Exclusivity and fail JQM, so the converse of Henson's theorem is false.

The result needs those branching measurements. Without branching between the two experimenters' wings, the
identity does not hold (Table 1 of the paper), and neither does the 2.883 benchmark it improves on.

What it does not show: that JQM together with closure under composition singles out quantum theory in
general (that question of Dowker, Henson and Wallden stays open), or anything about three or more copies or
other Bell inequalities.

## How it was verified

There are three kinds of check:

- **Machine checks in Lean.** The identity, the bound and the non-existence half of the isotropic theorem
  are proved in Lean 4 with Mathlib, with no `sorry`, no extra axioms and no compiled code trusted.
- **Scripts.** Every count, rank and search in the paper is computed by a script filed here, re-run from
  scratch, with its output recorded beside it.
- **Reading by hand.** Each author read the arguments, and Wren also wrote independent code for the
  computations.

Wren's checks were made with the proof in view, not blind. Wren ruled the earlier draft claim by claim, and
on 20 November 2026 ruled the finished paper the same way. Every statement and every figure in the paper
text holds. Wren's tools were hand algebra, fresh code sharing nothing with Clawd's, and the cited papers
read directly. Wren did not open or rebuild the Lean files. The table below says, claim by claim, who did
what and what no one has checked.

### Claim by claim

Numbering is the paper's. "Clawd's code" and "Wren's code" mean separate programs written independently.
A Lean lane is a subagent that Clawd ran to write the formalisation. Clawd then read its statements cold
against the paper and rebuilt it from source. "Wren, by table" means Wren's claim-by-claim ruling on the
paper, with the proof in view.

| Claim | Who proved it | Who checked it, and how | Machine check | Left unchecked |
|---|---|---|---|---|
| Definition 1 (JQM on the two-copy scenario, with branching) | Dowker, Henson and Wallden; the premise that cross-wing branching pairs decohere is theirs | Clawd and Wren each read it against their paper | the definitions in `TwoCopy/Defs.lean` | — |
| Theorem 2 (the identity) | Clawd found it as an average over 16 relabelled sets. Wren found the four-set form and a proof without a computer, by averaging over a group of 64 relabellings. The paper's proof combines Wren's four-set form with Clawd's cancellation of the cross terms by character sums, which needs no group averaging | Wren, by table: the earlier draft and the paper's proof line by line, and the pointwise identity on all 256 points in Wren's code. Clawd re-ran every step with separate code. A blind refuter subagent tried to break it and failed | Lean proves it in a stronger form: D need not be symmetric, and P is any real function. The pointwise identity is also checked by `character_sets.py` | — |
| Remark 3 (the four sets are images of S₀ under the relabellings τⱼ; the group has order 64) | Clawd | Wren, by table, by hand and on all 256 points; Clawd's code | Lean: group order 64 and exactly four images. That the character definition gives the same sets is checked by `character_sets.py`, not in Lean | — |
| Corollary 4 (−8 ≤ ⟨C₁C₂⟩ ≤ 8; 2√2 for two copies of one box) | Clawd | Wren, by table | Lean, both bounds, for a behaviour | The sentence on products of two different boxes is not in Lean; it follows in one line from ⟨C₁C₂⟩ = S₁S₂ |
| Theorem 5, "only if" (q > 1/√2 has no JQM) | Clawd | Wren, by table | Lean | — |
| Theorem 5, "if" (an explicit quantum construction for q ≤ 1/√2) | Clawd | Wren, by table: read line by line, and coded independently, tested at q = 1/√2, 0.6, 0.3 and 0. Clawd's `theorem5_if.py` at q = 1/√2, 0.6 and 0.3. Both agree to errors near 10⁻¹⁶ | none; the numerical checks are floating point | Not in Lean |
| Theorem 6 (Consistent Exclusivity holds exactly when q ≤ (√10 − 1)/3) | The threshold and the ten-event set are due to Sainz et al., who report from computation that the value is optimal for two copies. The proof that the threshold is exact on the whole isotropic line is Clawd's | Wren audited it independently, in exact arithmetic in ℚ(√10) with no symmetry reduction, and ruled the paper's proof by table, rounding bound included. Clawd's search is in integers with a stated rounding bound | none | Not in Lean |
| Corollary 7 (the converse of Henson's theorem fails) | Theorems 5 and 6 together | as for those; Wren, by table | the Lean half is Theorem 5's "only if" | as for Theorems 5 "if" and 6 |
| Lemma 8 (a shared box drops out; the value of D(1,1)) | Clawd | Wren, by table | Lean: the first statement and D(1,1) | — |
| Lemma 9 (each M_k is ±4 everywhere) | Clawd | Wren, by table, by hand and on all 256 points; Clawd's code | Lean for M₀; `character_sets.py` for all four | — |
| Lemma 10 (which cross terms survive) | Clawd | Wren, by table: all 15 nonzero d by hand and in Wren's code | Not in Lean as stated (Lean proves Theorem 2 by its own route). The surviving d = 0111, 1110, 1111 and the vanishing character sums are checked by `lemma10_survivors.py` | — |
| Table 1 (the identity needs branching across wings) | Clawd | Wren, by table, recomputed independently over two prime fields; Clawd's ranks are exact over the integers (FLINT). All agree | `exact_ranks.py` | — |
| Section 6 counts (at least 18 of 45 pairs; at most four good-good events with product measurements; 3,072 and 6,144) | The count of 18 and the two incomparable classes (3,072 and 6,144) are Wren's. The at-most-four statement and the q = 1 statement are Clawd's | Clawd's and Wren's code. Wren audited the at-most-four statement in exact arithmetic | the scripts | — |
| Section 7 (the 640 five-event sets fall into two classes, 512 and 128, joined by no relabelling and no graph automorphism) | Wren found the split. Clawd added the graph invariant | Clawd re-ran it with the stabiliser computed inside the full relabelling group of 98,304; Wren re-ran it, invariant included, and checked both sets printed by Sainz et al. | `orbits.py`, `graph_inv.py` | — |

Clawd judged how the Lean statements correspond to the paper's, by reading the Lean files cold against the
text. Wren did not inspect the Lean. Section 8 of the paper lists what is formalised and what is not.

### Earlier work, and how it was searched

We looked for anyone who had already shown that a joint quantum measure on two copies gives Tsirelson's
bound, or that JQM is strictly stronger than Consistent Exclusivity:

- **Citing lists.** We took the full lists of papers citing Henson (32 in Semantic Scholar, 29 in
  OpenAlex), Dowker–Henson–Wallden (30) and Sainz et al. (40), fetched on 20 November 2026.
- **Nearest papers.** The nearest of these were searched in their TeX source for joint quantum measure,
  JQM, 2.883, two copies, third-order interference, quantum measure theory and decoherence functional, and
  every hit was read.
- **A plain-language web search.** One extended web search on the plain statement found nothing.

None of them states the two-copy bound. The nearest results are Boës and Navascués, who assume a product
rule for independent systems and closure for every number of copies, and Niestegge, who uses a different
definition of third-order interference with further axioms.

**Added 23 November 2026.** The two-copy problem itself was posed before us. Paul Boës, in his MSc
dissertation "Closure and Stability in Quantum Measure Theory" (Imperial College London, 2014), section 6.3,
sought the largest noisy PR box whose two-copy behaviour admits a decoherence functional, with branching, by
a symmetrised linear programme. He reported no improvement on the local-orthogonality value 0.721 and left
open whether the joint quantum measure condition recovers Tsirelson's bound and whether it coincides with
local orthogonality. Theorem 5 and Corollary 7 of the paper answer both for two copies on the isotropic line.
The dissertation is not on any of the citing lists above; Clawd found it on 22 November and read section 6.3
on 23 November, from the page images. The paper now cites it in "Earlier work".

The search has limits. A paper that proves this result without citing Henson, Dowker–Henson–Wallden or
Sainz et al. would not be found this way, and some titles on those lists were read by title only.

### References not read first-hand

Sorkin (1994), Dowker–Henson–Wallden, Henson, Sainz et al., Boës–Navascués and Niestegge were read from
their arXiv sources, at the lines the paper relies on. Popescu–Rohrlich (1994) and Cirel'son (1980) are cited
only for the PR box and for Tsirelson's bound, which are standard. We did not read them first-hand for this
paper.

## Licence

Text: CC BY 4.0. Code: Apache 2.0.
