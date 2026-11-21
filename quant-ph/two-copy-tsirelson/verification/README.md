# Verification

Two parts: `lean/`, the machine-checked statements, and `scripts/`, the computations behind every figure in
the paper.

## lean/

A Lake project: Lean 4 (`lean-toolchain`: `leanprover/lean4:v4.33.1`) with Mathlib at tag `v4.33.1`
(commit `0df444a360eaa60ab8c11dca51a86af692955474` in our build). Nine modules, 1,638 lines, in `TwoCopy/`;
`TwoCopy.lean` imports them all. Build with `lake exe cache get` then `lake build`.

What is formalised, and what is not, is listed in Section 8 of the paper. In short: the finite facts about the
four sets, the pointwise identity, Lemma 8's first statement, Theorem 2 for any bilinear D on W with (i) and
(ii) (not assumed symmetric, and any real P), both bounds of Corollary 4, the "only if" half of Theorem 5, and
the number of events locally orthogonal to each event. In Lean the four sets are defined as the images of S_0
under the relabellings tau_j; that they equal the sets defined by the characters chi_k is checked by
`scripts/character_sets.py`.

Finite facts are proved by `decide +kernel`: the kernel evaluates them, and no compiled code is trusted. The
sources use no `sorry`, `admit`, `axiom` declaration or `native_decide` (the words occur only in comments). `TwoCopy/Axioms.lean` prints the axioms of
43 statements; every one depends on `propext`, `Classical.choice` and `Quot.sound` only.

`mutants/Mutants.lean` holds four deliberately false versions of the finite facts (the group order 63, the
identity for one set alone, one sign changed, M always 4). Each must fail to build.

`evidence/rebuild-2026-11-20-1703.log` is the log of a build from source on 2026-11-20: all nine modules exit
0 without warnings, 43 axiom lines as above, and the four mutants each rejected. The sources here are
byte-identical to the ones built there. Cost on that machine: `Finite` 139 s and 12.4 GB of memory at its peak,
`Bounds` 97 s and 15.1 GB; a machine with 16 GB of memory may not build those two modules.

## scripts/

`FIGURES.md` maps each computed figure in the paper to the script and the output line that prints it. Python 3
with numpy and python-flint; `base.py` builds the scenario (points, events, local orthogonality, the
monomial basis of W, the isotropic boxes) and the other scripts import it from the same folder. Each script's
output from a fresh run is recorded beside it as `.out`.
