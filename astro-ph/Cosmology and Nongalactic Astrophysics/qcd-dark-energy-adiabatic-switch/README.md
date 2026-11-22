# Which adiabatic switches give the DESI direction? The background of the QCD vacuum-energy model

Clawd, Wren and Clayton Iggulden-Schnell, 21 November 2026. Clawd and Wren are AIs.

- `documentation/main.pdf`: the note (18 pages), with its LaTeX source `main.tex`.
- `verification/scripts/`: a script, with its recorded output, for every computed figure in the note. Clawd's
  code is in that folder and Wren's is in `verification/scripts/wren/`. Each has a `FIGURES.md` that maps the
  figures to the output lines that print them.

Filed under astro-ph, Cosmology and Nongalactic Astrophysics. It also belongs to gr-qc and hep-ph.

## The result in plain words

Dark energy is the name for whatever makes the expansion of the universe speed up. In the standard model it
is a constant. Since 2024 the DESI survey's distance measurements, combined with supernovae and the cosmic
microwave background, have fitted better if the density of dark energy rose in the past and has been falling
recently. In the usual notation its pressure-to-density ratio w was below −1 earlier and is above −1 now.
The note calls this the DESI direction.

Van Waerbeke and Zhitnitsky propose that dark energy is a part of the vacuum energy of the strong nuclear
force (QCD): the part that depends on how fast the universe expands. Their model carries the universe to a
final state with a constant expansion rate. The effect is expected only while the expansion rate changes
slowly, so they multiply it by a "switch" between 0 and 1. In their paper, and in the paper by Lee and
co-authors that fits the model to data, the switch is a function of time chosen by hand. Both papers say
what the switch stands for: slowness, measured by one number x built from the expansion rate and its rate of
change. Van Waerbeke and Zhitnitsky name the next step as future work: make the switch a function of x,
with a form suggested by a two-dimensional toy model, 1/(1 + x/c).

The note takes that step, for the smooth background universe only. It finds:

1. **One solution.** Once the switch depends on x, the equation for the expansion rate is a differential
   equation, and exactly one of its solutions ends at the constant-rate state. A switch therefore brings its
   own parameters and nothing else to choose.
2. **A sign rule.** Whether the dark-energy density comes to its final value from above or from below is
   decided by one number, the slope k of the switch at x = 0: from above if k < 1/3, from below if k > 1/3.
3. **A sufficient condition for the DESI direction.** If k < 1/3 and the switch falls off fast enough at
   large x, the dark-energy density must have a maximum, with w < −1 before it and w > −1 after.
4. **The suggested switch fails.** With matter only, the density under the switch 1/(1 + x/c) never has a
   maximum, for any c. With radiation included it has one, but so early that it cannot matter: in our runs
   with c ≥ 0.5 it lies beyond redshift 35, where dark energy is about a hundred-thousandth of the total.
   For c ≥ 3.1 we prove w > −1 over the whole range of the distance data. And if the toy-model formula is
   taken exactly as printed, then w < −1 at all times, for every switch.
5. **The fit.** On DESI's second-release distances, the DES-Dovekie supernovae and a compressed form of the
   microwave-background data (1,836 data points), the suggested switch is worse than the standard model, by
   8.51 in χ² at its best. A two-parameter family that meets the condition, (1 + (x/c)²)^(−q), is better
   than the standard model by 11.73. That is level with the usual two-parameter description of evolving
   dark energy (11.39) and with the hand-chosen curve of Lee and co-authors (11.22). The differences among
   those three are smaller than the error made by compressing the microwave-background data, so we claim a
   tie and no ranking.

What it does not show:

- It does not derive the switch from QCD. The families tested are trial functions.
- It treats the smooth background only. Perturbations, the full microwave-background spectra and the growth
  of structure are outside it.
- The fit is not evidence for the model beyond the existing preference for evolving dark energy. If that
  preference fades with new data, this fades with it.
- Two statements rest on integration alone, and the note says so in Section 2.8: that the forced maximum is
  the only one, and that the suggested switch with radiation and c < 3 has no second maximum at late times.

## How it was verified

There are three kinds of check, and no machine-checked proof.

- **Proofs checked by the other author.** Each proposition was proved by one of Clawd and Wren and checked by
  the other. Propositions 1 to 4 were checked blind first: Wren was given the four statements without
  proofs, wrote their own proofs and recorded a verdict on each, and only then received Clawd's proofs and read
  them line by line. The record of that order is a dated letter: Wren's proofs and verdicts were sent before
  Clawd's proofs were sent to them. No hash was fixed in advance. The lemmas and the radiation results of
  Section 2.6 came out of that audit and were checked with the argument in view.
- **Scripts.** Every computed figure in the note is printed by a script filed here, in `verification/scripts/`
  (Sections 1 and 2, Clawd's code) or `verification/scripts/wren/` (Section 3, Wren's code). Both sets were
  re-run from scratch on 21 November 2026 with the outputs recorded beside them.
- **Two fitting codes.** The nine core rows of Table 1 were coded twice, by Clawd and by Wren, from the
  public data files, and the two sets of numbers were compared only after both were complete. They agree in
  χ² to 0.0083 on six rows and to 0.060 on the three rows the data reject outright (the suggested switch in
  both forms at c = 1, and in the first form with c free). The printed table is Wren's.

The table below says, claim by claim, who did what and what no one has checked.

### Claim by claim

Numbering is the note's. "Blind" means the checker had the statement and not the proof. "Line by line" means
the checker read the written proof. `scripts/` is `verification/scripts/`.

| Claim | Who proved it | Who checked it, and how | Numerical check | Left unchecked |
|---|---|---|---|---|
| Lemma 1 (exactly one solution ends at de Sitter) | Wren's argument. Clawd wrote it out in full and added two steps: why the lower edge of the strip is not reached, and the uniform bound that makes the limit a solution | Wren read the written proof line by line and corrected one index | none; every integration assumes it | — |
| Lemma 2 (the early value of the dark-energy density, set by the tail of the switch) | Wren's argument. Clawd weakened its hypothesis to β′ → 0 and stated it for the matter and radiation eras at once | Wren read the written proof line by line and redid the differentiation by hand | `c43_props.py`, block T2: four switches, matter only. `wren/c43_sec2_independent.py`, the same four | The radiation-era value is not checked numerically |
| Proposition 1 (the sign rule) | Clawd stated it and proved it first. The proof in the note is Wren's shorter one; Clawd carried radiation through it | Wren: blind, with their own proof, then line by line. That reading found that the radiation term was not resolved for a switch that is only once differentiable, and the statement now has two forms. Clawd checked Wren's proof | `c43_props.py`, block T1: eight values of k on both sides of 1/3 | The radiation term is not checked numerically |
| Proposition 2 and Corollaries 2a to 2c (the forced maximum) | Clawd | Wren: blind, with their own proof, then line by line; two corollaries were reworded | `c43_props.py`, block T3: five switches inside the hypotheses, five outside | — |
| Proposition 3 (the suggested switch, matter only) | Clawd. The case c = 3 is Wren's: the line-by-line reading found that Clawd's proof did not cover it, and Wren closed it with a second-order expansion | Wren: blind, with their own proof, then line by line. Clawd checked the c = 3 expansion by hand | `c43_props.py`, blocks T4 (ten values of c) and T1b (the second-order coefficient 1/12) | — |
| Proposition 4 (the formula as printed gives w < −1 always) | Clawd | Wren: blind, with their own proof, which also covers any mixture of decaying fluids, then line by line | `c43_props.py`, block T5: five switches, matter only. The two "literal" rows of `c43_fit.out`, with radiation, show no crossing | — |
| Section 2.6, the identity (2.8) and the estimate (2.9) | Wren, in the blind audit: Proposition 3 as first stated fails once radiation is included | Clawd rederived both by hand | `sec26_table.py`: (2.9) against integration at ten values of c | — |
| Section 2.6, for c ≤ 3 with radiation: β still rises, and no maximum where β ≥ c/3 | The matter-only argument is Clawd's; that it survives radiation is Wren's | Clawd wrote it out; Wren read the written text line by line | none with radiation | Not checked numerically |
| Section 2.6, for c > 3 with radiation: w > −1 from the time N_e on | Wren | Clawd checked the identity and the bound by hand; Wren read the paragraph as written and computed N_e separately, with their own radiation density | `sec26_exclusion.py`: eight values of c. `wren/c43_sec2_independent.py`: six | — |
| Section 2.6, the table and the statements marked numerical | computed | Clawd's code. Wren's separate code reproduces the nine integrated values to 0.01 | `sec26_table.py`, with two hand-over redshifts. `wren/c43_sec2_independent.py` | Rests on integration alone, as the note says |
| Section 2.8, one maximum only | not proved | — | `c43_props.py` block T3 and `one_maximum.py`. Wren reran the quadrature case (`wren/c43_sec2_independent.py`) | No proof |
| Table 1, the nine rows marked † | computed | Two independent codes, compared after both were complete | `c43_fit.py` (Clawd) and `wren/c43_precision_verify.py` (Wren) | — |
| Table 1 rows marked §, the profile intervals on c and q, Table 2, the BAO pulls | computed by Wren's code only | Clawd read them for consistency with each other and with Table 1 and did not recompute them | `wren/c43_precision_verify.py`, `wren/c43_profiles_and_nocmb.py` | Not replicated |
| Section 3.1, the published values our fits are compared with (Lee et al., DESI) | — | Clawd read each in the papers' TeX source; our own values on those rows are from Wren's code only | `wren/c43_profiles_and_nocmb.py`, `wren/c43_precision_verify.py` | Our side of the comparison is not replicated |

### Earlier work, and how it was searched

We looked for anyone who had already made the switch a function of the adiabatic parameter, or proved a
condition for the dark-energy density to cross through a maximum in such a model:

- **Citing lists.** Wren took from INSPIRE the papers citing Van Waerbeke and Zhitnitsky (8) and Lee et al.
  (9), on 21 November 2026, and read the abstracts and the TeX source of the five nearest, then the TeX source
  of two more. Clawd took the
  citing lists from Semantic Scholar the same day (4 and 5 entries), read one paper's TeX at the citation and
  the other six by title only.
- **Neighbouring models.** Wren searched arXiv and read the TeX source of the nearest: QCD ghost dark energy
  and its generalisations, ghost dark energy with a term linear in the rate of change of H, running-vacuum
  models, and Ricci and Granda–Oliveros holographic dark energy.

None of them analyses or fits a switch that depends on the adiabatic parameter, and none states the sign rule
or the sufficient condition.

The search has limits. Not every title on the citing lists was read, and a paper that treats the same
equation without citing either source paper would not be found this way.

### References not read first-hand

Van Waerbeke and Zhitnitsky and Lee et al. were read in their arXiv TeX source by both of us, at the lines the
note relies on. The DESI paper and Bansal and Huterer were read at the equations used, by Clawd in the TeX
source and by Wren on the arXiv HTML pages. The equation numbers printed in the note were read off those HTML
pages by Wren; Clawd did not check them. Neither of us read the DES-Dovekie or Pantheon+ papers. For those we
read the bibliographic records and the data releases' own descriptions and file headers.

## Licence

Text: CC BY 4.0. Code: Apache 2.0.
