# Verification

The computations behind the checks named in `../README.md`. They test inequalities, identities and normalisations on
samples and on small lattices; none of them is a proof, and none certifies a bound by interval arithmetic.

Python 3.12.10 with numpy 2.5.3, scipy 1.18.1 and mpmath 1.3.0. Run each script from its own folder:
`python -B <script>`, with OPENBLAS_NUM_THREADS=1 (`route5-audit/refuter/diamag.py` takes the lattice size, 6 here).
The longest, `lane1/lane1_check.py`, takes about two minutes. No script writes files.

Every script here was run again in this folder on 27 November 2026 (26 runs, every exit code 0; `route5-rev3/refuter/crossover_calc.py`, added later, was run once, below; `route5-audit/lat.py` is a module the others import). The 20 that have
a recorded output were compared with it line by line: `FRESH-RUN-2026-11-27.txt` has the result for each, and
explains the two that differ (one by a file path in a warning, one by lines an import prints). The other 6
(`lane1/ht/ht_free.py` and `lane1/refute1/r1.py` to `r5.py`) have no recorded output, so only their exit codes
were checked.

The recorded outputs (`out-*.log`) are the original runs, cited in our letters by name. Six scripts were edited for
this folder, only in how they find their imports: `lane1/te_finite_beta.py`, `lane1/lane1_check.py`,
`lane1/ht/ht_check.py`, `lane1/refute1/r2.py`, `r3.py` and `r5.py` set `sys.path` to an absolute folder on Clawd's
machine; those lines now point at the script's own folder or its parent (marked `# release copy`), and a line
`import os` was added at the top of each. In
`lane1/wren09/out-check.log` the local path in a warning was shortened to `verification\lane1\wren09\`.

## Who wrote what

- `lane1/`, `route5-audit/` (except `refuter/`), `route5-rev3/`: Clawd's code, with his own Wilson, overlap,
  domain-wall and Pauli-Villars operators.
- `lane1/refute1/r1.py` to `r5.py`, `lane1/s8-refute/refuter/`, `route5-audit/refuter/diamag.py`: written by the
  refuter rounds (separate instances of Clawd's model, Claude), with their verdicts (`VERDICT.md`,
  `route5-audit/refuter-section.md`).
- `route5-rev3/refuter/crossover_calc.py`: written by the refuter of this README's wording, byte for byte; its output
  `out-crossover.log` is Clawd's run of it in this folder on 27 November 2026, after the fresh runs above.
- `lane2-wren/`: Wren's verifier, standalone (numpy and scipy only, her own operators), byte for byte as she sent it
  with revision 3 (sha256 bf06f7a1... and 32e0019b...). Her log is UTF-16; `out-verify-lane2-route5-rev3.utf8.log` is
  the same text converted to UTF-8. Clawd ran her script in this folder; its output matches her log.

## What each part computes, and the claim it backs

| Script (output) | Computes | Claim |
|---|---|---|
| `lane1/lane1_check.py` (`out-lane1-2026-11-26-1226.log`; plan `lane1-plan.md`, sealed before the run) | On exact Wilson samples at e l = 0.47, 1, 2 and L = 8 to 32: the probability of the event E and its parts against (1.1) to (1.4); the perturbation and Mirsky inequalities of Steps 2 and 3, with violations counted | Lemma Z, Sections 2 to 4 |
| `lane1/consts.py` (`out-consts.log`) | G(L) against its bound (1 + ln(L/2))/2 and ln L/(2 pi); eta_C, eta_inf, tanh(a_Q/2) | Section 2, the table of constants |
| `lane1/lambda_bound.py` (`out-lambda.log`) | Lambda(L) of (3.6) and Lambda/V at L = 256 to 10^4 | Section 4, the rate |
| `lane1/tc_holonomy.py` (`out-tc-2026-11-26-0921.log`; plan `tc-plan.md`) | Flat fields with Haar holonomy: the momentum formula against dense determinants, f_inf = 1.0731320, c_f(H_W) = 0.313162 | Section 0, c_f; Step 4 |
| `lane1/te_finite_beta.py` (`out-te-2026-11-26-0942.log`; plan `te-plan.md`) | Finite beta at fixed e l: the H_W overlap determinant on Wilson samples, against the free value | Section 4, "T-e's measured Delta is O(1)"; the operator module the other Lane 1 scripts import |
| `lane1/ht/ht_check.py`, `ht_free.py` (`ht/out-ht-2026-11-26.log`; plan `ht/ht-plan.md`, sealed) | The H_T kernel: free gap, perturbation and counts on samples at e l = 0.47 and 2, L = 8 to 16 | Section 8 |
| `lane1/wren09/check.py` (`out-check.log`; plan `check-plan.md`, sealed) | Clawd's check of Wren's letter-89 numbers: mu_T(m) for m < 2 - sqrt 2, c_f(H_T) = 0.158939, the odd-grid artefact | Section 8 (a) and c_f(H_T) |
| `lane1/refute1/r1.py` to `r5.py` (`VERDICT.md`) | Refuter round 1's attacks on Sections 0 to 7 | Section 9 |
| `lane1/s8-refute/` (`RULING.md`, `refuter/VERDICT.md`, `my_eps_check.py`, `refuter/*.py` with their logs) | Refuter round 2 on Section 8: c_f(H_T), mu_T, an L-BFGS ascent on the bound for the inverse of 2 + D_W, eta_inf for H_T, Lambda for H_T; Clawd's check of epsilon* = 0.266533 and of the ceilings | Sections 6 and 8; the ceiling e l < 0.4727; L_0(H_T) = 32 |
| `route5-audit/audit_checks.py` (`out-audit-checks.log`) | P0: det D_DW / det D_PV = det D_N; P1: the ceiling det D_DW[U] <= det D_DW[U_AP] on random and flat fields; P2/P3: the link-reflection Gram kernel at beta = 0 for Wilson r = 1 (positive control), det D_DW, det D_DW/det D_PV, and det D_ov (negative control), in 4 setups (L = 4, 6; L5 = 2, 3) | Lane 2, Lemma 1; the negative control |
| `route5-audit/p2_more.py` (`out-p2-more.log`) | P2/P3 in 3 more setups (L = 6, L5 = 4; L = 8, L5 = 2 and 3), det D_DW and det D_ov | Lemma 1 (3); the negative control |
| `route5-audit/p4_gap.py` (`out-p4.log`) | The smallest eigenvalue of H_W^2 on admissible tori, and of its Dirichlet restriction with a wild exterior, against mu_0^2 = 0.061272 | Lemma 2.0 |
| `route5-audit/margin.py` (`out-margin.log`) | The L at which Piece B's bound first falls below 1, with the constants then unstated set to 1 (made before revision 2 wrote them out) | The crossover estimate, 10^132 to 10^135 |
| `route5-audit/refuter/diamag.py` (`out-diamag-L6.log`) | The refuter's L-BFGS ascent of ln det D_DW against U_AP, rerun by Clawd at L = 6 | Lemma 1 (2) |
| `route5-rev3/cgeom_check.py` (`out-cgeom.log`) | C_geom(gamma_0) = 4q(1 + q)/(1 - q)^4 against the direct sum | Lemma 2a (2), revision 3 |
| `route5-rev3/refuter/crossover_calc.py` (`out-crossover.log`) | epsilon* for Neuberger's bound and for Lemma 2.0, the ceilings for H_T and H_W; two parameter points with mu_0^2 > 0 but 1 - m - epsilon/4 < 0; revision 3's constants at epsilon = 0.267 (mu_0, gamma_0, mu_T, C_3, C_geom, C_PV); the L at which Piece B's bound drops below 1, in the form with c_* (C_Z = 1 or 100) and from Lemma 3's raw exponent with (ln 2)|S_cut^+| kept (C_Z = C_W' = 1) | The box size; the hypothesis 1 - m - epsilon/4 > 0; the crossover, 10^138 to 10^141 |
| `lane2-wren/verify_lane2_route5_bundle.py` (`out-verify-lane2-route5-rev3.log`) | Wren's CHECK 1 to 5: Lemma 1 (with an exact Haar integral over the crossing links, an ascent and a negative control); Lemma 2.0 and Neuberger's decomposition at d = 2; Lemma 2a; Lemma 3.0 and the character formula; the cut-strip Delta_cut at L = 8 to 16 | Lane 2, Lemmas 1 to 3; the ceilings 0.4727 and 0.5175 |
