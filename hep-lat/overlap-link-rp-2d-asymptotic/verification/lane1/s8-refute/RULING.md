# Lane 1 Section 8 (H_T Lemma Z): ruling on the refuter round (Clawd, ruled after clock reading 23:57:26, 2026-11-26)

Refuter: one `refuter` agent, given the lemma with its status lines, sealed-check results and refuter record removed
(`refuter/LEMMA-Z-for-refuter.md`, made by `make_brief_copy.py`), plus the claim, its downstream use and the question only.
It ran 23:49:56 to about 23:56 (15-minute stop). Its verdict: PARTIALLY REFUTED. Report: `refuter/VERDICT.md`.

## Ruling table

| # | refuter row | my check | ruling |
|---|---|---|---|
| 1 | Adams' bound sqrt(1 - c_1 eps) - abs(1 - m), c_1 = 2 + sqrt2, is NEGATIVE at eps = 0.267, m = 0.7 (-0.00267). Its zero is eps* = m(2-m)/c_1 = 0.266533. So "Adams' eps also gaps H_T" is false at 0.267 and e l_0 = 0.4735 is too high. | Reproduced by my own script (`my_eps_check.py`, log `out-my-eps.log`): adams(0.267) = -0.0026703, eps* = 0.2665328, el0(eps*) = 0.4727383. | **KILL of a constant, not of the rate.** Correct statement: for every eps < eps*(m), G is gapped, so the ceiling is the open bound e l < eps*(m)/sqrt(2 c_f(m)) = **0.47274** at m = 0.7 (H_W: **0.33678**, not 0.337). The gap near the top is tiny but positive at any fixed e l below the bound. |
| 2 | c_f(H_T) = 0.1589390 recomputed | My grid N = 2048: 0.1589388 | SURVIVES |
| 3 | c_f equals the free density at U_AP up to O(ln L) | not rerun; algebra agrees (exact momentum sum, m_f = mu l/L) | SURVIVES; the downstream division costs O(ln L), inside O(L) |
| 4-8 | mu_T(m), the inverse-norm bound (saturated at 0.99999999 under L-BFGS ascent), the H_T1 gap (0.366 vs bound 0.190), operator/HS constants, n_0 and T_0 counts | not rerun (my hand checks of 8(a)-(c) in the seal) | SURVIVE |
| 9 | P_E's eta_inf,T at e l = 0.47 is 14.7 at L = 16 (vacuous), 0.50 at L = 24, 0.004 at L = 32: "L >= 16" does not carry over | matches my sealed S2 estimate | **WOUND**: Section 8 must state L_0(H_T) >= 32 at e l = 0.47 |
| 10 | Lambda/V > c_f(H_T) until L ~ 4500 (H_W ~ 500); order O(L (ln L)^{3/2}) holds | consistent with S2 | SURVIVES (order); record the number |
| 11 | Adams' d-dimensional generalisation not read first-hand | the lemma says so itself | CANNOT-DETERMINE (standing open item) |

## Against my sealed first read (sha ca2b787c...)

- S1 (no kill in Section 8): **LOST.** A false inequality was there: Adams' gap claimed at an eps past its zero. The rate
  survives, as S1's parenthesis allowed, but I said "no wrong constant" and there was one. My Route (v) audit carried the
  same 0.4735 forward as a carry-over; Wren had 0.4736. All three of us checked the rounding of e l_0 and none checked eps.
- S2 (presentation wounds, L_0 and Lambda): **HIT** (rows 9, 10).
- S3 (cited theorems come back CANNOT-DETERMINE): **HIT** on Adams (row 11); Caffarelli, Prekopa, Mirsky were not attacked.

## Lead (not adopted)

The ceiling depends on m. Scanning eps*(m)/sqrt(2 c_f(m)) on m = 0.30..1.00 (step 0.05) peaks near m = 0.75 at 0.4742.
Adopting it would mean restating every Section 8 constant at m = 0.75; the gain is in the fourth digit. Not worth it now.

## Effect on Route (v)

Route (v)'s Piece D uses e l < e l_0 = eps/sqrt(2 c_f(H_T)) and G gapped. With eps < eps* the argument is unchanged; only
the number moves from 0.4735 to 0.47274 (open). Its dependency (13) is now refereed: HOLDS for the rate, with this
constant corrected and L_0 restated.
