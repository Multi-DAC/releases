# VERDICT: PARTIALLY REFUTED

Start 23:49:56, stop 23:57 PST (2026-11-26). Scripts and logs are in this folder.

The inequality ln E_W[det D^{H_T,(m_f)}] >= -c_f(H_T) V - O(L (ln L)^{3/2}) survived every attack I ran. What dies is a
stated Section 8 constant that the downstream proof uses: e l_0 = 0.4735, and the line-286 claim that Adams' eps = 0.267
gaps H_W and H_T at m = 0.7.

| # | section/line | attack | method | result | ruling |
|---|---|---|---|---|---|
| 1 | 8, l.284-286 and 6, l.248-250 (eps = 0.267 "also gaps H_T"; e l_0 = 0.4735) | Is Adams' bound positive at the eps used? | algebra plus cf_check.py, out-cf.log | Using the artefact's own quoted bound, sqrt(1 - c_1 eps) - abs(1 - m) with c_1 = 2 + sqrt2 at d = 2, the value at m = 0.7, eps = 0.267 is **-0.00267 < 0**. A positive gap needs eps < eps* = m(2-m)/c_1 = **0.266533**. So "|H_T| >= abs(H_W)/(6-m) > 0 on abs(1 - U_p) < 0.267" is not established, and e l_0 = eps/sqrt(2 c_f) must be **< 0.47274** (strict), not 0.4735. The note "rounded down" is false: 0.267 is 0.26653 rounded UP past the zero of the bound. The same defect gives H_W's 0.337, which should be < 0.33678. | **KILL** (of the constant; the order is unchanged) |
| 2 | 8, l.279 (c_f(H_T) = 0.158939) | Recompute it independently | cf_check.py: midpoint N = 512 to 4096 plus Richardson | 0.1589389597. H_W gives 0.3131624, which also matches. | SURVIVES |
| 3 | Use: c_f equals the free density at U_AP | Is ln det D_{m_f}[U_AP] + c_f V within O(L)? | free_check.py (exact momentum sum, AP in time, m_f = 1/L) | 0.61, 0.71, 0.80, 0.89, 0.99, 1.08 at L = 8 to 256: O(ln L) from m_f^2 sum 1/abs(lambda0)^2. Massless it is 0.3466, which is O(1). The sign is unfavourable (the denominator is larger) but it is only O(ln L). | SURVIVES |
| 4 | 8(a) mu_T(m), every holonomy, 0 < m <= 1 | Free min abs(H_T) against the formula | cf_check.py, grid 2001^2, m = 0.1 to 1.0 | Holds at every m tested. It is exact (m/(2-m)) for m <= 0.586. At m = 0.7: 0.3939 against 0.37928. | SURVIVES |
| 5 | 8(b), norm(inverse of (2 + D_W(-m))) <= 1/(2-m) for every U | Adversarial ascent over the link angles | op_check.py, L = 4, L-BFGS, 4 starts | The maximum of (2-m) times the norm is 0.99999999: the bound is saturated but never exceeded. | SURVIVES |
| 6 | 8(b), the gap of H_T1 >= mu_T/2 under the sup cap eta_T | Minimise the gap over a (abs(a) <= eta_T) and the holonomy phi | op_check.py, L-BFGS, 4 starts | The minimum reached was 0.366, against the bound 0.190. | SURVIVES |
| 7 | 8(b) op/HS constants 4/(2-m)^2 and 2sqrt2/(2-m)^2 | Maximise the ratios | op_check.py | The maxima reached were 0.745 (op, bound 2.367) and 0.502 (HS, bound 1.674). Both are loose and valid. h_T = 4.41265, which matches. | SURVIVES |
| 8 | 8(c) n_0 and T_0 bounds with c = 1/(4pi^2) | Exact free counts, worst of 144 holonomies | free_check.py, L = 32, 64, 128 | Every case satisfies the bound, by a factor of 10 to 100. The algebra (b_T < 0 implies n < 4 and w_T <= 1) checks. | SURVIVES |
| 9 | 2, l.146 ("P_E >= 2/3 - 2e-6 for L >= 16 at e l <= 0.47"), used under eta_T | Is it re-derived for the H_T cap? | etainf_T.py | No. With eta_T = 0.0801, eta_inf,T at e l = 0.47 is 14.7 (vacuous) at L = 16, 0.50 at L = 24, and 0.004 at L = 32. L_0 rises, and Section 8 gives no numbers for it. | WOUND (minor) |
| 10 | Lambda for H_T ("same form, bigger constants") | Evaluate it | lambda_T.py, out-lambdaT.log | At e l = 0.47, Lambda/V is 0.51 at L = 1024, 0.164 at L = 4096 and 0.078 at L = 10^4. It stays above c_f(H_T) = 0.159 until L is about 4500 (about 500 for H_W). The order O(L (ln L)^{3/2}) holds. | SURVIVES (order) |
| 11 | Adams' d-dimensional generalisation (hypotheses, normalisation of H(m)) | Are its hypotheses met at d = 2, r = 1? | not read first-hand (the artefact says so too) | Not checked. | CANNOT-DETERMINE |

**Overall: PARTIALLY REFUTED.** The lemma's inequality and its use of c_f survive. The Section 8 ceiling e l_0 = 0.4735,
and the gap claim at eps = 0.267 behind it, are false as stated.

**Repair (a lead; not verified).** The root cause is rounding a threshold in the unsafe direction: 0.26653 was written as
0.267. Replace eps with any eps < eps*(m) = m(2-m)/(2+sqrt2), for example eps = 0.2665. The H_T gap on the admissible set is
then (sqrt(1 - c_1 eps) - abs(1 - m))/(6 - m) > 0, a checkable positive number (about 3.5e-5 at eps = 0.2665, which is tiny).
The ceiling becomes e l_0 = eps/sqrt(2 c_f(H_T)) < 0.47274 (0.4726 at eps = 0.2665). Any downstream use that needs a
non-negligible gap on G has to pick a smaller eps and live with a lower ceiling. This adds a checkable condition; it does
not merely protect the claim. Separately, Section 8 should restate L_0 under eta_T (L >= 32 rather than 16 at e l = 0.47).
