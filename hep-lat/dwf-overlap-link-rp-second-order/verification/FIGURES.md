# Every computed figure in the paper, and the output line that prints it

All scripts are in this folder; see `README.md` for the environment and for how each `.out` was made. An output
line is named by a string to grep for in the `.out` file. "Wren's code" means a figure computed by Wren with her own
code. Six of them are recomputed by her standalone `wren_release10_checks.py` (her letter 81), whose output
`wren_release10_checks.out` is in this folder; the CHECK block that prints each is named. Clawd has not run this
code.

## Section 1 and Section 3 (Part (a))

| Paper | Figure | Source, output line |
|---|---|---|
| s1 (a); s3.2 "Size of the constant"; s6 | tr N~_PV linear in S, between 0.046 S and 0.076 S for (m, p1) in {0.7, 1.0} x {0, pi} and S in {4, 8, 16}, finite a5 -> 0 limit at fixed S | Wren's code (first in her letter 72, item B): `wren_release10_checks.out`, CHECK 1, "Range of tr N~_PV(2x2)/S across S in {4,8,16}: [0.04567, 0.07642]" and the "Fixed-S a5 -> 0 convergence" rows |
| s3.2 "Size of the constant" | with the other normalisation, 39.7 against 0.248 at L5 = 8, and 1346 against 0.255 at L5 = 64, S = 4 | Wren's code (first in her letter 73, Defect 5): `wren_release10_checks.out`, CHECK 2, "L5= 8 a5=0.5000 (S=4)": sum W~_j=0.24839, sum W~_j/mu_j^2=39.68; "L5=64 a5=0.0625 (S=4)": 0.25544 and 1345.84 |
| s3.2 Remark on the range of a5 | at a5 = 2/m, E_PV(0, pi) = 0 and rho_PV = 1 | `check_letter74_breaks.out`, "v2-1: a5 = 2/m, PV": bulk sigma_min 0.000e+00 at L5 = 1 and 3 (Wren's finding, letter 74; this is Clawd's check of it) |
| s3.2 Remark on the range of a5; s6 | at a5 > 1/m, L5 = 1, D_DW has an open-boundary zero mode; present for L5 = 1 and 3, absent for L5 = 2 and 4 at a5 = 1.45, 1.6 and 2.0 | `check_letter75_vmp.out`, "(b) min over k1 of sigma_min(V_{-,+})": rows "DW L5=1" and "DW L5=3" fall to 1e-14 or below from a5 = 1.44; rows "DW L5=2" and "DW L5=4" stay at 0.19 or above. `check_letter75_h4.out`, block "B. threshold 1/m": the plane-propagator error by 2L, including "L5=4 k1=0.0 DW a5=2.00" |
| s3.2 Example | rho_*(0.7) = 0.189036 | `pv_gap_check.out`, "m=0.7000: grid min ... rho_*(m) = exp(-2 Emin) = 0.189036" |
| s3.2 Example | a5_max(0.7) = 0.659341 | by hand: 2(0.3)/(0.7 x 1.3) = 0.659341 |
| s3.2 Remark (m >= 1); s6 | at a5 = 1 and m > 1, rho_PV(m, 1, p1) can exceed rho_*(m, p1) at intermediate p1 while staying <= rho_*(m) | Wren's code (first in her letter 73): `wren_release10_checks.out`, CHECK 3, "m=1.3 a5=1": max_p1 rho_PV=0.346021 = rho_*(m), max_p1 [rho_PV(p1) - rho_*(p1)]=+9.1305e-04 at p1=0.800*pi (controls m = 0.7, 1.0: +0.0000e+00) |
| s3.3 Lemma 1, numerical check | mirror identity exactly 0 for DW and PV; K Hermitian to 1e-14; min/max eigenvalue +3.1e-6 (DW), +9.6e-8 (PV) | `check_letter77.out`, "(C1) |D_-[theta U] ... = 0.0e+00" (mf=0 and mf=1); "mf=0 flat crossing: Hermitian dev 1.0e-14; eigenvalues min/max = +3.149e-06"; "mf=1 flat crossing: Hermitian dev 1.7e-14; ... +9.554e-08" |
| s3.3 Lemma 1, controls | crossing phase (2.0, 0.0): -4.7e-3 (DW), -3.5e-3 (PV); no conjugation on U0: -3.4e-5, -1.9e-6 | `check_letter77.out`, "control, crossing (2.0, 0.0)": -4.675e-03 and -3.527e-03; "control, no conj on U0": -3.440e-05 and -1.929e-06 |
| s3.3 Lemma 2, (H2) | d_y log det D_DW = +0.3128 at twist pi/3 (open time, L = 6, Ls = 4, m = 0.7, L5 = 6, a5 = 0.25) | `check_fig4_twist.out` (Clawd's code), "time_bc open twist pi/3 DW: d_y log det = +3.1277e-01", with twist 0 and pi at most 2.3e-08 (finite-difference noise). Wren's CHECK 4 in `wren_release10_checks.out`, block "time_bc = open", gives +3.1277e-01 too. Her letter 73 figure, +0.3156, was computed with periodic time (her letter 81); both codes give +0.3109 with antiperiodic time |
| s3.3 Lemma 2, proof of (H2) | gamma5 = -i gamma0 gamma1; C = sigma2 relations (exact) | `check_letter77.out`, "(A) gamma algebra": every line 0.0e+00 |
| s3.3 (H4) remark | (H4) fails above 1/m for L5 = 1 and 3, not for 2 and 4 | as for the remark on the range of a5, above |
| s3.3 (H4) remark | at a5 = 2, L5 = 1 the bulk is gapped, yet an exponentially small open-boundary mode and an O(1) error at the plane | `check_letter74_breaks.out`, "v2-2: a5 = 2, L5 = 1": "bulk sigma_min (p0 = 0) 2.000e-01"; sigma_min 2.519e-02 ... 1.790e-10 at 2L = 12 ... 96 and "max|G_open - G_inf| 1.281e+00" throughout; the control at a5 = 1 has 1.277e-15. (Wren's counterexample, letter 74; this is Clawd's check of it) |
| s3.3 Lemma 4, numerical check | closed form against a direct root to 7.66e-15 over 400 random points; r = 2 control agrees at 0 of 200 | `pv_gap_check.out`, "max |formula - root| = 7.66e-15"; "r=2 control: 0/200 random points agree" |
| s3.3 Lemma 4, numerical check | on 721 x 721 grids, min over k5 of E_PV minus E_cut(k1) exactly 0 where (i)'s condition holds, negative where it fails; m = 0.7, a5 = 1: -0.371 at k1 = 0 | `pv_gap_check.out`, block "(ii)": "+0.000e+00 ... True" rows and negative "... False" rows; "m=0.7 a5=1.0  min(E_PV - E_cut) = -3.711e-01" |
| s3.3 Lemma 4, numerical check | min over k1 of E_cut matches the formula to 10 digits at m = 0.1 ... 0.9 | `pv_gap_check.out`, block "(iii)": "grid min" equals "predicted" to the 10 printed digits on each row |
| s3.3 Lemma 4, numerical check | on a grid m = 0.05 ... 0.95, min over (k1, k5) of E_PV minus min_q E_cut(q) >= -1.1e-15 (rounding) for every a5 <= 1, and -4.8e-2 at a5 = 1.05, m = 0.95 | `lemma4_grid.out`, "min over m, a5 <= 1, (k1, k5) of E_PV - min_q E_cut: -1.05e-15" and "a5 = 1.05  m = 0.95  min diff = -4.757e-02". The script was written for this release; the original grid (25 November, 09:49) was not saved, and its bound read -1.0e-15 |

## Section 4 (Part (b))

| Paper | Figure | Source, output line |
|---|---|---|
| Cor 3.1(2) | the two sums are at most 0.354 and 0.173 (120 random points) | Wren's code (first in her letter 77, `attack_merged_t12.py`): `wren_release10_checks.out`, CHECK 5, "max sum_j |A1(1,j)|^2 = 0.353793" and "max sum_j |A0(1,j)|^2 = 0.172474" |
| Remark after Prop 4 | 240 grid cases, t up to 10^4, L5 16 to 128; maximum 1.078 on the grid | `prop4_large_t.out`, "WORST (np.float64(1.0780713023503659), (0.3, 1.0, 0.0, 16, 0))"; `prop4_parse.out`, "rows: 240" and "worst sup_t: 1.0781". (In the original run 8 of the 240 cases stopped on the pairing assertion and were run by `prop4_rerun.py`; the module `l4_assembled.py` here has the clustering fix, so the fresh `prop4_large_t.out` runs all 240, and `prop4_rerun.out` repeats those 8.) |
| Remark after Prop 4 | 2.044 beyond the grid, at m = 0.3, a5 = 0.999 a5max = 2.742 | `prop4_a5hi.out`, "m=0.3 a5=2.7424 k1=0.000 L5= 16 S=  43.9: max 2.0442" (the same at L5 = 32 and 64) |
| Remark after Prop 4 | the slowest nine to t = 10^7; at the exact tie 0.359, flat to t = 10^7 | `prop4_tail.out`: the nine rows with t to 10000000; "EXACT TIE m=0.7 a5=a5max L5=16: ... r = 0.359 0.359 ..." |
| Remark after Prop 4 | the supremum sits at t = 0 in 203 of the 240 cases; saturation in L5 | `prop4_parse.out`, "cases with argmax t > 0: 37 of 240" and "sup_t by L5 (saturation)" |
| Remark after Prop 4 | the supremum divided by a5 depends on S alone and saturates by S of about 13 | `prop4_large_t.out` and `prop4_a5hi.out` rows (columns a5, S, max); read off the table, not printed as one line |
| Remark after Prop 4 | the cancellation is collective: pole groups paired from the top have trace norms growing like L5 while the total stays fixed | `prop4_pairs.out` |
| Remark after Prop 4 | at k1 = pi, entries fall like t^-1/2 while the top singular vector spreads like sqrt(t); the trace norm stays flat | `letter78_tests.recorded.out`, block "T1  m=0.3 a5=1.0 k1=pi L5=128": "trace norm 0.0866 ... 0.0844", "sqrt(t)*max", "<wall dist>/sqrt(t)" |
| Lemma L4, numerical check | exact pole sums match position space to 1.4e-16 | `l4_assembled-validate.out`, "twist=0.00 mf=0: max |C_mine - C_pos| = 1.39e-16" |
| Lemma L4, numerical check | max ‖M_in(t)‖_op / rho_*^t = 0.120 over S = 3.2 to 51.2; / ((t + 2)^4 rho_*^t) = 7.5e-3, maximum at t = 0 | `l4_sdep.out`, "m=0.3 twist=3.14 L5= 64 a5=0.20 S= 12.8: max_t ||M_in||/rho^t = 0.1197"; 7.5e-3 = 0.1197 / 16, since the maximum is at t = 0 (the per-t rows of `l4_sdep.out`); `l4_assembled-main.out` for S = 3.2 and 6.4 ("/((t+2)^4 rho^t)") |
| Lemma L4, numerical check | at fixed S the value converges as a5 halves; at fixed a5 it rises with S and levels off by S of about 13 | `l4_assembled-main.out`, "M" rows (S = 3.2, a5 = 0.2, 0.1, 0.05); `l4_sdep.out`, rows S = 3.2, 6.4, 12.8, 25.6, 51.2 |
| Lemma L4, numerical check | a5 <= a5max is sharp at k1 = 0: just above it the whole DW bulk band at k1 = 0 crosses r_*(0); at theta_s = pi it is not necessary | `l4_sdep.out`, "C(ii) m=0.9 (a5max 0.202)" rows: "DW poles above r_* per k = [...]", L5 poles at k1 = 0 for a5 = 0.5, 1.0, 0.25 |

## Section 5 (Part (c))

| Paper | Figure | Source, output line |
|---|---|---|
| Thm T1-2, the domain, numerical check | lambda_0^ov <= e^{-c0 |k1|} in all 7 (m, k1) rows with a pole, at m = 0.3 and 0.7 | `letter78_tests.recorded.out`, block "T3": "overlap pole ... (r_* ...; exp(-c0|k1|) ...)" on the 7 rows with a pole |
| Thm T1-2, numerical check | the DW surface pole converges to lambda_0^ov exponentially in S; at m = 0.3, k1 = pi/16, a5 = 0.2 the gap is 5.1e-3, 3.6e-5, 1.8e-9, 4.8e-15 at S = 8, 16, 32, 64 | `letter78_t3c.recorded.out`, "k1=0.196 a5=0.20 overlap pole 0.823762": gaps 5.11e-03, 3.56e-05, 1.81e-09, 4.77e-15 |
| Order of limits, numerical check | the antiperiodic-time and IR entries agree at p1 = pi/2 and pi, the difference halving with a5; at p1 = 0 they differ | Wren's code (first in her letter 77, `verify_p5_and_catches.py`): `wren_release10_checks.out`, CHECK 6, the "theta_s=0 Ls=4 m=0.7" rows: max|diff| 1.092e-04 -> 5.225e-05 (p1 = pi/2), 1.581e-04 -> 8.436e-05 (p1 = pi), 1.253e-01 -> 1.251e-01 (p1 = 0). NOT checked with Clawd's code |
| Thm T1-2, p1 = 0 example | the k1 = 0 mode at large tau: +1.0 (IR) against +0.5 (AP) in the (1,1) entry, 0.0 (IR) against +0.5 (AP) in the (0,0) entry | Wren's statement (letter 78). Clawd's check: `oolimits00_test.recorded.out`, "AP, exact overlap ... Q_11 +0.500001 ... Q_00 +0.500000"; IR rows "Q_00 +0.000000"; "IR extrapolated to eps = 0 (quadratic in eps), Q_11: +0.993006"; the fit A e^{-2 eps tau} (A = 0.99970, tau = 8.98) is in Clawd's note `work/drives/2026-11-25/oolimits00-result.md`, not printed by the script; twist-pi control rows ~0 |
| Prop 5, D_eff | at m = 0.7, S = 8, |det| differences 5.0e-2, 2.6e-2, 1.35e-2, 6.8e-3 at a5 = 0.2, 0.1, 0.05, 0.025; S = 4 nearly the same | `check_prop5_deff.out`, "m=0.7 S=8": 5.03e-02, 2.63e-02, 1.35e-02, 6.81e-03; "m=0.7 S=4": 4.92e-02 ... 6.50e-03 |
| Prop 5, D_eff | Sylvester's identity to 2e-11; min w = m | `check_prop5_deff.out`, largest "Sylvester 2e-11" (m=0.3 S=8 a5=0.025); "min over torus of w: m=0.3 0.300000, m=0.7 0.700000" |
| Prop 5, numerical check (Steps 1-2) | T5 Hermitian to 1.1e-16 and positive definite; Mercator identity to 7.8e-16; det D_N = det D_DW / det D_PV to 5.2e-14 | `check_letter77.out`, block "(D)": largest "T5 Hermitian 1.1e-16", "min eig" all positive, largest "T5 - (I - a5 H_W + a5^2 R) 7.8e-16", largest "|det D_N / (det DW / det PV) - 1| = 5.2e-14" |
| Prop 5, numerical check | the Hessian of log det D_N equals the 5d one to 7e-6 (L5 = 6, a5 = 0.5) | `check_letter77.out`, block "(E2)": "max relative difference 7.0e-06" |
| Prop 5, numerical check (joint limit) | max |H_ov| = 3.73e-2; |H_N - H_S| = 8.36, 4.31, 2.19, 1.11 e-3 at S = 8, the same within 0.5 % at S = 16; |H_S - H_ov| = 1.65e-2, 8.05e-4, 1.77e-6 at S = 4, 8, 16 | `check_letter77.out`, block "(E)": "|H_ov| max 3.7343e-02"; rows "S=   8", "S=  16" (8.32, 4.30, 2.19, 1.11 e-3), "S=   4" |
| Prop 5, numerical check | about e^{-0.76 S}; the a5 constant is 22 to 27 % of |H_ov| at a5 = 0.2 | arithmetic on the block (E) lines: ln(1.65e-2 / 8.05e-4) / 4 = 0.755 and ln(8.05e-4 / 1.77e-6) / 8 = 0.765; 8.32e-3 / 3.7343e-2 = 22.3 % and 1.02e-2 / 3.7343e-2 = 27.3 % |
| Prop 5, k1 = 0, the pole z0 | 1 - z0 = 7.4901e-6 against 2m(1 - a5 m)^L5 = 8.0539e-6, ratio 0.9300; at a5 = 0.1 ratio 0.9650; 2m e^{-mS} = 1.91e-5 | `check_letter77.out`, block "(B)": "a5=0.20 L5=80: 1 - z0 = 7.4901e-06; 2m(1 - a5 m)^L5 = 8.0539e-06 (ratio 0.9300 ...); 2m e^-mS = 1.9144e-05" and "a5=0.10 L5=160: ... (ratio 0.9650 ...)" |
| Prop 5, k1 = 0, the residue | D_N(1, 0) = q/(1+q) I | No figure: the text now gives it exactly (at z = 1, T5 = diag((1 - a5 m)^-1, 1 - a5 m)) |
| Prop 5, k1 = 0, the residue | one zero in (0.9, 1), in the + channel only; agrees with the 5d pencil's z0 to 8 digits | `residue_z0_test.out`, "zeros in (0.9, 1): + ['0.9959409881'], - []" (m = 0.7, a5 = 0.1, S = 8) against `letter78_tests.recorded.out` T2, "m=0.7 a5=0.10 S= 8.0 L5= 80: DW inside poles >= r_mid: 1 ['0.99594099']"; likewise 0.9954096290 against 0.99540963 at a5 = 0.05 |
| Prop 5, k1 = 0, the residue | |c/m - 1| = 3.6, 1.9, 1.1 e-2 at a5 = 0.2, 0.1, 0.05 (m = 0.7, S = 8); 1.75e-2 at a5 = 0.1, S = 16; 3.8e-3 at m = 0.3, a5 = 0.05, S = 32 | `residue_z0_test.out`, "/(2m P0^+) dev" on the rows: 3.64e-02, 1.94e-02, 1.10e-02; 1.75e-02; 3.77e-03 |
| Prop 5 (a) | c1(0.1) = 7.63e-3, c1(0.3) = 7.15e-2; equality at k1 = 0, z = r_mid | `rew2_min.out`, "m=0.100 ... = 7.6316e-03" and "m=0.300 ... = 7.1471e-02" (the minimum of Re det X on |z| = r_mid equals the closed form c1(m)); `letter78_T4.out`, "c1 0.07147" |
| Prop 5, numerical check | min Re w^2 / c1 = 1.0000 at m = 0.05 and 0.3, larger above; max |b| = 3.95 | `letter78_T4.out`, "m=0.05: ... min Re(w^2)/c1 = 1.0000 ; max |b| = 3.950", "m=0.30: ... 1.0000", "m=0.70: ... 2.4413" and above |
| Prop 5, numerical check | at k1 = 0, a5 in {0.1, 0.05}: exactly one DW pole at or above r_mid once S >= S0(m); S0(0.3) in (4, 8], S0(0.7) in (2, 4], S0(0.9) <= 2 | `letter78_tests.recorded.out`, block "T2": "DW inside poles >= r_mid: 0" or ": 1" by (m, a5, S) |
| Prop 5, numerical check | the pole tends to 1: 0.9485 and 0.9955 at m = 0.3, S = 8, 16 | `letter78_tests.recorded.out`, T2, "m=0.3 a5=0.10 S= 8.0 ... ['0.94848937']" and "S=16.0 ... ['0.99548934']" |
| Prop 5, numerical check | every PV pole below 0.771 r_mid; nothing on |lambda| = 1 | `letter78_tests.recorded.out`, T2: largest "max PV/r_mid 0.7708"; "unit groups 0" on every row |

## Figures computed by hand, not by a script

- a5_max(0.7) = 0.659341, a5max(m) < 1 exactly when m > 2 - sqrt 2 = 0.586, r_mid >= 1/sqrt 2, |b| < 5,
  |B_++^-1| <= 4/3: closed forms in the text.
- C_s = 8, 2 C_s = 16 (Corollary 3.1(2)): bounds proved in the text.
- Hilbert's inequality constant pi: cited.
