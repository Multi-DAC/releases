# Every computed figure in the note, and the output line that prints it

Two codes stand behind the numbers. Clawd's is in this folder. Wren's, which produced Table 1, Table 2 and the
figures of Section 3 as printed, is in `wren/` (see the note at the end). The two were written independently
from the public data files and compared only after both were complete.

Clawd's scripts are Python 3 with numpy and scipy. `c43_fit.py` is the model and the likelihood; `sec26_table.py`,
`sec26_exclusion.py`, `one_maximum.py`, `c43_controls.py` and `c43_scan.py` import it from this folder.
`c43_props.py` is stand-alone and needs no data. The data files are not shipped: `DATA.md` says where they come
from. `rerun.sh` runs everything from scratch; `rerun.log` is the log of the run of 2026-11-21 (Python 3.12.10,
numpy 2.5.3, scipy 1.18.1), and each `X.out` is the output of that run.

## Sections 1 and 2 (Clawd's scripts)

| Note | Figure | Script, output line |
|---|---|---|
| s2.4 | fitted quadrature switch, c = 1.641: the only maximum of u on 0 <= z <= 1090 is at z = 0.47 | `one_maximum.out`: "quad: c 1.641 q 0.500 ... extrema of u on the whole grid (z, kind, u): [(0.473, 'max', 1.2827)]", both hand-overs |
| s2.6 | mu_m0 = 0.46942, mu_r0 = 1.3827e-4 | `sec26_table.out`, first line |
| s2.6, table, row "integrated" | 1 + z_max = 4.03, 7.29, 14.37, 36.17, 47.60, 68.16, 73.67, 82.66, 119.72 for c = 0.003 ... 10 | `sec26_table.out`: "c=... early maximum 1+z_max ..." (one line per c and hand-over; the two hand-overs agree to the digits printed) |
| s2.6, table, row "(2.9)" | 4.57, 7.39, 14.08, 35.36, 46.65, 66.97, 72.40, 81.23, 117.19 | the same lines: "(2.9) gives ..." |
| s2.6 | (2.9) within 2.3% of the integrated value for 0.05 <= c <= 10; off by 13% at c = 0.003 | `sec26_table.out`, summary: "largest \|(2.9)/integrated - 1\| for 0.05 <= c <= 10: 0.0224 (at c = 0.5)"; "(2.9)/integrated - 1 at c = 0.003: +0.1333" |
| s2.6 | dark energy between 0.8e-5 and 6.5e-5 of the total density at the maximum | `sec26_table.out`, summary: "7.83e-06 to 6.47e-05" |
| s2.6 | at c = 0.002 no extremum at all | `sec26_table.out`: "c=0.002 ... extrema of u on 0<=z<=1090: none" |
| s2.6 | N_e at 1 + z = 18.0, 46.3, 101.4 for c = 3.1, 4, 10 | `sec26_exclusion.out`: "E = c-3 at 1+z = 18.03", "46.31", "101.37" |
| s2.6 | the early maximum above N_e, at 1 + z = 74.6, 82.7, 119.7 | `sec26_exclusion.out`: "the first at 1+z = 74.64", "82.66", "119.72"; and w > -1 at every grid point below the bound (the two counts on each line are equal) |
| s2.6 | N_e -> infinity as c -> 3 | `sec26_exclusion.out`: c = 3.01 and 3.001 give 1 + z = 6.56 and 1.80 |
| s2.6, s2.8 | for c >= 0.5 the only maximum is at z >= 35; at lower redshift at most one minimum | `sec26_table.out`: the full list of extrema on each line with c >= 0.5 (lowest maximum: c = 0.5, 1 + z = 36.17) |
| s2.7, Prop. 1 | (u - 1)/mu at N = 3: +0.99959, +0.24968, +0.05249, +0.02028, -0.01001, -0.09098, -0.20003, -0.49998, against (1 - 3k)/(1 + 3k) | `c43_props.out`, block T1, eight lines |
| s2.7, Prop. 1 | linear switch, c = 3: (u - 1)/mu^2 = 0.08333 at mu = 1e-4, against 1/12 | `c43_props.out`, block T1b: "0.083327" from three starting points |
| s2.7, Lemma 2 | u at N = -4: 0.66645, 0.33355, 1.99147, 0.66623; x/h = 1.5000 | `c43_props.out`, block T2 |
| s2.7, Prop. 2 | five rational switches inside the hypotheses: a maximum above 1 in each; k = 0.34, 0.5: a minimum and no maximum; c_t = 3: no extremum | `c43_props.out`, block T3: the lines "inside", "mutant k", and "mutant ct (rational) k=0.20 ct=3.00" |
| s2.7, Prop. 3 | c = 0.3, 0.5, 1, 2, 2.9: one minimum, beta there 0.0973, 0.1567, 0.2881, 0.5800, 0.9508; c = 3, 3.1, 4, 6, 20: no extremum; beta rises in time in all ten | `c43_props.out`, block T4 |
| s2.7, Prop. 4 | literal form, five switches: x falls in time, u has no extremum | `c43_props.out`, block T5 |
| s2.8 | exactly one maximum of u: the five rational cases (matter only) | `c43_props.out`, block T3, the five lines "inside" |
| s2.8 | exactly one maximum of u on 0 <= z <= 1090, with radiation, at the fitted quadrature, Lorentzian and general-q points | `one_maximum.out`, six lines (z = 0.473, 0.516, 0.511) |
| s2.8 | "false for c below about 0.01, where the early maximum itself reaches z < 7" | `sec26_table.out`: c = 0.01, 1 + z_max = 7.29 |

Block T3 also prints a rational switch with c_t = 1.6, outside the hypotheses of Corollary 2a, that still has a
maximum. The condition is sufficient, not necessary, and the note does not use that line.

## Section 3, the parts Clawd's code replicates

Table 1 as printed is from Wren's code. Clawd's independent values for its nine rows marked with a dagger:

| Row of Table 1 | chi^2 in Table 1 | Clawd, `c43_fit.out` ("chi2") | z_c in Table 1 | Clawd ("w_crossings_z_lt_3", the grid point below the crossing) |
|---|---|---|---|---|
| LCDM | 1650.79 | 1650.7962 | none | none |
| generalized Lorentzian | 1639.06 | 1639.0548 (c = 2.620, q = 0.916) | 0.511 | 0.5096 |
| w0waCDM | 1639.40 | 1639.3983 (w0 = -0.818, wa = -0.666) | 0.376 | 0.3747 |
| Lorentzian | 1639.10 | 1639.0990 (c = 2.787, Hbar = 50.02) | 0.516 | 0.5149 |
| quadrature | 1642.57 | 1642.5742 (c = 1.641) | 0.473 | 0.4718 |
| linear H-form, free c | 1659.30 | 1659.3215 (c = 2.468) | 0.166 | 0.1664 |
| linear H-form, c = 1 | 1762.12 | 1762.0937 | 1.137 | 1.1360 |
| literal (1.1), free c | 1650.72 | 1650.7191 (c = 185.7; the direction is flat) | none | none |
| literal (1.1), c = 1 | 1912.49 | 1912.4262 | none | none |

Other figures of Section 3 that this folder prints:

| Note | Figure | Script, output line |
|---|---|---|
| s3.1 | matter-only thresholds c_min at Omega_m = 0.31: 1.07267 (quadrature), 1.59423 (Lorentzian) | `c43_props.out`, block T6 |
| s3.2, item 1 | linear H-form worse than LCDM along the profile in c | `c43_scan_lin.out`: twelve values of c from 1 to 30, least chi^2 1659.336 at c = 2.5, against 1650.796. Wren's scan covers 0.05 to 1000 |
| s3.2, item 1 | early maximum at z about 67 with a dark-energy share about 1.1e-5 | `sec26_table.out`, c = 2.469: 1 + z_max 68.16, share 1.07e-05 (at the parameters of Section 2.6, not at the fitted point) |
| s3.3 | x-evolution code: starting point 2.5, 3.0, 3.5 or tolerance 1e-10 to 1e-12 moves H(z) by less than 4e-11; hand-over at z = 8, 20, 60 by less than 3e-11 | `c43_controls.out`, last six lines |
| s3.3, s4.3 | dark-energy fraction at z = 1090: 9.7e-14 (Lorentzian), 1.1e-9 (quadrature) | `c43_fit.out`: "ODE_at_zstar" 9.70e-14 and 1.09e-09 |

## Section 3 as printed (Wren's code)

Table 1 as printed, Table 2, the rows marked with a section sign (the logistic and Gaussian rows, the profile
intervals on c and q, the five secondary dataset combinations), the BAO pulls, the comparison with the published
fits of Lee et al. and of DESI, and the numerical checks of the BDF code are from Wren's code alone. Those
scripts and their outputs are in `wren/`, and `wren/FIGURES.md` is Wren's own map of those figures to output
lines. `wren/c43_sec2_independent.py` also repeats, in a separate code, four of the checks listed above: the four
Lemma 2 values, the Section 2.6 table, the exclusion bound, and the single maximum at the fitted quadrature point.

The agreement between the two codes on the nine rows of the table above (0.0083 at most on six rows, 0.060 at most
on all nine) is the difference between the "chi2" values of `c43_fit.out` and the four-decimal values in
parentheses in `wren/c43_precision_verify.out`, lines 5 to 15.
