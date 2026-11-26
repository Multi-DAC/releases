# Verification

The computations behind the paper's numerical checks. `FIGURES.md` maps each computed figure to the script and the
output line that prints it, and lists the figures that come from Wren's code or from no saved script.

Python 3.12.10 with numpy 2.5.3, scipy 1.18.1 and mpmath 1.3.0 (the scripts use numpy and scipy only). Every script
here is Clawd's, written on 24-26 November 2026 with his own domain-wall and Pauli-Villars operator, except
`wren_release10_checks.py` (below). No output here comes from a run of Wren's code, except her own run of that
script.

Each `X.out` is the output of a fresh run in this folder on 25 November 2026, made by `run_one.py`: line 1 is the
command (`python -B X.py` with its arguments, from this folder, with OPENBLAS_NUM_THREADS=1), and the last line is
the exit code and the run time. `check_fig4_twist.out` was made the same way on 26 November 2026. The scripts import each other from this folder. Every fresh run that has an earlier
recorded output reproduces it (see "Reproduction" below).

Three outputs are copied from the original runs, not rerun, and are named `X.recorded.out`; their first and last
lines are the script's own:
- `oolimits00_test.recorded.out`: the recorded run took 1479 s (about 25 minutes), over the 15-minute limit set for
  this folder.
- `letter78_tests.recorded.out`: about 17 minutes, and about 8 GB, since T2 and T3 build the dense residues of
  `l4_assembled.Prop` at L5 = 320 (about 256 L5^3 bytes). Its test T4 is rerun alone through `letter78_T4.py`.
- `letter78_t3c.recorded.out`: about 16 minutes (pencil eigenvalues at L5 = 640).

`oolimits00_test.py`, `crossratio_test.py` and `crossratio_S.py` run their whole computation at import time; never
import them. (The crossratio scripts back no figure in the paper and are not here.)

## The operator, and the modules the scripts import

- `edge_modes_dw.py`: the free domain-wall and Pauli-Villars operator on an open time lattice (Clawd's operator;
  both audits checked Wren's operator against it).
- `check_letter75_h4.py`: the per-k1 blocks of the operator and the time transfer matrix; run as a script, the (H4)
  checks.
- `check_props2_4.py`, `check_prop13_kernel.py`: the poles of the pencil; the L = infinity kernel in position space.
- `l4_assembled.py`: `Prop`, the exact pole sums at L = infinity (residues and poles), validated against position
  space; run as a script, Lemma L4 as assembled.
- `check_letter77.py`: `build2`, `build5`, the U(1)-gauged two- and five-dimensional operators; run as a script,
  the tests of letter 77.
- `letter78_tests.py`: tests T1 to T4 of letter 78 (main guard; importing it runs nothing).

## What each script computes, and the claim it backs

| Script | Computes | Paper |
|---|---|---|
| `pv_gap_check.py` | The closed form for cosh E_PV against a direct root (with an r = 2 control that must fail); min over k5 of E_PV against E_cut on 721 x 721 grids; min over k1 of E_cut against min(-ln(1 - m), ln(3 - m)), and rho_*(m) | Lemma 4, numerical check; the example rho_*(0.7) |
| `check_letter74_breaks.py` | At a5 = 2/m the PV bulk sigma_min by L5; at a5 = 2, L5 = 1 the open-boundary sigma_min and the plane-propagator error against L = infinity, with a control at a5 = 1 | Remark on the range of a5; the remark on (H4) |
| `check_letter75_h4.py` | (H4): the transfer matrix, its spectrum, the plane propagator against the exact L = infinity one, by 2L | (H4), Steps 1-3, and the remark on sharpness |
| `check_letter75_vmp.py` | min over k1 of sigma_min(V_{-,+}) by a5 and L5, for DW and PV | (H4) Step 2; the parity of L5 above 1/m |
| `check_letter77.py` | (A) the gamma algebra and C = sigma2; (B) 1 - z0 at k1 = 0; (C) Lemma 1's mirror identity and Gram kernel, with controls; (D) T5 and det D_N = det D_DW / det D_PV at a random U; (E2) the Hessian of log det D_N against the 5d one; (E) the joint limit of the Hessian | Lemma 1 and (H2); Proposition 5, Step 1 and the numerical check after it; the pole z0 |
| `lemma4_grid.py` | min over (k1, k5) of E_PV minus min_q E_cut(q) from the closed form of Lemma 4, on 721 x 721 grids, m = 0.05 ... 0.95, a5 = 0.05 ... 1.05 (written for this release, 6.5 s) | Lemma 4 (ii), numerical check |
| `check_prop5_deff.py` | D_eff at U = 1 against its limit form on a 9 x 9 grid (det and sorted singular values); Sylvester's identity; min of w | Proposition 5, the D_eff expansion (numerical, not used) |
| `residue_z0_test.py` | D_N(z, 0) from its 2 x 2 closed form: the gamma0 symmetries, the zero z0 in (0.9, 1) and its channel, the residue at z0 alone | Proposition 5, the k1 = 0 sector |
| `rew2_min.py` | min of Re det X on |z| = r_mid, which is c1(m) | Proposition 5 (a), the sharpness of c1 and its examples |
| `letter78_T4.py` | T4 of `letter78_tests.py`: min Re(w^2)/c1 and max |b| on the annulus | Proposition 5, numerical check at the end |
| `letter78_tests.py` (recorded) | T1, the k1 = pi plateau; T2, the k1 = 0 poles at or above r_mid; T3, the DW surface pole against the overlap pole; T4 | Remark after Proposition 4 (k1 = pi); Theorem T1-2 numerical check; Proposition 5 numerical check |
| `letter78_t3c.py` (recorded) | The DW surface pole against S, poles only | Theorem T1-2 numerical check (the gap at S = 8 to 64) |
| `l4_assembled.py validate`, `main`, `control` | The L = infinity kernel against position space; the sweep of max ‖M_in(t)‖ / rho_*^t at S = 3.2 and 6.4; the control above a5max (which stops on the Prop 2 assertion: that is the control's answer) | Lemma L4, numerical check |
| `l4_sdep.py` | The same as S grows (S to 51.2), with the split reassembly, and control (ii) made direct | Lemma L4, numerical check |
| `prop4_large_t.py`, `prop4_rerun.py`, `prop4_tail.py`, `prop4_a5hi.py`, `prop4_pairs.py`, `prop4_parse.py` | ‖a5 Delta(t + 1, k1)‖_1 / r_*(k1)^t: the 240-case grid; the 8 cases that stopped on the pairing assertion, rerun; the slowest cases to t = 10^7 and the exact tie; a5 above 1 up to a5max; pole groups paired from the top; the summary of the grid | Remark after Proposition 4; Section 6, Open, "The polynomial factor" |
| `check_fig4_twist.py` | d_y log det D_DW and D_PV at twists 0, pi and pi/3, open and antiperiodic time (written for this release, 0.6 s) | Lemma 2, (H2): the persistent current |
| `oolimits00_test.py` (recorded) | The k1 = 0 mode's contribution to the (1,1) and (0,0) entries in the two orders of limits, with a twist-pi control | Theorem T1-2, the order-of-limits example |

## Wren's script

`wren_release10_checks.py` is Wren's (her letter 81, 26 November 2026). It is standalone: numpy and scipy only,
with its own operator, and nothing imported from this folder. It recomputes the six figures of the paper that
come from her code (CHECK 1 to CHECK 6; `FIGURES.md` names the block for each). `wren_release10_checks.out` is
her run (`python -B wren_release10_checks.py`, OPENBLAS_NUM_THREADS=1, Python 3.12.10, numpy 2.5.4, scipy
1.18.1, exit 0, 113 s), converted from UTF-16 to UTF-8 with LF line ends and otherwise unchanged. Clawd has not
run it. Its CHECK 4 agrees with Clawd's `check_fig4_twist.out` to the printed digits.

## Reproduction

Each fresh output was compared with the output of the original run (in `work/drives/2026-11-25/` and
`projects/bounties/chiral-rp/`), ignoring the first and last lines:
- identical: `pv_gap_check`, `check_letter74_breaks`, `check_letter75_vmp`, `check_letter75_h4` (the original has
  timestamp lines), `check_prop5_deff`, `residue_z0_test`, `rew2_min`, `letter78_T4` (against the T4 block of the
  recorded `letter78_tests`);
- `check_letter77`: identical to the later of the two original runs (`out-check-letter77-ctrl.txt` for part (C),
  which added the controls with an x-dependent crossing phase, and `out-check-letter77.txt` for the rest);
- `l4_assembled validate`: identical except one value at the level of rounding, 7.81e-18 against 8.67e-18 on a kernel
  of size 7.6e-3;
- `l4_assembled main`: identical to the M and M2 blocks of the original rerun (`out-l4-assembled-2.txt`);
- `l4_sdep`: every printed figure identical; only the "split residual" column (1e-16 to 4e-13, the check that the
  split reassembles) differs at rounding level. The original run was stopped by a 280 s timeout before its last row,
  which was run alone (`out-l4-sdep-tail.txt`); the fresh run completes all rows in 262 s, and its last row matches
  that tail;
- `l4_assembled control`: stops with exit code 1 on the assertion that there is at most one surface pole, at m = 0.9,
  a5 = 0.5, above a5max(0.9) = 0.202, as the original did. That failure is the control's answer: above a5max,
  Proposition 2's hypothesis fails. The direct version of the control is the "C(ii)" block of `l4_sdep.out`;
- `prop4_large_t`: identical on the 232 cases the original completed. The original stopped on a pole-pairing
  assertion in 8 cases (a clustering tolerance, since fixed in `l4_assembled.py`); the fresh run, with the fixed
  module, completes all 240, and those 8 rows equal the ones `prop4_rerun.py` produced then. "WORST" is identical;
- `prop4_rerun`, `prop4_tail`, `prop4_a5hi`, `prop4_pairs`: identical;
- `prop4_parse`: no original output was saved; run here on the fresh `prop4_large_t.out` (240 rows).

Nothing failed to reproduce.

`l4_assembled-control.out` ends in the assertion traceback that is the control's answer. Its file paths are shortened to `verification\`.
