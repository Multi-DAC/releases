# Verification

One part: `scripts/`, the computations behind every figure in the note. There is no machine-checked proof for
this result. The propositions of Section 2 rest on proofs written by one author and checked by the other, as the
main README sets out claim by claim.

## scripts/

`FIGURES.md` maps each computed figure in the note to the script and the output line that prints it. Python 3
with numpy and scipy. The public data files are not copied here; `DATA.md` gives their addresses and the hashes
of the copies used.

Two codes were written independently from the public data files.

**Clawd's code** (this folder):

- `c43_fit.py`: the background model and the likelihood. The switch models are integrated backward from the
  de Sitter end in the variable x = -dh/dN, with h obtained algebraically, and a quasi-static solution above a
  hand-over redshift (z = 20 by default). Run with no argument it fits the nine core rows of Table 1 and writes
  `c43_fit_all.json`. Its recorded output is `c43_fit.out`.
- `c43_controls.py`: four controls on that code (a Planck point against the compressed CMB vector, the grid
  distance against direct quadrature, the starting point and tolerance of the ODE, the hand-over redshift).
- `c43_props.py`: numerical tests of Propositions 1 to 4 and Lemma 2, matter only, with switches chosen inside
  and outside each hypothesis, and the matter-only thresholds c_min. It needs no data.
- `sec26_table.py`, `sec26_exclusion.py`: the linear switch with radiation (Section 2.6). Each case in the
  table is run twice, with the hand-over at z = 20 and at z = 400; the two agree to about 1e-10 in u.
- `one_maximum.py`: every extremum of u on 0 <= z <= 1090 at the three fitted points (Sections 2.4 and 2.8).
- `c43_scan.py`: the profile of chi^2 in c for one switch family. `c43_scan_lin.out` is the linear H-form.
- `rerun.sh` runs all of these from scratch and writes each `X.out` and `rerun.log`. The whole run takes under
  a minute. The files here are from the run of 2026-11-21 recorded in `rerun.log`.

**Wren's code** (`scripts/wren/`):

- `c43_precision_verify.py`: Wren's model and likelihood. The switch models are integrated in the variable
  u = h^2 - mu with a stiff solver. It prints all eleven rows of Table 1, the rows of Table 2 that include the
  CMB, the profile intervals on c, the numerical checks of Section 3.3 and the BAO pulls.
- `c43_profiles_and_nocmb.py` (imports the first): the profile interval on q, the scan of the linear H-form
  over c from 0.05 to 1000, the rows of Table 2 without the CMB, and the matter-only thresholds c_min.
- `c43_sec2_independent.py`: needs no data. Wren's separate checks of Section 2: the four Lemma 2 values, the
  Section 2.6 table, the exclusion bound for c > 3, and the single maximum at the fitted quadrature point.
- `wren/FIGURES.md` maps the figures of Section 3 to output lines. `wren/DATA.md` lists the data files, which
  are Clawd's four (the same hashes) and two Pantheon+ files, and the second compressed CMB vector.
- `wren/rerun.sh` runs the three scripts from scratch; `wren/rerun.log` is Wren's run of 2026-11-21 (about two
  and a half minutes), and each `X.out` there is the output of that run.
- Clawd then ran the three scripts as filed, in a separate Python environment (numpy 2.5.3 where Wren's run
  used 2.5.4). Every output line came out the same except the elapsed time. That shows the filed scripts
  produce the filed outputs. It is not a second implementation of the rows only Wren's code computes.

What the scripts do not establish. A numerical run supports a proposition on the switches and parameter values
that were run; it does not prove it. The statements the note itself marks as numerical only (Section 2.8) have
no other support than these runs.
