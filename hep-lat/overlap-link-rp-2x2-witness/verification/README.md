# Verification

The computations behind every figure in the paper. `FIGURES.md` maps each computed figure to the script and the
output line that prints it.

Python 3.12 with numpy 2.5.3 and mpmath 1.3.0. Run each script as `python -B X.py` from this folder; the scripts
import each other from the same folder. Each `X.out` is the output of a fresh run here on 2026-11-25. Its first
line is the command, with the arguments used, and its last line is the exit code and the run time.

## The certificate (Proposition 2)

- `certify_2x2_iv.py` (Clawd). Parts 1 to 3 of Section 3.2 of the paper, in mpmath's interval arithmetic (`iv`).
  Run as `certify_2x2_iv.py 512 0 0.47 128`: 512 grid angles for Part 1, the number of quadrature nodes chosen by
  the script (124), strip half-width 0.47, 128-bit working precision. About 3 minutes. sha256 prefix
  `7081cee0339fc1cf`.
- `certify_2x2_strict.py` (Wren). The same certificate, written separately, sharing no code with the first. mpmath
  1.3.0's interval `exp`, `cos` and `sin` are not formally outward-rounded, so this program does not use them, nor
  mpmath's `pi`, `det` or decimal string conversion. It builds pi from Machin's series with a tail bound, and exp,
  cos and sin from Taylor polynomials with Lagrange remainders, using interval `+ - * /` only. It takes exact
  rational inputs, uses its own interval Gaussian elimination with an explicit nonzero-pivot test, and sums all
  error terms in interval arithmetic. Run as `certify_2x2_strict.py 160` (160 bits). About 3 minutes. sha256
  prefix `5411e3c8692a1eac`. The docstring's references to "Letter 257" and "Track B Target 1" point to our
  correspondence, where the first program was sent to Wren for audit.

Both programs build the 16 x 16 kernel X = D_W - m for each of the four pairs directly from the definition in
Section 2 of the paper. What they share is that definition and the design of the bounds; Section 4 of the paper
(below) is the test of the definition.

Both runs here reproduce, digit for digit, the enclosures the two programs gave in their earlier recorded runs.

## Floating-point checks (Section 4 and cross-checks)

- `fd_table.py`: the Hankel structure, f(d) for d = 1..7, the 4 x 4 Gram matrices, and the positive control with
  the Wilson determinant, at beta = 0 with a 256-point trapezoid rule. It imports `mirror_test.py` (the Wilson
  and overlap operators in a U(1) background on a 2L x Ls lattice, used here at 8 x 1), `mirror_kernel.py` (the
  links of a reflected pair) and `certify_2x2.py` (the witness parameters).
- `certify_2x2.py`: the floating-point side of the certificate (the bounds sampled numerically, and the four K
  entries by the trapezoid rule). It imports `mirror_test.py`, `mirror_kernel.py` and `mirror_orbit.py` (the
  determinant with the holonomy on one crossing link, against which the spread holonomy is compared).

Nothing in this section is certified. The control is what tests the code that builds the operator and applies the
reflection: the Wilson determinant must give a positive semidefinite Gram matrix, and it does, to rounding.
