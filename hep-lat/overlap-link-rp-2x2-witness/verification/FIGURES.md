# Every computed figure in the paper, and the output line that prints it

All scripts are in this folder. Python 3.12 with numpy 2.5.3 and mpmath 1.3.0. Each `X.py` was re-run fresh in
this folder on 2026-11-25, and its output is `X.out` beside it. The first line of each `.out` is the command; the
last line is the exit code and the run time. Line numbers below count that first line as line 1.

| Paper | Figure | Script, output line |
|---|---|---|
| Abstract, Prop. 2, s3.2 "Result" | v^T K v in [-8.264675148e-9, -8.232272695e-9]; imaginary part in [-1.62e-11, 1.62e-11]; total radius 1.62e-11 | `certify_2x2_iv.out` line 30: "RESULT v^T K v in Re [-8.264675148e-9, -8.232272695e-9], Im [-1.62e-11, 1.62e-11]; total radius 1.62e-11 ..."; line 31 "CERTIFIED NEGATIVE" |
| Prop. 2 | centres K11 = K22 = 3.784508295226e-4, K12 = K21 = 3.784549537596e-4 | `certify_2x2_iv.out` lines 26-29, "K11 center", "K12 center", "K21 center", "K22 center" (printed to 40 digits; rounded in the paper) |
| s3.2 Part 1 | sigma_min >= 0.6840330 for every Phi; lambda_0 = 0.4679011 | `certify_2x2_iv.out` line 11: "PART1 certified: lambda_min >= [0.4679011582, ...] on the whole circle (sigma floor [0.6840330096, ...])". The grid threshold 0.6848 is an input (`certify_2x2_iv.py`, Part 1) |
| s3.2 Part 2 | delta <= 0.0606, Delta <= 0.4031, c_0 >= 0.06487, M_D <= 7.0970, M_D^16 <= 4.1414e13, n = 124, Q <= 4.0504e-12 | `certify_2x2_iv.out` line 12: "PART2 eta=0.47 delta<=[0.06051008, ...] Delta<=[0.403028, ...] c0>=[0.06487316, ...] M_D<=[7.0969422, ...] M_D^16<=[4.14133e+13, ...] n=124 Q<=4.05031e-12" |
| s3.2 Part 2 | 4Q <= 1.621e-11 on v^T K v | `certify_2x2_iv.out` line 30: "quadrature 1.62e-11"; `certify_2x2_strict.out`: "quad <= [1.6201e-11, ...]" |
| s3.2 Part 3 | determinant error at most 2.7e-28 per node and pair | `certify_2x2_iv.out` line 30: "worst per-point err 2.631e-28" |
| s3.2 "Second computation" | v^T K v in [-8.26467514775e-9, -8.23227269469e-9], at 160 bits | `certify_2x2_strict.out`: "STRICT RESULT v^T K v in Re [-8.26467514775e-9, -8.23227269469e-9]" and "STRICTLY CERTIFIED NEGATIVE"; its own Part 1 and Part 2 lines agree with the first program's |
| s3.4 | 2 sin(0.1) = 0.19967 | by hand (2 sin 0.1 = 0.1996668) |
| s4 table | f(1), ..., f(4) for the overlap, det(D_W - 0.7) and det(D_W + 0.1) | `fd_table.out` lines 5, 13 and 21: "f(d), d = 1..7: ..." |
| s4 | 2f(3) - 2f(4) = -8.2484739e-9, inside the certified enclosure | `fd_table.out` line 9: "= 2 f(3) - 2 f(4) = -8.2484739213e-09"; line 8 gives v^T K v from K directly, -8.2484739214e-09 |
| s4 | spread of K_st at fixed s+t-1 is 2.2e-19; abs(f(d) - f(8-d)) <= 1.1e-19 | `fd_table.out` line 4: "Hankel residual ... 2.17e-19"; line 6: "f(d) - f(8-d), d = 1..3: -5.42e-20  1.08e-19  0.00e+00" |
| s4 | overlap 4 x 4 eigenvalues -6.794442e-6, -3.009393e-6 (largest 1.510052e-3) | `fd_table.out` line 7: "eigenvalues of the 4x4 K: -6.794442e-06  -3.009393e-06  1.699110e-08  1.510052e-03" |
| s4 | "about 800 times the 2 x 2 value" | ratio of line 7's -6.794442e-06 to line 9's -8.2484739213e-09: 823.7 |
| s4 control | det(D_W - 0.7): f(3) - f(4) = +2.6e-4; eigenvalues 6.1e-17, 1.598e-2, 3.968e-2, 4.024 | `fd_table.out` line 13 (f(3) - f(4) = 1.0003223032582 - 1.0000575572252 = 2.647e-4); line 15: "eigenvalues of the 4x4 K: 6.135330e-17  1.598308e-02  3.967672e-02  4.023924e+00" |
| s4 control | det(D_W + 0.1): smallest eigenvalue -2.1e-15, largest 25.6 | `fd_table.out` line 23: "eigenvalues of the 4x4 K: -2.124549e-15 ... 2.559655e+01" |

Not a figure in the paper, filed as a cross-check: `certify_2x2.out` (floating point). Its block (a) compares the
determinant with the holonomy spread over all eight temporal links against the same determinant with the
holonomy on one crossing link (`mirror_orbit.py`): largest relative difference 9.81e-15, which checks Lemma 3(b)
in code. Block (b) samples sigma_min on 4,096 angles (smallest 0.684908, above the certified 0.6840330); block (d)
samples the strip at y = +-eta and finds every sampled quantity inside its bound; block (e) gives the four K
entries by the plain trapezoid rule at n = 123 and 246, equal to the certified centres in every printed digit.
