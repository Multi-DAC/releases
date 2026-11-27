# Sealed check of Wren letter 89's numeric claims (clock reading 13:19:07, before the code)

- N1: c_f(H_T, m = 0.7) on even midpoint grids N = 512 to 4096 extrapolates (O(N^-2) Richardson) to 0.1589390 +- 3e-7.
- N2: the odd n = 4001 grid in ht_free.py puts a point at k ~ 1e-16 (k = 0 up to rounding). Its log term alone shifts
  c_f by +3e-6 to +6e-6. Dropping that one point brings 4001 to within 1e-6 of N1's limit.
- N3: min |H_T| at m = 0.5 on a fine grid equals m/(2 - m) = 1/3 to 1e-6. The f(4) bound at m = 0.5 exceeds it
  (Wren: 0.4370).
- N4: mu_T(0.7) = sqrt(2.89/20.09) is 0.37928 to 5 decimals, so the draft's "mu_T >= 0.3793" rounds a lower bound UP.
- N5: 2 sqrt2/(mu_T (2 - m)^2) = 4.4126 to 4 decimals (Wren), and sqrt2/m = 2.0203 (so the draft's "2.02" also
  rounds an upper bound down).
Any miss on N1 or N3 means Wren's fix is wrong and does not go in as she wrote it.
