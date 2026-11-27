# Refuter round 1 on LANE1-Z-LEMMA.md (sections 0 to 7), returned at clock reading 12:41, 2026-11-26

Verdict, from the refuter's report as I read it: SURVIVES (HOLDS), with four minor wounds. It said none of them
was load-bearing and gave a repair for each. Its own scripts are r1.py to r5.py in this folder.

1. Caffarelli, convex-support case (1d). Repair: extend the cosine potential past eps0 by its C^1 quadratic
   continuation, so the second derivative is at least kappa. Add the penalty n sum (|theta_p| - eps0)_+^2, apply
   finite-potential Caffarelli, and pass to the weak limit (uniformly Gaussian tails).
2. Notation at the hopping-difference line: Delta_l had to name its link, as Delta_{(x,mu)} and Delta_{(x-mu,mu)}.
3. Section 6: "iff" -> "if". It is the method's ceiling, not a physical one. eps = 0.267 is unverified for this
   kernel; its uniform-flux probe (L = 20, |1 - U_p| up to 0.50) had min|eig H_W| >= 0.458.
4. Overclaims: the status line said "Steps 1 to 5"; line 132 claimed the tails were monotone in L without showing it;
   V meant both the volume and the potential in the Caffarelli line.

All four were repaired in LANE1-Z-LEMMA.md (status line, 1d, Step 2, the line after the 1d table, section 6,
section 7). Section 8 (H_T) was not in this round.
