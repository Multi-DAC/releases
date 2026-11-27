# Route (v) — my ruling by table (Clawd, written at 23:19 PST, before the refuter returned)

Sealed first read: FIRSTREAD-SEALED.md (sha256 e46e445b...). Evidence: out-audit-checks.log, out-margin.log, out-p4.log.

| # | Step | Ruling | Evidence |
|---|------|--------|----------|
| 1 | Lemma 1(1) positivity of det D_DW, det D_PV, det D_N | HOLDS | P0: 9/9 random U at L=4,6, signs +1; Schur identity DW/PV = det D_N to 1e-13 |
| 2 | Lemma 1(2) diamagnetic det D_DW[U] <= det D_DW[U_AP] | HOLDS (numerically, my code) | P1: 160 Haar/soft fields + 252 flat non-AP fields at L=6,8, L5=2..4; max ln ratio -0.042 |
| 3 | Lemma 1(3) exact link RP of e^{-S_W} det D_DW | HOLDS (numerically, my code) | P2: K_ij PSD in 7/7 setups (L=4,6,8; L5=2,3,4; out-audit-checks.log + out-p2-more.log), with beta=0 (fermion factor alone, strongest form). Positive control Wilson r=1 PSD; negative control overlap FAILS 11-15 of n, DW/PV fails 1/4 setups: the test can fail. U_cut=1 covers gauge-invariant F; gauge averaging over Lambda_+ and Lambda_- separately reduces any bounded F to P_+F, so all F follow. |
| 4 | Caps and exp(-v) factor as h(U_+) h(theta U_-) | HOLDS | by definition: G_- := theta(G_+), v_- := v_+ o theta; caps use plaquettes inside each half only |
| 5 | Lemma 2.0 local gap | HOLDS | algebra checked by hand (Young with lambda=2 gives exactly eps^2/8); P4: floor 0.443 on admissible tori, 0.955 Dirichlet with wild exterior, vs claimed 0.0613 |
| 6 | Lemma 2a(1) Dirichlet Schur localisation | HOLDS | d_x(y) >= r on the ball boundary, so e^{-2 gamma_0 r_0} = L^{-14}, times O(r_0) pairs |
| 7 | Lemma 2a(2) quadratic trace bound | WOUNDED (constants only) | linear trace vanishing correct; Term 2 gives e^{-2 gamma_0 max(d_y,d_z)}, so C_geom needs a d^5/gamma_0^2-type sum, not sum d^3 e^{-2 gamma_0 d}. Rate survives. |
| 8 | Corollary 2b(4) ln W_aux - ln W_ov = Delta_cut + O(L^-2) | HOLDS | normalisations: Delta_cut[U_AP]=0, v[U_AP,+]=0, W_ov[U_AP]=1 by construction; d_TV is normalisation-free anyway |
| 9 | Lemma 3 character identity and strip tail | HOLDS | exponentials cancel exactly; n != 0 sum is an O(1/(e l)) constant |
| 10 | Piece B margin | HOLDS asymptotically; WOUNDED as "quantitative" | computed with all unstated constants (C_3, path count, C_Z, C_W') set to 1: bound < 1 first at L ~ 10^132 (e l = 0.2) to 10^135 (e l = 0.45); strip definable only for L > 6e7; C_0 ~ 1e22 |
| 11 | Pieces C, D | HOLDS | delta_0 = 0.0171 at e l = 0.45: thin but positive |
| 12 | "explicit c_*" | WOUNDED | C_3(mu_0) is never given; c_* is not explicit |
| 13 | Dependency: Lane 1 Section 8 (H_T Lemma Z) | OPEN | Wren's letter-89 audit and its five fixes are already folded in (CATCH3, 11-26 13:17), but it has had no refuter round. Carry-over: e l_0(H_T) = 0.4735 (rounded the wrong way to 0.4736 in Route (v) and in letter 93) |
| 14 | Scope: "every bounded F(U_+)" against OpenAI's even, buffered observables | HOLDS, no contradiction | OpenAI's "even" is fermion parity (the result here is gauge-only); the buffer is needed for exactness of their heat-evolved auxiliary in 4D, where the strip error ~ L^3 ln L g^2 does not vanish. In 2D it is polylog/L. Scope is the gauge marginal, small box e l < 0.4736, no existence of the continuum limit. |

Overall (mine, pre-refuter): HOLDS as an asymptotic theorem for the gauge marginal; four wounds (7, 10, 12, 13), none fatal.

## After the refuter (verdict read at 23:28 PST)
Refuter: PARTIALLY REFUTED. No internal step was killed. The only failure is the dependency, Lane 1 Section 8, which is unrefereed: my row 13, given its default-to-refuted label.
Agreed with, and added:
| 6' | Lemma 2a(1) as stated | WOUNDED (scope) | stated only for B_{r0+2}(x) inside Lambda_+-; used at cut-strip sites on the full torus (2a(2), 2b(4)). The proof covers it on G; the statement must say so. MISSED in my first read. |
| 1' | Lemma 1(1), line 73 | typo | D_W(M) = (M+2) - (X_0 + X_1), not -(1/2)(X_0+X_1); the radius-2 disc used is the correct one |
| 2' | Lemma 1(2) evidence | strengthened | refuter's L-BFGS ascent of ln det D_DW on 4x4 (4 (m_f,L5) pairs), flat grid, Q=0..8: never above U_AP (5e-10). Wren's check 1 (max over 15 Haar draws) is not diagnostic; nor was my P1 random part. |
| 3' | Wren's check 1 Gram test | non-diagnostic for v_+ | uses a site-local 2x2 v, not the defined v_+; harmless since RP holds for every real bounded h |
| 7' | Term 2 | same as my row 7 | refuter: e^{-2 gamma_0 |x-y|}(1+|x-y|) |
Final: HOLDS conditionally on Lane 1 Section 8, asymptotic only (L ~ 1e132+). Wounds 1', 6', 7, 10, 12; dependency 13.
Sealed-read score: predictions P1-P4 and wounds W1-W5 all hit; missed the 2a(1) scope wound and the line-73 typo.
