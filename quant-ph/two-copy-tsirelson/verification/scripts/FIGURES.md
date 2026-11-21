# Every computed figure in the paper, and the output line that prints it

All scripts are in this folder and are stand-alone apart from `base.py` (points, events, local orthogonality,
the monomial basis of W, PR_q (x) PR_q) and `character_sets.py` (the four sets from the paper's definition),
both here. Python 3 with numpy and python-flint. Each `X.py` was re-run fresh on 2026-11-20 and its output
is `X.out` beside it. Run as `python -B X.py` from any directory.

| Paper | Figure | Script, output line |
|---|---|---|
| s2 | 256 events; 22,400 locally orthogonal pairs (175 per event) | `character_sets.out`: "events: 256 ; locally orthogonal pairs: 22400 ; orthogonal partners per event: [175]" |
| s2, s6 | dim W = 81 | `character_sets.out`: "dim W = exact rank of the 256 event indicators over the integers: 81" |
| s3, Lemma 9 | M_k = +-4 on every point; each S_k has 128 points | `character_sets.out`: "M_k takes values [-4, 4]" (k = 0..3); "|S_k| = [128, 128, 128, 128]" |
| Remark 3 | M_1 = M_0 o tau_A1, M_2 = M_0 o tau_B1 = M_0 o tau_A2, M_3 = M_0 o tau_B2 | `character_sets.out`: four lines "... on all 256 points: True", with the sets; control "(expected False): False" |
| Remark 3 | the group of the four tau_j has order 64; S_0 has exactly four images, the S_k | `character_sets.out`: "order ...: 64", "distinct images of S_0 under it: 4", "they are exactly S_0..S_3: True"; also `w09_fourset.out` "(a) \|H\| = 64" |
| Thm 2 | pointwise: sum_k 1_{S_k} = 2 - C1 C2 / 4 | `character_sets.out`: "4 * sum_k 1_{S_k} == 8 - C1 C2 on all 256 points: True" |
| s4 | sum_k c_k chi_k(s) = 2(-1)^{ab}; phi + ab = xy + x'y' | `character_sets.out`: the two lines "for all 16 s: True" |
| Thm 5 (if) | the construction satisfies (i), (ii), (iii) on all of Omega and reproduces PR_q (x) PR_q | `theorem5_if.out`: one line per q (1/sqrt 2, 0.6, 0.3); floating point, errors at 1e-16 |
| Thm 6 | 18,944 sets of type (5,2,3) | `part_e.out`: "(5,2,3) sets in all: 18944"; `orbits.out`: "(5,2,3) sets in all: 18944" |
| Thm 6 | 18,160 and 640 sets of 4 and 5 pairwise orthogonal good-good events; none of 6 | `part_e.out`: "(f) ... {..., 4: 18160, 5: 640}"; `orbits.out`: "sets of five ...: 640 ; of six: 0" |
| Thm 6 | at most 16 pairwise orthogonal events | `part_e.out`: "largest set of pairwise orthogonal events: 16" |
| Thm 6 (if, at q*) | integer search, scale 1e9, maximum = 144e9 exactly; types (5,2,3), (4,8,4) have excess 0 | `part_end.out`: "checks: (5,2,3) total - one = 0 ; (4,8,4) total - one = 0"; "=> maximum 144000000000 against one = 144000000000 ; difference 0"; "exact maximum <= 1 at q* ...: True" |
| Table 1 | 22,400 / 3,049 / fixed; 8,064 / 2,153 / not fixed (twice); 1,920 / 1,241 / not fixed | `exact_ranks.out`, four lines, exact integer ranks (FLINT); the floating ranks in `w09_fourset.out` agree |
| s6 | at least 18 of the 45 pairs of each (5,2,3) set are not outcomes of one product measurement per wing | `part_e.out`: "pairs (of 45) that are not within-wing, histogram: {18: 256, ...}" (least key 18) |
| s6 | with product measurements only, at most four good-good events pairwise exclusive; CE excludes no isotropic box | `part_res.out`: "restricted to within-wing pairs: GG cliques by size {..., 4: 144}"; q = 1 maximum 0 against 16; same for within-copy; `orbits.out`: "GG sets pairwise within-wing, by size 0..6: [1, 64, 544, 576, 144, 0, 0]" |
| s6 | 3,072 and 6,144 pairs | `orbits.out`: "within-wing not bipartite 3072 ; bipartite not within-wing 6144" |
| s7 | the 640 five-sets fall into classes of 512 and 128 | `orbits.out`: "orbits of the 640 under the stabiliser, sizes: [512, 128]" |
| s7 | no relabelling maps one class to the other | `orbits.out`: "relabellings (of all 98304 ) carrying a set of orbit 0 onto a set of another orbit: 0" |
| s7 | no automorphism of the orthogonality graph does: neighbour counts outside the set | `graph_inv.out`: "invariant separates the classes: True" (twice); class "disjoint" has outside vertices with exactly 1 neighbour, class "share" none |
| s7 | both printed sets of Sainz et al. are in the 512-class | `orbits.out`: "eq. (PRviol) ... in orbit 0"; "the ten-term inequality ... in orbit 0" |
| s7 | extensions to a (5,2,3) set: 30 and 28 | `orbits.out`: "orbit 0 : extensions ... [30]", "orbit 1 : ... [28]" |
| Lemma 10, Thm 2 proof | the surviving d are 0111, 1110, 1111; sum_k chi_k(d) = 0 for each | `lemma10_survivors.out`: "surviving nonzero d: ['0111', '1110', '1111']" (also by brute force); "all three sums zero: True" |

`w09_mono.out` is a side check of the one-copy fact behind the "if, below q*" step: in a single CHSH
scenario every set of pairwise locally orthogonal events has total at most 1 under PR_q for every q in
[0, 1] (the paper gets this from Sainz et al.'s bipartite result instead).

Exact or floating: dim W, Table 1, the counts and the integer search are exact. `theorem5_if.py` and the
floating ranks in `w09_fourset.py` are floating point and are corroboration only; the proof of Theorem 5's
"if" half is the argument in the paper.
