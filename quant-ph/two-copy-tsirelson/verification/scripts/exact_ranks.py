# Table 1 of the paper, in exact arithmetic. Unknowns: a symmetric form D on W, as its 3,321 entries on the
# 81 monomials. Rows: D(1_e, 1_e') = 0 for each decohering pair of the measurement list, and the 256 diagonal
# rows D(1_e, 1_e) (whose right-hand sides are P(e)). The four-set sum is fixed by these conditions exactly
# when its functional lies in the row space, that is when appending it does not raise the rank. All entries
# are integers (monomial coefficients scaled by 256), and ranks are computed over the integers by FLINT, so
# they are exact. Pair lists as in w09_fourset.py: all locally orthogonal pairs; pairs that are outcomes of
# one product of a branching measurement per wing; the same per copy; pairs in the same context.
import os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import flint
from base import NE, EVENTS, ORTH_PAIRS, EV_C, NM, coeffs_int
from character_sets import S  # the four sets, from the character definition

which = sys.argv[1:] or ["all", "within-wing", "within-copy", "none"]

iu = [(a, b) for a in range(NM) for b in range(a, NM)]
IU0 = np.array([a for a, b in iu]); IU1 = np.array([b for a, b in iu])
offd = IU0 != IU1


def row(cF, cG):
    r = cF[IU0] * cG[IU1]
    return r + np.where(offd, cF[IU1] * cG[IU0], 0)


def wing_class(i, k, groups):
    (b1, o1), (b2, o2) = EVENTS[i], EVENTS[k]
    alleq = True
    for grp in groups:
        same = all(b1[j] == b2[j] and o1[j] == o2[j] for j in grp)
        clash = any(b1[j] == b2[j] and o1[j] != o2[j] for j in grp)
        if not (same or clash):
            return False
        alleq &= same
    return not alleq


WINGS = [(0, 2), (1, 3)]
COPIES = [(0, 1), (2, 3)]
LISTS = {
    "all": ORTH_PAIRS,
    "within-wing": [(i, k) for (i, k) in ORTH_PAIRS if wing_class(i, k, WINGS)],
    "within-copy": [(i, k) for (i, k) in ORTH_PAIRS if wing_class(i, k, COPIES)],
    "none": [(i, k) for (i, k) in ORTH_PAIRS if EVENTS[i][0] == EVENTS[k][0]],
}
EC = EV_C.astype(object)
diag = [row(EC[i], EC[i]) for i in range(NE)]
cS = [coeffs_int(x.astype(np.int64)).astype(object) for x in S]
f4 = sum(row(c, c) for c in cS)
print("unknowns:", len(iu))
for label in which:
    pairs = LISTS[label]
    t0 = time.time()
    rows = [row(EC[i], EC[k]) for (i, k) in pairs] + diag
    A = flint.fmpz_mat([[int(v) for v in r] for r in rows])
    rA = A.rank()
    Af = flint.fmpz_mat([[int(v) for v in r] for r in rows + [f4]])
    rAf = Af.rank()
    print("%-12s decohering pairs %5d ; exact rank with the 256 diagonal rows %d ; with the four-set functional %d ; four-set sum %s (%.0f s)"
          % (label, len(pairs), rA, rAf, "FIXED" if rAf == rA else "NOT FIXED", time.time() - t0), flush=True)
print("DONE")
