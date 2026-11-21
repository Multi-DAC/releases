# Wren's letter 09, claims (a), (b), (c) and the forcedness of the 4-set identity, run with my own code.
# Scenario objects (points, events, orthogonality, M, C) come from refute-2c/base.py; the group H, the images,
# the pair classes and the linear system are built here.
#
# tau_j = swap party j's two boxes and flip the outcome of its partner's box 1.
# Parties 0..3 = A1, B1, A2, B2; partner of j is j ^ 1; wings: Alice {0, 2}, Bob {1, 3}; copies {0, 1}, {2, 3}.
import os, sys, itertools
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from base import *

M = M_func()
C1, C2 = C_func(0), C_func(1)
assert set(M.tolist()) == {-4, 4}
S0 = frozenset(np.nonzero(M == -4)[0].tolist())
print("M takes values", sorted(set(M.tolist())), "; |S| =", len(S0))


def tau(j):
    s = [0] * 4
    s[j] = 1
    f = [[0, 0] for _ in range(4)]
    f[j ^ 1][1] = 1
    return ((0, 1, 2, 3), tuple(s), f)


GEN = {name: tau(j) for j, name in enumerate(["A1", "B1", "A2", "B2"])}
PM = {name: tuple(point_map(g).tolist()) for name, g in GEN.items()}
EM = {name: event_map(g) for name, g in GEN.items()}

# each generator: a permutation of points, of events, keeping orthogonality and PR (x) PR
P1 = pr_q(1)
for name in GEN:
    assert sorted(PM[name]) == list(range(NP)) and sorted(EM[name].tolist()) == list(range(NE))
    em = EM[name]
    assert all(ORTH[em[i], em[k]] == ORTH[i, k] for i in range(NE) for k in range(NE))
    assert all(P1[em[i]] == P1[i] for i in range(NE))
    # the event map and the point map agree: the image of an event's point set is the mapped event's point set
    pm = np.array(PM[name])
    for i in range(0, NE, 7):
        img = np.zeros(NP, dtype=np.int64)
        img[pm[np.nonzero(EV_IND[i])[0]]] = 1
        assert np.array_equal(img, EV_IND[em[i]])
print("each tau_j: permutes points and events, keeps local orthogonality and PR (x) PR: True")


def closure(gens):
    ident = tuple(range(NP))
    seen = {ident}
    frontier = [ident]
    while frontier:
        new = []
        for p in frontier:
            for g in gens:
                q = tuple(g[i] for i in p)  # g after p
                if q not in seen:
                    seen.add(q)
                    new.append(q)
        frontier = new
    return seen


H = closure(list(PM.values()))
print("(a) |H| =", len(H))
D4_1 = closure([PM["A1"], PM["B1"]])
D4_2 = closure([PM["A2"], PM["B2"]])
print("    |<tau_A1, tau_B1>| =", len(D4_1), "; |<tau_A2, tau_B2>| =", len(D4_2))

img = lambda p, S: frozenset(p[w] for w in S)
images = sorted({img(p, S0) for p in H}, key=lambda s: sorted(s))
print("    distinct images of S under H:", len(images))
stab = sum(1 for p in H if img(p, S0) == S0)
print("    stabiliser of S in H:", stab, "; each image hit", {sum(1 for p in H if img(p, S0) == I) for I in images}, "times")
idx = {I: k for k, I in enumerate(images)}
k0 = idx[S0]
perm = {name: tuple(idx[img(PM[name], I)] for I in images) for name in GEN}
print("    S is image number", k0, "; action of the generators on the images:", perm)
print("    tau_B1 and tau_A2 act alike on the images:", perm["B1"] == perm["A2"])
acts = set(perm.values())
ident4 = tuple(range(len(images)))
invol = all(tuple(p[p[i]] for i in range(len(p))) == ident4 for p in acts)
fixfree = all(all(p[i] != i for i in range(len(p))) for p in acts if p != ident4)
print("    distinct actions:", len(acts), "; all involutions:", invol, "; none fixes an image:", fixfree,
      "(three fixed-point-free involutions on four objects that generate a transitive group = Klein four-group)")

# (b) pointwise identities
Mk = [np.array([-4 if w in I else 4 for w in range(NP)], dtype=np.int64) for I in images]
ind = [np.array([1 if w in I else 0 for w in range(NP)], dtype=np.int64) for I in images]
print("(b) sum of the image functions M_k == 2 C1 C2 on all 256 points:", np.array_equal(sum(Mk), 2 * C1 * C2))
print("    4 * (sum of the image indicators) == 8 - C1 C2 on all 256 points:", np.array_equal(4 * sum(ind), 8 - C1 * C2))
print("    control, the single set alone, 16 * 1_S == 8 - C1 C2 (expected False):",
      np.array_equal(16 * ind[k0], 8 - C1 * C2))
print("    values taken by C1 C2:", sorted(set((C1 * C2).tolist())))

# (c) the hand step
a = lambda j, b: ALPHA[:, j, b]
mono = {(x, y): a(0, x) * a(1, y) for x in range(2) for y in range(2)}
signed = {(x, y): (-1) ** (x * y) * mono[(x, y)] for x in range(2) for y in range(2)}
orbit = {tuple(mono[(0, 0)][list(p)].tolist()) for p in D4_1}
print("(c) orbit of a0 b0 under <tau_A1, tau_B1>:", len(orbit), "functions; equal to {a0b0, a0b1, a1b0, -a1b1}:",
      orbit == {tuple(v.tolist()) for v in signed.values()})
for (x, y), f in mono.items():
    avg4 = sum(f[list(p)] for p in D4_1)  # 8 * average
    assert np.array_equal(avg4, 2 * (-1) ** (x * y) * C1), (x, y)
print("    average over the 8 elements of a_x b_y == (-1)^(xy) C1 / 4 for all four (x, y): True")
ssum = sum((-1) ** ((x + x * y + x * xp + x * yp + y * yp + xp * yp) % 2) for x, y, xp, yp in itertools.product(range(2), repeat=4))
print("    sum over (x, y, x', y') of (-1)^(x + xy + xx' + xy' + yy' + x'y') =", ssum)
avgM = sum(M[list(p)] for p in H)
print("    sum over H of M o g == 32 C1 C2 (average C1 C2 / 2) on all 256 points:", np.array_equal(avgM, 32 * C1 * C2))


# The identity for every D: is the functional D -> sum_k D(1_Sk, 1_Sk) fixed by the constraints?
# Unknowns: symmetric D on the 81 monomials (3,321). Row of D(F, G) from monomial coefficients.
iu = [(m, n) for m in range(NM) for n in range(m, NM)]
IU0 = np.array([m for m, n in iu]); IU1 = np.array([n for m, n in iu])
offd = IU0 != IU1


def row(cF, cG):
    r = cF[IU0] * cG[IU1]
    r = r + np.where(offd, cF[IU1] * cG[IU0], 0)
    return r


EC = EV_C.astype(np.float64) / 256.0


def wing_class(i, k, groups):
    """(equal on every group or orthogonal inside it, not all equal): outcomes of one product measurement."""
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
pairs_all = ORTH_PAIRS
pairs_wing = [(i, k) for (i, k) in pairs_all if wing_class(i, k, WINGS)]
pairs_copy = [(i, k) for (i, k) in pairs_all if wing_class(i, k, COPIES)]
pairs_none = [(i, k) for (i, k) in pairs_all if EVENTS[i][0] == EVENTS[k][0]]
print("pair counts: all", len(pairs_all), "; within-wing", len(pairs_wing), "; within-copy", len(pairs_copy), "; same context", len(pairs_none))

cS = [coeffs_int(v).astype(np.float64) / 256.0 for v in ind]
f4 = sum(row(c, c) for c in cS)
f1 = row(cS[k0], cS[k0])
diag = np.array([row(EC[i], EC[i]) for i in range(NE)])


def forced(pairs, label):
    Ap = np.array([row(EC[i], EC[k]) for (i, k) in pairs])
    wp = np.linalg.eigvalsh(Ap.T @ Ap)
    rank_pairs = int((wp >= wp[-1] * 1e-10).sum())
    A = np.vstack([Ap, diag])
    G = A.T @ A
    w, V = np.linalg.eigh(G)
    tol = w[-1] * 1e-10
    null = V[:, w < tol]
    rank = int((w >= tol).sum())
    r4 = np.linalg.norm(null.T @ f4)
    r1 = np.linalg.norm(null.T @ f1)
    # where forced, read the value off: f4 = y^T A, value = y . rhs with rhs = (0 on pairs, P(e) on the diagonal)
    out = "%-12s pairs %d, rank of the pair rows %d ; with the 256 diagonal rows rank %d of %d ; component of the functional outside the row space: 4-set %.3e, single S %.3e (norms of the functionals %.3f, %.3f)" % (
        label, len(pairs), rank_pairs, rank, len(iu), r4, r1, np.linalg.norm(f4), np.linalg.norm(f1))
    if r4 < 1e-8:
        keep = w >= tol
        x = V[:, keep] @ ((V[:, keep].T @ f4) / w[keep])
        y = A @ x
        assert np.linalg.norm(A.T @ y - f4) < 1e-8
        ydiag = y[len(pairs):]
        for q in (1.0, 0.72, 0.5):
            P = np.array([float(p) for p in pr_q(Fraction(q).limit_denominator(1000))])
            val = float(ydiag @ P)
            out += "\n             forced value at PR_q (x) PR_q, q = %.2f: %.12f ; 2 - 4 q^2 = %.12f" % (q, val, 2 - 4 * q * q)
    print(out)


for pairs, label in [(pairs_all, "all"), (pairs_wing, "within-wing"), (pairs_copy, "within-copy"), (pairs_none, "none")]:
    forced(pairs, label)
