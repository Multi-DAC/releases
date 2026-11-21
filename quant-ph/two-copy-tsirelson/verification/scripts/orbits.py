# Wren's letter 10, section 3 (the 640 five-sets split into two classes) and the bipartite-LO counts of
# its row 41, run with my own code. Events and orthogonality from refute-2c/base.py; nothing of Wren's is read.
#
# The symmetry group is NOT taken as 16 x 16 x 2: every relabelling of the four-party scenario is built
# (a party permutation, and per party: swap the boxes, flip box 0's outcome, flip box 1's outcome;
# 24 * 8^4 = 98,304 maps of the 256 events) and the ones that keep the support of PR x PR are kept.
# If that stabiliser were larger than 512 the two classes could merge, so the check can fail.
import os, sys, time, itertools
from fractions import Fraction
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from base import *

t0 = time.time()


def bits(x):
    while x:
        lb = x & -x
        yield lb.bit_length() - 1
        x ^= lb


def adjacency(pred):
    adj = [0] * NE
    for i in range(NE):
        m = 0
        for k in range(NE):
            if i != k and pred(i, k):
                m |= 1 << k
        adj[i] = m
    return adj


def cliques(adj, cand, k):
    """All k-cliques inside the vertex mask cand, each once, as increasing tuples."""
    out = []

    def rec(cur, cand):
        if len(cur) == k:
            out.append(tuple(cur))
            return
        for v in bits(cand):
            rec(cur + [v], cand & adj[v] & ~((1 << (v + 1)) - 1))

    rec([], cand)
    return out


ADJ = adjacency(lambda i, k: bool(ORTH[i, k]))
P1 = pr_q(1)
Phalf = pr_q(Fraction(1, 2))
TYPE = [{Fraction(9, 64): 0, Fraction(3, 64): 1, Fraction(1, 64): 2}[Phalf[i]] for i in range(NE)]  # GG, mixed, BB
GG = [i for i in range(NE) if TYPE[i] == 0]
assert GG == [i for i in range(NE) if P1[i] != 0] and len(GG) == 64
ALL = (1 << NE) - 1
GGmask = sum(1 << i for i in GG)

# ---- 1. the five-sets and their clash pattern
fives = cliques(ADJ, GGmask, 5)
print("sets of five pairwise orthogonal GG events:", len(fives), "; of six:", len(cliques(ADJ, GGmask, 6)))


def clash_parties(i, k):
    (b1, o1), (b2, o2) = EVENTS[i], EVENTS[k]
    return tuple(j for j in range(4) if b1[j] == b2[j] and o1[j] != o2[j])


def shape(S):
    mult = sorted(len(clash_parties(i, k)) for i, k in itertools.combinations(S, 2))
    dbl = [(i, k) for i, k in itertools.combinations(S, 2) if len(clash_parties(i, k)) == 2]
    where = sorted(clash_parties(i, k) for i, k in dbl)
    if len(dbl) == 2:
        g = "share a vertex" if set(dbl[0]) & set(dbl[1]) else "disjoint"
    else:
        g = "other"
    return tuple(mult), g, tuple(where)


SHAPE = {S: shape(S) for S in fives}
hist = {}
for S in fives:
    key = SHAPE[S][:2]
    hist[key] = hist.get(key, 0) + 1
print("clash multiplicities over the 10 pairs, and how the two double-clash pairs sit:")
for key, n in sorted(hist.items()):
    print("   ", key, ":", n)
wh = {}
for S in fives:
    key = (SHAPE[S][1], SHAPE[S][2])
    wh[key] = wh.get(key, 0) + 1
print("which parties the double clashes are on:")
for key, n in sorted(wh.items()):
    print("   ", key, ":", n)

# ---- 2. every relabelling of the scenario, as a map of event indices
# index of event (b, o): b_0 is bit 7 ... b_3 bit 4, o_0 bit 3 ... o_3 bit 0 (the order of itertools.product)
B = np.array([e[0] for e in EVENTS])
O = np.array([e[1] for e in EVENTS])
assert all(EV_INDEX[EVENTS[i]] == sum(int(B[i, j]) << (7 - j) for j in range(4)) + sum(int(O[i, j]) << (3 - j) for j in range(4)) for i in range(NE))
LOCAL = list(itertools.product(range(2), repeat=3))  # (swap boxes, flip outcome of box 0, flip outcome of box 1)
contrib = {}
for j in range(4):
    for t, (s, f0, f1) in enumerate(LOCAL):
        nb = B[:, j] ^ s
        no = O[:, j] ^ np.where(B[:, j] == 0, f0, f1)
        for p in range(4):
            contrib[j, t, p] = (nb << (7 - p)) + (no << (3 - p))
isGG = np.zeros(NE, dtype=bool)
isGG[GG] = True
GGidx = np.array(GG)
TYPEarr = np.array(TYPE)
group = []
stab = []
for pi in itertools.permutations(range(4)):
    for ts in itertools.product(range(8), repeat=4):
        g = contrib[0, ts[0], pi[0]] + contrib[1, ts[1], pi[1]] + contrib[2, ts[2], pi[2]] + contrib[3, ts[3], pi[3]]
        group.append(g)
        if isGG[g[GGidx]].all():
            stab.append(g)
assert all(len(set(g.tolist())) == NE for g in group[::997])
assert len({g.tobytes() for g in group}) == len(group)
print("relabellings of the scenario:", len(group), "; keeping the support of PR x PR:", len(stab))
print("   each of those keeps the kind (GG, mixed, BB) of every event:", all((TYPEarr[g] == TYPEarr).all() for g in stab))
orth_ok = all(all(ORTH[g[i], g[k]] for i, k in ORTH_PAIRS[::37]) for g in group[::499])
print("   relabellings keep orthogonality (sampled):", orth_ok)

# ---- 3. orbits of the five-sets under the stabiliser
fiveset = {frozenset(S) for S in fives}
seen = {}
orbits = []
for S in fives:
    fs = frozenset(S)
    if fs in seen:
        continue
    orb = set()
    for g in stab:
        img = frozenset(int(g[v]) for v in S)
        assert img in fiveset
        orb.add(img)
    # an orbit of a group is closed; check on the images
    for T in list(orb)[:8]:
        for g in stab[:64]:
            assert frozenset(int(g[v]) for v in T) in orb
    for T in orb:
        seen[T] = len(orbits)
    orbits.append(orb)
print("orbits of the 640 under the stabiliser, sizes:", [len(o) for o in orbits])
for n, orb in enumerate(orbits):
    shapes = {SHAPE[tuple(sorted(T))][:2] for T in orb}
    print("   orbit", n, "size", len(orb), "clash pattern:", shapes)

# ---- 4. the whole scenario group: does any relabelling carry a set of one orbit onto a set of another?
if len(orbits) > 1:
    rep = min(orbits[0], key=lambda T: sorted(T))
    others = set().union(*orbits[1:])
    hits = sum(1 for g in group if frozenset(int(g[v]) for v in rep) in others)
    print("relabellings (of all", len(group), ") carrying a set of orbit 0 onto a set of another orbit:", hits)

# ---- 5. extensions to (5,2,3) sets, by orbit
NONGG = ALL & ~GGmask
ext = {}
for S in fives:
    cand = NONGG
    for v in S:
        cand &= ADJ[v]
    n = sum(1 for T in cliques(ADJ, cand, 5) if sorted(TYPE[v] for v in T) == [1, 1, 2, 2, 2])
    ext[frozenset(S)] = n
tot = 0
for n, orb in enumerate(orbits):
    vals = sorted({ext[T] for T in orb})
    tot += sum(ext[T] for T in orb)
    print("   orbit", n, ": extensions to a (5,2,3) set per five-set:", vals)
print("(5,2,3) sets in all:", tot)


# ---- 6. the sets printed in Sainz et al. (arXiv:1311.6699), P(a1 b1 a2 b2 | x1 y1 x2 y2)
def ev(s):
    o, b = s.split("|")
    return EV_INDEX[(tuple(int(c) for c in b), tuple(int(c) for c in o))]


prviol = [ev(s) for s in ["0000|0000", "1110|0011", "0011|0110", "1101|1011", "0111|1101"]]
ten = [ev(s) for s in ["1111|0000", "1100|1010", "0100|1100", "0011|0001", "0010|0111",
                       "1011|0000", "0101|1100", "1101|1100", "1010|0110", "1001|0100"]]
for name, S in [("eq. (PRviol), five events", prviol), ("the ten-term inequality", ten)]:
    clique = all(ORTH[i, k] for i, k in itertools.combinations(S, 2))
    kinds = tuple(sum(1 for v in S if TYPE[v] == t) for t in range(3))
    gg = frozenset(v for v in S if TYPE[v] == 0)
    where = seen.get(gg, None)
    print(name, ": pairwise orthogonal", clique, "; kinds (GG, mixed, BB)", kinds, "; its GG part is in orbit", where,
          "; clash pattern", SHAPE[tuple(sorted(gg))][:2] if gg in fiveset else None)

# ---- 7. bipartite orthogonality (each wing one party with four boxes and four outcomes) against within-wing
WINGS = [(0, 2), (1, 3)]


def bip(i, k):
    (b1, o1), (b2, o2) = EVENTS[i], EVENTS[k]
    return any(all(b1[j] == b2[j] for j in w) and any(o1[j] != o2[j] for j in w) for w in WINGS)


def product_pair(i, k, groups):
    (b1, o1), (b2, o2) = EVENTS[i], EVENTS[k]
    alleq = True
    for grp in groups:
        same = all(b1[j] == b2[j] and o1[j] == o2[j] for j in grp)
        clash = any(b1[j] == b2[j] and o1[j] != o2[j] for j in grp)
        if not (same or clash):
            return False
        alleq &= same
    return not alleq


ADJ_BIP = adjacency(bip)
ADJ_WING = adjacency(lambda i, k: product_pair(i, k, WINGS))
npairs = lambda a: sum(bin(m).count("1") for m in a) // 2
print("pairs: bipartite", npairs(ADJ_BIP), "; within-wing", npairs(ADJ_WING),
      "; within-wing not bipartite", npairs([ADJ_WING[i] & ~ADJ_BIP[i] for i in range(NE)]),
      "; bipartite not within-wing", npairs([ADJ_BIP[i] & ~ADJ_WING[i] for i in range(NE)]),
      "; bipartite inside all orthogonal pairs:", all(ADJ_BIP[i] & ~ADJ[i] == 0 for i in range(NE)))
print("GG sets pairwise bipartite-orthogonal, by size 0..6:", [len(cliques(ADJ_BIP, GGmask, k)) for k in range(7)])
print("GG sets pairwise within-wing, by size 0..6:", [len(cliques(ADJ_WING, GGmask, k)) for k in range(7)])
print("(%d s)" % (time.time() - t0))


def max_clique(adj, cand, W):
    best = [0, []]

    def expand(cur_w, cur, cand):
        order = sorted(bits(cand), key=lambda v: -W[v])
        classes = []
        for v in order:
            for cl in classes:
                if not (adj[v] & cl[0]):
                    cl[0] |= 1 << v
                    break
            else:
                classes.append([1 << v, W[v]])
        seq = []
        cum = 0
        for cl in classes:
            cum += cl[1]
            for v in bits(cl[0]):
                seq.append((v, cum))
        for v, cb in reversed(seq):
            if cur_w + cb <= best[0]:
                return
            nw = cur_w + W[v]
            if nw > best[0]:
                best[0] = nw
                best[1] = cur + [v]
            nc = cand & adj[v]
            if nc:
                expand(nw, cur + [v], nc)
            cand &= ~(1 << v)

    expand(0, [], cand)
    return best[0], best[1]


# q = 73/100: weights times 160000 are 173^2, 173*27, 27^2; a total of 1 is 160000
Wq = [{0: 173 * 173, 1: 173 * 27, 2: 27 * 27}[t] for t in TYPE]
assert all(Fraction(Wq[i], 160000) == pr_q(Fraction(73, 100))[i] for i in range(0, NE, 7))
sys.stdout.flush()
w, S = max_clique(ADJ_BIP, ALL, Wq)
print("bipartite orthogonality, q = 73/100: largest total", Fraction(w, 160000), "kinds",
      tuple(sum(1 for v in S if TYPE[v] == t) for t in range(3)), "(%d s)" % (time.time() - t0))
