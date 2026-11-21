# Wren's letter 09, claims (d), (e), (f), and the endpoint of Theorem 4, run with my own code.
# Events and orthogonality from refute-2c/base.py; the branch and bound is the one in s5_clique_exact.py,
# taking the adjacency as an argument.
#
# Weights. Event types by copy: good on both (GG), on one (mixed), on neither (BB); probabilities
# (1+q)^2/16, (1-q^2)/16, (1-q)^2/16. A set of type (a, b, c) has total [a(1+q)^2 + b(1-q^2) + c(1-q)^2]/16.
#
# Endpoint q* = (sqrt(10) - 1)/3, the positive root of 3q^2 + 2q - 3. There, times 144:
#   GG 14 + 4 r, mixed -2 + 2 r, BB 26 - 8 r, with r = sqrt(10); total 1 <-> 144.
# The search runs on integers W = round(1e9 * weight). Rigour by a gap: for a set of type (a, b, c) with
# a + b + c <= 16 the exact excess over 144 is p + s r with p = 14a - 2b + 26c - 144, s = 4a + 2b - 8c, and if
# it is not zero then |p + s r| = |p^2 - 10 s^2| / |p - s r| >= 1 / (|p| + |s| r) >= 1/965 (|p| <= 560,
# |s| <= 128). Rounding moves a 16-event total by at most 8. So "integer maximum <= 144e9 + 8" proves "exact
# maximum <= 144", because an exact excess would show as at least 1e9/965 - 8 > 1e6.
import os, sys, time, itertools
from math import isqrt
from fractions import Fraction
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from base import *


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


WINGS = [(0, 2), (1, 3)]
COPIES = [(0, 1), (2, 3)]
ADJ = adjacency(lambda i, k: bool(ORTH[i, k]))
ADJ_WING = adjacency(lambda i, k: product_pair(i, k, WINGS))
ADJ_COPY = adjacency(lambda i, k: product_pair(i, k, COPIES))
for name, a in [("all", ADJ), ("within-wing", ADJ_WING), ("within-copy", ADJ_COPY)]:
    print("pairs,", name, ":", sum(bin(m).count("1") for m in a) // 2)
assert all(ADJ_WING[i] & ~ADJ[i] == 0 and ADJ_COPY[i] & ~ADJ[i] == 0 for i in range(NE))


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
                best[0] = nw; best[1] = cur + [v]
            nc = cand & adj[v]
            if nc:
                expand(nw, cur + [v], nc)
            cand &= ~(1 << v)

    expand(0, [], cand)
    return best[0], best[1]


P1 = pr_q(1)
Phalf = pr_q(Fraction(1, 2))
TYPE = []  # 0 GG, 1 mixed, 2 BB
for i in range(NE):
    TYPE.append({Fraction(9, 64): 0, Fraction(3, 64): 1, Fraction(1, 64): 2}[Phalf[i]])
GG = [i for i in range(NE) if TYPE[i] == 0]
assert GG == [i for i in range(NE) if P1[i] != 0]
print("events by type (GG, mixed, BB):", [TYPE.count(t) for t in range(3)])
ALL = (1 << NE) - 1
GGmask = sum(1 << i for i in GG)
NONGG = ALL & ~GGmask


def triple(S):
    return tuple(sum(1 for i in S if TYPE[i] == t) for t in range(3))


def gg_cliques(adj):
    out = {}

    def enum(cur, cand):
        out.setdefault(len(cur), []).append(list(cur))
        for v in list(bits(cand)):
            cand &= ~(1 << v)
            enum(cur + [v], cand & adj[v])

    enum([], GGmask)
    return out


CL = gg_cliques(ADJ)
print("(f) sets of pairwise orthogonal GG events by size:", {k: len(v) for k, v in sorted(CL.items())})
sz, _ = max_clique(ADJ, ALL, [1] * NE)
print("    largest set of pairwise orthogonal events:", sz)


def decide(adj, cl, wtype, one, label):
    """Maximum total weight of a clique, by the number k of GG members. wtype = (GG, mixed, BB) integer weights."""
    W = [wtype[TYPE[i]] for i in range(NE)]
    best = (0, [])
    lines = []
    for k in sorted(cl):
        cap = k * wtype[0] + (16 - k) * wtype[1]
        if cap <= one:
            lines.append("k=%d: at most %d GG + %d mixed = %d <= %d, no search" % (k, k, 16 - k, cap, one))
            continue
        kb, kset = 0, None
        for K in cl[k]:
            cand = NONGG
            for v in K:
                cand &= adj[v]
            ew, es = max_clique(adj, cand, W) if cand else (0, [])
            if k * wtype[0] + ew > kb:
                kb, kset = k * wtype[0] + ew, K + es
        lines.append("k=%d: search over %d GG sets: maximum %d, type %s" % (k, len(cl[k]), kb, triple(kset)))
        if kb > best[0]:
            best = (kb, kset)
    print(label)
    for l in lines:
        print("   " + l)
    print("   => maximum %d against one = %d ; difference %d ; type of a maximiser %s" % (best[0], one, best[0] - one, triple(best[1])))
    return best


def rational(q):
    q = Fraction(q)
    den = q.denominator
    one = 16 * den * den
    w = tuple(int(x * one) for x in ((1 + q) ** 2 / 16, (1 - q * q) / 16, (1 - q) ** 2 / 16))
    assert Fraction(w[0], one) == (1 + q) ** 2 / 16
    return w, one


# (e) every (5,2,3) set: how many of its 45 pairs are not within-wing, not within-copy
t0 = time.time()
n523 = 0
minw, minc = 99, 99
histw, histc = {}, {}
per5 = []
for K in CL[5]:
    cand = NONGG
    for v in K:
        cand &= ADJ[v]
    mixed = [v for v in bits(cand) if TYPE[v] == 1]
    bb = [v for v in bits(cand) if TYPE[v] == 2]
    found = 0
    for m1, m2 in itertools.combinations(mixed, 2):
        if not (ADJ[m1] >> m2) & 1:
            continue
        cb = [v for v in bb if (ADJ[m1] >> v) & 1 and (ADJ[m2] >> v) & 1]
        for trio in itertools.combinations(cb, 3):
            if all((ADJ[x] >> y) & 1 for x, y in itertools.combinations(trio, 2)):
                S = K + [m1, m2] + list(trio)
                nw = sum(1 for x, y in itertools.combinations(S, 2) if not (ADJ_WING[x] >> y) & 1)
                nc = sum(1 for x, y in itertools.combinations(S, 2) if not (ADJ_COPY[x] >> y) & 1)
                histw[nw] = histw.get(nw, 0) + 1
                histc[nc] = histc.get(nc, 0) + 1
                found += 1
    per5.append(found)
    n523 += found
print("(e) GG sets of size 5:", len(CL[5]), "; each extends to a (5,2,3) set:", all(f > 0 for f in per5),
      "; extensions per GG set:", sorted(set(per5)), "; (5,2,3) sets in all:", n523)
print("    pairs (of 45) that are not within-wing, histogram:", dict(sorted(histw.items())))
print("    pairs (of 45) that are not within-copy, histogram:", dict(sorted(histc.items())))
k5w = sorted({sum(1 for x, y in itertools.combinations(K, 2) if not (ADJ_WING[x] >> y) & 1) for K in CL[5]})
k5c = sorted({sum(1 for x, y in itertools.combinations(K, 2) if not (ADJ_COPY[x] >> y) & 1) for K in CL[5]})
print("    among the 10 pairs of the 5 GG events alone: not within-wing", k5w, "; not within-copy", k5c, "(%.0f s)" % (time.time() - t0))

