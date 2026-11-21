# Is the two-class split of the 640 five-sets visible to the orthogonality GRAPH alone?
# Sainz et al.'s "up to symmetry" could mean automorphisms of the graph induced on the 64 possible events,
# a group that may be larger than the relabellings. An invariant computed from adjacency only decides it
# one way: if it differs between the classes, no graph automorphism joins them.
import os, sys, itertools
from fractions import Fraction
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from base import *

P1 = pr_q(1)
GG = [i for i in range(NE) if P1[i] != 0]
ADJ = {i: {k for k in GG if k != i and ORTH[i, k]} for i in GG}
ADJALL = {i: {k for k in range(NE) if k != i and ORTH[i, k]} for i in range(NE)}
print("GG graph: 64 vertices, degrees", sorted({len(ADJ[i]) for i in GG}), "edges", sum(len(ADJ[i]) for i in GG) // 2)


def fives():
    out = []
    def rec(cur, cand):
        if len(cur) == 5:
            out.append(tuple(cur)); return
        for v in sorted(cand):
            rec(cur + [v], {w for w in cand if w > v and w in ADJ[v]})
    rec([], set(GG))
    return out


def nclash(i, k):
    (b1, o1), (b2, o2) = EVENTS[i], EVENTS[k]
    return sum(1 for j in range(4) if b1[j] == b2[j] and o1[j] != o2[j])


def cls(S):
    dbl = [(i, k) for i, k in itertools.combinations(S, 2) if nclash(i, k) == 2]
    assert len(dbl) == 2
    return "share" if set(dbl[0]) & set(dbl[1]) else "disjoint"


def inv(S, adj, verts):
    Sset = set(S)
    outside = sorted(len(adj[v] & Sset) for v in verts if v not in Sset)
    pairs = sorted(len(adj[i] & adj[k]) for i, k in itertools.combinations(S, 2))
    triples = sorted(len(adj[i] & adj[k] & adj[l]) for i, k, l in itertools.combinations(S, 3))
    return (tuple((n, outside.count(n)) for n in sorted(set(outside))), tuple(pairs), tuple(triples))


F = fives()
assert len(F) == 640
for name, adj, verts in [("graph on the 64 GG events", ADJ, GG), ("graph on all 256 events", ADJALL, range(NE))]:
    table = {}
    for S in F:
        table.setdefault(cls(S), set()).add(inv(S, adj, verts))
    print(name)
    for c in sorted(table):
        for v in sorted(table[c]):
            print("   class", c, ": outside vertices by number of neighbours in the set", v[0])
            print("        common neighbours of pairs", v[1])
            print("        common neighbours of triples", v[2])
    print("   invariant separates the classes:", table["share"].isdisjoint(table["disjoint"]))
