# The one-copy fact the monotonicity lemma needs: in the two-party, two-box scenario every set of pairwise
# locally orthogonal events has total probability at most 1 under PR_q, for every q in [0, 1].
# An event (x, y; a, b) is "good" if a xor b = x*y, with probability (1+q)/4, else (1-q)/4.
# A set with g good and b bad events has total [g(1+q) + b(1-q)]/4, linear in q, so it is at most 1 on [0, 1]
# exactly when g + b <= 4 (q = 0) and g <= 2 (q = 1). All cliques are enumerated.
import itertools

EV = [(x, y, a, b) for x in range(2) for y in range(2) for a in range(2) for b in range(2)]
orth = lambda e, f: (e[0] == f[0] and e[2] != f[2]) or (e[1] == f[1] and e[3] != f[3])
good = lambda e: (e[2] ^ e[3]) == (e[0] & e[1])
types = set()
ncl = 0
for mask in range(1, 1 << 16):
    S = [EV[i] for i in range(16) if (mask >> i) & 1]
    if all(orth(e, f) for e, f in itertools.combinations(S, 2)):
        ncl += 1
        types.add((sum(1 for e in S if good(e)), sum(1 for e in S if not good(e))))
print("nonempty sets of pairwise orthogonal events:", ncl)
print("(good, bad) types:", sorted(types))
print("every type has good + bad <= 4 and good <= 2:", all(g + b <= 4 and g <= 2 for g, b in types))
print("control, the count of orthogonal pairs (expected 16 * 7 / 2 = 56):", sum(1 for e, f in itertools.combinations(EV, 2) if orth(e, f)))
