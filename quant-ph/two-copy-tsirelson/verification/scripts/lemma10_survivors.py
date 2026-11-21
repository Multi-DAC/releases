"""Lemma 10: which nonzero d in F_2^4 survive, and the character sums of Theorem 2's quadratic terms.

d = (d_x, d_y, d_x', d_y'). phi(s) = x + x x' + x y' + y y'; B(s, d) = phi(s+d) + phi(s) + phi(d).
A nonzero d contributes only if B(e_j, d) = 0 for every coordinate j with d_j = 0.
Stand-alone: Python 3, no imports beyond itertools.
"""
from itertools import product


def phi(s):
    x, y, xp, yp = s
    return (x + x * xp + x * yp + y * yp) % 2


def add(s, t):
    return tuple((a + b) % 2 for a, b in zip(s, t))


def B(s, d):
    return (phi(add(s, d)) + phi(s) + phi(d)) % 2


E = [tuple(1 if i == j else 0 for i in range(4)) for j in range(4)]
V = list(product((0, 1), repeat=4))

# B is bilinear: check it, so the survivor test below is the lemma's test
bilinear = all(B(add(s, t), d) == (B(s, d) + B(t, d)) % 2 for s in V for t in V for d in V)
print("B additive in s for all s, t, d:", bilinear)
nondeg = [d for d in V if all(B(s, d) == 0 for s in V)]
print("d with B(., d) = 0 identically (non-degenerate iff only 0000):", ["".join(map(str, d)) for d in nondeg])

# direct form of the lemma: summing (-1)^B(s,d) over the box of a party j outside d gives zero unless B(e_j, d) = 0
surv = [d for d in V if any(d) and all(B(E[j], d) == 0 for j in range(4) if d[j] == 0)]
print("surviving nonzero d:", ["".join(map(str, d)) for d in surv])

# brute force of the same statement: the s-sum over the coordinates outside d of (-1)^B(s, d) is nonzero
brute = []
for d in V:
    if not any(d):
        continue
    out = [j for j in range(4) if d[j] == 0]
    tot = 0
    for bits in product((0, 1), repeat=len(out)):
        s = [0, 0, 0, 0]
        for j, b in zip(out, bits):
            s[j] = b
        tot += (-1) ** B(tuple(s), d)
    if tot != 0:
        brute.append(d)
print("same by summing over the boxes outside d:", ["".join(map(str, d)) for d in brute], brute == surv)


def chi(k, s):
    x, y, xp, yp = s
    e = [0, y + xp + yp, x + yp, x + y + xp][k]
    return (-1) ** (e % 2)


for d in surv:
    print("d =", "".join(map(str, d)), ": chi_k(d) for k = 0..3:", [chi(k, d) for k in range(4)],
          "; sum:", sum(chi(k, d) for k in range(4)))
print("all three sums zero:", all(sum(chi(k, d) for k in range(4)) == 0 for d in surv))
