# Independent construction of the two-copy CHSH scenario, from letter 150's text only.
# Parties j = 0..3 are A1, B1, A2, B2. Copy of j is j // 2, partner is j ^ 1.
# A point w in {0,1}^8: value of party j, box b is bit (2j + b).
# alpha_{j,b}(w) = 1 - 2 * value  (so value 0 <-> +1).
# Event: boxes (b_0..b_3) and outcomes (o_0..o_3); indicator of points with v[j][b_j] = o_j for all j.
# Monomial basis of W: t_j in {0,1,2}: 0 -> 1, 1 -> alpha_{j,0}, 2 -> alpha_{j,1}; 81 monomials.
import itertools
from fractions import Fraction
import numpy as np

NP = 256


def val(w, j, b):
    return (w >> (2 * j + b)) & 1


VAL = np.array([[[val(w, j, b) for b in range(2)] for j in range(4)] for w in range(NP)])  # [w, j, b]
ALPHA = 1 - 2 * VAL  # [w, j, b]

EVENTS = []
for bs in itertools.product(range(2), repeat=4):
    for os_ in itertools.product(range(2), repeat=4):
        EVENTS.append((bs, os_))
EV_INDEX = {e: i for i, e in enumerate(EVENTS)}
NE = len(EVENTS)
EV_IND = np.zeros((NE, NP), dtype=np.int64)
for i, (bs, os_) in enumerate(EVENTS):
    m = np.ones(NP, dtype=bool)
    for j in range(4):
        m &= VAL[:, j, bs[j]] == os_[j]
    EV_IND[i] = m


def orthogonal(e1, e2):
    (b1, o1), (b2, o2) = e1, e2
    return any(b1[j] == b2[j] and o1[j] != o2[j] for j in range(4))


ORTH = np.zeros((NE, NE), dtype=bool)
for i in range(NE):
    for k in range(NE):
        ORTH[i, k] = orthogonal(EVENTS[i], EVENTS[k])
ORTH_PAIRS = [(i, k) for i in range(NE) for k in range(i + 1, NE) if ORTH[i, k]]

MONOS = list(itertools.product(range(3), repeat=4))
MONO_INDEX = {t: i for i, t in enumerate(MONOS)}
NM = len(MONOS)
MONO_VEC = np.ones((NM, NP), dtype=np.int64)
for mi, t in enumerate(MONOS):
    for j in range(4):
        if t[j] > 0:
            MONO_VEC[mi] *= ALPHA[:, j, t[j] - 1]


def coeffs_int(F):
    """Exact monomial coefficients (times 256) of an integer function F; asserts F lies in W."""
    F = np.asarray(F, dtype=np.int64)
    c = MONO_VEC @ F
    rec = c @ MONO_VEC
    assert np.array_equal(rec, 256 * F), "function not in W"
    return c


EV_C = np.array([coeffs_int(EV_IND[i]) for i in range(NE)])  # scaled by 256


def pr_q(q):
    """P(e) for PR_q (x) PR_q as Fractions. PR: P(a,b|x,y) = 1/2 if a xor b = x*y."""
    q = Fraction(q)
    P = []
    for bs, os_ in EVENTS:
        p = Fraction(1)
        for c in range(2):
            x, y = bs[2 * c], bs[2 * c + 1]
            a, b = os_[2 * c], os_[2 * c + 1]
            pr = Fraction(1, 2) if (a ^ b) == (x & y) else Fraction(0)
            p *= q * pr + (1 - q) * Fraction(1, 4)
        P.append(p)
    return P


def M_func():
    M = np.zeros(NP, dtype=np.int64)
    for x, y, xp, yp in itertools.product(range(2), repeat=4):
        sgn = (-1) ** ((x + x * xp + x * yp + y * yp) % 2)
        M += sgn * ALPHA[:, 0, x] * ALPHA[:, 1, y] * ALPHA[:, 2, xp] * ALPHA[:, 3, yp]
    return M


def C_func(c):
    C = np.zeros(NP, dtype=np.int64)
    for x, y in itertools.product(range(2), repeat=2):
        C += (-1) ** (x * y) * ALPHA[:, 2 * c, x] * ALPHA[:, 2 * c + 1, y]
    return C


# relabeling g = (sigma, s, f): party j -> sigma[j]; box b -> b ^ s[j]; outcome of (j,b) xor f[j][b]
def point_map(g):
    sigma, s, f = g
    out = np.zeros(NP, dtype=np.int64)
    for w in range(NP):
        nw = 0
        for j in range(4):
            for b in range(2):
                v = val(w, j, b) ^ f[j][b]
                nw |= v << (2 * sigma[j] + (b ^ s[j]))
        out[w] = nw
    return out


def event_map(g):
    sigma, s, f = g
    out = np.zeros(NE, dtype=np.int64)
    for i, (bs, os_) in enumerate(EVENTS):
        nb = [0] * 4
        no = [0] * 4
        for j in range(4):
            nb[sigma[j]] = bs[j] ^ s[j]
            no[sigma[j]] = os_[j] ^ f[j][bs[j]]
        out[i] = EV_INDEX[(tuple(nb), tuple(no))]
    return out


def all_relabelings():
    for sigma in itertools.permutations(range(4)):
        for s in itertools.product(range(2), repeat=4):
            for fl in itertools.product(range(2), repeat=8):
                f = [[fl[2 * j], fl[2 * j + 1]] for j in range(4)]
                yield (sigma, s, f)


def pr_group():
    """All relabelings (24 party perms x box swaps x outcome flips = 98304) that keep PR (x) PR."""
    P1 = pr_q(1)
    nz = np.array([p != 0 for p in P1])
    # fast event map: vectorised
    B = np.array([e[0] for e in EVENTS])
    O = np.array([e[1] for e in EVENTS])
    keep = []
    for g in all_relabelings():
        sigma, s, f = g
        idx = np.zeros(NE, dtype=np.int64)
        fa = np.array(f)
        for j in range(4):
            nb = B[:, j] ^ s[j]
            no = O[:, j] ^ fa[j][B[:, j]]
            idx += (nb * 16 + no) * 0  # placeholder to keep shape
            idx += 0
        # EV index layout: boxes major (b0 b1 b2 b3), then outcomes (o0..o3), each in product order
        nbs = np.zeros((NE, 4), dtype=np.int64)
        nos = np.zeros((NE, 4), dtype=np.int64)
        for j in range(4):
            nbs[:, sigma[j]] = B[:, j] ^ s[j]
            nos[:, sigma[j]] = O[:, j] ^ fa[j][B[:, j]]
        idx = ((nbs[:, 0] * 8 + nbs[:, 1] * 4 + nbs[:, 2] * 2 + nbs[:, 3]) * 16
               + nos[:, 0] * 8 + nos[:, 1] * 4 + nos[:, 2] * 2 + nos[:, 3])
        if np.array_equal(nz[idx], nz):
            keep.append(g)
    return keep
