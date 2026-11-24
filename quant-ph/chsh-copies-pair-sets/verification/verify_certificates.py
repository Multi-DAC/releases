"""One verifier for the four identity files of Theorem 9 of the paper (documentation/paper.pdf). Self-contained:
python 3 with numpy, no other file of code.

The P line checks, besides the polynomial itself: p(-1) > 0, p(0) > 0, p(1) < 0, that p is square-free, and that
p has exactly one root in (-1, 1] ; the CHSH bound is printed as the six-decimal ROUND-UP b of 4 q*, certified
exactly by p(b / 4) < 0 and p((b - 10^-6) / 4) > 0.

What a certificate is (the uniform statement). n copies, N = 16^n points, W_n spanned by the 9^n Walsh monomials.
G = the 16^n n! point maps (in each copy one of the 16 CHSH relabellings that fix PR_q for every q ; times the
copy orders). Data: sets f_j, rationals y_j > 0, and a pair set (nb or mid). The claim of a file:

    for every real bilinear K on W_n with (i) K(1_e, 1_e) = P_q(e) on the product events of PR_q^(x n) and
    K(1_e, 1_e') = 0 on the pair set:      sum_j y_j avg_{g in G} K(1_{g f_j}, 1_{g f_j}) = p(q),

with p the polynomial stored in the file. If every f_j is spanned, a K that is also >= 0 on spanned sets needs
p(q) >= 0. With c(f) = the integer Walsh coordinates of 1_f, the claim follows from one of two matrix facts:

  nb  (files without `pairs`):  Chat = sum_j y_j sum_g c(g f_j) c(g f_j)^T vanishes on every INCOMPATIBLE pair of
      monomials (some wing of some copy named at two different settings). On a compatible pair K(m, m') is the
      expectation <m m'> in PR_q^(x n): both monomials are combinations of events of one common setting, so only
      (i) and the non-branching zeros are used. p(q) = sum over compatible pairs of Chat <m m'> / (|G| N^2).
  mid (files with `pairs`):  2 sum_j y_j sum_g c(g f_j) c(g f_j)^T
                               = sum_c mu_c sum_g [ c(g e_c) c(g e'_c)^T + c(g e'_c) c(g e_c)^T ],
      each (e_c, e'_c) a pair of product events that is either equal or a middle pair (same setting pair and
      different outcome pair in some copy) ; p(q) = sum over the equal pairs of mu_c P_q(e_c).
      mu = lambda as stored (two copies) ; mu = lambda / 4096 (three copies ; the stored normalisation).

What is checked, per file (stated before any run):
  S   every stored set is spanned (B B^T f = N f) ; weights > 0 ;
  G   the group: 128 one-copy relabellings, 16 fix PR_q ; each of the 16 maps events to events, keeps the sign of
      every event and the relation "same setting pair" ; |G| equals the file's group_size ; the coordinate rule
      c(g f) = signed permutation of c(f) equals the coordinates of the LITERAL image set, for the first set of
      the file (all g at two copies ; 2,004 seeded draws of g, with replacement, at three copies) ;
  I   the matrix fact, entry by entry, in integers (python integers for the weighted sums) ;
  P   the polynomial read off equals the file's `poly` exactly ; p(1) < 0 ; p > 0 at the lower end of the
      bracket and < 0 at the upper end ; by a Sturm count p has no root in (upper end, 1] ;
  E   (mid) every stored pair is equal or a middle pair, by the definition.
Controls that can fail, per file:
  K1  the first weight + 1/1000: I or P must fail ;
  K2  the first set with point 0 flipped: S must fail ;
  K3  (mid) the first non-equal pair with e' replaced by e with the setting of both wings changed in every copy:
      E must fail ;
  K4  (two copies) the coordinate rule with one sign of one relabelling changed: G must fail.
Exactness: Gram matrices are computed in float64 on integer arrays ; every partial sum is an integer below
16^n n! 16^(2n) < 2^53, so the products are exact ; they are converted to integers and asserted integral.

File conventions read (not changed): two copies, point L = u_A1 + 4 u_B1 + 16 u_A2 + 64 u_B2, u = 2 n_0 + n_1,
events [x, a] over the four wings ; three copies, point 256 k1 + 16 k2 + k3, k = 8 a0 + 4 a1 + 2 b0 + b1, events
(x, a, y, b) per copy.

Usage: python verify_certificates.py nb nb_n2.json mid mid_n2.json nb nb_n3.json mid mid_n3.json
"""
import itertools
import json
import math
import random
import sys
import time
from fractions import Fraction as Fr

import numpy as np

PTS1 = list(itertools.product((0, 1), repeat=4))                    # (a0, a1, b0, b1) ; index 8a0+4a1+2b0+b1
MON1 = [(sa, sb) for sa in (None, 0, 1) for sb in (None, 0, 1)]     # alpha_sa * beta_sb
EV1 = [(x, a, y, b) for x in (0, 1) for a in (0, 1) for y in (0, 1) for b in (0, 1)]


def mval(m, p):
    v = 1
    if m[0] is not None:
        v *= (-1) ** p[m[0]]
    if m[1] is not None:
        v *= (-1) ** p[2 + m[1]]
    return v


B1 = np.array([[mval(m, p) for m in MON1] for p in PTS1], dtype=np.int64)   # 16 x 9


def ev_points(e):
    x, a, y, b = e
    return frozenset(i for i, p in enumerate(PTS1) if p[x] == a and p[2 + y] == b)


def ev_sign(e):
    x, a, y, b = e
    return (-1) ** (a + b + x * y)


def one_copy_group():
    idx = {p: i for i, p in enumerate(PTS1)}
    maps = set()
    for sw, ta, tb in itertools.product((0, 1), repeat=3):
        for fl in itertools.product((0, 1), repeat=4):
            img = []
            for p in PTS1:
                a, b = (p[0], p[1]), (p[2], p[3])
                if sw:
                    a, b = b, a
                if ta:
                    a = (a[1], a[0])
                if tb:
                    b = (b[1], b[0])
                img.append(idx[(a[0] ^ fl[0], a[1] ^ fl[1], b[0] ^ fl[2], b[1] ^ fl[3])])
            maps.add(tuple(img))
    assert len(maps) == 128
    where = {ev_points(e): e for e in EV1}
    assert len(where) == 16
    keep = []
    for pi in sorted(maps):
        ok = True
        imgs = {}
        for e in EV1:
            s = frozenset(pi[i] for i in ev_points(e))
            if s not in where or ev_sign(where[s]) != ev_sign(e):
                ok = False
                break
            imgs[e] = where[s]
        if not ok:
            continue
        # keeps "same setting pair"
        for e in EV1:
            for f in EV1:
                same = (e[0], e[2]) == (f[0], f[2])
                same_img = (imgs[e][0], imgs[e][2]) == (imgs[f][0], imgs[f][2])
                assert same == same_img
        keep.append(pi)
    assert len(keep) == 16
    sig = np.zeros((16, 9), dtype=np.int64)
    sgn = np.zeros((16, 9), dtype=np.int64)
    for g, pi in enumerate(keep):
        for i in range(9):
            v = B1[list(pi), i]
            hit = [(j, s) for j in range(9) for s in (1, -1) if np.array_equal(v, s * B1[:, j])]
            assert len(hit) == 1
            sig[g, i], sgn[g, i] = hit[0]
    return np.array(keep, dtype=np.int64), sig, sgn


G1, SIG, SGN = one_copy_group()


def coords(F):
    n = F.ndim
    c = F.astype(np.int64)
    for k in range(n):
        c = np.moveaxis(np.tensordot(c, B1, axes=([k], [0])), -1, k)
    return c


def back(c):
    n = c.ndim
    f = c
    for k in range(n):
        f = np.moveaxis(np.tensordot(f, B1.T, axes=([k], [0])), -1, k)
    return f


def is_spanned(F):
    return bool(set(np.unique(F).tolist()) <= {0, 1}) and np.array_equal(back(coords(F)), 16 ** F.ndim * F)


def block(ct, sgn=SGN):
    """rows c(g f) for g = (pi_1 .. pi_n) at one copy order ; ct = c(f) with its axes already permuted"""
    n = ct.ndim
    idx, sg = [], 1
    for k in range(n):
        shape = [1] * (2 * n)
        shape[k], shape[n + k] = 16, 9
        idx.append(SIG.reshape(shape))
        sg = sg * sgn.reshape(shape)
    return (ct[tuple(idx)] * sg).reshape(16 ** n, 9 ** n)


def gram(F, F2=None):
    """sum over G of c(g f) c(g f2)^T, exact, as int64"""
    n = F.ndim
    c, c2 = coords(F), None if F2 is None else coords(F2)
    out = np.zeros((9 ** n, 9 ** n))
    for tau in itertools.permutations(range(n)):
        X = block(c.transpose(tau)).astype(np.float64)
        Y = X if c2 is None else block(c2.transpose(tau)).astype(np.float64)
        out += X.T @ Y
    r = np.rint(out)
    assert np.array_equal(r, out) and np.abs(r).max() < 2 ** 53
    return r.astype(np.int64)


def literal_rows(F, tau, gs):
    """coordinates of the literal image sets: F' [pi_1 k_1, .., pi_n k_n] = F.transpose(tau)[k_1, .., k_n]"""
    n = F.ndim
    Ft = F.transpose(tau)
    rows = []
    for g in gs:
        Fi = np.zeros_like(Ft)
        Fi[np.ix_(*[G1[g[k]] for k in range(n)])] = Ft
        rows.append(coords(Fi).reshape(-1))
    return np.array(rows)


def check_group(F, group_size, sgn=SGN, sample=None):
    n = F.ndim
    size = 16 ** n * math.factorial(n)
    c = coords(F)
    all_g = list(itertools.product(range(16), repeat=n))
    bad = total = 0
    rng = random.Random(5801)
    for tau in itertools.permutations(range(n)):
        gs = all_g if sample is None else [all_g[rng.randrange(len(all_g))] for _ in range(sample // math.factorial(n) + 1)]
        rule = block(c.transpose(tau), sgn)
        lit = literal_rows(F, tau, gs)
        pos = [sum(g[k] * 16 ** (n - 1 - k) for k in range(n)) for g in gs]
        bad += int((rule[pos] != lit).any(axis=1).sum())
        total += len(gs)
    return size == group_size, bad, total


# ---------- polynomials in q, ascending coefficients, Fractions ----------
def padd(a, b):
    m = max(len(a), len(b))
    return [(a[i] if i < len(a) else Fr(0)) + (b[i] if i < len(b) else Fr(0)) for i in range(m)]


def pmul(a, b):
    out = [Fr(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return out


def peval(p, q):
    v = Fr(0)
    for c in reversed(p):
        v = v * q + c
    return v


def trim(p):
    p = list(p)
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    return p


def prem(a, b):
    a, b = trim(a), trim(b)
    while len(a) >= len(b) and any(a):
        f = a[-1] / b[-1]
        s = len(a) - len(b)
        for i, y in enumerate(b):
            a[s + i] -= f * y
        a.pop()
        if not a:
            a = [Fr(0)]
    return trim(a)


def sturm_roots(p, lo, hi):
    """distinct real roots of p in (lo, hi]"""
    seq = [trim(p), trim([i * c for i, c in enumerate(p)][1:])]
    while any(seq[-1]) and len(seq[-1]) > 1:
        r = prem(seq[-2], seq[-1])
        seq.append([-c for c in r])
    def var(x):
        s = [v for v in (peval(t, x) for t in seq) if v != 0]
        return sum(1 for u, v in zip(s, s[1:]) if (u > 0) != (v > 0))
    return var(lo) - var(hi)


def squarefree(p):
    """True if p and p' have no common root: the Euclidean sequence of (p, p') ends in a non-zero constant"""
    seq = [trim(p), trim([i * c for i, c in enumerate(p)][1:])]
    while any(seq[-1]) and len(seq[-1]) > 1:
        seq.append(prem(seq[-2], seq[-1]))
    return len(seq[-1]) == 1 and seq[-1][0] != 0


# ---------- loading ----------
def event_array(ev):
    """ev: list of per-copy (x, a, y, b)"""
    vecs = []
    for e in ev:
        v = np.zeros(16, dtype=np.int64)
        v[list(ev_points(tuple(e)))] = 1
        vecs.append(v)
    F = vecs[0]
    for v in vecs[1:]:
        F = np.multiply.outer(F, v)
    return F


def is_mid(e, f):
    return any((a[0], a[2]) == (b[0], b[2]) and tuple(a) != tuple(b) for a, b in zip(e, f))


def load(path):
    d = json.load(open(path))
    raw = d["sets"] if "sets" in d else d["sets_hex"]
    assert isinstance(raw, dict) and set(raw) == set(d["weights"])
    length = len(next(iter(raw.values())))
    n = {256: 2, 4096: 3}[length]
    sets = {}
    for k, s in raw.items():
        assert len(s) == length and set(s) <= {"0", "1"}
        v = np.array([int(ch) for ch in s], dtype=np.int64)
        if n == 2:      # L order -> 16 k1 + k2
            F = np.zeros((16, 16), dtype=np.int64)
            for k1 in range(16):
                for k2 in range(16):
                    F[k1, k2] = v[(k1 >> 2) + 4 * (k1 & 3) + 16 * ((k2 >> 2) + 4 * (k2 & 3))]
        else:
            F = v.reshape(16, 16, 16)
        sets[k] = F
    pairs = None
    if "pairs" in d:
        pairs = []
        for pr in d["pairs"]:
            if n == 2:
                def conv(e):
                    x, a = e
                    return ((x[0], a[0], x[1], a[1]), (x[2], a[2], x[3], a[3]))
                e, e2, mu = conv(pr["e"]), conv(pr["e2"]), Fr(pr["lambda"])
            else:
                e, e2 = tuple(tuple(t) for t in pr["e"]), tuple(tuple(t) for t in pr["e2"])
                mu = Fr(pr["lambda"]) / 4096
            pairs.append((e, e2, mu))
    if "root_lo" in d:
        lo, hi = Fr(d["root_lo"]), Fr(d["root_hi"])
    else:
        lo, hi = Fr(d["root_float"]) - Fr(1, 10 ** 9), Fr(d["root_float"]) + Fr(1, 10 ** 9)
    return {"n": n, "sets": sets, "weights": {k: Fr(v) for k, v in d["weights"].items()}, "pairs": pairs,
            "poly": [Fr(c) for c in d["poly"]], "group_size": d["group_size"], "lo": lo, "hi": hi}


def lcm_den(vals):
    den = 1
    for v in vals:
        den = den * v.denominator // math.gcd(den, v.denominator)
    return den


def one_copy_tables():
    comp = np.zeros((9, 9), dtype=bool)
    coef = np.zeros((9, 9), dtype=np.int64)
    deg = np.zeros((9, 9), dtype=np.int64)
    for i, m in enumerate(MON1):
        for j, mm in enumerate(MON1):
            ok = all(m[w] is None or mm[w] is None or m[w] == mm[w] for w in (0, 1))
            comp[i, j] = ok
            if not ok:
                continue
            prod = tuple(None if m[w] == mm[w] else (m[w] if m[w] is not None else mm[w]) for w in (0, 1))
            if prod == (None, None):
                coef[i, j], deg[i, j] = 1, 0
            elif prod[0] is None or prod[1] is None:
                coef[i, j], deg[i, j] = 0, 0
            else:
                coef[i, j], deg[i, j] = (-1) ** (prod[0] * prod[1]), 1
    return comp, coef, deg


def n_copy_tables(n):
    comp1, coef1, deg1 = one_copy_tables()
    comp, coef, deg = comp1, coef1, deg1
    for _ in range(n - 1):
        m = comp.shape[0]
        comp = (comp[:, None, :, None] & comp1[None, :, None, :]).reshape(m * 9, m * 9)
        coef = (coef[:, None, :, None] * coef1[None, :, None, :]).reshape(m * 9, m * 9)
        deg = (deg[:, None, :, None] + deg1[None, :, None, :]).reshape(m * 9, m * 9)
    return comp, coef, deg


def event_prob(e):
    p = [Fr(1)]
    for c in e:
        p = pmul(p, [Fr(1, 4), Fr(ev_sign(tuple(c)), 4)])
    return p


def verify(C, kind, grams, tag, weights=None, sets=None, pairs=None):
    """returns dict of results ; raises nothing on a failed check"""
    n = C["n"]
    N, size = 16 ** n, 16 ** n * math.factorial(n)
    weights = C["weights"] if weights is None else weights
    sets = C["sets"] if sets is None else sets
    pairs = C["pairs"] if pairs is None else pairs
    res = {}
    res["S"] = all(is_spanned(F) for F in sets.values()) and all(y > 0 for y in weights.values())
    if not res["S"]:
        return res
    Dy = lcm_den(weights.values())
    M = np.zeros((9 ** n, 9 ** n), dtype=object)
    for k, F in sets.items():
        if k not in grams:
            grams[k] = gram(F)
        M = M + grams[k].astype(object) * int(weights[k] * Dy)
    if kind == "nb":
        comp, coef, deg = n_copy_tables(n)
        res["I_nonzero"] = int(sum(1 for v in M[~comp] if v != 0))
        p = []
        for dd in range(n + 1):
            mask = comp & (deg == dd)
            p.append(Fr(int(sum(int(a) * int(b) for a, b in zip(M[mask], coef[mask]))), Dy * size * N * N))
    else:
        res["E"] = all(e == e2 or is_mid(e, e2) for (e, e2, mu) in pairs)
        if not res["E"]:
            return res
        Dl = lcm_den([mu for (_, _, mu) in pairs])
        R = np.zeros((9 ** n, 9 ** n), dtype=object)
        p = [Fr(0)]
        for (e, e2, mu) in pairs:
            key = ("pair", e, e2)
            if key not in grams:
                X = gram(event_array(e), event_array(e2))
                grams[key] = X + X.T
            R = R + grams[key].astype(object) * int(mu * Dl)
            if e == e2:
                p = padd(p, [mu * c for c in event_prob(e)])
        diff = 2 * Dl * M - Dy * R
        res["I_nonzero"] = int(sum(1 for v in diff.ravel() if v != 0))
    p = (p + [Fr(0)] * (n + 1))[:n + 1]
    res["poly"] = p
    res["P_equal"] = p == (C["poly"] + [Fr(0)] * (n + 1))[:n + 1]
    return res


def run_file(kind, path):
    t0 = time.time()
    C = load(path)
    n = C["n"]
    tag = "%s n=%d" % (kind, n)
    print("FILE %s (%s): %d sets, %s pairs, weights sum %s" % (path, tag, len(C["sets"]),
          "no" if C["pairs"] is None else len(C["pairs"]), sum(C["weights"].values())), flush=True)
    assert (kind == "mid") == (C["pairs"] is not None)
    k0 = sorted(C["sets"], key=int)[0]
    ok_size, bad, total = check_group(C["sets"][k0], C["group_size"], sample=None if n == 2 else 2000)
    print("  G  |G| = %d equals the file's group_size: %s ; coordinate rule against literal images: %d of %d "
          "differ" % (16 ** n * math.factorial(n), ok_size, bad, total), flush=True)
    g_ok = ok_size and bad == 0
    k4 = None
    if n == 2:
        sgn_bad = SGN.copy()
        sgn_bad[3, 4] *= -1
        _, badk, totk = check_group(C["sets"][k0], C["group_size"], sgn=sgn_bad)
        k4 = badk > 0
        print("  K4 one sign of one relabelling changed: %d of %d rows differ (must be > 0)" % (badk, totk), flush=True)
    grams = {}
    res = verify(C, kind, grams, tag)
    p = res["poly"]
    print("  S  sets spanned and weights > 0: %s" % res["S"], flush=True)
    if kind == "mid":
        print("  E  every stored pair is equal or a middle pair: %s" % res["E"], flush=True)
    print("  I  non-zero entries (%s): %d" % ("incompatible pairs of Chat" if kind == "nb" else
          "left side minus right side", res["I_nonzero"]), flush=True)
    den = lcm_den(p)
    print("  P  p(q) = (%s) / %d ; equals the file's poly: %s" % (" , ".join(str(int(c * den)) for c in p), den,
          res["P_equal"]), flush=True)
    s_lo, s_hi, s_1 = peval(p, C["lo"]), peval(p, C["hi"]), peval(p, Fr(1))
    roots = sturm_roots(p, C["hi"], Fr(1))
    s_m1, s_0 = peval(p, Fr(-1)), peval(p, Fr(0))
    sqf = squarefree(p)
    roots_all = sturm_roots(p, Fr(-1), Fr(1))
    step = Fr(1, 10 ** 6)
    bound = Fr(math.ceil(4 * C["hi"] / step)) * step
    b_ok = peval(p, bound / 4) < 0 and peval(p, (bound - step) / 4) > 0
    p_ok = (res["P_equal"] and s_lo > 0 and s_hi < 0 and s_1 < 0 and roots == 0 and s_m1 > 0 and s_0 > 0 and sqf
            and roots_all == 1 and b_ok)
    print("     bracket %.12f .. %.12f ; p(lower) > 0: %s ; p(upper) < 0: %s ; p(1) < 0: %s ; roots "
          "in (upper, 1]: %d" % (float(C["lo"]), float(C["hi"]), s_lo > 0, s_hi < 0, s_1 < 0, roots), flush=True)
    print("     p(-1) > 0: %s ; p(0) > 0: %s ; square-free: %s ; roots in (-1, 1]: %d" % (s_m1 > 0, s_0 > 0, sqf,
                                                                                         roots_all), flush=True)
    print("     CHSH bound, 4 q* rounded UP to six decimals: %d.%06d ; p(bound / 4) < 0 and p((bound - 1e-6) / 4) "
          "> 0: %s" % (bound.numerator * 10 ** 6 // bound.denominator // 10 ** 6,
                       bound.numerator * 10 ** 6 // bound.denominator % 10 ** 6, b_ok), flush=True)
    main_ok = g_ok and res["S"] and res.get("E", True) and res["I_nonzero"] == 0 and p_ok
    # controls
    w1 = dict(C["weights"])
    w1[k0] += Fr(1, 1000)
    r1 = verify(C, kind, grams, tag, weights=w1)
    k1 = r1["I_nonzero"] > 0 or not r1["P_equal"]
    print("  K1 first weight + 1/1000: non-zero entries %d ; poly equal %s (one must fail): control %s"
          % (r1["I_nonzero"], r1["P_equal"], "FAILS AS IT MUST" if k1 else "DID NOT FAIL"), flush=True)
    s2 = dict(C["sets"])
    F = s2[k0].copy()
    F.flat[0] ^= 1
    s2[k0] = F
    r2 = verify(C, kind, {}, tag, sets=s2)
    k2 = not r2["S"]
    print("  K2 first set with point 0 flipped: spanned check %s: control %s" % (r2["S"],
          "FAILS AS IT MUST" if k2 else "DID NOT FAIL"), flush=True)
    k3 = True
    if kind == "mid":
        j = next(i for i, (e, e2, mu) in enumerate(C["pairs"]) if e != e2)
        e, e2, mu = C["pairs"][j]
        e_bad = tuple((1 - c[0], c[1], 1 - c[2], c[3]) for c in e)
        pb = list(C["pairs"])
        pb[j] = (e, e_bad, mu)
        r3 = verify(C, kind, grams, tag, pairs=pb)
        k3 = not r3["E"]
        print("  K3 one pair replaced by a pair outside the pair set: pair check %s: control %s" % (r3["E"],
              "FAILS AS IT MUST" if k3 else "DID NOT FAIL"), flush=True)
    ok = main_ok and k1 and k2 and k3 and (k4 is None or k4)
    print("  RESULT %s: %s ; %.0f s" % (tag, "VERIFIED, controls fail as they must" if ok else "NOT VERIFIED",
                                        time.time() - t0), flush=True)
    return ok


def main():
    args = sys.argv[1:]
    assert len(args) % 2 == 0 and args
    print("one-copy relabellings: 128 ; fixing PR_q: %d" % len(G1), flush=True)
    oks = [run_file(args[i], args[i + 1]) for i in range(0, len(args), 2)]
    print("AS STATED BEFORE THE RUN: %s" % all(oks))
    print("DONE-VERIFY")


if __name__ == "__main__":
    main()
