"""My own check of refuter row 1 (Adams' bound at eps = 0.267) and a scan of the ceiling over m.
Written independently of the refuter's cf_check.py: H_T free eigenvalue from the lemma's 8 formulas."""
import numpy as np

c1 = 2 + np.sqrt(2)  # (2 + sqrt2) d(d-1)/2 at d = 2


def adams(eps, m):
    return np.sqrt(1 - c1 * eps) - abs(1 - m)


def cf_T(m, N=1024):
    k = (np.arange(N) + 0.5) * 2 * np.pi / N - np.pi
    K1, K2 = np.meshgrid(k, k)
    U = (1 - np.cos(K1)) + (1 - np.cos(K2))
    b = U - m
    s2 = np.sin(K1) ** 2 + np.sin(K2) ** 2
    n = (2 + b) ** 2 + s2
    bT = (b * (2 + b) + s2) / n
    wT = np.sqrt(bT ** 2 + 4 * s2 / n ** 2)
    lam = (1 + bT / wT) / 2
    return -np.mean(np.log(lam))


m = 0.7
es = m * (2 - m) / c1
print('m=0.7: adams(0.267) =', adams(0.267, m), ' eps* =', es, ' adams(eps*) =', adams(es, m))
cf = cf_T(m, 2048)
print('c_f_T(0.7) N=2048:', cf, ' el0(0.267) =', 0.267 / np.sqrt(2 * cf), ' el0(eps*) =', es / np.sqrt(2 * cf))
print('scan: m, eps*(m), c_f_T(m), el0(m) = eps*/sqrt(2 c_f)')
best = None
for mm in np.arange(0.30, 1.0001, 0.05):
    e = mm * (2 - mm) / c1
    c = cf_T(mm)
    el = e / np.sqrt(2 * c)
    print(round(mm, 2), round(e, 5), round(c, 6), round(el, 5))
    if best is None or el > best[1]:
        best = (round(mm, 2), el)
print('best on grid:', best)
