"""Direct lattice sum for Route (v) rev 3's C_geom (line 201): fix a plaquette p = (s_p, y_p); sum over sites x and
spatial links l = ((t_l, y_l), (t_l, y_l + 1)) with y_l = y_p and p strictly between t_x and t_l (t_x <= s_p < t_l or
t_l <= s_p < t_x), of |t_l - t_x| exp(-2 g d(x, l)), d(x, l) = |t_l - t_x| + min(|x1 - y_l|, |x1 - y_l - 1|).
Compare with 4 q (1 + q)/(1 - q)^4, q = exp(-2 g). Written independently of Wren's verifier."""
import math

def direct(g, R):
    s_p, y_p = 0, 0
    tot = 0.0
    for tx in range(-R, R + 1):
        for tl in range(-R, R + 1):
            between = (tx <= s_p < tl) or (tl <= s_p < tx)
            if not between:
                continue
            n = abs(tl - tx)
            for x1 in range(-R, R + 1):
                k = min(abs(x1 - y_p), abs(x1 - y_p - 1))
                tot += n * math.exp(-2 * g * (n + k))
    return tot

for g in (0.2, 0.5, 1.0):
    q = math.exp(-2 * g)
    closed = 4 * q * (1 + q) / (1 - q) ** 4
    R = int(40 / g)
    d = direct(g, R)
    print('gamma', g, 'R', R, 'direct', d, 'closed', closed, 'rel diff', d / closed - 1)
