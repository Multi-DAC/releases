"""Letter 75, Step 2 consequences, on my transfer matrix (check_letter75_h4.transfer):
(a) the spectrum of T is closed under lambda -> 1/lambda (her Q_t T Q_t^H = T^-1), with exactly L5 eigenvalues < 1;
(b) sigma_min of V_{-,+} (u+ rows of the L5 eigenvectors with lambda < 1) as a5 crosses 1/m, m = 0.7,
    minimised over a k1 grid of 401 points in [-pi, pi], for DW and PV, L5 = 1..4.
Expect (b) bounded away from 0 for a5 < 1/m; at k1 = 0 it should vanish for a5 > 1/m if the edge mode is there."""
import sys
import numpy as np

sys.path.insert(0, ".")
from check_letter75_h4 import transfer  # noqa: E402

m = 0.7
ks = np.linspace(-np.pi, np.pi, 401)
print("(a) inversion symmetry of spec(T) and count below 1")
worst, badcount = 0.0, 0
for mf in (0.0, 1.0):
    for L5 in (1, 2, 3, 4):
        for k1 in ks[::20]:
            for a5 in (0.3, 1.0, 1.4, 1.6, 2.0):
                T = transfer(L5, k1, m, a5, mf)[0]
                ev = np.sort(np.linalg.eigvalsh((T + T.conj().T) / 2))
                worst = max(worst, np.abs(ev * ev[::-1] - 1).max())
                badcount += int((ev < 1).sum() != L5)
print("   max |lambda_j lambda_(2L5+1-j) - 1| = %.2e ; cases without exactly L5 below 1: %d" % (worst, badcount))
print("(b) min over k1 of sigma_min(V_{-,+}); k1 of the minimum")
for name, mf in (("DW", 0.0), ("PV", 1.0)):
    for L5 in (1, 2, 3, 4):
        row = []
        for a5 in (0.3, 1.0, 1.3, 1.4, 1.42, 1.44, 1.5, 1.6, 2.0):
            best = (9.0, 0.0)
            for k1 in ks:
                T = transfer(L5, k1, m, a5, mf)[0]
                ev, V = np.linalg.eigh((T + T.conj().T) / 2)
                Vm = V[:, :L5]
                s = np.linalg.svd(Vm[:L5, :], compute_uv=False)[-1]
                if s < best[0]:
                    best = (s, k1)
            row.append("%.2f:%.1e@%+.2f" % (a5, best[0], best[1]))
        print("   %s L5=%d  %s" % (name, L5, "  ".join(row)))
