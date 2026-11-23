"""Independent check (Clawd) of Proposition 11, float: smallest eigenvalue of K_np on W_2 in the orthonormal Walsh basis, against -8/65536 (Reduction B)."""
import json
import os

import numpy as np

import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import knp_check as k  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
d = json.load(open(os.path.join(HERE, "data_208a.json")))
q = int(d["q_t"])
B = np.array(k.B, dtype=float) / 4.0
Dt = np.array([[int(v) / q for v in row] for row in d["Z_t"]])
Z = np.array(k.read_Z(os.path.join(HERE, "..", "documentation", "paper.tex")), dtype=float)
H = np.zeros((9, 9))
H[6, 7] = H[7, 6] = 1
H[5, 8] = H[8, 5] = -1
Qt = B.T @ Dt @ B
QD = Z / (4096 * 16)   # (B/4)^T (B Z B^T / (256*4096)) (B/4) = 16 Z 16 / (16 * 256 * 4096)
QH = H / 16
K = np.kron(Qt, Qt) + np.kron(QD, QH)
ev = np.linalg.eigvalsh(K)
print("Q_t eigenvalues min %.4e max %.4e ; Q_Delta min %.4e max %.4e" % (np.linalg.eigvalsh(Qt).min(), np.linalg.eigvalsh(Qt).max(),
                                                                       np.linalg.eigvalsh(QD).min(), np.linalg.eigvalsh(QD).max()))
print("K_np on W_2: smallest eigenvalue %.6e ; -8/65536 = %.6e ; -6/65536 = %.6e ; margin to lambda = 8: %.3e"
      % (ev.min(), -8 / 65536, -6 / 65536, ev.min() + 8 / 65536))
