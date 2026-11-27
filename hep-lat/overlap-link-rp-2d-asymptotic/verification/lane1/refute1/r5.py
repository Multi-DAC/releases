import os  # release copy
# Probe of Section 6's imported admissibility eps = 0.267: uniform-flux fields, min |eig H_W| vs |1 - U_p|
import numpy as np, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))  # release copy
from te_finite_beta import build_DW_links, links_from_plaquettes, s3
L = 20; V = L*L; G5 = np.kron(np.eye(V), s3)
for Q in (1, 4, 8, 12, 16, 18, 20, 24, 28, 32):
    th = np.full((L, L), 2*np.pi*Q/V)
    U0, U1 = links_from_plaquettes(th, Q, (0.3, 1.1))
    w = np.linalg.eigvalsh(G5 @ build_DW_links(L, U0, U1))
    print(f"Q={Q:2d} theta_p={2*np.pi*Q/V:.4f} |1-U_p|={2*np.sin(np.pi*Q/V):.4f}  min|eig H_W| = {np.abs(w).min():.4f}", flush=True)
