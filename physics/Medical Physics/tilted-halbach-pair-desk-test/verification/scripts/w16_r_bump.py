"""S3 and S4 with the pocket circle enlarged by the same amount in both (Wren's suggestion, letter 11): centre field,
ratio, smallest gap between cubes. z = +-12.3 mm kept, and also z scaled with R (12.3/29.3 R).
Usage: python -I -B w16_r_bump.py
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import w16_prism as wp
from w16_wren_checks import build, gap, bx, S, N  # noqa: E402  (prints their tables on import)

print()
print("--- R enlarged identically in S3 and S4 ---")
for R in (0.0293, 0.0303, 0.0313):
    for Z, tag in ((0.0123, "z 12.3 mm"), (0.0123 / 0.0293 * R, "z scaled")):
        s3 = build(S, N, R, (Z, -Z), "halbach"); s4 = build(S, N, R, (Z, -Z), "max")
        print("R %.1f mm, %s (%.2f mm): S3 %.2f mT, S4 %.2f mT, ratio %.4f; gap S3 %.3f mm, S4 %.3f mm"
              % (R * 1e3, tag, Z * 1e3, bx(s3) * 1e3, bx(s4) * 1e3, bx(s4) / bx(s3), gap(s3)[0] * 1e3, gap(s4)[0] * 1e3))
