"""Cold audit of Wren's PROOF-T1-2-L4-and-overlap.md, Propositions 2, 3, 4 (Lemmas 4A-4C), with MY operator.

D(z) = B0 + z B+ + B-/z from edge_modes_dw.op's row blocks (my conventions); Wren's a5-scaled D_X is a5 times this,
so her a5 R_j is the residue R of H(z)^-1 with H(z) := gamma0 D(z) (mine), and her a5 G(tau) is mine.
Poles: finite nonzero eigenvalues of the quadratic pencil z^2 B+ + z B0 + B- (companion form), cross-checked against
transfer() from check_letter75_h4. Residues: R = V (V^H H'(lam) V)^-1 V^H on ker H(lam), clustered for degeneracy;
cross-checked by a unit-circle contour where it converges.

Tests, per (m, a5, L5, k1):
 P2  number of DW poles in (r_*(k1), 1) <= 1; max PV pole <= r_*(k1); interlacing x_{j-1}^PV <= x_j^DW <= x_{j+1}^PV;
     eigenvalues of gamma0 (a5 (D_DW - D_PV)) (her claim +-1/2).
 P3  P0- (sum_all R) P0- = P0- ; P0+ (sum_all R / mu^2) P0+ = P0+ ; spec(sum_{|lam|<1} R) in [0, 2].
 4B  det D_DW(i) / det D_PV(i) against her closed form; > 1/8.
 4C  min over |z| = r of |det M2(z)| / eps, with r, eps chosen as in her Prop 4.
 P4  ||a5 Delta(t + 1)||_1 / ((t + 2)^4 r_*^t), and / r_*^t, for t = 0 .. 80.
Controls: a5 above a5max(m) (pointwise Lemma 4(i) fails there) to see which claims fail with the hypothesis."""
import sys
import numpy as np
import scipy.linalg as sla

sys.path.insert(0, ".")
from check_letter75_h4 import blocks, transfer  # noqa: E402
from edge_modes_dw import g0  # noqa: E402


TIE = 1e-12


def l5scan():
    """Uniformity in L5 of ||a5 Delta(t + 1)||_1 / r_*^t (max over t <= 80), L5 up to 128; and Prop 3's identity
    residual at the same points (to see where the 2e-7 comes from)."""
    for m, a5, k1 in ((0.7, 0.2, 0.0), (0.7, a5max(0.7), 0.0), (0.3, 1.0, 0.0), (0.3, 0.05, 0.0),
                      (0.9, a5max(0.9), np.pi / 2), (0.7, 0.5, np.pi)):
        line = []
        for L5 in (8, 16, 32, 64, 128):
            r = case(m, a5, L5, k1, tmax=80, do_contour=False)
            line.append("L5=%d: %.3f (n_surf %d, idp %.0e)" % (L5, max(r["r0"]), r["n_surf"],
                                                               max(r["res"]["DW"]["idp"], r["res"]["PV"]["idp"])))
        print("  max_t ||a5 Delta||_1 / r_*^t, m=%.2f a5=%.4f k1=%.2f: %s" % (m, a5, k1, "; ".join(line)))
        sys.stdout.flush()


def a5max(m):
    return 2 * (1 - m) / (m * (2 - m))


def rstar(k1, m):
    c = 2 - m - np.cos(k1)
    ch = (c * c + 1 + np.sin(k1) ** 2) / (2 * c)
    return np.exp(-np.arccosh(ch))


def Dz(B, z):
    B0, Bp, Bm = B
    return B0 + z * Bp + Bm / z


def poles(B):
    B0, Bp, Bm = B
    n = B0.shape[0]
    Z, I = np.zeros((n, n)), np.eye(n)
    A = np.block([[Z, I], [-Bm, -B0]])
    Bmat = np.block([[I, Z], [Z, Bp]])
    w = sla.eigvals(A, Bmat)
    w = w[np.isfinite(w)]
    w = w[np.abs(w) > 1e-9]
    return np.sort(w.real), np.abs(w.imag).max() if len(w) else 0.0


def cluster(vals, tol=1e-6):
    groups = []
    for v in np.sort(vals):
        if groups and abs(v - groups[-1][-1]) <= tol * max(1.0, abs(v)):
            groups[-1].append(v)
        else:
            groups.append([v])
    return [(np.mean(g), len(g)) for g in groups]


def residues(B, lams, G0):
    out = []
    B0, Bp, Bm = B
    herm = 0.0
    for lam, d in cluster(lams):
        H = G0 @ Dz(B, lam)
        herm = max(herm, np.abs(H - H.conj().T).max())
        H = (H + H.conj().T) / 2
        w, U = np.linalg.eigh(H)
        idx = np.argsort(np.abs(w))[:d]
        V = U[:, idx]
        Hp = G0 @ (Bp - Bm / lam ** 2)
        Hp = (Hp + Hp.conj().T) / 2
        R = V @ np.linalg.inv(V.conj().T @ Hp @ V) @ V.conj().T
        out.append((lam, d, (R + R.conj().T) / 2, np.sort(np.abs(w))[d] if d < len(w) else np.inf))
    return out, herm


def contour_sum(B, G0, N):
    th = 2 * np.pi * (np.arange(N) + 0.5) / N
    z = np.exp(1j * th)
    B0, Bp, Bm = B
    S = np.zeros_like(B0)
    for ch in np.array_split(np.arange(N), max(1, N // 2048)):
        Dk = B0[None] + z[ch, None, None] * Bp[None] + Bm[None] / z[ch, None, None]
        Hk = G0[None] @ Dk
        S += (np.linalg.inv(Hk) * z[ch, None, None]).sum(axis=0)
    return S / N


def tracenorm(A):
    return np.linalg.svd(A, compute_uv=False).sum()


def case(m, a5, L5, k1, tmax=80, do_contour=True):
    G0 = np.kron(np.eye(L5), g0)
    P0p = np.kron(np.eye(L5), 0.5 * (np.eye(2) + g0))
    P0m = np.kron(np.eye(L5), 0.5 * (np.eye(2) - g0))
    rs = rstar(k1, m)
    res = {}
    data = {}
    for X, mf in (("DW", 0.0), ("PV", 1.0)):
        B = blocks(L5, k1, m, a5, mf)
        lam, im = poles(B)
        Tm = transfer(L5, k1, m, a5, mf)[0]
        lt = np.sort(np.linalg.eigvals(Tm).real)
        cross = np.abs(np.sort(lam) - lt).max() / max(1, np.abs(lt).max()) if len(lt) == len(lam) else np.inf
        inside = np.sort(lam)[:L5]   # poles pair as (lam, 1/lam): the L5 smallest are inside (robust at 1 - 1e-16)
        n_strict = int(np.sum(lam < 1))
        R1, h1 = residues(B, inside, G0)                 # cluster inside and outside separately, never across 1
        R2, h2 = residues(B, np.sort(lam)[L5:], G0)
        R, herm = R1 + R2, max(h1, h2)
        Sall = sum(r[2] for r in R)
        Smu2 = sum(r[2] / r[0] ** 2 for r in R)
        Sin = sum(r[2] for r in R if r[0] <= inside.max())
        idm = np.abs(P0m @ Sall @ P0m - P0m).max()
        idp = np.abs(P0p @ Smu2 @ P0p - P0p).max()
        ev = np.linalg.eigvalsh(Sin)
        cont = np.nan
        if do_contour and inside.max() ** 4096 < 1e-13:
            Sc = contour_sum(B, G0, 4096)
            cont = np.abs(Sc - Sin).max()
        data[X] = dict(B=B, lam=lam, inside=np.sort(inside)[::-1], R=R)
        res[X] = dict(n=len(lam), n_in=n_strict, im=im, cross=cross, herm=herm, idm=idm, idp=idp,
                      evmin=ev.min(), evmax=ev.max(), cont=cont, gap=min(r[3] for r in R))
    # P2
    dw, pv = data["DW"]["inside"], data["PV"]["inside"]
    n_surf = int(np.sum(dw > rs + TIE))   # TIE: at a5 = a5max(m), k1 = 0 every bulk pole equals r_* exactly
    pvmax = pv.max()
    xdw = np.sort(dw + 1 / dw)
    xpv = np.sort(pv + 1 / pv)
    il = 0.0
    for j in range(L5):
        lo = xpv[j - 1] if j >= 1 else 2.0
        hi = xpv[j + 1] if j + 1 < L5 else np.inf
        il = max(il, lo - xdw[j], xdw[j] - hi)
    E = data["DW"]["B"][0] - data["PV"]["B"][0]
    e2 = np.linalg.eigvalsh(a5 * G0 @ (data["DW"]["B"][0] - data["PV"]["B"][0]) + 0j * E)
    e2 = e2[np.abs(e2) > 1e-12]
    # 4B
    dD = np.linalg.slogdet(Dz(data["DW"]["B"], 1j))
    dP = np.linalg.slogdet(Dz(data["PV"]["B"], 1j))
    M2i = dD[0] / dP[0] * np.exp(dD[1] - dP[1])
    c = 2 - m - np.cos(k1)
    detX = c * c + 1 + np.sin(k1) ** 2
    tr = 2 + a5 ** 2 * detX / (1 + a5 * c)
    mu = (tr + np.sqrt(tr * tr - 4)) / 2
    Tpp = 1 / (1 + a5 * c)
    closed = ((mu - Tpp) + mu ** (-2 * L5) * (Tpp - 1 / mu)) / ((mu - 1 / mu) * (1 + mu ** (-L5)) ** 2)
    # P4 and 4C
    ratios4, ratios0, ratios4c = [], [], []
    lam1 = dw[0] if n_surf else None
    for t in range(tmax + 1):
        A = np.zeros((2 * L5, 2 * L5), complex)
        for lam, d, Rm, _ in data["DW"]["R"]:
            if lam <= rs + TIE:
                A += lam ** t * Rm @ G0
        for lam, d, Rm, _ in data["PV"]["R"]:
            if lam <= data["PV"]["inside"].max():
                A -= lam ** t * Rm @ G0
        tn = tracenorm(A)
        ratios4.append(tn / ((t + 2) ** 4 * rs ** t))
        ratios0.append(tn / rs ** t)
        if t in (0, 5, 20, 80):
            dlt = (1 - rs) / (4 * (t + 2))
            r = rs + dlt
            if lam1 is not None and abs(lam1 - r) < dlt / 3:
                r = rs + 2 * dlt
            eps = dlt / 3
            th = 2 * np.pi * np.arange(256) / 256
            vals = []
            for zz in r * np.exp(1j * th):
                a = np.linalg.slogdet(Dz(data["DW"]["B"], zz))
                b = np.linalg.slogdet(Dz(data["PV"]["B"], zz))
                vals.append(np.exp(a[1] - b[1]))
            ratios4c.append(min(vals) / eps)
    return dict(res=res, rs=rs, n_surf=n_surf, lam1=lam1, pvmax=pvmax, il=il, e2=e2, M2i=M2i, closed=closed,
                r4=max(ratios4), r0=ratios0, r4c=min(ratios4c))


def main():
    big = dict(idm=0, idp=0, evmin=1, evmax=0, cross=0, im=0, herm=0, cont=0)
    worst = {}
    rows = []
    ms = [0.3, 0.7, 0.9]
    for m in ms:
        amax = a5max(m)
        a5s = sorted(set([0.05, 0.2, min(0.5, amax), amax] + ([1.0] if amax > 1 else [])))
        ctrl = [1.0, 1.3] if amax < 1 else [amax * 1.2]
        for a5 in a5s + ctrl:
            is_ctrl = a5 > amax + 1e-12
            for L5 in (1, 2, 3, 4, 6, 8, 12, 16):
                for k1 in np.linspace(0, np.pi, 7):
                    r = case(m, a5, L5, k1, do_contour=(L5 <= 6))
                    tag = "CTRL" if is_ctrl else "in"
                    rows.append((m, a5, L5, k1, tag, r))
                    if not is_ctrl:
                        for X in ("DW", "PV"):
                            q = r["res"][X]
                            big["idm"] = max(big["idm"], q["idm"]); big["idp"] = max(big["idp"], q["idp"])
                            big["evmin"] = min(big["evmin"], q["evmin"]); big["evmax"] = max(big["evmax"], q["evmax"])
                            big["cross"] = max(big["cross"], q["cross"]); big["im"] = max(big["im"], q["im"])
                            big["herm"] = max(big["herm"], q["herm"])
                            if not np.isnan(q["cont"]):
                                big["cont"] = max(big["cont"], q["cont"])
                            assert q["n"] == 2 * L5 and q["n_in"] == L5, (m, a5, L5, k1, X, q)
    print("in-range cases (a5 <= a5max(m)), m in", ms, ", L5 in 1..16, k1 in [0, pi] (7):")
    print("  pole count 2 L5 with L5 inside, every case; max |Im pole| %.1e; pencil vs transfer() %.1e; "
          "max |H - H^H| at poles %.1e" % (big["im"], big["cross"], big["herm"]))
    print("  P3 identities: max |P0- Sall P0- - P0-| %.1e, max |P0+ S/mu^2 P0+ - P0+| %.1e" % (big["idm"], big["idp"]))
    print("  P3 spec(sum_{|lam|<1} R) in [%.6f, %.6f]  (claim [0, 2]);  contour vs residues max %.1e"
          % (big["evmin"], big["evmax"], big["cont"]))
    for tag in ("in", "CTRL"):
        sel = [x for x in rows if x[4] == tag]
        nsurf = max(x[5]["n_surf"] for x in sel)
        pvx = max(x[5]["pvmax"] - x[5]["rs"] for x in sel)
        il = max(x[5]["il"] for x in sel)
        m2 = min(x[5]["M2i"] for x in sel)
        m2err = max(abs(x[5]["M2i"] - x[5]["closed"]) / abs(x[5]["closed"]) for x in sel)
        r4 = max(x[5]["r4"] for x in sel)
        r4c = min(x[5]["r4c"] for x in sel)
        e2 = sorted(set(np.round(np.concatenate([x[5]["e2"] for x in sel]), 9)))
        print("%s: max # DW poles in (r_*, 1) = %d; max (PV pole - r_*) = %.3e; max interlacing violation %.2e"
              % (tag, nsurf, pvx, il))
        print("    det M2(i): min %.6f (claim > 1/8), closed-form rel err max %.1e; nonzero eigs of gamma0 a5(D_DW-D_PV): %s"
              % (m2, m2err, e2))
        print("    P4 max_t ||a5 Delta(t+1)||_1 / ((t+2)^4 r_*^t) = %.3e; 4C min_t min_|z|=r |det M2| / eps = %.3e"
              % (r4, r4c))
        if tag == "CTRL":
            bad = [x for x in sel if x[5]["n_surf"] > 1 or x[5]["il"] > 1e-9]
            for x in bad[:8]:
                print("    CTRL failure m=%.2f a5=%.3f L5=%d k1=%.3f: n_surf %d, PV-r* %.3e, interlace %.2e"
                      % (x[0], x[1], x[2], x[3], x[5]["n_surf"], x[5]["pvmax"] - x[5]["rs"], x[5]["il"]))
    # growth of ||a5 Delta||_1 / r_*^t in t: worst over in-range cases at t = 0, 10, 20, 40, 80
    sel = [x for x in rows if x[4] == "in"]
    for t in (0, 10, 20, 40, 80):
        w = max(sel, key=lambda x: x[5]["r0"][t])
        print("  t=%2d: max ||a5 Delta(t+1)||_1 / r_*^t = %.3e at m=%.2f a5=%.3f L5=%d k1=%.3f"
              % (t, w[5]["r0"][t], w[0], w[1], w[2], w[3]))
    # L5 dependence at fixed (m, a5, k1): does the ratio grow with L5?
    for m, a5, k1 in ((0.7, 0.2, 0.0), (0.3, 1.0, 0.0), (0.9, a5max(0.9), np.pi / 2)):
        line = []
        for L5 in (1, 2, 4, 8, 16):
            x = [y for y in sel if abs(y[0] - m) < 1e-9 and abs(y[1] - a5) < 1e-9 and y[2] == L5
                 and abs(y[3] - k1) < 1e-9]
            if x:
                line.append("L5=%d: %.2e" % (L5, max(x[0][5]["r0"])))
        print("  max_t ||a5 Delta||_1 / r_*^t at m=%.2f a5=%.3f k1=%.2f: %s" % (m, a5, k1, "; ".join(line)))


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "l5scan":
        l5scan()
    else:
        main()
