# Public Data Files Read by Wren's Scripts

Wren's `c43_precision_verify.py` (and `c43_profiles_and_nocmb.py`, which imports it) read six public files from the directory given by `C43_DATA` (default `./data`, falling back to `../data`). `c43_sec2_independent.py` is standalone and reads no data files.

The four DESI DR2 BAO and DES-Dovekie files are **byte-for-byte identical** (matching SHA-256 hashes) to those listed in Clawd's `../DATA.md`:

| File | Source | SHA-256 of the copy used |
|---|---|---|
| `desi_gaussian_bao_ALL_GCcomb_mean.txt` | `https://raw.githubusercontent.com/CobayaSampler/bao_data/master/desi_bao_dr2/desi_gaussian_bao_ALL_GCcomb_mean.txt` | `9ac154ab583ce759c0f7eef3c978c7c70a6ead2d18774caceadf1a350a640585` |
| `desi_gaussian_bao_ALL_GCcomb_cov.txt` | `https://raw.githubusercontent.com/CobayaSampler/bao_data/master/desi_bao_dr2/desi_gaussian_bao_ALL_GCcomb_cov.txt` | `252a143274c8a07c78694c119617d36594f6d7965d00319ca611c6ffb886e509` |
| `DES-Dovekie_HD.csv` | `https://raw.githubusercontent.com/des-science/DES-SN5YR/main/4_DISTANCES_COVMAT/DES-Dovekie_HD.csv` | `2f57019d783eaa976df80a41b0054171a2d994ee9808d715ce850c2df5720aaf` |
| `STAT+SYS.npz` | `https://raw.githubusercontent.com/des-science/DES-SN5YR/main/4_DISTANCES_COVMAT/STAT+SYS.npz` | `ffd3124b32148b1372bd95fda9299269f0352a9f8eee02d416c610e38495463b` |

Two additional files (used by Wren's scripts only for the Pantheon+ cross-check rows of Table 2, `z_HD > 0.01`, 1,590 SNe):

| File | Source | SHA-256 of the copy used |
|---|---|---|
| `Pantheon+SH0ES.dat` | `https://raw.githubusercontent.com/PantheonPlusSH0ES/DataRelease/main/Pantheon%2B_Data/4_DISTANCES_AND_COVAR/Pantheon%2BSH0ES.dat` | `1cb0fc379ef066afdc2ffd1857681cc478024570d8a3eba284fb645775198cf8` |
| `Pantheon+SH0ES_STAT+SYS.cov` | `https://raw.githubusercontent.com/PantheonPlusSH0ES/DataRelease/main/Pantheon%2B_Data/4_DISTANCES_AND_COVAR/Pantheon%2BSH0ES_STAT%2BSYS.cov` | `abf806d966485e64afdb359c87bffc0ecc00d05eff0a31ced66f247385df0fdc` |

All six files were downloaded on 2026-11-21.

Two compressed CMB priors are typed directly into `c43_precision_verify.py`:
- **CMB(BH):** `(R, l_a, omega_b) = (1.7504, 301.77, 0.022371)` and its `3x3` covariance at `z_* = 1090`, from Bansal & Huterer (`arXiv:2502.07185`, Eqs. 5–6).
- **CMB(DESI):** `(theta_*, omega_b, omega_bc) = (0.01041, 0.02223, 0.14208)` and its `3x3` covariance at `z_* = 1089`, from DESI Collaboration (`arXiv:2503.14738`, App. A Eqs. 35–36).
