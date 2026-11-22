# The public data files the scripts read

The data are not copied into this repository. `c43_fit.py` reads four files from the directory named by the
environment variable `C43_DATA` (default: a folder `data/` beside the script). The other scripts that import
`c43_fit.py` need the files only because it loads them when imported.

| File | Source | SHA-256 of the copy used |
|---|---|---|
| `desi_gaussian_bao_ALL_GCcomb_mean.txt` | DESI DR2 BAO, 13 points: `https://raw.githubusercontent.com/CobayaSampler/bao_data/master/desi_bao_dr2/desi_gaussian_bao_ALL_GCcomb_mean.txt` | `9ac154ab583ce759c0f7eef3c978c7c70a6ead2d18774caceadf1a350a640585` |
| `desi_gaussian_bao_ALL_GCcomb_cov.txt` | the same folder, `desi_gaussian_bao_ALL_GCcomb_cov.txt` | `252a143274c8a07c78694c119617d36594f6d7965d00319ca611c6ffb886e509` |
| `DES-Dovekie_HD.csv` | DES-Dovekie Hubble diagram, 1,820 supernovae: `https://raw.githubusercontent.com/des-science/DES-SN5YR/main/4_DISTANCES_COVMAT/DES-Dovekie_HD.csv` | `2f57019d783eaa976df80a41b0054171a2d994ee9808d715ce850c2df5720aaf` |
| `STAT+SYS.npz` | the same folder, `STAT+SYS.npz`: the inverse covariance, statistical plus systematic, upper triangle | `ffd3124b32148b1372bd95fda9299269f0352a9f8eee02d416c610e38495463b` |

All four were downloaded on 2026-11-21 from the `master` and `main` branches named in the addresses. A later
change to either repository would change the hashes; `rerun.log` prints the hashes of the files a run used.

The compressed CMB likelihood is not a file. Its data vector (R, l_a, omega_b) = (1.7504, 301.77, 0.022371) and its
covariance are typed into `c43_fit.py` from Bansal and Huterer, arXiv:2502.07185, Eqs. 5 and 6, read in the
TeX source.

Constants: T_CMB = 2.7255 K, N_eff = 3.044, one massive neutrino of 0.06 eV (the fitting formula of Bansal and
Huterer's Appendix B), z_* = 1090, the DESI DR2 fitting formula for the sound horizon r_d.
