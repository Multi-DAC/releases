# Wren's Section 3 (and Independent Section 2) Figures Map

Maps every computed figure in Section 3 of the note (and Wren's independent cross-checks referenced in the main README) to the exact script and output line in `verification/scripts/wren/`.

Three scripts run from scratch via `rerun.sh` (`rerun.log` records the Python/NumPy/SciPy versions, SHA-256 hashes of all six public data files, exit codes, and timings):
- `c43_precision_verify.py` -> `c43_precision_verify.out` (Table 1, CMB rows of Table 2, 1-sigma profile CIs on `c`, Section 3.3 BDF numerical robustness, DESI DR2 diagonal BAO pulls).
- `c43_profiles_and_nocmb.py` -> `c43_profiles_and_nocmb.out` (1-sigma profile CI on `q`, linear H-form wide `c` scan `0.05..1000`, no-CMB rows of Table 2 `DD+DESI` and `PP+DESI`, matter-only `c_min`, Section 2.6 `1+z(N_e)` exclusion bound).
- `c43_sec2_independent.py` -> `c43_sec2_independent.out` (standalone, needs no data files: Lemma 2 four switches, Section 2.6 9-row table + `c=0.002`, `1+z(N_e)` exclusion bound, fitted quadrature one-maximum check, matter-only `c_min`).

## Section 3.1 (Datasets, Replication, and External Calibration)

| Note | Figure | Script, output line |
|---|---|---|
| s3.1 | Dataset counts: 13 DESI BAO, 1,820 DES-Dovekie, 1,590 Pantheon+ (`z_HD > 0.01`) | `c43_precision_verify.out`, line 1: `Sample counts: DESI BAO = 13, DES-Dovekie = 1820, Pantheon+ (zHD > 0.01) = 1590` |
| s3.1 | Unrounded `chi^2` & `z_c` for all 9 replicated rows (`<= 0.0083` viable, `<= 0.060` all nine; `z_c <= 0.002`) | `c43_precision_verify.out`, lines 5–15 (4-decimal `chi^2` in parentheses and `zc=...`), compared against Clawd's `c43_fit.out` |
| s3.1 | Matter-only `c_min(Omega_m = 0.31)`: `1.07267` (quadrature), `1.59423` (Lorentzian) | `c43_profiles_and_nocmb.out`, lines 54–55; also `c43_sec2_independent.out`, lines 32–33 |
| s3.1 | No-CMB external calibration vs [L+]: `DD+DESI` `-7.29` (`w0wa`), `-7.28` (logistic); `PP+DESI` `-5.01` (`w0wa`), `-6.14` (logistic) | `c43_profiles_and_nocmb.out`, lines 37–38 (`DD+DESI`) and 46–47 (`PP+DESI`) |
| s3.1 | `CMB(BH)+DESI` `w0wa` calibration: `w_0 = -0.44, w_a = -1.67, H_0 = 63.74` | `c43_precision_verify.out`, line 36 (`CMB(BH)+DESI`, `w0wa`) |

## Section 3.2 (Table 1 and Conclusions 1–4)

| Note | Figure | Script, output line |
|---|---|---|
| s3.2, Table 1 | All 11 rows on `CMB(BH)+DD+DESI` (`N = 1836`): `chi^2`, `dchi^2`, `dAIC`, `dBIC`, `dBAO/dSN/dCMB`, `H_0`, `Omega_m`, switch/DE params, `w_0`, `z_c` | `c43_precision_verify.out`, lines 5–15 |
| s3.2, Table 1 & item 3 | 1-sigma profile CI on `q`: `q = 0.915 (+0.584 / -0.272)` (`0.643`–`1.499`); `dchi^2(q=1) = +0.045`, `dchi^2(q=1/2) = +3.51`, Gaussian `+6.78` (`1645.84 - 1639.06`) | `c43_profiles_and_nocmb.out`, lines 2, 4, 9, 15; Gaussian `chi^2 = 1645.84` in `c43_precision_verify.out`, line 14 |
| s3.2, Table 1 & item 4 | 1-sigma profile CIs on `c`: Lorentzian `2.787 (+0.237 / -0.213)` (`Hbar = 50.02`), quadrature `1.641 (+0.155 / -0.138)` | `c43_precision_verify.out`, lines 12–13, 16–17 |
| s3.2, item 1 | Linear H-form wide `c` profile scan over `c in [0.05, 1000]`: minimum `+8.51` at `c = 2.469`, `w(1) = -0.904`, `w_0 = -1.011`, `z_c = 0.166*` | `c43_profiles_and_nocmb.out`, lines 18–31; `c43_precision_verify.out`, line 11 |

## Section 3.3 (Table 2, BAO Pulls, and BDF Numerical Robustness)

| Note | Figure | Script, output line |
|---|---|---|
| s3.3, Table 2 | `CMB(BH)+PP+DESI` (`N = 1606`) row + 1-sigma CIs on `c` (`2.794 +0.257/-0.229`, `1.622 +0.166/-0.147`) | `c43_precision_verify.out`, lines 20–32 |
| s3.3, Table 2 | `CMB(DESI)+DD+DESI` (`N = 1836`) row | `c43_precision_verify.out`, lines 48–58 |
| s3.3, Table 2 | `CMB(BH)+DESI` (no SN, `N = 16`) row | `c43_precision_verify.out`, lines 35–45 |
| s3.3, Table 2 | `DD+DESI` (no CMB, `N = 1833`) and `PP+DESI` (no CMB, `N = 1603`) rows, including linear H-form `H_0 ~ 59`, `c ~ 4.4` (`4.417` and `4.352`) | `c43_profiles_and_nocmb.out`, lines 36–42 (`DD+DESI`) and 45–51 (`PP+DESI`) |
| s3.3 | Lorentzian `w_0 = -0.899, w(0.5) = -0.9945, w(1) = -1.217, z_c = 0.516` | `c43_precision_verify.out`, line 13 |
| s3.3 | DESI DR2 diagonal BAO pulls `(data - model)/sigma_i` for `LCDM` vs Lorentzian across all 13 measurements | `c43_precision_verify.out`, lines 66–78 |
| s3.3 | `u`-evolution BDF numerical robustness: `N_start in {3.5, 4.5, 5.5}` (`<= 8.1e-12` Lor., `<= 4.3e-12` quad.), `rtol=1e-12, atol=1e-14` (`<= 2.2e-11`), `D_M(1090)` direct vs `z = 10` splice (`1.20e-6` Lor., `2.64e-8` quad.) | `c43_precision_verify.out`, lines 61–62 |

## Independent Checks of Section 2 (Referenced in Main README)

| README / Note | Figure | Script, output line |
|---|---|---|
| s2.7, Lemma 2 | Four matter-only switches `u(-4)` and `x/h(-4)` | `c43_sec2_independent.out`, lines 2–5 |
| s2.6, table | Nine integrated `1 + z_max` values (agreeing with Clawd to `0.01`), `(2.9)` row, and `c = 0.002` no-extremum check | `c43_sec2_independent.out`, lines 9–18 |
| s2.6, exclusion | Independent `1 + z(N_e)` for `c = 3.1, 3.5, 4, 6, 10, 20` (`46.3123` at `c = 4`) | `c43_sec2_independent.out`, lines 21–26 (also `c43_profiles_and_nocmb.out`, lines 58–64) |
| s2.4, s2.8 | Fitted quadrature `c = 1.641` unique maximum on `0 <= z <= 1090` at `z = 0.473` (`u = 1.2827`) | `c43_sec2_independent.out`, line 29 |
