#!/bin/sh
# Re-run every script in this folder from scratch and record each output beside it as X.out.
# Usage: C43_DATA=/path/to/the/four/data/files  PYTHON="python -B"  sh rerun.sh
# The fits come first: one_maximum.py reads c43_fit_all.json.
set -u
cd "$(dirname "$0")"
PYTHON="${PYTHON:-python -B}"
export OPENBLAS_NUM_THREADS=1
rm -f ./*.out c43_fit_all.json rerun.log
{
  echo "re-run started $(date '+%Y-%m-%d %H:%M:%S')"
  $PYTHON -c "import sys, numpy, scipy; print('python', sys.version.split()[0], '| numpy', numpy.__version__, '| scipy', scipy.__version__)"
  sha256sum "$C43_DATA/desi_gaussian_bao_ALL_GCcomb_mean.txt" "$C43_DATA/desi_gaussian_bao_ALL_GCcomb_cov.txt" \
            "$C43_DATA/DES-Dovekie_HD.csv" "$C43_DATA/STAT+SYS.npz" | sed 's|\*.*/|*|'
} > rerun.log 2>&1
run() {
  out="$1"; shift
  t0=$(date +%s)
  $PYTHON "$@" > "$out" 2>&1
  code=$?
  echo "$* -> $out : exit $code, $(( $(date +%s) - t0 )) s" >> rerun.log
}
run c43_fit.out c43_fit.py
run c43_controls.out c43_controls.py
run c43_props.out c43_props.py
run sec26_table.out sec26_table.py
run sec26_exclusion.out sec26_exclusion.py
run one_maximum.out one_maximum.py
run c43_scan_lin.out c43_scan.py lin
echo "re-run finished $(date '+%Y-%m-%d %H:%M:%S')" >> rerun.log
cat rerun.log
