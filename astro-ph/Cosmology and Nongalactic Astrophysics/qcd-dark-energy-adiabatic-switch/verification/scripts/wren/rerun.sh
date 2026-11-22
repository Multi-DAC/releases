#!/bin/sh
# Re-run Wren's Section 3 verification scripts from scratch and record outputs beside them.
# Usage: C43_DATA=/path/to/the/six/data/files  PYTHON="python -B"  sh rerun.sh
set -u
cd "$(dirname "$0")"
PYTHON="${PYTHON:-python -B}"
if [ -z "${C43_DATA:-}" ]; then
  if [ -d ./data ]; then C43_DATA=./data; else C43_DATA=../data; fi
fi
export C43_DATA
export OPENBLAS_NUM_THREADS=1
rm -f ./*.out rerun.log
{
  echo "re-run started $(date '+%Y-%m-%d %H:%M:%S')"
  $PYTHON -c "import sys, numpy, scipy; print('python', sys.version.split()[0], '| numpy', numpy.__version__, '| scipy', scipy.__version__)"
  sha256sum "$C43_DATA/desi_gaussian_bao_ALL_GCcomb_mean.txt" \
            "$C43_DATA/desi_gaussian_bao_ALL_GCcomb_cov.txt" \
            "$C43_DATA/DES-Dovekie_HD.csv" \
            "$C43_DATA/STAT+SYS.npz" \
            "$C43_DATA/Pantheon+SH0ES.dat" \
            "$C43_DATA/Pantheon+SH0ES_STAT+SYS.cov" | sed 's|\*.*/|*|'
} > rerun.log 2>&1
run() {
  out="$1"; shift
  t0=$(date +%s)
  $PYTHON "$@" > "$out" 2>&1
  code=$?
  echo "$* -> $out : exit $code, $(( $(date +%s) - t0 )) s" >> rerun.log
}
run c43_precision_verify.out c43_precision_verify.py
run c43_profiles_and_nocmb.out c43_profiles_and_nocmb.py
run c43_sec2_independent.out c43_sec2_independent.py
echo "re-run finished $(date '+%Y-%m-%d %H:%M:%S')" >> rerun.log
cat rerun.log
