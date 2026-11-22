#!/bin/sh
# Run every script from scratch into a fresh directory and compare with the outputs filed here.
# Usage: sh rerun.sh <empty_or_new_directory>      (PYTHON=<interpreter> to choose the Python; default "python")
# A line that begins "wrote " (it names the output directory) is ignored in the comparison.
# Needs: numpy, scipy, manifold3d, trimesh, matplotlib. Nothing filed here is overwritten.
set -u
HERE=$(cd "$(dirname "$0")" && pwd)
OUT=${1:?give a directory for the new outputs}
PY=${PYTHON:-python}
export OPENBLAS_NUM_THREADS=1
mkdir -p "$OUT/print/out" "$OUT/print/out-holder" "$OUT/print/out-seat" "$OUT/read"
same=0; diff=0

cmp_txt() {   # <new> <filed>
    if diff -q -I '^wrote ' --strip-trailing-cr "$1" "$2" >/dev/null 2>&1; then same=$((same+1)); echo "SAME   ${2#$HERE/}"
    else diff=$((diff+1)); echo "DIFFER ${2#$HERE/} ($(diff --strip-trailing-cr "$1" "$2" | grep -c '^[<>]') lines)"; fi
}
cmp_bin() {   # <new> <filed>
    if cmp -s "$1" "$2"; then same=$((same+1)); echo "SAME   ${2#$HERE/}"; else diff=$((diff+1)); echo "DIFFER ${2#$HERE/}"; fi
}
run() {       # <script> <output name> [args...]
    s=$1; o=$2; shift 2
    start=$(date +%s)
    "$PY" -I -B -X utf8 "$HERE/$s" "$@" > "$OUT/$o" 2> "$OUT/$o.err"
    echo "ran    $s (exit $?, $(( $(date +%s) - start )) s)"
    cmp_txt "$OUT/$o" "$HERE/$o"
}

run w16_prism_test.py     out-prism-test.txt
run w16_desk.py           out-desk.txt "$OUT/desk-envelope-1121.json"
cmp_txt "$OUT/desk-envelope-1121.json" "$HERE/desk-envelope-1121.json"
run w16_desk_cmp.py       out-desk-cmp.txt
run w16_wren_checks.py    out-wren-checks.txt
run w16_r_bump.py         out-r-bump.txt
run w16_forces.py         out-forces.txt
run w16_loading_order.py  out-loading-order.txt
run w16_single_ring.py    out-single-ring.txt
run w16_cradle_side.py    out-cradle-side.txt
run w16_envelope_dims.py  out-envelope-dims.txt
run w16_seating.py        out-seating.txt
run w16_seat_envelope.py  out-seat-envelope.txt
run w16_rehberg_cmp.py    out-rehberg-cmp.txt

run print/gauge_strip.py  print/out-gauge-strip.txt "$OUT/print/out"
run print/gauge_cap.py    print/out-gauge-cap.txt   "$OUT/print/out"
run print/holder.py       print/out-holder.txt      "$OUT/print/out-holder"
run print/pole_faces.py   print/out-pole-faces.txt
run print/sensor_seat.py  print/out-sensor-seat.txt "$OUT/print/out-seat"
for f in out/w16-gauge-strip.stl out/w16-gauge-cap.stl \
         out-holder/w16-S3-outer-cradle.stl out-holder/w16-S3-inner-lid.stl \
         out-holder/w16-S4-outer-cradle.stl out-holder/w16-S4-inner-lid.stl \
         out-seat/w16-sensor-seat.stl out-seat/w16-seat-shim-0.6mm.stl out-seat/w16-seat-shim-1.0mm.stl \
         out-seat/w16-seat-shim-2.0mm.stl out-seat/w16-seat-shim-3.0mm.stl; do
    cmp_bin "$OUT/print/$f" "$HERE/print/$f"
done

run read/test_w16_read.py read/out-test.txt

echo "same $same, different $diff"
