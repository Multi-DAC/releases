"""Run one verification script fresh and record its output.

Usage: python -B run_one.py OUTNAME SCRIPT.py [ARGS ...]

Runs `python -B SCRIPT.py ARGS` from this folder with OPENBLAS_NUM_THREADS=1 and writes OUTNAME.out:
line 1 is the command, then the script's stdout and stderr, and the last line is the exit code and the run time.
A run is stopped after 900 s (15 minutes); the last line then says so.
"""
import os
import subprocess
import sys
import time

here = os.path.dirname(os.path.abspath(__file__))
outname, script, args = sys.argv[1], sys.argv[2], sys.argv[3:]
env = dict(os.environ, OPENBLAS_NUM_THREADS="1", PYTHONDONTWRITEBYTECODE="1", PYTHONIOENCODING="utf-8")
cmd = [sys.executable, "-B", script] + args
shown = "python -B " + " ".join([script] + args)
t0 = time.time()
with open(os.path.join(here, outname + ".out"), "w", encoding="utf-8", newline="\n") as f:
    f.write(shown + "\n")
    f.flush()
    try:
        p = subprocess.run(cmd, cwd=here, env=env, stdout=f, stderr=subprocess.STDOUT, timeout=900)
        code = p.returncode
        note = ""
    except subprocess.TimeoutExpired:
        code = "none"
        note = " (stopped at the 900 s limit)"
    f.flush()
with open(os.path.join(here, outname + ".out"), "a", encoding="utf-8", newline="\n") as f:
    f.write("exit code %s, run time %.1f s%s\n" % (code, time.time() - t0, note))
