# Verification

One part: `scripts/`. There is no machine-checked proof for this release and no measurement. What is here is the
code behind every computed figure in the write-up, with its recorded output, the print files with the code that
builds and tests them, and the reading program with its self-test. The main README says, claim by claim, which
of these were checked by a second code or a second person and which were not.

Python 3 with numpy and scipy; the print files also need manifold3d and trimesh, the pictures matplotlib, and
the reading program on the computer pyserial. All of this is Clawd's code. Wren's own codes (the second field
code, the three force routes, the ring-gap table) are not filed here. Comments in the scripts refer to numbered
letters between Clawd and Wren; the letters are not part of the release. W16 is the project's working label, S3
the plain Halbach holder and S4 the tilted one.

## scripts/

**The field of a block.**

- `w16_prism.py`: the exact field of uniformly magnetised rectangular blocks, in closed form. Every other field
  script imports it.
- `w16_prism_test.py` -> `out-prism-test.txt`: tests of that kernel, each with a wrong kernel that must fail.

**Choosing the design** (10 mm cubes; the build is the same design with every length times 0.9525).

- `w16_desk.py` -> `out-desk.txt`, `desk-envelope-1121.json`: plain Halbach rings of catalogue cubes over a
  range of sizes.
- `w16_desk_cmp.py` -> `out-desk-cmp.txt`: the same cube centres under the two rules, plain and tilted.
- `w16_wren_checks.py` -> `out-wren-checks.txt`, `w16_r_bump.py` -> `out-r-bump.txt`: three of Wren's numbers
  for an earlier geometry reproduced on this kernel, and both holders with the circle of cube centres enlarged
  by the same amount.
- `w16_envelope_dims.py` -> `out-envelope-dims.txt`: where the cube corners reach, for the holder's outline.

**The prediction and what moves it.**

- `w16_seat_envelope.py` -> `out-seat-envelope.txt`: at the build scale (3/8 inch, 1.30 T): centre fields and
  the ratio, parity of the field components, curvatures, the sensor's allowed box, the field along the bore
  axis, the height fit, the room field.
- `w16_seating.py` -> `out-seating.txt`: cube strength, turning and shifting in the pocket, radius errors,
  random position and direction (the first seven rows of the table in Section 3).
- `w16_rehberg_cmp.py` -> `out-rehberg-cmp.txt`: the tilted rule against the two printed formulas of Rehberg and
  Bluemler at our twelve angles, the same ratio for point dipoles, and a wrong formula that must disagree.

**Forces and loading.**

- `w16_forces.py` -> `out-forces.txt`: uniformity, and the force and torque on every seated cube, with controls
  on the force code. This is the output that was fixed by hash before Wren's blind check.
- `w16_loading_order.py` -> `out-loading-order.txt`: an earlier study, every stage of two loading orders with
  both rings in place (round the ring, and opposite pockets in turn). The write-up uses neither order.
- `w16_single_ring.py` -> `out-single-ring.txt`: one ring loaded by itself in the order used, then two loaded
  rings brought together.
- `w16_cradle_side.py` -> `out-cradle-side.txt`: which half of a ring should be the cradle, for two orders.

**Print files** (`scripts/print/`). Each builder writes its STL files into the directory it is given and prints
its checks; the filed STL files are in `out/`, `out-holder/` and `out-seat/`.

- `holder.py <out_dir> [pocket_mm]` -> `out-holder.txt`, four STL files (cradle and lid for S3 and S4). The
  pocket size defaults to the 9.60 mm placeholder. **Regenerate with the gauge's result before printing.**
- `pole_faces.py` -> `out-pole-faces.txt`: which cube faces show from above in each pocket of a cradle.
- `gauge_strip.py <out_dir>` -> `out-gauge-strip.txt`, `gauge_cap.py <out_dir>` -> `out-gauge-cap.txt`: the
  pocket gauge and the cap for its tilted nests.
- `sensor_seat.py <out_dir>` -> `out-sensor-seat.txt` (24 checks): the plug and four shims.
- `render_holder.py`, `render_gauge.py`, `render_seat.py`: the pictures, drawn from the STL files or from
  slices of the solids.
- `brd_read.py` -> `out-brd-6490.txt`: outline, holes and part positions read from Adafruit's Eagle board file
  for product 6490. The board file is Adafruit's and is not copied here (their repository
  `Adafruit-TMAG5273-3D-Hall-Effect-Magnetometer-Breakout-PCB`, commit 543c3c9).

**Reading program** (`scripts/read/`).

- `code.py`: runs on the Trinkey QT2040 under CircuitPython. It needs Adafruit's `adafruit_tmag5273` driver
  and its two support libraries on the board.
- `w16_read.py`: runs on the computer. `collect` takes one seating at a time, `analyse` gives the ratio on each
  range, `zfit` fits the height scan.
- `test_w16_read.py` -> `out-test.txt`: 18 checks on made-up readings, no hardware. They are: exact recovery
  of a known ratio over 200 random sensor orientations and offsets on both ranges; two wrong formulas that
  must miss (half-sum, single seating); a half turn that is off by 2 and by 5 degrees; noise alone; re-seating
  scatter; a sensor axis that reads 5% high, at 0 and 45 degrees; the warning for a field not on one axis; the
  two ranges kept apart; the line parser accepting a data line and refusing eight others; the collector
  dropping readings after a range switch and reporting a stream that ends early; a round trip through the
  data file; the warnings for a missing seating, a reading near full scale and a temperature difference; and
  the height fit, inside and outside the heights tried.

**Rerun.** `rerun.sh <new_directory>` runs every script above except the three picture scripts and
`brd_read.py` into a fresh directory and compares each output, and each STL file byte for byte, with the filed
one. `rerun.log` is the run of 22 November 2026. That run shows the filed scripts produce the filed outputs and
print files on one machine. It is not a second implementation of anything.

## What the scripts do not establish

- A field script evaluates a model: rigid, uniform magnetisation, exact cube positions. It says nothing about
  the bought cubes or the printed holders.
- The checks in `print/` are tests of solids in memory and of the files written from them. They do not show
  that a part prints, fits, or holds a cube.
- The reading program's self-test checks its arithmetic and its file handling against simulated readings. The
  program has never talked to the sensor.
