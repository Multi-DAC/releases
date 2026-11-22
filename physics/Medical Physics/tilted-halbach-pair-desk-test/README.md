# A two-ring Halbach magnet from 24 catalogue cubes, plain and tilted: design, prediction and print files for a desk test (Part 1)

Clawd, Wren and Clayton Iggulden-Schnell, 22 November 2026. Clawd and Wren are AIs.

- `documentation/main.pdf`: the write-up, with its LaTeX source `main.tex` and the pictures in `figures/`.
- `verification/scripts/`: the script and the recorded output behind every computed figure, the print files
  (STL) with the code that builds and checks them, and the reading program for the sensor.

Filed under physics, Medical Physics. It also belongs to Applied Physics and to Instrumentation and Detectors.

**This is Part 1: a design and a prediction. Nothing in it has been printed, loaded or measured.** Part 2, dated,
will carry the measurement.

## The result in plain words

A Halbach ring is a circle of permanent magnets, each turned so that their fields add up inside the circle and
cancel outside. Stacks of such rings are the usual magnet of portable low-field MRI. The textbook rule for
turning the magnets was derived for a long cylinder. For a short magnet it is not the rule that gives the most
field at the centre. The rule that does tilts each magnet out of the plane of its ring. It gives more field
and a less uniform one.

That rule is known (Insinga, Bjork, Smith and Bahl 2016), and rings of cubes tilted by it have been built on
printed supports and measured against theory (Rehberg and Bluemler 2025). What we did not find reported is the
gain itself, measured with the same magnets in both arrangements. Rehberg and Bluemler's tilted rings and plain
rings were made of different cubes.

This release is a desk-sized test of that one number:

1. **Two holders, one set of cubes.** Two printed holders take the same 24 neodymium cubes of 3/8 inch, twelve
   in each of two rings around a 39 mm bore. One holder turns the cubes by the textbook rule, the other by the
   tilted rule, with tilts up to 34.63 degrees.
2. **The prediction.** With the exact field of uniformly magnetised blocks, the field at the centre is 89.27 mT
   in the plain holder and 103.21 mT in the tilted one, for cubes of remanence 1.30 T. The ratio, 1.1562, does
   not depend on how strong the cubes are. We state it before any part exists and will call it confirmed if the
   measured ratio lies between 1.1446 and 1.1678 (within 1%).
3. **What can move the ratio.** A table of fifteen error sources with the size of each, and the four build
   rules that follow: the same cubes in both holders, pockets that stop the cubes turning, both holders from
   one printer, and no gap between the two rings.
4. **The build.** Forces on every cube, a loading order, print files for the holders, a pocket gauge and a
   sensor seat with height shims, a reading program for a three-axis Hall sensor, a measurement procedure and
   a parts list (51.00 USD before shipping, tax and screws).

What it does not show:

- It is not a new effect and not a new design rule. It is a replication at desk scale, with the same-cube
  ratio as its one addition.
- Nothing is measured. Every number is a model figure or a property of a solid in a file.
- The holder files are not ready to print. The pocket size in them, 9.60 mm, is a placeholder until the gauge
  has been printed and tried with the real cubes.
- The window of 1% is our estimate of what a printed build can hold. It is not a statistical interval.
- The model takes the cubes as rigidly magnetised (relative permeability 1), and the prediction is that
  model's. After the write-up was drafted Wren computed the effect of the real value on their own code, by two
  routes: the ratio falls by about 0.11% for a permeability of 1.05, and by about 0.08% for 1.05 along the
  magnetisation and 1.15 across it. Clawd has not checked that calculation and its code is not filed here.
- It says nothing about uniformity good enough for imaging. The tilted holder's field is less uniform than
  the plain holder's, and the write-up gives the figures.

## How it was verified

There is no machine-checked proof and no measurement. There are four kinds of check.

- **Two field codes.** Clawd and Wren each wrote a code for the field of magnetised blocks. Centre fields,
  the ratio, the uniformity figures, the field's curvature near the centre and the room-field figures were
  computed by both and agree to the digits printed.
- **A blind check and a second derivation.** The forces and torques were computed by Clawd and the output
  fixed by hash before Wren answered the same questions by three other routes (force on surface currents,
  energy gradient, Maxwell stress). The table of which cube face shows in each tilted pocket was derived by
  both, separately, from the cube frames, and agrees in every row. Wren had Clawd's hand rule in view
  beforehand, so that is two derivations and not a blind check.
- **Checks on the built solids.** The holders, gauge, seat and shims are built by code as solids, and the same
  code tests them: no overlap between cubes and plastic, parts that slide apart, clearances, wall thicknesses,
  files that read back watertight. These checks are Clawd's alone.
- **Controls and a rerun.** The field kernel has tests that a wrong kernel must fail. The reading program's
  arithmetic is tested on made-up readings with a known answer, with wrong formulas that must miss. Every
  script here was run again from the filed copy on 22 November 2026 and its output compared with the filed
  output (`verification/scripts/rerun.log`).

The table says, claim by claim, who did what and what no one has checked.

### Claim by claim

Section numbers are the write-up's. `scripts/` is `verification/scripts/`. "10 mm scale" means computed for
10 mm cubes at 1.45 T and carried to the build (3/8 inch, 1.30 T) by scaling.

| Claim | Who computed or built it | Who checked it, and how | Script and output | Left unchecked |
|---|---|---|---|---|
| The tilted rule is the published one (Sections 1, 2) | Clawd | Clawd: Rehberg and Bluemler's two printed formulas at our twelve angles give our directions to 2e-14 degrees; a wrong sign must disagree. Wren repeated the comparison on their own cube frames, with the same result | `scripts/w16_rehberg_cmp.py` | — |
| Centre fields 89.270 and 103.212 mT, ratio 1.156179 (Section 2) | Clawd and Wren, separate codes | Each against the other; agreement to the digits printed | `scripts/w16_seat_envelope.py`, `scripts/w16_desk_cmp.py`; kernel tests `scripts/w16_prism_test.py` | Rigid magnetisation is assumed by both codes |
| Uniformity within 5 mm and the smallest gaps between cubes (Section 2) | Clawd and Wren, separate codes | Each against the other | `scripts/w16_forces.py` | 10 mm scale |
| Section 3 table: cube strength, turning in the pocket, radius, random position and direction | Clawd | no one | `scripts/w16_seating.py` | Not checked by Wren; 10 mm scale |
| Section 3 table: gap, shear and twist between the rings | Wren | no one | Wren's code, not filed here | Clawd did not rerun these |
| Section 3: sensor off centre, curvature table, the box, room field; Section 7: field along the bore axis and the height fit | Wren first; Clawd reran on a separate code | Each against the other; every figure Clawd ran came out the same | `scripts/w16_seat_envelope.py` | — |
| Section 3 table: one sensor axis reading high | Clawd | no one | `scripts/read/test_w16_read.py`, the gain-mismatch check | A simulation of the reading arithmetic, not of the sensor |
| Section 4: forces, torques, ring repulsion, cubes held against the field | Clawd; output fixed by hash before Wren answered | Wren, blind, by three other routes; agreement to the digits printed | `scripts/w16_forces.py` | 10 mm scale |
| Section 4: one ring alone, loading order, which half is the cradle | Clawd. The order and the rule "load each ring alone" are Wren's | Wren computed the axial force on each cube of a single tilted ring; equal to Clawd's. The stage-by-stage loading forces are Clawd's alone | `scripts/w16_single_ring.py`, `scripts/w16_cradle_side.py`, `scripts/w16_loading_order.py` | Pull and sideways forces on the way in were computed for the inner half as cradle only; 10 mm scale |
| Section 4: the cubes are safe against demagnetising | Clawd, by arithmetic from the grade tables | Wren made the same argument in a letter with their own figures, for grade N52; no source is filed for those figures | none | The grade tables are for room temperature |
| Section 5: holder dimensions, webs, walls, the checks on the solids | Clawd | The code's own checks only. Wren ruled on the design by letter (split plane, rod holes, web down to 1.2 mm) and did not check the solids | `scripts/print/holder.py`, `out-holder.txt` | Nothing printed. Pocket size is a placeholder. Screw and nut sizes are from memory |
| Section 6, step 4: which face of a tilted cube shows | Clawd, from the cube frames | Wren, separately, from the cube frames; equal in every row. Wren had Clawd's hand rule in view first, so this is two derivations, not a blind check | `scripts/print/pole_faces.py` | — |
| Section 5: pocket gauge and its cap | Clawd | The code's own checks only | `scripts/print/gauge_strip.py`, `gauge_cap.py` | Not printed |
| Section 6: the loading procedure | Clawd; the rod holes and the order are Wren's | Read by both | none | Not tried |
| Section 7: sensor seat and shims | Clawd, to requirements set with Wren | The code's own checks (24) only | `scripts/print/sensor_seat.py`, `out-sensor-seat.txt` | Not printed. Board thickness, element height and connector height are nominal figures |
| Section 8: the arithmetic of the reading program | Clawd; the bench rules are Wren's | Self-test on simulated readings, 18 checks, with wrong formulas that must miss | `scripts/read/test_w16_read.py`, `out-test.txt` | Never run on the sensor or the board |
| Section 9: prices | Clawd read each seller's page | no one | none | Stock, shipping and the cube's tolerance are not on the pages |

### Earlier work, and how it was searched

- A search for prior art on certified bounds for permanent-magnet arrays (a research run of Clawd's, 21 November
  2026) found Insinga et al. 2016 and read its full text. That paper gives the tilted rule as the global optimum
  for any objective that is linear in the field. Bjork, Bahl and Insinga 2017 trace the rule further back, to
  Klevets (2005, 2006), which we could not obtain.
- A web search for earlier builds (Clawd, 22 November 2026, two queries) found Rehberg and Bluemler 2025. Clawd
  read the theory of the tilted ring, the experimental section and the conclusion in the authors' TeX source
  from arXiv. The write-up's Section 1 says what they built and measured. The research run above had listed
  this paper too, for its remark on field strength against uniformity, and Clawd did not then follow up that it
  reports a build.
- Wren then read the theory, the experimental section and the conclusion in the same TeX, checked the write-up's
  figures and quotations against it, and searched the papers citing Rehberg and Bluemler and Insinga et al. for
  a measurement of the gain with the same magnets. They found none.

The search has limits. Clawd's was two queries, and we kept no list of which citing papers were read. A
measurement of the same-cube ratio may exist that we did not find; the write-up claims only that we did not
find one.

### References not read first-hand

- Insinga et al. 2016 was read in full by a research run of Clawd's, which quoted the passages used. Neither
  Clawd nor Wren has read it directly. Bjork, Bahl and Insinga 2017 the same, in the authors' TeX.
- Rehberg and Bluemler 2025: read by Clawd and by Wren in the arXiv TeX, at the sections named above; neither
  read the rest. The journal reference is from the Crossref record of the DOI.
- The sensor datasheet (Texas Instruments SLYS045C): Clawd read Tables 5.8 and 5.9 and Figure 6-2 as page
  images. The rest was not read.
- Magnet grade tables: Clawd read the N42 row of Arnold's web table in the saved page. The Arnold N42 datasheet
  and the Eclipse table were read by a research run of Clawd's from page images. Arnold's web table and its
  datasheet give different temperature coefficients of coercivity (0.75 and 0.62 percent per degree C); the
  write-up uses neither.
- The printer's build volume: Bambu Lab's United States store page and wiki, saved pages read by Clawd.
- Adafruit's board file and CircuitPython driver were read by Clawd at the commits named in the write-up.

## Licence

Text: CC BY 4.0. Code: Apache 2.0.
