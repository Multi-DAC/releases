/-
Mutants (lane lean2c, evidence only; not part of the library, not built by lake).
Each statement below is a false variant of a fact of `Finite.lean`; `decide +kernel` must
REJECT each one. Run from the worktree:
  lake env lean <this file>
Expected: one error per theorem ("decide failed ... proposition ... is false"), 4 errors.
-/
import TwoCopy.Defs

namespace TwoCopy.Mutants

open Finset TwoCopy

set_option maxRecDepth 100000 in
/-- mutant of `coeff`: 63 instead of 64. -/
theorem coeff_63 : ∀ d t : Boxes,
    (∑ s with s * d = t, ∑ k, epsk k s * epsk k (s + d)) = if d = 0 ∧ t = 0 then 63 else 0 := by
  decide +kernel

set_option maxRecDepth 100000 in
/-- mutant of `coeff`: the single set (k = 0 only) instead of the four. -/
theorem coeff_single : ∀ d t : Boxes,
    (∑ s with s * d = t, eps s * eps (s + d)) = if d = 0 ∧ t = 0 then 16 else 0 := by
  decide +kernel

set_option maxRecDepth 100000 in
/-- mutant of `sum_epsk`: the sign of the third image flipped. -/
theorem sum_epsk_mut : ∀ s : Boxes,
    epsk 0 s + epsk 1 s - epsk 2 s + epsk 3 s
      = 2 * (-1) ^ ((s 0).val * (s 1).val + (s 2).val * (s 3).val) := by
  decide +kernel

set_option maxRecDepth 100000 in
/-- mutant of `M_pm`: M is always 4. -/
theorem M_always_4 : ∀ ω : Pt, M ω = 4 := by decide +kernel

end TwoCopy.Mutants
