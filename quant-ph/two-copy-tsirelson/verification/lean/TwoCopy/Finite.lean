/-
Two-copy identity, file 2 of 6: the closed computations (the finite core of target (a)), and
the tests of the definitions. The statements about the group H, the images of S as sets, and
the pointwise identity over ℝ are derived from these in `Images.lean`.

What this file proves, in plain words.

  * `M_pm`: M = 4 or M = −4 at each of the 256 points. `card_S`: S = [M = −4] has 128 points.
  * `rep_invol`: each of the four named relabelings (identity, τ₀, τ₁, τ₃) is an involution.
  * `stab`: each τ_j maps each of the four sets [M ∘ rep k = −4] onto one of the four.
    `distinct`: the four sets are pairwise distinct.
  * `pointwiseZ`: at every point 4(1_{S₀} + 1_{S₁} + 1_{S₂} + 1_{S₃}) = 8 − C₁C₂, over ℤ.
  * The computations the proof of Theorem 1 stands on: `fourier` (M ∘ rep k is the signed sum
    of monomials with signs `epsk k`), `sum_epsk` (Σ_k epsk k s = 2(−1)^(xy + x′y′)), `coeff`
    (for each difference pattern d and each restriction t, the coefficient
    Σ_{s : s·d = t} Σ_k epsk k s · epsk k (s + d) is 64 at d = 0, t = 0 and 0 elsewhere).
  * Tests of the definitions: 256 points, 256 events, 16 points per event, 22,400 unordered
    locally orthogonal pairs (the draft's count), LO is symmetric and irreflexive, values of
    the variables at named points, the action of τ_j on the variables.

Every fact marked `decide +kernel` is a closed computation checked by the Lean kernel; no
`native_decide` is used.
-/
import TwoCopy.Defs

namespace TwoCopy

open Finset

/-! ### Closed computations -/

set_option maxRecDepth 100000 in
/-- M = ±4 on Ω. -/
theorem M_pm : ∀ ω : Pt, M ω = 4 ∨ M ω = -4 := by decide +kernel

set_option maxRecDepth 100000 in
/-- S has 128 points. -/
theorem card_S : S.card = 128 := by decide +kernel

set_option maxRecDepth 100000 in
theorem rep_invol : ∀ k ω, rep k (rep k ω) = ω := by decide +kernel

set_option maxRecDepth 100000 in
/-- M ∘ rep k = Σ_s epsk k s · m_s. -/
theorem fourier : ∀ k ω, M (rep k ω) = ∑ s, epsk k s * mono s ω := by decide +kernel

set_option maxRecDepth 100000 in
/-- Each τ_j maps each of the four sets [M ∘ rep k = −4] onto one of the four. -/
theorem stab : ∀ j k, ∃ k', ∀ ω, (M (rep k (tau j ω)) = -4 ↔ M (rep k' ω) = -4) := by
  decide +kernel

set_option maxRecDepth 100000 in
/-- The four sets are pairwise distinct: a point of one outside the other. -/
theorem distinct : ∀ k k', k ≠ k' → ∃ ω, M (rep k ω) = -4 ∧ M (rep k' ω) ≠ -4 := by
  decide +kernel

set_option maxRecDepth 100000 in
/-- The pointwise identity over ℤ: 4 (1_{S₀} + 1_{S₁} + 1_{S₂} + 1_{S₃}) = 8 − C₁C₂. -/
theorem pointwiseZ :
    ∀ ω, 4 * ∑ k, (if M (rep k ω) = -4 then (1 : ℤ) else 0) = 8 - C 0 ω * C 1 ω := by
  decide +kernel

set_option maxRecDepth 100000 in
/-- Σ_k epsk k s = 2 (−1)^(x y + x′ y′). -/
theorem sum_epsk :
    ∀ s : Boxes, ∑ k, epsk k s = 2 * (-1) ^ ((s 0).val * (s 1).val + (s 2).val * (s 3).val) := by
  decide +kernel

set_option maxRecDepth 100000 in
/-- The coefficient fact. For a difference pattern d and a restriction t, the total coefficient
of the monomial pairs (s, s + d) with s·d = t in Σ_k (Σ_s epsk k s m_s)² is 64 for d = 0, t = 0
and 0 otherwise. -/
theorem coeff : ∀ d t : Boxes,
    (∑ s with s * d = t, ∑ k, epsk k s * epsk k (s + d)) = if d = 0 ∧ t = 0 then 64 else 0 := by
  decide +kernel

/-! ### Tests of the definitions -/

theorem card_Pt : Fintype.card Pt = 256 := by simp [Pt]

theorem card_Event : Fintype.card Event = 256 := by simp [Event, Boxes]

set_option maxRecDepth 100000 in
/-- Each event is 16 points of Ω. -/
theorem card_event_points : ∀ e : Event, (univ.filter fun ω : Pt => Holds e ω).card = 16 := by
  decide +kernel

set_option maxRecDepth 100000 in
theorem LO_symm : ∀ e e' : Event, LO e e' → LO e' e := by
  rintro e e' ⟨j, h1, h2⟩; exact ⟨j, h1.symm, fun h => h2 h.symm⟩

theorem LO_irrefl : ∀ e : Event, ¬ LO e e := by
  rintro e ⟨j, _, h2⟩; exact h2 rfl

set_option maxRecDepth 100000 in
/-- Locally orthogonal events are disjoint. -/
theorem LO_disjoint : ∀ e e' : Event, LO e e' → ∀ ω, ¬ (Holds e ω ∧ Holds e' ω) := by
  rintro e e' ⟨j, h1, h2⟩ ω ⟨h, h'⟩
  exact h2 (by rw [← h j, ← h' j, h1])

set_option maxRecDepth 100000 in
/-- Each event is locally orthogonal to exactly 175 events; 256 · 175 = 44,800 ordered pairs,
22,400 unordered (the draft's count). -/
theorem card_LO_row : ∀ e : Event, (univ.filter fun e' : Event => LO e e').card = 175 := by
  decide +kernel

set_option maxRecDepth 100000 in
/-- Values at named points: at the all-zero point every variable is 1, M = 4, C₁ = C₂ = 2. -/
theorem values_at_zero :
    (∀ j b, var j b (0 : Pt) = 1) ∧ M (0 : Pt) = 4 ∧ C 0 (0 : Pt) = 2 ∧ C 1 (0 : Pt) = 2 := by
  decide +kernel

set_option maxRecDepth 100000 in
/-- The variable of a box showing outcome 1 is −1. -/
theorem var_at_one : ∀ j b, var j b (fun _ _ => 1) = -1 := by decide +kernel

set_option maxRecDepth 100000 in
/-- What τ_j does to the variables: it exchanges party j's two variables, negates the variable
of box 1 of j's partner, and fixes the rest. -/
theorem var_tau : ∀ j i b ω, var i b (tau j ω) =
    if i = j then var i (b + 1) ω else if i = partner j ∧ b = 1 then -var i b ω else var i b ω := by
  decide +kernel

set_option maxRecDepth 100000 in
/-- τ_j is not the identity. -/
theorem tau_ne_id : ∀ j, ∃ ω, tau j ω ≠ ω := by decide +kernel

end TwoCopy
