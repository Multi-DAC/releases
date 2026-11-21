/-
Two-copy identity, file 4 of 6: Step 1 of the draft's proof, as lemmas (target (b)).

What this file proves, in plain words. D is any real bilinear form on the functions on Ω
(`Form`), and the only hypothesis is (ii), `Decoheres D`: D(1_e, 1_e′) = 0 for all locally
orthogonal events e, e′. No symmetry of D, no behaviour P, no normalisation.

  * `step1`: for a party j, a box b, α = var(j, b), and F, G in W_{−j} (the span of the
    indicators of events of the other three parties): D(αF, αG) = D(F, G).
  * `D_one_ind`, `D_ind_one`: D(1, 1_e) = D(1_e, 1_e) = D(1_e, 1) for every event e.
  * `D_one_one`: D(1, 1) = Σ_o D(1_(s,o), 1_(s,o)) for every choice of boxes s. So under (i) the
    behaviour has the same total in every context, namely D(1, 1).
  * Membership lemmas: a function of the four outcomes read in the boxes s lies in W
    (`boxfun_mem_W`); if it does not depend on party j's outcome it lies in W_{−j}
    (`boxfun_mem_Wexc`).
-/
import TwoCopy.Images

namespace TwoCopy

open Finset

/-- A real bilinear form on the functions on Ω. -/
abbrev Form := (Pt → ℝ) →ₗ[ℝ] (Pt → ℝ) →ₗ[ℝ] ℝ

/-- Condition (ii): D vanishes on every locally orthogonal pair of events. -/
def Decoheres (D : Form) : Prop := ∀ e e', LO e e' → D (ind e) (ind e') = 0

/-- The indicator of "box b of party j shows outcome a". -/
def piR (j : Fin 4) (b a : Fin 2) : Pt → ℝ := fun ω => if ω j b = a then 1 else 0

theorem holds_update (j : Fin 4) (b a : Fin 2) (e : Event) (ω : Pt) :
    Holds (Function.update e.1 j b, Function.update e.2 j a) ω ↔
      (ω j b = a ∧ ∀ i, i ≠ j → ω i (e.1 i) = e.2 i) := by
  unfold Holds
  constructor
  · intro h
    refine ⟨by simpa using h j, fun i hi => ?_⟩
    have := h i
    simpa [Function.update_of_ne hi] using this
  · rintro ⟨h1, h2⟩ i
    by_cases hi : i = j
    · subst hi; simpa using h1
    · simpa [Function.update_of_ne hi] using h2 i hi

/-- (outcome a in box b of party j) · (event of the other three parties) is an event. -/
theorem piR_mul_indExcept (j : Fin 4) (b a : Fin 2) (e : Event) :
    piR j b a * indExcept j e = ind (Function.update e.1 j b, Function.update e.2 j a) := by
  funext ω
  simp only [Pi.mul_apply, piR, indExcept, ind, holds_update]
  by_cases h1 : ω j b = a <;> by_cases h2 : (∀ i, i ≠ j → ω i (e.1 i) = e.2 i) <;> simp [h1, h2]

theorem varR_eq_piR (j : Fin 4) (b : Fin 2) : varR j b = piR j b 0 - piR j b 1 := by
  funext ω
  simp only [varR, var, piR, Pi.sub_apply]
  generalize ω j b = a
  fin_cases a <;> simp [sgn]

theorem one_eq_piR (j : Fin 4) (b : Fin 2) : (1 : Pt → ℝ) = piR j b 0 + piR j b 1 := by
  funext ω
  simp only [piR, Pi.add_apply, Pi.one_apply]
  generalize ω j b = a
  fin_cases a <;> simp

/-- The heart of Step 1: opposite outcomes of one box of party j decohere, whatever the other
three parties' events. -/
theorem step1_cross (D : Form) (hii : Decoheres D) (j : Fin 4) (b a a' : Fin 2) (haa : a ≠ a')
    {F G : Pt → ℝ} (hF : F ∈ Wexc j) (hG : G ∈ Wexc j) :
    D (piR j b a * F) (piR j b a' * G) = 0 := by
  unfold Wexc at hF hG
  induction hF using Submodule.span_induction with
  | mem x hx =>
    obtain ⟨e, rfl⟩ := hx
    induction hG using Submodule.span_induction with
    | mem y hy =>
      obtain ⟨e', rfl⟩ := hy
      rw [piR_mul_indExcept, piR_mul_indExcept]
      exact hii _ _ ⟨j, by simp, by simpa using haa⟩
    | zero => simp
    | add y z _ _ hy hz => rw [mul_add, map_add, hy, hz, add_zero]
    | smul c y _ hy => rw [mul_smul_comm, map_smul, hy, smul_zero]
  | zero => simp
  | add x y _ _ hx hy => rw [mul_add, map_add, LinearMap.add_apply, hx, hy, add_zero]
  | smul c x _ hx => rw [mul_smul_comm, map_smul, LinearMap.smul_apply, hx, smul_zero]

/-- **Step 1 (rule).** For α = var(j, b) and F, G in W_{−j}: D(αF, αG) = D(F, G). -/
theorem step1 (D : Form) (hii : Decoheres D) (j : Fin 4) (b : Fin 2)
    {F G : Pt → ℝ} (hF : F ∈ Wexc j) (hG : G ∈ Wexc j) :
    D (varR j b * F) (varR j b * G) = D F G := by
  have h01 := step1_cross D hii j b 0 1 (by decide) hF hG
  have h10 := step1_cross D hii j b 1 0 (by decide) hF hG
  have eF : F = piR j b 0 * F + piR j b 1 * F := by rw [← add_mul, ← one_eq_piR, one_mul]
  have eG : G = piR j b 0 * G + piR j b 1 * G := by rw [← add_mul, ← one_eq_piR, one_mul]
  have key : D F G = D (piR j b 0 * F + piR j b 1 * F) (piR j b 0 * G + piR j b 1 * G) := by
    rw [← eF, ← eG]
  rw [key, varR_eq_piR, sub_mul, sub_mul]
  simp only [map_add, map_sub, LinearMap.add_apply, LinearMap.sub_apply, h01, h10]
  ring

/-- The outcomes of one choice of boxes partition Ω. -/
theorem sum_ind (s : Boxes) : ∑ o : Boxes, ind (s, o) = 1 := by
  funext ω
  rw [Finset.sum_apply, Finset.sum_eq_single (fun j => ω j (s j))]
  · simp [ind, Holds]
  · intro o _ ho
    have : ¬ Holds (s, o) ω := fun h => ho (funext fun j => (h j).symm)
    simp [ind, this]
  · simp

/-- D(1, 1_e) = D(1_e, 1_e), from (ii) alone. -/
theorem D_one_ind (D : Form) (hii : Decoheres D) (e : Event) : D 1 (ind e) = D (ind e) (ind e) := by
  obtain ⟨s, o⟩ := e
  rw [← sum_ind s, map_sum, LinearMap.sum_apply, Finset.sum_eq_single o]
  · intro o' _ ho
    obtain ⟨j, hj⟩ := Function.ne_iff.1 ho
    exact hii _ _ ⟨j, rfl, hj⟩
  · simp

/-- D(1_e, 1) = D(1_e, 1_e), from (ii) alone. -/
theorem D_ind_one (D : Form) (hii : Decoheres D) (e : Event) : D (ind e) 1 = D (ind e) (ind e) := by
  obtain ⟨s, o⟩ := e
  rw [← sum_ind s, map_sum, Finset.sum_eq_single o]
  · intro o' _ ho
    obtain ⟨j, hj⟩ := Function.ne_iff.1 ho
    exact hii _ _ ⟨j, rfl, fun h => hj h.symm⟩
  · simp

/-- D(1, 1) is the total of the diagonal over the outcomes of any one choice of boxes. -/
theorem D_one_one (D : Form) (hii : Decoheres D) (s : Boxes) :
    D 1 1 = ∑ o : Boxes, D (ind (s, o)) (ind (s, o)) := by
  calc D 1 1 = D 1 (∑ o : Boxes, ind (s, o)) := by rw [sum_ind]
    _ = ∑ o : Boxes, D 1 (ind (s, o)) := map_sum _ _ _
    _ = ∑ o : Boxes, D (ind (s, o)) (ind (s, o)) :=
        Finset.sum_congr rfl fun o _ => D_one_ind D hii (s, o)

/-- A function of the four outcomes read in the boxes s is a combination of the 16 events with
boxes s. -/
theorem boxfun_eq (s : Boxes) (g : Boxes → ℝ) :
    (fun ω : Pt => g (fun i => ω i (s i))) = ∑ o : Boxes, g o • ind (s, o) := by
  funext ω
  rw [Finset.sum_apply, Finset.sum_eq_single (fun i => ω i (s i))]
  · simp [ind, Holds]
  · intro o _ ho
    have : ¬ Holds (s, o) ω := fun h => ho (funext fun j => (h j).symm)
    simp [ind, this]
  · simp

theorem boxfun_mem_W (s : Boxes) (g : Boxes → ℝ) : (fun ω : Pt => g (fun i => ω i (s i))) ∈ W := by
  rw [boxfun_eq]
  exact Submodule.sum_mem _ fun o _ => Submodule.smul_mem _ _ (Submodule.subset_span ⟨_, rfl⟩)

/-- A function of the outcomes read in the boxes s that does not depend on party j's outcome
lies in W_{−j}. -/
theorem boxfun_mem_Wexc (j : Fin 4) (s : Boxes) (g : Boxes → ℝ)
    (hg : ∀ o a, g (Function.update o j a) = g o) :
    (fun ω : Pt => g (fun i => ω i (s i))) ∈ Wexc j := by
  have key : (fun ω : Pt => g (fun i => ω i (s i))) =
      ∑ o : Boxes with o j = 0, g o • indExcept j (s, o) := by
    funext ω
    rw [Finset.sum_apply, Finset.sum_eq_single (Function.update (fun i => ω i (s i)) j 0)]
    · have h : ∀ i, i ≠ j → ω i (s i) = Function.update (fun i => ω i (s i)) j 0 i :=
        fun i hi => by rw [Function.update_of_ne hi]
      have h2 : g (Function.update (fun i => ω i (s i)) j 0) = g (fun i => ω i (s i)) := hg _ _
      simp only [indExcept, Pi.smul_apply, smul_eq_mul]
      rw [if_pos h, h2, mul_one]
    · intro o ho hne
      have : ¬ ∀ i, i ≠ j → ω i (s i) = o i := by
        intro h; apply hne; funext i
        by_cases hi : i = j
        · subst hi; simpa using (Finset.mem_filter.1 ho).2
        · rw [Function.update_of_ne hi]; exact (h i hi).symm
      simp [indExcept, this]
    · intro h; exact absurd (by simp) h
  rw [key]
  exact Submodule.sum_mem _ fun o _ => Submodule.smul_mem _ _ (Submodule.subset_span ⟨_, rfl⟩)

end TwoCopy
