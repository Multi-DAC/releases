/-
Two-copy identity, file 6 of 6: Corollary 2 and the "only if" half of Theorem 3
(targets (d) and (e)), and a model that shows Theorem 1's hypotheses can be met.

What this file proves, in plain words.

  * `corollary2`: if P sums to 1 over the outcomes of one choice of boxes, and D is a bilinear
    form with (i), (ii) and (iii) D(1_T, 1_T) ≥ 0 for every subset T of Ω, then
    −8 ≤ ⟨C₁C₂⟩_P ≤ 8. D is not assumed symmetric. The sharper forms say which sets are used:
    `corollary2_upper` needs (iii) only on the four sets S_k, `corollary2_lower` only on their
    four complements.
  * The lower bound is proved here from the complements (`theorem1_compl`), not by the
    relabeling the draft alludes to. `corollary2_lower_by_relabeling` then carries out the
    relabeling argument as well (flip both outcomes of party A₁, which changes the sign of C₁),
    so that route goes through too.
  * `corr_PRq`: for the product of two isotropic boxes, ⟨C₁C₂⟩ = 16 q².
    `PRq_norm`: that behaviour sums to 1 in every context.
  * `theorem3_only_if`: for q² > 1/2 there is no bilinear D (symmetric or not) with (i), (ii),
    (iii) for PR_q ⊗ PR_q. `theorem3_only_if_sqrt`: the same stated as q > 1/√2.
    `theorem3_JQM`: if PR_q ⊗ PR_q has a JQM (symmetric D with (i), (ii), (iii)) then
    −1/√2 ≤ q ≤ 1/√2.
  * Tests: `corr_pointmass` (for the behaviour concentrated at a point ω, ⟨C₁C₂⟩ is
    C₁(ω)C₂(ω): ties `corr` and `E` to `C`), `PRq_at_zero`, and `evalForm_hyps` (the form
    D(F, G) = F(ω)G(ω) is symmetric and satisfies (i), (ii), (iii) for the point-mass behaviour,
    so the hypotheses of Theorem 1 and Corollary 2 are satisfiable, and Theorem 1 for it is the
    pointwise identity).
-/
import TwoCopy.Identity

namespace TwoCopy

open Finset

/-! ### Corollary 2 -/

/-- Upper bound: (iii) is needed only on the four sets S_k. -/
theorem corollary2_upper (P : Event → ℝ) (D : Form) (hi : ∀ e, D (ind e) (ind e) = P e)
    (hii : ∀ e e', LO e e' → D (ind e) (ind e') = 0)
    (s₀ : Boxes) (hnorm : ∑ o : Boxes, P (s₀, o) = 1)
    (hpos : ∀ k, 0 ≤ D (indF (Sset k)) (indF (Sset k))) : corr P ≤ 8 := by
  have h := theorem1 P D hi hii s₀ hnorm
  have h0 : 0 ≤ ∑ k, D (indF (Sset k)) (indF (Sset k)) := Finset.sum_nonneg fun k _ => hpos k
  linarith

/-- Lower bound: (iii) is needed only on the four complements of the S_k. -/
theorem corollary2_lower (P : Event → ℝ) (D : Form) (hi : ∀ e, D (ind e) (ind e) = P e)
    (hii : ∀ e e', LO e e' → D (ind e) (ind e') = 0)
    (s₀ : Boxes) (hnorm : ∑ o : Boxes, P (s₀, o) = 1)
    (hpos : ∀ k, 0 ≤ D (indF (Sset k)ᶜ) (indF (Sset k)ᶜ)) : -8 ≤ corr P := by
  have h := theorem1_compl P D hi hii
  rw [D_one_one_eq P D hi hii s₀, hnorm] at h
  have h0 : 0 ≤ ∑ k, D (indF (Sset k)ᶜ) (indF (Sset k)ᶜ) := Finset.sum_nonneg fun k _ => hpos k
  linarith

/-- **Corollary 2.** A normalised P with a bilinear D satisfying (i), (ii) and
(iii) D(1_T, 1_T) ≥ 0 for every T ⊆ Ω has −8 ≤ ⟨C₁C₂⟩_P ≤ 8. -/
theorem corollary2 (P : Event → ℝ) (D : Form) (hi : ∀ e, D (ind e) (ind e) = P e)
    (hii : ∀ e e', LO e e' → D (ind e) (ind e') = 0)
    (s₀ : Boxes) (hnorm : ∑ o : Boxes, P (s₀, o) = 1)
    (hiii : ∀ T : Finset Pt, 0 ≤ D (indF T) (indF T)) : -8 ≤ corr P ∧ corr P ≤ 8 :=
  ⟨corollary2_lower P D hi hii s₀ hnorm fun _ => hiii _,
    corollary2_upper P D hi hii s₀ hnorm fun _ => hiii _⟩

/-! ### The lower bound again, by the draft's relabeling -/

/-- The relabeling of Ω that flips both outcomes of party 0 (A₁); it changes the sign of C₁. -/
def flip0 (ω : Pt) : Pt := fun i b => if i = 0 then ω i b + 1 else ω i b

/-- The same relabeling on events: flip the outcome of party 0. -/
def flipE (e : Event) : Event := (e.1, fun i => if i = 0 then e.2 i + 1 else e.2 i)

set_option maxRecDepth 100000 in
theorem flip0_invol : ∀ ω, flip0 (flip0 ω) = ω := by decide +kernel

set_option maxRecDepth 100000 in
theorem flipE_invol : ∀ e, flipE (flipE e) = e := by decide +kernel

set_option maxRecDepth 100000 in
theorem C_flip0 : ∀ ω, C 0 (flip0 ω) = -C 0 ω ∧ C 1 (flip0 ω) = C 1 ω := by decide +kernel

set_option maxRecDepth 100000 in
theorem holds_flip : ∀ (e : Event) (ω : Pt), Holds e (flip0 ω) ↔ Holds (flipE e) ω := by
  decide +kernel

set_option maxRecDepth 100000 in
theorem LO_flip : ∀ e e' : Event, LO e e' → LO (flipE e) (flipE e') := by
  rintro e e' ⟨j, h1, h2⟩
  refine ⟨j, h1, ?_⟩
  simp only [flipE]
  by_cases hj : j = 0 <;> simp [hj] <;> simpa [hj] using h2

theorem ind_comp_flip (e : Event) : ind e ∘ flip0 = ind (flipE e) := by
  funext ω; simp only [Function.comp, ind, holds_flip]

theorem indF_comp_flip (T : Finset Pt) : indF T ∘ flip0 = indF (T.image flip0) := by
  funext ω
  have : flip0 ω ∈ T ↔ ω ∈ T.image flip0 := by
    rw [Finset.mem_image]
    constructor
    · intro h; exact ⟨flip0 ω, h, flip0_invol ω⟩
    · rintro ⟨ω', h, rfl⟩; rw [flip0_invol]; exact h
  simp only [Function.comp, indF, this]

/-- The outcome flip, on outcome tuples, as a bijection. -/
def flipO : Boxes ≃ Boxes where
  toFun o := fun i => if i = 0 then o i + 1 else o i
  invFun o := fun i => if i = 0 then o i + 1 else o i
  left_inv := by decide
  right_inv := by decide

set_option maxRecDepth 100000 in
theorem sign_flipO : ∀ o : Boxes,
    (-1 : ℤ) ^ ((flipO o 0).val + (flipO o 1).val + (flipO o 2).val + (flipO o 3).val) =
      -(-1 : ℤ) ^ ((o 0).val + (o 1).val + (o 2).val + (o 3).val) := by decide +kernel

theorem neg_one_pow_flipO (o : Boxes) :
    (-1 : ℝ) ^ ((flipO o 0).val + (flipO o 1).val + (flipO o 2).val + (flipO o 3).val) =
      -(-1 : ℝ) ^ ((o 0).val + (o 1).val + (o 2).val + (o 3).val) := by
  exact_mod_cast congrArg (fun z : ℤ => (z : ℝ)) (sign_flipO o)

/-- The relabeled behaviour has the opposite correlator. -/
theorem corr_flip (P : Event → ℝ) : corr (P ∘ flipE) = -corr P := by
  unfold corr
  rw [← Finset.sum_neg_distrib]
  refine Finset.sum_congr rfl fun s _ => ?_
  have hE : E (P ∘ flipE) s = -E P s := by
    unfold E
    rw [← Finset.sum_neg_distrib, ← Equiv.sum_comp flipO
      (fun o => -((-1 : ℝ) ^ ((o 0).val + (o 1).val + (o 2).val + (o 3).val) * P (s, o)))]
    refine Finset.sum_congr rfl fun o _ => ?_
    rw [neg_one_pow_flipO]
    have : (P ∘ flipE) (s, o) = P (s, flipO o) := rfl
    rw [this]; ring
  rw [hE]; ring

theorem norm_flip (P : Event → ℝ) (s : Boxes) :
    ∑ o : Boxes, (P ∘ flipE) (s, o) = ∑ o : Boxes, P (s, o) := by
  rw [← Equiv.sum_comp flipO (fun o => P (s, o))]
  exact Finset.sum_congr rfl fun o _ => rfl

/-- **The draft's route to the lower bound.** Relabel by flipping the outcomes of party A₁:
D′(F, G) = D(F ∘ flip, G ∘ flip) satisfies (i) for the relabeled behaviour, (ii) and (iii), and
the relabeled behaviour has correlator −⟨C₁C₂⟩_P; the upper bound for it is the lower bound
for P. The relabeling argument goes through as the draft states it. -/
theorem corollary2_lower_by_relabeling (P : Event → ℝ) (D : Form)
    (hi : ∀ e, D (ind e) (ind e) = P e)
    (hii : ∀ e e', LO e e' → D (ind e) (ind e') = 0)
    (s₀ : Boxes) (hnorm : ∑ o : Boxes, P (s₀, o) = 1)
    (hiii : ∀ T : Finset Pt, 0 ≤ D (indF T) (indF T)) : -8 ≤ corr P := by
  let L : (Pt → ℝ) →ₗ[ℝ] (Pt → ℝ) := LinearMap.funLeft ℝ ℝ flip0
  have hL : ∀ F : Pt → ℝ, L F = F ∘ flip0 := fun F => rfl
  have h := corollary2_upper (P ∘ flipE) (D.compl₁₂ L L)
    (fun e => by
      rw [LinearMap.compl₁₂_apply, hL, ind_comp_flip, hi]; rfl)
    (fun e e' he => by
      rw [LinearMap.compl₁₂_apply, hL, hL, ind_comp_flip, ind_comp_flip]
      exact hii _ _ (LO_flip e e' he))
    s₀ (by rw [norm_flip, hnorm])
    (fun k => by
      rw [LinearMap.compl₁₂_apply, hL, indF_comp_flip]; exact hiii _)
  rw [corr_flip] at h
  linarith

/-! ### Theorem 3, "only if" -/

/-- +1 if the copy-1 part of the event wins the PR game (a xor b = x·y), −1 otherwise. -/
def good1 (e : Event) : ℤ := if e.2 0 + e.2 1 = e.1 0 * e.1 1 then 1 else -1

/-- The same for the copy-2 part. -/
def good2 (e : Event) : ℤ := if e.2 2 + e.2 3 = e.1 2 * e.1 3 then 1 else -1

/-- The sign with which P(e) enters ⟨C₁C₂⟩_P. -/
def sig (e : Event) : ℤ :=
  (-1) ^ ((e.1 0).val * (e.1 1).val + (e.1 2).val * (e.1 3).val) *
    (-1) ^ ((e.2 0).val + (e.2 1).val + (e.2 2).val + (e.2 3).val)

theorem PRq_eq (q : ℝ) (e : Event) :
    PRq q e = (1 + q * (good1 e : ℝ)) * (1 + q * (good2 e : ℝ)) / 16 := by
  by_cases h1 : e.2 0 + e.2 1 = e.1 0 * e.1 1 <;> by_cases h2 : e.2 2 + e.2 3 = e.1 2 * e.1 3 <;>
    simp [PRq, pq, good1, good2, h1, h2] <;> ring

theorem corr_eq_sum (P : Event → ℝ) : corr P = ∑ e : Event, (sig e : ℝ) * P e := by
  unfold corr E
  rw [Fintype.sum_prod_type]
  refine Finset.sum_congr rfl fun s _ => ?_
  rw [Finset.mul_sum]
  refine Finset.sum_congr rfl fun o _ => ?_
  simp only [sig]; push_cast; ring

set_option maxRecDepth 100000 in
theorem sig_sums :
    (∑ e : Event, sig e) = 0 ∧ (∑ e : Event, sig e * good1 e) = 0 ∧
      (∑ e : Event, sig e * good2 e) = 0 ∧ (∑ e : Event, sig e * good1 e * good2 e) = 256 := by
  decide +kernel

set_option maxRecDepth 100000 in
theorem good_sums : ∀ s : Boxes,
    (∑ o : Boxes, good1 (s, o)) = 0 ∧ (∑ o : Boxes, good2 (s, o)) = 0 ∧
      (∑ o : Boxes, good1 (s, o) * good2 (s, o)) = 0 := by
  decide +kernel

/-- For the product of two isotropic boxes, ⟨C₁C₂⟩ = 16 q². -/
theorem corr_PRq (q : ℝ) : corr (PRq q) = 16 * q ^ 2 := by
  obtain ⟨h0, h1, h2, h12⟩ := sig_sums
  have r0 : ∑ e : Event, (sig e : ℝ) = 0 := by exact_mod_cast congrArg (fun z : ℤ => (z : ℝ)) h0
  have r1 : ∑ e : Event, (sig e : ℝ) * (good1 e : ℝ) = 0 := by
    exact_mod_cast congrArg (fun z : ℤ => (z : ℝ)) h1
  have r2 : ∑ e : Event, (sig e : ℝ) * (good2 e : ℝ) = 0 := by
    exact_mod_cast congrArg (fun z : ℤ => (z : ℝ)) h2
  have r12 : ∑ e : Event, (sig e : ℝ) * (good1 e : ℝ) * (good2 e : ℝ) = 256 := by
    exact_mod_cast congrArg (fun z : ℤ => (z : ℝ)) h12
  have e : ∀ ev : Event, (sig ev : ℝ) * PRq q ev =
      (1 / 16) * (sig ev : ℝ) + (q / 16) * ((sig ev : ℝ) * (good1 ev : ℝ)) +
        (q / 16) * ((sig ev : ℝ) * (good2 ev : ℝ)) +
        (q ^ 2 / 16) * ((sig ev : ℝ) * (good1 ev : ℝ) * (good2 ev : ℝ)) := by
    intro ev; rw [PRq_eq]; ring
  rw [corr_eq_sum]
  simp only [e]
  rw [Finset.sum_add_distrib, Finset.sum_add_distrib, Finset.sum_add_distrib,
    ← Finset.mul_sum, ← Finset.mul_sum, ← Finset.mul_sum, ← Finset.mul_sum, r0, r1, r2, r12]
  ring

/-- The product of two isotropic boxes sums to 1 over the outcomes of every choice of boxes. -/
theorem PRq_norm (q : ℝ) (s : Boxes) : ∑ o : Boxes, PRq q (s, o) = 1 := by
  obtain ⟨h1, h2, h12⟩ := good_sums s
  have r0 : ∑ _o : Boxes, (1 : ℝ) = 16 := by simp
  have r1 : ∑ o : Boxes, (good1 (s, o) : ℝ) = 0 := by
    exact_mod_cast congrArg (fun z : ℤ => (z : ℝ)) h1
  have r2 : ∑ o : Boxes, (good2 (s, o) : ℝ) = 0 := by
    exact_mod_cast congrArg (fun z : ℤ => (z : ℝ)) h2
  have r12 : ∑ o : Boxes, (good1 (s, o) : ℝ) * (good2 (s, o) : ℝ) = 0 := by
    exact_mod_cast congrArg (fun z : ℤ => (z : ℝ)) h12
  have e : ∀ o : Boxes, PRq q (s, o) =
      (1 / 16) * (1 : ℝ) + (q / 16) * (good1 (s, o) : ℝ) + (q / 16) * (good2 (s, o) : ℝ) +
        (q ^ 2 / 16) * ((good1 (s, o) : ℝ) * (good2 (s, o) : ℝ)) := by
    intro o; rw [PRq_eq]; ring
  simp only [e]
  rw [Finset.sum_add_distrib, Finset.sum_add_distrib, Finset.sum_add_distrib,
    ← Finset.mul_sum, ← Finset.mul_sum, ← Finset.mul_sum, ← Finset.mul_sum, r0, r1, r2, r12]
  ring

/-- **Theorem 3, "only if".** For q² > 1/2 there is no bilinear form D, symmetric or not, with
(i) D(1_e, 1_e) = PR_q ⊗ PR_q (e), (ii) D(1_e, 1_e′) = 0 on locally orthogonal pairs, and
(iii) D(1_T, 1_T) ≥ 0 for all T ⊆ Ω. -/
theorem theorem3_only_if (q : ℝ) (hq : 1 / 2 < q ^ 2) :
    ¬ ∃ D : Form, (∀ e, D (ind e) (ind e) = PRq q e) ∧
      (∀ e e', LO e e' → D (ind e) (ind e') = 0) ∧
      (∀ T : Finset Pt, 0 ≤ D (indF T) (indF T)) := by
  rintro ⟨D, hi, hii, hiii⟩
  have h := (corollary2 (PRq q) D hi hii 0 (PRq_norm q 0) hiii).2
  rw [corr_PRq] at h
  linarith

/-- The same with the draft's threshold: q > 1/√2. -/
theorem theorem3_only_if_sqrt (q : ℝ) (hq : 1 / Real.sqrt 2 < q) :
    ¬ ∃ D : Form, (∀ e, D (ind e) (ind e) = PRq q e) ∧
      (∀ e e', LO e e' → D (ind e) (ind e') = 0) ∧
      (∀ T : Finset Pt, 0 ≤ D (indF T) (indF T)) := by
  apply theorem3_only_if
  have h2 : (1 / Real.sqrt 2) ^ 2 = 1 / 2 := by
    rw [div_pow, Real.sq_sqrt (by norm_num)]; norm_num
  have hpos : 0 < 1 / Real.sqrt 2 := by positivity
  nlinarith

/-- A joint quantum measure for P, as the draft defines it: a symmetric bilinear form with
(i), (ii), (iii). -/
def IsJQM (P : Event → ℝ) (D : Form) : Prop :=
  (∀ F G, D F G = D G F) ∧ (∀ e, D (ind e) (ind e) = P e) ∧
    (∀ e e', LO e e' → D (ind e) (ind e') = 0) ∧ (∀ T : Finset Pt, 0 ≤ D (indF T) (indF T))

/-- **Theorem 3, "only if", in the draft's words.** If PR_q ⊗ PR_q has a JQM then q ≤ 1/√2
(and q ≥ −1/√2). -/
theorem theorem3_JQM (q : ℝ) (h : ∃ D : Form, IsJQM (PRq q) D) :
    -(1 / Real.sqrt 2) ≤ q ∧ q ≤ 1 / Real.sqrt 2 := by
  obtain ⟨D, -, hi, hii, hiii⟩ := h
  have hq : q ^ 2 ≤ 1 / 2 := by
    by_contra hlt
    exact theorem3_only_if q (not_le.1 hlt) ⟨D, hi, hii, hiii⟩
  have h2 : (1 / Real.sqrt 2) ^ 2 = 1 / 2 := by
    rw [div_pow, Real.sq_sqrt (by norm_num)]; norm_num
  have hpos : 0 < 1 / Real.sqrt 2 := by positivity
  constructor <;> nlinarith

/-! ### Tests and a model -/

/-- The behaviour concentrated at the point ω. -/
def pointP (ω : Pt) (e : Event) : ℝ := if Holds e ω then 1 else 0

set_option maxRecDepth 100000 in
theorem corr_pointmassZ : ∀ ω : Pt,
    (∑ e : Event, sig e * (if Holds e ω then 1 else 0)) = C 0 ω * C 1 ω := by
  decide +kernel

/-- Test of `corr`, `E` against `C`: for the behaviour concentrated at ω, ⟨C₁C₂⟩ = C₁(ω)C₂(ω). -/
theorem corr_pointmass (ω : Pt) : corr (pointP ω) = (C 0 ω : ℝ) * (C 1 ω : ℝ) := by
  rw [corr_eq_sum]
  have h := congrArg (fun z : ℤ => (z : ℝ)) (corr_pointmassZ ω)
  simp only [pointP]
  push_cast at h
  exact h

/-- Test of `PRq`: the event "all boxes 0, all outcomes 0" is good on both copies. -/
theorem PRq_at_zero (q : ℝ) : PRq q (0, 0) = ((1 + q) / 4) ^ 2 := by
  simp [PRq, pq]; ring

/-- The form D(F, G) = F(ω)G(ω). -/
def evalForm (ω : Pt) : Form :=
  LinearMap.mk₂ ℝ (fun F G => F ω * G ω)
    (fun F F' G => by simp [add_mul]) (fun c F G => by simp [mul_assoc])
    (fun F G G' => by simp [mul_add]) (fun c F G => by simp [mul_left_comm])

/-- The hypotheses of Theorem 1 and Corollary 2 are satisfiable: evaluation at a point is a
symmetric form with (i) for the point-mass behaviour, (ii) and (iii); the behaviour is
normalised. So `evalForm ω` is a JQM of `pointP ω`. -/
theorem evalForm_hyps (ω : Pt) :
    IsJQM (pointP ω) (evalForm ω) ∧ ∀ s : Boxes, ∑ o : Boxes, pointP ω (s, o) = 1 := by
  refine ⟨⟨fun F G => by simp [evalForm, mul_comm], fun e => ?_, fun e e' he => ?_, fun T => ?_⟩,
    fun s => ?_⟩
  · by_cases h : Holds e ω <;> simp [evalForm, ind, pointP, h]
  · have := LO_disjoint e e' he ω
    by_cases h : Holds e ω <;> by_cases h' : Holds e' ω <;> simp [evalForm, ind, h, h']
    exact this ⟨h, h'⟩
  · simp only [evalForm, LinearMap.mk₂_apply]; exact mul_self_nonneg _
  · have := congrFun (sum_ind s) ω
    simpa [Finset.sum_apply, ind, pointP] using this

end TwoCopy
