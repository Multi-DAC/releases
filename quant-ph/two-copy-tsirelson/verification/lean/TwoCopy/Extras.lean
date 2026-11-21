/-
Two-copy identity: extras. Checks of side claims of the draft, a test of the definition of W,
and Theorem 3 "only if" in its sharpest form. Nothing in the other files depends on this one.

What this file proves, in plain words.

  * `image_tau2`: the image of S under τ₂ is the same set as its image under τ₁ (so the four
    named images S, τ₀S, τ₁S, τ₃S are the images under the identity and all four generators).
  * `step3_survivors`, `step3_survivors_nonzero`: the draft's Step 3 for the single set S. In
    D(M, M) the difference patterns d whose coefficient does not vanish are among: none,
    {B₁, A₂, B₂}, {A₁, B₁, A₂}, all four; and each of those three does survive. So for a single
    set the monomial reduction leaves unknowns (the draft removes them by invariance, Step 4;
    the four-set sum removes them by cancellation, `Finite.coeff`).
  * `same_party_product_not_mem_W`: the product of the two box variables of one party is not
    in W. A test of the definition of W: it is not the whole function space (the draft:
    dim W = 81 < 256; the dimension itself is not proved here).
  * `theorem3_only_if_W`: for q² > 1/2 there is no bilinear form on W alone with (i), (ii) for
    PR_q ⊗ PR_q and non-negative on the four sets S_k. This implies `theorem3_only_if` and uses
    the least: a form on W, not symmetric, positivity on four sets.
-/
import TwoCopy.Bounds

namespace TwoCopy

open Finset

/-! ### τ₂ S = τ₁ S -/

set_option maxRecDepth 100000 in
theorem tau2_tau1 : ∀ ω, (M (tau 2 ω) = -4 ↔ M (tau 1 ω) = -4) := by decide +kernel

theorem image_tau2 : S.image (tau 2) = Sset 2 := by
  ext ω
  rw [mem_Sset, Finset.mem_image]
  constructor
  · rintro ⟨ω', h, rfl⟩
    exact (tau2_tau1 (tau 2 ω')).1 (by rw [tau_invol]; exact (mem_S _).1 h)
  · intro h
    exact ⟨tau 2 ω, (mem_S _).2 ((tau2_tau1 ω).2 h), tau_invol 2 ω⟩

/-! ### The draft's Step 3, for the single set -/

set_option maxRecDepth 100000 in
/-- In D(M, M) = Σ_{s,s′} eps s · eps s′ · D(m_s, m_s′), the difference patterns with a
non-vanishing coefficient are among ∅, {B₁, A₂, B₂}, {A₁, B₁, A₂} and all four parties. -/
theorem step3_survivors : ∀ d t : Boxes, (∑ s : Boxes with s * d = t, eps s * eps (s + d)) ≠ 0 →
    (d = 0 ∨ d = ![0, 1, 1, 1] ∨ d = ![1, 1, 1, 0] ∨ d = ![1, 1, 1, 1]) := by
  decide +kernel

set_option maxRecDepth 100000 in
/-- Each of the three non-empty patterns does survive, for some restriction t. -/
theorem step3_survivors_nonzero :
    ∀ d : Boxes, (d = ![0, 1, 1, 1] ∨ d = ![1, 1, 1, 0] ∨ d = ![1, 1, 1, 1]) →
      ∃ t : Boxes, (∑ s : Boxes with s * d = t, eps s * eps (s + d)) ≠ 0 := by
  decide +kernel

/-! ### A test of W -/

/-- The product of the two box variables of party 0, as an integer function. -/
def cc (ω : Pt) : ℤ := var 0 0 ω * var 0 1 ω

/-- The linear functional F ↦ Σ_ω cc(ω) F(ω). -/
def phiW : (Pt → ℝ) →ₗ[ℝ] ℝ := ∑ ω : Pt, (cc ω : ℝ) • LinearMap.proj ω

theorem phiW_apply (F : Pt → ℝ) : phiW F = ∑ ω : Pt, (cc ω : ℝ) * F ω := by
  simp [phiW, LinearMap.sum_apply]

set_option maxRecDepth 100000 in
theorem cc_ind : ∀ e : Event, (∑ ω : Pt, cc ω * (if Holds e ω then 1 else 0)) = 0 := by
  decide +kernel

set_option maxRecDepth 100000 in
theorem cc_sq : ∀ ω : Pt, cc ω * cc ω = 1 := by decide +kernel

theorem W_le_ker : W ≤ LinearMap.ker phiW := by
  unfold W
  rw [Submodule.span_le]
  rintro _ ⟨e, rfl⟩
  rw [SetLike.mem_coe, LinearMap.mem_ker, phiW_apply]
  have h := congrArg (fun z : ℤ => (z : ℝ)) (cc_ind e)
  simp only [ind]
  push_cast at h
  exact h

/-- The product of the two box variables of one party is not in W: W is a proper subspace. -/
theorem same_party_product_not_mem_W : varR 0 0 * varR 0 1 ∉ W := by
  intro h
  have h0 := LinearMap.mem_ker.1 (W_le_ker h)
  rw [phiW_apply] at h0
  have h1 : ∀ ω : Pt, (cc ω : ℝ) * (varR 0 0 * varR 0 1) ω = 1 := by
    intro ω
    have := congrArg (fun z : ℤ => (z : ℝ)) (cc_sq ω)
    simp only [Pi.mul_apply, varR]
    simp only [cc] at this ⊢
    push_cast at this ⊢
    exact this
  simp only [h1] at h0
  simp at h0

/-! ### Theorem 3 "only if", sharpest form -/

/-- For q² > 1/2 there is no bilinear form on W, symmetric or not, with (i) and (ii) for
PR_q ⊗ PR_q that is non-negative on the four sets S_k. -/
theorem theorem3_only_if_W (q : ℝ) (hq : 1 / 2 < q ^ 2) :
    ¬ ∃ D : W →ₗ[ℝ] W →ₗ[ℝ] ℝ,
      (∀ e, D ⟨ind e, ind_mem_W e⟩ ⟨ind e, ind_mem_W e⟩ = PRq q e) ∧
      (∀ e e', LO e e' → D ⟨ind e, ind_mem_W e⟩ ⟨ind e', ind_mem_W e'⟩ = 0) ∧
      (∀ k, 0 ≤ D ⟨indF (Sset k), indF_Sset_mem_W k⟩ ⟨indF (Sset k), indF_Sset_mem_W k⟩) := by
  rintro ⟨D, hi, hii, hpos⟩
  have h := theorem1_W (PRq q) D hi hii 0 (PRq_norm q 0)
  have h0 : 0 ≤ ∑ k, D ⟨indF (Sset k), indF_Sset_mem_W k⟩ ⟨indF (Sset k), indF_Sset_mem_W k⟩ :=
    Finset.sum_nonneg fun k _ => hpos k
  rw [corr_PRq] at h
  linarith

#print axioms image_tau2
#print axioms step3_survivors
#print axioms step3_survivors_nonzero
#print axioms same_party_product_not_mem_W
#print axioms theorem3_only_if_W

end TwoCopy
