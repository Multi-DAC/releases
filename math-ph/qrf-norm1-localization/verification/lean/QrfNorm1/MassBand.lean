import QrfNorm1.Basic

/-!
# QrfNorm1.MassBand

Formal verification of the exact vertical-slice inequality at the heart of
Proposition 5.1(i) (`prop:mass-band`): for all `r` and `0 ≤ m₁ ≤ m₂`,
  `√(r² + m₂²) - √(r² + m₁²) ≤ m₂ - m₁`,
with equality if and only if `r = 0 ∨ m₁ = m₂`.

Proved purely algebraically over any `OrderedCommRing` (covering `Int`, `Rat`, `Real`).
-/

namespace QrfNorm1

section MassBand
variable {R : Type} [OrderedCommRing R]

/-- Algebraic cross-product identity:
    `(b m₂)² - (a m₁)² = r² (m₂ - m₁)(m₂ + m₁)` whenever `a² = r² + m₂²`, `b² = r² + m₁²`. -/
private theorem slice_cross_sq_identity (r m1 m2 a b : R)
    (ha2 : a * a = r * r + m2 * m2) (hb2 : b * b = r * r + m1 * m1) :
    (b * m2) * (b * m2) - (a * m1) * (a * m1) =
      (r * r) * ((m2 - m1) * (m2 + m1)) := by
  have h1 : (b * m2) * (b * m2) = (b * b) * (m2 * m2) := by ac_rfl
  have h2 : (a * m1) * (a * m1) = (a * a) * (m1 * m1) := by ac_rfl
  rw [h1, h2, hb2, ha2, right_distrib, right_distrib]
  have hcancel :
      (r * r * (m2 * m2) + m1 * m1 * (m2 * m2)) - (r * r * (m1 * m1) + m2 * m2 * (m1 * m1)) =
      r * r * (m2 * m2) - r * r * (m1 * m1) := by
    have hcomm : m1 * m1 * (m2 * m2) = m2 * m2 * (m1 * m1) := OrderedCommRing.mul_comm _ _
    rw [hcomm, OrderedCommRing.sub_eq_add_neg, OrderedCommRing.sub_eq_add_neg, neg_add]
    calc (r * r * (m2 * m2) + m2 * m2 * (m1 * m1)) +
            (-(r * r * (m1 * m1)) + -(m2 * m2 * (m1 * m1)))
        = (r * r * (m2 * m2) + -(r * r * (m1 * m1))) +
            (m2 * m2 * (m1 * m1) + -(m2 * m2 * (m1 * m1))) := by ac_rfl
      _ = (r * r * (m2 * m2) + -(r * r * (m1 * m1))) + 0 := by rw [add_neg_cancel]
      _ = r * r * (m2 * m2) + -(r * r * (m1 * m1)) := add_zero _
  rw [hcancel, ← left_distrib_sub]
  congr 1
  -- m2*m2 - m1*m1 = (m2 - m1)*(m2 + m1)
  rw [right_distrib_sub, OrderedCommRing.left_distrib, OrderedCommRing.left_distrib,
      OrderedCommRing.mul_comm m1 m2, OrderedCommRing.sub_eq_add_neg,
      OrderedCommRing.sub_eq_add_neg, neg_add]
  calc m2 * m2 + -(m1 * m1)
      = m2 * m2 + -(m1 * m1) + 0 := by rw [add_zero]
    _ = m2 * m2 + -(m1 * m1) + (m2 * m1 + -(m2 * m1)) := by rw [add_neg_cancel]
    _ = (m2 * m2 + m2 * m1) + (-(m2 * m1) + -(m1 * m1)) := by ac_rfl

/-- Algebraic square identity:
    `(b + m₂)² - (a + m₁)² = 2 (b m₂ - a m₁)` whenever `a² = r² + m₂²`, `b² = r² + m₁²`. -/
private theorem slice_sum_sq_identity (r m1 m2 a b : R)
    (ha2 : a * a = r * r + m2 * m2) (hb2 : b * b = r * r + m1 * m1) :
    (b + m2) * (b + m2) - (a + m1) * (a + m1) =
      (b * m2 + b * m2) - (a * m1 + a * m1) := by
  have hexp : ∀ x y : R, (x + y) * (x + y) = (x * x + y * y) + (x * y + x * y) := by
    intro x y
    rw [right_distrib, OrderedCommRing.left_distrib, OrderedCommRing.left_distrib,
        OrderedCommRing.mul_comm y x]
    ac_rfl
  rw [hexp b m2, hexp a m1, hb2, ha2]
  have hC : (r * r + m1 * m1) + m2 * m2 = (r * r + m2 * m2) + m1 * m1 := by ac_rfl
  rw [hC, add_sub_add_left]

/-- Proposition 5.1(i) sharp vertical-slice bound:
    if `0 ≤ m₁ ≤ m₂`, `0 ≤ a`, `0 ≤ b`, `a² = r² + m₂²`, `b² = r² + m₁²`,
    then `a - b ≤ m₂ - m₁`. -/
theorem vertical_slice_bound (r m1 m2 a b : R)
    (hm1 : (0 : R) ≤ m1) (hm : m1 ≤ m2) (ha : (0 : R) ≤ a) (hb : (0 : R) ≤ b)
    (ha2 : a * a = r * r + m2 * m2) (hb2 : b * b = r * r + m1 * m1) :
    a - b ≤ m2 - m1 := by
  have hm2 : (0 : R) ≤ m2 := OrderedCommRing.le_trans hm1 hm
  have hcross_sq : (a * m1) * (a * m1) ≤ (b * m2) * (b * m2) := by
    apply le_of_sub_nonneg
    rw [slice_cross_sq_identity r m1 m2 a b ha2 hb2]
    have hsum0 : (0 : R) ≤ m2 + m1 := by
      have := add_le_add hm2 hm1
      rwa [add_zero] at this
    exact mul_nonneg (sq_nonneg r)
      (mul_nonneg (sub_nonneg_of_le hm) hsum0)
  have hcross : a * m1 ≤ b * m2 :=
    le_of_sq_le_sq (mul_nonneg ha hm1) (mul_nonneg hb hm2) hcross_sq
  have hsum_sq : (a + m1) * (a + m1) ≤ (b + m2) * (b + m2) := by
    apply le_of_sub_nonneg
    rw [slice_sum_sq_identity r m1 m2 a b ha2 hb2]
    exact sub_nonneg_of_le (add_le_add hcross hcross)
  have ham1 : (0 : R) ≤ a + m1 := by
    have := add_le_add ha hm1; rwa [add_zero] at this
  have hbm2 : (0 : R) ≤ b + m2 := by
    have := add_le_add hb hm2; rwa [add_zero] at this
  have hsum : a + m1 ≤ b + m2 := le_of_sq_le_sq ham1 hbm2 hsum_sq
  -- a + m1 ≤ b + m2 ↔ a - b ≤ m2 - m1
  apply le_of_sub_nonneg
  have h0 : (0 : R) ≤ (b + m2) - (a + m1) := sub_nonneg_of_le hsum
  have hrew : (b + m2) - (a + m1) = (m2 - m1) - (a - b) := by
    rw [OrderedCommRing.sub_eq_add_neg, OrderedCommRing.sub_eq_add_neg,
        OrderedCommRing.sub_eq_add_neg, OrderedCommRing.sub_eq_add_neg,
        neg_add, neg_add, neg_neg]
    ac_rfl
  rwa [hrew] at h0

/-- Proposition 5.1(i) sharpness characterization:
    equality `a - b = m₂ - m₁` holds if and only if `r = 0 ∨ m₁ = m₂`. -/
theorem vertical_slice_eq_iff (r m1 m2 a b : R)
    (hm1 : (0 : R) ≤ m1) (hm : m1 ≤ m2) (ha : (0 : R) ≤ a) (hb : (0 : R) ≤ b)
    (ha2 : a * a = r * r + m2 * m2) (hb2 : b * b = r * r + m1 * m1) :
    a - b = m2 - m1 ↔ r = 0 ∨ m1 = m2 := by
  have hm2 : (0 : R) ≤ m2 := OrderedCommRing.le_trans hm1 hm
  constructor
  · intro heq
    have hsum : a + m1 = b + m2 := by
      have h0 : (m2 - m1) - (a - b) = 0 := by rw [heq, sub_self]
      have hrew : (m2 - m1) - (a - b) = (b + m2) - (a + m1) := by
        rw [OrderedCommRing.sub_eq_add_neg, OrderedCommRing.sub_eq_add_neg,
            OrderedCommRing.sub_eq_add_neg, OrderedCommRing.sub_eq_add_neg,
            neg_add, neg_add, neg_neg]
        ac_rfl
      rw [hrew] at h0
      exact (sub_eq_zero_iff.mp h0).symm
    have hsum_sq_zero : (b + m2) * (b + m2) - (a + m1) * (a + m1) = 0 := by
      rw [hsum, sub_self]
    rw [slice_sum_sq_identity r m1 m2 a b ha2 hb2] at hsum_sq_zero
    have h2eq : a * m1 + a * m1 = b * m2 + b * m2 := (sub_eq_zero_iff.mp hsum_sq_zero).symm
    have hcross : a * m1 = b * m2 := eq_of_two_mul_eq_two_mul h2eq
    have hcross_sq_zero : (b * m2) * (b * m2) - (a * m1) * (a * m1) = 0 := by
      rw [hcross, sub_self]
    rw [slice_cross_sq_identity r m1 m2 a b ha2 hb2] at hcross_sq_zero
    rcases mul_eq_zero hcross_sq_zero with hr2 | hprod
    · rcases mul_eq_zero hr2 with hr | hr
      · exact Or.inl hr
      · exact Or.inl hr
    · rcases mul_eq_zero hprod with hdiff | hadd
      · exact Or.inr (sub_eq_zero_iff.mp hdiff).symm
      · have hm2_le0 : m2 ≤ 0 := by
          have := OrderedCommRing.add_le_add_left hm1 m2
          rw [add_zero, hadd] at this
          exact this
        have hm2_0 : m2 = 0 := OrderedCommRing.le_antisymm hm2_le0 hm2
        have hm1_0 : m1 = 0 := OrderedCommRing.le_antisymm (by rwa [hm2_0] at hm) hm1
        exact Or.inr (by rw [hm1_0, hm2_0])
  · rintro (hr | hm_eq)
    · rw [hr, mul_zero, zero_add] at ha2 hb2
      have hae : a = m2 := eq_of_sq_eq_sq ha hm2 ha2
      have hbe : b = m1 := eq_of_sq_eq_sq hb hm1 hb2
      rw [hae, hbe]
    · subst hm_eq
      have hab : a = b := eq_of_sq_eq_sq ha hb (by rw [ha2, hb2])
      rw [hab, sub_self, sub_self]

end MassBand
end QrfNorm1
