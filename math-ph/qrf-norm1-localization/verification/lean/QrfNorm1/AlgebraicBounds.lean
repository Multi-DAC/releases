import QrfNorm1.Basic

/-!
# QrfNorm1.AlgebraicBounds

Formal verification of the elementary algebraic and finite-dimensional matrix
identities used in Paper 3:
- Theorem 4.1(i) (`thm:op-tail`): 2D Cauchy–Schwarz split inequality
- Theorem 3.1(i) (`thm:open-borel`): final Tonelli floor step `(1 - η/2)(β - η/2) ≥ β - η`
- Theorem 4.2(i) (`thm:lsk`): operator-contraction complement factor `1 - (1 - δ)² = δ(2 - δ)`
- Proposition 4.3(c) (`prop:examples`): exact `2×2` block characteristic polynomial & top eigenvalue
- Corollary 3.4 (`cor:phase`): `N = 2` Gram-block eigenvalues and necessary condition `|⟨ξ_n, ξ_{n+1}⟩| = 1`.
-/

namespace QrfNorm1

section AlgebraicBounds
variable {R : Type} [OrderedCommRing R]

/-- Theorem 4.1(i) 2D Cauchy–Schwarz inequality:
    `(A x + B y)² ≤ (A² + B²)(x² + y²)`. -/
theorem cauchy_schwarz_2d (A B x y : R) :
    (A * x + B * y) * (A * x + B * y) ≤ (A * A + B * B) * (x * x + y * y) := by
  apply le_of_sub_nonneg
  have hid : (A * A + B * B) * (x * x + y * y) - (A * x + B * y) * (A * x + B * y) =
      (A * y - B * x) * (A * y - B * x) := by
    simp only [OrderedCommRing.sub_eq_add_neg, OrderedCommRing.left_distrib, right_distrib,
               mul_neg, neg_mul, neg_add, neg_neg]
    have h1 : A * A * (x * x) = A * x * (A * x) := by ac_rfl
    have h2 : B * B * (y * y) = B * y * (B * y) := by ac_rfl
    have h3 : A * A * (y * y) = A * y * (A * y) := by ac_rfl
    have h4 : B * B * (x * x) = B * x * (B * x) := by ac_rfl
    have h5 : B * y * (A * x) = A * x * (B * y) := by ac_rfl
    have h6 : A * y * (B * x) = A * x * (B * y) := by ac_rfl
    have h7 : B * x * (A * y) = A * x * (B * y) := by ac_rfl
    rw [h1, h2, h3, h4, h5, h6, h7]
    calc A * x * (A * x) + B * x * (B * x) + (A * y * (A * y) + B * y * (B * y)) +
            (-(A * x * (A * x)) + -(A * x * (B * y)) + (-(A * x * (B * y)) + -(B * y * (B * y))))
        = (A * x * (A * x) + -(A * x * (A * x))) + (B * y * (B * y) + -(B * y * (B * y))) +
            (A * y * (A * y) + -(A * x * (B * y)) + (-(A * x * (B * y)) + B * x * (B * x))) := by ac_rfl
      _ = 0 + 0 + (A * y * (A * y) + -(A * x * (B * y)) + (-(A * x * (B * y)) + B * x * (B * x))) := by
            rw [add_neg_cancel, add_neg_cancel]
      _ = A * y * (A * y) + -(A * x * (B * y)) + (-(A * x * (B * y)) + B * x * (B * x)) := by
            rw [add_zero, zero_add]
  rw [hid]
  exact sq_nonneg (A * y - B * x)

/-- Theorem 3.1(i) final Tonelli floor algebra (scaled by 4 to avoid denominators):
    for `β ≤ 1` and `0 ≤ η`, `(2 - η)(2β - η) ≥ 4(β - η)`, i.e.
    `(1 - η/2)(β - η/2) ≥ β - η`. -/
theorem tonelli_floor_step (beta eta : R) (hb1 : beta ≤ 1) (he0 : (0 : R) ≤ eta) :
    (beta - eta) + (beta - eta) + ((beta - eta) + (beta - eta)) ≤
      ((1 + 1) - eta) * ((beta + beta) - eta) := by
  apply le_of_sub_nonneg
  have hid :
      ((1 + 1) - eta) * ((beta + beta) - eta) -
        ((beta - eta) + (beta - eta) + ((beta - eta) + (beta - eta))) =
      eta * (((1 - beta) + (1 - beta)) + eta) := by
    simp only [OrderedCommRing.sub_eq_add_neg, OrderedCommRing.left_distrib, right_distrib,
               mul_neg, neg_mul, neg_add, neg_neg, one_mul, mul_one]
    calc beta + beta + -(eta * beta) + (beta + beta + -(eta * beta)) +
            (-(eta) + -(eta) + eta * eta) +
            (-(beta) + eta + (-(beta) + eta) + (-(beta) + eta + (-(beta) + eta)))
        = ((beta + -(beta)) + (beta + -(beta)) + ((beta + -(beta)) + (beta + -(beta)))) +
            ((eta + -(eta)) + (eta + -(eta))) +
            (eta + -(eta * beta) + (eta + -(eta * beta)) + eta * eta) := by ac_rfl
      _ = ((0 + 0) + (0 + 0)) + (0 + 0) +
            (eta + -(eta * beta) + (eta + -(eta * beta)) + eta * eta) := by
            rw [add_neg_cancel, add_neg_cancel]
      _ = eta + -(eta * beta) + (eta + -(eta * beta)) + eta * eta := by
            simp only [add_zero, zero_add]
  rw [hid]
  have h1b : (0 : R) ≤ 1 - beta := sub_nonneg_of_le hb1
  have h2b : (0 : R) ≤ (1 - beta) + (1 - beta) := by
    have := add_le_add h1b h1b; rwa [add_zero] at this
  have hpos : (0 : R) ≤ ((1 - beta) + (1 - beta)) + eta := by
    have := add_le_add h2b he0; rwa [add_zero] at this
  exact mul_nonneg he0 hpos

/-- Theorem 4.2(i) Logvinenko–Sereda–Kovrijkine contraction complement factor:
    `1 - (1 - δ)² = δ (2 - δ)`. -/
theorem lsk_factor_identity (delta : R) :
    (1 : R) - (1 - delta) * (1 - delta) = delta * ((1 + 1) - delta) := by
  simp only [OrderedCommRing.sub_eq_add_neg, OrderedCommRing.left_distrib, right_distrib,
             mul_neg, neg_mul, neg_add, neg_neg, one_mul, mul_one]
  calc (1 : R) + (-(1) + delta + (delta + -(delta * delta)))
      = (1 + -(1)) + (delta + delta + -(delta * delta)) := by ac_rfl
    _ = 0 + (delta + delta + -(delta * delta)) := by rw [add_neg_cancel]
    _ = delta + delta + -(delta * delta) := zero_add _

/-- Proposition 4.3(c) exact `2×2` block characteristic polynomial:
    if `s² = a₂ d` (where `d = |β₀|² / M`), then
    `(λ - a₂)(λ - d) - s² = λ (λ - (a₂ + d))`,
    so the exact eigenvalues are `0` and `a₂ + d`. -/
theorem ex5c_char_poly (a2 d s lam : R) (hs : s * s = a2 * d) :
    (lam - a2) * (lam - d) - s * s = lam * (lam - (a2 + d)) := by
  rw [hs]
  simp only [OrderedCommRing.sub_eq_add_neg, OrderedCommRing.left_distrib, right_distrib,
             mul_neg, neg_mul, neg_add, neg_neg]
  rw [OrderedCommRing.mul_comm a2 lam]
  calc lam * lam + -(lam * a2) + (-(lam * d) + a2 * d) + -(a2 * d)
      = (lam * lam + (-(lam * a2) + -(lam * d))) + (a2 * d + -(a2 * d)) := by ac_rfl
    _ = (lam * lam + (-(lam * a2) + -(lam * d))) + 0 := by rw [add_neg_cancel]
    _ = lam * lam + (-(lam * a2) + -(lam * d)) := add_zero _

/-- Corollary 3.4 (`N = 2`) Gram-block characteristic polynomial:
    `(λ - 1)² - c² = (λ - (1 + c))(λ - (1 - c))`, so the eigenvalues of
    `[[1, c], [c, 1]]` are `1 + c` and `1 - c`. -/
theorem phase_N2_char_poly (c lam : R) :
    (lam - 1) * (lam - 1) - c * c = (lam - (1 + c)) * (lam - (1 - c)) := by
  have hdiff : ∀ x y : R, x * x - y * y = (x - y) * (x + y) := by
    intro x y
    rw [right_distrib_sub, OrderedCommRing.left_distrib, OrderedCommRing.left_distrib,
        OrderedCommRing.mul_comm y x, OrderedCommRing.sub_eq_add_neg,
        OrderedCommRing.sub_eq_add_neg, neg_add]
    calc x * x + -(y * y)
        = x * x + -(y * y) + 0 := by rw [add_zero]
      _ = x * x + -(y * y) + (x * y + -(x * y)) := by rw [add_neg_cancel]
      _ = (x * x + x * y) + (-(x * y) + -(y * y)) := by ac_rfl
  rw [hdiff (lam - 1) c]
  have h1 : (lam - 1) - c = lam - (1 + c) := by
    rw [OrderedCommRing.sub_eq_add_neg, OrderedCommRing.sub_eq_add_neg,
        OrderedCommRing.sub_eq_add_neg, neg_add, OrderedCommRing.add_assoc]
  have h2 : (lam - 1) + c = lam - (1 - c) := by
    rw [OrderedCommRing.sub_eq_add_neg, OrderedCommRing.sub_eq_add_neg,
        OrderedCommRing.sub_eq_add_neg, neg_add, neg_neg, OrderedCommRing.add_assoc]
  rw [h1, h2]

/-- Corollary 3.4 (`N = 2`) necessary condition: the normalized `N = 2` Gram block
    attains top eigenvalue `2` (`(1 + c)/2 = 1`) iff `c = 1`. -/
theorem phase_N2_necessary (c : R) :
    (1 : R) + c = 1 + 1 ↔ c = 1 := by
  constructor
  · exact add_left_cancel
  · intro h; rw [h]

end AlgebraicBounds
end QrfNorm1
