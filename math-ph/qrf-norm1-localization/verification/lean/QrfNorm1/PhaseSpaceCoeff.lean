import QrfNorm1.Basic

/-!
# QrfNorm1.PhaseSpaceCoeff

Formal verification of:
1. Kaiser's `2+1`D `n`-body phase-space coefficient recursion and closed form
   (Lemma 5.2(ii), matching `arXiv:2510.08704v2` Eqs. (9)–(10)).
2. The exact Källén-discriminant monotonicity identity in `μ²` used inside the
   inductive step of Lemma 5.2(iii).
-/

namespace QrfNorm1

/-- Standard factorial on `Nat`. -/
def fact : Nat → Nat
  | 0 => 1
  | k + 1 => (k + 1) * fact k

/-- Rational denominator of Kaiser's `2+1`D coefficient `α_n` (modulo `π^{n-2}`),
    defined recursively by `D_2 = 4` and `D_{n+1} = 4 (n - 1) D_n` for `n ≥ 2`. -/
def kaiserDenom : Nat → Nat
  | 0 => 0
  | 1 => 0
  | 2 => 4
  | k + 3 => 4 * (k + 1) * kaiserDenom (k + 2)

/-- Exact closed form `D_n = 4 * 4^(n - 2) * (n - 2)!` for all `n ≥ 2`. -/
theorem kaiserDenom_closed_form (n : Nat) (hn : 2 ≤ n) :
    kaiserDenom n = 4 * 4 ^ (n - 2) * fact (n - 2) := by
  obtain ⟨k, rfl⟩ : ∃ k, n = k + 2 := ⟨n - 2, by omega⟩
  induction k with
  | zero => rfl
  | succ k ih =>
    have hstep : kaiserDenom (k + 3) = 4 * (k + 1) * kaiserDenom (k + 2) := rfl
    rw [hstep, ih (by omega)]
    have hsub1 : k + 3 - 2 = k + 1 := by omega
    have hsub2 : k + 2 - 2 = k := by omega
    rw [hsub1, hsub2, fact, Nat.pow_succ]
    -- 4 * (k + 1) * (4 * 4^k * fact k) = 4 * (4^k * 4) * ((k + 1) * fact k)
    calc 4 * (k + 1) * (4 * 4 ^ k * fact k)
        = (4 * (4 ^ k * 4)) * ((k + 1) * fact k) := by
            ac_rfl

/-- Concrete values matching `Φ_2 = 1/(4M)`, `Φ_3 = (M-3m)/(16π M)`,
    `Φ_4 = (M-4m)²/(128π² M)`, `Φ_5 = (M-5m)³/(1536π³ M)`. -/
theorem kaiserDenom_values :
    kaiserDenom 2 = 4 ∧ kaiserDenom 3 = 16 ∧ kaiserDenom 4 = 128 ∧ kaiserDenom 5 = 1536 := by
  decide

/-- Källén kinematic triangle function `λ(a, b, c) = a² + b² + c² - 2ab - 2ac - 2bc`. -/
def kallen {R : Type} [OrderedCommRing R] (a b c : R) : R :=
  a * a + b * b + c * c - (a * b + a * b) - (a * c + a * c) - (b * c + b * c)

/-- Exact algebraic difference identity for Lemma 5.2(iii):
    `λ(1, t, u₁) - λ(1, t, u₂) = (u₂ - u₁) * (2(1 + t) - (u₁ + u₂))`. -/
theorem kallen_diff_identity {R : Type} [OrderedCommRing R] (t u1 u2 : R) :
    kallen (1 : R) t u1 - kallen (1 : R) t u2 =
      (u2 - u1) * ((1 + t) + (1 + t) - (u1 + u2)) := by
  unfold kallen
  have h1 : ∀ u : R,
      (1 : R) * 1 + t * t + u * u - (1 * t + 1 * t) - (1 * u + 1 * u) - (t * u + t * u) =
      (1 + t * t - (t + t)) + (u * u - ((1 + t) + (1 + t)) * u) := by
    intro u
    simp only [one_mul, right_distrib, OrderedCommRing.sub_eq_add_neg, neg_add]
    ac_rfl
  rw [h1 u1, h1 u2, add_sub_add_left]
  simp only [OrderedCommRing.sub_eq_add_neg, OrderedCommRing.left_distrib, right_distrib,
             mul_neg, neg_mul, neg_add, neg_neg, one_mul, mul_one]
  have hm1 : u2 * t = t * u2 := OrderedCommRing.mul_comm _ _
  have hm2 : u1 * t = t * u1 := OrderedCommRing.mul_comm _ _
  have hm3 : u2 * u1 = u1 * u2 := OrderedCommRing.mul_comm _ _
  rw [hm1, hm2, hm3]
  have hadd0 :
      u1 * u1 + (-(u1) + -(t * u1) + (-(u1) + -(t * u1))) +
        (-(u2 * u2) + (u2 + t * u2 + (u2 + t * u2))) =
      u1 * u1 + (-(u1) + -(t * u1) + (-(u1) + -(t * u1))) +
        (-(u2 * u2) + (u2 + t * u2 + (u2 + t * u2))) + (u1 * u2 + -(u1 * u2)) := by
    rw [add_neg_cancel, add_zero]
  rw [hadd0]
  ac_rfl

/-- Lemma 5.2(iii) kinematic monotonicity: on the physical domain
    `0 ≤ s`, `0 ≤ μ₁ ≤ μ₂`, `μ₂ + s ≤ 1`, decreasing `μ` increases `λ(1, s², μ²)`. -/
theorem kallen_monotone_in_mu_sq {R : Type} [OrderedCommRing R]
    (s mu1 mu2 : R) (hs0 : (0 : R) ≤ s) (hmu1 : (0 : R) ≤ mu1)
    (hle : mu1 ≤ mu2) (hkin : mu2 + s ≤ 1) :
    kallen (1 : R) (s * s) (mu2 * mu2) ≤ kallen (1 : R) (s * s) (mu1 * mu1) := by
  apply le_of_sub_nonneg
  rw [kallen_diff_identity (s * s) (mu1 * mu1) (mu2 * mu2)]
  have hmu2 : (0 : R) ≤ mu2 := OrderedCommRing.le_trans hmu1 hle
  have hmu2_le1 : mu2 ≤ 1 := by
    have h1 : mu2 + 0 ≤ mu2 + s := OrderedCommRing.add_le_add_left hs0 mu2
    rw [add_zero] at h1
    exact OrderedCommRing.le_trans h1 hkin
  have hmu1_le1 : mu1 ≤ 1 := OrderedCommRing.le_trans hle hmu2_le1
  have hu_le : mu1 * mu1 ≤ mu2 * mu2 := mul_le_mul_nonneg hmu1 hmu1 hle hle
  have hu2_le1 : mu2 * mu2 ≤ 1 := by
    have := mul_le_mul_nonneg hmu2 hmu2 hmu2_le1 hmu2_le1
    rwa [mul_one] at this
  have hu1_le1 : mu1 * mu1 ≤ 1 := OrderedCommRing.le_trans hu_le hu2_le1
  have hsum_le2 : mu1 * mu1 + mu2 * mu2 ≤ (1 : R) + 1 := add_le_add hu1_le1 hu2_le1
  have h2_le_S : (1 : R) + 1 ≤ (1 + s * s) + (1 + s * s) := by
    have hs2 : (0 : R) ≤ s * s := sq_nonneg s
    have h1s : (1 : R) ≤ 1 + s * s := by
      have := OrderedCommRing.add_le_add_left hs2 (1 : R)
      rwa [add_zero] at this
    exact add_le_add h1s h1s
  have hfac1 : (0 : R) ≤ mu2 * mu2 - mu1 * mu1 := sub_nonneg_of_le hu_le
  have hfac2 : (0 : R) ≤ (1 + s * s) + (1 + s * s) - (mu1 * mu1 + mu2 * mu2) :=
    sub_nonneg_of_le (OrderedCommRing.le_trans hsum_le2 h2_le_S)
  exact mul_nonneg hfac1 hfac2

end QrfNorm1
