/-!
# QrfNorm1.ScalingLaw

Formal verification of the Diophantine scaling-dimension classification
in Theorem 5.3 (`thm:wick`) and the integrability exponent lower bound
in Lemma 5.2(iii) (`lem:phase-space`) of Paper 3.
-/

namespace QrfNorm1

/-- Exact exponent factorization: `n(d - 2) - d = (n - 1)(d - 2) - 2`. -/
theorem wick_exponent_identity (d n : Int) :
    n * (d - 2) - d = (n - 1) * (d - 2) - 2 := by
  rw [Int.sub_mul, Int.one_mul]
  omega

/-- Auxiliary product classification on positive integers: `m * k = 2`. -/
private theorem mul_eq_two_iff (m k : Nat) (hm : 1 ≤ m) (hk : 1 ≤ k) :
    m * k = 2 ↔ (m = 1 ∧ k = 2) ∨ (m = 2 ∧ k = 1) := by
  constructor
  · intro h
    have hm_le : m ≤ 2 := by
      have : m * 1 ≤ m * k := Nat.mul_le_mul_left m hk
      omega
    have hk_le : k ≤ 2 := by
      have : 1 * k ≤ m * k := Nat.mul_le_mul_right k hm
      omega
    match m, k, hm, hk, hm_le, hk_le, h with
    | 1, 1, _, _, _, _, h => simp at h
    | 1, 2, _, _, _, _, _ => exact Or.inl ⟨rfl, rfl⟩
    | 2, 1, _, _, _, _, _ => exact Or.inr ⟨rfl, rfl⟩
    | 2, 2, _, _, _, _, h => simp at h
  · rintro (⟨rfl, rfl⟩ | ⟨rfl, rfl⟩) <;> rfl

/-- Auxiliary product classification on positive integers: `m * k < 2`. -/
private theorem mul_lt_two_iff (m k : Nat) (hm : 1 ≤ m) (hk : 1 ≤ k) :
    m * k < 2 ↔ m = 1 ∧ k = 1 := by
  constructor
  · intro h
    have hm_le : m ≤ 1 := by
      have : m * 1 ≤ m * k := Nat.mul_le_mul_left m hk
      omega
    have hk_le : k ≤ 1 := by
      have : 1 * k ≤ m * k := Nat.mul_le_mul_right k hm
      omega
    omega
  · rintro ⟨rfl, rfl⟩; decide

/-- Theorem 5.3(i) — Critical scaling dimension `n(d - 2) = d` for `d ≥ 3, n ≥ 2`:
    holds iff `(d, n) ∈ {(4, 2), (3, 3)}` (i.e. `:φ²:` in `3+1`D or `:φ³:` in `2+1`D). -/
theorem wick_critical_iff (d n : Nat) (hd : 3 ≤ d) (hn : 2 ≤ n) :
    n * (d - 2) = d ↔ (d = 4 ∧ n = 2) ∨ (d = 3 ∧ n = 3) := by
  have hrel : n * (d - 2) = (n - 1) * (d - 2) + (d - 2) := by
    have : n = (n - 1) + 1 := by omega
    conv => lhs; rw [this]
    rw [Nat.add_mul, Nat.one_mul]
  have hiff : n * (d - 2) = d ↔ (n - 1) * (d - 2) = 2 := by omega
  rw [hiff, mul_eq_two_iff (n - 1) (d - 2) (by omega) (by omega)]
  omega

/-- Theorem 5.3(ii) — Subcritical scaling dimension `n(d - 2) < d` for `d ≥ 3, n ≥ 2`:
    unique solution `(d, n) = (3, 2)` (i.e. `:φ²:` in `2+1`D). -/
theorem wick_subcritical_iff (d n : Nat) (hd : 3 ≤ d) (hn : 2 ≤ n) :
    n * (d - 2) < d ↔ d = 3 ∧ n = 2 := by
  have hrel : n * (d - 2) = (n - 1) * (d - 2) + (d - 2) := by
    have : n = (n - 1) + 1 := by omega
    conv => lhs; rw [this]
    rw [Nat.add_mul, Nat.one_mul]
  have hiff : n * (d - 2) < d ↔ (n - 1) * (d - 2) < 2 := by omega
  rw [hiff, mul_lt_two_iff (n - 1) (d - 2) (by omega) (by omega)]
  omega

/-- Theorem 5.3(iii) — Supercritical scaling dimension `n(d - 2) > d` for `d ≥ 3, n ≥ 2`:
    every remaining pair `(d, n)`. -/
theorem wick_supercritical_iff (d n : Nat) (hd : 3 ≤ d) (hn : 2 ≤ n) :
    d < n * (d - 2) ↔ ¬((d = 4 ∧ n = 2) ∨ (d = 3 ∧ n = 3) ∨ (d = 3 ∧ n = 2)) := by
  have hcrit := wick_critical_iff d n hd hn
  have hsub := wick_subcritical_iff d n hd hn
  omega

/-- Lemma 5.2(iii) integrability exponent bound: for all `d ≥ 3` and `n ≥ 3`,
    `2 * α_{n,d} = (n - 1)(d - 2) - d ≥ -1` (hence `α_{n,d} ≥ -1/2 > -1`),
    with equality iff `(d, n) = (3, 3)`. -/
theorem alpha_integrability_bound (d n : Nat) (hd : 3 ≤ d) (hn : 3 ≤ n) :
    d ≤ (n - 1) * (d - 2) + 1 := by
  have h1 : 2 * (d - 2) ≤ (n - 1) * (d - 2) := Nat.mul_le_mul_right (d - 2) (by omega)
  omega

theorem alpha_integrability_eq_iff (d n : Nat) (hd : 3 ≤ d) (hn : 3 ≤ n) :
    (n - 1) * (d - 2) + 1 = d ↔ d = 3 ∧ n = 3 := by
  have h1 : 2 * (d - 2) ≤ (n - 1) * (d - 2) := Nat.mul_le_mul_right (d - 2) (by omega)
  constructor
  · intro h
    have hd3 : d = 3 := by omega
    subst hd3
    simp at h
    omega
  · rintro ⟨rfl, rfl⟩; decide

/-- Remark 5.4: in `d = 1 + 1 = 2` dimensions, every `n ≥ 2` has engineering
    exponent `n(d - 2) - d = -2 < 0` (always subcritical). -/
theorem dim_1p1_always_subcritical (n : Int) :
    n * ((2 : Int) - 2) - 2 = -2 := by
  simp

end QrfNorm1
