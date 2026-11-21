/-!
# QrfNorm1.Basic

Minimal zero-dependency ordered commutative ring interface (with a concrete
non-vacuousness witness on `Int`) used to formalize the algebraic inequalities
and finite-dimensional eigenvalue identities of Paper 3 (`qrf_norm1`).
-/

namespace QrfNorm1

/-- A linearly ordered commutative ring (sufficient for all polynomial and
    square-root-by-squaring inequalities in Paper 3, covering `Int`, `Rat`, `Real`). -/
class OrderedCommRing (R : Type) extends Add R, Sub R, Mul R, Neg R, LE R, LT R where
  zero : R
  one : R
  add_comm : ∀ a b : R, a + b = b + a
  add_assoc : ∀ a b c : R, (a + b) + c = a + (b + c)
  zero_add : ∀ a : R, zero + a = a
  add_neg_cancel : ∀ a : R, a + (-a) = zero
  sub_eq_add_neg : ∀ a b : R, a - b = a + (-b)
  mul_comm : ∀ a b : R, a * b = b * a
  mul_assoc : ∀ a b c : R, (a * b) * c = a * (b * c)
  one_mul : ∀ a : R, one * a = a
  left_distrib : ∀ a b c : R, a * (b + c) = a * b + a * c
  le_refl : ∀ a : R, a ≤ a
  le_trans : ∀ {a b c : R}, a ≤ b → b ≤ c → a ≤ c
  le_antisymm : ∀ {a b : R}, a ≤ b → b ≤ a → a = b
  le_total : ∀ a b : R, a ≤ b ∨ b ≤ a
  lt_iff_le_not_le : ∀ a b : R, a < b ↔ a ≤ b ∧ ¬(b ≤ a)
  add_le_add_left : ∀ {a b : R}, a ≤ b → ∀ c : R, c + a ≤ c + b
  mul_nonneg : ∀ {a b : R}, zero ≤ a → zero ≤ b → zero ≤ a * b
  zero_lt_one : zero < one
  mul_eq_zero : ∀ {a b : R}, a * b = zero → a = zero ∨ b = zero

instance (R : Type) [OrderedCommRing R] : OfNat R 0 := ⟨OrderedCommRing.zero⟩
instance (R : Type) [OrderedCommRing R] : OfNat R 1 := ⟨OrderedCommRing.one⟩

instance (R : Type) [OrderedCommRing R] : Std.Associative (α := R) (· + ·) :=
  ⟨OrderedCommRing.add_assoc⟩

instance (R : Type) [OrderedCommRing R] : Std.Commutative (α := R) (· + ·) :=
  ⟨OrderedCommRing.add_comm⟩

instance (R : Type) [OrderedCommRing R] : Std.Associative (α := R) (· * ·) :=
  ⟨OrderedCommRing.mul_assoc⟩

instance (R : Type) [OrderedCommRing R] : Std.Commutative (α := R) (· * ·) :=
  ⟨OrderedCommRing.mul_comm⟩

instance : OrderedCommRing Int where
  zero := 0
  one := 1
  add_comm := Int.add_comm
  add_assoc := Int.add_assoc
  zero_add := Int.zero_add
  add_neg_cancel := Int.add_right_neg
  sub_eq_add_neg := fun a b => Int.sub_eq_add_neg
  mul_comm := Int.mul_comm
  mul_assoc := Int.mul_assoc
  one_mul := Int.one_mul
  left_distrib := Int.mul_add
  le_refl := Int.le_refl
  le_trans := Int.le_trans
  le_antisymm := Int.le_antisymm
  le_total := Int.le_total
  lt_iff_le_not_le := fun a b => by omega
  add_le_add_left := fun h c => Int.add_le_add_left h c
  mul_nonneg := fun ha hb => Int.mul_nonneg ha hb
  zero_lt_one := by omega
  mul_eq_zero := fun h => Int.mul_eq_zero.mp h

section RingLemmas
variable {R : Type} [OrderedCommRing R]

theorem zero_add (a : R) : (0 : R) + a = a := OrderedCommRing.zero_add a
theorem add_neg_cancel (a : R) : a + (-a) = (0 : R) := OrderedCommRing.add_neg_cancel a
theorem one_mul (a : R) : (1 : R) * a = a := OrderedCommRing.one_mul a
theorem mul_nonneg {a b : R} (ha : (0 : R) ≤ a) (hb : (0 : R) ≤ b) : (0 : R) ≤ a * b :=
  OrderedCommRing.mul_nonneg ha hb
theorem zero_lt_one : (0 : R) < (1 : R) := OrderedCommRing.zero_lt_one
theorem mul_eq_zero {a b : R} (h : a * b = (0 : R)) : a = (0 : R) ∨ b = (0 : R) :=
  OrderedCommRing.mul_eq_zero h

theorem add_zero (a : R) : a + 0 = a := by
  rw [OrderedCommRing.add_comm, zero_add]

theorem neg_add_cancel (a : R) : (-a) + a = 0 := by
  rw [OrderedCommRing.add_comm, add_neg_cancel]

theorem mul_one (a : R) : a * 1 = a := by
  rw [OrderedCommRing.mul_comm, one_mul]

theorem right_distrib (a b c : R) : (a + b) * c = a * c + b * c := by
  rw [OrderedCommRing.mul_comm, OrderedCommRing.left_distrib,
      OrderedCommRing.mul_comm c a, OrderedCommRing.mul_comm c b]

theorem add_left_cancel {a b c : R} (h : c + a = c + b) : a = b := by
  have h1 : (-c) + (c + a) = (-c) + (c + b) := congrArg (fun x => (-c) + x) h
  rw [← OrderedCommRing.add_assoc, ← OrderedCommRing.add_assoc,
      neg_add_cancel, zero_add, zero_add] at h1
  exact h1

theorem add_right_cancel {a b c : R} (h : a + c = b + c) : a = b := by
  rw [OrderedCommRing.add_comm a c, OrderedCommRing.add_comm b c] at h
  exact add_left_cancel h

theorem mul_zero (a : R) : a * 0 = 0 := by
  have h : a * 0 + a * 0 = a * 0 + 0 := by
    rw [← OrderedCommRing.left_distrib, add_zero, add_zero]
  exact add_left_cancel h

theorem zero_mul (a : R) : 0 * a = 0 := by
  rw [OrderedCommRing.mul_comm, mul_zero]

theorem neg_eq_of_add_eq_zero {a b : R} (h : a + b = 0) : -a = b := by
  have h1 : a + (-a) = a + b := by rw [add_neg_cancel, h]
  exact add_left_cancel h1

theorem neg_zero : -(0 : R) = 0 :=
  neg_eq_of_add_eq_zero (add_zero 0)

theorem neg_neg (a : R) : -(-a) = a :=
  neg_eq_of_add_eq_zero (neg_add_cancel a)

theorem neg_add (a b : R) : -(a + b) = (-a) + (-b) := by
  apply neg_eq_of_add_eq_zero
  calc (a + b) + ((-a) + (-b))
      = a + (b + ((-b) + (-a))) := by
          rw [OrderedCommRing.add_assoc, OrderedCommRing.add_comm (-a) (-b)]
    _ = a + ((b + (-b)) + (-a)) := by rw [← OrderedCommRing.add_assoc b (-b) (-a)]
    _ = a + (0 + (-a)) := by rw [add_neg_cancel]
    _ = a + (-a) := by rw [zero_add]
    _ = 0 := add_neg_cancel a

theorem mul_neg (a b : R) : a * (-b) = -(a * b) := by
  symm
  apply neg_eq_of_add_eq_zero
  rw [← OrderedCommRing.left_distrib, add_neg_cancel, mul_zero]

theorem neg_mul (a b : R) : (-a) * b = -(a * b) := by
  rw [OrderedCommRing.mul_comm, mul_neg, OrderedCommRing.mul_comm]

theorem neg_mul_neg (a b : R) : (-a) * (-b) = a * b := by
  rw [neg_mul, mul_neg, neg_neg]

theorem sub_self (a : R) : a - a = 0 := by
  rw [OrderedCommRing.sub_eq_add_neg, add_neg_cancel]

theorem sub_zero (a : R) : a - 0 = a := by
  rw [OrderedCommRing.sub_eq_add_neg, neg_zero, add_zero]

theorem left_distrib_sub (a b c : R) : a * (b - c) = a * b - a * c := by
  rw [OrderedCommRing.sub_eq_add_neg, OrderedCommRing.left_distrib, mul_neg,
      ← OrderedCommRing.sub_eq_add_neg]

theorem right_distrib_sub (a b c : R) : (a - b) * c = a * c - b * c := by
  rw [OrderedCommRing.sub_eq_add_neg, right_distrib, neg_mul,
      ← OrderedCommRing.sub_eq_add_neg]

theorem add_le_add_right {a b : R} (h : a ≤ b) (c : R) : a + c ≤ b + c := by
  rw [OrderedCommRing.add_comm a c, OrderedCommRing.add_comm b c]
  exact OrderedCommRing.add_le_add_left h c

theorem add_le_add {a b c d : R} (h1 : a ≤ b) (h2 : c ≤ d) : a + c ≤ b + d :=
  OrderedCommRing.le_trans (add_le_add_right h1 c) (OrderedCommRing.add_le_add_left h2 b)

theorem le_of_add_le_add_left {a b c : R} (h : c + a ≤ c + b) : a ≤ b := by
  have h1 := OrderedCommRing.add_le_add_left h (-c)
  rw [← OrderedCommRing.add_assoc, ← OrderedCommRing.add_assoc,
      neg_add_cancel, zero_add, zero_add] at h1
  exact h1

theorem le_of_add_le_add_right {a b c : R} (h : a + c ≤ b + c) : a ≤ b := by
  rw [OrderedCommRing.add_comm a c, OrderedCommRing.add_comm b c] at h
  exact le_of_add_le_add_left h

theorem sub_nonneg_of_le {a b : R} (h : a ≤ b) : (0 : R) ≤ b - a := by
  rw [OrderedCommRing.sub_eq_add_neg]
  have h1 := add_le_add_right h (-a)
  rw [add_neg_cancel] at h1
  exact h1

theorem le_of_sub_nonneg {a b : R} (h : (0 : R) ≤ b - a) : a ≤ b := by
  rw [OrderedCommRing.sub_eq_add_neg] at h
  have h1 := add_le_add_right h a
  rw [zero_add, OrderedCommRing.add_assoc, neg_add_cancel, add_zero] at h1
  exact h1

theorem sub_eq_zero_iff {a b : R} : a - b = 0 ↔ a = b := by
  constructor
  · intro h
    rw [OrderedCommRing.sub_eq_add_neg] at h
    have h1 := congrArg (fun x => x + b) h
    rw [zero_add, OrderedCommRing.add_assoc, neg_add_cancel, add_zero] at h1
    exact h1
  · intro h; rw [h, sub_self]

theorem sq_nonneg (a : R) : (0 : R) ≤ a * a := by
  rcases OrderedCommRing.le_total 0 a with h | h
  · exact mul_nonneg h h
  · have hn : (0 : R) ≤ -a := by
      have := sub_nonneg_of_le h
      rw [OrderedCommRing.sub_eq_add_neg, zero_add] at this
      exact this
    have hmul := mul_nonneg hn hn
    rw [neg_mul_neg] at hmul
    exact hmul

theorem le_of_sq_le_sq {a b : R} (ha : (0 : R) ≤ a) (hb : (0 : R) ≤ b)
    (h : a * a ≤ b * b) : a ≤ b := by
  rcases OrderedCommRing.le_total a b with hab | hba
  · exact hab
  · have hdiff : (0 : R) ≤ (a - b) * (a + b) := by
      exact mul_nonneg (sub_nonneg_of_le hba) (add_le_add ha hb |> fun x => by
        rw [add_zero] at x; exact x)
    have hid : (a - b) * (a + b) = a * a - b * b := by
      rw [right_distrib_sub, OrderedCommRing.left_distrib, OrderedCommRing.left_distrib,
          OrderedCommRing.mul_comm b a]
      rw [OrderedCommRing.sub_eq_add_neg, neg_add, OrderedCommRing.add_assoc,
          ← OrderedCommRing.add_assoc (a * b) (-(a * b)) (-(b * b)),
          add_neg_cancel, zero_add,
          ← OrderedCommRing.sub_eq_add_neg]
    rw [hid] at hdiff
    have hge : b * b ≤ a * a := le_of_sub_nonneg hdiff
    have heq : a * a = b * b := OrderedCommRing.le_antisymm h hge
    have hprod : (a - b) * (a + b) = 0 := by rw [hid, heq, sub_self]
    rcases mul_eq_zero hprod with h1 | h2
    · rw [sub_eq_zero_iff.mp h1]; exact OrderedCommRing.le_refl b
    · have ha0 : a ≤ 0 := by
        have h' : a + 0 ≤ a + b := OrderedCommRing.add_le_add_left hb a
        rw [add_zero, h2] at h'
        exact h'
      have hae : a = 0 := OrderedCommRing.le_antisymm ha0 ha
      rw [hae]; exact hb

theorem eq_of_sq_eq_sq {a b : R} (ha : (0 : R) ≤ a) (hb : (0 : R) ≤ b)
    (h : a * a = b * b) : a = b := by
  apply OrderedCommRing.le_antisymm
  · exact le_of_sq_le_sq ha hb (by rw [h]; exact OrderedCommRing.le_refl _)
  · exact le_of_sq_le_sq hb ha (by rw [h]; exact OrderedCommRing.le_refl _)

theorem mul_le_mul_nonneg {a b c d : R}
    (ha : (0 : R) ≤ a) (hc : (0 : R) ≤ c) (hab : a ≤ b) (hcd : c ≤ d) :
    a * c ≤ b * d := by
  have hb : (0 : R) ≤ b := OrderedCommRing.le_trans ha hab
  have h1 : (0 : R) ≤ (b - a) * c := mul_nonneg (sub_nonneg_of_le hab) hc
  rw [right_distrib_sub] at h1
  have h1' : a * c ≤ b * c := le_of_sub_nonneg h1
  have h2 : (0 : R) ≤ b * (d - c) := mul_nonneg hb (sub_nonneg_of_le hcd)
  rw [left_distrib_sub] at h2
  have h2' : b * c ≤ b * d := le_of_sub_nonneg h2
  exact OrderedCommRing.le_trans h1' h2'

theorem le_of_two_mul_le_two_mul {a b : R} (h : a + a ≤ b + b) : a ≤ b := by
  rcases OrderedCommRing.le_total a b with hab | hba
  · exact hab
  · have hge : b + b ≤ a + a := add_le_add hba hba
    have heq : a + a = b + b := OrderedCommRing.le_antisymm h hge
    have hdiff : ((1 : R) + 1) * (a - b) = 0 := by
      rw [right_distrib, one_mul, OrderedCommRing.sub_eq_add_neg]
      calc (a + (-b)) + (a + (-b))
          = (a + a) + ((-b) + (-b)) := by
              rw [OrderedCommRing.add_assoc, ← OrderedCommRing.add_assoc (-b) a (-b),
                  OrderedCommRing.add_comm (-b) a, OrderedCommRing.add_assoc a (-b) (-b),
                  ← OrderedCommRing.add_assoc]
        _ = (b + b) + (-(b + b)) := by rw [heq, ← neg_add]
        _ = 0 := add_neg_cancel (b + b)
    rcases mul_eq_zero hdiff with h20 | hab0
    · exfalso
      have h01 : (0 : R) ≤ 1 := ((OrderedCommRing.lt_iff_le_not_le 0 1).mp zero_lt_one).1
      have h1_le_2 : (1 : R) ≤ 1 + 1 := by
        have := OrderedCommRing.add_le_add_left h01 1
        rwa [add_zero] at this
      rw [h20] at h1_le_2
      exact ((OrderedCommRing.lt_iff_le_not_le 0 1).mp zero_lt_one).2 h1_le_2
    · rw [sub_eq_zero_iff.mp hab0]
      exact OrderedCommRing.le_refl b

theorem eq_of_two_mul_eq_two_mul {a b : R} (h : a + a = b + b) : a = b :=
  OrderedCommRing.le_antisymm
    (le_of_two_mul_le_two_mul (by rw [h]; exact OrderedCommRing.le_refl _))
    (le_of_two_mul_le_two_mul (by rw [h]; exact OrderedCommRing.le_refl _))

theorem add_sub_add_left (c a b : R) : (c + a) - (c + b) = a - b := by
  simp only [OrderedCommRing.sub_eq_add_neg, neg_add]
  calc (c + a) + (-c + -b)
      = (c + -c) + (a + -b) := by ac_rfl
    _ = 0 + (a + -b) := by rw [add_neg_cancel]
    _ = a + -b := zero_add _

end RingLemmas
end QrfNorm1
