/-
  Theorem101.Basic — Zero-dependency algebraic & inner-product structures for Paper 1.
  All hypotheses are carried explicitly as typeclass fields; ZERO custom `axiom`
  declarations are introduced anywhere in this project.
-/

namespace Theorem101

/-- An abstract commutative additive group (for vector-space & operator identities). -/
class AddCommGrp (V : Type _) where
  zero : V
  add  : V → V → V
  neg  : V → V
  add_assoc : ∀ a b c : V, add (add a b) c = add a (add b c)
  add_comm  : ∀ a b : V, add a b = add b a
  zero_add  : ∀ a : V, add zero a = a
  neg_add   : ∀ a : V, add (neg a) a = zero

namespace AddCommGrp
variable {V : Type _} [AddCommGrp V]

def sub (a b : V) : V := add a (neg b)

theorem add_zero (a : V) : add a zero = a := by
  rw [add_comm, zero_add]

theorem add_neg (a : V) : add a (neg a) = zero := by
  rw [add_comm, neg_add]

theorem neg_zero : neg (zero : V) = zero := by
  have h := add_neg (zero : V)
  rw [zero_add] at h
  exact h

theorem sub_self (a : V) : sub a a = zero := add_neg a

theorem eq_of_sub_eq_zero {a b : V} (h : sub a b = zero) : a = b := by
  have h1 : add (sub a b) b = add zero b := congrArg (fun x => add x b) h
  unfold sub at h1
  rw [add_assoc, neg_add, add_zero, zero_add] at h1
  exact h1

theorem sub_eq_zero_of_eq {a b : V} (h : a = b) : sub a b = zero := by
  rw [h, sub_self]

theorem add_sub_cancel (a b : V) : sub (add a b) b = a := by
  unfold sub
  rw [add_assoc, add_neg, add_zero]

theorem sub_add_cancel (a b : V) : add (sub a b) b = a := by
  unfold sub
  rw [add_assoc, neg_add, add_zero]

theorem neg_neg (a : V) : neg (neg a) = a := by
  have h1 : add (neg (neg a)) (neg a) = zero := neg_add (neg a)
  have h2 := congrArg (fun x => add x a) h1
  rw [add_assoc, neg_add, add_zero, zero_add] at h2
  exact h2

theorem neg_add_rev (a b : V) : neg (add a b) = add (neg a) (neg b) := by
  have h1 : add (add a b) (add (neg a) (neg b)) = zero := by
    rw [add_assoc a b _, ← add_assoc b (neg a) _, add_comm b (neg a),
        add_assoc (neg a) b _, add_neg, add_zero, add_neg]
  have h2 := congrArg (fun x => add (neg (add a b)) x) h1
  rw [← add_assoc, neg_add, zero_add, add_zero] at h2
  exact h2.symm

end AddCommGrp

/-- A strictly ordered commutative ring (instantiated by `ℝ`, `ℚ`, `ℤ`). -/
class OrderedCommRing (R : Type _) where
  zero : R
  one  : R
  add  : R → R → R
  neg  : R → R
  mul  : R → R → R
  le   : R → R → Prop
  add_assoc : ∀ a b c : R, add (add a b) c = add a (add b c)
  add_comm  : ∀ a b : R, add a b = add b a
  zero_add  : ∀ a : R, add zero a = a
  neg_add   : ∀ a : R, add (neg a) a = zero
  mul_assoc : ∀ a b c : R, mul (mul a b) c = mul a (mul b c)
  mul_comm  : ∀ a b : R, mul a b = mul b a
  one_mul   : ∀ a : R, mul one a = a
  left_distrib : ∀ a b c : R, mul a (add b c) = add (mul a b) (mul a c)
  le_refl   : ∀ a : R, le a a
  le_trans  : ∀ a b c : R, le a b → le b c → le a c
  le_antisymm : ∀ a b : R, le a b → le b a → a = b
  le_total  : ∀ a b : R, le a b ∨ le b a
  add_le_add_left : ∀ a b c : R, le a b → le (add c a) (add c b)
  mul_nonneg : ∀ a b : R, le zero a → le zero b → le zero (mul a b)
  mul_pos    : ∀ a b : R, (le zero a ∧ zero ≠ a) → (le zero b ∧ zero ≠ b) →
                 (le zero (mul a b) ∧ zero ≠ mul a b)
  zero_lt_one : ¬ (le one zero)

namespace OrderedCommRing
variable {R : Type _} [OrderedCommRing R]

def sub (a b : R) : R := add a (neg b)
def lt (a b : R) : Prop := le a b ∧ a ≠ b
def two : R := add one one
def four : R := add two two
def eight : R := add four four
def sq (a : R) : R := mul a a

theorem add_zero (a : R) : add a zero = a := by rw [add_comm, zero_add]
theorem add_neg (a : R) : add a (neg a) = zero := by rw [add_comm, neg_add]
theorem mul_one (a : R) : mul a one = a := by rw [mul_comm, one_mul]
theorem right_distrib (a b c : R) : mul (add a b) c = add (mul a c) (mul b c) := by
  rw [mul_comm, left_distrib, mul_comm c a, mul_comm c b]

theorem eq_of_sub_eq_zero {a b : R} (h : sub a b = zero) : a = b := by
  have h1 := congrArg (fun x => add x b) h
  unfold sub at h1
  rw [add_assoc, neg_add, add_zero, zero_add] at h1
  exact h1

theorem sub_self (a : R) : sub a a = zero := add_neg a

theorem mul_zero (a : R) : mul a zero = zero := by
  have h : add (mul a zero) (mul a zero) = mul a zero := by
    rw [← left_distrib, zero_add]
  have h2 := congrArg (fun x => add x (neg (mul a zero))) h
  rw [add_assoc, add_neg, add_zero] at h2
  exact h2

theorem zero_mul (a : R) : mul zero a = zero := by rw [mul_comm, mul_zero]

theorem neg_neg (a : R) : neg (neg a) = a := by
  have h1 : add (neg (neg a)) (neg a) = zero := neg_add (neg a)
  have h2 := congrArg (fun x => add x a) h1
  rw [add_assoc, neg_add, add_zero, zero_add] at h2
  exact h2

theorem neg_zero : neg (zero : R) = zero := by
  have h := add_neg (zero : R)
  rw [zero_add] at h
  exact h

theorem mul_neg (a b : R) : mul a (neg b) = neg (mul a b) := by
  have h1 : add (mul a (neg b)) (mul a b) = zero := by
    rw [← left_distrib, neg_add, mul_zero]
  have h2 := congrArg (fun x => add x (neg (mul a b))) h1
  rw [add_assoc, add_neg, add_zero, zero_add] at h2
  exact h2

theorem neg_mul (a b : R) : mul (neg a) b = neg (mul a b) := by
  rw [mul_comm, mul_neg, mul_comm]

theorem neg_add_rev (a b : R) : neg (add a b) = add (neg a) (neg b) := by
  have h1 : add (add a b) (add (neg a) (neg b)) = zero := by
    rw [add_assoc a b _, ← add_assoc b (neg a) _, add_comm b (neg a),
        add_assoc (neg a) b _, add_neg, add_zero, add_neg]
  have h2 := congrArg (fun x => add (neg (add a b)) x) h1
  rw [← add_assoc, neg_add, zero_add, add_zero] at h2
  exact h2.symm

theorem add_le_add_right {a b : R} (h : le a b) (c : R) : le (add a c) (add b c) := by
  rw [add_comm a c, add_comm b c]
  exact add_le_add_left a b c h

theorem add_le_add {a b c d : R} (h1 : le a b) (h2 : le c d) : le (add a c) (add b d) :=
  le_trans _ _ _ (add_le_add_right h1 c) (add_le_add_left c d b h2)

theorem sub_nonneg_of_le {a b : R} (h : le a b) : le zero (sub b a) := by
  have h1 := add_le_add_right h (neg a)
  rw [add_neg] at h1
  exact h1

theorem le_of_sub_nonneg {a b : R} (h : le zero (sub b a)) : le a b := by
  have h1 := add_le_add_right h a
  unfold sub at h1
  rw [zero_add, add_assoc, neg_add, add_zero] at h1
  exact h1

theorem neg_nonneg_of_nonpos {a : R} (h : le a zero) : le zero (neg a) := by
  have h1 := add_le_add_right h (neg a)
  rw [add_neg, zero_add] at h1
  exact h1

theorem nonpos_of_neg_nonneg {a : R} (h : le zero (neg a)) : le a zero := by
  have h1 := add_le_add_right h a
  rw [zero_add, neg_add] at h1
  exact h1

theorem sq_nonneg (a : R) : le zero (sq a) := by
  unfold sq
  rcases le_total zero a with h | h
  · exact mul_nonneg a a h h
  · have hn := neg_nonneg_of_nonpos h
    have h2 := mul_nonneg (neg a) (neg a) hn hn
    rw [mul_neg, neg_mul, neg_neg] at h2
    exact h2

theorem zero_le_one : le (zero : R) one := by
  rcases le_total (zero : R) one with h | h
  · exact h
  · exact False.elim (zero_lt_one h)

theorem zero_lt_one' : lt (zero : R) one :=
  ⟨zero_le_one, fun h => zero_lt_one (by rw [← h]; exact le_refl _)⟩

theorem lt_of_le_of_lt {a b c : R} (h1 : le a b) (h2 : lt b c) : lt a c := by
  refine ⟨le_trans _ _ _ h1 h2.1, fun heq => h2.2 ?_⟩
  rw [heq] at h1
  exact le_antisymm _ _ h2.1 h1

theorem lt_of_lt_of_le {a b c : R} (h1 : lt a b) (h2 : le b c) : lt a c := by
  refine ⟨le_trans _ _ _ h1.1 h2, fun heq => h1.2 ?_⟩
  rw [← heq] at h2
  exact le_antisymm _ _ h1.1 h2

theorem add_pos_of_pos_of_nonneg {a b : R} (ha : lt zero a) (hb : le zero b) :
    lt zero (add a b) := by
  have h1 : le a (add a b) := by
    have h := add_le_add_left zero b a hb
    rw [add_zero] at h
    exact h
  exact lt_of_lt_of_le ha h1

theorem zero_lt_two : lt (zero : R) two :=
  add_pos_of_pos_of_nonneg zero_lt_one' zero_le_one

theorem zero_lt_four : lt (zero : R) four :=
  add_pos_of_pos_of_nonneg zero_lt_two zero_lt_two.1

theorem zero_lt_eight : lt (zero : R) eight :=
  add_pos_of_pos_of_nonneg zero_lt_four zero_lt_four.1

theorem mul_le_mul_of_nonneg_left {a b c : R} (h : le a b) (hc : le zero c) :
    le (mul c a) (mul c b) := by
  apply le_of_sub_nonneg
  have hsub : sub (mul c b) (mul c a) = mul c (sub b a) := by
    unfold sub; rw [left_distrib, mul_neg]
  rw [hsub]
  exact mul_nonneg c (sub b a) hc (sub_nonneg_of_le h)

theorem mul_le_mul_of_nonneg_right {a b c : R} (h : le a b) (hc : le zero c) :
    le (mul a c) (mul b c) := by
  rw [mul_comm a c, mul_comm b c]
  exact mul_le_mul_of_nonneg_left h hc

theorem mul_le_mul {a b c d : R} (hab : le a b) (hcd : le c d)
    (hc : le zero c) (hb : le zero b) : le (mul a c) (mul b d) :=
  le_trans _ _ _ (mul_le_mul_of_nonneg_right hab hc) (mul_le_mul_of_nonneg_left hcd hb)

theorem nonneg_of_mul_nonneg_right {a b : R} (h : le zero (mul a b)) (hb : lt zero b) :
    le zero a := by
  rcases le_total zero a with ha | ha
  · exact ha
  · by_cases haz : a = zero
    · rw [haz]; exact le_refl _
    · have hna : lt zero (neg a) := by
        refine ⟨neg_nonneg_of_nonpos ha, fun heq => haz ?_⟩
        have h2 := congrArg neg heq
        rw [neg_zero, neg_neg] at h2
        exact h2.symm
      have hpos := mul_pos (neg a) b hna hb
      rw [neg_mul] at hpos
      have hle0 := nonpos_of_neg_nonneg hpos.1
      have hzero : mul a b = zero := le_antisymm _ _ hle0 h
      rw [hzero, neg_zero] at hpos
      exact False.elim (hpos.2 rfl)

theorem nonneg_of_mul_nonneg_left {a b : R} (h : le zero (mul a b)) (ha : lt zero a) :
    le zero b := by
  rw [mul_comm] at h
  exact nonneg_of_mul_nonneg_right h ha

theorem le_of_mul_le_mul_left {a b c : R} (h : le (mul c a) (mul c b)) (hc : lt zero c) :
    le a b := by
  apply le_of_sub_nonneg
  have h1 := sub_nonneg_of_le h
  have h2 : sub (mul c b) (mul c a) = mul c (sub b a) := by
    unfold sub; rw [left_distrib, mul_neg]
  rw [h2] at h1
  exact nonneg_of_mul_nonneg_left h1 hc

theorem le_of_mul_le_mul_right {a b c : R} (h : le (mul a c) (mul b c)) (hc : lt zero c) :
    le a b := by
  rw [mul_comm a c, mul_comm b c] at h
  exact le_of_mul_le_mul_left h hc

theorem eq_zero_of_sq_eq_zero {a : R} (h : sq a = zero) : a = zero := by
  by_cases ha : a = zero
  · exact ha
  · rcases le_total zero a with hle | hle
    · have hp := mul_pos a a ⟨hle, fun e => ha e.symm⟩ ⟨hle, fun e => ha e.symm⟩
      exact False.elim (hp.2 h.symm)
    · have hna : lt zero (neg a) := by
        refine ⟨neg_nonneg_of_nonpos hle, fun e => ha ?_⟩
        have h2 := congrArg neg e; rw [neg_zero, neg_neg] at h2; exact h2.symm
      have hp := mul_pos (neg a) (neg a) hna hna
      rw [mul_neg, neg_mul, neg_neg] at hp
      exact False.elim (hp.2 h.symm)

theorem le_of_sq_le_sq {a b : R} (_ha : le zero a) (hb : le zero b) (h : le (sq a) (sq b)) :
    le a b := by
  rcases le_total a b with h1 | h1
  · exact h1
  · by_cases hab : a = b
    · rw [hab]; exact le_refl b
    · have hdiff : lt zero (sub a b) := by
        refine ⟨sub_nonneg_of_le h1, fun e => hab ?_⟩
        exact eq_of_sub_eq_zero e.symm
      have hsum : lt zero (add a b) := by
        have ha_pos : lt zero a := lt_of_le_of_lt hb ⟨h1, fun e => hab e.symm⟩
        exact add_pos_of_pos_of_nonneg ha_pos hb
      have hprod := mul_pos (sub a b) (add a b) hdiff hsum
      have hexp : mul (sub a b) (add a b) = sub (sq a) (sq b) := by
        unfold sub sq
        rw [right_distrib, left_distrib, left_distrib, neg_mul, neg_mul,
            mul_comm b a, add_assoc (mul a a) (mul a b),
            ← add_assoc (mul a b) (neg (mul a b)), add_neg, zero_add]
      rw [hexp] at hprod
      have hle := sub_nonneg_of_le h
      have hsub_neg : sub (sq a) (sq b) = neg (sub (sq b) (sq a)) := by
        unfold sub; rw [neg_add_rev, neg_neg, add_comm]
      rw [hsub_neg] at hprod
      have hle0 := nonpos_of_neg_nonneg hprod.1
      have hz := le_antisymm _ _ hle0 hle
      rw [hz, neg_zero] at hprod
      exact False.elim (hprod.2 rfl)

theorem sq_add (a b : R) :
    sq (add a b) = add (add (sq a) (mul two (mul a b))) (sq b) := by
  unfold sq two
  rw [right_distrib, left_distrib, left_distrib, right_distrib, one_mul, mul_comm b a,
      add_assoc (mul a a) (mul a b), ← add_assoc (mul a b) (mul a b),
      ← add_assoc (mul a a)]

theorem sq_sub (a b : R) :
    sq (sub a b) = sub (add (sq a) (sq b)) (mul two (mul a b)) := by
  unfold sub
  rw [sq_add, mul_neg, mul_neg]
  have hnb : sq (neg b) = sq b := by
    unfold sq; rw [mul_neg, neg_mul, neg_neg]
  rw [hnb, add_assoc (sq a), add_comm (neg (mul two (mul a b))) (sq b), ← add_assoc]

/-- Elementary AM–GM: `4 a b ≤ (a + b)²`, with equality iff `a = b`. -/
theorem four_mul_le_sq_add (a b : R) : le (mul four (mul a b)) (sq (add a b)) := by
  apply le_of_sub_nonneg
  have hid : sub (sq (add a b)) (mul four (mul a b)) = sq (sub a b) := by
    rw [sq_add, sq_sub]
    unfold sub four
    rw [right_distrib, neg_add_rev]
    -- (a² + 2ab) + b² + (-2ab + -2ab) = (a² + b²) + -2ab
    rw [add_assoc (add (sq a) (mul two (mul a b))),
        add_comm (sq b) (add (neg (mul two (mul a b))) (neg (mul two (mul a b)))),
        ← add_assoc (add (sq a) (mul two (mul a b))),
        add_assoc (sq a) (mul two (mul a b)),
        ← add_assoc (mul two (mul a b)) (neg (mul two (mul a b))),
        add_neg, zero_add,
        add_assoc (sq a), add_comm (neg (mul two (mul a b))) (sq b), ← add_assoc]
  rw [hid]
  exact sq_nonneg (sub a b)

theorem eq_of_four_mul_eq_sq_add {a b : R} (h : mul four (mul a b) = sq (add a b)) :
    a = b := by
  have hid : sub (sq (add a b)) (mul four (mul a b)) = sq (sub a b) := by
    rw [sq_add, sq_sub]
    unfold sub four
    rw [right_distrib, neg_add_rev]
    rw [add_assoc (add (sq a) (mul two (mul a b))),
        add_comm (sq b) (add (neg (mul two (mul a b))) (neg (mul two (mul a b)))),
        ← add_assoc (add (sq a) (mul two (mul a b))),
        add_assoc (sq a) (mul two (mul a b)),
        ← add_assoc (mul two (mul a b)) (neg (mul two (mul a b))),
        add_neg, zero_add,
        add_assoc (sq a), add_comm (neg (mul two (mul a b))) (sq b), ← add_assoc]
  rw [h, sub_self] at hid
  exact eq_of_sub_eq_zero (eq_zero_of_sq_eq_zero hid.symm)

theorem two_mul_eq_add (a : R) : mul two a = add a a := by
  unfold two; rw [right_distrib, one_mul]

theorem four_mul_eq (a : R) : mul four a = mul two (mul two a) := by
  unfold four; rw [right_distrib, ← two_mul_eq_add]

end OrderedCommRing

/-- A real inner-product space over a strictly ordered commutative ring `R`. -/
class RealInnerProd (R : outParam (Type _)) (V : Type _) [OrderedCommRing R] [AddCommGrp V] where
  smul  : R → V → V
  inner : V → V → R
  inner_symm : ∀ u v : V, inner u v = inner v u
  inner_add_left : ∀ u v w : V, inner (AddCommGrp.add u v) w =
                                 OrderedCommRing.add (inner u w) (inner v w)
  inner_smul_left : ∀ (c : R) (u v : V), inner (smul c u) v = OrderedCommRing.mul c (inner u v)
  inner_self_nonneg : ∀ u : V, OrderedCommRing.le OrderedCommRing.zero (inner u u)
  inner_self_eq_zero : ∀ u : V, inner u u = OrderedCommRing.zero → u = AddCommGrp.zero

namespace RealInnerProd
open OrderedCommRing AddCommGrp
variable {R V : Type _} [OrderedCommRing R] [AddCommGrp V] [RealInnerProd R V]

def normSq (u : V) : R := inner (R := R) u u

theorem inner_add_right (u v w : V) :
    inner (R := R) u (AddCommGrp.add v w) = OrderedCommRing.add (inner u v) (inner u w) := by
  rw [inner_symm, inner_add_left, inner_symm v u, inner_symm w u]

theorem inner_smul_right (c : R) (u v : V) :
    inner (R := R) u (smul c v) = OrderedCommRing.mul c (inner u v) := by
  rw [inner_symm, inner_smul_left, inner_symm]

theorem inner_zero_left (v : V) : inner (R := R) (AddCommGrp.zero : V) v = OrderedCommRing.zero := by
  have h : inner (R := R) (AddCommGrp.add AddCommGrp.zero AddCommGrp.zero) v =
           OrderedCommRing.add (inner AddCommGrp.zero v) (inner AddCommGrp.zero v) :=
    inner_add_left _ _ _
  rw [AddCommGrp.zero_add] at h
  have h2 := congrArg (fun x => OrderedCommRing.add x (OrderedCommRing.neg (inner (R := R) AddCommGrp.zero v))) h
  rw [OrderedCommRing.add_neg, OrderedCommRing.add_assoc, OrderedCommRing.add_neg,
      OrderedCommRing.add_zero] at h2
  exact h2.symm

theorem inner_zero_right (u : V) : inner (R := R) u (AddCommGrp.zero : V) = OrderedCommRing.zero := by
  rw [inner_symm, inner_zero_left]

theorem inner_neg_left (u v : V) : inner (R := R) (AddCommGrp.neg u) v = OrderedCommRing.neg (inner u v) := by
  have h : OrderedCommRing.add (inner (R := R) (AddCommGrp.neg u) v) (inner u v) = OrderedCommRing.zero := by
    rw [← inner_add_left, AddCommGrp.neg_add, inner_zero_left]
  have h2 := congrArg (fun x => OrderedCommRing.add x (OrderedCommRing.neg (inner (R := R) u v))) h
  rw [OrderedCommRing.add_assoc, OrderedCommRing.add_neg, OrderedCommRing.add_zero,
      OrderedCommRing.zero_add] at h2
  exact h2

theorem inner_neg_right (u v : V) : inner (R := R) u (AddCommGrp.neg v) = OrderedCommRing.neg (inner u v) := by
  rw [inner_symm, inner_neg_left, inner_symm]

theorem inner_sub_left (u v w : V) :
    inner (R := R) (AddCommGrp.sub u v) w = OrderedCommRing.sub (inner u w) (inner v w) := by
  unfold AddCommGrp.sub OrderedCommRing.sub
  rw [inner_add_left, inner_neg_left]

theorem inner_sub_right (u v w : V) :
    inner (R := R) u (AddCommGrp.sub v w) = OrderedCommRing.sub (inner u v) (inner u w) := by
  unfold AddCommGrp.sub OrderedCommRing.sub
  rw [inner_add_right, inner_neg_right]

theorem normSq_add (u v : V) :
    normSq (R := R) (AddCommGrp.add u v) =
    OrderedCommRing.add (OrderedCommRing.add (normSq u) (OrderedCommRing.mul two (inner u v)))
                        (normSq v) := by
  unfold normSq
  rw [inner_add_left, inner_add_right, inner_add_right, inner_symm v u, two_mul_eq_add]
  rw [OrderedCommRing.add_assoc (inner u u), ← OrderedCommRing.add_assoc (inner u v),
      ← OrderedCommRing.add_assoc (inner u u)]

theorem normSq_sub (u v : V) :
    normSq (R := R) (AddCommGrp.sub u v) =
    OrderedCommRing.sub (OrderedCommRing.add (normSq u) (normSq v))
                        (OrderedCommRing.mul two (inner u v)) := by
  unfold AddCommGrp.sub
  rw [normSq_add, inner_neg_right, OrderedCommRing.mul_neg]
  have hnv : normSq (R := R) (AddCommGrp.neg v) = normSq v := by
    unfold normSq; rw [inner_neg_left, inner_neg_right, OrderedCommRing.neg_neg]
  rw [hnv]
  unfold OrderedCommRing.sub
  rw [OrderedCommRing.add_assoc (normSq u),
      OrderedCommRing.add_comm (OrderedCommRing.neg (OrderedCommRing.mul two (inner u v))) (normSq v),
      ← OrderedCommRing.add_assoc]

theorem pythagoras {u v : V} (h : inner (R := R) u v = OrderedCommRing.zero) :
    normSq (R := R) (AddCommGrp.add u v) = OrderedCommRing.add (normSq u) (normSq v) := by
  rw [normSq_add, h, OrderedCommRing.mul_zero, OrderedCommRing.add_zero]

/-- Unit vectors with inner product `1` must be equal (`‖u - v‖² = 1 + 1 - 2 = 0`). -/
theorem eq_of_unit_inner_one {u v : V}
    (hu : normSq (R := R) u = OrderedCommRing.one)
    (hv : normSq (R := R) v = OrderedCommRing.one)
    (huv : inner (R := R) u v = OrderedCommRing.one) : u = v := by
  apply AddCommGrp.eq_of_sub_eq_zero
  apply inner_self_eq_zero (R := R)
  change normSq (R := R) (AddCommGrp.sub u v) = OrderedCommRing.zero
  rw [normSq_sub, hu, hv, huv, OrderedCommRing.mul_one]
  change OrderedCommRing.sub two two = OrderedCommRing.zero
  exact OrderedCommRing.sub_self two

/-- Cauchy–Schwarz over any strictly ordered commutative ring (no division needed!):
    `⟨u, v⟩² ≤ ‖u‖² · ‖v‖²`. -/
theorem cauchy_schwarz (u v : V) :
    OrderedCommRing.le (sq (inner (R := R) u v)) (OrderedCommRing.mul (normSq u) (normSq v)) := by
  by_cases hv : normSq (R := R) v = OrderedCommRing.zero
  · have hv0 : v = AddCommGrp.zero := inner_self_eq_zero (R := R) v hv
    rw [hv, hv0, inner_zero_right, OrderedCommRing.mul_zero]
    unfold sq
    rw [OrderedCommRing.mul_zero]
    exact OrderedCommRing.le_refl _
  · have hv_pos : OrderedCommRing.lt OrderedCommRing.zero (normSq (R := R) v) :=
      ⟨inner_self_nonneg (R := R) v, fun e => hv e.symm⟩
    let c : R := normSq v
    let b : R := inner u v
    let w : V := AddCommGrp.sub (smul c u) (smul b v)
    have hw : OrderedCommRing.le OrderedCommRing.zero (normSq (R := R) w) := inner_self_nonneg w
    have hexp : normSq (R := R) w =
        OrderedCommRing.mul c (OrderedCommRing.sub (OrderedCommRing.mul (normSq u) c) (sq b)) := by
      unfold w normSq
      rw [inner_sub_left, inner_sub_right, inner_sub_right,
          inner_smul_left, inner_smul_right,
          inner_smul_left, inner_smul_right,
          inner_smul_left, inner_smul_right,
          inner_smul_left, inner_smul_right,
          inner_symm v u]
      change OrderedCommRing.sub
        (OrderedCommRing.sub (OrderedCommRing.mul c (OrderedCommRing.mul c (normSq u)))
                             (OrderedCommRing.mul c (OrderedCommRing.mul b b)))
        (OrderedCommRing.sub (OrderedCommRing.mul b (OrderedCommRing.mul c b))
                             (OrderedCommRing.mul b (OrderedCommRing.mul b c))) = _
      have h0 : OrderedCommRing.sub (OrderedCommRing.mul b (OrderedCommRing.mul c b))
                                    (OrderedCommRing.mul b (OrderedCommRing.mul b c)) =
                OrderedCommRing.zero := by
        rw [OrderedCommRing.mul_comm c b]
        exact OrderedCommRing.sub_self _
      rw [h0]
      unfold OrderedCommRing.sub sq
      rw [OrderedCommRing.neg_zero, OrderedCommRing.add_zero,
          OrderedCommRing.left_distrib, OrderedCommRing.mul_neg,
          OrderedCommRing.mul_comm c (normSq u)]
      rfl
    rw [hexp] at hw
    have hsub := nonneg_of_mul_nonneg_left hw hv_pos
    exact le_of_sub_nonneg hsub

/-- Squared triangle inequality: if `su² = ‖u‖²`, `sv² = ‖v‖²` with `su, sv ≥ 0`,
    then `‖u + v‖² ≤ (su + sv)²`. -/
theorem normSq_add_le_sq_add {u v : V} {su sv : R}
    (hsu0 : OrderedCommRing.le OrderedCommRing.zero su)
    (hsv0 : OrderedCommRing.le OrderedCommRing.zero sv)
    (hsu : sq su = normSq (R := R) u) (hsv : sq sv = normSq (R := R) v) :
    OrderedCommRing.le (normSq (R := R) (AddCommGrp.add u v)) (sq (OrderedCommRing.add su sv)) := by
  rw [normSq_add, sq_add, hsu, hsv]
  have hcs : OrderedCommRing.le (inner (R := R) u v) (OrderedCommRing.mul su sv) := by
    rcases OrderedCommRing.le_total OrderedCommRing.zero (inner (R := R) u v) with hpos | hneg
    · apply le_of_sq_le_sq hpos (OrderedCommRing.mul_nonneg su sv hsu0 hsv0)
      have hsq : sq (OrderedCommRing.mul su sv) = OrderedCommRing.mul (normSq u) (normSq v) := by
        rw [← hsu, ← hsv]
        unfold sq
        rw [OrderedCommRing.mul_assoc, ← OrderedCommRing.mul_assoc sv su sv,
            OrderedCommRing.mul_comm sv su, OrderedCommRing.mul_assoc su sv sv,
            ← OrderedCommRing.mul_assoc]
      rw [hsq]
      exact cauchy_schwarz u v
    · exact OrderedCommRing.le_trans _ _ _ hneg (OrderedCommRing.mul_nonneg su sv hsu0 hsv0)
  have h2 : OrderedCommRing.le (OrderedCommRing.mul two (inner (R := R) u v))
                               (OrderedCommRing.mul two (OrderedCommRing.mul su sv)) :=
    mul_le_mul_of_nonneg_left hcs zero_lt_two.1
  exact add_le_add_right (OrderedCommRing.add_le_add_left _ _ (normSq u) h2) (normSq v)

theorem normSq_neg (u : V) : normSq (R := R) (AddCommGrp.neg u) = normSq u := by
  unfold normSq
  rw [inner_neg_left, inner_neg_right, OrderedCommRing.neg_neg]

end RealInnerProd

end Theorem101
