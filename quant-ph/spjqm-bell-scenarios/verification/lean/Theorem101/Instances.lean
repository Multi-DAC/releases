/-
  Theorem101.Instances — Non-vacuousness witnesses instantiating `AddCommGrp`,
  `OrderedCommRing`, and `RealInnerProd` on Lean's built-in `Int`.
  This guarantees that our parameterized typeclasses are consistent and inhabited.
-/
import Theorem101.Basic

namespace Theorem101.Instances
open Theorem101

instance instAddCommGrpInt : AddCommGrp Int where
  zero := 0
  add := (· + ·)
  neg := (- ·)
  add_assoc := Int.add_assoc
  add_comm := Int.add_comm
  zero_add := Int.zero_add
  neg_add := Int.add_left_neg

instance instOrderedCommRingInt : OrderedCommRing Int where
  zero := 0
  one := 1
  add := (· + ·)
  neg := (- ·)
  mul := (· * ·)
  le := (· ≤ ·)
  add_assoc := Int.add_assoc
  add_comm := Int.add_comm
  zero_add := Int.zero_add
  neg_add := Int.add_left_neg
  mul_assoc := Int.mul_assoc
  mul_comm := Int.mul_comm
  one_mul := Int.one_mul
  left_distrib := Int.mul_add
  le_refl := Int.le_refl
  le_trans := fun _ _ _ => Int.le_trans
  le_antisymm := fun _ _ => Int.le_antisymm
  le_total := Int.le_total
  add_le_add_left := fun a b c h => by omega
  mul_nonneg := fun a b ha hb => Int.mul_nonneg ha hb
  mul_pos := fun a b ha hb => by
    have ha' : (0 : Int) < a := by omega
    have hb' : (0 : Int) < b := by omega
    have hp : (0 : Int) < a * b := Int.mul_pos ha' hb'
    exact ⟨by omega, by omega⟩
  zero_lt_one := by omega

instance instRealInnerProdInt : RealInnerProd Int Int where
  smul := (· * ·)
  inner := (· * ·)
  inner_symm := Int.mul_comm
  inner_add_left := Int.add_mul
  inner_smul_left := Int.mul_assoc
  inner_self_nonneg := fun u => by
    rcases Int.le_total 0 u with h | h
    · exact Int.mul_nonneg h h
    · have hn : 0 ≤ -u := by omega
      have h2 : 0 ≤ (-u) * (-u) := Int.mul_nonneg hn hn
      rw [Int.neg_mul, Int.mul_neg, Int.neg_neg] at h2
      exact h2
  inner_self_eq_zero := fun u h => by
    rcases Int.mul_eq_zero.mp h with h1 | h1 <;> exact h1

/-- Explicit sanity witness: `OrderedCommRing Int` and `RealInnerProd Int Int` are inhabited. -/
theorem typeclasses_nonvacuous :
    Nonempty (OrderedCommRing Int) ∧ Nonempty (RealInnerProd Int Int) :=
  ⟨⟨instOrderedCommRingInt⟩, ⟨instRealInnerProdInt⟩⟩

end Theorem101.Instances
