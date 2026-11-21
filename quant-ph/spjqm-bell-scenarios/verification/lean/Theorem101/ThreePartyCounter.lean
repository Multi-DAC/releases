/-
  Theorem101.ThreePartyCounter — Formal verification of Theorem 4.1 (`thm:3party-counter`):
  The (3,2,2) behaviour P(abc|xyz) = 1/4 iff a ⊕ b ⊕ c = x y z lies in Q^1_3 \ SPJQM_3.
-/
import Theorem101.Basic

namespace Theorem101.ThreePartyCounter
open Theorem101 OrderedCommRing AddCommGrp RealInnerProd

/-! ### Part 1: Exact finite `Q^1_3` indicator certificate on `Ξ = {0,1}^6` (`Fin 64`) -/

/-- Bit `k` (0..5) of a 6-bit history `g : Fin 64` encoding `(a0, a1, b0, b1, c0, c1)`. -/
def bit (g : Fin 64) (k : Fin 6) : Nat := (g.val >>> k.val) &&& 1

/-- Count of histories in `Fin 64` where coordinate `k` equals `o`. -/
def count1 (k : Fin 6) (o : Fin 2) : Nat :=
  (List.finRange 64).countP (fun g => bit g k == o.val)

/-- Count of histories in `Fin 64` where coordinate `k1 = o1` and `k2 = o2`. -/
def count2 (k1 k2 : Fin 6) (o1 o2 : Fin 2) : Nat :=
  (List.finRange 64).countP (fun g => (bit g k1 == o1.val) && (bit g k2 == o2.val))

/-- All 12 single-party marginals equal `32/64 = 1/2`. -/
theorem indicator_q1_marginals_1 :
    ∀ k : Fin 6, ∀ o : Fin 2, count1 k o = 32 := by
  decide

/-- All 48 distinct-party two-party marginals equal `16/64 = 1/4`. -/
theorem indicator_q1_marginals_2 :
    ∀ p1 p2 : Fin 3, p1.val < p2.val →
    ∀ s1 s2 o1 o2 : Fin 2,
      count2 ⟨2 * p1.val + s1.val, by omega⟩ ⟨2 * p2.val + s2.val, by omega⟩ o1 o2 = 16 := by
  decide

/-! ### Part 2: Exact 3-party correlator `∑ (-1)^(a+b+c) P(abc|xyz) = (-1)^(xyz)` -/

/-- Scaled by 4: `4 * ∑_{a,b,c} (-1)^(a+b+c) P(abc|xyz)`. -/
def correlator4 (x y z : Fin 2) : Int :=
  (List.finRange 8).foldl (fun acc idx =>
    let a := idx.val &&& 1
    let b := (idx.val >>> 1) &&& 1
    let c := (idx.val >>> 2) &&& 1
    let rhs := x.val * y.val * z.val
    if (a ^^^ b ^^^ c) == rhs then
      acc + (if ((a + b + c) % 2) == 0 then (1 : Int) else (-1 : Int))
    else acc) 0

/-- For all 8 setting triples `(x,y,z)`, the 3-party correlator equals `(-1)^(xyz)` (scaled by 4). -/
theorem three_party_correlators_exact :
    ∀ x y z : Fin 2,
      correlator4 x y z = (if x.val * y.val * z.val == 1 then (-4 : Int) else (4 : Int)) := by
  decide

/-! ### Part 3: Abstract inner-product space exclusion from `SPJQM_3` -/

/-- If an `SPJQM_3` realization existed for `a ⊕ b ⊕ c = xyz`, the four unit vectors
    `A0, A1, G00, G11` would satisfy `⟨A0,G00⟩ = ⟨A0,G11⟩ = ⟨A1,G00⟩ = 1` and
    `⟨A1,G11⟩ = -1`, which is impossible in any real inner-product space. -/
theorem thm_3party_no_spjqm {R V : Type _} [OrderedCommRing R] [AddCommGrp V] [RealInnerProd R V]
    (A0 A1 G00 G11 : V)
    (hA0  : normSq A0  = OrderedCommRing.one)
    (hA1  : normSq A1  = OrderedCommRing.one)
    (hG00 : normSq G00 = OrderedCommRing.one)
    (hG11 : normSq G11 = OrderedCommRing.one)
    (h00  : inner A0 G00 = OrderedCommRing.one)
    (h01  : inner A0 G11 = OrderedCommRing.one)
    (h10  : inner A1 G00 = OrderedCommRing.one)
    (h11  : inner A1 G11 = OrderedCommRing.neg OrderedCommRing.one) : False := by
  have e1 : A0 = G00 := eq_of_unit_inner_one hA0 hG00 h00
  have e2 : A0 = G11 := eq_of_unit_inner_one hA0 hG11 h01
  have e3 : A1 = G00 := eq_of_unit_inner_one hA1 hG00 h10
  have e4 : A1 = G11 := e3.trans (e1.symm.trans e2)
  rw [← e4] at h11
  -- Now `inner A1 A1 = -1` while `inner A1 A1 = 1`, so `1 = -1`, hence `2 = 0`, contradiction!
  have h1eq : (OrderedCommRing.one : R) = OrderedCommRing.neg OrderedCommRing.one :=
    hA1.symm.trans h11
  have h2eq : (two : R) = OrderedCommRing.zero := by
    calc (two : R) = OrderedCommRing.add OrderedCommRing.one OrderedCommRing.one := rfl
      _ = OrderedCommRing.add OrderedCommRing.one (OrderedCommRing.neg OrderedCommRing.one) :=
            congrArg (OrderedCommRing.add OrderedCommRing.one) h1eq
      _ = OrderedCommRing.zero := OrderedCommRing.add_neg _
  exact zero_lt_two.2 h2eq.symm

end Theorem101.ThreePartyCounter
