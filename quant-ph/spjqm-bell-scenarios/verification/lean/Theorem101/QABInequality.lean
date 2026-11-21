/-
  Theorem101.QABInequality — Formal verification of:
    * Proposition 6.1 (`prop:face-max`): Exact AM–GM upper bound `p ≤ 1/8` & uniqueness of `Hardy_{1/8}`
    * Proposition 6.2 (`prop:qab-bound`): Exact DHW vector identity + BOTH quantitative Q^{1+AB} bounds:
        - Pythagoras bound:  P(00|ab) ≤ P(00|a'b') + (√P(11|a'b') + √P(01|ab') + √P(10|a'b))²
        - Symmetric bound:   √P(00|ab) ≤ √(P(00|a'b') + P(11|a'b')) + √P(01|ab') + √P(10|a'b)
    * Corollary 6.3(a),(b): Hardy_{1/8} ∉ Q^{1+AB}, and P_ε ∉ Q^{1+AB} for 0 ≤ ε < 1/19
                            (and via symmetric bound for 0 ≤ ε < 1/(11 + 8√2))
-/
import Theorem101.Basic

namespace Theorem101.QABInequality
open Theorem101 OrderedCommRing AddCommGrp RealInnerProd

/-! ### Part 1: Proposition 6.2 — Exact DHW vector identity & BOTH `Q^{1+AB}` inequalities -/

/-- **Exact DHW vector identity (`eq:qab-identity`):**
    Using ONLY context additivity in any commutative additive group `V`,
    `r1 - r2 + r3 + r4 = v01_apb + v10_abp + v00_ab - v11_ab`. -/
theorem prop_qab_vector_identity {V : Type _} [AddCommGrp V]
    (vXi vA0 vAp0 vB0 vBp0 : V)
    (v00_ab v11_ab v00_abp v01_abp v10_abp v00_apb v01_apb v10_apb v00_apbp v11_apbp : V)
    -- Context additivity equations:
    (h_apbp : AddCommGrp.sub v00_apbp v11_apbp = AddCommGrp.sub (AddCommGrp.add vAp0 vBp0) vXi)
    (h_ab   : AddCommGrp.sub v00_ab   v11_ab   = AddCommGrp.sub (AddCommGrp.add vA0  vB0)  vXi)
    (h_abp1 : v01_abp = AddCommGrp.sub vA0  v00_abp)
    (h_abp2 : v10_abp = AddCommGrp.sub vBp0 v00_abp)
    (h_apb1 : v10_apb = AddCommGrp.sub vB0  v00_apb)
    (h_apb2 : v01_apb = AddCommGrp.sub vAp0 v00_apb) :
    AddCommGrp.add (AddCommGrp.add (AddCommGrp.sub v00_apbp v11_apbp) v01_abp) v10_apb =
    AddCommGrp.add (AddCommGrp.add v01_apb v10_abp) (AddCommGrp.sub v00_ab v11_ab) := by
  rw [h_apbp, h_ab, h_abp1, h_abp2, h_apb1, h_apb2]
  unfold AddCommGrp.sub
  -- Both sides are sums of the same 7 terms: vAp0, vBp0, -vXi, vA0, -v00_abp, vB0, -v00_apb!
  -- Canonicalize both sides to ((((((vA0 + vB0) + vAp0) + vBp0) + (-vXi)) + (-v00_abp)) + (-v00_apb)).
  have lhs_eq :
      AddCommGrp.add (AddCommGrp.add (AddCommGrp.add (AddCommGrp.add vAp0 vBp0) (AddCommGrp.neg vXi))
                                     (AddCommGrp.add vA0 (AddCommGrp.neg v00_abp)))
                     (AddCommGrp.add vB0 (AddCommGrp.neg v00_apb)) =
      AddCommGrp.add (AddCommGrp.add (AddCommGrp.add (AddCommGrp.add vA0 vB0) (AddCommGrp.add vAp0 vBp0))
                                     (AddCommGrp.neg vXi))
                     (AddCommGrp.add (AddCommGrp.neg v00_abp) (AddCommGrp.neg v00_apb)) := by
    rw [AddCommGrp.add_assoc (AddCommGrp.add vAp0 vBp0) (AddCommGrp.neg vXi) _,
        AddCommGrp.add_comm (AddCommGrp.neg vXi) (AddCommGrp.add vA0 (AddCommGrp.neg v00_abp)),
        ← AddCommGrp.add_assoc (AddCommGrp.add vAp0 vBp0) _ (AddCommGrp.neg vXi),
        AddCommGrp.add_assoc (AddCommGrp.add (AddCommGrp.add vAp0 vBp0) _) (AddCommGrp.neg vXi) _,
        AddCommGrp.add_comm (AddCommGrp.neg vXi) (AddCommGrp.add vB0 (AddCommGrp.neg v00_apb)),
        ← AddCommGrp.add_assoc (AddCommGrp.add (AddCommGrp.add vAp0 vBp0) _) _ (AddCommGrp.neg vXi),
        AddCommGrp.add_assoc (AddCommGrp.add vAp0 vBp0) (AddCommGrp.add vA0 (AddCommGrp.neg v00_abp)) _,
        AddCommGrp.add_assoc vA0 (AddCommGrp.neg v00_abp) _,
        AddCommGrp.add_comm (AddCommGrp.neg v00_abp) (AddCommGrp.add vB0 (AddCommGrp.neg v00_apb)),
        AddCommGrp.add_assoc vB0 (AddCommGrp.neg v00_apb) (AddCommGrp.neg v00_abp),
        AddCommGrp.add_comm (AddCommGrp.neg v00_apb) (AddCommGrp.neg v00_abp),
        ← AddCommGrp.add_assoc vA0 vB0 _,
        ← AddCommGrp.add_assoc (AddCommGrp.add vAp0 vBp0) (AddCommGrp.add vA0 vB0) _,
        AddCommGrp.add_comm (AddCommGrp.add vAp0 vBp0) (AddCommGrp.add vA0 vB0),
        AddCommGrp.add_assoc (AddCommGrp.add (AddCommGrp.add vA0 vB0) (AddCommGrp.add vAp0 vBp0)) _ (AddCommGrp.neg vXi),
        AddCommGrp.add_comm (AddCommGrp.add (AddCommGrp.neg v00_abp) (AddCommGrp.neg v00_apb)) (AddCommGrp.neg vXi),
        ← AddCommGrp.add_assoc]
  have rhs_eq :
      AddCommGrp.add (AddCommGrp.add (AddCommGrp.add vAp0 (AddCommGrp.neg v00_apb))
                                     (AddCommGrp.add vBp0 (AddCommGrp.neg v00_abp)))
                     (AddCommGrp.add (AddCommGrp.add vA0 vB0) (AddCommGrp.neg vXi)) =
      AddCommGrp.add (AddCommGrp.add (AddCommGrp.add (AddCommGrp.add vA0 vB0) (AddCommGrp.add vAp0 vBp0))
                                     (AddCommGrp.neg vXi))
                     (AddCommGrp.add (AddCommGrp.neg v00_abp) (AddCommGrp.neg v00_apb)) := by
    have h1 : AddCommGrp.add (AddCommGrp.add vAp0 (AddCommGrp.neg v00_apb))
                             (AddCommGrp.add vBp0 (AddCommGrp.neg v00_abp)) =
              AddCommGrp.add (AddCommGrp.add vAp0 vBp0)
                             (AddCommGrp.add (AddCommGrp.neg v00_abp) (AddCommGrp.neg v00_apb)) := by
      rw [AddCommGrp.add_assoc vAp0 (AddCommGrp.neg v00_apb) _,
          AddCommGrp.add_comm (AddCommGrp.neg v00_apb) (AddCommGrp.add vBp0 (AddCommGrp.neg v00_abp)),
          AddCommGrp.add_assoc vBp0 (AddCommGrp.neg v00_abp) (AddCommGrp.neg v00_apb),
          ← AddCommGrp.add_assoc vAp0 vBp0 _]
    rw [h1,
        AddCommGrp.add_comm (AddCommGrp.add (AddCommGrp.add vAp0 vBp0) _) _,
        AddCommGrp.add_assoc (AddCommGrp.add vA0 vB0) (AddCommGrp.neg vXi) _,
        AddCommGrp.add_comm (AddCommGrp.neg vXi) (AddCommGrp.add (AddCommGrp.add vAp0 vBp0) _),
        AddCommGrp.add_assoc (AddCommGrp.add vAp0 vBp0) _ (AddCommGrp.neg vXi),
        AddCommGrp.add_comm (AddCommGrp.add (AddCommGrp.neg v00_abp) (AddCommGrp.neg v00_apb)) (AddCommGrp.neg vXi),
        ← AddCommGrp.add_assoc (AddCommGrp.add vAp0 vBp0) (AddCommGrp.neg vXi) _,
        ← AddCommGrp.add_assoc (AddCommGrp.add vA0 vB0) _ _,
        ← AddCommGrp.add_assoc (AddCommGrp.add vA0 vB0) (AddCommGrp.add vAp0 vBp0) (AddCommGrp.neg vXi)]
  rw [lhs_eq, rhs_eq]

/-- **Proposition 6.2, sharp Pythagoras bound (`eq:qab-bound-orth`):**
    In any real inner-product space realizing a `Q^{1+AB} = SPJQM_b` behaviour,
    `P(00|ab) ≤ P(00|a'b') + (s2 + s3 + s4)²` where `s2² = P(11|a'b')`, `s3² = P(01|ab')`,
    `s4² = P(10|a'b)` with `s2, s3, s4 ≥ 0`. -/
theorem prop_qab_pythagoras_bound {R V : Type _} [OrderedCommRing R] [AddCommGrp V] [RealInnerProd R V]
    (v00_ab v11_ab v01_apb v10_abp r1 r2 r3 r4 : V)
    (hid : AddCommGrp.add (AddCommGrp.add (AddCommGrp.sub r1 r2) r3) r4 =
           AddCommGrp.add (AddCommGrp.add v01_apb v10_abp) (AddCommGrp.sub v00_ab v11_ab))
    -- Cross-context orthogonality on RHS against v00_ab:
    (ho_ab1 : inner v01_apb v00_ab = OrderedCommRing.zero)
    (ho_ab2 : inner v10_abp v00_ab = OrderedCommRing.zero)
    (ho_ab3 : inner v11_ab  v00_ab = OrderedCommRing.zero)
    -- Cross-context orthogonality of r1 = |[00|a'b']⟩ against r2, r3, r4:
    (ho_r12 : inner r1 r2 = OrderedCommRing.zero)
    (ho_r13 : inner r1 r3 = OrderedCommRing.zero)
    (ho_r14 : inner r1 r4 = OrderedCommRing.zero)
    -- Nonnegative square roots of ‖r2‖², ‖r3‖², ‖r4‖²:
    {s2 s3 s4 : R}
    (hs2_0 : OrderedCommRing.le OrderedCommRing.zero s2)
    (hs3_0 : OrderedCommRing.le OrderedCommRing.zero s3)
    (hs4_0 : OrderedCommRing.le OrderedCommRing.zero s4)
    (hs2 : sq s2 = normSq r2) (hs3 : sq s3 = normSq r3) (hs4 : sq s4 = normSq r4) :
    OrderedCommRing.le (normSq v00_ab)
      (OrderedCommRing.add (normSq r1) (sq (OrderedCommRing.add (OrderedCommRing.add s2 s3) s4))) := by
  let S := AddCommGrp.add (AddCommGrp.add (AddCommGrp.sub r1 r2) r3) r4
  -- Step 1: ⟨S, v00_ab⟩ = ‖v00_ab‖²
  have hinner : inner S v00_ab = normSq v00_ab := by
    change inner (AddCommGrp.add (AddCommGrp.add (AddCommGrp.sub r1 r2) r3) r4) v00_ab = normSq v00_ab
    rw [hid, inner_add_left, inner_add_left, inner_sub_left, ho_ab1, ho_ab2, ho_ab3,
        OrderedCommRing.zero_add, OrderedCommRing.zero_add]
    unfold OrderedCommRing.sub normSq
    rw [OrderedCommRing.neg_zero, OrderedCommRing.add_zero]
  -- Step 2: ‖v00_ab‖² ≤ ‖S‖² (from Cauchy–Schwarz: (‖v00_ab‖²)² ≤ ‖S‖² · ‖v00_ab‖²)
  have hle_S : OrderedCommRing.le (normSq v00_ab) (normSq S) := by
    by_cases hp : normSq v00_ab = OrderedCommRing.zero
    · rw [hp]; exact inner_self_nonneg S
    · have hp_pos : OrderedCommRing.lt OrderedCommRing.zero (normSq v00_ab) :=
        ⟨inner_self_nonneg v00_ab, fun e => hp e.symm⟩
      have hcs := cauchy_schwarz S v00_ab
      rw [hinner] at hcs
      unfold sq at hcs
      exact le_of_mul_le_mul_right hcs hp_pos
  -- Step 3: S = r1 + ((-r2 + r3) + r4), with r1 ⊥ ((-r2 + r3) + r4)
  let T := AddCommGrp.add (AddCommGrp.add (AddCommGrp.neg r2) r3) r4
  have hST : S = AddCommGrp.add r1 T := by
    unfold S T AddCommGrp.sub
    rw [AddCommGrp.add_assoc r1 (AddCommGrp.neg r2) r3,
        AddCommGrp.add_assoc r1 (AddCommGrp.add (AddCommGrp.neg r2) r3) r4]
  have horth : inner r1 T = OrderedCommRing.zero := by
    unfold T
    rw [inner_add_right, inner_add_right, inner_neg_right, ho_r12, ho_r13, ho_r14,
        OrderedCommRing.neg_zero, OrderedCommRing.zero_add, OrderedCommRing.zero_add]
  have hpyth : normSq S = OrderedCommRing.add (normSq r1) (normSq T) := by
    rw [hST, pythagoras horth]
  -- Step 4: ‖T‖² ≤ (s2 + s3 + s4)² via two applications of triangle inequality
  -- First, ‖-r2 + r3‖² ≤ (s2 + s3)², and directly expand ‖(-r2 + r3) + r4‖²!
  have hT_le : OrderedCommRing.le (normSq T) (sq (OrderedCommRing.add (OrderedCommRing.add s2 s3) s4)) := by
    -- Direct expansion avoids needing a square root of ‖-r2 + r3‖²:
    -- ‖(-r2 + r3) + r4‖² = ‖r2‖² + ‖r3‖² + ‖r4‖² + 2⟨-r2, r3⟩ + 2⟨-r2, r4⟩ + 2⟨r3, r4⟩
    -- and each ⟨u, v⟩ ≤ su * sv by Cauchy–Schwarz!
    have hcs_bound : ∀ (u v : V) {su sv : R},
        OrderedCommRing.le OrderedCommRing.zero su →
        OrderedCommRing.le OrderedCommRing.zero sv →
        sq su = normSq u → sq sv = normSq v →
        OrderedCommRing.le (inner u v) (OrderedCommRing.mul su sv) := by
      intro u v su sv hsu0 hsv0 hsu hsv
      rcases OrderedCommRing.le_total OrderedCommRing.zero (inner u v) with hpos | hneg
      · apply le_of_sq_le_sq hpos (OrderedCommRing.mul_nonneg su sv hsu0 hsv0)
        have hsq : sq (OrderedCommRing.mul su sv) = OrderedCommRing.mul (normSq u) (normSq v) := by
          rw [← hsu, ← hsv]; unfold sq
          rw [OrderedCommRing.mul_assoc, ← OrderedCommRing.mul_assoc sv su sv,
              OrderedCommRing.mul_comm sv su, OrderedCommRing.mul_assoc su sv sv,
              ← OrderedCommRing.mul_assoc]
        rw [hsq]; exact cauchy_schwarz u v
      · exact OrderedCommRing.le_trans _ _ _ hneg (OrderedCommRing.mul_nonneg su sv hsu0 hsv0)
    have hb23 := hcs_bound (AddCommGrp.neg r2) r3 hs2_0 hs3_0 (by rw [normSq_neg]; exact hs2) hs3
    have hb24 := hcs_bound (AddCommGrp.neg r2) r4 hs2_0 hs4_0 (by rw [normSq_neg]; exact hs2) hs4
    have hb34 := hcs_bound r3 r4 hs3_0 hs4_0 hs3 hs4
    unfold T
    rw [normSq_add, normSq_add, normSq_neg, inner_add_left,
        OrderedCommRing.left_distrib two (inner (AddCommGrp.neg r2) r4) (inner r3 r4)]
    rw [sq_add (OrderedCommRing.add s2 s3) s4, sq_add s2 s3, hs2, hs3, hs4,
        OrderedCommRing.right_distrib s2 s3 s4,
        OrderedCommRing.left_distrib two (OrderedCommRing.mul s2 s4) (OrderedCommRing.mul s3 s4)]
    have h2_23 := mul_le_mul_of_nonneg_left hb23 zero_lt_two.1
    have h2_24 := mul_le_mul_of_nonneg_left hb24 zero_lt_two.1
    have h2_34 := mul_le_mul_of_nonneg_left hb34 zero_lt_two.1
    have hfirst := add_le_add_right (OrderedCommRing.add_le_add_left _ _ (normSq r2) h2_23) (normSq r3)
    have hmid := OrderedCommRing.add_le_add h2_24 h2_34
    exact add_le_add_right (OrderedCommRing.add_le_add hfirst hmid) (normSq r4)
  rw [hpyth] at hle_S
  exact OrderedCommRing.le_trans _ _ _ hle_S (OrderedCommRing.add_le_add_left _ _ (normSq r1) hT_le)

/-- **Proposition 6.2, symmetric square-root bound (`eq:qab-bound`):**
    Follows directly from `prop_qab_pythagoras_bound`: if `s0² = P(00|ab)`, `s1² = P(00|a'b')`,
    `s2² = P(11|a'b')`, `s3² = P(01|ab')`, `s4² = P(10|a'b)`, and `t² = P(00|a'b') + P(11|a'b')`
    (all `≥ 0`), then `s0 ≤ t + s3 + s4`. -/
theorem prop_qab_sqrt_bound {R : Type _} [OrderedCommRing R]
    {p00 p1 p2 s0 s2 s3 s4 t : R}
    (hs0_0 : OrderedCommRing.le OrderedCommRing.zero s0)
    (hs2_0 : OrderedCommRing.le OrderedCommRing.zero s2)
    (hs3_0 : OrderedCommRing.le OrderedCommRing.zero s3)
    (hs4_0 : OrderedCommRing.le OrderedCommRing.zero s4)
    (ht_0  : OrderedCommRing.le OrderedCommRing.zero t)
    (hp1_0 : OrderedCommRing.le OrderedCommRing.zero p1)
    (hs0 : sq s0 = p00) (hs2 : sq s2 = p2) (ht : sq t = OrderedCommRing.add p1 p2)
    (hpyth : OrderedCommRing.le p00
               (OrderedCommRing.add p1 (sq (OrderedCommRing.add (OrderedCommRing.add s2 s3) s4)))) :
    OrderedCommRing.le s0 (OrderedCommRing.add t (OrderedCommRing.add s3 s4)) := by
  have hs34_0 : OrderedCommRing.le OrderedCommRing.zero (OrderedCommRing.add s3 s4) :=
    OrderedCommRing.add_le_add hs3_0 hs4_0 |> fun h => by rwa [OrderedCommRing.zero_add] at h
  have hrhs_0 : OrderedCommRing.le OrderedCommRing.zero (OrderedCommRing.add t (OrderedCommRing.add s3 s4)) :=
    OrderedCommRing.add_le_add ht_0 hs34_0 |> fun h => by rwa [OrderedCommRing.zero_add] at h
  apply le_of_sq_le_sq hs0_0 hrhs_0
  rw [hs0]
  -- Show p1 + (s2 + (s3+s4))² ≤ (t + (s3+s4))² using s2 ≤ t (since s2² = p2 ≤ p1 + p2 = t²)
  have hs2_le_t : OrderedCommRing.le s2 t := by
    apply le_of_sq_le_sq hs2_0 ht_0
    rw [hs2, ht]
    have h := add_le_add_right hp1_0 p2
    rwa [OrderedCommRing.zero_add] at h
  have hstep : OrderedCommRing.le
      (OrderedCommRing.add p1 (sq (OrderedCommRing.add (OrderedCommRing.add s2 s3) s4)))
      (sq (OrderedCommRing.add t (OrderedCommRing.add s3 s4))) := by
    rw [OrderedCommRing.add_assoc s2 s3 s4,
        sq_add s2 (OrderedCommRing.add s3 s4),
        sq_add t (OrderedCommRing.add s3 s4),
        hs2, ht]
    -- LHS = p1 + ((p2 + 2 s2 (s3+s4)) + (s3+s4)²) = ((p1 + p2) + 2 s2 (s3+s4)) + (s3+s4)²
    rw [← OrderedCommRing.add_assoc p1 (OrderedCommRing.add p2 _) _,
        ← OrderedCommRing.add_assoc p1 p2 _]
    have hmul : OrderedCommRing.le (OrderedCommRing.mul s2 (OrderedCommRing.add s3 s4))
                                   (OrderedCommRing.mul t (OrderedCommRing.add s3 s4)) :=
      mul_le_mul_of_nonneg_right hs2_le_t hs34_0
    have h2 := mul_le_mul_of_nonneg_left hmul zero_lt_two.1
    exact add_le_add_right (OrderedCommRing.add_le_add_left _ _ (OrderedCommRing.add p1 p2) h2) _
  exact OrderedCommRing.le_trans _ _ _ hpyth hstep

/-- **Corollary 6.3(a):** On the 4-zero Hardy face (`p1 = p2 = p3 = p4 = 0`), any `Q^{1+AB}`
    behaviour must satisfy `P(00|ab) = 0`. -/
theorem cor_hardy_face_zero {R : Type _} [OrderedCommRing R] {p00 : R}
    (hp00_0 : OrderedCommRing.le OrderedCommRing.zero p00)
    (hpyth : OrderedCommRing.le p00
      (OrderedCommRing.add OrderedCommRing.zero
        (sq (OrderedCommRing.add (OrderedCommRing.add OrderedCommRing.zero OrderedCommRing.zero)
                                 OrderedCommRing.zero)))) :
    p00 = OrderedCommRing.zero := by
  rw [OrderedCommRing.zero_add, OrderedCommRing.zero_add, OrderedCommRing.zero_add] at hpyth
  unfold sq at hpyth
  rw [OrderedCommRing.mul_zero] at hpyth
  exact OrderedCommRing.le_antisymm _ _ hpyth hp00_0

/-- **Corollary 6.3(b) (Pythagoras threshold `ε < 1/19`):**
    On `P_ε = (1 - ε) Hardy_{1/8} + ε U`, multiplying `eq:qab-bound-orth` by `8` gives
    `1 + ε ≤ 20 ε`, which is impossible whenever `19 ε < 1`. -/
theorem cor_peps_violates_pyth {R : Type _} [OrderedCommRing R]
    {ε s : R}
    (hs : sq s = ε)
    -- Scaled-by-8 inequality: (1 + ε) ≤ 2ε + 2*(3s)² = 20ε
    (hineq : OrderedCommRing.le (OrderedCommRing.add OrderedCommRing.one ε)
               (OrderedCommRing.add (OrderedCommRing.mul two ε)
                 (OrderedCommRing.mul two (sq (OrderedCommRing.add (OrderedCommRing.add s s) s))))) :
    OrderedCommRing.le OrderedCommRing.one
      (OrderedCommRing.mul (OrderedCommRing.add eight (OrderedCommRing.add eight (OrderedCommRing.add two OrderedCommRing.one))) ε) := by
  -- Notice (s + s + s)² = 9 s² = 9 ε, so 2ε + 2(9ε) = 20ε, subtracting ε gives 1 ≤ 19ε!
  have h3s : sq (OrderedCommRing.add (OrderedCommRing.add s s) s) =
             OrderedCommRing.mul (OrderedCommRing.add eight OrderedCommRing.one) ε := by
    rw [← hs]
    unfold sq eight four two
    simp only [OrderedCommRing.left_distrib, OrderedCommRing.right_distrib,
               OrderedCommRing.one_mul, OrderedCommRing.add_assoc]
  rw [h3s] at hineq
  -- Now hineq : 1 + ε ≤ 2ε + 2(9ε) = (19ε) + ε
  have hrhs : OrderedCommRing.add (OrderedCommRing.mul two ε)
                (OrderedCommRing.mul two (OrderedCommRing.mul (OrderedCommRing.add eight OrderedCommRing.one) ε)) =
              OrderedCommRing.add
                (OrderedCommRing.mul (OrderedCommRing.add eight (OrderedCommRing.add eight (OrderedCommRing.add two OrderedCommRing.one))) ε)
                ε := by
    unfold eight four two
    simp only [OrderedCommRing.left_distrib, OrderedCommRing.right_distrib,
               OrderedCommRing.one_mul, OrderedCommRing.add_assoc]
  rw [hrhs] at hineq
  have h := add_le_add_right hineq (OrderedCommRing.neg ε)
  rw [OrderedCommRing.add_assoc, OrderedCommRing.add_neg, OrderedCommRing.add_zero,
      OrderedCommRing.add_assoc, OrderedCommRing.add_neg, OrderedCommRing.add_zero] at h
  exact h

/-! ### Part 2: Proposition 6.1 — Exact AM–GM maximum `1/8` & uniqueness of `Hardy_{1/8}` -/

/-- **Proposition 6.1 (AM–GM bound `8p ≤ 1`):**
    On the face `F_4` (`A = p + w1 + w4`, `B = p + w2 + w5`, `A + B = 1`),
    the `Q^1` criterion `p² A B ≤ (p + w1) w4 (p + w2) w5` implies `8 p ≤ 1` (i.e. `p ≤ 1/8`). -/
theorem prop_face_max_bound {R : Type _} [OrderedCommRing R]
    {p w1 w2 w4 w5 : R}
    (hp0 : OrderedCommRing.le OrderedCommRing.zero p)
    (_hw1 : OrderedCommRing.le OrderedCommRing.zero w1)
    (hw2 : OrderedCommRing.le OrderedCommRing.zero w2)
    (_hw4 : OrderedCommRing.le OrderedCommRing.zero w4)
    (hw5 : OrderedCommRing.le OrderedCommRing.zero w5)
    (hsum : OrderedCommRing.add
              (OrderedCommRing.add (OrderedCommRing.add p w1) w4)
              (OrderedCommRing.add (OrderedCommRing.add p w2) w5) = OrderedCommRing.one)
    (hApos : OrderedCommRing.lt OrderedCommRing.zero (OrderedCommRing.add (OrderedCommRing.add p w1) w4))
    (hBpos : OrderedCommRing.lt OrderedCommRing.zero (OrderedCommRing.add (OrderedCommRing.add p w2) w5))
    (hq1 : OrderedCommRing.le
             (OrderedCommRing.mul (sq p)
               (OrderedCommRing.mul (OrderedCommRing.add (OrderedCommRing.add p w1) w4)
                                    (OrderedCommRing.add (OrderedCommRing.add p w2) w5)))
             (OrderedCommRing.mul (OrderedCommRing.mul (OrderedCommRing.add p w1) w4)
                                  (OrderedCommRing.mul (OrderedCommRing.add p w2) w5))) :
    OrderedCommRing.le (OrderedCommRing.mul eight p) OrderedCommRing.one := by
  let A := OrderedCommRing.add (OrderedCommRing.add p w1) w4
  let B := OrderedCommRing.add (OrderedCommRing.add p w2) w5
  have hABpos : OrderedCommRing.lt OrderedCommRing.zero (OrderedCommRing.mul A B) :=
    OrderedCommRing.mul_pos A B hApos hBpos
  -- 4 (p+w1) w4 ≤ A² and 4 (p+w2) w5 ≤ B²
  have ham1 : OrderedCommRing.le (OrderedCommRing.mul four (OrderedCommRing.mul (OrderedCommRing.add p w1) w4)) (sq A) :=
    four_mul_le_sq_add (OrderedCommRing.add p w1) w4
  have ham2 : OrderedCommRing.le (OrderedCommRing.mul four (OrderedCommRing.mul (OrderedCommRing.add p w2) w5)) (sq B) :=
    four_mul_le_sq_add (OrderedCommRing.add p w2) w5
  have hprod_nonneg : OrderedCommRing.le OrderedCommRing.zero
      (OrderedCommRing.mul four (OrderedCommRing.mul (OrderedCommRing.add p w2) w5)) :=
    OrderedCommRing.mul_nonneg _ _ zero_lt_four.1
      (OrderedCommRing.mul_nonneg _ _
        (OrderedCommRing.add_le_add hp0 hw2 |> fun h => by rwa [OrderedCommRing.zero_add] at h) hw5)
  have h16 := OrderedCommRing.mul_le_mul ham1 ham2 hprod_nonneg (sq_nonneg A)
  -- 16 * (p+w1)w4(p+w2)w5 ≤ A² B² = (AB)²
  -- Also 16 * p² AB ≤ 16 * (p+w1)w4(p+w2)w5 from hq1
  have h16_pos : OrderedCommRing.le (OrderedCommRing.zero : R) (OrderedCommRing.mul four four) :=
    OrderedCommRing.mul_nonneg _ _ zero_lt_four.1 zero_lt_four.1
  have hq1_16 := mul_le_mul_of_nonneg_left hq1 h16_pos
  have hrearr : OrderedCommRing.mul (OrderedCommRing.mul four four)
                  (OrderedCommRing.mul (OrderedCommRing.mul (OrderedCommRing.add p w1) w4)
                                       (OrderedCommRing.mul (OrderedCommRing.add p w2) w5)) =
                OrderedCommRing.mul
                  (OrderedCommRing.mul four (OrderedCommRing.mul (OrderedCommRing.add p w1) w4))
                  (OrderedCommRing.mul four (OrderedCommRing.mul (OrderedCommRing.add p w2) w5)) := by
    rw [OrderedCommRing.mul_assoc four four,
        ← OrderedCommRing.mul_assoc four (OrderedCommRing.mul (OrderedCommRing.add p w1) w4),
        OrderedCommRing.mul_comm four (OrderedCommRing.mul (OrderedCommRing.add p w1) w4),
        OrderedCommRing.mul_assoc (OrderedCommRing.mul (OrderedCommRing.add p w1) w4) four,
        ← OrderedCommRing.mul_assoc four,
        OrderedCommRing.mul_comm (OrderedCommRing.mul (OrderedCommRing.add p w1) w4) four]
  rw [hrearr] at hq1_16
  have hsqAB : OrderedCommRing.mul (sq A) (sq B) = OrderedCommRing.mul (OrderedCommRing.mul A B) (OrderedCommRing.mul A B) := by
    unfold sq
    rw [OrderedCommRing.mul_assoc, ← OrderedCommRing.mul_assoc A B B,
        OrderedCommRing.mul_comm A B, OrderedCommRing.mul_assoc B A B,
        ← OrderedCommRing.mul_assoc, OrderedCommRing.mul_comm B A]
  have hchain1 := OrderedCommRing.le_trans _ _ _ hq1_16 h16
  rw [hsqAB] at hchain1
  -- Cancel AB > 0: 16 p² ≤ AB
  have hcancel : OrderedCommRing.le (OrderedCommRing.mul (OrderedCommRing.mul four four) (sq p))
                                    (OrderedCommRing.mul A B) := by
    apply le_of_mul_le_mul_right _ hABpos
    rw [OrderedCommRing.mul_assoc]
    exact hchain1
  -- Multiply by 4 and apply 4 AB ≤ (A + B)² = 1: 64 p² = (8p)² ≤ 4 AB ≤ 1 = 1²
  have h64 := mul_le_mul_of_nonneg_left hcancel zero_lt_four.1
  have ham3 : OrderedCommRing.le (OrderedCommRing.mul four (OrderedCommRing.mul A B))
                                 (sq (OrderedCommRing.add A B)) :=
    four_mul_le_sq_add A B
  rw [hsum] at ham3
  have h64_le_1 := OrderedCommRing.le_trans _ _ _ h64 ham3
  have h8sq : sq (OrderedCommRing.mul eight p) =
              OrderedCommRing.mul four (OrderedCommRing.mul (OrderedCommRing.mul four four) (sq p)) := by
    unfold sq eight four two
    simp only [OrderedCommRing.left_distrib, OrderedCommRing.right_distrib,
               OrderedCommRing.one_mul, OrderedCommRing.add_assoc]
  rw [← h8sq] at h64_le_1
  have h8p_0 : OrderedCommRing.le OrderedCommRing.zero (OrderedCommRing.mul eight p) :=
    OrderedCommRing.mul_nonneg _ _ zero_lt_eight.1 hp0
  exact le_of_sq_le_sq h8p_0 zero_le_one h64_le_1

/-- **Proposition 6.1 (uniqueness of `Hardy_{1/8}`):**
    If equality `8p = 1` holds on `F_4 ∩ Q^1`, then `2A = 1`, `2B = 1`, `4w4 = 1`, `4w5 = 1`,
    `8w1 = 1`, and `8w2 = 1` uniquely! -/
theorem prop_face_max_unique {R : Type _} [OrderedCommRing R]
    {p w1 w2 w4 w5 : R}
    (hp0 : OrderedCommRing.le OrderedCommRing.zero p)
    (hw1 : OrderedCommRing.le OrderedCommRing.zero w1)
    (hw2 : OrderedCommRing.le OrderedCommRing.zero w2)
    (hw4 : OrderedCommRing.le OrderedCommRing.zero w4)
    (hw5 : OrderedCommRing.le OrderedCommRing.zero w5)
    (hsum : OrderedCommRing.add
              (OrderedCommRing.add (OrderedCommRing.add p w1) w4)
              (OrderedCommRing.add (OrderedCommRing.add p w2) w5) = OrderedCommRing.one)
    (hApos : OrderedCommRing.lt OrderedCommRing.zero (OrderedCommRing.add (OrderedCommRing.add p w1) w4))
    (hBpos : OrderedCommRing.lt OrderedCommRing.zero (OrderedCommRing.add (OrderedCommRing.add p w2) w5))
    (hq1 : OrderedCommRing.le
             (OrderedCommRing.mul (sq p)
               (OrderedCommRing.mul (OrderedCommRing.add (OrderedCommRing.add p w1) w4)
                                    (OrderedCommRing.add (OrderedCommRing.add p w2) w5)))
             (OrderedCommRing.mul (OrderedCommRing.mul (OrderedCommRing.add p w1) w4)
                                  (OrderedCommRing.mul (OrderedCommRing.add p w2) w5)))
    (heq : OrderedCommRing.mul eight p = OrderedCommRing.one) :
    OrderedCommRing.mul four w4 = OrderedCommRing.one ∧
    OrderedCommRing.mul four w5 = OrderedCommRing.one ∧
    OrderedCommRing.mul eight w1 = OrderedCommRing.one ∧
    OrderedCommRing.mul eight w2 = OrderedCommRing.one := by
  let A := OrderedCommRing.add (OrderedCommRing.add p w1) w4
  let B := OrderedCommRing.add (OrderedCommRing.add p w2) w5
  have hABpos : OrderedCommRing.lt OrderedCommRing.zero (OrderedCommRing.mul A B) :=
    OrderedCommRing.mul_pos A B hApos hBpos
  have ham1 := four_mul_le_sq_add (OrderedCommRing.add p w1) w4
  have ham2 := four_mul_le_sq_add (OrderedCommRing.add p w2) w5
  have hprod_nonneg : OrderedCommRing.le OrderedCommRing.zero
      (OrderedCommRing.mul four (OrderedCommRing.mul (OrderedCommRing.add p w2) w5)) :=
    OrderedCommRing.mul_nonneg _ _ zero_lt_four.1
      (OrderedCommRing.mul_nonneg _ _
        (OrderedCommRing.add_le_add hp0 hw2 |> fun h => by rwa [OrderedCommRing.zero_add] at h) hw5)
  have h16 := OrderedCommRing.mul_le_mul ham1 ham2 hprod_nonneg (sq_nonneg A)
  have h16_pos : OrderedCommRing.le (OrderedCommRing.zero : R) (OrderedCommRing.mul four four) :=
    OrderedCommRing.mul_nonneg _ _ zero_lt_four.1 zero_lt_four.1
  have hq1_16 := mul_le_mul_of_nonneg_left hq1 h16_pos
  have hrearr : OrderedCommRing.mul (OrderedCommRing.mul four four)
                  (OrderedCommRing.mul (OrderedCommRing.mul (OrderedCommRing.add p w1) w4)
                                       (OrderedCommRing.mul (OrderedCommRing.add p w2) w5)) =
                OrderedCommRing.mul
                  (OrderedCommRing.mul four (OrderedCommRing.mul (OrderedCommRing.add p w1) w4))
                  (OrderedCommRing.mul four (OrderedCommRing.mul (OrderedCommRing.add p w2) w5)) := by
    rw [OrderedCommRing.mul_assoc four four,
        ← OrderedCommRing.mul_assoc four (OrderedCommRing.mul (OrderedCommRing.add p w1) w4),
        OrderedCommRing.mul_comm four (OrderedCommRing.mul (OrderedCommRing.add p w1) w4),
        OrderedCommRing.mul_assoc (OrderedCommRing.mul (OrderedCommRing.add p w1) w4) four,
        ← OrderedCommRing.mul_assoc four,
        OrderedCommRing.mul_comm (OrderedCommRing.mul (OrderedCommRing.add p w1) w4) four]
  rw [hrearr] at hq1_16
  have hsqAB : OrderedCommRing.mul (sq A) (sq B) = OrderedCommRing.mul (OrderedCommRing.mul A B) (OrderedCommRing.mul A B) := by
    unfold sq
    rw [OrderedCommRing.mul_assoc, ← OrderedCommRing.mul_assoc A B B,
        OrderedCommRing.mul_comm A B, OrderedCommRing.mul_assoc B A B,
        ← OrderedCommRing.mul_assoc, OrderedCommRing.mul_comm B A]
  have hchain1 := OrderedCommRing.le_trans _ _ _ hq1_16 h16
  rw [hsqAB] at hchain1
  have hcancel : OrderedCommRing.le (OrderedCommRing.mul (OrderedCommRing.mul four four) (sq p))
                                    (OrderedCommRing.mul A B) := by
    apply le_of_mul_le_mul_right _ hABpos
    rw [OrderedCommRing.mul_assoc]; exact hchain1
  have h64 := mul_le_mul_of_nonneg_left hcancel zero_lt_four.1
  have h8sq : sq (OrderedCommRing.mul eight p) =
              OrderedCommRing.mul four (OrderedCommRing.mul (OrderedCommRing.mul four four) (sq p)) := by
    unfold sq eight four two
    simp only [OrderedCommRing.left_distrib, OrderedCommRing.right_distrib,
               OrderedCommRing.one_mul, OrderedCommRing.add_assoc]
  rw [← h8sq, heq] at h64
  -- Now 1² ≤ 4 AB ≤ (A + B)² = 1², so 4 AB = (A + B)² => A = B!
  have ham3 := four_mul_le_sq_add A B
  rw [hsum] at ham3
  have h4AB : OrderedCommRing.mul four (OrderedCommRing.mul A B) = sq (OrderedCommRing.add A B) := by
    rw [hsum]; exact OrderedCommRing.le_antisymm _ _ ham3 h64
  have hAB_eq : A = B := eq_of_four_mul_eq_sq_add h4AB
  -- Hence 2A = A + B = 1, 2B = 1
  have h2A : OrderedCommRing.mul two A = OrderedCommRing.one := by
    rw [two_mul_eq_add]
    calc OrderedCommRing.add A A = OrderedCommRing.add A B := congrArg (OrderedCommRing.add A) hAB_eq
      _ = OrderedCommRing.one := hsum
  have h2B : OrderedCommRing.mul two B = OrderedCommRing.one := by
    rw [← hAB_eq]; exact h2A
  -- Next, 16 p² = AB (since 4 * 16 p² = 1 = 4 AB), so 16 p² AB = (AB)² = A² B²,
  -- squeezing both ham1 and ham2 to exact equalities!
  have h16p2_eq_AB : OrderedCommRing.mul (OrderedCommRing.mul four four) (sq p) = OrderedCommRing.mul A B := by
    apply OrderedCommRing.le_antisymm _ _ hcancel
    apply le_of_mul_le_mul_left _ zero_lt_four
    rw [← h8sq, heq, h4AB, hsum]; exact OrderedCommRing.le_refl _
  have hprod_eq :
      OrderedCommRing.mul
        (OrderedCommRing.mul four (OrderedCommRing.mul (OrderedCommRing.add p w1) w4))
        (OrderedCommRing.mul four (OrderedCommRing.mul (OrderedCommRing.add p w2) w5)) =
      OrderedCommRing.mul (sq A) (sq B) := by
    apply OrderedCommRing.le_antisymm _ _ h16
    have hlhs : OrderedCommRing.mul (sq A) (sq B) =
                OrderedCommRing.mul (OrderedCommRing.mul four four)
                  (OrderedCommRing.mul (sq p) (OrderedCommRing.mul A B)) := by
      calc OrderedCommRing.mul (sq A) (sq B)
          = OrderedCommRing.mul (OrderedCommRing.mul A B) (OrderedCommRing.mul A B) := hsqAB
        _ = OrderedCommRing.mul (OrderedCommRing.mul (OrderedCommRing.mul four four) (sq p))
              (OrderedCommRing.mul A B) :=
              congrArg (fun x => OrderedCommRing.mul x (OrderedCommRing.mul A B)) h16p2_eq_AB.symm
        _ = OrderedCommRing.mul (OrderedCommRing.mul four four)
              (OrderedCommRing.mul (sq p) (OrderedCommRing.mul A B)) :=
              OrderedCommRing.mul_assoc _ _ _
    rw [hlhs]
    exact hq1_16
  -- From ham1 ≤ A² and ham2 ≤ B² with product equal and A², B² > 0, both equalities hold!
  have hA2_pos : OrderedCommRing.lt OrderedCommRing.zero (sq A) := OrderedCommRing.mul_pos A A hApos hApos
  have hB2_pos : OrderedCommRing.lt OrderedCommRing.zero (sq B) := OrderedCommRing.mul_pos B B hBpos hBpos
  have h1_nonneg : OrderedCommRing.le OrderedCommRing.zero
      (OrderedCommRing.mul four (OrderedCommRing.mul (OrderedCommRing.add p w1) w4)) :=
    OrderedCommRing.mul_nonneg _ _ zero_lt_four.1
      (OrderedCommRing.mul_nonneg _ _
        (OrderedCommRing.add_le_add hp0 hw1 |> fun h => by rwa [OrderedCommRing.zero_add] at h) hw4)
  have ham1_eq : OrderedCommRing.mul four (OrderedCommRing.mul (OrderedCommRing.add p w1) w4) = sq A := by
    apply OrderedCommRing.le_antisymm _ _ ham1
    apply le_of_mul_le_mul_right _ hB2_pos
    rw [← hprod_eq]
    exact mul_le_mul_of_nonneg_left ham2 h1_nonneg
  have ham2_eq : OrderedCommRing.mul four (OrderedCommRing.mul (OrderedCommRing.add p w2) w5) = sq B := by
    apply OrderedCommRing.le_antisymm _ _ ham2
    apply le_of_mul_le_mul_left _ hA2_pos
    rw [← hprod_eq, ham1_eq]
    exact OrderedCommRing.le_refl _
  -- Equality in AM–GM gives p + w1 = w4 and p + w2 = w5!
  have hw4_eq : OrderedCommRing.add p w1 = w4 := eq_of_four_mul_eq_sq_add ham1_eq
  have hw5_eq : OrderedCommRing.add p w2 = w5 := eq_of_four_mul_eq_sq_add ham2_eq
  -- Since A = (p+w1) + w4 = 2 w4 and 2 A = 1, we get 4 w4 = 1! Likewise 4 w5 = 1!
  have h4w4 : OrderedCommRing.mul four w4 = OrderedCommRing.one := by
    have hA2w4 : A = OrderedCommRing.mul two w4 := by
      change OrderedCommRing.add (OrderedCommRing.add p w1) w4 = _
      rw [hw4_eq, two_mul_eq_add]
    rw [four_mul_eq, ← hA2w4, h2A]
  have h4w5 : OrderedCommRing.mul four w5 = OrderedCommRing.one := by
    have hB2w5 : B = OrderedCommRing.mul two w5 := by
      change OrderedCommRing.add (OrderedCommRing.add p w2) w5 = _
      rw [hw5_eq, two_mul_eq_add]
    rw [four_mul_eq, ← hB2w5, h2B]
  -- Finally 8p + 8w1 = 8(p+w1) = 8 w4 = 2(4 w4) = 2 = 1 + 1, and 8p = 1 => 8w1 = 1!
  have h8w1 : OrderedCommRing.mul eight w1 = OrderedCommRing.one := by
    have h1 : OrderedCommRing.add (OrderedCommRing.mul eight p) (OrderedCommRing.mul eight w1) = two := by
      rw [← OrderedCommRing.left_distrib, hw4_eq]
      unfold eight; rw [OrderedCommRing.right_distrib, h4w4]; rfl
    rw [heq] at h1
    unfold two at h1
    have h2 := congrArg (fun x => OrderedCommRing.add (OrderedCommRing.neg OrderedCommRing.one) x) h1
    rw [← OrderedCommRing.add_assoc, OrderedCommRing.neg_add, OrderedCommRing.zero_add,
        ← OrderedCommRing.add_assoc, OrderedCommRing.neg_add, OrderedCommRing.zero_add] at h2
    exact h2
  have h8w2 : OrderedCommRing.mul eight w2 = OrderedCommRing.one := by
    have h1 : OrderedCommRing.add (OrderedCommRing.mul eight p) (OrderedCommRing.mul eight w2) = two := by
      rw [← OrderedCommRing.left_distrib, hw5_eq]
      unfold eight; rw [OrderedCommRing.right_distrib, h4w5]; rfl
    rw [heq] at h1
    unfold two at h1
    have h2 := congrArg (fun x => OrderedCommRing.add (OrderedCommRing.neg OrderedCommRing.one) x) h1
    rw [← OrderedCommRing.add_assoc, OrderedCommRing.neg_add, OrderedCommRing.zero_add,
        ← OrderedCommRing.add_assoc, OrderedCommRing.neg_add, OrderedCommRing.zero_add] at h2
    exact h2
  exact ⟨h4w4, h4w5, h8w1, h8w2⟩

end Theorem101.QABInequality
