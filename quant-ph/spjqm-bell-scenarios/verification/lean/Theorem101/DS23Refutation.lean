/-
  Theorem101.DS23Refutation — Formal verification of:
    * Theorem 5.2(a): Exact rank-4 rational Gram matrix placing Hardy_{1/8} in SPJQM
    * Theorem 5.2(b): Abstract 4-line operator impossibility (no PoZ+LoN model for Hardy_{1/8})
    * Corollary 6.5 (left): Exact 64-atom rank-4 Gram matrix placing Hardy_{1/8} ⊗ det_C in SPJQM_3
    * Theorem 5.1(d): Exact J3 perturbation Δ vanishes on all 64 context blocks & branching = -2/125
-/
import Theorem101.Basic

namespace Theorem101.DS23Refutation
open Theorem101 AddCommGrp

/-! ### Part 1: Theorem 5.2(b) — 4-line operator impossibility (zero axioms, pure algebra!) -/

/-- An additive endomorphism on `H` (event operator in a DS23 PoZ+LoN model).
    Note: neither Hermiticity, idempotence, nor any norm bound is assumed! -/
structure EventOp (H : Type _) [AddCommGrp H] where
  fn : H → H
  map_zero : fn zero = zero
  map_add  : ∀ u v : H, fn (add u v) = add (fn u) (fn v)

theorem EventOp.map_neg {H : Type _} [AddCommGrp H] (E : EventOp H) (u : H) :
    E.fn (neg u) = neg (E.fn u) := by
  have h : add (E.fn (neg u)) (E.fn u) = zero := by
    rw [← E.map_add, neg_add, E.map_zero]
  have h2 := congrArg (fun x => add x (neg (E.fn u))) h
  rw [add_assoc, add_neg, add_zero, zero_add] at h2
  exact h2

theorem EventOp.map_sub {H : Type _} [AddCommGrp H] (E : EventOp H) (u v : H) :
    E.fn (sub u v) = sub (E.fn u) (E.fn v) := by
  unfold sub; rw [E.map_add, E.map_neg]

/-- **Theorem 5.2(b) (4-line refutation of DS23's Converse Theorem):**
    In any DS23 PoZ+LoN model with commuting linear event operators summing to `I`,
    the four Hardy zeros `Ap0 Bp0 Ω = 0`, `Ap1 Bp1 Ω = 0`, `A0 Bp1 Ω = 0`, `Ap1 B0 Ω = 0`
    force `A0 B0 Ω = 0` (hence `P(00|ab) = 0 ≠ 1/8`). -/
theorem thm_ds23_no_model {H : Type _} [AddCommGrp H]
    (A0 Ap0 Ap1 B0 Bp0 Bp1 : EventOp H) (Ω : H)
    (hAp : ∀ v, add (Ap0.fn v) (Ap1.fn v) = v)
    (hBp : ∀ v, add (Bp0.fn v) (Bp1.fn v) = v)
    (hcomm_A0_B0  : ∀ v, A0.fn (B0.fn v) = B0.fn (A0.fn v))
    (hcomm_Ap1_B0 : ∀ v, Ap1.fn (B0.fn v) = B0.fn (Ap1.fn v))
    (hz1 : Ap0.fn (Bp0.fn Ω) = zero)
    (hz2 : Ap1.fn (Bp1.fn Ω) = zero)
    (hz3 : A0.fn  (Bp1.fn Ω) = zero)
    (hz4 : Ap1.fn (B0.fn Ω)  = zero) :
    A0.fn (B0.fn Ω) = zero := by
  -- Line 1: Ap0 Ω = Bp1 Ω
  have hstep1 : Ap0.fn Ω = Bp1.fn Ω := by
    calc Ap0.fn Ω = Ap0.fn (add (Bp0.fn Ω) (Bp1.fn Ω)) := by rw [hBp]
      _ = add (Ap0.fn (Bp0.fn Ω)) (Ap0.fn (Bp1.fn Ω)) := Ap0.map_add _ _
      _ = Ap0.fn (Bp1.fn Ω) := by rw [hz1, zero_add]
      _ = add (Ap0.fn (Bp1.fn Ω)) (Ap1.fn (Bp1.fn Ω)) := by rw [hz2, add_zero]
      _ = Bp1.fn Ω := hAp (Bp1.fn Ω)
  -- Line 2: A0 (Ap0 Ω) = 0
  have hstep2 : A0.fn (Ap0.fn Ω) = zero := by
    rw [hstep1, hz3]
  -- Line 4: A0 (B0 Ω) = B0 (A0 (Ap0 Ω)) + A0 (Ap1 (B0 Ω)) = 0 + 0 = 0
  have hstep4 : A0.fn (B0.fn Ω) = add (B0.fn (A0.fn (Ap0.fn Ω))) (A0.fn (Ap1.fn (B0.fn Ω))) := by
    calc A0.fn (B0.fn Ω)
        = A0.fn (B0.fn (add (Ap0.fn Ω) (Ap1.fn Ω))) := by rw [hAp]
      _ = add (A0.fn (B0.fn (Ap0.fn Ω))) (A0.fn (B0.fn (Ap1.fn Ω))) := by
            rw [B0.map_add, A0.map_add]
      _ = add (B0.fn (A0.fn (Ap0.fn Ω))) (A0.fn (Ap1.fn (B0.fn Ω))) := by
            rw [hcomm_A0_B0 (Ap0.fn Ω), ← hcomm_Ap1_B0 Ω]
  rw [hstep4, hstep2, hz4, B0.map_zero, A0.map_zero, zero_add]

/-! ### Part 2: Theorem 5.2(a) — Exact 16×16 rank-4 rational Gram matrix for `Hardy_{1/8}` -/

/-- 4D integer vectors (scaled by 64 so every atom vector has exact integer coordinates).
    The original inner product has weights `w = (1/8, 1/8, 1/8, 5/128) = (16, 16, 16, 5) / 128`.
    Scaling vectors by `64` and weights by `128` multiplies every inner product by
    `128 * 64^2 = 524288`. -/
abbrev Vec4 := Int × Int × Int × Int

def v4Add (u v : Vec4) : Vec4 :=
  (u.1 + v.1, u.2.1 + v.2.1, u.2.2.1 + v.2.2.1, u.2.2.2 + v.2.2.2)

def v4Sub (u v : Vec4) : Vec4 :=
  (u.1 - v.1, u.2.1 - v.2.1, u.2.2.1 - v.2.2.1, u.2.2.2 - v.2.2.2)

def v4Scale (c : Int) (u : Vec4) : Vec4 :=
  (c * u.1, c * u.2.1, c * u.2.2.1, c * u.2.2.2)

def ip4 (u v : Vec4) : Int :=
  16 * u.1 * v.1 + 16 * u.2.1 * v.2.1 + 16 * u.2.2.1 * v.2.2.1 + 5 * u.2.2.2 * v.2.2.2

/-- All four diagonal Gram weights `(16, 16, 16, 5)` are strictly positive, so `D_joint`
    is positive semidefinite of rank ≤ 4 by construction. -/
theorem hardy18_weights_pos : (16 : Int) > 0 ∧ (16 : Int) > 0 ∧ (16 : Int) > 0 ∧ (5 : Int) > 0 := by
  decide

-- Exact scaled-by-4 marginal/context vectors (multiplied by 16 again in the S4 formula -> scale 64):
def psi4 : Vec4 := (8, 8, 0, 0)
def a0_4 : Vec4 := (0, 4, 4, 0)
def a1_4 : Vec4 := v4Sub psi4 a0_4
def ap0_4 : Vec4 := (8, 0, 0, 0)
def ap1_4 : Vec4 := v4Sub psi4 ap0_4
def b0_4 : Vec4 := (4, 0, 4, 0)
def b1_4 : Vec4 := v4Sub psi4 b0_4
def bp0_4 : Vec4 := (0, 8, 0, 0)
def bp1_4 : Vec4 := v4Sub psi4 bp0_4

def A_vec (x o : Fin 2) : Vec4 :=
  match x.val, o.val with
  | 0, 0 => a0_4
  | 0, _ => a1_4
  | _, 0 => ap0_4
  | _, _ => ap1_4

def B_vec (y o : Fin 2) : Vec4 :=
  match y.val, o.val with
  | 0, 0 => b0_4
  | 0, _ => b1_4
  | _, 0 => bp0_4
  | _, _ => bp1_4

def W_vec (x y a b : Fin 2) : Vec4 :=
  let w00 : Vec4 :=
    match x.val, y.val with
    | 0, 0 => (1, 1, 3, 4)   -- 4 * (1/4, 1/4, 3/4, 1)
    | 0, _ => a0_4
    | _, 0 => b0_4
    | _, _ => (0, 0, 0, 0)
  let w01 := v4Sub (A_vec x 0) w00
  let w10 := v4Sub (B_vec y 0) w00
  let w11 := v4Sub (B_vec y 1) w01
  match a.val, b.val with
  | 0, 0 => w00
  | 0, _ => w01
  | _, 0 => w10
  | _, _ => w11

/-- S4 atom vector `64 * |i, i', j, j'⟩ = 4*(W00+W01+W10+W11) - 2*(A0+A1+B0+B1) + psi4`. -/
def hardyAtom (g : Fin 16) : Vec4 :=
  let i   : Fin 2 := ⟨(g.val / 8) % 2, by omega⟩
  let ip_ : Fin 2 := ⟨(g.val / 4) % 2, by omega⟩
  let j   : Fin 2 := ⟨(g.val / 2) % 2, by omega⟩
  let jp  : Fin 2 := ⟨g.val % 2, by omega⟩
  let sW := v4Add (v4Add (W_vec 0 0 i j) (W_vec 0 1 i jp))
                  (v4Add (W_vec 1 0 ip_ j) (W_vec 1 1 ip_ jp))
  let sM := v4Add (v4Add (A_vec 0 i) (A_vec 1 ip_))
                  (v4Add (B_vec 0 j) (B_vec 1 jp))
  v4Add (v4Sub (v4Scale 4 sW) (v4Scale 2 sM)) psi4

/-- Sum of `hardyAtom` over all atoms in cylinder `[a, b | x, y]`. -/
def hardyCellVec (x y a b : Fin 2) : Vec4 :=
  (List.finRange 16).foldl (fun acc g =>
    let ai := if x.val == 0 then (g.val >>> 3) &&& 1 else (g.val >>> 2) &&& 1
    let bj := if y.val == 0 then (g.val >>> 1) &&& 1 else g.val &&& 1
    if ai == a.val && bj == b.val then v4Add acc (hardyAtom g) else acc) (0, 0, 0, 0)

/-- Target `Hardy_{1/8}` probabilities scaled by 8 (exact integers). -/
def hardyTarget8 (x y a b : Fin 2) : Int :=
  match x.val, y.val, a.val, b.val with
  | 0, 0, 0, 0 => 1
  | 0, 0, 0, 1 => 1
  | 0, 0, 1, 0 => 1
  | 0, 0, 1, 1 => 5
  | 0, 1, 0, 0 => 2
  | 0, 1, 0, 1 => 0
  | 0, 1, 1, 0 => 2
  | 0, 1, 1, 1 => 4
  | 1, 0, 0, 0 => 2
  | 1, 0, 0, 1 => 2
  | 1, 0, 1, 0 => 0
  | 1, 0, 1, 1 => 4
  | 1, 1, 0, 0 => 0
  | 1, 1, 0, 1 => 4
  | 1, 1, 1, 0 => 4
  | _, _, _, _ => 0

/-- **Theorem 5.2(a) (kernel-verified):** All 64 context decoherence entries
    `D_joint([a1 b1 | x y], [a2 b2 | x y])` equal `δ_{a1,a2} δ_{b1,b2} Hardy_{1/8}(a1 b1 | x y)`
    (scaled by `524288 = 8 * 65536`). -/
theorem hardy18_spjqm_exact :
    ∀ x y a1 b1 a2 b2 : Fin 2,
      ip4 (hardyCellVec x y a1 b1) (hardyCellVec x y a2 b2) =
        (if a1 == a2 && b1 == b2 then 65536 * hardyTarget8 x y a1 b1 else 0) := by
  decide

/-! ### Part 3: Corollary 6.5 (left) — 64-atom 3-party `Hardy_{1/8} ⊗ det_C ∈ SPJQM_3` -/

/-- 64-atom vector for `Hardy_{1/8} ⊗ det_C`: `hardyAtom` when `c0 = c1 = 0`, else `0`. -/
def hardy3Atom (g : Fin 64) : Vec4 :=
  if (g.val &&& 3) == 0 then hardyAtom ⟨g.val >>> 2, by omega⟩ else (0, 0, 0, 0)

def hardy3CellVec (x y z a b c : Fin 2) : Vec4 :=
  (List.finRange 64).foldl (fun acc g =>
    let ai := if x.val == 0 then (g.val >>> 5) &&& 1 else (g.val >>> 4) &&& 1
    let bj := if y.val == 0 then (g.val >>> 3) &&& 1 else (g.val >>> 2) &&& 1
    let ck := if z.val == 0 then (g.val >>> 1) &&& 1 else g.val &&& 1
    if ai == a.val && bj == b.val && ck == c.val then v4Add acc (hardy3Atom g) else acc) (0, 0, 0, 0)

/-- **Corollary 6.5 (left, kernel-verified):** All 512 three-party context decoherence entries
    of `Hardy_{1/8} ⊗ det_C` match exactly. -/
theorem hardy18_3party_spjqm_exact :
    ∀ x y z a1 b1 c1 a2 b2 c2 : Fin 2,
      ip4 (hardy3CellVec x y z a1 b1 c1) (hardy3CellVec x y z a2 b2 c2) =
        (if a1 == a2 && b1 == b2 && c1 == c2 && c1.val == 0
         then 65536 * hardyTarget8 x y a1 b1 else 0) := by
  decide

/-! ### Part 4: Theorem 5.1(d) — Exact `J3` perturbation `Δ` & branching entry `-2/125` -/

def s13 (g : Fin 16) : Int :=
  if (((g.val >>> 3) &&& 1) + ((g.val >>> 1) &&& 1)) % 2 == 0 then 1 else -1

def s14 (g : Fin 16) : Int :=
  if (((g.val >>> 3) &&& 1) + (g.val &&& 1)) % 2 == 0 then 1 else -1

def deltaEntry (g1 g2 : Fin 16) : Int :=
  s13 g1 * s14 g2 + s14 g1 * s13 g2

def inCell (g : Fin 16) (x y a b : Fin 2) : Bool :=
  let ai := if x.val == 0 then (g.val >>> 3) &&& 1 else (g.val >>> 2) &&& 1
  let bj := if y.val == 0 then (g.val >>> 1) &&& 1 else g.val &&& 1
  ai == a.val && bj == b.val

/-- `Δ` sums to `0` on every one of the 64 context blocks `(x, y, a1, b1, a2, b2)`. -/
theorem j3_delta_context_zero :
    ∀ x y a1 b1 a2 b2 : Fin 2,
      (List.finRange 16).foldl (fun acc1 g1 =>
        if inCell g1 x y a1 b1 then
          acc1 + (List.finRange 16).foldl (fun acc2 g2 =>
            if inCell g2 x y a2 b2 then acc2 + deltaEntry g1 g2 else acc2) 0
        else acc1) 0 = 0 := by
  decide

/-- The branching off-diagonal sum of `Δ` on `(a=0 b=0, a=1 b'=0)` equals `-16`,
    so `Δ / 1000` contributes `-16 / 1000 = -2 / 125` (while `J2` vanishes on disjoint `a=0, a=1`). -/
theorem j3_delta_branching_sum :
    (List.finRange 16).foldl (fun acc1 g1 =>
      if ((g1.val >>> 3) &&& 1) == 0 && ((g1.val >>> 1) &&& 1) == 0 then
        acc1 + (List.finRange 16).foldl (fun acc2 g2 =>
          if ((g2.val >>> 3) &&& 1) == 1 && (g2.val &&& 1) == 0 then
            acc2 + deltaEntry g1 g2
          else acc2) 0
      else acc1) 0 = -16 := by
  decide

end Theorem101.DS23Refutation
