/-
  Theorem101.BipartiteLift — Formal verification of the marginal-recovery identity
  underlying Lemma 4.2 (Step (a) for arbitrary number of settings `m ≥ 1` and arbitrary
  finite outcome counts `dims : List Nat`) and Theorem 3.2 (Step 4, Claim A).
  Works over ANY `AddCommGrp V` (using natural-number / integer scalar multiplication `nsmul`,
  avoiding division by clearing the denominator `N₁ = ∏ dims`)!
-/
import Theorem101.Basic

namespace Theorem101.BipartiteLift
open Theorem101 AddCommGrp

/-- Natural-number scalar multiplication on any `AddCommGrp V`. -/
def nsmul {V : Type _} [AddCommGrp V] : Nat → V → V
  | 0, _     => zero
  | n + 1, v => add (nsmul n v) v

theorem nsmul_zero {V : Type _} [AddCommGrp V] (n : Nat) : nsmul n (zero : V) = zero := by
  induction n with
  | zero => rfl
  | succ k ih => simp [nsmul, ih, zero_add]

theorem nsmul_add {V : Type _} [AddCommGrp V] (n : Nat) (u v : V) :
    nsmul n (add u v) = add (nsmul n u) (nsmul n v) := by
  induction n with
  | zero => simp [nsmul, zero_add]
  | succ k ih =>
    simp only [nsmul, ih]
    rw [add_assoc (nsmul k u) (nsmul k v) (add u v),
        ← add_assoc (nsmul k v) u v,
        add_comm (nsmul k v) u,
        add_assoc u (nsmul k v) v,
        ← add_assoc (nsmul k u) u _]

theorem nsmul_add_nat {V : Type _} [AddCommGrp V] (m n : Nat) (v : V) :
    nsmul (m + n) v = add (nsmul m v) (nsmul n v) := by
  induction n with
  | zero => simp [nsmul, add_zero]
  | succ k ih =>
    rw [Nat.add_succ, nsmul, ih, nsmul, add_assoc]

theorem nsmul_mul {V : Type _} [AddCommGrp V] (m n : Nat) (v : V) :
    nsmul (m * n) v = nsmul m (nsmul n v) := by
  induction m with
  | zero => simp [nsmul]
  | succ k ih =>
    rw [Nat.succ_mul, nsmul_add_nat, ih, nsmul]

theorem nsmul_neg {V : Type _} [AddCommGrp V] (n : Nat) (v : V) :
    nsmul n (neg v) = neg (nsmul n v) := by
  induction n with
  | zero => simp [nsmul, neg_zero]
  | succ k ih => simp only [nsmul, ih, neg_add_rev]

theorem nsmul_sub {V : Type _} [AddCommGrp V] (n : Nat) (u v : V) :
    nsmul n (sub u v) = sub (nsmul n u) (nsmul n v) := by
  unfold sub; rw [nsmul_add, nsmul_neg]

/-- Finite sum of `f : Fin n → V` in an `AddCommGrp V`. -/
def finSum {V : Type _} [AddCommGrp V] : {n : Nat} → (Fin n → V) → V
  | 0, _     => zero
  | _ + 1, f => add (finSum (fun i => f i.castSucc)) (f (Fin.last _))

theorem finSum_const {V : Type _} [AddCommGrp V] (n : Nat) (v : V) :
    finSum (fun (_ : Fin n) => v) = nsmul n v := by
  induction n with
  | zero => rfl
  | succ k ih => simp only [finSum, ih, nsmul]

theorem finSum_add {V : Type _} [AddCommGrp V] {n : Nat} (f g : Fin n → V) :
    finSum (fun i => add (f i) (g i)) = add (finSum f) (finSum g) := by
  induction n with
  | zero => simp [finSum, zero_add]
  | succ k ih =>
    simp only [finSum, ih]
    rw [add_assoc (finSum fun i => f i.castSucc),
        ← add_assoc (finSum fun i => g i.castSucc) (f (Fin.last k)),
        add_comm (finSum fun i => g i.castSucc) (f (Fin.last k)),
        add_assoc (f (Fin.last k)),
        ← add_assoc (finSum fun i => f i.castSucc)]

theorem finSum_neg {V : Type _} [AddCommGrp V] {n : Nat} (f : Fin n → V) :
    finSum (fun i => neg (f i)) = neg (finSum f) := by
  induction n with
  | zero => simp [finSum, neg_zero]
  | succ k ih => simp only [finSum, ih, neg_add_rev]

theorem finSum_sub {V : Type _} [AddCommGrp V] {n : Nat} (f g : Fin n → V) :
    finSum (fun i => sub (f i) (g i)) = sub (finSum f) (finSum g) := by
  unfold sub; rw [finSum_add, finSum_neg]

theorem finSum_nsmul {V : Type _} [AddCommGrp V] {n : Nat} (c : Nat) (f : Fin n → V) :
    finSum (fun i => nsmul c (f i)) = nsmul c (finSum f) := by
  induction n with
  | zero => simp [finSum, nsmul_zero]
  | succ k ih => simp only [finSum, ih, nsmul_add]

theorem finSum_swap {V : Type _} [AddCommGrp V] {m n : Nat} (f : Fin m → Fin n → V) :
    finSum (fun i => finSum (fun j => f i j)) = finSum (fun j => finSum (fun i => f i j)) := by
  induction m with
  | zero =>
    simp only [finSum]
    induction n with
    | zero => rfl
    | succ k ih =>
      rw [finSum]
      have h0 := ih (fun i j => f i j.castSucc)
      change zero = add (finSum (fun j => zero)) zero
      rw [← h0, add_zero]
  | succ k ih =>
    simp only [finSum, ih, finSum_add]

/-- A tuple of outcomes across settings with outcome counts `dims : List Nat`. -/
def Tuple : List Nat → Type
  | []      => Unit
  | d :: ds => Fin d × Tuple ds

/-- Product of all outcome counts `N₁ = ∏ d`. -/
def prodDims : List Nat → Nat
  | []      => 1
  | d :: ds => d * prodDims ds

/-- Sum a function `F : Tuple dims → V` over ALL tuples in `Tuple dims`. -/
def sumTuple {V : Type _} [AddCommGrp V] : {dims : List Nat} → (Tuple dims → V) → V
  | [],      F => F ()
  | _ :: _,  F => finSum (fun a => sumTuple (fun rest => F (a, rest)))

/-- Sum a function `F : Tuple dims → V` over all tuples whose coordinate `x0` is fixed to `a0`. -/
def sumTupleFix {V : Type _} [AddCommGrp V] :
    {dims : List Nat} → (x0 : Fin dims.length) → Fin (dims.get x0) → (Tuple dims → V) → V
  | [],     x0, _, _ => Fin.elim0 x0
  | _ :: _, ⟨0, _⟩,     a0, F => sumTuple (fun rest => F (a0, rest))
  | _ :: _, ⟨k + 1, h⟩, a0, F =>
      finSum (fun a => sumTupleFix ⟨k, Nat.lt_of_succ_lt_succ h⟩ a0 (fun rest => F (a, rest)))

/-- Coordinate projection `proj x t : Fin (dims.get x)` from `t : Tuple dims`. -/
def proj : {dims : List Nat} → (x : Fin dims.length) → Tuple dims → Fin (dims.get x)
  | [],     x, _         => Fin.elim0 x
  | _ :: _, ⟨0, _⟩,     (a, _)    => a
  | _ :: _, ⟨k + 1, h⟩, (_, rest) => proj ⟨k, Nat.lt_of_succ_lt_succ h⟩ rest

theorem sumTuple_const {V : Type _} [AddCommGrp V] (dims : List Nat) (v : V) :
    sumTuple (dims := dims) (fun _ => v) = nsmul (prodDims dims) v := by
  induction dims with
  | nil => simp [sumTuple, prodDims, nsmul, zero_add]
  | cons d ds ih =>
    simp only [sumTuple, prodDims, ih, finSum_const, nsmul_mul]

theorem sumTuple_add {V : Type _} [AddCommGrp V] {dims : List Nat} (F G : Tuple dims → V) :
    sumTuple (fun t => add (F t) (G t)) = add (sumTuple F) (sumTuple G) := by
  induction dims with
  | nil => rfl
  | cons d ds ih => simp only [sumTuple, ih, finSum_add]

theorem sumTupleFix_add {V : Type _} [AddCommGrp V] {dims : List Nat}
    (x0 : Fin dims.length) (a0 : Fin (dims.get x0)) (F G : Tuple dims → V) :
    sumTupleFix x0 a0 (fun t => add (F t) (G t)) =
    add (sumTupleFix x0 a0 F) (sumTupleFix x0 a0 G) := by
  induction dims with
  | nil => exact Fin.elim0 x0
  | cons d ds ih =>
    cases x0 with
    | mk idx h =>
      cases idx with
      | zero => exact sumTuple_add _ _
      | succ k => simp only [sumTupleFix, ih, finSum_add]

theorem sumTupleFix_sub {V : Type _} [AddCommGrp V] {dims : List Nat}
    (x0 : Fin dims.length) (a0 : Fin (dims.get x0)) (F G : Tuple dims → V) :
    sumTupleFix x0 a0 (fun t => sub (F t) (G t)) =
    sub (sumTupleFix x0 a0 F) (sumTupleFix x0 a0 G) := by
  induction dims with
  | nil => exact Fin.elim0 x0
  | cons d ds ih =>
    cases x0 with
    | mk idx h =>
      cases idx with
      | zero =>
        unfold sumTupleFix sub
        have hneg : ∀ {ds' : List Nat} (H : Tuple ds' → V),
            sumTuple (fun t => neg (H t)) = neg (sumTuple H) := by
          intro ds'
          induction ds' with
          | nil => intro _; rfl
          | cons d' ds'' ih' => intro H; simp only [sumTuple, ih', finSum_neg]
        rw [sumTuple_add, hneg]
      | succ k => simp only [sumTupleFix, ih, finSum_sub]

/-- Helper: product of `dims` divided by `dims.get x0` (computed cleanly as product of other entries). -/
def prodExcept : {dims : List Nat} → Fin dims.length → Nat
  | [],     x0          => Fin.elim0 x0
  | _ :: ds, ⟨0, _⟩     => prodDims ds
  | d :: _, ⟨k + 1, h⟩  => d * prodExcept ⟨k, Nat.lt_of_succ_lt_succ h⟩

theorem dims_mul_prodExcept {dims : List Nat} (x0 : Fin dims.length) :
    dims.get x0 * prodExcept x0 = prodDims dims := by
  induction dims with
  | nil => exact Fin.elim0 x0
  | cons d ds ih =>
    cases x0 with
    | mk idx h =>
      cases idx with
      | zero => rfl
      | succ k =>
        simp only [List.get, prodExcept, prodDims]
        rw [← Nat.mul_assoc, Nat.mul_comm (ds.get ⟨k, _⟩) d, Nat.mul_assoc, ih]

theorem sumTupleFix_const {V : Type _} [AddCommGrp V] {dims : List Nat}
    (x0 : Fin dims.length) (a0 : Fin (dims.get x0)) (v : V) :
    sumTupleFix x0 a0 (fun _ => v) = nsmul (prodExcept x0) v := by
  induction dims with
  | nil => exact Fin.elim0 x0
  | cons d ds ih =>
    cases x0 with
    | mk idx h =>
      cases idx with
      | zero => exact sumTuple_const ds v
      | succ k =>
        simp only [sumTupleFix, prodExcept, ih, finSum_const, nsmul_mul]

/-- Summing `f(proj x t)` over all `t` equals `(prodExcept x) • (∑_a f a)`. -/
theorem sumTuple_proj {V : Type _} [AddCommGrp V] {dims : List Nat}
    (x : Fin dims.length) (f : Fin (dims.get x) → V) :
    sumTuple (fun t => f (proj x t)) = nsmul (prodExcept x) (finSum f) := by
  induction dims with
  | nil => exact Fin.elim0 x
  | cons d ds ih =>
    cases x with
    | mk idx h =>
      cases idx with
      | zero =>
        simp only [sumTuple, proj, prodExcept, sumTuple_const, finSum_nsmul]
      | succ k =>
        simp only [sumTuple, proj, prodExcept, ih, finSum_const, nsmul_mul]

/-- Summing `f(proj x t)` over all `t` with `proj x0 t = a0`:
    - if `x = x0`, equals `(prodExcept x0) • f a0`;
    - if `x ≠ x0`, multiplying by `dims.get x` equals `(prodExcept x0) • (∑_a f a)`! -/
theorem sumTupleFix_proj_eq {V : Type _} [AddCommGrp V] {dims : List Nat}
    (x0 : Fin dims.length) (a0 : Fin (dims.get x0)) (f : Fin (dims.get x0) → V) :
    sumTupleFix x0 a0 (fun t => f (proj x0 t)) = nsmul (prodExcept x0) (f a0) := by
  induction dims with
  | nil => exact Fin.elim0 x0
  | cons d ds ih =>
    cases x0 with
    | mk idx h =>
      cases idx with
      | zero => exact sumTuple_const ds (f a0)
      | succ k =>
        simp only [sumTupleFix, proj, prodExcept, ih, finSum_const, nsmul_mul]

theorem sumTupleFix_proj_ne_scaled {V : Type _} [AddCommGrp V] {dims : List Nat}
    (x0 x : Fin dims.length) (hne : x ≠ x0)
    (a0 : Fin (dims.get x0)) (f : Fin (dims.get x) → V) :
    nsmul (dims.get x) (sumTupleFix x0 a0 (fun t => f (proj x t))) =
    nsmul (prodExcept x0) (finSum f) := by
  induction dims with
  | nil => exact Fin.elim0 x0
  | cons d ds ih =>
    cases x0 with
    | mk i0 h0 =>
      cases x with
      | mk ix hx =>
        cases i0 with
        | zero =>
          cases ix with
          | zero => exact False.elim (hne rfl)
          | succ kx =>
            simp only [sumTupleFix, proj, prodExcept, List.get]
            rw [sumTuple_proj ⟨kx, Nat.lt_of_succ_lt_succ hx⟩ f, ← nsmul_mul,
                dims_mul_prodExcept ⟨kx, Nat.lt_of_succ_lt_succ hx⟩]
        | succ k0 =>
          cases ix with
          | zero =>
            simp only [sumTupleFix, proj, prodExcept, List.get, sumTupleFix_const, finSum_nsmul,
                       ← nsmul_mul]
          | succ kx =>
            have hne' : (⟨kx, Nat.lt_of_succ_lt_succ hx⟩ : Fin ds.length) ≠
                        ⟨k0, Nat.lt_of_succ_lt_succ h0⟩ := by
              intro e; apply hne
              cases e; rfl
            simp only [sumTupleFix, proj, prodExcept, List.get]
            rw [← finSum_nsmul, ih ⟨k0, _⟩ ⟨kx, _⟩ hne' a0 f, finSum_const, nsmul_mul]

theorem sumTupleFix_finSum {V : Type _} [AddCommGrp V] {dims : List Nat} {m : Nat}
    (x0 : Fin dims.length) (a0 : Fin (dims.get x0)) (F : Fin m → Tuple dims → V) :
    sumTupleFix x0 a0 (fun t => finSum (fun i => F i t)) =
    finSum (fun i => sumTupleFix x0 a0 (F i)) := by
  induction m with
  | zero =>
    simp only [finSum]
    rw [sumTupleFix_const, nsmul_zero]
  | succ k ih =>
    simp only [finSum]
    rw [sumTupleFix_add, ih]

/-- Helper: splitting a `finSum` at one index `x0 : Fin m` into `g x0` plus the sum of the
    remaining `m - 1` terms when all `x ≠ x0` equal a constant `c`. -/
theorem finSum_split_const {V : Type _} [AddCommGrp V] {m : Nat}
    (x0 : Fin m) (g : Fin m → V) (c : V)
    (hother : ∀ x : Fin m, x ≠ x0 → g x = c) :
    finSum g = add (g x0) (nsmul (m - 1) c) := by
  induction m with
  | zero => exact Fin.elim0 x0
  | succ k ih =>
    cases k with
    | zero =>
      cases x0 with
      | mk i hi =>
        have hi0 : i = 0 := by omega
        subst hi0
        simp [finSum, nsmul, zero_add, add_zero]
    | succ k' =>
      by_cases hlast : x0 = Fin.last (k' + 1)
      · subst hlast
        rw [finSum]
        have hconst : (fun i : Fin (k' + 1) => g i.castSucc) = (fun _ => c) := by
          funext i
          apply hother
          intro heq
          have hval := congrArg Fin.val heq
          simp [Fin.castSucc, Fin.last] at hval
          omega
        rw [hconst, finSum_const, add_comm]
        rfl
      · have hlt : x0.val < k' + 1 := by
          have := x0.isLt
          have hne : x0.val ≠ k' + 1 := fun e => hlast (Fin.ext e)
          omega
        let x0' : Fin (k' + 1) := ⟨x0.val, hlt⟩
        have hx0_eq : x0'.castSucc = x0 := Fin.ext rfl
        have hlast_eq : g (Fin.last (k' + 1)) = c := hother _ (fun e => hlast e.symm)
        rw [finSum]
        rw [hlast_eq]
        have ih_ap := ih x0' (fun i => g i.castSucc) (fun x hx => by
          apply hother
          intro e
          apply hx
          apply Fin.ext
          simpa [x0', Fin.castSucc] using congrArg Fin.val e)
        rw [ih_ap, hx0_eq, add_assoc]
        change add (g x0) (add (nsmul k' c) c) = add (g x0) (nsmul (k' + 1) c)
        rfl

/-- **Lemma 4.2, Step (a) (arbitrary number of settings `m ≥ 1` & arbitrary outcome counts `dims`):**
    Let `dims` be any non-empty list of outcome counts with product `N₁ = prodDims dims`.
    Given any consistent family `f x : Fin (dims.get x) → V` with common marginal
    `h = ∑_a f x a`, define the cleared-denominator lift on `Tuple dims`:
      `G(t) := ∑_x (dims.get x) • f x (proj x t) - (m - 1) • h`.
    Then for EVERY setting `x0` and EVERY outcome `a0`,
      `∑_{t : proj x0 t = a0} G(t) = N₁ • f x0 a0`! -/
theorem lem_lifting_step_a {V : Type _} [AddCommGrp V]
    (dims : List Nat) (f : (x : Fin dims.length) → Fin (dims.get x) → V) (h : V)
    (hcons : ∀ x : Fin dims.length, finSum (f x) = h)
    (x0 : Fin dims.length) (a0 : Fin (dims.get x0)) :
    sumTupleFix x0 a0
      (fun t => sub (finSum (fun x => nsmul (dims.get x) (f x (proj x t))))
                    (nsmul (dims.length - 1) h))
    = nsmul (prodDims dims) (f x0 a0) := by
  rw [sumTupleFix_sub, sumTupleFix_finSum, sumTupleFix_const]
  -- Evaluate the `finSum` term-by-term:
  -- For `x = x0`: `(dims.get x0) • (prodExcept x0) • f x0 a0 = N₁ • f x0 a0`
  -- For `x ≠ x0`: `(prodExcept x0) • h`
  have hterm : ∀ x : Fin dims.length, x ≠ x0 →
      sumTupleFix x0 a0 (fun t => nsmul (dims.get x) (f x (proj x t))) =
      nsmul (prodExcept x0) h := by
    intro x hne
    -- Pull `nsmul (dims.get x)` out of `sumTupleFix`:
    have hpull : ∀ {ds : List Nat} (c : Nat) (y0 : Fin ds.length) (b0 : Fin (ds.get y0)) (F : Tuple ds → V),
        sumTupleFix y0 b0 (fun t => nsmul c (F t)) = nsmul c (sumTupleFix y0 b0 F) := by
      intro ds
      induction ds with
      | nil => intro _ y0; exact Fin.elim0 y0
      | cons d ds' ih =>
        intro c y0 b0 F
        cases y0 with
        | mk idx hy =>
          cases idx with
          | zero =>
            simp only [sumTupleFix]
            have hst : ∀ {ds'' : List Nat} (H : Tuple ds'' → V),
                sumTuple (fun t => nsmul c (H t)) = nsmul c (sumTuple H) := by
              intro ds''
              induction ds'' with
              | nil => intro _; rfl
              | cons d'' ds''' ih'' => intro H; simp only [sumTuple, ih'', finSum_nsmul]
            exact hst _
          | succ k => simp only [sumTupleFix, ih, finSum_nsmul]
    rw [hpull, sumTupleFix_proj_ne_scaled x0 x hne a0 (f x), hcons x]
  rw [finSum_split_const x0 _ (nsmul (prodExcept x0) h) hterm]
  -- Evaluate the `x = x0` term:
  have hpull0 : sumTupleFix x0 a0 (fun t => nsmul (dims.get x0) (f x0 (proj x0 t))) =
                nsmul (dims.get x0) (sumTupleFix x0 a0 (fun t => f x0 (proj x0 t))) := by
    have hpull : ∀ {ds : List Nat} (c : Nat) (y0 : Fin ds.length) (b0 : Fin (ds.get y0)) (F : Tuple ds → V),
        sumTupleFix y0 b0 (fun t => nsmul c (F t)) = nsmul c (sumTupleFix y0 b0 F) := by
      intro ds
      induction ds with
      | nil => intro _ y0; exact Fin.elim0 y0
      | cons d ds' ih =>
        intro c y0 b0 F
        cases y0 with
        | mk idx hy =>
          cases idx with
          | zero =>
            simp only [sumTupleFix]
            have hst : ∀ {ds'' : List Nat} (H : Tuple ds'' → V),
                sumTuple (fun t => nsmul c (H t)) = nsmul c (sumTuple H) := by
              intro ds''
              induction ds'' with
              | nil => intro _; rfl
              | cons d'' ds''' ih'' => intro H; simp only [sumTuple, ih'', finSum_nsmul]
            exact hst _
          | succ k => simp only [sumTupleFix, ih, finSum_nsmul]
    exact hpull _ _ _ _
  rw [hpull0, sumTupleFix_proj_eq x0 a0 (f x0), ← nsmul_mul, dims_mul_prodExcept x0]
  -- Now cancel `nsmul (m - 1) (nsmul (prodExcept x0) h)` with `nsmul (prodExcept x0) (nsmul (m - 1) h)`!
  rw [← nsmul_mul (dims.length - 1) (prodExcept x0) h,
      ← nsmul_mul (prodExcept x0) (dims.length - 1) h,
      Nat.mul_comm (prodExcept x0) (dims.length - 1),
      add_sub_cancel]

end Theorem101.BipartiteLift
