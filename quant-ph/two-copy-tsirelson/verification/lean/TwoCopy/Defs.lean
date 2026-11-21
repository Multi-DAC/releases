/-
Two-copy identity, file 1 of 6: the objects.

What this file does, in plain words. It defines, and proves nothing about, the objects of the
draft "A joint quantum measure on two copies gives Tsirelson's bound" (Sections 1 and 2):

  * the sample space Ω (`Pt`): an outcome in {0, 1} for each of the 4 parties' 2 boxes;
    parties 0, 1, 2, 3 are A₁, B₁, A₂, B₂, the partner of j is j xor 1;
  * the box variables `var j b`, with values ±1 in ℤ;
  * events (`Event`: a box and an outcome for each party), their indicators `ind`, and local
    orthogonality `LO`;
  * M (`M`), the set S = [M = −4] (`S`), the relabelings τ_j (`tau`, `tauPerm`), the group H
    they generate (`H`), the four named images `Sset k` = S, τ₀S, τ₁S, τ₃S;
  * the CHSH expressions C₁, C₂ (`C 0`, `C 1`);
  * for a function P on events, E_P (`E`) and ⟨C₁C₂⟩_P (`corr`);
  * the isotropic box p_q (`pq`) and the product of two of them (`PRq`).

Integer-valued objects are ℤ-valued so that the finite facts of `Finite.lean` are closed
computations; the real functions on Ω that a bilinear form eats are their casts (`varR`, `monoR`,
`ind`, `indF`). Every definition here is tested by a fact in `Finite.lean` that would be false if
the definition were wrong (counts, values at named points).
-/
import Mathlib

namespace TwoCopy

/-- A point of Ω: an outcome in {0, 1} for each (party, box). 256 points. -/
abbrev Pt := Fin 4 → Fin 2 → Fin 2

/-- One element of {0, 1} per party: a choice of boxes, or a tuple of outcomes. -/
abbrev Boxes := Fin 4 → Fin 2

/-- An event: a box (first component) and an outcome (second component) for each party. -/
abbrev Event := Boxes × Boxes

/-- The partner of party j is j xor 1: A₁ ↔ B₁, A₂ ↔ B₂. -/
def partner (j : Fin 4) : Fin 4 := ![1, 0, 3, 2] j

/-- (−1)^a for a ∈ {0, 1}. -/
def sgn (a : Fin 2) : ℤ := (-1) ^ a.val

/-- The box variable: var(j, b)(ω) = (−1)^ω(j, b). -/
def var (j : Fin 4) (b : Fin 2) (ω : Pt) : ℤ := sgn (ω j b)

/-- The event e happens at ω: every party's chosen box shows the event's outcome. -/
def Holds (e : Event) (ω : Pt) : Prop := ∀ j, ω j (e.1 j) = e.2 j

instance (e : Event) : DecidablePred (Holds e) := fun _ => by unfold Holds; infer_instance

/-- Local orthogonality: some party has the same box in both events and different outcomes. -/
def LO (e e' : Event) : Prop := ∃ j, e.1 j = e'.1 j ∧ e.2 j ≠ e'.2 j

instance : DecidableRel LO := fun _ _ => by unfold LO; infer_instance

/-- M = Σ_{x,y,x′,y′} (−1)^(x + x x′ + x y′ + y y′) α_x β_y α′_x′ β′_y′. -/
def M (ω : Pt) : ℤ :=
  ∑ x : Fin 2, ∑ y : Fin 2, ∑ x' : Fin 2, ∑ y' : Fin 2,
    (-1) ^ (x.val + x.val * x'.val + x.val * y'.val + y.val * y'.val)
      * (var 0 x ω * var 1 y ω * var 2 x' ω * var 3 y' ω)

/-- S = [M = −4]. -/
def S : Finset Pt := Finset.univ.filter fun ω => M ω = -4

/-- τ_j: swap party j's two boxes, and flip the outcome of box 1 of j's partner. -/
def tau (j : Fin 4) (ω : Pt) : Pt := fun i b =>
  if i = j then ω i (b + 1)
  else if i = partner j ∧ b = 1 then ω i b + 1
  else ω i b

/-- The four relabelings that name the images of S: the identity, τ₀, τ₁, τ₃. -/
def rep (k : Fin 4) : Pt → Pt := ![id, tau 0, tau 1, tau 3] k

/-- The four named images of S: S₀ = S, S₁ = τ₀ S, S₂ = τ₁ S, S₃ = τ₃ S. -/
def Sset (k : Fin 4) : Finset Pt := S.image (rep k)

/-- The CHSH expression of copy c: C = α₀β₀ + α₀β₁ + α₁β₀ − α₁β₁, on parties 2c, 2c + 1. -/
def C (c : Fin 2) (ω : Pt) : ℤ :=
  let a : Fin 4 := ![0, 2] c
  let b : Fin 4 := ![1, 3] c
  var a 0 ω * var b 0 ω + var a 0 ω * var b 1 ω + var a 1 ω * var b 0 ω - var a 1 ω * var b 1 ω

/-- The monomial m_s = α_x β_y α′_x′ β′_y′ for the box choice s = (x, y, x′, y′). -/
def mono (s : Boxes) (ω : Pt) : ℤ := ∏ j, var j (s j) ω

/-- The sign of m_s in M: (−1)^(x + x x′ + x y′ + y y′). -/
def eps (s : Boxes) : ℤ :=
  (-1) ^ ((s 0).val + (s 0).val * (s 2).val + (s 0).val * (s 3).val + (s 1).val * (s 3).val)

/-- The signs of the monomials in M ∘ rep k (proved in `Finite.fourier`). -/
def epsk (k : Fin 4) (s : Boxes) : ℤ :=
  ![eps s, -(eps s * sgn (s 1) * sgn (s 2) * sgn (s 3)), eps s * sgn (s 0) * sgn (s 3),
    eps s * sgn (s 0) * sgn (s 1) * sgn (s 2)] k

/-! ### Real functions on Ω -/

/-- The box variable as a real function on Ω. -/
def varR (j : Fin 4) (b : Fin 2) : Pt → ℝ := fun ω => (var j b ω : ℝ)

/-- The monomial m_s as a real function on Ω. -/
def monoR (s : Boxes) : Pt → ℝ := fun ω => (mono s ω : ℝ)

/-- The indicator of an event. -/
def ind (e : Event) : Pt → ℝ := fun ω => if Holds e ω then 1 else 0

/-- The indicator of a subset T of Ω. -/
def indF (T : Finset Pt) : Pt → ℝ := fun ω => if ω ∈ T then 1 else 0

/-- W, the span of the event indicators. -/
def W : Submodule ℝ (Pt → ℝ) := Submodule.span ℝ (Set.range ind)

/-- The indicator of "every party other than j shows the event's outcome in the event's box":
an event of the other three parties. -/
def indExcept (j : Fin 4) (e : Event) : Pt → ℝ :=
  fun ω => if ∀ i, i ≠ j → ω i (e.1 i) = e.2 i then 1 else 0

/-- W_{−j}, the span of the indicators of events of the three parties other than j. -/
def Wexc (j : Fin 4) : Submodule ℝ (Pt → ℝ) := Submodule.span ℝ (Set.range (indExcept j))

/-! ### The relabeling group -/

theorem tau_invol : ∀ j ω, tau j (tau j ω) = ω := by decide +kernel

/-- τ_j as a permutation of Ω. -/
def tauPerm (j : Fin 4) : Equiv.Perm Pt := ⟨tau j, tau j, tau_invol j, tau_invol j⟩

/-- H, the group the four τ_j generate. -/
def H : Subgroup (Equiv.Perm Pt) := Subgroup.closure (Set.range tauPerm)

/-! ### Behaviours and correlators -/

/-- E_P(x, y, x′, y′) = Σ over outcomes o of (−1)^(o₀+o₁+o₂+o₃) P(event with boxes s, outcomes o). -/
def E (P : Event → ℝ) (s : Boxes) : ℝ :=
  ∑ o : Boxes, (-1 : ℝ) ^ ((o 0).val + (o 1).val + (o 2).val + (o 3).val) * P (s, o)

/-- ⟨C₁C₂⟩_P = Σ_{x,y,x′,y′} (−1)^(x y + x′ y′) E_P(x, y, x′, y′). -/
def corr (P : Event → ℝ) : ℝ :=
  ∑ s : Boxes, (-1 : ℝ) ^ ((s 0).val * (s 1).val + (s 2).val * (s 3).val) * E P s

/-- The isotropic box: p_q(x, y; a, b) = (1 + q)/4 if a xor b = x·y, and (1 − q)/4 otherwise. -/
noncomputable def pq (q : ℝ) (x y a b : Fin 2) : ℝ := if a + b = x * y then (1 + q) / 4 else (1 - q) / 4

/-- The product of two isotropic boxes: P(e) = p_q(copy-1 part of e) · p_q(copy-2 part of e). -/
noncomputable def PRq (q : ℝ) (e : Event) : ℝ :=
  pq q (e.1 0) (e.1 1) (e.2 0) (e.2 1) * pq q (e.1 2) (e.1 3) (e.2 2) (e.2 3)

end TwoCopy
