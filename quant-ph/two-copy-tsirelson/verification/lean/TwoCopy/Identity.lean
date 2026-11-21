/-
Two-copy identity, file 5 of 6: Theorem 1 (target (c)).

What this file proves, in plain words. P is any real function on the 256 events, D any real
bilinear form on the functions on Ω, with
  (i)  D(1_e, 1_e) = P(e) for every event e,
  (ii) D(1_e, 1_e′) = 0 for every locally orthogonal pair e, e′.
Nothing else is assumed in `theorem1_general`: D need not be symmetric, P need not be
normalised, non-negative or no-signalling, D need not be invariant under the τ_j.

  * `theorem1_general`: Σ_k D(1_{S_k}, 1_{S_k}) = 2·D(1, 1) − ⟨C₁C₂⟩_P / 4.
  * `D_one_one_eq`: D(1, 1) = Σ_o P(s, o) for every choice of boxes s. So the constant "2" of
    the draft is 2·D(1, 1), and it is 2 as soon as P sums to 1 over the outcomes of ONE choice
    of boxes (and then, by (i) and (ii), of every choice).
  * `theorem1`: with Σ_o P(s₀, o) = 1 for one s₀: Σ_k D(1_{S_k}, 1_{S_k}) = 2 − ⟨C₁C₂⟩_P / 4.
  * `theorem1_compl`: the same for the four complements, with the sign of the correlator
    reversed: Σ_k D(1 − 1_{S_k}, 1 − 1_{S_k}) = 2·D(1, 1) + ⟨C₁C₂⟩_P / 4.
  * `indF_Sset_mem_W`: each 1_{S_k} lies in W, the span of the event indicators; and
    `theorem1_W`: Theorem 1 for a form defined on W only (the draft's setting).

Route. Not the draft's Steps 4 and 5 (invariant D, then averaging over H). By Step 1, for two
monomials m_s, m_s′ the value D(m_s, m_s′) depends only on the pattern d of parties where s and
s′ differ and on s restricted to d (`D_mono`); so Σ_k D(M_k, M_k), M_k = M ∘ rep k, is an
integer combination of the unknowns R(d, t), and `Finite.coeff` says every coefficient vanishes
except that of R(0, 0) = D(1, 1), which is 64 (`sum_D_MkR`). The linear terms give the
correlator through `Finite.sum_epsk`.
-/
import TwoCopy.Step1

namespace TwoCopy

open Finset

/-! ### Partial monomials and the reduction of D(m_s, m_s′) -/

/-- The partial monomial: the product of the variables var(j, s j) over the parties j in A. -/
def pm (A : Finset (Fin 4)) (s : Boxes) : Pt → ℝ := ∏ j ∈ A, varR j (s j)

theorem pm_apply (A : Finset (Fin 4)) (s : Boxes) (ω : Pt) :
    pm A s ω = ∏ j ∈ A, (sgn (ω j (s j)) : ℝ) := by
  simp [pm, Finset.prod_apply, varR, var]

theorem pm_univ (s : Boxes) : pm univ s = monoR s := by
  funext ω; simp [pm_apply, monoR, mono, var]

theorem pm_congr (A : Finset (Fin 4)) (s s' : Boxes) (h : ∀ j ∈ A, s j = s' j) :
    pm A s = pm A s' :=
  Finset.prod_congr rfl fun j hj => by rw [h j hj]

theorem pm_mem_Wexc (A : Finset (Fin 4)) (s : Boxes) (j : Fin 4) (hj : j ∉ A) :
    pm A s ∈ Wexc j := by
  have h : pm A s =
      fun ω : Pt => (fun o : Boxes => ∏ i ∈ A, (sgn (o i) : ℝ)) (fun i => ω i (s i)) := by
    funext ω; exact pm_apply A s ω
  rw [h]
  refine boxfun_mem_Wexc j s (fun o : Boxes => ∏ i ∈ A, (sgn (o i) : ℝ)) ?_
  intro o a
  refine Finset.prod_congr rfl fun i hi => ?_
  have hij : i ≠ j := by rintro rfl; exact hj hi
  rw [Function.update_of_ne hij]

/-- A party with the same box in s and s′ drops out (Step 1). -/
theorem pm_erase (D : Form) (hii : Decoheres D) (A : Finset (Fin 4)) (s s' : Boxes) (j : Fin 4)
    (hj : j ∈ A) (h : s j = s' j) :
    D (pm A s) (pm A s') = D (pm (A.erase j) s) (pm (A.erase j) s') := by
  have e1 : pm A s = varR j (s j) * pm (A.erase j) s :=
    (Finset.mul_prod_erase A (fun i => varR i (s i)) hj).symm
  have e2 : pm A s' = varR j (s j) * pm (A.erase j) s' := by
    rw [h]; exact (Finset.mul_prod_erase A (fun i => varR i (s' i)) hj).symm
  rw [e1, e2]
  exact step1 D hii j (s j) (pm_mem_Wexc _ _ _ (Finset.notMem_erase j A))
    (pm_mem_Wexc _ _ _ (Finset.notMem_erase j A))

theorem pm_reduce (D : Form) (hii : Decoheres D) (s s' : Boxes) (B : Finset (Fin 4))
    (hB : ∀ j ∈ B, s j = s' j) :
    D (pm univ s) (pm univ s') = D (pm (univ \ B) s) (pm (univ \ B) s') := by
  induction B using Finset.induction_on with
  | empty => simp
  | insert j B hj ih =>
    rw [ih (fun i hi => hB i (Finset.mem_insert_of_mem hi)),
      pm_erase D hii (univ \ B) s s' j (by simp [hj]) (hB j (Finset.mem_insert_self j B)),
      Finset.sdiff_insert]

/-- The parties where the pattern d is 1. -/
def supp (d : Boxes) : Finset (Fin 4) := univ.filter fun j => d j = 1

/-- The unknowns: for a difference pattern d and boxes t, the value of D on the two partial
monomials over the parties of d, with boxes t and the opposite boxes. -/
def R (D : Form) (d t : Boxes) : ℝ := D (pm (supp d) t) (pm (supp d) (t + d))

theorem R_zero (D : Form) : R D 0 0 = D 1 1 := by
  have : supp 0 = ∅ := by ext j; simp [supp]
  simp [R, this, pm]

/-- D(m_s, m_{s+d}) depends only on d and on s restricted to d. -/
theorem D_mono (D : Form) (hii : Decoheres D) (s d : Boxes) :
    D (monoR s) (monoR (s + d)) = R D d (s * d) := by
  have h0 : ∀ a : Fin 2, a ≠ 1 → a = 0 := by decide
  have hB : ∀ j ∈ univ \ supp d, s j = (s + d) j := by
    intro j hj
    have : d j = 0 := h0 _ (by simpa [supp] using hj)
    simp [this]
  have hs : univ \ (univ \ supp d) = supp d := by ext j; simp
  rw [← pm_univ, ← pm_univ, pm_reduce D hii s (s + d) _ hB, hs]
  unfold R
  have h1 : ∀ j ∈ supp d, d j = 1 := fun j hj => by simpa [supp] using hj
  rw [pm_congr (supp d) s (s * d) (fun j hj => by simp [h1 j hj]),
    pm_congr (supp d) (s + d) (s * d + d) (fun j hj => by simp [h1 j hj])]

/-! ### The functions M_k = M ∘ rep k -/

/-- M ∘ rep k as a real function, written as its signed sum of monomials. -/
def MkR (k : Fin 4) : Pt → ℝ := ∑ s : Boxes, (epsk k s : ℝ) • monoR s

theorem MkR_apply (k : Fin 4) (ω : Pt) : MkR k ω = (M (rep k ω) : ℝ) := by
  rw [fourier]; simp [MkR, Finset.sum_apply, monoR]

/-- 1_{S_k} = 1/2 − M_k/8, as functions. -/
theorem indF_Sset_eq (k : Fin 4) :
    indF (Sset k) = (1 / 2 : ℝ) • (1 : Pt → ℝ) - (1 / 8 : ℝ) • MkR k := by
  funext ω
  simp only [Pi.sub_apply, Pi.smul_apply, Pi.one_apply, smul_eq_mul, indF_Sset, MkR_apply]
  ring

theorem D_MkR (D : Form) (k : Fin 4) :
    D (MkR k) (MkR k) =
      ∑ s : Boxes, ∑ s' : Boxes, ((epsk k s * epsk k s' : ℤ) : ℝ) * D (monoR s) (monoR s') := by
  simp only [MkR, map_sum, map_smul, LinearMap.sum_apply, LinearMap.smul_apply, smul_eq_mul,
    Finset.mul_sum]
  rw [Finset.sum_comm]
  exact Finset.sum_congr rfl fun s _ => Finset.sum_congr rfl fun s' _ => by push_cast; ring

/-- **The quadratic terms.** Σ_k D(M_k, M_k) = 64·D(1, 1), from (ii) alone. -/
theorem sum_D_MkR (D : Form) (hii : Decoheres D) : ∑ k, D (MkR k) (MkR k) = 64 * D 1 1 := by
  have hc : ∀ d t : Boxes,
      (∑ s : Boxes with s * d = t, ∑ k, ((epsk k s * epsk k (s + d) : ℤ) : ℝ)) =
        if d = 0 ∧ t = 0 then 64 else 0 := by
    intro d t
    exact_mod_cast congrArg (fun z : ℤ => (z : ℝ)) (coeff d t)
  calc ∑ k, D (MkR k) (MkR k)
      = ∑ k, ∑ s : Boxes, ∑ s' : Boxes,
          ((epsk k s * epsk k s' : ℤ) : ℝ) * D (monoR s) (monoR s') :=
        Finset.sum_congr rfl fun k _ => D_MkR D k
    _ = ∑ s : Boxes, ∑ s' : Boxes, ∑ k,
          ((epsk k s * epsk k s' : ℤ) : ℝ) * D (monoR s) (monoR s') := by
        rw [Finset.sum_comm]; exact Finset.sum_congr rfl fun s _ => Finset.sum_comm
    _ = ∑ s : Boxes, ∑ d : Boxes, ∑ k,
          ((epsk k s * epsk k (s + d) : ℤ) : ℝ) * R D d (s * d) := by
        refine Finset.sum_congr rfl fun s _ => ?_
        rw [← Equiv.sum_comp (Equiv.addLeft s)
          (fun s' => ∑ k, ((epsk k s * epsk k s' : ℤ) : ℝ) * D (monoR s) (monoR s'))]
        refine Finset.sum_congr rfl fun d _ => Finset.sum_congr rfl fun k _ => ?_
        simp only [Equiv.coe_addLeft, D_mono D hii]
    _ = ∑ d : Boxes, ∑ s : Boxes,
          (∑ k, ((epsk k s * epsk k (s + d) : ℤ) : ℝ)) * R D d (s * d) := by
        rw [Finset.sum_comm]
        exact Finset.sum_congr rfl fun d _ => Finset.sum_congr rfl fun s _ =>
          (Finset.sum_mul _ _ _).symm
    _ = ∑ d : Boxes, ∑ t : Boxes,
          (∑ s : Boxes with s * d = t, ∑ k, ((epsk k s * epsk k (s + d) : ℤ) : ℝ)) * R D d t := by
        refine Finset.sum_congr rfl fun d _ => ?_
        rw [← Finset.sum_fiberwise univ (fun s : Boxes => s * d)]
        refine Finset.sum_congr rfl fun t _ => ?_
        rw [Finset.sum_mul]
        refine Finset.sum_congr rfl fun s hs => ?_
        rw [(Finset.mem_filter.1 hs).2]
    _ = ∑ d : Boxes, ∑ t : Boxes, (if d = 0 ∧ t = 0 then (64 : ℝ) else 0) * R D d t := by
        simp only [hc]
    _ = 64 * D 1 1 := by
        rw [Finset.sum_eq_single 0, Finset.sum_eq_single 0]
        · simp [R_zero]
        · intro t _ ht; simp [ht]
        · simp
        · intro d _ hd; simp [hd]
        · simp

/-! ### The linear terms -/

theorem monoR_eq (s : Boxes) :
    monoR s = ∑ o : Boxes,
      ((-1 : ℝ) ^ ((o 0).val + (o 1).val + (o 2).val + (o 3).val)) • ind (s, o) := by
  rw [← boxfun_eq s (fun o => (-1 : ℝ) ^ ((o 0).val + (o 1).val + (o 2).val + (o 3).val))]
  funext ω
  simp [monoR, mono, var, sgn, Fin.prod_univ_four, pow_add]

theorem D_one_mono (P : Event → ℝ) (D : Form) (hi : ∀ e, D (ind e) (ind e) = P e)
    (hii : Decoheres D) (s : Boxes) : D 1 (monoR s) = E P s := by
  rw [monoR_eq, map_sum]; unfold E
  refine Finset.sum_congr rfl fun o _ => ?_
  rw [map_smul, smul_eq_mul, D_one_ind D hii, hi]

theorem D_mono_one (P : Event → ℝ) (D : Form) (hi : ∀ e, D (ind e) (ind e) = P e)
    (hii : Decoheres D) (s : Boxes) : D (monoR s) 1 = E P s := by
  rw [monoR_eq, map_sum, LinearMap.sum_apply]; unfold E
  refine Finset.sum_congr rfl fun o _ => ?_
  rw [map_smul, LinearMap.smul_apply, smul_eq_mul, D_ind_one D hii, hi]

theorem sum_epsk_E (P : Event → ℝ) : ∑ k, ∑ s : Boxes, (epsk k s : ℝ) * E P s = 2 * corr P := by
  rw [Finset.sum_comm]; unfold corr; rw [Finset.mul_sum]
  refine Finset.sum_congr rfl fun s _ => ?_
  rw [← Finset.sum_mul]
  have h2 : (∑ k, (epsk k s : ℝ)) =
      2 * (-1) ^ ((s 0).val * (s 1).val + (s 2).val * (s 3).val) := by
    exact_mod_cast congrArg (fun z : ℤ => (z : ℝ)) (sum_epsk s)
  rw [h2]; ring

theorem sum_D_one_MkR (P : Event → ℝ) (D : Form) (hi : ∀ e, D (ind e) (ind e) = P e)
    (hii : Decoheres D) : ∑ k, D 1 (MkR k) = 2 * corr P := by
  simp only [MkR, map_sum, map_smul, smul_eq_mul, D_one_mono P D hi hii]
  exact sum_epsk_E P

theorem sum_D_MkR_one (P : Event → ℝ) (D : Form) (hi : ∀ e, D (ind e) (ind e) = P e)
    (hii : Decoheres D) : ∑ k, D (MkR k) 1 = 2 * corr P := by
  simp only [MkR, map_sum, map_smul, LinearMap.sum_apply, LinearMap.smul_apply, smul_eq_mul,
    D_mono_one P D hi hii]
  exact sum_epsk_E P

/-! ### Theorem 1 -/

/-- D(1, 1) is the total of P over the outcomes of any one choice of boxes. -/
theorem D_one_one_eq (P : Event → ℝ) (D : Form) (hi : ∀ e, D (ind e) (ind e) = P e)
    (hii : Decoheres D) (s : Boxes) : D 1 1 = ∑ o : Boxes, P (s, o) := by
  rw [D_one_one D hii s]; exact Finset.sum_congr rfl fun o _ => hi _

theorem D_indS (D : Form) (k : Fin 4) :
    D (indF (Sset k)) (indF (Sset k)) =
      D 1 1 / 4 - (D 1 (MkR k) + D (MkR k) 1) / 16 + D (MkR k) (MkR k) / 64 := by
  rw [indF_Sset_eq]
  simp only [map_sub, map_smul, LinearMap.sub_apply, LinearMap.smul_apply, smul_eq_mul]
  ring

/-- **Theorem 1, general form.** For every P and every bilinear D with (i) and (ii), and
nothing else: Σ_k D(1_{S_k}, 1_{S_k}) = 2·D(1, 1) − ⟨C₁C₂⟩_P / 4. -/
theorem theorem1_general (P : Event → ℝ) (D : Form) (hi : ∀ e, D (ind e) (ind e) = P e)
    (hii : ∀ e e', LO e e' → D (ind e) (ind e') = 0) :
    ∑ k, D (indF (Sset k)) (indF (Sset k)) = 2 * D 1 1 - corr P / 4 := by
  have h1 := sum_D_MkR D hii
  have h2 := sum_D_one_MkR P D hi hii
  have h3 := sum_D_MkR_one P D hi hii
  simp only [Fin.sum_univ_four] at h1 h2 h3 ⊢
  rw [D_indS, D_indS, D_indS, D_indS]
  linarith

/-- **Theorem 1.** For every P : events → ℝ that sums to 1 over the outcomes of one choice of
boxes s₀, and every bilinear D with (i) D(1_e, 1_e) = P(e) and (ii) D(1_e, 1_e′) = 0 for
locally orthogonal e, e′: Σ_k D(1_{S_k}, 1_{S_k}) = 2 − ⟨C₁C₂⟩_P / 4. -/
theorem theorem1 (P : Event → ℝ) (D : Form) (hi : ∀ e, D (ind e) (ind e) = P e)
    (hii : ∀ e e', LO e e' → D (ind e) (ind e') = 0)
    (s₀ : Boxes) (hnorm : ∑ o : Boxes, P (s₀, o) = 1) :
    ∑ k, D (indF (Sset k)) (indF (Sset k)) = 2 - corr P / 4 := by
  rw [theorem1_general P D hi hii, D_one_one_eq P D hi hii s₀, hnorm, mul_one]

/-- The indicator of the complement of S_k is 1 − 1_{S_k}. -/
theorem indF_compl (T : Finset Pt) : indF Tᶜ = 1 - indF T := by
  funext ω
  by_cases h : ω ∈ T <;> simp [indF, h]

/-- **Theorem 1 for the complements.** Σ_k D(1_{S_kᶜ}, 1_{S_kᶜ}) = 2·D(1, 1) + ⟨C₁C₂⟩_P / 4. -/
theorem theorem1_compl (P : Event → ℝ) (D : Form) (hi : ∀ e, D (ind e) (ind e) = P e)
    (hii : ∀ e e', LO e e' → D (ind e) (ind e') = 0) :
    ∑ k, D (indF (Sset k)ᶜ) (indF (Sset k)ᶜ) = 2 * D 1 1 + corr P / 4 := by
  have h1 := sum_D_MkR D hii
  have h2 := sum_D_one_MkR P D hi hii
  have h3 := sum_D_MkR_one P D hi hii
  have hk : ∀ k, D (indF (Sset k)ᶜ) (indF (Sset k)ᶜ) =
      D 1 1 / 4 + (D 1 (MkR k) + D (MkR k) 1) / 16 + D (MkR k) (MkR k) / 64 := by
    intro k
    rw [indF_compl, indF_Sset_eq]
    simp only [map_sub, map_smul, LinearMap.sub_apply, LinearMap.smul_apply, smul_eq_mul]
    ring
  simp only [Fin.sum_univ_four] at h1 h2 h3 ⊢
  rw [hk, hk, hk, hk]
  linarith

/-! ### The draft's setting: forms on W -/

theorem ind_mem_W (e : Event) : ind e ∈ W := Submodule.subset_span ⟨e, rfl⟩

theorem one_mem_W : (1 : Pt → ℝ) ∈ W := by
  rw [← sum_ind 0]; exact Submodule.sum_mem _ fun o _ => ind_mem_W _

theorem monoR_mem_W (s : Boxes) : monoR s ∈ W := by
  rw [monoR_eq]; exact Submodule.sum_mem _ fun o _ => Submodule.smul_mem _ _ (ind_mem_W _)

/-- Each 1_{S_k} lies in W, the span of the event indicators. -/
theorem indF_Sset_mem_W (k : Fin 4) : indF (Sset k) ∈ W := by
  rw [indF_Sset_eq]
  refine Submodule.sub_mem _ (Submodule.smul_mem _ _ one_mem_W) (Submodule.smul_mem _ _ ?_)
  exact Submodule.sum_mem _ fun s _ => Submodule.smul_mem _ _ (monoR_mem_W s)

/-- **Theorem 1 for a form on W only**, as in the draft. -/
theorem theorem1_W (P : Event → ℝ) (D : W →ₗ[ℝ] W →ₗ[ℝ] ℝ)
    (hi : ∀ e, D ⟨ind e, ind_mem_W e⟩ ⟨ind e, ind_mem_W e⟩ = P e)
    (hii : ∀ e e', LO e e' → D ⟨ind e, ind_mem_W e⟩ ⟨ind e', ind_mem_W e'⟩ = 0)
    (s₀ : Boxes) (hnorm : ∑ o : Boxes, P (s₀, o) = 1) :
    ∑ k, D ⟨indF (Sset k), indF_Sset_mem_W k⟩ ⟨indF (Sset k), indF_Sset_mem_W k⟩ =
      2 - corr P / 4 := by
  obtain ⟨p, hp⟩ := LinearMap.exists_leftInverse_of_injective W.subtype (Submodule.ker_subtype W)
  have hp' : ∀ (F : Pt → ℝ) (hF : F ∈ W), p F = ⟨F, hF⟩ := fun F hF =>
    LinearMap.congr_fun hp ⟨F, hF⟩
  have h := theorem1 P (D.compl₁₂ p p)
    (fun e => by rw [LinearMap.compl₁₂_apply, hp' _ (ind_mem_W e)]; exact hi e)
    (fun e e' he => by
      rw [LinearMap.compl₁₂_apply, hp' _ (ind_mem_W e), hp' _ (ind_mem_W e')]; exact hii e e' he)
    s₀ hnorm
  rw [← h]
  refine Finset.sum_congr rfl fun k _ => ?_
  rw [LinearMap.compl₁₂_apply, hp' _ (indF_Sset_mem_W k)]

end TwoCopy
