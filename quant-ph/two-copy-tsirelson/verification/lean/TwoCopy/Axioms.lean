/-
Two-copy identity: the axiom audit.

What this file does, in plain words. It proves nothing new. It imports the six files and
prints, for every final theorem, the axioms its proof depends on. Building it
(`lake build TwoCopy.Axioms`) shows the list. The expected output for each line is
`[propext, Classical.choice, Quot.sound]` (Lean's three standard axioms) or a subset: no
`sorryAx`, and no `Lean.ofReduceBool` / `Lean.trustCompiler` (no `native_decide` is used;
the closed computations are `decide +kernel`, checked by the kernel).
-/
import TwoCopy.Bounds

namespace TwoCopy

-- target (a)
#print axioms M_pm
#print axioms card_S
#print axioms image_mem_H
#print axioms Sset_is_image
#print axioms Sset_injective
#print axioms images_eq
#print axioms card_images
#print axioms pointwise
#print axioms pointwise_const_unique
#print axioms single_set_false
-- the closed computations Theorem 1 stands on
#print axioms fourier
#print axioms sum_epsk
#print axioms coeff
-- target (b)
#print axioms step1
#print axioms D_one_ind
#print axioms D_ind_one
-- target (c)
#print axioms theorem1_general
#print axioms D_one_one_eq
#print axioms theorem1
#print axioms theorem1_compl
#print axioms indF_Sset_mem_W
#print axioms theorem1_W
-- target (d)
#print axioms corollary2
#print axioms corollary2_upper
#print axioms corollary2_lower
#print axioms corollary2_lower_by_relabeling
-- target (e)
#print axioms corr_PRq
#print axioms PRq_norm
#print axioms theorem3_only_if
#print axioms theorem3_only_if_sqrt
#print axioms theorem3_JQM
-- tests
#print axioms card_event_points
#print axioms card_LO_row
#print axioms var_tau
#print axioms corr_pointmass
#print axioms evalForm_hyps

end TwoCopy
