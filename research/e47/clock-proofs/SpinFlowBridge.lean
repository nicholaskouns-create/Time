import Lean
theorem spin_lock (s d z : Nat) (h : 2*z=d+z) (v : d+z=10)
    (c : d=2*s+1) : s=2 ∧ d=5 ∧ z=5 := by omega
theorem walk_identity (q d z : Nat) (h : q*z=d+z) :
    q*(2*z)=2*(d+z) := by
  rw [Nat.mul_left_comm, h]
def sizes := [1,3,6,10,15,18,19,18,15,10,6,3,1]
example : (sizes.map (fun n => n^3)).sum = 27761 := by decide
#print axioms spin_lock
