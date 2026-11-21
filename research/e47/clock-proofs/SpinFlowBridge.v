From Coq Require Import Arith Lia List.
Import ListNotations.
Theorem spin_lock (s d z:nat) :
  2*z=d+z -> d+z=10 -> d=2*s+1 -> s=2 /\ d=5 /\ z=5.
Proof. lia. Qed.
Theorem walk_identity (q d z:nat) :
  q*z=d+z -> q*(2*z)=2*(d+z).
Proof. nia. Qed.
Print Assumptions spin_lock.
