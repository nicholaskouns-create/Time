import Lean
set_option maxRecDepth 20000
set_option maxHeartbeats 0
abbrev T := Nat × Nat × Nat
def ts : List T := (List.range 5).flatMap fun a => (List.range 5).flatMap fun b => (List.range 5).map fun c => (a,b,c)
def entries (t : T) := [t.1,t.2.1,t.2.2]
def basis (w : Nat) := ts.filter fun t => (entries t).sum == w
def cas (w : Nat) (u v : T) : Int := Id.run do
  let m : Int := (w:Int)-6
  let mut a := if u == v then m*(m+1) else 0
  let l := entries v
  for i in List.range 3 do
    let p := l[i]!
    if p < 4 then
      let r := l.set i (p+1)
      for j in List.range 3 do
        let q := r[j]!
        if q > 0 && r.set j (q-1) == entries u then
          a := a + ((4-p)*q : Nat)
  return a
def diagData (w : Nat) := ((basis w).length, ((basis w).map fun t => cas w t t).sum)
example : (List.range 13).map diagData =
[(1, 42), (3, 102), (6, 162), (10, 210), (15, 240), (18, 246), (19, 246), (18, 246), (15, 240), (10, 210), (6, 162), (3, 102), (1, 42)] := by decide
