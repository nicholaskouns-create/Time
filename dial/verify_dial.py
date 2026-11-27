#!/usr/bin/env python3
"""A Dial of Forty-Seven Hours: independent NumPy replay of the proof's finite claims.

Builds V2 (x) V2 (x) V2 from the spin-2 ladder, then checks the grading S, the index of every
pairing, the symmetry obstruction, the symmetric dial's seven hours (commutant M5 + M2) and the
two clocks. Writes dial_data.json with --write; by default checks the committed file.
Usage: python dial/verify_dial.py [--write]      Exit 0 = every check passes.
"""
import json, sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
TOL = 1e-9
rng = np.random.default_rng(47)
results, failures = [], []


def check(cid, title, ok, detail=""):
    results.append({"id": cid, "title": title, "pass": bool(ok), "detail": detail})
    (failures.append(cid) if not ok else None)
    print(("PASS " if ok else "FAIL ") + cid + "  " + title + (("  [" + detail + "]") if detail else ""))


# ---------------------------------------------------------------- carrier
m = np.array([2, 1, 0, -1, -2.0])
Jp1 = np.zeros((5, 5))
for k in range(1, 5):
    Jp1[k - 1, k] = np.sqrt(6 - m[k] * (m[k] + 1))
single = [(Jp1 + Jp1.T) / 2, (Jp1 - Jp1.T) / 2j, np.diag(m).astype(complex)]
I5, I = np.eye(5), np.eye(125)


def site(A, k):
    f = [I5, I5, I5]; f[k] = A
    return np.kron(np.kron(f[0], f[1]), f[2])


Js = [[site(A, k) for A in single] for k in range(3)]
Jt = [sum(Js[k][a] for k in range(3)) for a in range(3)]
C = sum(J @ J for J in Jt)
X = sum(sum(Js[k][a] @ Js[l][a] for a in range(3)) for k, l in [(0, 1), (1, 2), (2, 0)])
w, V = np.linalg.eigh(C)
lams = [0, 2, 6, 12, 20, 30, 42]
states = [int(np.sum(np.abs(w - lam) < 1e-8)) for lam in lams]
spins = [int(round((-1 + np.sqrt(1 + 4 * lam)) / 2)) for lam in lams]
copies = [s // (2 * J + 1) for s, J in zip(states, spins)]
check("D1", "sector sizes 1, 9, 25, 28, 27, 22, 13; kernel 25 + 22 = 47",
      states == [1, 9, 25, 28, 27, 22, 13] and copies == [1, 3, 5, 4, 3, 2, 1])

rX = np.abs(C - 18 * I - 2 * X).max()
check("D2", "C = 18I + 2X (total pair exchange)", rX < TOL, f"{rX:.1e}")
K = (C - 6 * I) @ (C - 30 * I)
rK = np.abs(K - 4 * (X @ X - 36 * I)).max()
check("D3", "K = (C - 6I)(C - 30I) = 4(X^2 - 36I)", rK < TOL, f"{rK:.1e}")
mu = float(np.trace(C).real / 125); sigma2 = float(np.trace((C - mu * I) @ (C - mu * I)).real / 125)
check("D3b", "moments mu = 18, sigma^2 = 144; mu -+ sigma = 6, 30", abs(mu - 18) < TOL and abs(sigma2 - 144) < TOL)

Hp, Hm = V[:, np.abs(w - 6) < 1e-8], V[:, np.abs(w - 30) < 1e-8]
E = np.hstack([Hp, Hm]); P = E @ E.conj().T
S = (18 * I - C) @ P / 12
rS = max(np.abs(S + X @ P / 6).max(), np.abs(S @ S - P).max(), np.abs(S - S.conj().T).max())
trS = float(np.trace(S).real)
check("D4", "grading S = (18I - C)P/12 = -XP/6, S^2 = P, Tr S = 3", rS < TOL and abs(trS - 3) < TOL, f"{rS:.1e}")
Pp, Pm = (P + S) / 2, (P - S) / 2
rP = max(np.abs(Pp - (30 * I - C) @ P / 24).max(), np.abs(Pm - (C - 6 * I) @ P / 24).max(), np.abs(Pp @ Pm).max())
check("D5", "P+ = (30I - C)P/24, P- = (C - 6I)P/24, P+P- = 0, ranks 25 and 22",
      rP < TOL and round(np.trace(Pp).real) == 25 and round(np.trace(Pm).real) == 22, f"{rP:.1e}")

# ---------------------------------------------------------------- pairing and index (in the basis [H+ | H-])
Sloc = np.diag([1.0] * 25 + [-1.0] * 22)


def pairing_stats(r):
    Q, _ = np.linalg.qr(rng.normal(size=(25, 22)) + 1j * rng.normal(size=(25, 22)))
    A = (Q.conj().T)[:22].copy()
    A[r:] = 0                                            # rank r
    D = np.block([[np.zeros((25, 25)), A.conj().T], [A, np.zeros((22, 22))]])
    ev, vec = np.linalg.eigh(D)
    z = vec[:, np.abs(ev) < 1e-9]
    return D, ev, float(np.trace(z.conj().T @ Sloc @ z).real), np.abs(Sloc @ D + D @ Sloc).max()


D, ev, chir, anti = pairing_stats(22)
counts = [int(np.sum(np.abs(ev - t) < 1e-9)) for t in (-1, 0, 1)]
check("D6", "full pairing: {S, D} = 0, spectrum (-1, 0, +1) with multiplicities (22, 3, 22)",
      anti < TOL and counts == [22, 3, 22], f"{counts}")
index_by_rank = {}
for r in (0, 7, 15, 22):
    _, evr, chr_, _ = pairing_stats(r)
    index_by_rank[r] = {"zero_modes": int(np.sum(np.abs(evr) < 1e-9)), "trace_S_on_kernel": round(chr_, 9)}
check("D7", "every rank r leaves 47 - 2r = 3 + 2(22 - r) zero modes with Tr S = 3",
      all(v["zero_modes"] == 47 - 2 * r and abs(v["trace_S_on_kernel"] - 3) < 1e-8 for r, v in index_by_rank.items()))

# ---------------------------------------------------------------- symmetry: commutant M5 + M2, no invariant odd operator
Je = [E.conj().T @ J @ E for J in Jt]
S47 = E.conj().T @ S @ E
Lc = np.vstack([np.kron(np.eye(47), Ja) - np.kron(Ja.T, np.eye(47)) for Ja in Je])
sv = np.linalg.svd(Lc, compute_uv=False)
commutant_dim = int(np.sum(sv < 1e-8))
Lo = np.vstack([Lc, np.kron(np.eye(47), S47) + np.kron(S47.T, np.eye(47))])
odd_invariant_dim = int(np.sum(np.linalg.svd(Lo, compute_uv=False) < 1e-8))
check("D8", "SU(2)-invariant operators on E47 form M5 + M2 (dim 29); invariant odd operators: none",
      commutant_dim == 29 and odd_invariant_dim == 0, f"commutant {commutant_dim}, odd {odd_invariant_dim}")

# symmetric dial with 7 distinct hours: phases on the 5 + 2 copies (isotypic multiplets via the lowering operator)
Jm_t = Jt[0] - 1j * Jt[1]


def ladders(Hsec, J):
    Mz = Hsec.conj().T @ Jt[2] @ Hsec
    mz, U = np.linalg.eigh(Mz)
    tops = (Hsec @ U)[:, np.abs(mz - J) < 1e-8]
    out = []
    for t in tops.T:
        lad, v = [t], t
        for _ in range(2 * J):
            v = Jm_t @ v; v = v / np.linalg.norm(v); lad.append(v)
        out.append(np.array(lad).T)
    return out


blocks = ladders(Hp, 2) + ladders(Hm, 5)
thetas = 2 * np.pi * np.array([0, 9, 19, 28, 38, 14, 33]) / 47
Usym = sum(np.exp(1j * th) * (Bk @ Bk.conj().T) for th, Bk in zip(thetas, blocks)) + (I - P)
rinv = max(np.abs(Usym @ J - J @ Usym).max() for J in Jt)
reps = []                                                    # distinct eigenvalues on the unit circle
for z in np.linalg.eigvals(E.conj().T @ Usym @ E):
    if all(abs(z - q) > 1e-6 for q in reps):
        reps.append(z)
hours = len(reps)
check("D9", "a symmetric dial with 7 distinct hours exists (5 copies x 5 states, 2 copies x 11)",
      rinv < TOL and hours == 7 and [b.shape[1] for b in blocks] == [5] * 5 + [11] * 2, f"{rinv:.1e}")

# ---------------------------------------------------------------- the two clocks
zeta = np.exp(2j * np.pi / 47)
U0 = zeta * Pp + Pm + (I - P)
order = next(k for k in range(1, 100) if np.abs(np.linalg.matrix_power(U0, k) - I).max() < 1e-9)
r0 = max(np.abs(U0 @ J - J @ U0).max() for J in Jt)
f_diag = np.diag(E.conj().T @ (C / 24) @ E).real           # the linear clock f = omega C/24, omega = 1
step = np.exp(-1j * 2 * np.pi / 47 * f_diag)                 # one 1/47 turn, on [H+ | H-]
ratio = step[0] / step[-1]
check("D10", "U0 = zeta on H+, 1 on H-: order 47, commutes with SU(2); one 1/47 turn of f = C/24 is U0 up to phase",
      order == 47 and r0 < TOL and abs(ratio - zeta) < 1e-12, f"order {order}")
Qb, _ = np.linalg.qr(rng.normal(size=(47, 47)) + 1j * rng.normal(size=(47, 47)))
U1 = Qb @ np.diag(zeta ** np.arange(47)) @ Qb.conj().T
U1t = E @ U1 @ E.conj().T + (I - P)
r1 = max(np.abs(U1t @ J - J @ U1t).max() for J in Jt)
check("D11", "a 47-hour clock U1 (+) I78 keeps E47 but fails to commute with SU(2)",
      np.abs(U1t @ P - P @ U1t).max() < TOL and r1 > 1e-3, f"commutator {r1:.3f}")

# ---------------------------------------------------------------- index along the Fibonacci-Lucas tower (exact)
F, Lu = [0, 1], [2, 1]
for i in range(2, 30):
    F.append(F[-1] + F[-2]); Lu.append(Lu[-1] + Lu[-2])
tower = [{"k": k, "two_j_plus_1": F[k], "two_J_plus_1": Lu[k], "dim_H_plus": F[k] ** 2, "dim_H_minus": F[k - 2] * Lu[k],
          "kernel": F[k] ** 2 + F[k - 2] * Lu[k], "carrier": F[k] ** 3, "index": F[k] ** 2 - F[k - 2] * Lu[k]} for k in range(3, 14, 2)]
check("D12", "index F_k^2 - F_(k-2) L_k = F_(k-3) F_(k-1): 0, 3, 24, 168, 1155, 7920",
      all(t["index"] == F[t["k"] - 3] * F[t["k"] - 1] for t in tower) and [t["index"] for t in tower] == [0, 3, 24, 168, 1155, 7920])

data = {
    "schema": "e47-dial-of-47-hours/1", "certificate": "MC-E47-DIAL-20261126-001",
    "carrier": {"spin": 2, "dim": 125, "single_site_casimir": 6},
    "sectors": [{"lambda": l, "J": J, "copies": c, "states": s} for l, J, c, s in zip(lams, spins, copies, states)],
    "moments": {"mu": 18, "sigma_squared": 144, "selected": [6, 30]},
    "exchange_X_by_sector": {str(l): (l - 18) // 2 for l in lams},
    "kernel": {"dim": 47, "H_plus": {"lambda": 6, "dim": 25, "structure": "C^5 (x) V2", "S": 1},
               "H_minus": {"lambda": 30, "dim": 22, "structure": "C^2 (x) V5", "S": -1}, "trace_S": 3},
    "pairing": {"full_rank_spectrum": {"-1": 22, "0": 3, "+1": 22}, "index": 3, "zero_modes_by_rank": {str(r): v["zero_modes"] for r, v in index_by_rank.items()}},
    "symmetry": {"commutant_dim": commutant_dim, "commutant": "M5 (+) M2", "invariant_odd_operators": odd_invariant_dim,
                 "max_hours_symmetric_dial": 7, "hour_weights": [5, 5, 5, 5, 5, 11, 11]},
    "clocks": {"U0": {"hours": 2, "order": order, "commutes_with_SU2": True},
               "U1": {"hours": 47, "commutes_with_SU2": False},
               "linear_clock_on_E47": "f(C)P = omega(3/4 P - 1/2 S)"},
    "tower": tower,
}
path = HERE / "dial_data.json"
if "--write" in sys.argv:
    path.write_text(json.dumps(data, indent=1, ensure_ascii=False) + "\n", encoding="utf-8")
    print("wrote", path.relative_to(HERE.parent))
else:
    committed = json.loads(path.read_text(encoding="utf-8"))
    check("DATA", "dial_data.json matches this replay", committed == data)
print(f"RESULT: {len(results) - len(failures)}/{len(results)} PASS")
sys.exit(1 if failures else 0)
