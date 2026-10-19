"""Temporal Quantum Graph Neural Network on the E47 carrier: validation runner.

Run from the repository root:   python validation/tqgnn_validation.py
Writes results/results.json; exits non-zero if any check fails.
"""
import json
import os
import sys
import time

import numpy as np
from scipy.integrate import trapezoid

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from tqgnn.carrier import Carrier, SPECTRUM, components, distances           # noqa: E402
from tqgnn.clocks import (click_time_precision, cycle_unitary, matrix_order,  # noqa: E402
                          phase_time)
from tqgnn.learning import jacobian_rank, residuals, simulate, train           # noqa: E402
from tqgnn.lindblad import counting_rates, evolve, liouvillian, steady_state   # noqa: E402
from tqgnn.network import TQGNN                                                # noqa: E402
from tqgnn.node import (ClockMemoryNode, adiabatic_coherence_rate,             # noqa: E402
                        adiabatic_dephasing, resonant_fano, stationary_leak)
from tqgnn.spacetime import METRIC, Worldline, lapse, pullback_metric          # noqa: E402
from tqgnn.stochastic import conditional_density, jump_trajectories, thin      # noqa: E402

T_START = time.time()
checks, values = {}, {}


def check(name, cond, value=None):
    checks[name] = bool(cond)
    if value is not None:
        values[name] = value
    tag = "PASS" if cond else "FAIL"
    suffix = "" if value is None else f" | {value}"
    print(f"{tag:4s}  {name}{suffix}", flush=True)


def section(title):
    print(f"\n# {title}", flush=True)


OMEGA = 2.0

# =============================================================================
section("A. Carrier as a graph: E47 is a spectral band, Gamma is a 4-local filter")
car = Carrier()
spec = {int(k): int(v) for k, v in zip(*np.unique(car.labels, return_counts=True))}
check("A1 Casimir spectrum and eigenspace dimensions", spec == SPECTRUM, spec)
check("A2 rank P = 47", round(np.trace(car.P)) == 47, f"tr P = {np.trace(car.P):.12f}")

adj = car.adjacency()
mm = car.site_m()
steps = sorted({tuple(int(x) for x in mm[u] - mm[v]) for u, v in zip(*np.nonzero(adj))})
face = sorted({tuple(int(x) for x in s) for s in
               [np.eye(3, dtype=int)[a] - np.eye(3, dtype=int)[b] for a in range(3) for b in range(3) if a != b]})
check("A3 C is a graph operator whose edges are face diagonals e_a - e_b", steps == face, f"{len(steps)} edge steps")

comps = sorted(components(adj), key=lambda c: mm[c[0]].sum())
sizes = [len(c) for c in comps]
check("A4 graph has 13 components (planes m1+m2+m3 = M)", sizes == [1, 3, 6, 10, 15, 18, 19, 18, 15, 10, 6, 3, 1], sizes)

D = distances(adj)
far = float(np.abs(car.Gamma[D > 4]).max())
check("A5 Gamma = I - K^2/99144 vanishes beyond graph distance 4", far == 0.0, f"max |Gamma_uv|, d > 4: {far}")

gaps = [abs(np.linalg.norm(np.linalg.matrix_power(car.Gamma, n) - car.P, 2) - (15 / 17) ** n) for n in range(1, 11)]
check("A6 ||Gamma^n - P||_2 = (15/17)^n, n = 1..10", max(gaps) < 1e-12,
      f"max |err| = {max(gaps):.2e}; (15/17)^220 = {(15 / 17) ** 220:.4e}")

pairs24 = [(x, y) for x in SPECTRUM for y in SPECTRUM if y - x == 24]
check("A7 only (6, 30) has Casimir gap 24", pairs24 == [(6, 30)], pairs24)

# =============================================================================
section("B. Time formalism: metric, lapse, time-chain identity")
g, phi_res, u, Wcopy = pullback_metric(car)
check("B1 pullback metric B4^dag (P6 - P0) B4 = diag(-1, 2, 8, 2)", np.allclose(g, METRIC, atol=1e-12),
      f"diag = {np.round(np.diag(g), 12).tolist()}, intertwiner residual {phi_res:.2e}")

N1 = lapse([0.3, 0, 0])
try:
    lapse([0.5, 0.2, 0.3])
    spacelike_rejected = False
except ValueError:
    spacelike_rejected = True
check("B2 lapse N = sqrt(1 - 2v1^2 - 8v2^2 - 2v3^2); spacelike worldline rejected",
      abs(N1 - np.sqrt(0.82)) < 1e-15 and spacelike_rejected, f"N(0.3,0,0) = {N1:.12f}")

node = ClockMemoryNode(car)
Hn, Rn = node.reduced(OMEGA, 0.6, 2.5)
wl_t = Worldline(lambda t: np.array([0.3 * np.sin(0.5 * t), 0.0, 0.2]))
psi_n = np.array([1, 0, 1, 0], dtype=complex) / np.sqrt(2)
rho_n = np.outer(psi_n, psi_n.conj())
tg = np.linspace(0, 8, 17)
coord = evolve(lambda t: (wl_t.N(t) * Hn, [np.sqrt(wl_t.N(t)) * Rn]), rho_n, tg)
taus = np.array([wl_t.tau(t) for t in tg])
proper = evolve(lambda s: (Hn, [Rn]), rho_n, taus)
tc_err = float(np.abs(coord - proper).max())
check("B3 time-chain identity: d rho/dt = N(t) L_tau rho equals proper-time evolution at tau(t)",
      tc_err < 1e-9, f"max |delta rho| = {tc_err:.2e}")

K2u = car.K @ car.K @ u
gu = car.Gamma @ u
check("B4 time direction u lies outside E47: K^2 u = 32400 u, Gamma u = 0.6732026144 u",
      np.linalg.norm(K2u - 32400 * u) < 1e-7 and abs(np.vdot(u, gu).real - (1 - 32400 / 99144)) < 1e-12,
      f"<u|Gamma|u> = {np.vdot(u, gu).real:.10f}, P u = {np.linalg.norm(car.P @ u):.1e}")

# =============================================================================
section("C. Node: E47 clock-memory unit")
check("C1 V inside E47, W outside, W^dag V = 0, [V W] orthonormal",
      np.linalg.norm(car.P @ node.V - node.V) < 1e-12 and np.linalg.norm(car.P @ node.W) < 1e-12
      and np.linalg.norm(node.U.T @ node.U - np.eye(4)) < 1e-12)

Hf, Rf = node.full(OMEGA, 0.7, 3.0)
Hr, Rr = node.reduced(OMEGA, 0.7, 3.0)
psi_r = np.array([1, 0, 1j, 0]) / np.sqrt(2)
rr = np.outer(psi_r, psi_r.conj())
ts = np.linspace(0, 6, 13)
red = evolve(lambda t: (Hr, [Rr]), rr, ts)
ful = evolve(lambda t: (Hf, [Rf]), node.U @ rr @ node.U.T, ts, rtol=1e-10, atol=1e-12)
dev = max(np.abs(ful[k] - node.U @ red[k] @ node.U.T).max() for k in range(len(ts)))
rdev = max(abs(np.trace(car.P @ ful[k]).real - (red[k][0, 0] + red[k][2, 2]).real) for k in range(len(ts)))
check("C2 exact reduction: 125-dim carrier dynamics = 4-dim logical (x) flag dynamics",
      dev < 1e-9 and rdev < 1e-9, f"max |delta rho| = {dev:.2e}, retention {rdev:.2e}")

H0, _ = node.reduced(OMEGA, 0.0, 0.0)
tsc = np.linspace(0, 10, 41)
r0 = evolve(lambda t: (H0, []), np.outer([1, 0, 1, 0], [1, 0, 1, 0]).astype(complex) / 2, tsc)
Xl = np.zeros((4, 4)); Xl[0, 2] = Xl[2, 0] = Xl[1, 3] = Xl[3, 1] = 1
cerr = max(abs(np.trace(Xl @ r).real - np.cos(OMEGA * t)) for r, t in zip(r0, tsc))
rerr = max(abs((r[0, 0] + r[2, 2]).real - 1) for r in r0)
check("C3 clock H = omega C/24: <O> = cos(omega tau) with Tr(P rho) = 1", cerr < 1e-9 and rerr < 1e-12,
      f"cos error {cerr:.1e}")

leak_rows = {}
ok = True
for gam, ret_ref, clk_ref in ((20.0, 99.019608, 0.19607843), (40.0, 99.751244, 0.09950249)):
    Hb, Rb = node.branch(0, 0.0, 1.0, gam)
    q = steady_state(Hb, [Rb])[1, 1].real
    leak_rows[gam] = (round(100 * (1 - q), 6), round(gam * q, 8))
    ok &= abs(100 * (1 - q) - ret_ref) < 5e-7 and abs(gam * q - clk_ref) < 5e-9
    ok &= abs(q - stationary_leak(1.0, gam)) < 1e-12
check("C4 matched recovery, c = 1: retention 99.019608 % / 99.751244 %, clicks 0.19607843 / 0.09950249",
      ok, leak_rows)

Da, Db = node.detunings(OMEGA)
errs = []
for l, dl in ((0, Da), (1, Db)):
    for cc, gg in ((1.0, 4.0), (0.5, 1.0)):
        Hb, Rb = node.branch(l, OMEGA, cc, gg)
        errs.append(abs(steady_state(Hb, [Rb])[1, 1].real - stationary_leak(cc, gg, dl)))
check("C5 detuned leak q = 4c^2/(gamma^2 + 8c^2 + 4 Delta^2) per clock branch", max(errs) < 1e-12,
      f"Delta_a = {Da}, Delta_b = {Db}; max err {max(errs):.1e}")

ferr = []
for cc, gg in ((1.0, 20.0), (1.0, 40.0), (1.0, 4.0)):
    Hb, Rb = node.branch(0, 0.0, cc, gg)
    mean, F = counting_rates(Hb, [Rb], [0])
    ferr.append(max(abs(F - resonant_fano(cc, gg)), abs(mean - gg * stationary_leak(cc, gg))))
check("C6 recovery clicks are sub-Poissonian: F = 1 - 24c^2 gamma^2/(gamma^2 + 8c^2)^2", max(ferr) < 1e-7,
      f"F(c=1, gamma=20) = {resonant_fano(1, 20):.6f}, F(c=1, gamma=4) = {resonant_fano(1, 4):.6f}")


def clock_mode(H, R, omega):
    ev = np.linalg.eigvals(liouvillian(H, [R]))
    cand = [e for e in ev if abs(abs(e.imag) - omega) < 0.3 * omega]
    e = max(cand, key=lambda e: e.real)
    return -e.real, abs(e.imag)


rel = []
for cc in (0.1, 0.05):
    Hc, Rc = node.reduced(OMEGA, cc, 4.0)
    rate, freq = clock_mode(Hc, Rc, OMEGA)
    lam = adiabatic_coherence_rate(cc, 4.0, Da, Db, OMEGA)
    rel.append((rate / -lam.real - 1, freq - lam.imag))
check("C7 linear clock + recovery: clock dephases at gamma c^2 |z_a - z_b|^2 / 2 (second order in c/gamma)",
      abs(rel[1][0]) < 5e-3 and abs(rel[0][0]) / abs(rel[1][0]) > 3.5 and abs(rel[1][1]) < 1e-6,
      f"rate rel. error {rel[0][0]:.2e} (c=0.1), {rel[1][0]:.2e} (c=0.05); freq error {rel[1][1]:.1e}")

Hm, Rm = node.reduced(OMEGA, 0.7, 3.0, "matched")
mrate, mfreq = clock_mode(Hm, Rm, OMEGA)
rm = evolve(lambda t: (Hm, [Rm]), np.outer([1, 0, 1, 0], [1, 0, 1, 0]).astype(complex) / 2, tsc)
cohm = max(abs(abs(r[0, 2] + r[1, 3]) - 0.5) for r in rm)
check("C8 matched spectral clock (equal detunings): clock coherence exactly preserved under leak + recovery",
      abs(mrate) < 1e-12 and cohm < 1e-9 and abs(mfreq - OMEGA) < 1e-12,
      f"dephasing {mrate:.1e}, |coherence| - 1/2 = {cohm:.1e}, f(42) moved to {node.clock_function(OMEGA, 'matched')[42]}")

# =============================================================================
section("D. Network: temporal quantum graph neural network")
wls = [Worldline([0, 0, 0]), Worldline([0.3, 0, 0]),
       Worldline(lambda t: np.array([0.0, 0.15 * np.sin(0.4 * t), 0.2]))]
edges = [(0, 1), (1, 2)]
psi3 = TQGNN.product_state([TQGNN.clock_state() for _ in range(3)])
rho3 = np.outer(psi3, psi3.conj())
tn = np.linspace(0, 12, 61)

free = TQGNN(3, edges, omega=OMEGA, worldlines=wls)
rf = free.run_dynamic(rho3, tn)
coh = np.array([free.readout(r, t)[1] for r, t in zip(rf, tn)])
taus3 = np.array([[w.tau(t) for w in wls] for t in tn])
est = np.column_stack([phase_time(coh[:, v], OMEGA) for v in range(3)])
check("D1 each node's clock phase reads its own proper time (time dilation)", np.abs(est - taus3).max() < 1e-9,
      f"tau(12) = {np.round(taus3[-1], 6).tolist()}, max err {np.abs(est - taus3).max():.1e}")

net = TQGNN(3, edges, omega=OMEGA, c=0.3, gamma=2.0, J=0.25, Jzz=0.05, worldlines=wls)
rn = net.run_dynamic(rho3, tn, rtol=1e-9, atol=1e-11)
tr = max(abs(np.trace(r) - 1) for r in rn)
me = min(np.linalg.eigvalsh(r).min() for r in rn)
check("D2 full network (leak, recovery, exchange, Ising, node lapses) stays a density matrix",
      tr < 1e-12 and me > -1e-8, f"trace err {tr:.1e}, min eigenvalue {me:.1e}")

wl_same = [Worldline(lambda t: np.array([0.25 * np.cos(0.3 * t), 0.1, 0.0]))] * 3
uni = TQGNN(3, edges, omega=OMEGA, c=0.3, gamma=2.0, J=0.25, worldlines=wl_same)
tu = np.linspace(0, 6, 7)
ru = uni.run_dynamic(rho3, tu, rtol=1e-10, atol=1e-12)
Hp, Lp = uni.proper_generator()
rp = evolve(lambda s: (Hp, Lp), rho3, np.array([wl_same[0].tau(t) for t in tu]), rtol=1e-10, atol=1e-12)
check("D3 equal lapses: network obeys the time-chain identity", np.abs(ru - rp).max() < 1e-8,
      f"max |delta rho| = {np.abs(ru - rp).max():.1e}")

w2 = [Worldline([0, 0, 0]), Worldline([0.3, 0, 0])]
J2 = 0.05
pair = TQGNN(2, [(0, 1)], omega=OMEGA, J=J2, worldlines=w2)
Nn = pair.lapses(0.0)
Je = 0.5 * (Nn[0] + Nn[1]) * J2
delta = (Nn[0] - Nn[1]) * OMEGA
OmR = np.sqrt(4 * Je ** 2 + delta ** 2)
b_a = TQGNN.product_state([TQGNN.clock_state(np.pi), TQGNN.clock_state(0.0)])
a_b = TQGNN.product_state([TQGNN.clock_state(0.0), TQGNN.clock_state(np.pi)])
rr2 = pair.run_dynamic(np.outer(b_a, b_a.conj()), np.array([0.0, np.pi / OmR]))[-1]
Ptr = np.vdot(a_b, rr2 @ a_b).real
Ppred = 4 * Je ** 2 / (4 * Je ** 2 + delta ** 2)
check("D4 message passing between clocks on different worldlines is detuned by (N_j - N_k) omega",
      abs(Ptr - Ppred) < 1e-9, f"max transfer {Ptr:.10f} = 4J^2/(4J^2 + delta^2) = {Ppred:.10f}, delta = {delta:.6f}")

# =============================================================================
section("S. From stochastic to dynamic: detection efficiency eta and recovery strength gamma")
gen1 = lambda t: (wl_t.N(t) * Hn, [np.sqrt(wl_t.N(t)) * Rn])
te = np.linspace(0, 8, 9)
exact1 = evolve(gen1, rho_n, te)
st0, ck0 = conditional_density(gen1, rho_n, te, 4, 0.0, 0.004, np.random.default_rng(1))
d0 = float(np.abs(st0 - exact1[:, None]).max())
check("S1 eta = 0: no record, conditional state = GKSL solution (dynamic system)",
      d0 < 1e-8 and sum(len(c) for c in ck0) == 0, f"max |delta rho| = {d0:.1e}")

Pn = np.diag([1, 0, 1, 0]).astype(complex)
Cn = np.zeros((4, 4), complex); Cn[0, 2] = Cn[1, 3] = 1
er = np.array([np.trace(Pn @ r).real for r in exact1])
ec = np.array([np.trace(Cn @ r) for r in exact1])
tf = np.linspace(0, 8, 401)
ef = evolve(gen1, rho_n, tf)
lam_t = np.array([wl_t.N(t) * 2.5 * (1 - np.trace(Pn @ r).real) for t, r in zip(tf, ef)])
mean_clicks = trapezoid(lam_t, tf)
zmax, var_by_eta, click_z = 0.0, {0.0: 0.0}, []
for k, eta in enumerate((0.25, 0.5, 0.75, 1.0)):
    st, ck = conditional_density(gen1, rho_n, te, 1500, eta, 0.005, np.random.default_rng(100 + k))
    ret = np.einsum("ij,tmji->tm", Pn, st).real
    cc = np.einsum("ij,tmji->tm", Cn, st)
    se = lambda x: x.std(1) / np.sqrt(x.shape[1]) + 1e-15
    z = max((np.abs(ret.mean(1) - er) / se(ret))[1:].max(),
            (np.abs(cc.real.mean(1) - ec.real) / se(cc.real))[1:].max(),
            (np.abs(cc.imag.mean(1) - ec.imag) / se(cc.imag))[1:].max())
    zmax = max(zmax, z)
    var_by_eta[eta] = float(ret[-1].var())
    n = np.array([len(c) for c in ck])
    click_z.append(abs(n.mean() - eta * mean_clicks) / (n.std() / np.sqrt(len(n))))
check("S2 unravelling invariance: E[rho_c] = rho for eta = 0.25, 0.5, 0.75, 1", zmax < 4.5,
      f"max |z| = {zmax:.2f} over retention and clock coherence")
check("S3 detected clicks have mean eta * integral of lambda(t)", max(click_z) < 4.0,
      f"exact mean {mean_clicks:.4f} at eta = 1; max |z| = {max(click_z):.2f}")
vs = [var_by_eta[e] for e in sorted(var_by_eta)]
check("S4 conditional variance rises with eta from 0 (stochastic) and vanishes at eta = 0 (dynamic)",
      all(np.diff(vs) > 0), {str(k): round(v, 6) for k, v in sorted(var_by_eta.items())})

traj_z = 0.0
st3, ck3 = jump_trajectories(net.generator, psi3, np.linspace(0, 6, 7), 2000, 0.005, np.random.default_rng(12))
exact3 = net.run_dynamic(rho3, np.linspace(0, 6, 7), rtol=1e-9, atol=1e-11)
for v in range(3):
    Pv, Cv = net.P_node[v], net.coh_node[v]
    ret = np.einsum("tmi,ij,tmj->tm", st3.conj(), Pv, st3).real
    cv = np.einsum("tmi,ij,tmj->tm", st3.conj(), Cv, st3)
    erv = np.array([np.trace(Pv @ r).real for r in exact3])
    ecv = np.array([np.trace(Cv @ r) for r in exact3])
    se = lambda x: x.std(1) / np.sqrt(x.shape[1]) + 1e-15
    traj_z = max(traj_z, (np.abs(ret.mean(1) - erv) / se(ret))[1:].max(),
                 (np.abs(cv.real.mean(1) - ecv.real) / se(cv.real))[1:].max(),
                 (np.abs(cv.imag.mean(1) - ecv.imag) / se(cv.imag))[1:].max())
check("S5 network click trajectories (eta = 1) average to the network GKSL flow", traj_z < 4.5,
      f"2000 trajectories, 64-dim, max |z| = {traj_z:.2f}")

fano_ok = []
for eta in (0.25, 0.5, 0.75):
    Hb, Rb = node.branch(0, 0.0, 1.0, 20.0)
    m1, F1 = counting_rates(Hb, [Rb], [0], 1.0)
    me_, Fe = counting_rates(Hb, [Rb], [0], eta)
    fano_ok.append(max(abs(Fe - (1 + eta * (F1 - 1))), abs(me_ - eta * m1)))
check("S6 click statistics thin exactly: rate eta*lambda, Fano 1 + eta(F - 1)", max(fano_ok) < 1e-7,
      f"max err {max(fano_ok):.1e}")

Hb, Rb = node.branch(0, 0.0, 1.0, 4.0)
lam4, F4 = 4.0 * stationary_leak(1, 4), resonant_fano(1, 4)
tau = 300.0
_, ckc = jump_trajectories(lambda t: (Hb, [Rb]), np.array([1, 0], complex), np.array([0.0, tau]), 1000, 0.005,
                           np.random.default_rng(5))
n1 = np.array([len(c) for c in ckc])
n5 = np.array([len(c) for c in thin(ckc, 0.5, np.random.default_rng(9))])
r1, p1 = (n1 / (lam4 * tau)).std(), click_time_precision(F4, lam4, tau)
r5, p5 = (n5 / (0.5 * lam4 * tau)).std(), click_time_precision(1 + 0.5 * (F4 - 1), 0.5 * lam4, tau)
check("S7 click clock: relative error of tau from counts = sqrt(F_eta/(eta lambda tau))",
      abs(r1 / p1 - 1) < 0.1 and abs(r5 / p5 - 1) < 0.1,
      f"eta=1: {r1:.4f} vs {p1:.4f}; eta=0.5: {r5:.4f} vs {p5:.4f} (tau = 300, 1000 trajectories)")

gs = np.array([16.0, 32.0, 64.0, 128.0])
rates = [clock_mode(*node.reduced(OMEGA, 0.5, g_), OMEGA)[0] for g_ in gs]
slope = np.polyfit(np.log(gs), np.log(rates), 1)[0]
lam_g = [g_ * stationary_leak(1, g_) for g_ in (4.0, 8.0, 16.0, 32.0)]
check("S8 Zeno route: stronger recovery makes clicks rarer and the clock deterministic (gamma_phi ~ gamma^-3)",
      abs(slope + 3) < 0.05 and all(np.diff(lam_g) < 0),
      f"log-log slope {slope:.4f}; click rate {np.round(lam_g, 5).tolist()} for gamma = 4, 8, 16, 32")

# =============================================================================
section("L. Learning: the network recovers each node's lapse from clock data")
n3, e3 = 3, [(0, 1), (1, 2)]
N_true = np.array([lapse([0, 0, 0]), lapse([0.3, 0, 0]), lapse([0, 0.15, 0.2])])
J_true = np.array([0.15, 0.10])
s = lambda th, ph: np.array([np.cos(th / 2), np.exp(1j * ph) * np.sin(th / 2)])
psi_l = np.kron(np.kron(s(np.pi / 2, 0.0), s(0.6, np.pi / 3)), s(2.2, -np.pi / 4))
times = np.linspace(0, 30, 121)
data = simulate(n3, e3, N_true, [OMEGA] * 3, J_true, psi_l, times)
rank, sv = jacobian_rank(n3, e3, N_true, np.array([OMEGA] * 3), J_true, psi_l, times)
check("L1 uncalibrated clocks: (N_v, omega_v, J_e) has 3 flat directions (only N_v omega_v, N_e J_e seen)",
      rank == 5, f"rank {rank} of 8; smallest singular values {np.round(sv[-3:], 12).tolist()}")
fit = train(data, n3, e3, OMEGA, psi_l, times)
lerr = float(np.abs(fit.x[:3] - N_true).max())
check("L2 calibrated clocks: trained lapses equal the metric lapses", lerr < 1e-8,
      f"fitted {np.round(fit.x[:3], 10).tolist()} vs {np.round(N_true, 10).tolist()}; J {np.round(fit.x[3:], 10).tolist()}")
rng = np.random.default_rng(3)
noisy = data + 0.002 * (rng.standard_normal(data.shape) + 1j * rng.standard_normal(data.shape))
fitn = train(noisy, n3, e3, OMEGA, psi_l, times)
check("L3 with readout noise 0.002 the lapses are still recovered to < 1e-3",
      np.abs(fitn.x[:3] - N_true).max() < 1e-3, f"max err {np.abs(fitn.x[:3] - N_true).max():.1e}")

# =============================================================================
section("P. Prime clock")
B47 = np.hstack([car.basis(6), car.basis(30)])
U47 = cycle_unitary(B47)
check("P1 47-cycle on E47 has order 47 and preserves E47",
      matrix_order(U47, 60) == 47 and np.linalg.norm(U47 @ car.P - car.P @ U47) < 1e-12)
X = np.kron(U47, np.roll(np.eye(5), 1, axis=0))
I = np.eye(X.shape[0])
mp = np.linalg.matrix_power
check("P2 combined with a 5-cycle the period is 235 = lcm(47, 5)",
      np.linalg.norm(mp(X, 235) - I) < 1e-9 and np.linalg.norm(mp(X, 47) - I) > 1 and np.linalg.norm(mp(X, 5) - I) > 1)

# =============================================================================
n_pass = sum(checks.values())
all_pass = n_pass == len(checks)
print("\n" + ("TQGNN_TIME_PASS" if all_pass else "VALIDATION_FAIL"))
print(f"{n_pass}/{len(checks)} checks PASS   ({time.time() - T_START:.0f} s)")

os.makedirs(os.path.join(os.path.dirname(__file__), "..", "results"), exist_ok=True)
with open(os.path.join(os.path.dirname(__file__), "..", "results", "results.json"), "w") as f:
    json.dump({"status": "TQGNN_TIME_PASS" if all_pass else "VALIDATION_FAIL",
               "passed": n_pass, "total": len(checks), "checks": checks,
               "values": {k: (v if isinstance(v, (int, float, str, list, dict)) else str(v)) for k, v in values.items()}},
              f, indent=2, default=str)
sys.exit(0 if all_pass else 1)
