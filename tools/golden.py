#!/usr/bin/env python3
"""Reference values for the Quantum Clock engine, computed by the Python TQGNN package.

The JavaScript engine in ../engine is an independent re-implementation. Its test suite and its
in-browser precision certificate compare against the numbers written here (data/golden.json).
Propagation here is exact (scipy.linalg.expm of the Liouvillian), so agreement is limited only by
floating point.

    python tools/golden.py            # write data/golden.json from the tqgnn package in this repository
    python tools/golden.py --check    # recompute and compare with the committed data/golden.json
"""
import argparse
import hashlib
import json
import math
import platform
import sys
from pathlib import Path

import numpy as np
import scipy
from scipy.linalg import expm

ROOT = Path(__file__).resolve().parent.parent
ap = argparse.ArgumentParser()
ap.add_argument("--tqgnn", default=str(ROOT), help="folder that contains tqgnn/ (default: this repository)")
ap.add_argument("--out", default=str(ROOT / "data" / "golden.json"), help="file to write, or with --check, to compare against")
ap.add_argument("--check", action="store_true", help="compare with --out instead of writing it; exit 1 on any difference")
args = ap.parse_args()
sys.path.insert(0, args.tqgnn)

from tqgnn.carrier import Carrier, SPECTRUM  # noqa: E402
from tqgnn.lindblad import counting_rates, liouvillian  # noqa: E402
from tqgnn.network import TQGNN  # noqa: E402
from tqgnn.node import ClockMemoryNode, adiabatic_coherence_rate  # noqa: E402
from tqgnn.spacetime import Worldline, pullback_metric  # noqa: E402

car = Carrier()
node = ClockMemoryNode(car)
out = {"schema": "QCLOCK-GOLDEN-1.0", "generator": "tools/golden.py",
       "environment": {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__}}

g, resid, _, _ = pullback_metric(car)
out["metric"] = {"g": g.tolist(), "intertwiner_residual": float(resid)}
out["multiplicities"] = {str(k): int(np.sum(car.labels == k)) for k in SPECTRUM}


def cplx(z):
    return [float(np.real(z)), float(np.imag(z))]


def rho0(theta=np.pi / 2, phase=0.0):
    psi = np.kron(np.array([np.cos(theta / 2), np.exp(1j * phase) * np.sin(theta / 2)]), np.array([1, 0]))
    return np.outer(psi, psi.conj())


def readout(rho):
    return {"retention": float(np.real(rho[0, 0] + rho[2, 2])),
            "coherence": cplx(rho[2, 0] + rho[3, 1])}


node_cases = []
for clock in ("linear", "matched"):
    for omega, c, gamma in ((2.0, 1.0, 20.0), (2 * np.pi / 8, 1.0, 20.0), (1.0, 2.0, 5.0)):
        H, R = node.reduced(omega, c, gamma, clock)
        Lv = liouvillian(H, [R])
        r0 = rho0(np.pi / 2, 0.3)
        for tau in (0.5, 2.0, 10.0):
            r = (expm(Lv * tau) @ r0.reshape(-1)).reshape(4, 4)
            node_cases.append({"clock": clock, "omega": omega, "c": c, "gamma": gamma, "theta": np.pi / 2, "phase": 0.3,
                               "tau": tau, **readout(r), "rho": [[cplx(r[i, j]) for j in range(4)] for i in range(4)]})
out["node"] = node_cases

fcs = []
for gamma in (5.0, 20.0, 40.0):
    for delta in (0.0, 0.5, 1.0):
        H = np.array([[0, 1.0], [1.0, delta]], dtype=complex)
        R = np.sqrt(gamma) * np.array([[0, 1], [0, 0]], dtype=complex)
        mean, fano = counting_rates(H, [R], [0], eta=1.0, h=1e-3)
        fcs.append({"c": 1.0, "gamma": gamma, "delta": delta, "rate_fd": float(mean), "fano_fd": float(fano)})
out["fcs"] = fcs

modes = []
for omega, c, gamma in ((2.0, 1.0, 20.0), (2.0, 0.2, 20.0)):
    Ha, Ra = node.branch(0, omega, c, gamma, "linear")
    Hb, _ = node.branch(1, omega, c, gamma, "linear")
    I = np.eye(2)
    RdR = Ra.conj().T @ Ra
    Lba = -1j * (np.kron(Hb, I) - np.kron(I, Ha.T)) + np.kron(Ra, Ra.conj()) - 0.5 * np.kron(RdR, I) - 0.5 * np.kron(I, RdR.T)
    ev = np.linalg.eigvals(Lba)
    k = int(np.argmax(ev.real))
    da, db = node.detunings(omega, "linear")
    modes.append({"omega": omega, "c": c, "gamma": gamma, "slowest_ba": cplx(ev[k]),
                  "adiabatic_ab": cplx(adiabatic_coherence_rate(c, gamma, da, db, omega))})
out["coherence_modes"] = modes

net_cases = []
for clock in ("matched", "linear"):
    net = TQGNN(2, [(0, 1)], omega=2.0, c=0.5, gamma=10.0, J=0.3, Jzz=0.2,
                worldlines=[Worldline([0.0, 0.0, 0.0]), Worldline([0.3, 0.0, 0.0])], clock=clock)
    H, Ls = net.generator(0.0)
    Lv = liouvillian(H, Ls)
    psi = TQGNN.product_state([TQGNN.clock_state(np.pi / 2, 0.0), TQGNN.clock_state(np.pi / 3, 0.7)])
    r0 = np.outer(psi, psi.conj())
    for t in (0.5, 1.5):
        r = (expm(Lv * t) @ r0.reshape(-1)).reshape(16, 16)
        ret, coh, clicks = net.readout(r, t)
        net_cases.append({"clock": clock, "t": t, "lapses": net.lapses(0.0).tolist(),
                          "retention": ret.tolist(), "coherence": [cplx(z) for z in coh], "clicks": clicks.tolist()})
out["network"] = {"n": 2, "edges": [[0, 1]], "omega": 2.0, "c": 0.5, "gamma": 10.0, "J": 0.3, "Jzz": 0.2,
                  "velocities": [[0, 0, 0], [0.3, 0, 0]], "initial": [[np.pi / 2, 0.0], [np.pi / 3, 0.7]], "cases": net_cases}

src = Path(__file__).read_bytes()
out["generator_sha256"] = hashlib.sha256(src).hexdigest()

if not args.check:
    Path(args.out).write_text(json.dumps(out, indent=1) + "\n")
    print(f"wrote {args.out}: {len(node_cases)} node, {len(fcs)} fcs, {len(modes)} mode, {len(net_cases)} network cases")
    sys.exit(0)

# --check: every number must agree to floating point (the finite-difference counting statistics
# amplify rounding by about 1/h^2, hence rel_tol 1e-8); strings, keys and lengths must match exactly.
ref, new, bad = json.loads(Path(args.out).read_text()), json.loads(json.dumps(out)), []
worst = 0.0


def walk(a, b, path):
    global worst
    if path in ("environment", "generator_sha256"):
        return
    if isinstance(a, dict) and isinstance(b, dict):
        if a.keys() != b.keys():
            bad.append(f"{path or '.'}: keys {sorted(a)} != {sorted(b)}")
            return
        for k in a:
            walk(a[k], b[k], f"{path}.{k}" if path else k)
    elif isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b):
            bad.append(f"{path}: length {len(a)} != {len(b)}")
            return
        for i, (x, y) in enumerate(zip(a, b)):
            walk(x, y, f"{path}[{i}]")
    elif isinstance(a, (int, float)) and isinstance(b, (int, float)) and not isinstance(a, bool) and not isinstance(b, bool):
        worst = max(worst, abs(a - b))
        if not math.isclose(a, b, rel_tol=1e-8, abs_tol=1e-11):
            bad.append(f"{path}: committed {a!r}, computed {b!r}")
    elif a != b:
        bad.append(f"{path}: committed {a!r}, computed {b!r}")


walk(ref, new, "")
if ref.get("generator_sha256") != new["generator_sha256"]:
    print("note: tools/golden.py has changed since the committed file was written")
if bad:
    print(f"GOLDEN_MISMATCH: {len(bad)} value(s) differ from {args.out}")
    print("\n".join(bad[:20]))
    sys.exit(1)
print(f"GOLDEN_MATCH: {len(node_cases)} node, {len(fcs)} fcs, {len(modes)} mode, {len(net_cases)} network cases "
      f"agree with {Path(args.out).name} (largest difference {worst:.2e})")
