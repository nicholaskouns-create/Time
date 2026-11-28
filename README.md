
> **[Universal Recursive Identity Theorem](https://nicholaskouns-create.github.io/Identity/leri/universal-invariant/)** — Complete continuous invariant factorization, quotient identity and 125→47→5 spectral descent. Analytic proof and computational exhibit; new formal compiler status separate.

> ⚛️ **[E₄₇ Φ₄ LIVE QUANTUM CIRCUIT](https://nicholaskouns-create.github.io/qsvt-phi4/)** · [Full 125D → 250D → 500D circuit](https://nicholaskouns-create.github.io/E47-Kartekeya/interfaces/qsvt-phi4-circuit/) · [Exact complementary factorization + reproducible NumPy certificate](https://github.com/nicholaskouns-create/E47-Kartekeya/tree/main/research/e47/qsp_phi4). Shared derivation; no independent application benchmark inferred.

# E47 Quantum Clock

> **Current license:** [All Rights Reserved](LICENSE) for owned original material. No new reuse permission. Earlier licenses and third-party/public-domain rights remain unaffected.


> **Current licensing:** [All Rights Reserved](LICENSE) for original material owned by Nicholas Shane Kouns. No new reuse permission is granted. Previously granted licenses and third-party/public-domain rights remain unaffected.

Phase readout, protected clock memory and networked quantum dynamics

**[Open the Quantum Clock](https://nicholaskouns-create.github.io/Time/)**

A qubit held inside the protected 47-dimensional sector of V₂⊗V₂⊗V₂. Its hand
is the phase between the Casimir sectors C = 6 and C = 30. It reads its own
proper time, keeps that time exactly while recovery pulls every leak back, and
clocks on different worldlines drift apart by exactly (N₁ − N₂)ω. Networked,
the clocks form a Temporal Quantum Graph Neural Network (TQGNN). The same
network runs as a stochastic click process and as a deterministic flow.

| check | result |
|---|---|
| Precision certificate `MC-QCLOCK-PRECISION-20261019-001` | **48/48 PASS**, in your browser on every page load, in the test suite and in Actions |
| TQGNN validator | **36/36 PASS** |
| Canonical spacetime runner | **26/26 PASS** |
| JavaScript engine against the Python TQGNN package | parity to 10⁻¹⁴, through `data/golden.json` |
| Dial certificate `MC-E47-DIAL-20261126-001` | **12/12 PASS** in your browser and in Actions; NumPy replay **13/13** against `dial/dial_data.json` |

## What it is built from

| formalism | where |
|---|---|
| V₂⊗³, C, K = (C − 6I)(C − 30I), E47, P47, Γ = I − K²/99144, ‖Γⁿ − P‖ = (15/17)ⁿ | E47-Kartekeya canonical core |
| Kernel-aligned decoherence-free Lindblad theorem, protected readout, gap 11664 | Drive: *Clock Math — First-Principles Lindblad Proof* |
| Prime 47-cycle, torus clock, joint period 235, synchronised spatial kernel | Notion: *Prime Clock and Synchronized Kernel Dynamics* |
| z ↦ z² mod 5875 (period 220), U = Γ ⊗ T | MC-E47-MODSYNC-20261015-001 |
| dτ = ψ_C dt ⇒ d/dτ = ψ_C⁻¹ d/dt | ARR-C005 Time-Chain Identity (C0-TIME-CHAIN) |
| g = B4†(P6 − P0)B4 = diag(−1, 2, 8, 2), lapse N | MC-E47-SPACETIME-PRIMA-FACIE-26-20260924 |
| Clock-memory node, matched and linear clocks, recovery statistics, TQGNN | *Time, Protected Memory, and E47: A Prima Facie Review*, and `tqgnn/` here |
| QGNN graph Hamiltonians | Verdon et al., arXiv:1909.12264 |
| Radix dialects of one phase | Drive: *Clock-Based Mathematics* |
| City state receipts | THE MATRIX — Quantum Simulator |

## The clock

- **Node.** |a⟩ = pair singlet ⊗ |m₃ = 0⟩ (C = 6), |b⟩ = P₃₀|38⟩ (C = 30),
  syndromes |a′⟩ = P₁₂|38⟩, |b′⟩ = P₄₂|38⟩. H = f(C) + c(WV† + VW†),
  R = √γ VW†. The 125-dimensional generator restricts exactly to 4 dimensions,
  logical ⊗ flag.
- **Clock.** The linear clock f = ωC/24 ticks at ω inside E47, but only (6, 30)
  has Casimir gap 24, so leak and recovery see different detunings for |a⟩ and
  |b⟩ and dephase it at γ_φ = ½γc²|z_a − z_b|². The matched clock moves f(42)
  to f(12) + ω; then the coherence obeys κ(τ) = e^(−iωτ)κ(0) exactly, through
  every click.
- **Ticks.** Recovery clicks are sub-Poissonian, F = 1 − 24c²γ²/(γ² + 8c²)².
  Efficiency η thins them to F_η = 1 + η(F − 1). Counting them estimates proper
  time with relative error √(F_η/(ηλτ)), which fails as η → 0, where the phase
  clock takes over.
- **Matrix bridge.** |0⟩ ↦ |b⟩, |1⟩ ↦ |a⟩. One forward tick ωτ = π/4 is the
  Matrix T gate, exactly, under leak and recovery.

## The network

- **Time.** Each node moves on its own worldline in the pullback metric
  g = diag(−1, 2, 8, 2), which gives it a lapse N_v = dτ_v/dt.
- **Graph.** QGNN exchange and Ising message passing, with recovery, run in each
  node's proper time: H_t = Σ N_v H_v + Σ N_e H_e and R_{v,t} = √N_v R_v. With
  equal lapses this is the time-chain identity dρ/dt = N(t) L_τ ρ.
- **Stochastic to dynamic.** Detection efficiency η unravels the same generator.
  At η = 1 the network is a click-resolved stochastic process. At η = 0 it is
  the deterministic lapse-scaled flow. In between, the conditional variance
  rises monotonically with η.

What `validation/tqgnn_validation.py` establishes:

- **Spectral GCN.** The Casimir is a graph operator on the 125 sites, E47 is its
  spectral band {6, 30}, Γ is a 4-hop graph filter, and stacking it converges
  to P at exactly (15/17)ⁿ.
- **Proper time.** Each node's clock phase reads that node's proper time, and
  exchange between clocks on different worldlines is detuned by (N_j − N_k)ω.
- **Zeno route.** Stronger recovery makes clicks rarer and the clock
  deterministic, with γ_φ ∝ γ⁻³.
- **Learning.** Lapses are not identifiable from uncalibrated clocks (one flat
  direction per node). With calibrated clocks, QGRNN-style training recovers the
  metric lapses 1, √0.82, √0.74.
- **Prime clock.** The periods 47 and 235 are confirmed.

Derivations are in [`docs/tqgnn-formalism.md`](docs/tqgnn-formalism.md). The
status of every claim, the modelling choices and what is not established are in
[`docs/scope.md`](docs/scope.md). The dynamic formalism (clock generator, E47
contraction, temporal QGNN, click time) is in [`docs/FORMALISM.md`](docs/FORMALISM.md),
with its integrator in [`src/time_dynamics.py`](src/time_dynamics.py).

## Run

The clock needs Node 22 and a browser, nothing else:

```
node --test tests/engine.test.js        # 10 tests, includes the certificate
node tools/qclock-cli.js certificate    # rewrites data/certificate.json
node tools/qclock-cli.js sweep --param gamma --from 2 --to 60 --steps 30 --csv
node tools/qclock-cli.js run --tau 50 --clock linear --gamma 10
node tools/qclock-cli.js verify-receipt matrix-quantum-receipt.json
node audit/site_matrix_audit.js         # Matrix precision audit, on its current source
python -m http.server                   # then open http://localhost:8000/
```

The TQGNN package needs Python 3.10+:

```
pip install -r requirements.txt
python validation/e47_spacetime_prima_facie.py   # 26 checks
python validation/tqgnn_validation.py            # 36 checks, about a minute
python -m unittest tests.test_time_dynamics      # src/time_dynamics.py
python tools/golden.py --check                   # the clock's reference values, recomputed
python tools/golden.py                           # rewrite data/golden.json
```

GitHub Actions runs both sets on every push and pull request.

## Layout

```
index.html, app.js, ui.js     the Quantum Clock page
worker.js                     the certificate, off the main thread
engine/linalg.js              dense real and complex linear algebra, Padé-13 expm, Jacobi eigen
engine/carrier.js             V2⊗3, Casimir graph, projector, pullback metric, lapse
engine/node.js                the clock-memory node, Liouvillian, counting statistics, closed forms
engine/sim.js                 the running clock (exact propagators, monitored runs)
engine/network.js             TQGNN networks, Taylor propagation, trajectories
engine/clocks.js              prime clock, torus clock, synchronisation clock, radix dialects
engine/receipt.js             Matrix receipts: digest, replay in both orders, clock-tick export
engine/certificate.js         the 48 checks
engine/automate.js            sweeps, CSV, permalinks
tools/                        qclock-cli.js (the engine, headless), golden.py (reference values)
data/golden.json              reference values from the Python TQGNN package (exact expm propagation)
data/certificate.json         the committed certificate
dial/                         A Dial of Forty-Seven Hours: page, certificate, NumPy replay, data, plates
audit/                        precision audit of both Matrix implementations
tqgnn/                        the TQGNN in Python: carrier, spacetime, node, network,
                              stochastic and jump unravellings, Lindblad, clocks, learning
validation/                   canonical spacetime runner (unchanged) and TQGNN validator
src/time_dynamics.py          Lindblad generator, projected Hamiltonian, RK4 step
tests/                        engine.test.js (the clock), test_time_dynamics.py (src/)
docs/FORMALISM.md             the dynamic formalism
docs/tqgnn-formalism.md       TQGNN derivations
docs/scope.md                 status of every claim
results/                      last validator output
```

## A Dial of Forty-Seven Hours

**[Open the dial](https://nicholaskouns-create.github.io/Time/dial/)** · [source](dial/)

The grading S = (18I − C)P/12 = −XP/6 splits E47 into 25 states of total spin 2 and 22 of
total spin 5, with Tr S = 3. Every pairing of the two sides leaves at least three states
unpaired, no nonzero pairing commutes with SU(2), and a clock that commutes with SU(2) shows
at most seven distinct hours. The clock's own hand is this grading: f(C)P = ω(¾P − ½S), and
one 1/47 turn is U₀ = ζI₂₅ ⊕ I₂₂ up to a global phase. The prime 47-cycle shows all 47
hours, so it cannot commute with SU(2) in any basis.

```
node --test tests/dial.test.js     # 12 checks, on this engine
python dial/verify_dial.py          # independent NumPy replay against dial/dial_data.json
```

## Credits

E47 framework and kernel: Nicholas Shane Kouns, AIMS Research Institute. Clock
construction, engine, instrument and TQGNN: Claude (Anthropic), with Nicholas
Shane Kouns, 19 October 2026. Quantum clock construction with phase-based elapsed-time readout, validated in mathematics and numerical quantum simulation. Laboratory hardware realization and calibration to physical seconds remain open.

[All Rights Reserved](LICENSE)


## E47 spectral and noiseless-subsystem certificate (2026-10-24)

[Read the cross-repository evidence note](E47_SPECTRAL_ADVANTAGES_2026-10-24.md). Documents the exact 125→47 Casimir construction, the (15/17)^n spectral convergence certificate, and collective-SU(2) noiseless multiplicity encoding, separately from channel-dependent QEC and unmeasured hardware advantages.

[The E₄₇ spectral kernel · formalism and verification](e47-spectral-kernel/README.md)

[Spin, weight blocks and the matched-clock illustration](research/e47/clock-proofs/README.md)


<!-- FUN-WITH-EINSTEIN-2026 -->
### 🪐 Fun With Einstein · Interactive spectral & quantum lab

[Launch Fun With Einstein](https://nicholaskouns-create.github.io/fun-with-einstein/) · [Source & Python/Lean/Coq proof modules](https://github.com/nicholaskouns-create/E47-Kartekeya/tree/main/docs/fun-with-einstein) · [CI compiler status](https://github.com/nicholaskouns-create/E47-Kartekeya/actions/workflows/fun-with-einstein.yml).

Fivefold C₅ eigenvalue ratio = φ⁻²; 125D Casimir E₄₇ projector; 250D quantum-ancilla probability α=47/(125φ²); stationary Einstein-equation recovery. Spectral simulation and gravitational statements have different evidence classes.
<!-- /FUN-WITH-EINSTEIN-2026 -->


---

### Interactive E₄₇ spectral geometry (computed)

[**Open the standalone E₄₇ 3D visualizer**](./e47-geometry/index.html) · [**Open on GitHub Pages**](https://nicholaskouns-create.github.io/Time/e47-geometry/)

An additive, mobile-friendly research instrument showing the 125-site product-spin lattice, the exact 47-dimensional projector-density field (125 rational samples with denominator 315), five orthogonal axial slices, the seven SU(2) spectral sectors, optimal contraction Γ* = I − K²/99144 and the five-stage exact polynomial projector. The volumetric embedding is a coordinate visualization, not a claim that E₄₇ is a physical three-dimensional solid. No existing instruments are replaced.


<!-- ΩBA-PROOF-LINK-20261124 -->
### ΩBA · First-Principles Triadic Contraction and E₄₇ Continuity

[**Open the light-mode interactive ΩBA proof plate**](https://nicholaskouns-create.github.io/omega-ba/) · [Local Unicode derivation](omega-ba/README.md) · [Executable NumPy validation](omega-ba/omega_ba_validation.py) · [Local interactive page](omega-ba/index.html).

A mathematical three-stage alignment → binding → commitment map (𝒪=A∘B∘Ω), Banach strict-contraction proof with explicit assumptions, and a distinct E₄₇ spectral complement contraction (125→47, rate 15/17). Python tests and evidence boundaries are supplied; acoustic, biological and retrocausal effects are not asserted.
<!-- /ΩBA-PROOF-LINK-20261124 -->
