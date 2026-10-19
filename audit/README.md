# The Matrix Quantum Simulator: precision audit

`MC-MATRIX-PRECISION-AUDIT-20261019-001` · 19 October 2026 · evidence class E1

The Quantum Clock imports and exports Matrix state receipts, so both Matrix
implementations were checked by running their own code, unmodified. Every
finding below is reproduced by the scripts in this folder; the numbers are in
[`MC-MATRIX-PRECISION-AUDIT-20261019-001.json`](MC-MATRIX-PRECISION-AUDIT-20261019-001.json).

```
node audit/site_matrix_audit.js [app folder]   # site app, nicholaskouns-create.github.io
node audit/city_matrix_audit.js [worker]       # City MATRIX, E47-Kartekeya
```

Given no path, each script fetches the app's current source from GitHub, so the
site audit's F1–F5 read `resolved` as soon as a fix is merged there. It also
reports the SHA-256 of the files it ran, to compare with the ones pinned in the
JSON record.

## The site app (`matrix-quantum-simulator/`)

| id | finding | evidence | fix |
|---|---|---|---|
| F1 | The importer rejects every receipt the app exports with at least one gate. | `valid()` checks `o.g`; `receipt()` writes `o.gate`. Bell, GHZ and superposition receipts: "Invalid circuit". Only an empty circuit imports. | In `receipt.js` `valid()`, read `o.gate` (two places). |
| F2 | `run()` applies gates in the order they were placed, not the column order the circuit shows. | Re-placing H in the Bell preset turns (\|000⟩+\|110⟩)/√2 into (\|000⟩+\|100⟩)/√2 with the same diagram. H at column 0 placed after S at column 1 computes S then H. | Iterate over `[...ops].sort((x, y) => x.col - y.col)`. |
| F3 | Placing a two-qubit gate deletes every gate in its column, including on the third wire. | X on q2 at column 3 disappears when CX q0→q1 is placed at column 3. | Remove only gates touching the two targets. |
| F4 | Placing a gate does not recompute the state, so EXPORT can pair a circuit with an older statevector. | Circuit [H on q0], statevector \|000⟩. | Call `run()` at the start of `receipt()`. |
| F5 | MEASURE falls through to \|111⟩ when round-off leaves the cumulative probability below the draw. | The loop stops at k = 7 regardless of its probability (≈ 10⁻¹⁶ chance). | Stop at the last state with nonzero probability. |
| F6 | The SHA-256 digest detects accidental edits, not deliberate ones. | Edited amplitudes with a recomputed digest import cleanly. | Replay the circuit on import. |

## The City MATRIX two-site MPS worker (`E47-Kartekeya/website/interfaces/matrix/`)

The certified untruncated run (8 qubits, 24 layers, χ = 16) is exact and is not
affected. Truncated runs are:

| χ | worker fidelity | canonical TEBD fidelity, same χ | truncating steps worse than the Schmidt optimum | worst error ratio | actual ÷ reported norm loss |
|---|---|---|---|---|---|
| 8 | 0.946930 | 0.994436 | 8 of 8 | 54.3× | 0.52 – 1.09 |
| 4 | 0.811701 | 0.944280 | 27 of 27 | 15.3× | 0.40 – 1.73 |

| id | finding | fix |
|---|---|---|
| M1 | Two-site updates are not made at the orthogonality centre, so the SVD weights are not the state's Schmidt coefficients: every truncation is suboptimal, and the reported discarded weight is not the norm actually lost. At χ = 8 the final norm² 0.92537 lies outside the window [1 − Σε, e^(−Σε)] = [0.90225, 0.90688] that canonical truncation guarantees. The χ = 4 stress-test figures measure this gauge choice, not the best χ = 4 MPS. | Keep Vidal Γλ form, or QR-sweep the orthogonality centre to the active bond before each SVD. |
| M2 | The circuit drawing shows a T1/Tφ noise box after every layer; the worker applies the channel once, after the last layer. | Apply it per layer, or draw it once. |
| M3 | In controlled-truncation mode the validation passes whenever the error is finite. | After M1, compare the fidelity with Π(1 − ε_k). |
| M4 | Only the brickwork entangler has a dense reference; the chain entangler is never validated. | Add the chain ordering to `exact()`. |
| M5 | `hermitianEig` is defined but unused. | Remove or use. |

## What the Quantum Clock does about it

- Replays every imported receipt in diagram order and in placement order and
  says which one the statevector matches (F2, F4, F6).
- Writes each exported gate as both `gate` and `g`, so the deployed importer
  accepts it (F1). Verified end to end: a GHZ receipt exported by the Matrix,
  seven clock ticks on q2, re-imported by the Matrix ("Imported and SHA-256
  verified"), and the Matrix's RUN reproduces the state to 0.
- Appends ticks after the last used column, where both orders agree (F2, F3).

Software audit of browser simulators only; no hardware claim.
