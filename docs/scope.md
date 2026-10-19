# Scope and evidence

Evidence classes follow the E47 convention: **E0** exact proof, **E1**
executable reconstruction, **E2** simulation, **E3** observation, **E4**
experiment or hardware, **OPEN** unresolved. Check IDs refer to
`validation/tqgnn_validation.py` unless marked *canonical runner*.

## Established

| Claim | Class | Check |
|---|---|---|
| Casimir spectrum, rank P = 47 | E0 | A1, A2 |
| C is a graph operator on the 125 sites; 13 planar components | E0 / E1 | A3, A4 |
| Gamma = I - K^2/99144 is 4-local on that graph; ‖Γⁿ − P‖₂ = (15/17)ⁿ | E0 | A5, A6 |
| Only (6, 30) has Casimir gap 24 | E0 | A7 |
| g = B4†(P6 − P0)B4 = diag(−1, 2, 8, 2) and the curvature chain | E1 | B1, *canonical runner 26/26* |
| Lapse N = sqrt(1 − 2v1² − 8v2² − 2v3²) for timelike worldlines | E0 | B2 |
| Time-chain identity dρ/dt = N(t) L_τ ρ | E0 | B3, D3 |
| Frame direction u lies outside E47 (Γu = 0.6732026144 u) | E1 | B4 |
| 125-dim node dynamics reduces exactly to 4 dims | E0 | C1, C2 |
| Clock ⟨O⟩ = cos ωτ inside E47 | E0 | C3 |
| Matched recovery retention and click rates (review table) | E0 | C4 |
| Stationary leak 4c²/(γ² + 8c² + 4Δ²) per branch | E0 | C5 |
| Fano factor of recovery clicks 1 − 24c²γ²/(γ² + 8c²)² | E1 | C6 |
| Linear clock + recovery dephases the clock (second-order rate) | E1 | C7 |
| Equal detunings (spectral clock) preserve the clock exactly | E0 | C8 |
| Network layers are CPTP for any non-negative lapses | E0 | D2 |
| Each node's clock reads its own proper time | E0 | D1 |
| Exchange between nodes is detuned by (N_j − N_k)ω | E0 | D4 |
| η = 0 conditional state equals the GKSL solution | E0 | S1 |
| Unravelling invariance E[ρ_c] = ρ | E0 (proof), E2 (Monte Carlo) | S2, S5 |
| Detected click mean and thinning law F_η = 1 + η(F − 1) | E0 | S3, S6 |
| Conditional variance non-decreasing in η | E0 (proof), E2 | S4 |
| Click-time precision sqrt(F_η/(ηλτ)) | E2 | S7 |
| Zeno route: γ_φ ∝ γ⁻³, click rate falls with γ | E1 | S8 |
| Lapses unidentifiable from uncalibrated clocks; recovered with calibrated clocks | E0 / E2 | L1–L3 |
| Prime clock periods 47 and 235 | E0 | P1, P2 |

## Modelling choices (definitions, not results)

- Clock and syndrome states: |a⟩ is the canonical runner's C = 6 copy with
  m₃ = 0; |b⟩, |a′⟩, |b′⟩ project the site m = (1, 0, −1) onto C = 30, 12, 42.
  Results C4–C8 depend on the syndrome sectors only through Δ_a, Δ_b.
- Clocks: the review's linear clock ω C/24, and a spectral clock that moves
  the C = 42 level so the detunings match.
- Edges act on the logical registers regardless of the flag (exchange and
  Ising terms).
- Edge lapse N_e = (N_j + N_k)/2. The geometric mean is available; every
  non-negative rule keeps the generator GKSL.
- Worldlines are curves in the fixed background metric; nothing back-reacts
  on the geometry.
- Units are dimensionless (ħ = 1; ω = 2 in validation).

## Not established

- **OPEN.** That the E47 selection rule determines physical time. The
  metric's indefinite form and the geometric embedding are additional
  choices.
- **OPEN.** Physical realization of the nodes, the edge couplings and the
  recovery reservoir.
- **OPEN.** Calibration of any clock in seconds: phase clock, click clock or
  prime cycle.
- **E3, E4.** No observation, experiment or hardware result is claimed.

## Numerical and statistical checks

Deterministic checks use DOP853 with rtol 1e-11 and atol 1e-13, or exact
diagonalisation. Monte Carlo checks use fixed seeds. They compare sample
means with exact GKSL values in units of the sample standard error and pass
below 4.5. The click-precision check passes within 10 % of the asymptotic
formula.
