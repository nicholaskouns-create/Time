# Quantum Clock: dynamic formalism

## Clock generator
For a finite-dimensional clock, take a Hermitian Hamiltonian H and density operator ρ. Define physical-time evolution by

dρ/dt = -i[H,ρ] + Σ_j γ_j(L_j ρ L_j† - ½{L_j†L_j,ρ}).

The Hamiltonian contribution is deterministic and trace preserving. The Lindblad terms define a deterministic evolution of the ensemble density matrix even when a stochastic quantum-trajectory unraveling is available.

## E47 contraction
On H₁₂₅ = V₂⊗V₂⊗V₂, set C = J_tot², K = (C-6I)(C-30I), P = 1_{C∈{6,30}}. Then rank(P)=47.

Γ = I-K²/99144, Γⁿ→P and ‖Γⁿ-P‖₂=(15/17)ⁿ.

For continuous contraction S(t)=exp(-tK²), S(t)→P as t→∞.

These are spectral selection operators; their contraction parameter is not automatically laboratory clock time.

## Temporal QGNN
Let G(t)=(V,E(t)) be a time-varying graph with Hermitian coupling H_G(t), and let ρ(t) be its encoded density state. Define

dρ/dt = -i[H_G(t),ρ] + D_t(ρ).

For exact preservation of E47, require H_G(t)=PH_G(t)P + (I-P)H_G(t)(I-P) and each Lindblad operator L_j(t) to preserve the P block (e.g. [L_j(t),P]=0). An initially P-supported ρ remains P-supported.

A projected generator H_P(t)=PH_G(t)P gives a direct invariant implementation. Time ordering is required for noncommuting H_P(t).

## Click time
A click count N(t) is an event observable. A specified counting process or detector model is needed to map click index to physical time; an integer count alone does not define a duration.

## Evidence
The E47 rank and contraction identities are established algebraic claims. The time-dependent graph evolution above is a mathematically specified construction, not a claim of experimental implementation.
