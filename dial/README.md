# A Dial of Forty-Seven Hours

**[Open the dial](https://nicholaskouns-create.github.io/Time/dial/)** · Symmetry, pairing, and the three that remain · 26 November 2026

The 47 protected states of the Quantum Clock split into 25 and 22. Every pairing of the two
sides leaves at least three states unpaired, no nonzero pairing commutes with SU(2), and a
clock that commutes with SU(2) shows at most seven of its forty-seven hours.

| file | what it is |
|---|---|
| `index.html`, `dial.js` | the page: the interactive dial (grading, pairing, symmetric dial, U₀, prime clock) and the live certificate |
| `dial-certificate.js` | `MC-E47-DIAL-20261126-001`: 12 checks on the Quantum Clock's own engine (`../engine`), in the browser and in `tests/dial.test.js` |
| `verify_dial.py` | independent NumPy replay of the same claims, plus the commutant M₅ ⊕ M₂ of dimension 29; checks `dial_data.json` |
| `dial_data.json` | the verified numbers: sectors, moments, grading, pairing, symmetry, clocks, tower |
| `plates/` | the proof plate and the E47 Machine Certificate (Rev. B), light house style: full resolution PNG, web preview JPG |

## The proof

1. **Carrier.** j = 2, c = j(j + 1) = 6, ℋ = V₂ ⊗ V₂ ⊗ V₂, dim ℋ = 125.
2. **Exchange and selector.** X = J₁·J₂ + J₂·J₃ + J₃·J₁, C = 18I + 2X, μ = 18, σ² = 144,
   K = (C − μI)² − σ²I = (C − 6I)(C − 30I) = 4(X² − 36I). The moment rule is chosen; the sectors it lands on are forced.
3. **Kernel.** ℋ₊ = ker(C − 6I) ≅ ℂ⁵ ⊗ V₂ (25), ℋ₋ = ker(C − 30I) ≅ ℂ² ⊗ V₅ (22), E₄₇ = ℋ₊ ⊕ ℋ₋.
4. **Grading.** S = (18I − C)P/12 = −XP/6, S² = P, P± = (P ± S)/2, Tr S = 25 − 22 = 3.
5. **Pairing and index.** {S, D} = 0 ⇔ D(u, v) = (A†v, Au); ind A = 3; dim ker D = 47 − 2r ≥ 3;
   AA† = I₂₂ gives spec D = {−1, 0, +1} with multiplicities 22, 3, 22.
6. **Symmetry forbids pairing.** [D, ρ(g)] = 0 for all g ⇒ [D, S] = 0; with {D, S} = 0, D = DS² = 0.
7. **Symmetric dial.** [U, ρ(g)] = 0 for all g ⇔ U = (B ⊗ I₅) ⊕ (F ⊗ I₁₁): at most 5 + 2 = 7 hours.
8. **Two clocks.** U₀ = ζI₂₅ ⊕ I₂₂ (ζ = e^(2πi/47)) has order 47 and commutes with SU(2);
   U₁ = diag(1, ζ, …, ζ⁴⁶) keeps E₄₇ but shows 47 > 7 hours, so it cannot commute with SU(2). ∎

## How it meets the Quantum Clock

- **The hand is the grading.** On E₄₇ the linear clock f = ωC/24 is f(C)P = ω(¾P − ½S); one 1/47
  turn is U₀ up to a global phase. The Quantum Clock is a symmetric dial with two hours (check D10).
- **The prime 47-cycle** (`engine/clocks.js`) shows all 47 hours, so it cannot commute with SU(2) in
  any basis, and it crosses the grading (check D11).
- **The tower.** On every Fibonacci–Lucas carrier the grading's index is F_k² − F_{k−2}L_k =
  F_{k−3}F_{k−1}: 0, 3, 24, 168, 1155, 7920 (check D12).

## Evidence

- **E0 exact:** the grading and the index (rank–nullity), the symmetry obstruction and the seven-hour bound (Schur's lemma), the tower identity (integer arithmetic).
- **E1 executable:** the 12 engine checks (residuals near 10⁻¹⁴) and the NumPy replay, including the commutant of dimension 29.
- **Not claimed:** the arcs draw one pairing in a basis adapted to it; "hours" are distinct eigenvalues of a unitary on E₄₇; no physical identification of the dial is claimed.

```
node --test tests/dial.test.js     # 12 checks, on the clock's engine
python dial/verify_dial.py          # NumPy replay against dial_data.json
python dial/verify_dial.py --write  # rewrite dial_data.json
```

Proof: Nicholas Shane Kouns. Plates, instrument and certificate with Claude (Anthropic), 26 November 2026.
