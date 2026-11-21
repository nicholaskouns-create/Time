# E₄₇ Φ₄ · complete QSP / LCU circuit cross-reference

Canonical shared result (E0 exact algebra, E1 locally executed NumPy certificate, 14/14 PASS):

Q₄(s) = 1 − 2176s²/353 + 2048s⁴/353

P(s) = Q₄(s) + 272is²/353

R(s) = s(2048s² − 1152 + 16i√817)/353

|P(s)|² + (1−s²)|R(s)|² = 1.

Φ₄ ≈ (0, 0.618892044578, 0.439842582816, −0.179049461762, 0) radians.

A = |(C−6I)(C−30I)|/432 on V₂⊗³ (125D).

U₂₅₀ = Z(φ₀)W(A)Z(φ₁)W(A)Z(φ₂)W(A)Z(φ₃)W(A)Z(φ₄), with W(A) = [[A,i√(I−A²)],[i√(I−A²),A]].

V₅₀₀ = (H⊗I₂₅₀)diag(U₂₅₀,U₂₅₀†)(H⊗I₂₅₀).

⟨00|V₅₀₀|00⟩ = Q₄(A), and ‖Q₄(A)−P₄₇‖₂ = 225/353.

The principal block is postselected; no hardware run, physical resource-count proof, or independent current-repository matrix replay is asserted.

## Source and run

- [Interactive full 250D and 500D circuit](https://nicholaskouns-create.github.io/E47-Kartekeya/interfaces/qsvt-phi4-circuit/)
- [Full-screen embedded circuit](https://nicholaskouns-create.github.io/qsvt-phi4/)
- [First-principles proof](https://github.com/nicholaskouns-create/E47-Kartekeya/blob/main/research/e47/qsp_phi4/README.md)
- [Full product-basis CG Python certificate](https://github.com/nicholaskouns-create/E47-Kartekeya/blob/main/research/e47/qsp_phi4/e47_phi4_cg_lcu_certificate.py)
- [JSON numeric receipt](https://github.com/nicholaskouns-create/E47-Kartekeya/blob/main/research/e47/qsp_phi4/e47_phi4_cg_lcu_receipt.json)

**Scope for this repository:** Clock error/recovery performance requires separate temporal noise and readout tests.
