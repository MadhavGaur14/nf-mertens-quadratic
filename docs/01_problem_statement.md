# 01 — Problem statement and logic of the disproof

## Definitions (following Hu–Kaneko–Martin–Schildkraut, arXiv:2109.06665, §1–3)

- μ_K: Möbius function on integral ideals of K. M_K(x) = Σ_{N(𝔞)≤x} μ_K(𝔞), with the term N(𝔞) = x halved when x is an integer.
- M_K^± = limsup / liminf of M_K(x)/x^{1/2}.
- **Naive Mertens-type conjecture over K:** −1 ≤ M_K^− ≤ M_K^+ ≤ 1.
- For quadratic K of discriminant D: ζ_K(s) = ζ(s) L(s, χ_D), χ_D the Kronecker symbol.
  - Q(√−3): D = −3, χ = χ_{−3} (odd), python-flint `dirichlet_char(3, 2)`.
  - Q(√5): D = 5, χ = χ_5 (even), `dirichlet_char(5, 4)`.

## The test quantity (paper, eq. (3.6))

h*_{K,T}(t) = Σ_{|γ|<T} f(|γ|/T) e^{iγt} / (ρ ζ_K′(ρ)),  f(y) = (1−y)cos πy + (1/π) sin πy (Jurkat–Peyerimhoff kernel),

summed over the nontrivial zeros ρ = 1/2 + iγ of ζ_K. Since χ_D is real, zeros come in conjugate pairs and
h* = 2 Re Σ_{0<γ<T} f(γ/T) e^{iγt}/(ρ ζ_K′(ρ)).

At a zero of ζ: ζ_K′(ρ) = ζ′(ρ) L(ρ, χ). At a zero of L: ζ_K′(ρ) = ζ(ρ) L′(ρ, χ).

## Logic (paper, Lemma 3.4 proof and §3.2 end)

1. If RH fails for ζ_K, or ζ_K has a multiple zero, then M_K^+ = ∞ or M_K^− = −∞ (paper, Thm 2.10): conjecture false.
2. Otherwise (RH + simple zeros), if 1/|ρ ζ_K′(ρ)| is unbounded the same conclusion holds (Prop. 2.8); else Jurkat's
   theorem applies and M_K^− ≤ liminf_t h*_{K,T}(t) ≤ h*_{K,T}(t) ≤ limsup_t h*_{K,T}(t) ≤ M_K^+ for every T > 0, t ∈ ℝ.
3. Hence **one pair (T, t) with h*_{K,T}(t) > 1 (or < −1) disproves the conjecture for K, unconditionally.**

For the computed value to *be* h*_{K,T}(t) we must know, rigorously, every zero of ζ_K with |γ| < T. In case 2 those
are all on the critical line and simple; if our certification finds otherwise, case 1 applies and the conjecture is false
anyway. So the certificate must establish: (a) the complete list of zeros of ζ_K with 0 < γ < T, each on the line and
simple; (b) enclosures of ρ and ζ_K′(ρ); (c) an enclosure of h*_{K,T}(t) lying entirely above 1 (or below −1).

## Starting point from the paper

M*_K(1) ≈ −0.4851 (Q(√−3)) and −0.4857 (Q(√5)) give M_K^+ ≥ 0.515 and M_K^− ≤ −0.485: both fields start about halfway
to the threshold. The paper's Remark after Theorem 1.1: "It is, however, likely that computations analogous to those used
in [OtR85] to falsify the original Mertens conjecture would suffice for both of these quadratic fields as well."
