# 04 — Certification plan

Goal: a computer-assisted proof, in Arb ball arithmetic, that h*_{K,T}(t) > 1 or < −1 for an explicit (T, t), for
K = Q(√−3) and K = Q(√5). By 01_problem_statement.md this disproves the naive conjecture for K unconditionally.

All steps are implemented in `src/certify/`, run through `src/certify/run_certificate.py`, and write a JSON certificate
plus a log to `certificates/`.

## Step 1 — Rigorous zero count N(T) (argument principle) — `count_zeros.py`

For a real primitive character χ mod q with parity a (root number 1),
Λ(s) = (q/π)^{(s+a)/2} Γ((s+a)/2) L(s,χ) is entire, Λ(s) = Λ(1−s), Λ(s̄) = conj Λ(s). For T not a zero ordinate,

  N(T,χ) := #{ρ : 0 < γ < T} = (1/π) [ θ(T) + Δ_C arg L(s,χ) ],   θ(T) = (T/2) log(q/π) + Im log Γ((1/2+a+iT)/2),

where C is the path 2 → 2+iT → 1/2+iT and arg L(2) = 0. Requires: L(σ,χ) ≠ 0 for σ ∈ [1/2, 2] (no zeros on the bottom
edge; certified by bisection with ball enclosures) and L(1/2+iT,χ) ≠ 0.
For ζ (q = 1) use ξ(s) = ½ s(s−1) π^{−s/2} Γ(s/2) ζ(s): the factor s(s−1) contributes +π, so N(T) = (θ(T)+Δ_C arg ζ)/π + 1;
ζ(σ) < 0 on (0,1) classically (ζ = η/(1−2^{1−σ}), η(σ) > 0).

- Vertical segment σ = 2: Re L(2+it) ≥ 2 − ζ(2) > 0, so Δ arg = principal arg L(2+iT) exactly.
- Horizontal segment: adaptive bisection; on each piece the ball enclosure of L over the whole piece must exclude 0
  (then the arg variation on the piece is < π and equals the principal arg of L(end)/L(start)).
- Output: an integer, certified (the enclosure of the right-hand side contains exactly one integer).

## Step 2 — Isolation of every zero on the critical line — `isolate_zeros.py`

For each approximate zero γ (from data/zeros), certify that the Hardy Z-function has strictly opposite signs at γ ± ε
(ε = 10^-28, balls excluding 0). Each interval then contains a zero of odd multiplicity on the critical line.
If the number of disjoint intervals in (0, T) equals N(T) from Step 1, then: every zero with 0 < γ < T is on the
critical line, is simple, and lies in exactly one of the intervals. Done separately for ζ and for L(s,χ); then check
that no ζ-interval meets an L-interval (so ζ_K has only simple zeros below T).

## Step 3 — Enclosure of h*_{K,T}(t) — `certify_hstar.py`

For each zero, ρ is the ball 1/2 + i[γ−ε, γ+ε]. ζ_K′(ρ) is enclosed using power-series (Taylor) arithmetic
(`acb_series.dirichlet_l`), which encloses the derivative over the whole ball, hence at the true zero:
- ζ-zero: ζ_K′(ρ) = ζ′(ρ) L(ρ,χ);  L-zero: ζ_K′(ρ) = ζ(ρ) L′(ρ,χ).
Then h* = Σ 2 f(γ/T) Re( e^{iγt} / (ρ ζ_K′(ρ)) ) in ball arithmetic. Certified if the final ball lies in (1, ∞) or (−∞, −1).

## Step 4 — Independent checks (not part of the proof, but reported)

- ζ zeros also compared with Arb's `acb.zeta_zeros` (itself rigorous).
- h* recomputed in mpmath from midpoints (hiprec_check.py).
- Sensitivity: certificate re-run at a second working precision gives overlapping enclosures.

## Witnesses to certify

| Field | T | t | numerical h* |
|---|---|---|---|
| Q(√−3) | 200 | 1560.51100774 | −1.2755 |
| Q(√−3) | 200 | 18558.8477306 | +1.2079 |
| Q(√5) | 1000 | 1016.24860331 | +1.1302 |
| Q(√5) | 1000 | 722.734868673 | −1.1108 |

## Known gaps / to discuss in the paper
- The step from "h*_{K,T}(t) > 1" to "M_K^+ > 1" is the paper's argument (Jurkat; their Prop. 3.3 and §3.2 end); we rely
  on it as published and should re-derive it in our write-up.
- Arb/FLINT correctness is assumed (standard for computer-assisted proofs); record versions in each certificate.
