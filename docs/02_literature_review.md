# 02 — Literature review (as of 2026-10-08)

Compiled from three web searches run on 2026-10-08. Items marked *(unverified)* were seen only second-hand.

## The direct predecessor

**Hu, Kaneko, Martin (Spencer), Schildkraut**, "On a Mertens-type conjecture for number fields", MPCPS 177 (2024)
481–523; arXiv:2109.06665 (v2, 16 Mar 2023; local copy in `references/`).
- Thm 1.1: naive conjecture false for all imaginary quadratic fields except Q(√−3), all real quadratic except Q(√5).
- Thm 1.2: false for all fields of fixed signature with large discriminant; Cor. 2.11: false for non-abelian Galois K.
- Explicit witnesses: Q(i) h*_{600}(72.85) ≤ −1.008, h*_{600}(−85.15) ≥ 1.029; Q(√3) h*_{200}(17.32) > 1.027;
  Q(√2) h*_{150}(−24.64) > 1.049. **We could not reproduce three of these four** (see 03_findings_log.md).
- Open questions §6–7: Conj. 6.1 (finitely many sign changes for infinitely many K), Q. 6.2, Q. 6.3, Q. 7.1–7.4
  (limiting distributions, ν_K([−1,1]), tail bounds).

## Competition check for Q(√−3), Q(√5)

- Citing works (OpenAlex, Google Scholar): one — the annotated bibliography for comparative prime number theory
  (Martin, Yang et al., Expo. Math. 2025) *(entry not read)*. Semantic Scholar was rate-limited.
- No paper found settling either field, or applying Odlyzko–te Riele / LLL over number fields.
- Related but different: Simonič, arXiv:2201.10771 (GRH upper bounds for Mertens-type functions incl. Dedekind zeta);
  Kim–Nguyen, Res. Number Theory 11:3 (2025), arXiv:2502.21021, and Rozmarynowycz–Kim, arXiv:2305.00345 (both over Q).

## Deep-research prior-work check (2026-10-10)

Full report: `05_prior_work_deep_research_2026-10-10.json` (98 agents, 3-vote adversarial verification).
**No source resolves Q(√−3) or Q(√5).** Verified (3-0 votes):
- HKMS leave both open (Thm 1.1; Thms 3.7/3.10; remark suggesting OtR85 computations) — arXiv and Cambridge Core text.
- The annotated bibliography (Expo. Math. 43 (2025) 125644, arXiv:2309.08729v3, entry [356]/[A108]) repeats the
  exclusion, records no resolution, and no other entry builds on HKMS; it claims completeness only through 2023.
- CPNT Problem List (Hamieh–Kadiri–Martin–Ng, June 2024): nothing on M_K(x)/√x; HKMS not cited.
- Rozmarynowycz–Kim (2305.00345), Kim–Nguyen (2502.21021), Johnston–Trudgian (2511.09978 v4, Jul 2026) work over
  Q only and do not cite HKMS; Humphries (function fields), Mossinghoff–Trudgian (over Z) do not apply.
Remaining gaps: no full citation-index sweep (Google Scholar/zbMATH/MathSciNet); 2025–26 talks/theses not excluded. Note: Omega-results for M(x, χ) would not
settle ζ_K = ζ·L(s,χ), since μ_K = μ * (μχ) is a different function.

## Records over Q (context)

| Function | Best proven | Source |
|---|---|---|
| M(x)/√x | limsup > 1.826054, liminf < −1.837625 | Hurst, Math. Comp. 87 (2018) 1013–1028, arXiv:1610.08551 |
| earlier | 1.06/−1.009 (OtR 1985); 1.218/−1.229 (Kotnik–te Riele 2006); ±1.6383 (Best–Trudgian, Math. Comp. 2015, arXiv:1209.3843) | |
| first counterexample | < exp(1.96·10^19) | Kim–Nguyen 2025 |
| L(x)/√x (Pólya) | limsup ≥ 1 | Mossinghoff–Trudgian, in *Exploring the Riemann Zeta Function* (Springer 2017) |
| Σλ(n)/n (Turán) | √x·T(x) > 2.37 i.o. | Mossinghoff–Trudgian, IJNT 2021, arXiv:2007.14537 |
| Σμ(n)/n | no numerical Ω-record found | — |

Hurst's cost: N = 800 zeros, LLL (fplll) 2 × 35 days, zeros to ~8151 digits. Under the linear independence hypothesis
limsup = +∞ (Ingham 1942), so constants can only be improved, never made unbounded, by this method.

## Other gaps identified (possible follow-ups)

1. Explicit Ω-constants for M(x; q, a) and Σ μ(n)χ(n) (q = 3, 4, 5): none found.
2. Explicit constants for "fake μ's" (Martin–Mossinghoff–Trudgian PAMS 2023, arXiv:2112.05227; Martin–Yip CJM 2025,
   arXiv:2411.06610 — constants left implicit).
3. Σ μ(n)/n · √x: same zero amplitudes as M(x), shifted phases; no record found.

## Zero data sources

- LMFDB ζ zeros (Platt): first ~1.04·10^11 zeros, ±2^-102.
- Odlyzko tables: first 2·10^6 zeros to ~4·10^-9; first 100 to 1000+ digits.
- Dirichlet L-zeros: LMFDB pages give ~50 zeros per character; this project computes its own with Arb.
