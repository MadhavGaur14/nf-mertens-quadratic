# 03 — Findings log

Newest first. "Numerical" = floating/multiprecision but not interval-certified. "Certified" = Arb ball arithmetic with
a complete proof chain (see 04_certification_plan.md).

---


---

## 2026-10-10 — Independent verification in PARI/GP 2.19.0

Code: `src/independent/pari_check*.gp`; logs: `results/2026-10-10_pari_independent/`. No code shared with `src/`:
zeta_K built from the field polynomial (`lfuncreate(x^2+x+1)` etc.), zeros from `lfunzeros`, zeta_K'(rho) from
`lfun(L, rho, 1)`.

- **Pitfall found:** `lfunzeros(zeta_K, T)` at its default subdivision **misses close zeta/L zero pairs**
  (Q(√−3), T=200: 185 instead of 193; missed pairs 37.5518/37.5862, 69.5130/69.5464, 138.0954/138.1160,
  182.1921/182.2071, gaps 0.015–0.034; similar for Q(i), Q(√2), Q(√3)). With `divz=64`, or with zeta and L(chi)
  separately, PARI finds exactly our certified counts (193, 814, 154, 237, 1553). This is an argument for the
  rigorous argument-principle count in our certificates.
- **With complete zero lists, every value agrees with our Arb certificates to all 15 printed digits:** all four
  first witnesses (Q(√−3) ±1.2755/1.2079, Q(√5) 1.1302/−1.1108), all six replacement witnesses, and the four HKMS
  witnesses (−0.93672, 1.20022, −0.79792, −0.30039). The HKMS discrepancy is therefore independently confirmed.
- **LLL witnesses (Table 2), checked 2026-10-15** (`src/independent/pari_check_lll.gp`, \p 130; log
  `results/2026-10-10_pari_independent/run4_lll_witnesses.log`): PARI gives 2.683073787077661771,
  −2.784942445945935525 (Q(√−3)) and 1.764536680198591170, −1.755074060289126858 (Q(√5)), agreeing with the Arb
  certificates to all 18 printed digits. PARI's zeros below T = 1000 (1472 and 1553) agree with our certified
  100-digit zeros to ≥ 102.3 decimal places each (`src/independent/compare_zeros.py`).

---

## 2026-10-08 — Stronger bounds by lattice reduction (rigorous)

Zeros below T = 1000 refined to 100 digits (`src/refine_zeros.py`, sign changes certified at radius 1e-102);
Odlyzko–te Riele-type simultaneous-approximation lattices (`src/lll_search.py`, FLINT LLL, n = 20–40 terms,
ν = 100–280 bits) swept over 96 parameter sets (`src/lll_sweep.py`). Best witnesses certified with ε = 1e-100, 600 bits:

| Field | M_K^+ > | M_K^− < |
|---|---|---|
| Q(√5) | **1.7645** | **−1.7550** |
| Q(√−3) | **2.6830** | **−2.7849** |

(Compare: direct scan at T = 1000 gave 1.130 / −1.111 and 1.683 / −1.770.) Best runs used n = 40, weight √c.

**Incident.** The first sweep parsed zeros at mpmath's default 15 digits, so every h*(t) with |t| > 10^15 was noise;
the reported "best" values (up to ±2.30) were selection-biased noise. Certification of four of them gave |h*| < 0.5,
exposing the bug. Fixed (precision set before parsing, plus an assertion); the folder is kept and marked INVALID.

---

## 2026-10-08 — Certification (rigorous)

Code: `src/certify/`. Certificates: `certificates/` (summary table in certificates/README.md).

**Theorem (computer-assisted).** The naive Mertens-type conjecture fails for K = Q(√−3) and K = Q(√5):
- Q(√−3): h*_{K,200}(1560.51100774) ∈ −1.27551112261409151 ± 6.4·10^-21 and h*_{K,200}(18558.8477306) ∈
  +1.20792583877255878 ± 2.9·10^-20; hence M_K^− < −1.2755 and M_K^+ > 1.2079.
- Q(√5): h*_{K,1000}(1016.24860331) ∈ +1.13020150032987766 ± 4.2·10^-22 and h*_{K,1000}(722.734868673) ∈
  −1.11075627303412053 ± 1.4·10^-20; hence M_K^+ > 1.1302 and M_K^− < −1.1107.
(Inequalities M^− ≤ h* ≤ M^+ via HKMS §3.2 / Jurkat; if RH fails or zeros are multiple the conjecture fails anyway.)

Proof ingredients: certified zero counts by the argument principle (Q(√−3), T=200: 79 + 114; Q(√5), T=1000: 649 + 904),
sign-change isolation of every zero, derivative enclosures by Taylor-series ball arithmetic. Real-axis positivity
L(σ,χ) > 0 (σ > 0) by elementary grouping/convexity lemmas (Arb cannot evaluate L on balls containing the Hurwitz pole at
σ = 1). Robust under 300 → 450 bits and ε 1e-28 → 1e-33; negative controls behave correctly.

**Correction to HKMS (rigorous).** Certified enclosures give h*_{Q(√2),150}(−24.64) = −0.7979…,
h*_{Q(√3),200}(17.32) = −0.3004…, h*_{Q(i),600}(72.85) = −0.9367…, contradicting the published values (> 1.049, > 1.027,
≤ −1.008) under the paper's eq. (3.6). Replacement witnesses at the same T are certified (Q(√2): ±1.30/−1.26;
Q(√3): ±1.16; Q(i): ±1.57). HKMS's h*_{Q(i),600}(−85.15) ≥ 1.029 is confirmed (1.2002…). All their theorems stand.

---

## 2026-10-08 — Feasibility study (numerical)

Logs: `results/2026-10-08_feasibility/`. Data: `data/zeros/quadratic_T1000.json`.

### Zero data
All zeros of ζ (Arb, rigorous) and L(s,χ_D) (sign-change scan, step 0.02) with 0 < γ < 1000 for Q(i), Q(√2), Q(√3),
Q(√−3), Q(√5) (and Q(√−2) to 700). L-zero counts agree with θ(T)/π to within ±0.7 at T = 250, 500, 750, 1000.
Min |γ_ζ − γ_L| ≥ 1.6·10^-4; all derivatives bounded away from 0. First zeros agree with known values
(L(χ_{−4}): 6.0209489047; L(χ_{−3}): 8.0397371557).

### Ground-truth validation
M_K(x) computed directly to 10^6 (b = μ * μχ). Kernel-smoothed x^{-1/2}M_K vs h*_{1000}(log x): agreement within
0.01–0.06 for t ∈ [12, 13.8] in all five fields (smaller t differ by the omitted x^{-1/2}M*_K(x) term).
Jump at t = 0 reproduces the paper's M*_K(1): h(0^−) → −0.45, h(0^+) → +0.47 at T = 1000 (limits −0.485, +0.515).

### Main result (numerical)
| Field | T | t | h*_{K,T}(t) (40 digits) |
|---|---|---|---|
| Q(√−3) | 200 | 1560.51100774 | −1.27551112 |
| Q(√−3) | 200 | 18558.8477306 | +1.20792584 |
| Q(√−3) | 400 | 1560.51013516 | −1.46222247 |
| Q(√−3) | 1000 | 800.889 / 477.692 | −1.7700 / +1.6828 (float scan) |
| Q(√5) | 1000 | 1016.24860331 | +1.13020150 |
| Q(√5) | 1000 | 722.734868673 | −1.11075627 |
| Q(√5) | 400 | 1016.244 / 722.739 | +1.0083 / −1.0057 (float scan) |
| Q(√5) | 200 | (scan |t| ≤ 2·10^4) | max 0.867, min −0.836 — insufficient |

Amplitude budget U(T) = Σ_{|γ|<T} f(|γ|/T)/|ρζ_K′(ρ)| at T = 200 / 400 / 1000:
Q(√−3) 1.85 / 2.69 / 4.18; Q(√5) 1.42 / 2.06 / 2.98 (Q(√5) is the hard case).

### Irreproducible witnesses in the 2024 paper
Our implementation of eq. (3.6), validated against ground truth, gives:

| Claim (paper) | Our value | Best within ±0.006 |
|---|---|---|
| Q(i): h*_{600}(72.85) ≤ −1.008 | −0.9367 | −0.9447 |
| Q(i): h*_{600}(−85.15) ≥ 1.029 | +1.2002 ✔ | +1.2210 |
| Q(√2): h*_{150}(−24.64) > 1.049 | −0.7979 | −0.7820 |
| Q(√3): h*_{200}(17.32) > 1.027 | −0.3004 | −0.2662 |

Ruled out: rounding of t, t ↦ −t / conjugation, γ>0-only sums, missing 1/ρ, no kernel, T as zero count, character
mix-ups (all six fields incl. Q(√−2) evaluated at each claimed (T,t)). The theorems remain true: at the same T our scan
finds Q(√2) max 1.3034 (t ≈ 2607.18), min −1.2559 (t ≈ 195.91); Q(√3) ±1.16 (t ≈ 631.34, 1746.88); Q(i) ±1.57.
Wording for any write-up: "we were unable to reproduce", not "incorrect". Consider contacting the authors.
