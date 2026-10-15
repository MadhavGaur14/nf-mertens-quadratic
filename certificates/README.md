# Certificates

Produced by `src/certify/run_certificate.py` (Arb ball arithmetic via python-flint 0.9.0, Python 3.14). Each `.json`
contains the field, character, T, working precision, isolation radius ε, the certified zero counts, the midpoints of all
isolating intervals, and the enclosure of h*_{K,T}(t) at each requested t. Each `.log` is the console transcript.

Proof chain per certificate (docs/04_certification_plan.md):
1. N(T) for ζ and L(s,χ) by the argument principle (certified integer).
2. A sign change of the Hardy Z-function on [γ−ε, γ+ε] for every zero; #intervals = N(T) ⇒ all zeros with 0<γ<T are
   on the critical line, simple, and enclosed; ζ- and L-intervals pairwise disjoint.
3. Ball enclosure of h*_{K,T}(t).

## Results (2026-10-08)

| Certificate | Witness t | Enclosure of h*_{K,T}(t) | Verdict |
|---|---|---|---|
| Qsqrt-3_T200 | 1560.51100774 | −1.27551112261409151 ± 6.4e-21 | **< −1** |
| Qsqrt-3_T200 | 18558.8477306 | +1.20792583877255878 ± 2.9e-20 | **> 1** |
| Qsqrt5_T1000 | 1016.24860331 | +1.13020150032987766 ± 4.2e-22 | **> 1** |
| Qsqrt5_T1000 | 722.734868673 | −1.11075627303412053 ± 1.4e-20 | **< −1** |
| Qsqrt2_T150 | 2607.183 | +1.30416564176212412 ± 9.2e-21 | > 1 |
| Qsqrt2_T150 | 195.906 | −1.25583009973636698 ± 1.2e-20 | < −1 |
| Qsqrt2_T150 | −24.64 (HKMS witness) | −0.79792323924516169 ± 1.2e-21 | not certified (HKMS claim > 1.049) |
| Qsqrt3_T200 | 631.336 | +1.16443450961437372 ± 4.1e-21 | > 1 |
| Qsqrt3_T200 | 1746.876 | −1.16296223358856553 ± 1.8e-20 | < −1 |
| Qsqrt3_T200 | 17.32 (HKMS witness) | −0.30038943431042942 ± 2.9e-21 | not certified (HKMS claim > 1.027) |
| Qi_T600 | 158.854 | +1.57214733418839255 ± 1.7e-20 | > 1 |
| Qi_T600 | 160.454 | −1.55457646591422597 ± 2.2e-20 | < −1 |
| Qi_T600 | −85.15 (HKMS witness) | +1.20021579231940321 ± 4.6e-20 | > 1 (reproduces HKMS) |
| Qi_T600 | 72.85 (HKMS witness) | −0.93672311792431125 ± 2.7e-21 | not certified (HKMS claim ≤ −1.008) |

### Stronger bounds via lattice reduction (LLL), T = 1000, 100-digit zeros, ε = 1e-100, 600 bits

| Certificate | Witness t (integer) | Enclosure of h*_{K,1000}(t) | Bound |
|---|---|---|---|
| Qsqrt5_T1000_LLLv2 | 6701024002746072686211658750800192982264466383201014802804 | +1.76453668019859117 ± 1.3e-20 | **M^+ > 1.7645** |
| Qsqrt5_T1000_LLLv2 | −3826144185540628524340499373955955297176444083490968202937 | −1.75507406028912686 ± 1.5e-20 | **M^− < −1.7550** |
| Qsqrt-3_T1000_LLLv2 | 165562943829466520552313925650444335142025426569159622719799633443383327010312592 | +2.68307378707766177 ± 5.8e-20 | **M^+ > 2.6830** |
| Qsqrt-3_T1000_LLLv2 | 3431377126899343121084375092595828581527345096545547547377 | −2.78494244594593553 ± 4.4e-21 | **M^− < −2.7849** |
| Qsqrt5_T1000_LLL | 1253089510300830588081520614430270334673580 | +1.68580041679481615 ± 4.3e-20 | (superseded) |
| Qsqrt-3_T1000_LLL | −7627565930529569320173884736841401921126271 | −2.23328909155063425 ± 1.7e-20 | (superseded) |

Witnesses found by `src/lll_search.py` (sweep log: `results/2026-10-08_lll_v2/`). Note: an earlier sweep
(`results/2026-10-08_lll_INVALID_precision_bug/`) had a precision bug; its candidates failed certification and are void.

Zero counts: Q(√−3), T=200: N_ζ = 79, N_L = 114. Q(√5), T=1000: N_ζ = 649, N_L = 904. Q(√2), T=150: 52, 102.
Q(√3), T=200: 79, 158. Q(i), T=600: 341, 473.

`robustness/`: Q(√−3) and Q(√5) re-run at 450 bits, ε = 1e-33 (identical enclosures), with negative controls
(t = 1560.6 → +0.4396; t = 0 → +0.01367; both correctly "NOT CERTIFIED").
