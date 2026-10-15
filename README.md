# The naive Mertens-type conjecture for Q(√−3) and Q(√5)

Research project, started 2026-10-08.

**Question.** For a number field K let M_K(x) = Σ_{N(𝔞) ≤ x} μ_K(𝔞). Hu, Kaneko, Martin & Schildkraut
(*Math. Proc. Camb. Phil. Soc.* 177 (2024) 481–523, arXiv:2109.06665) disproved the "naive Mertens-type conjecture"
−1 ≤ liminf M_K(x)/√x ≤ limsup M_K(x)/√x ≤ 1 for every quadratic field **except Q(√−3) and Q(√5)**, which they left
open, suggesting an Odlyzko–te Riele-type computation. This project settles those two cases.

**Status (2026-10-08).** **Certified:** the conjecture fails for both fields (computer-assisted proof, Arb ball
arithmetic): Q(√−3) M^+ > 2.6830, M^− < −2.7849; Q(√5) M^+ > 1.7645, M^− < −1.7550 (LLL witnesses, T = 1000;
first disproofs at T = 200 / 1000 gave ±1.2 / ±1.1). Three of the four explicit witnesses
published by HKMS do not hold; certified replacements are provided. See [certificates/README.md](certificates/README.md)
and [docs/03_findings_log.md](docs/03_findings_log.md). Next: larger T / LLL for stronger bounds, write-up.

**Archived record.** Paper, code, zero data and certificates: Zenodo, [doi:10.5281/zenodo.23403118](https://doi.org/10.5281/zenodo.23403118).
Repository: <https://github.com/MadhavGuar14/nf-mertens-quadratic>

## Layout

| Path | Contents |
|---|---|
| `docs/01_problem_statement.md` | Precise statement, method, and the logic of the disproof |
| `docs/02_literature_review.md` | Prior art, current records, competition check (with URLs) |
| `docs/03_findings_log.md` | Dated log of results, including the irreproducible witnesses in the 2024 paper |
| `docs/04_certification_plan.md` | What must be proved rigorously, and how |
| `docs/references.md` | Bibliography (verified links only) |
| `references/` | Local copies of key papers (not included in the public repository) |
| `src/common.py` | Shared paths, field ↔ character table, kernel |
| `src/compute_zeros.py` | Zeros of ζ·L(s,χ_D) and coefficients 1/(ρ ζ_K′(ρ)) (fast, non-rigorous) |
| `src/scan.py` | Scan h*_{K,T}(t) over t for large values |
| `src/certify/` | Rigorous (Arb ball-arithmetic) certification code |
| `src/exploration/` | One-off checks used during the feasibility study (ground truth, paper controls, conventions) |
| `data/zeros/` | Computed zero data (JSON) |
| `results/<date>_<topic>/` | Raw logs of each computational run |
| `certificates/` | Machine-checkable certificates produced by `src/certify/` |
| `paper/` | Manuscript (LaTeX source to be added once the arXiv version is posted) and `lint_tex.py` |
| `src/independent/` | Independent PARI/GP verification (`pari_check*.gp`) |

## Environment

Python 3.14, python-flint 0.9.0 (Arb ball arithmetic), mpmath 1.3.0, numpy 2.4.
`pip install python-flint mpmath numpy`

## Reproducing

```
python src/compute_zeros.py 1000                    # ~8 min; writes data/zeros/quadratic_T1000.json
python src/exploration/feasibility.py data/zeros/quadratic_T1000.json
python src/exploration/groundtruth.py data/zeros/quadratic_T1000.json 1e6
python src/exploration/hiprec_check.py "Q(sqrt-3):200:1560.511:-1" "Q(sqrt5):1000:1016.248:1"
python src/certify/<step>.py ...                    # see docs/04_certification_plan.md
```

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23403118.svg)](https://doi.org/10.5281/zenodo.23403118)