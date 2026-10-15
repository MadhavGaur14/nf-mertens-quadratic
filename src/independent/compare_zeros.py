"""Digit-by-digit comparison of PARI/GP zeros (pari_check_lll.gp output) with the Arb-certified 100-digit zeros.
Reports the count match and the smallest number of agreeing decimals over all zeros."""
import json, sys
from pathlib import Path
import mpmath as mp

ROOT = Path(__file__).resolve().parent.parent.parent
RES = ROOT / 'results' / '2026-10-10_pari_independent'
mp.mp.dps = 130
for tag, field in [('Qsqrt-3', 'Q(sqrt-3)'), ('Qsqrt5', 'Q(sqrt5)')]:
    pf = RES / f'pari_zeros_{tag}_T1000.txt'
    if not pf.exists():
        print(f"{field}: {pf.name} not found (PARI run not finished?)")
        continue
    pari = sorted(mp.mpf(l.strip()) for l in pf.read_text().splitlines() if l.strip())
    arb = sorted(mp.mpf(r[1]) for r in json.load(open(ROOT / 'data' / 'zeros' / f'hp_{tag}_T1000_D100.json'))['fields'][field]['zeros'])
    print(f"{field}: PARI {len(pari)} zeros, Arb-certified {len(arb)} zeros -> {'counts agree' if len(pari) == len(arb) else 'COUNT MISMATCH'}")
    if len(pari) == len(arb):
        worst = min(-mp.log10(abs(a - b)) if a != b else 130 for a, b in zip(pari, arb))
        print(f"   every zero agrees to at least {float(worst):.1f} decimal places")
