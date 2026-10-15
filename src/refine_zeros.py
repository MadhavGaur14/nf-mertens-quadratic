"""Refine all zeros below T of zeta and L(s, chi_D) to D significant decimals (needed when the witness t is huge:
the phase gamma*t must be known to ~1e-10, so gamma needs ~log10(t)+12 correct decimals).
Coefficients a = 1/(rho zeta_K'(rho)) are kept from the source file (they enter h* linearly, 15 digits suffice).
usage: python src/refine_zeros.py FIELD T D [--src FILE]   -> data/zeros/hp_<FIELD>_T<T>_D<D>.json"""
import argparse, json, sys, time
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent / 'certify'))
from flint import arb, ctx, dirichlet_char
from common import FIELDS, DATA
from arb_tools import ZETA
from isolate_zeros import sign_change

ap = argparse.ArgumentParser()
ap.add_argument('field')
ap.add_argument('T', type=float)
ap.add_argument('D', type=int)
ap.add_argument('--src', default=str(DATA / 'quadratic_T2000.json'))
A = ap.parse_args()

q, l = FIELDS[A.field]
chi = dirichlet_char(q, l)
ctx.prec = int(3.33 * (A.D + 45)) + 64
src = json.load(open(A.src))['fields'][A.field]
recs = [r for r in src['zeros'] if float(r[1]) < A.T]
eps = arb(10) ** (-(A.D + 2))


def newton_secant(fn, g0):
    """Secant iteration from two nearby points; converges superlinearly from a 30-digit start."""
    h = arb(10) ** -28
    x0, x1 = g0 - h, g0 + h
    f0, f1 = fn(x0).mid(), fn(x1).mid()
    for _ in range(40):
        if f1 == f0:
            break
        x2 = (x1 - f1 * (x1 - x0) / (f1 - f0)).mid()
        x0, f0, x1, f1 = x1, f1, x2, fn(x2).mid()
        if abs(x1 - x0) < eps / 100:
            break
    return x1


t0 = time.time()
out = []
for i, r in enumerate(recs):
    c = ZETA if r[0] == 'zeta' else chi
    g = newton_secant(lambda t: c.hardy_z(t).real, arb(r[1]))
    if not sign_change(c, g, eps):
        raise RuntimeError(f"refinement failed near {r[1]}")
    out.append([r[0], g.str(A.D + 5, radius=False), r[2], r[3], r[4]])
    if (i + 1) % 250 == 0:
        print(f"  {i + 1}/{len(recs)} ({time.time() - t0:.0f}s)", flush=True)
path = DATA / f"hp_{A.field.replace('(', '').replace(')', '')}_T{int(A.T)}_D{A.D}.json"
json.dump({'TMAX': A.T, 'digits': A.D, 'fields': {A.field: {**{k: v for k, v in src.items() if k != 'zeros'},
                                                             'zeros': out}}}, open(path, 'w'))
print(f"{A.field}: {len(out)} zeros refined to {A.D} digits (sign change certified at radius 1e-{A.D + 2})"
      f" in {time.time() - t0:.0f}s -> {path.name}")
