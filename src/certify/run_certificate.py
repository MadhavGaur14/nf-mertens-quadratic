"""Produce a certificate that h*_{K,T}(t) lies outside [-1, 1] for given witnesses t.
usage: python src/certify/run_certificate.py FIELD T t1 [t2 ...] [--prec 300] [--eps 1e-28]
writes certificates/<FIELD>_T<T>.json and .log"""
import argparse, json, sys, time
from datetime import datetime
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from flint import arb, acb, ctx, dirichlet_char
from common import FIELDS, CERTS, load_fields
from arb_tools import ZETA, versions
from count_zeros import count_zeros
from isolate_zeros import isolate
from certify_hstar import coefficients, hstar

ap = argparse.ArgumentParser()
ap.add_argument('field')
ap.add_argument('T')
ap.add_argument('ts', nargs='+')
ap.add_argument('--prec', type=int, default=300)
ap.add_argument('--eps', default='1e-28')
ap.add_argument('--zeros', default=None, help='zero file (default: data/zeros/quadratic_T1000.json); use the hp_* '
                                             'files when |t| is large (need eps*|t| << 1)')
ap.add_argument('--tag', default='', help='suffix for the certificate file name')
A = ap.parse_args()

CERTS.mkdir(exist_ok=True)
tag = f"{A.field}_T{A.T}{A.tag}".replace('(', '').replace(')', '')
logf = open(CERTS / f"{tag}.log", 'w', encoding='utf-8')


def log(*a):
    s = ' '.join(str(x) for x in a)
    print(s, flush=True)
    logf.write(s + '\n')


q, l = FIELDS[A.field]
chi = dirichlet_char(q, l)
T, Tf = arb(A.T), float(A.T)
cert = {'field': A.field, 'character': [q, l], 'T': A.T, 'eps': A.eps, 'prec': A.prec, 'versions': versions(),
        'date': datetime.now().isoformat(timespec='seconds')}
log(f"== Certificate for the naive Mertens-type conjecture over {A.field}; chi = dirichlet_char{(q, l)}; T = {A.T}")
t0 = time.time()

# Step 1: rigorous zero counts
ctx.prec = 128
NL, infoL = count_zeros(chi, A.T)
NZ, infoZ = count_zeros(ZETA, A.T)
log(f"[1] N_L(T) = {NL}  ({infoL});  N_zeta(T) = {NZ}  ({infoZ})  [{time.time() - t0:.0f}s]")

# Step 2: isolation
ctx.prec = A.prec
rec = (load_fields(Path(A.zeros)) if A.zeros else load_fields())[A.field]['zeros']
for ts in A.ts:                                  # guard: phase uncertainty eps*|t| must be tiny
    assert (arb(A.eps) * abs(arb(ts))).upper() < arb(10) ** -15, f"eps too large for t = {ts}; use hp zeros"
gz = [r[1] for r in rec if r[0] == 'zeta' and float(r[1]) < Tf]
gl = [r[1] for r in rec if r[0] == 'L' and float(r[1]) < Tf]
iz, il = isolate(ZETA, gz, A.eps), isolate(chi, gl, A.eps)
eps = arb(A.eps)
for lst in (iz, il):
    assert (lst[0] - eps).lower() > 0 and (lst[-1] + eps).upper() < T, "interval outside (0, T)"
assert len(iz) == NZ, f"zeta: isolated {len(iz)} != N(T) = {NZ}"
assert len(il) == NL, f"L: isolated {len(il)} != N(T) = {NL}"
allz = sorted([(float(g), 'zeta', g) for g in iz] + [(float(g), 'L', g) for g in il])
for (fa, sa, a), (fb, sb, b) in zip(allz, allz[1:]):
    assert (b - a - 2 * eps).lower() > 0, f"zeta/L intervals not disjoint near {fa}"
# cross-check zeta zeros against Arb's own rigorous zeta_zeros
arbz = [z.imag for z in acb.zeta_zeros(1, NZ)]
assert all(abs(a - b).upper() < eps for a, b in zip(arbz, iz)), "mismatch vs acb.zeta_zeros"
log(f"[2] isolated {len(iz)} zeta zeros and {len(il)} L zeros in disjoint intervals of radius {A.eps} inside (0,T);"
    f" counts equal N(T) => all zeros below T are on the line and simple; zeta zeros match acb.zeta_zeros"
    f"  [{time.time() - t0:.0f}s]")
cert['N_zeta'], cert['N_L'] = NZ, NL
cert['zeta_zero_midpoints'] = [g.str(40, radius=False) for g in iz]
cert['L_zero_midpoints'] = [g.str(40, radius=False) for g in il]

# Step 3: enclosure of h*
co = coefficients(chi, iz, il, A.eps)
res = []
for ts in A.ts:
    h = hstar(co, A.T, ts)
    verdict = 'CERTIFIED > 1' if h.lower() > 1 else 'CERTIFIED < -1' if h.upper() < -1 else 'NOT CERTIFIED'
    log(f"[3] h*_(K,T)({ts}) in {h.str(20)}   -> {verdict}")
    res.append({'t': ts, 'enclosure': h.str(30), 'lower': h.lower().str(20, radius=False),
                'upper': h.upper().str(20, radius=False), 'verdict': verdict})
cert['witnesses'] = res
cert['seconds'] = round(time.time() - t0, 1)
json.dump(cert, open(CERTS / f"{tag}.json", 'w'), indent=1)
log(f"done in {time.time() - t0:.0f}s -> certificates/{tag}.json")
