"""High-precision (mpmath, 40 digits) re-evaluation of h*_{K,T}(t) at the scan's extremal points, from the 35-digit
zeros, plus refinement of the local extremum."""
import sys
from pathlib import Path
import mpmath as mp
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import load_fields

mp.mp.dps = 40
A = load_fields()
# usage: python src/exploration/hiprec_check.py "FIELD:T:t:sign" ...   (sign = +1 for max, -1 for min)


def h(name, T, t):
    T, t = mp.mpf(T), mp.mpf(t)
    s = mp.mpf(0)
    for kind, g, ar, ai, _ in A[name]['zeros']:
        g = mp.mpf(g)
        if g >= T:
            break
        y = g / T
        f = (1 - y) * mp.cos(mp.pi * y) + mp.sin(mp.pi * y) / mp.pi
        s += 2 * f * mp.re(mp.mpc(ar, ai) * mp.expj(g * t))
    return s


for name, T, t, sgn in [(a, int(b), c, int(d)) for a, b, c, d in (x.split(':') for x in sys.argv[1:])]:
    t = mp.mpf(t)
    v0 = h(name, T, t)
    # golden-section refinement of the local extremum within +-0.003
    lo, hi = t - mp.mpf('0.003'), t + mp.mpf('0.003')
    phi = (mp.sqrt(5) - 1) / 2
    for _ in range(60):
        a, b = hi - phi * (hi - lo), lo + phi * (hi - lo)
        if sgn * h(name, T, a) > sgn * h(name, T, b):
            hi = b
        else:
            lo = a
    tb = (lo + hi) / 2
    print(f"{name:10s} T={T}: h*({mp.nstr(t, 10)}) = {mp.nstr(v0, 12)};  local extremum h*({mp.nstr(tb, 12)}) = {mp.nstr(h(name, T, tb), 12)}")
