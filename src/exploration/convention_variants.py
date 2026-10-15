"""Which convention (if any) reproduces the four values claimed in Hu-Kaneko-Martin-Schildkraut?"""
import json, sys
import numpy as np
from flint import dirichlet_char, acb, arb, ctx

D = json.load(open(sys.argv[1]))
CLAIMS = [('Q(i)', 600, 72.85, '<= -1.008'), ('Q(i)', 600, -85.15, '>= 1.029'),
          ('Q(sqrt2)', 150, -24.64, '> 1.049'), ('Q(sqrt3)', 200, 17.32, '> 1.027')]


def kern(y):
    return (1 - y) * np.cos(np.pi * y) + np.sin(np.pi * y) / np.pi


def load(name):
    z = D['fields'][name]['zeros']
    g = np.array([float(r[1]) for r in z])
    a = np.array([float(r[2]) + 1j * float(r[3]) for r in z])
    src = np.array([r[0] for r in z])
    return g, a, src


VARIANTS = {
    'ours: 2Re sum f a e^{igt}': lambda g, a, s, t, T: np.real(2 * np.sum(kern(g / T) * a * np.exp(1j * g * t))),
    'conjugate e^{-igt}': lambda g, a, s, t, T: np.real(2 * np.sum(kern(g / T) * a * np.exp(-1j * g * t))),
    'gamma>0 only, Re': lambda g, a, s, t, T: np.real(np.sum(kern(g / T) * a * np.exp(1j * g * t))),
    'gamma>0 only, |.|': lambda g, a, s, t, T: np.abs(np.sum(kern(g / T) * a * np.exp(1j * g * t))),
    '1/zeta_K\'(rho) (no 1/rho)': lambda g, a, s, t, T: np.real(2 * np.sum(kern(g / T) * a * (0.5 + 1j * g) * np.exp(1j * g * t))),
    'zeta zeros only': lambda g, a, s, t, T: np.real(2 * np.sum((s == 'zeta') * kern(g / T) * a * np.exp(1j * g * t))),
    'L zeros only': lambda g, a, s, t, T: np.real(2 * np.sum((s == 'L') * kern(g / T) * a * np.exp(1j * g * t))),
    'no kernel (f=1)': lambda g, a, s, t, T: np.real(2 * np.sum((g < T) * a * np.exp(1j * g * t))),
    'kernel f(g/T) with T=#zeros?': None,
}

print(f"{'variant':32s}" + "".join(f"{c[0] + ' t=' + str(c[2]):>22s}" for c in CLAIMS))
print(f"{'paper claims':32s}" + "".join(f"{c[3]:>22s}" for c in CLAIMS))
for vn, fn in VARIANTS.items():
    if fn is None:
        continue
    row = []
    for name, T, t, _ in CLAIMS:
        g, a, s = load(name)
        m = g < T
        row.append(fn(g[m], a[m], s[m], t, T))
    print(f"{vn:32s}" + "".join(f"{v:22.4f}" for v in row))
# T interpreted as the number of zeros of zeta_K (in the upper half plane) rather than a height
row = []
for name, N, t, _ in CLAIMS:
    g, a, s = load(name)
    T = g[N - 1] + 1e-9 if N <= len(g) else g[-1]
    m = g < T
    row.append(np.real(2 * np.sum(kern(g[m] / T) * a[m] * np.exp(1j * g[m] * t))))
print(f"{'T := gamma_N (N = stated T)':32s}" + "".join(f"{v:22.4f}" for v in row))

# Wrong-character hypothesis: real fields computed with the other real character of the same modulus
print("\nOther real characters of the same modulus (in case of a character mix-up):")
for q in (8, 12):
    for l in range(1, q):
        try:
            c = dirichlet_char(q, l)
        except Exception:
            continue
        if c.is_real() and not c.is_principal():
            print(f"  ({q},{l}) primitive={c.is_primitive()} parity={c.parity()} conductor={c.conductor()}")
