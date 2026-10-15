"""Feasibility analysis for the naive Mertens-type conjecture over quadratic fields (Hu-Kaneko-Martin-Schildkraut).
  h*_{K,T}(t) = sum_{|gamma|<T} f(|gamma|/T) e^{i gamma t} / (rho zeta_K'(rho)),  f(y) = (1-y)cos(pi y) + sin(pi y)/pi
A value h* > 1 (or < -1) at any t implies M_K^+ > 1 (M_K^- < -1)  [their eq. (3.6) and the remark after it].
usage: python nf_scan.py nf_data_1000.json"""
import json, sys, time
import numpy as np

D = json.load(open(sys.argv[1]))
F = {}
for name, rec in D['fields'].items():
    z = rec['zeros']
    F[name] = (np.array([float(r[1]) for r in z]), np.array([float(r[2]) + 1j * float(r[3]) for r in z]),
               np.array([r[0] for r in z]))


def kern(y):
    return (1 - y) * np.cos(np.pi * y) + np.sin(np.pi * y) / np.pi


def coeffs(name, T):
    g, a, _ = F[name]
    m = g < T
    return g[m], 2 * kern(g[m] / T) * a[m]          # factor 2: conjugate zeros (chi real)


def h(name, T, t):
    g, c = coeffs(name, T)
    t = np.atleast_1d(np.asarray(t, dtype=float))
    return np.real(np.exp(1j * np.outer(t, g)) @ c)


print("=== Positive controls (values claimed in Hu-Kaneko-Martin-Schildkraut, Section 3)")
for name, T, t, claim in [('Q(i)', 600, 72.85, '<= -1.008'), ('Q(i)', 600, -85.15, '>= 1.029'),
                          ('Q(sqrt2)', 150, -24.64, '> 1.049'), ('Q(sqrt3)', 200, 17.32, '> 1.027')]:
    print(f"  {name:9s} T={T:4d} t={t:7.2f}: h* = {h(name, T, t)[0]: .4f}   (paper: {claim})")

print("\n=== Jump at t=0 (should straddle M*(1) ~ -0.485 and 1+M*(1) ~ 0.515 as T grows)")
for name in ['Q(sqrt-3)', 'Q(sqrt5)']:
    print(f"  {name:9s} " + "  ".join(f"T={T}: h(-0.02)={h(name, T, -0.02)[0]:+.3f} h(+0.02)={h(name, T, 0.02)[0]:+.3f}"
                                      for T in (200, 1000)))

print("\n=== Amplitude budget U(T) = sum_{|gamma|<T} f(|gamma|/T)/|rho zeta_K'(rho)|  (sup of |h*| over all t)")
Ts = [100, 200, 300, 400, 600, 800, 1000]
print("  field        " + "".join(f"{T:>8d}" for T in Ts) + "     #zeros<1000")
for name in F:
    us = [np.abs(coeffs(name, T)[1]).sum() for T in Ts]
    print(f"  {name:11s}  " + "".join(f"{u:8.3f}" for u in us) + f"   {len(F[name][0]):6d}")


def scan(name, T, R, pts_per_period=8, chunk=200_000):
    g, c = coeffs(name, T)
    step = 2 * np.pi / T / pts_per_period
    best_hi, best_lo = (-9, 0), (9, 0)
    t0 = time.time()
    for start in np.arange(-R, R, chunk * step):
        t = start + step * np.arange(chunk)
        v = np.real(np.exp(1j * np.outer(t, g)) @ c)
        i, j = v.argmax(), v.argmin()
        if v[i] > best_hi[0]:
            best_hi = (v[i], t[i])
        if v[j] < best_lo[0]:
            best_lo = (v[j], t[j])
    return best_hi, best_lo, len(g), time.time() - t0


print("\n=== Direct scan of h*_{K,T}(t) for t in [-R, R]")
for name in ['Q(sqrt-3)', 'Q(sqrt5)', 'Q(i)']:
    for T, R in [(200, 20000), (400, 5000), (1000, 1000)]:
        hi, lo, n, dt = scan(name, T, R)
        flag = '  <-- exceeds 1 !' if hi[0] > 1 or lo[0] < -1 else ''
        print(f"  {name:9s} T={T:5d} ({n:4d} zeros) R={R:6d}: max {hi[0]:.4f} at t={hi[1]:.3f};"
              f"  min {lo[0]:.4f} at t={lo[1]:.3f}   [{dt:.0f}s]{flag}", flush=True)
