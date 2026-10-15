"""Scan h*_{K,T}(t) on [-R, R] for chosen fields.   usage: python nf_scan_fields.py DATA.json FIELD:T:R [...]"""
import json, sys, time
import numpy as np

A = json.load(open(sys.argv[1]))['fields']


def kern(y):
    return (1 - y) * np.cos(np.pi * y) + np.sin(np.pi * y) / np.pi


for spec in sys.argv[2:]:
    name, T, R = spec.split(':')
    T, R = float(T), float(R)
    z = A[name]['zeros']
    g = np.array([float(r[1]) for r in z])
    a = np.array([float(r[2]) + 1j * float(r[3]) for r in z])
    m = g < T
    g, c = g[m], 2 * kern(g[m] / T) * a[m]
    step, chunk = 2 * np.pi / T / 8, 200_000
    hi, lo, t0 = (-9, 0), (9, 0), time.time()
    for start in np.arange(-R, R, chunk * step):
        t = start + step * np.arange(chunk)
        v = np.real(np.exp(1j * np.outer(t, g)) @ c)
        i, j = v.argmax(), v.argmin()
        hi = max(hi, (v[i], t[i]))
        lo = min(lo, (v[j], t[j]))
    print(f"{name:9s} T={T:5.0f} ({len(g)} zeros) R={R:6.0f}: max {hi[0]:.4f} at t={hi[1]:.3f};"
          f"  min {lo[0]:.4f} at t={lo[1]:.3f}   [{time.time() - t0:.0f}s]", flush=True)
