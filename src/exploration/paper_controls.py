"""Positive controls: local extrema of h*_{K,T} near the (2-decimal) t values reported by Hu-Kaneko-Martin-Schildkraut."""
import json, sys
import numpy as np

D = json.load(open(sys.argv[1]))


def kern(y):
    return (1 - y) * np.cos(np.pi * y) + np.sin(np.pi * y) / np.pi


def h(name, T, t):
    z = D['fields'][name]['zeros']
    g = np.array([float(r[1]) for r in z])
    a = np.array([float(r[2]) + 1j * float(r[3]) for r in z])
    m = g < T
    c = 2 * kern(g[m] / T) * a[m]
    return np.real(np.exp(1j * np.outer(np.atleast_1d(t), g[m])) @ c)


for name, T, t, s, claim in [('Q(i)', 600, 72.85, -1, '<= -1.008'), ('Q(i)', 600, -85.15, 1, '>= 1.029'),
                             ('Q(sqrt2)', 150, -24.64, 1, '> 1.049'), ('Q(sqrt3)', 200, 17.32, 1, '> 1.027')]:
    print(f"{name:9s} T={T} t={t} (paper: h* {claim}); value at t: {h(name, T, t)[0]:+.4f}")
    for w in [0.005, 0.01, 0.05]:
        tt = np.linspace(t - w, t + w, 4001)
        v = s * h(name, T, tt)
        i = v.argmax()
        print(f"     best within +-{w}: h* = {s * v[i]:+.4f} at t = {tt[i]:.4f}")
