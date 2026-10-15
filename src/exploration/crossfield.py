"""Evaluate each claimed (T, t) for every field we have (mix-up hypothesis), with local extremum within +-0.006."""
import sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import load_fields

A = load_fields()


def kern(y):
    return (1 - y) * np.cos(np.pi * y) + np.sin(np.pi * y) / np.pi


def h(name, T, t):
    z = A[name]['zeros']
    g = np.array([float(r[1]) for r in z])
    a = np.array([float(r[2]) + 1j * float(r[3]) for r in z])
    m = g < T
    return np.real(np.exp(1j * np.outer(np.atleast_1d(t), g[m])) @ (2 * kern(g[m] / T) * a[m]))


for lab, T, t, s in [('Q(i) claim <=-1.008', 600, 72.85, -1), ('Q(i) claim >=1.029', 600, -85.15, 1),
                     ('Q(sqrt2) claim >1.049', 150, -24.64, 1), ('Q(sqrt3) claim >1.027', 200, 17.32, 1)]:
    print(f"\n{lab}: T={T}, t={t}")
    for name in A:
        tt = np.linspace(t - 0.006, t + 0.006, 601)
        v = h(name, T, tt)
        best = v.max() if s > 0 else v.min()
        print(f"   {name:10s}  h*(t) = {h(name, T, t)[0]:+.4f}   best nearby = {best:+.4f}")
