"""Odlyzko-te Riele-type search for large |h*_{K,T}(t)| by lattice reduction.
h*_{K,T}(t) = sum_j c_j cos(gamma_j t + phi_j),  c_j = 2 f(gamma_j/T)|a_j|,  phi_j = arg a_j   (0 < gamma_j < T).
Pick the n terms with the largest c_j and look for t with gamma_j t + phi_j (+ pi for a minimum) = 0 mod 2pi for all of
them simultaneously: inhomogeneous simultaneous Diophantine approximation via LLL on the (n+2)-dimensional lattice
   row j (1..n):  W_j e_j                                     (W_j = round(2^nu w_j), w_j = weight of term j)
   row n+1:       (round(2^nu w_j gamma_j/(2 pi)) )_j , S , 0  (coefficient = t, an integer multiple of 1)
   row n+2:       (round(2^nu w_j beta_j))_j       , 0 , M  (embedding row; beta_j = phi_j/(2pi) [+1/2])
A short vector with last coordinate +-M gives t with all weighted errors small.
usage: python src/lll_search.py FIELD T n nu [--sign +1|-1] [--tmax-digits D] [--zeros FILE]"""
import argparse, json, math, sys, time
from pathlib import Path
import mpmath as mp
import numpy as np
from flint import fmpz_mat
sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import DATA, kern_np

ap = argparse.ArgumentParser()
ap.add_argument('field')
ap.add_argument('T', type=float)
ap.add_argument('n', type=int)
ap.add_argument('nu', type=int)
ap.add_argument('--sign', type=int, default=1)
ap.add_argument('--zeros', default=str(DATA / 'quadratic_T1000.json'))
ap.add_argument('--weight', default='c', choices=['c', 'sqrtc', 'one'])
ap.add_argument('--S', type=float, default=1.0, help='weight of the t coordinate')
ap.add_argument('--rows', type=int, default=999, help='number of reduced rows to try')
ap.add_argument('--mfac', type=float, default=0.1, help='embedding constant M = mfac * 2^nu (Kannan embedding)')
A = ap.parse_args()

rec = json.load(open(A.zeros))['fields'][A.field]['zeros']
digits_zero = min(len(r[1].split('.')[1]) for r in rec if float(r[1]) < A.T)
# BUGFIX 2026-10-08: the working precision must be set BEFORE the zeros are parsed; parsing at mpmath's default
# 15 digits silently truncated gamma and made every h*(t) with |t| > 1e15 meaningless.
mp.mp.dps = max(digits_zero + 10, A.nu // 3 + 40)
G = [(mp.mpf(r[1]), mp.mpc(r[2], r[3])) for r in rec if float(r[1]) < A.T]
assert abs(G[0][0] - mp.mpf(rec[0][1])) == 0 or len(rec[0][1]) < mp.mp.dps, "zero parsed at too low precision"
gam = np.array([float(g) for g, _ in G])
amp = np.array([2 * kern_np(float(g) / A.T) * abs(complex(a)) for g, a in G])
order = np.argsort(-amp)[:A.n]
U = amp.sum()
print(f"{A.field} T={A.T}: {len(G)} zeros, amplitude budget U = {U:.4f}; top-{A.n} terms carry {amp[order].sum():.4f};"
      f" zeros known to {digits_zero} digits", flush=True)

C = mp.mpf(2) ** A.nu
w = {'c': amp[order] / amp[order].max(), 'sqrtc': np.sqrt(amp[order] / amp[order].max()), 'one': np.ones(A.n)}[A.weight]
n = A.n
rows = []
for k, j in enumerate(order):
    r = [0] * (n + 2)
    r[k] = int(mp.nint(C * w[k]))
    rows.append(r)
x = [mp.mpf(G[j][0]) / (2 * mp.pi) for j in order]
beta = [mp.arg(G[j][1]) / (2 * mp.pi) + (mp.mpf(1) / 2 if A.sign < 0 else 0) for j in order]
S = int(round(A.S))
M = int(mp.nint(C * A.mfac))
rows.append([int(mp.nint(C * w[k] * x[k])) for k in range(n)] + [S, 0])
rows.append([int(mp.nint(C * w[k] * beta[k])) for k in range(n)] + [0, M])
t0 = time.time()
R = fmpz_mat(rows).lll()
print(f"LLL dim {n + 2}, nu={A.nu}: {time.time() - t0:.1f}s", flush=True)


def h(t, dps):
    with mp.workdps(dps):
        t = mp.mpf(t)
        s = mp.mpf(0)
        for g, a in G:
            y = g / A.T
            f = (1 - y) * mp.cos(mp.pi * y) + mp.sin(mp.pi * y) / mp.pi
            s += 2 * f * mp.re(a * mp.expj(g * t))
        return s


best = None
for i in range(min(A.rows, n + 2)):
    v = [int(R[i, j]) for j in range(n + 2)]
    if abs(v[-1]) != M:
        continue
    if v[-1] < 0:
        v = [-c for c in v]
    t = v[n] // S
    if t == 0:
        continue
    lt = math.log10(abs(t))
    if lt > digits_zero - 12:
        print(f"  row {i}: t ~ 1e{lt:.1f} too large for {digits_zero}-digit zeros; skipped")
        continue
    val = h(t, int(lt) + 30)
    print(f"  row {i}: t = {t}  (~1e{lt:.1f})  h* = {mp.nstr(val, 12)}", flush=True)
    if best is None or A.sign * val > A.sign * best[1]:
        best = (t, val)
if best:
    print(f"BEST: {A.field} T={A.T} sign={A.sign:+d}: h*({best[0]}) = {mp.nstr(best[1], 15)}")
