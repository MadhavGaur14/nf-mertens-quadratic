"""Ground truth: M_K(x) = sum_{N(a)<=x} mu_K(a) computed directly, compared with the zero sum h*_{K,T}(log x).
For quadratic K, 1/zeta_K = (1/zeta)(1/L(chi)), so the ideal-count coefficients are b = mu * (mu chi) and
M_K(x) = sum_{d<=x} mu(d) Mchi(x/d), Mchi(y) = sum_{m<=y} mu(m)chi(m).
h* is a smoothing of x^{-1/2}(M_K(x)+M*_K(x)) in u = log x with a kernel of width ~ 1/T."""
import json, sys
import numpy as np

D = json.load(open(sys.argv[1]))
N = int(float(sys.argv[2])) if len(sys.argv) > 2 else 10 ** 7


def mobius(N):
    mu = np.ones(N + 1, dtype=np.int8)
    isp = np.ones(N + 1, dtype=bool)
    isp[:2] = False
    for p in range(2, int(N ** 0.5) + 1):
        if isp[p]:
            isp[p * p::p] = False
    for p in np.nonzero(isp)[0]:
        p = int(p)
        mu[p::p] *= -1
        if p * p <= N:
            mu[p * p::p * p] = 0
    mu[0] = 0
    return mu


def kron(q, l, N):
    from flint import dirichlet_char
    c = dirichlet_char(q, l)
    per = np.array([int(round(float(c(n).real.mid()))) for n in range(q)], dtype=np.int8)
    return per[np.arange(N + 1) % q]


def kern(y):
    return (1 - y) * np.cos(np.pi * y) + np.sin(np.pi * y) / np.pi


mu = mobius(N)
print(f"sieved mu to {N:.0e}", flush=True)
for name in ['Q(i)', 'Q(sqrt2)', 'Q(sqrt3)', 'Q(sqrt-3)', 'Q(sqrt5)']:
    rec = D['fields'][name]
    chi = kron(rec['q'], rec['l'], N)
    Mchi = np.cumsum((mu * chi).astype(np.int64))            # Mchi[y]
    # b = mu * (mu chi) via Dirichlet convolution on the full range
    b = np.zeros(N + 1, dtype=np.int64)
    muchi = (mu * chi).astype(np.int64)
    for d in range(1, N + 1):
        if mu[d]:
            b[d::d][: N // d] += int(mu[d]) * muchi[1: N // d + 1]
    MK = np.cumsum(b)
    z = rec['zeros']
    g = np.array([float(r[1]) for r in z])
    a = np.array([float(r[2]) + 1j * float(r[3]) for r in z])
    T = 1000.0
    c = 2 * kern(g / T) * a
    print(f"\n{name}:  M_K(x) at a few x (first values {list(MK[1:12])})")
    print("    t=log x    smoothed x^-1/2 M_K (data)    h*_1000(t) (zeros)")
    for t in [t for t in (4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0) if t < np.log(N) - 0.06] + [np.log(N) - 0.06]:
        # Fourier transform of f(|gamma|/T) is the smoothing kernel; approximate by direct inversion on a grid
        u = np.linspace(t - 0.05, t + 0.05, 2001)
        x = np.exp(u)
        vals = MK[np.minimum(x.astype(np.int64), N)] / np.sqrt(x)
        # kernel K(v) = (1/2pi) int f(|w|/T) e^{-iwv} dw, evaluated numerically
        w = np.linspace(-T, T, 4001)
        Kv = np.real(np.exp(-1j * np.outer(t - u, w)) @ kern(np.abs(w) / T)) * (w[1] - w[0]) / (2 * np.pi)
        sm = np.sum(Kv * vals) * (u[1] - u[0])
        hz = np.real(np.exp(1j * t * g) @ c)
        print(f"    {t:7.3f}    {sm:+.4f}                      {hz:+.4f}")
