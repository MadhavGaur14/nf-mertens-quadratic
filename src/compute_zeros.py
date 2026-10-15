"""Zeros of zeta_K = zeta * L(s, chi_D) for quadratic fields K and the coefficients a_rho = 1/(rho zeta_K'(rho)).
Non-rigorous (sign-change scan + Illinois refinement); certified versions live in src/certify/.
usage: python src/compute_zeros.py TMAX [FIELD ...]   -> writes data/zeros/quadratic_T{TMAX}.json"""
import json, math, sys, time
from flint import dirichlet_char, acb, arb, ctx
from common import FIELDS as ALL_FIELDS, DATA

ZETA = dirichlet_char(1, 1)
TMAX = float(sys.argv[1]) if len(sys.argv) > 1 else 1000.0
FIELDS = {k: ALL_FIELDS[k] for k in (sys.argv[2:] or ['Q(i)', 'Q(sqrt2)', 'Q(sqrt3)', 'Q(sqrt-3)', 'Q(sqrt5)'])}
STEP = 0.02


def refine(chi, a, b):
    """Illinois iteration on Hardy Z of chi in [a, b] (sign change), 256-bit."""
    ctx.prec = 256
    a, b = arb(a), arb(b)
    fa, fb = chi.hardy_z(a).real.mid(), chi.hardy_z(b).real.mid()
    for _ in range(200):
        c = (b - fb * (b - a) / (fb - fa)).mid()
        fc = chi.hardy_z(c).real.mid()
        step = abs(c - b)
        if (fc > 0) != (fb > 0):
            a, fa = b, fb
        else:
            fa = fa / 2
        b, fb = c, fc
        if step < arb(10) ** -40 or fc == 0:
            break
    return b


def theta(chi, t):
    """Phase of the completed L-function: Z(t) = e^{i theta(t)} L(1/2+it) (real primitive chi)."""
    ctx.prec = 128
    q, par = chi.conductor(), chi.parity()
    t = arb(t)
    return (t / 2 * (arb(q) / arb.pi()).log() + acb(arb(0.5 + par) / 2, t / 2).lgamma().imag)


def L_zeros(chi, T):
    ctx.prec = 96
    zs, n, prev = [], int(T / STEP) + 1, None
    t_prev, z_prev = None, None
    for k in range(1, n + 1):
        t = k * STEP
        z = float(chi.hardy_z(arb(t)).real.mid())
        if z_prev is not None and (z > 0) != (z_prev > 0):
            zs.append(refine(chi, t_prev, t))
            ctx.prec = 96
        t_prev, z_prev = t, z
    return zs


def deriv(chi, g):
    """d/ds L(s, chi) at s = 1/2 + i g by a central difference along the critical line (400-bit)."""
    ctx.prec = 400
    h = arb(10) ** -30
    f = lambda tt: chi.l(acb(arb(1) / 2, tt))
    return (f(g + h) - f(g - h)) / (2 * h) / acb(0, 1)


def main():
    out = {'TMAX': TMAX, 'fields': {}}
    t0 = time.time()
    ctx.prec = 256
    nest = int(TMAX / (2 * math.pi) * math.log(TMAX / (2 * math.pi * math.e)) + 60)
    zz = [z.imag.mid() for z in acb.zeta_zeros(1, nest) if z.imag < TMAX]   # rigorous (Arb) zeta zeros
    print(f"zeta: {len(zz)} zeros below {TMAX} ({time.time() - t0:.0f}s)", flush=True)
    zeta_der = [deriv(ZETA, g) for g in zz]
    for name, (q, l) in FIELDS.items():
        chi = dirichlet_char(q, l)
        t1 = time.time()
        lz = L_zeros(chi, TMAX)
        # completeness heuristic: count vs theta(T)/pi
        checks = []
        for T in [TMAX * f for f in (0.25, 0.5, 0.75, 1.0)]:
            cnt = sum(1 for g in lz if g < T)
            checks.append((T, cnt, float((theta(chi, T) / arb.pi()).mid())))
        sep = min(abs(float(a.mid()) - float(b.mid())) for a in lz for b in zz)
        ctx.prec = 400
        Lhalf = chi.l(acb(0.5))
        recs = []
        for g, zd in zip(zz, zeta_der):          # zeros of zeta: zeta_K' = zeta' * L
            rho = acb(arb(1) / 2, g)
            dK = zd * chi.l(rho)
            a = 1 / (rho * dK)
            recs.append(('zeta', g.str(35, radius=False), a.real.str(30, radius=False), a.imag.str(30, radius=False),
                         float(abs(zd).mid())))
        for g in lz:                              # zeros of L: zeta_K' = zeta * L'
            rho = acb(arb(1) / 2, g)
            ld = deriv(chi, g)
            dK = ZETA.l(rho) * ld
            a = 1 / (rho * dK)
            recs.append(('L', g.str(35, radius=False), a.real.str(30, radius=False), a.imag.str(30, radius=False),
                         float(abs(ld).mid())))
        recs.sort(key=lambda r: float(r[1]))
        mind = min(r[4] for r in recs)
        out['fields'][name] = {'q': q, 'l': l, 'L_half': Lhalf.real.str(20, radius=False), 'n_zeta': len(zz),
                               'n_L': len(lz), 'min_sep_zeta_L': sep, 'min_abs_derivative': mind,
                               'count_checks': checks, 'zeros': recs}
        print(f"{name:10s} chi=({q},{l}) L(1/2)={Lhalf.real.str(8, radius=False)}  L-zeros:{len(lz)}"
              f"  min|gamma_zeta-gamma_L|={sep:.3g}  min|deriv|={mind:.3g}  ({time.time() - t1:.0f}s)", flush=True)
        for T, cnt, th in checks:
            print(f"     T={T:7.1f}: found {cnt:5d} sign changes, theta(T)/pi = {th:9.2f}, diff {cnt - th:+.2f}")
    json.dump(out, open(DATA / f'quadratic_T{int(TMAX)}.json', 'w'))
    print(f"done in {time.time() - t0:.0f}s")


if __name__ == '__main__':
    main()
