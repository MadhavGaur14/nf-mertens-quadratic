"""Ball-arithmetic helpers (python-flint / Arb). Every function returns enclosures; 'certain' predicates are True only
when the property is proved for every point of the enclosing balls."""
from flint import acb, arb, acb_series, dirichlet_char, ctx

ZETA = dirichlet_char(1, 1)


def is_zeta(chi):
    return chi.modulus() == 1


def positive(x):
    """x is an arb: True iff every point of the ball is > 0."""
    return x.lower() > 0


def negative(x):
    return x.upper() < 0


def nonzero(z):
    """z is an acb: True iff the ball excludes 0."""
    return abs(z).lower() > 0


def theta(chi, T):
    """theta(T) = (T/2) log(q/pi) + Im log Gamma((1/2 + a + iT)/2), q = conductor, a = parity (principal branch of
    log Gamma, continuous on the path from the real axis)."""
    q, a = chi.conductor(), chi.parity()
    T = arb(T)
    return T / 2 * (arb(q) / arb.pi()).log() + acb((arb(1) / 2 + a) / 2, T / 2).lgamma().imag


def L(chi, s):
    return chi.l(s)


def L_and_derivative(chi, rho):
    """Enclosures of L(rho) and L'(rho) over the whole ball rho, by Taylor-series arithmetic."""
    c = acb_series([rho, 1], prec=2).dirichlet_l(chi).coeffs()
    return c[0], c[1]


def hardy_z(chi, t):
    z = chi.hardy_z(t)
    return z.real


def kernel(y):
    """Jurkat-Peyerimhoff kernel f(y) = (1-y) cos(pi y) + sin(pi y)/pi, y an arb in [0,1]."""
    py = arb.pi() * y
    return (1 - y) * py.cos() + py.sin() / arb.pi()


def versions():
    import flint, sys
    return {'python': sys.version.split()[0], 'python-flint': flint.__version__}
