"""Step 2: certify a sign change of the Hardy Z-function across [gamma - eps, gamma + eps] for every approximate zero.
Combined with Step 1 (count equality) this proves all zeros below T are on the critical line, simple, and isolated."""
from flint import arb, ctx
from arb_tools import hardy_z, positive, negative


def refine(chi, g, eps):
    """Illinois refinement of an approximate zero (arb midpoint) to well below eps."""
    a, b = g - 100 * eps, g + 100 * eps
    fa, fb = hardy_z(chi, a).mid(), hardy_z(chi, b).mid()
    if (fa > 0) == (fb > 0):
        return None
    for _ in range(200):
        c = (b - fb * (b - a) / (fb - fa)).mid()
        fc = hardy_z(chi, c).mid()
        step = abs(c - b)
        if (fc > 0) != (fb > 0):
            a, fa = b, fb
        else:
            fa = fa / 2
        b, fb = c, fc
        if step < eps / 1000:
            break
    return b


def sign_change(chi, g, eps):
    zl, zr = hardy_z(chi, g - eps), hardy_z(chi, g + eps)
    return (positive(zl) and negative(zr)) or (negative(zl) and positive(zr))


def isolate(chi, gammas, eps):
    """gammas: approximate ordinates (decimal strings). Returns certified midpoints (arb, exact) such that Z changes sign
    on [g - eps, g + eps]; raises if any zero cannot be certified."""
    eps = arb(eps)
    out = []
    for s in gammas:
        g = arb(s).mid()
        if not sign_change(chi, g, eps):
            g2 = refine(chi, g, eps)
            if g2 is None or not sign_change(chi, g2, eps):
                raise RuntimeError(f"could not isolate zero near {s} (prec {ctx.prec})")
            g = g2
        out.append(g)
    out.sort(key=lambda x: float(x))
    for x, y in zip(out, out[1:]):                       # disjointness
        if not (y - x - 2 * eps).lower() > 0:
            raise RuntimeError(f"isolating intervals overlap near {x}")
    return out
