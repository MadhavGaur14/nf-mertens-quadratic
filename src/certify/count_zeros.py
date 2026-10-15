"""Step 1: rigorous count N(T) = #{rho : 0 < gamma < T} of nontrivial zeros of L(s, chi) (or zeta) by the argument
principle along 2 -> 2+iT -> 1/2+iT. See docs/04_certification_plan.md."""
from flint import acb, arb, ctx
from arb_tools import theta, L, nonzero, positive, is_zeta


# L(sigma, chi) > 0 for all sigma > 0 is proved analytically (ball evaluation of L via Hurwitz zeta returns NaN on any
# ball containing the Hurwitz pole sigma = 1, so a bisection certificate is not available there):
#   chi_{-3}: L = sum_k [(3k+1)^-s - (3k+2)^-s], every bracket > 0.
#   chi_5:    L = sum_k [f(1) - f(2) - f(3) + f(4)], f(u) = (5k+u)^-s decreasing and convex, so f(1)-f(2) > f(3)-f(4).
REAL_AXIS_LEMMA = {(3, 2): 'L(s,chi_-3) > 0 for s > 0: pairs (3k+1)^-s - (3k+2)^-s > 0',
                   (5, 4): 'L(s,chi_5) > 0 for s > 0: blocks f(1)-f(2)-f(3)+f(4) > 0 by convexity of (5k+u)^-s',
                   (4, 3): 'L(s,chi_-4) > 0 for s > 0: pairs (4k+1)^-s - (4k+3)^-s > 0',
                   (8, 5): 'L(s,chi_8) > 0 for s > 0: blocks f(1)-f(3)-f(5)+f(7) > 0 by convexity of (8k+u)^-s',
                   (12, 11): 'L(s,chi_12) > 0 for s > 0: blocks f(1)-f(5)-f(7)+f(11) > 0 by convexity of (12k+u)^-s'}


def real_axis_positive(chi):
    key = (chi.modulus(), chi.number())
    if key not in REAL_AXIS_LEMMA:
        raise RuntimeError(f"no real-axis positivity lemma recorded for character {key}")
    for sg in ('0.5', '0.75', '0.999', '1.001', '1.5', '2'):     # numerical sanity check only
        assert positive(L(chi, acb(arb(sg))).real)
    return REAL_AXIS_LEMMA[key]


def arg_change_horizontal(chi, T, s_from=2.0, s_to=0.5, max_depth=60):
    """Continuous change of arg L(sigma + iT) as sigma goes from s_from to s_to, certified piecewise."""
    T = arb(T)
    total, pieces = arb(0), 0
    stack = [(arb(s_from), arb(s_to), 0)]
    while stack:
        a, b, d = stack.pop()
        ball = acb(arb((a + b) / 2, abs(b - a) / 2), T)
        enc = L(chi, ball)
        if nonzero(enc):
            dv = (L(chi, acb(b, T)) / L(chi, acb(a, T))).arg()
            if dv.rad() < 0.25:
                total += dv
                pieces += 1
                continue
        if d > max_depth:
            raise RuntimeError(f"arg tracking failed near sigma={a}")
        m = (a + b) / 2
        stack += [(m, b, d + 1), (a, m, d + 1)]          # order irrelevant: increments are summed
    return total, pieces


def count_zeros(chi, T):
    """Returns (N, info) with N a certified integer."""
    T = arb(T)
    if not nonzero(L(chi, acb(arb(1) / 2, T))):
        raise RuntimeError("T may be a zero ordinate; choose another T")
    if not is_zeta(chi):
        real_pieces = real_axis_positive(chi)
    else:
        real_pieces = 'classical (zeta(sigma) < 0 on (0,1))'
    v = L(chi, acb(2, T))
    assert positive(v.real)                                # Re L(2+it) >= 2 - zeta(2) > 0
    vert = v.arg()
    horiz, pieces = arg_change_horizontal(chi, T)
    extra = 1 if is_zeta(chi) else 0                       # factor s(s-1) of xi(s) contributes +pi
    Nb = (theta(chi, T) + vert + horiz) / arb.pi() + extra
    n = int(round(float(Nb.mid())))
    if not (abs(Nb - n).upper() < 0.5):
        raise RuntimeError(f"count not certified: {Nb}")
    return n, {'N_ball': Nb.str(15), 'horizontal_pieces': pieces, 'real_axis_pieces': real_pieces,
               'prec': ctx.prec}


if __name__ == '__main__':
    import sys
    from flint import dirichlet_char
    ctx.prec = 128
    q, l, T = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
    print(count_zeros(dirichlet_char(q, l), T))
