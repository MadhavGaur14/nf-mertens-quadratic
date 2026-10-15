"""Step 3: ball enclosure of h*_{K,T}(t) = sum_{|gamma|<T} f(|gamma|/T) e^{i gamma t} / (rho zeta_K'(rho))
from certified isolating intervals of the zeros of zeta and of L(s, chi)."""
from flint import acb, arb
from arb_tools import ZETA, L, L_and_derivative, kernel, nonzero


def coefficients(chi, zeta_mid, L_mid, eps):
    """Per-zero enclosures (gamma ball, a = 1/(rho zeta_K'(rho))). zeta_K = zeta * L(chi)."""
    eps = arb(eps)
    out = []
    for g in zeta_mid:
        rho = acb(arb(1) / 2, arb(g, eps))
        _, zp = L_and_derivative(ZETA, rho)
        dK = zp * L(chi, rho)
        assert nonzero(dK)
        out.append((arb(g, eps), 1 / (rho * dK), 'zeta'))
    for g in L_mid:
        rho = acb(arb(1) / 2, arb(g, eps))
        _, lp = L_and_derivative(chi, rho)
        dK = L(ZETA, rho) * lp
        assert nonzero(dK)
        out.append((arb(g, eps), 1 / (rho * dK), 'L'))
    return out


def hstar(coeffs, T, t):
    """Enclosure of h*_{K,T}(t); coefficients must contain every zero with 0 < gamma < T (conjugates added via 2 Re)."""
    T, t = arb(T), arb(t)
    s = arb(0)
    for g, a, _ in coeffs:
        y = g / T
        s += 2 * kernel(y) * (a * acb(0, g * t).exp()).real
    return s
