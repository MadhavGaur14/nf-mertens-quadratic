"""Shared definitions for the nf_mertens project: paths, fields, kernel, data loading."""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / 'data' / 'zeros'
RESULTS = ROOT / 'results'
CERTS = ROOT / 'certificates'

# Quadratic field -> real primitive Dirichlet character (modulus, Conrey index) in python-flint's dirichlet_char.
# zeta_K(s) = zeta(s) L(s, chi_D) with chi_D the Kronecker symbol of the field discriminant D.
FIELDS = {'Q(i)': (4, 3), 'Q(sqrt2)': (8, 5), 'Q(sqrt3)': (12, 11), 'Q(sqrt-3)': (3, 2), 'Q(sqrt5)': (5, 4),
          'Q(sqrt-2)': (8, 3)}

DEFAULT_ZEROS = DATA / 'quadratic_T1000.json'


def load_fields(path=DEFAULT_ZEROS):
    """Zero records per field: (source 'zeta'|'L', gamma, Re a, Im a, |derivative|), a = 1/(rho zeta_K'(rho))."""
    fields = json.load(open(path))['fields']
    extra = DATA / 'Qsqrt-2_T700.json'
    if path == DEFAULT_ZEROS and extra.exists():
        fields.update(json.load(open(extra))['fields'])
    return fields


def kern_np(y):
    """Jurkat-Peyerimhoff kernel f(y) = (1-y)cos(pi y) + sin(pi y)/pi on [0,1] (numpy)."""
    import numpy as np
    return (1 - y) * np.cos(np.pi * y) + np.sin(np.pi * y) / np.pi
