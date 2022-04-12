"""Procedures defined in extension module ``dpv_th.deriv._deriv_f77``."""
# Local
from ._deriv_f77 import aura
from ._deriv_f77 import deriv
from ._deriv_f77 import div
from ._deriv_f77 import grad
from ._deriv_f77 import gridok
from ._deriv_f77 import rot

__all__: list[str] = [
    "aura",
    "deriv",
    "div",
    "grad",
    "gridok",
    "rot",
]
