"""Subpackage ``dpv_th.deriv``."""
# Local
from ._deriv_f77 import aura as aura_f77
from ._deriv_f77 import deriv as deriv_f77
from ._deriv_f77 import div as div_f77
from ._deriv_f77 import grad as grad_f77
from ._deriv_f77 import gridok as gridok_f77
from ._deriv_f77 import rot as rot_f77

__all__: list[str] = [
    "aura_f77",
    "deriv_f77",
    "div_f77",
    "grad_f77",
    "gridok_f77",
    "rot_f77",
]
