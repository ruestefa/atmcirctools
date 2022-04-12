"""Procedures defined in extension module ``dpv_th.deriv._deriv_f95``."""
# Local
from ._deriv_f90 import aura
from ._deriv_f90 import deformation
from ._deriv_f90 import deriv
from ._deriv_f90 import deriv_time
from ._deriv_f90 import div
from ._deriv_f90 import grad
from ._deriv_f90 import gradmat
from ._deriv_f90 import gradvec
from ._deriv_f90 import hadv
from ._deriv_f90 import rot

__all__: list[str] = [
    "aura",
    "deformation",
    "deriv",
    "deriv_time",
    "div",
    "grad",
    "gradmat",
    "gradvec",
    "hadv",
    "rot",
]
