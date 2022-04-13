"""Subpackage ``dpv_th.deriv``."""
# Local
from . import f77
from . import f90
from .f77 import aura as aura_f77
from .f77 import deriv as deriv_f77
from .f77 import div as div_f77
from .f77 import grad as grad_f77
from .f77 import gridok as gridok_f77
from .f77 import rot as rot_f77
from .f90 import aura as aura_f90
from .f90 import deformation as deformation_f90
from .f90 import deriv as deriv_f90
from .f90 import deriv_time as deriv_time_f90
from .f90 import div as div_f90
from .f90 import gradmat as gradmat_f90
from .f90 import gradvec as gradvec_f90
from .f90 import hadv as hadv_f90
from .f90 import rot as rot_f90

__all__: list[str] = [
    "f77",
    "f90",
    "aura_f77",
    "deriv_f77",
    "div_f77",
    "grad_f77",
    "gridok_f77",
    "rot_f77",
    "aura_f90",
    "deformation_f90",
    "deriv_f90",
    "deriv_time_f90",
    "div_f90",
    "gradmat_f90",
    "gradvec_f90",
    "hadv_f90",
    "rot_f90",
]
