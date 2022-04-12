"""Top-level package ``dpv_th``."""
# Standard library
import importlib.metadata

# Local
from . import deriv

metadata = importlib.metadata.metadata(__package__)

__author__ = "Stefan Ruedisuehli"
__email__ = "stefan.ruedisuehli@env.ethz.ch"
__version__ = importlib.metadata.version(__package__)

del importlib

__all__: list[str] = [
    "deriv",
]
