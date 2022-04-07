"""Interpolate 3D field to a 2D surface at a given leven.

This module uses gt4py, hence no ``from __future__ import annotations``.

"""

# Standard library
from enum import Enum
from typing import Any
from typing import Optional
from typing import Union

# Third-party
import gt4py.gtscript as gts
import gt4py.storage as gt_store
import numpy as np
import numpy.typing as npt
from gt4py.gtscript import BACKWARD
from gt4py.gtscript import computation
from gt4py.gtscript import FORWARD
from gt4py.gtscript import interval
from gt4py.gtscript import isnan

# from gt4py.gtscript import PARALLEL
# from gt4py.gtscript import sqrt


class StrEnum(str, Enum):
    """Enum that compares with string values."""

    def __str__(self) -> str:
        """Return value as string."""
        return str(self.value)


class InterpolationDirection(StrEnum):
    """Direction in which vertical interpolation is conducted."""

    up = "up"
    down = "down"


InterpolationDirectionLike_T = Union[InterpolationDirection, str]


class LevelInterpolatorError(Exception):
    """Base exception for error during level interpolation."""


class NonMonotonicGridError(LevelInterpolatorError):
    """Grid is not monotonically in-/decreasing in z."""


class LevelInterpolator:
    """Interpolate 3D fields to 2D surfaces at given levels."""

    def __init__(
        self,
        grid: npt.ArrayLike,
        *,
        direction: Optional[InterpolationDirectionLike_T] = None,
        dtype: npt.DTypeLike = np.float32,
        gt_backend: str = "gtc:numpy"
    ) -> None:
        """Create a new instance.

        Args:
            grid: 3D grid array, ideally with monotonically increasing or
                decreasing values.

            direction (optional): Direction in which vertical interpolation
                is conducted; typically 'forward'/'backward' for monotonically
                in-/decreasing grids in z.

            dtype (optional): Data type used internally during computations.

            gt_backend (optional): Backend used by GT4Py.

        """
        self.grid: npt.NDArray[np.float_] = np.asarray(grid, dtype)
        if direction is None:
            try:
                direction = self._derive_direction()
            except NonMonotonicGridError as e:
                raise ValueError(
                    "must specify direction for grid that is not non-monotonically"
                    " in-/decreasing grid in z"
                ) from e
        self.direction: InterpolationDirection = InterpolationDirection(direction)
        self.dtype: npt.DTypeLike = dtype
        self.gt_backend: str = gt_backend

    def to_level(self, fld: npt.ArrayLike, lvl: float) -> npt.NDArray[np.float_]:
        """Interpolate to a single level."""
        dtype_in = fld.dtype if isinstance(fld, np.ndarray) else self.dtype
        fld = np.asarray(fld, self.dtype)
        grid = self.grid[:, :, ::-1] if str(self.direction) == "down" else self.grid
        DTYPE = self.dtype

        @gts.stencil(backend=self.gt_backend)
        def to_level(
            fld: gts.Field[gts.IJK, DTYPE],
            grid: gts.Field[gts.IJK, DTYPE],
            lvl: float,
            intp: gts.Field[gts.IJ, DTYPE],
            fld_below: gts.Field[gts.IJ, DTYPE],
            fld_above: gts.Field[gts.IJ, DTYPE],
            d_below: gts.Field[gts.IJ, DTYPE],
            d_above: gts.Field[gts.IJ, DTYPE],
            vnan: float,
        ) -> None:
            """Interpolate ``fld`` to ``val``-surface of ``grid``."""
            with computation(FORWARD), interval(...):
                fld_below[...] = vnan
                d_below[...] = 0.0
            with computation(FORWARD), interval(...):
                if grid < lvl:
                    fld_below[...] = fld[0, 0, 0]
                    d_below[...] = lvl - grid[0, 0, 0]
            with computation(BACKWARD), interval(...):
                fld_above[...] = vnan
                d_above[...] = 0.0
            with computation(BACKWARD), interval(...):
                if grid > lvl:
                    fld_above[...] = fld[0, 0, 0]
                    d_above[...] = grid[0, 0, 0] - lvl
            with computation(FORWARD), interval(...):
                if (
                    isnan(vnan)
                    and (isnan(fld_below[0, 0]) == vnan or isnan(fld_above[0, 0]))
                ) or (fld_below[0, 0] == vnan or fld_above[0, 0] == vnan):
                    intp[...] = vnan
                else:
                    intp[...] = fld_below[0, 0] + d_below[0, 0] / (
                        d_below[0, 0] + d_above[0, 0]
                    )

        # Define stores
        kwargs3d: dict[str, Any] = {
            "default_origin": (0, 0, 0),
            "dtype": self.dtype,
            "backend": self.gt_backend,
        }
        kwargs2d: dict[str, Any] = {
            **kwargs3d,
            "default_origin": (0, 0),
            "mask": "IJ",
        }
        shape2d = fld.shape[:2]
        fld_store = gt_store.from_array(fld, **kwargs3d)
        grid_store = gt_store.from_array(grid, **kwargs3d)
        intp_store = gt_store.empty(shape=shape2d, **kwargs2d)
        wk_fld_below_store = gt_store.empty(shape=shape2d, **kwargs2d)
        wk_fld_above_store = gt_store.empty(shape=shape2d, **kwargs2d)
        wk_d_below_store = gt_store.empty(shape=shape2d, **kwargs2d)
        wk_d_above_store = gt_store.empty(shape=shape2d, **kwargs2d)

        to_level(
            fld_store,
            grid_store,
            lvl,
            intp_store,
            wk_fld_below_store,
            wk_fld_above_store,
            wk_d_below_store,
            wk_d_above_store,
            np.nan,
        )

        return np.asarray(intp_store.data, dtype_in)

    def _derive_direction(self) -> InterpolationDirection:
        """Derive interpol. direction from monotonically in-/decreasing grid."""
        deltas = self.grid[:, :, 1:] - self.grid[:, :, :-1]
        if (deltas > 0).all():
            return InterpolationDirection.up
        elif (deltas < 0).all():
            return InterpolationDirection.down
        else:
            raise NonMonotonicGridError(
                "cannot derive direction from non-monotonic grid in z"
            )
