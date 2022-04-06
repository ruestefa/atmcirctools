"""Test the basic functionality of GT4Py."""

# Third-party
import gt4py
import gt4py.gtscript as gts
import numpy as np
from gt4py.gtscript import BACKWARD
from gt4py.gtscript import FORWARD
from gt4py.gtscript import PARALLEL
from gt4py.gtscript import computation
from gt4py.gtscript import interval


BACKEND = "gtc:numpy"
DTYPE = np.float32
ORIGIN = (0, 0, 0)


def test_copy() -> None:
    """Copy a field."""

    @gts.stencil(backend=BACKEND)
    def copy(
        in_fld: gts.Field[DTYPE],
        out_fld: gts.Field[DTYPE],
    ) -> None:
        """Copy the input field to the output field."""
        with computation(PARALLEL), interval(...):
            out_fld = in_fld[0, 0, 0]

    # Define fields
    (nx, ny, nz) = shape3d = (4, 5, 6)
    in_data = np.zeros(shape3d, DTYPE)
    for i in range(nx):
        for k in range(nz):
            in_data[i, :, k] = i + k
    out_data = np.zeros(in_data.shape, DTYPE)

    # Define stores
    in_store = gt4py.storage.from_array(
        in_data, backend=BACKEND, default_origin=ORIGIN, dtype=DTYPE
    )
    out_store = gt4py.storage.from_array(
        out_data, backend=BACKEND, default_origin=ORIGIN, dtype=DTYPE
    )

    # Make sure the fields differ
    assert not np.equal(in_store, out_store).all()

    # Run the test: Copy one field to the other
    copy(in_store, out_store)

    # Now the fields are equal
    assert np.equal(in_store, out_store).all()


def test_vertical_shift_periodic() -> None:
    """Shift a field upward twice, whereby the vertical axis acts as periodic."""

    @gts.stencil(backend=BACKEND)
    def shift_up(
        in_fld: gts.Field[gts.IJK, DTYPE],
        wk2d_fld: gts.Field[gts.IJ, DTYPE],
        out_fld: gts.Field[gts.IJK, DTYPE],
    ) -> None:
        """Shift field upward by one, moving the uppermost layer to the bottom."""
        with computation(FORWARD), interval(-1, None):
            wk2d_fld = in_fld[0, 0, 0]
        with computation(FORWARD), interval(0, 1):
            out_fld = wk2d_fld[0, 0]
        with computation(FORWARD), interval(1, None):
            out_fld = in_fld[0, 0, -1]

    # Define fields
    (nx, ny, nz) = shape3d = *shape2d, _ = (4, 5, 6)
    in_data = np.zeros(shape3d, DTYPE)
    for i in range(nx):
        for k in range(nz):
            in_data[i, :, k] = i + k
    out_data = np.zeros(shape3d, DTYPE)

    # Define stores
    in_store = gt4py.storage.from_array(
        in_data, backend=BACKEND, default_origin=ORIGIN, dtype=DTYPE
    )
    wk2d_store = gt4py.storage.empty(
        shape=shape2d, mask="IJ", backend=BACKEND, default_origin=(0, 0), dtype=DTYPE
    )
    out_store = gt4py.storage.from_array(
        out_data, backend=BACKEND, default_origin=ORIGIN, dtype=DTYPE
    )

    # Create test reference and make sure the output array doesn't validate yet
    ref = np.dstack([in_data[:, :, -2:], in_data[:, :, :-2]])
    assert not np.equal(in_store, ref).all()

    # Run test: Shift field upward twice
    shift_up(in_store, wk2d_store, out_store)
    shift_up(out_store, wk2d_store, in_store)

    # Now the output field validates
    assert np.equal(in_store, ref).all()


def test_vert_intp_full() -> None:
    """Interpolate a 3D field in the vertical to a given surface.

    The surface value lies on the grid, making nearest neighbor interpolation
    sufficient.

    """

    @gts.stencil(backend=BACKEND)
    def interpolate(
        fld: gts.Field[gts.IJK, DTYPE],
        grid: gts.Field[gts.IJK, DTYPE],
        val: float,
        intp: gts.Field[gts.IJ, DTYPE],
        vnan: float,
    ) -> None:
        """Interpolate ``fld`` to ``val``-surface of ``grid``."""
        with computation(FORWARD), interval(-1, None):
            intp = vnan
        with computation(FORWARD), interval(...):
            if grid == val:
                intp = fld[0, 0, 0]

    # Define fields
    (nx, ny, nz) = shape3d = *shape2d, _ = (4, 5, 6)
    grid_data = np.zeros(shape3d, DTYPE)
    in_data = np.zeros(shape3d, DTYPE)
    for k in range(nz):
        in_data[:, :, k] = k
        for i in range(nx):
            grid_data[i, :, k] = i + k
    val = 6.0

    # Define stores
    in_store = gt4py.storage.from_array(
        in_data, backend=BACKEND, default_origin=ORIGIN, dtype=DTYPE
    )
    grid_store = gt4py.storage.from_array(
        grid_data, backend=BACKEND, default_origin=ORIGIN, dtype=DTYPE
    )
    intp_store = gt4py.storage.empty(
        shape=shape2d, mask="IJ", backend=BACKEND, default_origin=(0, 0), dtype=DTYPE
    )

    # Create test reference and ensure the output array doesn't validate yet
    idcs = np.where(grid_data < val, 1, 0).sum(axis=2)
    ref = np.full(shape2d, np.nan, DTYPE)
    for i in range(nx):
        for j in range(ny):
            if idcs[i, j] < nz:
                ref[i, j] = in_data[i, j, idcs[i, j]]
    assert not np.equal(intp_store, ref).all()

    # Run test: Interpolate field to surface
    interpolate(in_store, grid_store, val, intp_store, np.nan)

    # Now the output array validates
    assert np.allclose(intp_store.data, ref, equal_nan=True)


def test_vert_intp_between() -> None:
    """Interpolate a 3D field in the vertical to a given surface.

    The surface value lies between the grid, requiring interpolation between
    vertical grid points.

    """
    NAN = np.nan

    @gts.stencil(backend=BACKEND)
    def interpolate(
        fld: gts.Field[gts.IJK, DTYPE],
        grid: gts.Field[gts.IJK, DTYPE],
        val: float,
        intp: gts.Field[gts.IJ, DTYPE],
        wk_lower: gts.Field[gts.IJ, DTYPE],
        wk_upper: gts.Field[gts.IJ, DTYPE],
        d_lower: gts.Field[gts.IJ, DTYPE],
        d_upper: gts.Field[gts.IJ, DTYPE],
        vnan: float,
    ) -> None:
        """Interpolate ``fld`` to ``val``-surface of ``grid``."""
        with computation(FORWARD), interval(...):
            wk_lower = vnan
            d_lower = 0.0
        with computation(FORWARD), interval(...):
            if grid < val:
                wk_lower = fld[0, 0, 0]
                d_lower = val - grid[0, 0, 0]
        with computation(BACKWARD), interval(...):
            wk_upper = vnan
            d_upper = 0.0
        with computation(BACKWARD), interval(...):
            if grid > val:
                wk_upper = fld[0, 0, 0]
                d_upper = grid[0, 0, 0] - val
        with computation(FORWARD), interval(...):
            if wk_lower[0, 0] == vnan or wk_upper[0, 0] == vnan:
                intp = vnan
            else:
                intp = (
                    d_lower[0, 0] / (d_lower[0, 0] + d_upper[0, 0]) * wk_lower[0, 0]
                    + d_upper[0, 0] / (d_lower[0, 0] + d_upper[0, 0]) * wk_upper[0, 0]
                )

    # Define fields
    (nx, ny, nz) = shape3d = *shape2d, _ = (4, 5, 6)
    grid_data = np.zeros(shape3d, DTYPE)
    fld_data = np.zeros(shape3d, DTYPE)
    for k in range(nz):
        fld_data[:, :, k] = k
        for i in range(nx):
            grid_data[i, :, k] = i + k
    val = 5.2

    # Define stores
    in_store = gt4py.storage.from_array(
        fld_data, backend=BACKEND, default_origin=ORIGIN, dtype=DTYPE
    )
    grid_store = gt4py.storage.from_array(
        grid_data, backend=BACKEND, default_origin=ORIGIN, dtype=DTYPE
    )
    intp_store = gt4py.storage.empty(
        shape=shape2d, mask="IJ", backend=BACKEND, default_origin=(0, 0), dtype=DTYPE
    )
    wk1_store = gt4py.storage.empty(
        shape=shape2d, mask="IJ", backend=BACKEND, default_origin=(0, 0), dtype=DTYPE
    )
    wk2_store = gt4py.storage.empty(
        shape=shape2d, mask="IJ", backend=BACKEND, default_origin=(0, 0), dtype=DTYPE
    )
    wk3_store = gt4py.storage.empty(
        shape=shape2d, mask="IJ", backend=BACKEND, default_origin=(0, 0), dtype=DTYPE
    )
    wk4_store = gt4py.storage.empty(
        shape=shape2d, mask="IJ", backend=BACKEND, default_origin=(0, 0), dtype=DTYPE
    )

    # Create test reference and ensure the output array doesn't validate yet
    idcs_lower = (np.where(grid_data < val, 1, 0).sum(axis=2) - 1)
    idcs_upper = idcs_lower + 1
    idcs_lower[idcs_lower > nz - 1] = -1
    idcs_upper[idcs_upper > nz - 1] = -1
    # grid_lower = np.take(grid_data, idcs_lower)
    # grid_upper = np.take(grid_data, idcs_upper)
    # ref_lower = np.take(fld_data, idcs_lower)
    # ref_upper = np.take(fld_data, idcs_upper)
    grid_lower = np.full(shape2d, np.nan, DTYPE)
    grid_upper = np.full(shape2d, np.nan, DTYPE)
    ref_lower = np.full(shape2d, np.nan, DTYPE)
    ref_upper = np.full(shape2d, np.nan, DTYPE)
    for i in range(nx):
        for j in range(ny):
            if idcs_lower[i, j] > 0:
                grid_lower[i, j] = grid_data[i, j, idcs_lower[i, j]]
                ref_lower[i, j] = fld_data[i, j, idcs_lower[i, j]]
            if idcs_upper[i, j] > 0:
                grid_upper[i, j] = grid_data[i, j, idcs_upper[i, j]]
                ref_upper[i, j] = fld_data[i, j, idcs_upper[i, j]]
    d_lower = val - grid_lower
    d_upper = grid_upper - val
    d_tot = d_lower + d_upper
    ref = d_lower / d_tot * ref_lower + d_upper / d_tot * ref_upper
    assert not np.equal(intp_store, ref).all()

    # Run test: Interpolate field to surface
    interpolate(
        in_store,
        grid_store,
        val,
        intp_store,
        wk1_store,
        wk2_store,
        wk3_store,
        wk4_store,
        np.nan,
    )

    # Now the output array validates
    assert np.allclose(intp_store.data, ref, equal_nan=True)
