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
        """Copy the field."""
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
        shape=(nx, ny), mask="IJ", backend=BACKEND, default_origin=(0, 0), dtype=DTYPE
    )
    out_store = gt4py.storage.from_array(
        out_data, backend=BACKEND, default_origin=ORIGIN, dtype=DTYPE
    )

    # Create test reference and make sure it doesn't validate yet
    ref = np.dstack([in_data[:, :, -2:], in_data[:, :, :-2]])
    assert not np.equal(in_store, ref).all()

    # Run test: Shift field upward twice
    shift_up(in_store, wk2d_store, out_store)
    shift_up(out_store, wk2d_store, in_store)

    # Now the field validates against the reference
    assert np.equal(in_store, ref).all()


def test_vert_intp_full() -> None:
    """Interpolate a 3D field in the vertical to a given surface.

    The surface value lies on the grid, making nearest neighbor interpolation
    sufficient.

    """

    @gts.stencil(backend=BACKEND)
    def interpolate(
        in_fld: gts.Field[gts.IJK, DTYPE],
        grid_fld: gts.Field[gts.IJK, DTYPE],
        val: float,
        intp_fld: gts.Field[gts.IJ, DTYPE],
    ) -> None:
        """Shift field upward by one, moving the uppermost layer to the bottom."""
        with computation(FORWARD), interval(-1, None):
            intp_fld = in_fld[0, 0, 0]
        with computation(FORWARD), interval(...):
            if grid_fld == val:
                intp_fld = in_fld[0, 0, 0]

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
        shape=(nx, ny), mask="IJ", backend=BACKEND, default_origin=(0, 0), dtype=DTYPE
    )

    # Create test reference and ensure the output array doesn't validate yet
    ref = np.take(in_data, np.where(grid_data < val, 1, 0).sum(axis=2).clip(max=nz - 1))
    assert not np.equal(intp_store, ref).all()

    # Run test: Interpolate field to surface
    interpolate(in_store, grid_store, val, intp_store)

    # Now the output array validates
    assert np.equal(intp_store, ref).all()
