"""Test the basic functionality of GT4Py."""

# Third-party
import gt4py
import gt4py.gtscript as gts
import numpy as np
from gt4py.gtscript import computation
from gt4py.gtscript import interval
from gt4py.gtscript import PARALLEL


BACKEND = "gtc:numpy"
DTYPE = np.float32
ORIGIN = (0, 0, 0)


def test_copy() -> None:
    """Copy a field."""

    @gts.stencil(backend=BACKEND)
    def copy(
        in_fld: gts.Field[DTYPE],
        out_fld: gts.Field[DTYPE],
    ):
        with computation(PARALLEL), interval(...):
            out_fld = in_fld[0, 0, 0]

    (nx, ny, nz) = shape = (5, 5, 5)
    fld = np.zeros(shape, DTYPE)
    for i in range(nx):
        for k in range(nz):
            fld[i, :, k] = i + k

    in_data = fld
    out_data = np.zeros(fld.shape, DTYPE)

    in_store = gt4py.storage.from_array(
        in_data, BACKEND, default_origin=ORIGIN, dtype=DTYPE
    )
    out_store = gt4py.storage.from_array(
        out_data, BACKEND, default_origin=ORIGIN, dtype=DTYPE
    )

    assert not np.equal(in_store, out_store).all()
    copy(in_store, out_store)
    assert np.equal(in_store, out_store).all()
