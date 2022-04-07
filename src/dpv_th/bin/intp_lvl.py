"""Interpolate a 3D field to a 2D surface at a given level."""
from __future__ import annotations

# Standard library
import dataclasses as dc
import sys
from typing import Any

# Third-party
import click
import numpy as np
import numpy.typing as npt
import xarray as xr
from atmcirclib.typing import PathLike_T

# First-party
from dpv_th.intp import LevelInterpolator


def write(intp_flds: dict[str, npt.NDArray[np.float_]], config: Config) -> None:
    """Write interpolated field to disk."""
    ds_in = xr.open_dataset(config.infile)
    grid_var = ds_in.variables[config.grid_var]
    dims = grid_var.dims
    assert dims == ("time", "lev", "lat", "lon"), str(dims)
    da_lev = xr.DataArray(
        name="lev",
        data=[config.lvl],
        coords={"lev": [config.lvl]},
        attrs={
            "long_name": f"{config.grid_var} levels",
            "units": ds_in.variables[config.grid_var].attrs["units"],
        },
    )
    coords = {
        "time": ds_in.time,
        "lev": da_lev,
        "lat": ds_in.lat,
        "lon": ds_in.lon,
    }
    data_vars: dict[str, xr.DataArray] = {}
    fld_var = ds_in.variables[config.fld_var]
    assert fld_var.dims == dims, str(fld_var.dims)
    for name, intp_fld in intp_flds.items():
        assert len(intp_fld.shape) == 2, str(intp_fld.shape)
        da_fld = xr.DataArray(
            name=name,
            data=intp_fld[None, None, :, :],
            dims=dims,
            attrs=fld_var.attrs,
        )
        data_vars[name] = da_fld
    ds = xr.Dataset(
        coords=coords,
        data_vars=data_vars,
    )
    ds.to_netcdf(config.outfile)


@dc.dataclass
class Config:
    """Script configuration."""

    infile: PathLike_T
    outfile: PathLike_T
    grid_var: str
    fld_var: str
    lvl: float


@click.command(
    context_settings={"show_default": True, "help_option_names": ["-h", "--help"]},
    help="Interpolate a 3D field to a 2D surface at a given level",
)
@click.option(
    "-i",
    "--infile",
    help="Input file containing both grid and field variables.",
    type=click.Path(exists=True),
    required=True,
)
@click.option(
    "-o",
    "--outfile",
    help="Output file.",
    type=click.Path(exists=False),
    required=True,
)
@click.option(
    "-g",
    "--grid-var",
    help="Name of grid variable.",
    default="TH",
)
@click.option(
    "-f",
    "--fld-var",
    help="Name of field variable to be interpolated.",
    default="PV",
)
@click.option(
    "-l",
    "--lvl",
    help="Level to which to interpolate field.",
    type=float,
    default=320,
)
def cli(**config_kwargs: Any) -> int:
    """Entry point from command line."""
    config = Config(**config_kwargs)
    print(f"read {config.grid_var}, {config.fld_var} from {config.infile}")
    with xr.open_dataset(config.infile) as ds:
        dims = ds.variables[config.grid_var].dims
        assert dims == ("time", "lev", "lat", "lon"), f"unexpected_dims: {dims}"
        grid = ds.variables[config.grid_var].data[0, :, :, :]
        fld = ds.variables[config.fld_var].data[0, :, :, :]
    grid = np.moveaxis(grid, 0, 2)
    fld = np.moveaxis(fld, 0, 2)
    print(f"interpolate {config.fld_var} to {config.grid_var} at {config.lvl}")
    intp_fld_up = LevelInterpolator(grid, direction="up").to_level(fld, config.lvl)
    print(f"interpolate {config.fld_var} to {config.grid_var} at {config.lvl}")
    intp_fld_down = LevelInterpolator(grid, direction="down").to_level(fld, config.lvl)
    intp_flds = {
        f"{config.fld_var}_up": intp_fld_up,
        f"{config.fld_var}_down": intp_fld_down,
    }
    print(f"write {', '.join(intp_flds)} to {config.outfile}")
    write(intp_flds, config)
    return 0


if __name__ == "__main__":
    sys.exit(cli())
