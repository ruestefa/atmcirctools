"""Interpolate a 3D field to a 2D surface at a given level."""
from __future__ import annotations

# Standard library
import dataclasses as dc
import sys
from collections.abc import Sequence
from typing import Any

# Third-party
import click
import numpy as np
import numpy.typing as npt
import xarray as xr
from atmcirclib.intp import LevelInterpolator
from atmcirclib.typing import PathLike_T


def check_dims(dims: tuple[str, ...], xp_dims: tuple[str, ...]) -> None:
    """Check dimension names."""
    if dims != xp_dims:
        raise Exception(f"expected dims {xp_dims}, got {dims}")


def write(intp_flds: dict[str, npt.NDArray[np.float_]], config: Config) -> None:
    """Write interpolated field to disk."""
    ds_in = xr.open_dataset(config.infile)
    grid_var = ds_in.variables[config.grid_var]
    dims = grid_var.dims
    check_dims(dims, config.dim_names)
    lvls: npt.NDArray[np.float_] = np.array(config.lvls)
    da_vert = xr.DataArray(
        name=config.dim_names[1],
        data=lvls,
        coords={config.dim_names[1]: lvls},
        attrs={
            "long_name": f"{config.grid_var} levels",
            "units": ds_in.variables[config.grid_var].attrs["units"],
        },
    )
    coords = {
        config.dim_names[0]: ds_in.variables[config.dim_names[0]],
        config.dim_names[1]: da_vert,
        config.dim_names[2]: ds_in.variables[config.dim_names[2]],
        config.dim_names[3]: ds_in.variables[config.dim_names[3]],
    }
    data_vars: dict[str, xr.DataArray] = {}
    fld_var = ds_in.variables[config.fld_var]
    check_dims(fld_var.dims, dims)
    for name, intp_fld in intp_flds.items():
        assert len(intp_fld.shape) == 3, str(intp_fld.shape)
        da_fld = xr.DataArray(
            name=name,
            data=intp_fld[None, :, :, :],
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
    lvls: Sequence[float]
    dim_names: tuple[str, str, str, str]
    direction: str

    @classmethod
    def create(cls, **kwargs: Any) -> Config:
        """Create a new instance from raw arguments."""
        kwargs["dim_names"] = tuple(kwargs["dim_names"])
        return cls(**kwargs)


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
    "lvls",
    help="Level to which to interpolate field; may be repeated.",
    type=float,
    default=list(range(305, 371, 5)),
    multiple=True,
)
@click.option(
    "--dim-names",
    help="Dimension names.",
    default=["time", "lev", "lat", "lon"],
    nargs=4,
)
@click.option(
    "--direction",
    help="Direction in which vertical inpolation is performed.",
    type=click.Choice(["down", "up", "both"]),
    default="down",
)
def cli(**config_kwargs: Any) -> int:
    """Entry point from command line."""
    config = Config.create(**config_kwargs)
    print(f"read {config.grid_var}, {config.fld_var} from {config.infile}")
    with xr.open_dataset(config.infile) as ds:
        dims = ds.variables[config.grid_var].dims
        check_dims(dims, config.dim_names)
        grid = ds.variables[config.grid_var].data[0, :, :, :]
        fld = ds.variables[config.fld_var].data[0, :, :, :]
    grid = np.moveaxis(grid, 0, 2)
    fld = np.moveaxis(fld, 0, 2)
    intp_flds: dict[str, npt.NDArray[np.float_]] = {}
    for direction in ["down", "up"]:
        if config.direction in [direction, "both"]:
            print(
                f"interpolate ({direction}) {config.fld_var} to {config.grid_var} at "
                + ", ".join(map(str, config.lvls))
            )
            intp = LevelInterpolator(grid, direction=direction).to_levels(
                fld, config.lvls
            )
            name = config.fld_var
            if config.direction == "both":
                name += f"_{direction}"
            intp_flds[name] = intp
    if config.direction == "both":
        intp_flds[f"{config.fld_var}_down-up"] = (
            intp_flds[f"{config.fld_var}_down"] - intp_flds[f"{config.fld_var}_up"]
        )
    print(f"write {', '.join(intp_flds)} to {config.outfile}")
    write(intp_flds, config)
    return 0


if __name__ == "__main__":
    sys.exit(cli())
