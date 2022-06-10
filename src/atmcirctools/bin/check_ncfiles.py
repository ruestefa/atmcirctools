"""Sanity-check NetCDF files."""
from __future__ import annotations

# Standard library
import sys
from pathlib import Path
from random import shuffle
from typing import Optional

# Third-party
import click
import netCDF4 as nc


def check_ncfile(file: Path) -> None:
    """Check NetCDF file."""
    print(file)
    nc.Dataset(file)


@click.command(
    help="Sanity-check NetCDF files in DIRS.",
)
@click.argument(
    "dirs",
    type=click.Path(
        exists=True,
        dir_okay=True,
        file_okay=False,
        path_type=Path,
    ),
    nargs=-1,
)
@click.option(
    "--pattern",
    help="Glob pattern used to match NetCDF files of interest.",
    default="*.nc",
)
@click.option(
    "-s/",
    "--shuffle/--no-shuffle",
    "do_shuffle",
    help="Shuffle files.",
)
@click.option(
    "-n",
    "--n-max",
    help="Restrict maximum number of files to check.",
    type=int,
    default=None,
)
@click.option(
    "-p",
    "--percent-max",
    help="Restrict maximum number of files to check to a percentage.",
    type=float,
    default=None,
)
def cli(
    dirs: list[Path],
    pattern: str,
    do_shuffle: bool,
    n_max: Optional[int] = None,
    percent_max: Optional[float] = None,
) -> int:
    """Run command line interface."""
    files: list[Path] = [f for d in dirs for f in d.rglob(pattern)]
    if do_shuffle:
        shuffle(files)
    if percent_max is not None:
        files = files[: int(len(files) * max(0.0, min(100.0, percent_max)) / 100)]
    if n_max is not None:
        files = files[:n_max]
    for file in files:
        check_ncfile(file)
    return 0


if __name__ == "__main__":
    sys.exit(cli())
