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


def check_ncfile(  # noqa: max-complexity=13
    file: Path,
    style: str = "human",
    verbosity: int = 1,
) -> None:
    """Check NetCDF file."""
    if style in ["h", "human"]:
        plain = False
    elif style in ["p", "plain"]:
        plain = True
    else:
        raise ValueError(f"invalid style: '{style}'")
    if not plain and verbosity > 0:
        print(f"[      ] {file}", end="", flush=True)
    try:
        nc.Dataset(file)
    except Exception as e:
        if plain and verbosity == 0:
            print(file)
        elif plain and verbosity == 1:
            print(f"{file} {type(e).__name__}")
        elif plain and verbosity > 1:
            print(f"F {file} {type(e).__name__}")
        elif not plain:
            print(f"\r[ FAIL ] {file} ({type(e).__name__})")
    else:
        if plain and verbosity > 1:
            print(f"- {file} -")
        elif not plain and verbosity <= 1:
            print("\r" + (len("[      ] ") + len(str(file))) * " " + "\r", end="")
        elif not plain and verbosity > 1:
            print(f"\r[  OK  ] {file}")


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
@click.option(
    "-o",
    "--output-style",
    help="Style of output.",
    type=click.Choice(["h", "human", "p", "plain"]),
    default="human",
)
@click.option(
    "-v/",
    "--verbose",
    "verbosity",
    help="Increase verbosity; may be repeated.",
    count=True,
)
def cli(
    dirs: list[Path],
    pattern: str,
    do_shuffle: bool,
    n_max: Optional[int],
    percent_max: Optional[float],
    output_style: str,
    verbosity: int,
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
        check_ncfile(file, style=output_style, verbosity=verbosity)
    print("\r")
    return 0


if __name__ == "__main__":
    sys.exit(cli())
