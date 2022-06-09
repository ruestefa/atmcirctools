"""Sanity-check NetCDF files."""
from __future__ import annotations

# Standard library
import sys

# Third-party
import click


@click.command()
def cli() -> int:
    """Run command line interface."""
    return 0


if __name__ == "__main__":
    sys.exit(cli())
