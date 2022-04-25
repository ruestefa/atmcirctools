"""Main command line entry point to atmcirctools."""
from __future__ import annotations

# Standard library
import sys

# Third-party
import click
from atmcirclib.click import Group

# Local
from .intp_lvl import cli as intp_lvl_cli
from .shared import CONTEXT_SETTINGS


@click.group(
    cls=Group,
    context_settings=CONTEXT_SETTINGS,
    help="AtmCircTools: A collection of tools useful for atmospheric science",
)
def cli() -> int:
    """Run umbrella command ``act`` of ``AtmCircTools`` from command line."""
    return 0


cli.add_command(intp_lvl_cli, name="intp-lvl")


if __name__ == "__main__":
    sys.exit(cli())
