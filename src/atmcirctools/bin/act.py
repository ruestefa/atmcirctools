"""Main command line entry point to atmcirctools."""
from __future__ import annotations

# Standard library
import sys

# Third-party
import click
from atmcirclib.click import Group

# Local
from . import call_graph
from . import check_nc
from . import create_recipe
from . import intp_lvl
from .shared import CONTEXT_SETTINGS


@click.group(
    cls=Group,
    context_settings=CONTEXT_SETTINGS,
    help="AtmCircTools: A collection of tools useful for atmospheric science",
)
def cli() -> int:
    """Run umbrella command ``act`` of ``AtmCircTools`` from command line."""
    return 0


cli.add_command(call_graph.cli, name="call-graph")
cli.add_command(check_nc.cli, name="check-nc")
cli.add_command(create_recipe.cli, name="create-recipe")
cli.add_command(intp_lvl.cli, name="intp-lvl")


if __name__ == "__main__":
    sys.exit(cli())
