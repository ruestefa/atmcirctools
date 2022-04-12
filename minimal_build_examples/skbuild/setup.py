"""Set up the project."""
from __future__ import annotations

# Standard library
import os
import sys

# Third-party
from setuptools import find_packages
from skbuild import setup

python_requires: str = ">=3.9"
install_requires: list[str] = []

# Obtain currently active Python version for cmake
current_python: str = f"{sys.version_info.major}.{sys.version_info.minor}"

setup(
    name="hello",
    packages=find_packages("src"),
    package_dir={"": "src"},
    python_requires=python_requires,
    install_requires=install_requires,
    cmake_args=[
        f"-DCMAKE_PREFIX_PATH={os.environ['CONDA_PREFIX']}",
        f"-DCMAKE_PYTHON_VERSION={current_python}",
    ],
)
