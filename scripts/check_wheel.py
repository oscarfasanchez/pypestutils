"""Verify the installed wheel supplies Python code and integrated native APIs."""

from importlib.metadata import version
from pathlib import Path

import numpy as np

import pypestutils
from pypestutils.gplane import GPlane
from pypestutils.pestutilslib import PestUtilsLib
from pypestutils.version import __version__


def main():
    """Require a bundled library with working GPLANE and SVA2 computations."""
    package = Path(pypestutils.__file__).resolve().parent
    checkout = Path(__file__).resolve().parent.parent
    if package == checkout / "pypestutils":
        raise RuntimeError("wheel check imported the source checkout")
    if __version__ != version("pypestutils"):
        raise RuntimeError("wheel metadata and package versions differ")
    lib = PestUtilsLib()
    library = Path(lib.pestutils._name).resolve()
    if library.parent != package / "lib":
        raise RuntimeError(f"native library is not bundled in the wheel: {library}")
    plane = GPlane(
        [[0, 0, 0], [2, 0, 0], [0, 2, 0]],
        [[0, 0], [2, 0], [0, 2]],
        (2, 2),
        (1, 1),
    )
    factors = plane.calc_factors(
        [[0.5, 0.5, 0], [1.5, 1.5, 1]], conwidth=0, a=1, lib=lib
    )
    result = factors.apply([[1, 2], [3, 4]], 10)
    np.testing.assert_allclose(result, [1, 10 + np.exp(-1) * (4 - 10)])
    # A single contributing source has unit normalized weight: mean + iid*sd.
    field = lib.fieldgen3d_sva2(
        [0], [0], [0], 1.0, [2.0],
        [0, 0.5], [0, 0], [0, 0], 1, 10.0, 4.0,
        1.0, 1.0, 1.0, 0.0, 0.0, 0.0, "none", "exp", 1.0,
    )  # fmt: skip
    np.testing.assert_allclose(field, [14.0, 14.0])
    print(f"Installed wheel verified: pypestutils {__version__}; {library.name}")


if __name__ == "__main__":
    main()
