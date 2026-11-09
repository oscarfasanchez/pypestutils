# PyPestUtils

## GPLANE / MODFLOW-USG / DISU fork

This fork includes [3D GPLANE structural overlays](docs/gplane.md),
MODFLOW-USG and MODFLOW 6 DISU support. See [FORK_CHANGES.md](FORK_CHANGES.md)
and the [PLPROC comparisons and benchmarks](benchmarks/README.md).

Prebuilt wheels are available in the
[v0.3.2-usg.1 GitHub Release](https://github.com/oscarfasanchez/pypestutils/releases/tag/v0.3.2-usg.1)
with package version `0.3.2+usg.1`. No compiler or repository checkout is needed.

Install with pip (the requirements file selects the platform wheel):

```sh
python -m pip install -r https://github.com/oscarfasanchez/pypestutils/releases/download/v0.3.2-usg.1/requirements-usg.txt
```

For Conda, download [environment-usg.yml](environment-usg.yml) and run
`conda env create -f environment-usg.yml`, or put this in your own file:

```yaml
dependencies:
  - python=3.12
  - pip
  - pip:
      - "-r https://github.com/oscarfasanchez/pypestutils/releases/download/v0.3.2-usg.1/requirements-usg.txt"
```

The automatic requirements support Windows x64, Linux glibc 2.27+ on x64 or
ARM64, and macOS 15+ on Intel or Apple Silicon. Alpine/musl wheels are separate
release assets. See [RELEASE_NOTES.md](RELEASE_NOTES.md) for direct wheel URLs
and validation. The package and import name remain `pypestutils`; separate
environments can retain upstream and fork installations.

Suite of Python tools based on PEST utilities to support parameterization with pilot points, drawing stochastic realizations, and MODFLOW post-processing for structured and unstructured grids. This package consists of a (pre-)compiled shared fortran library, and a corresponding low-level python wrapper around the library functions.  There are also higher-level "helper" functions to further abstract the granular low-level function calls for common workflow elements.  

The low-level functions are relatively strict in their argument types - this is required to effectively pass the array-type data container references to the fortran library.  As such, the low-level python functions perform considerable type checking.  The higher-level helpers attempt to support a variety of argument types and will attempt coerce to the required type.  

This package is currently in pre-alpha development, and is not suitable for use, but early adopters are welcome to have a go.

## Examples

There are several jupyter notebook examples of using pypestutils for a structured and quadtree Freyberg model.  These notebooks rely on both [flopy](https://github.com/modflowpy/flopy) and [pyEMU](https://github.com/pypest/pyemu) to help with visualization and processing.

The use the low-level python interface to the shared fortran library, you create a `PESTUTILSLIB` instance and then can directly call the shared library routines:

```
from pypestutils.pestutilslib import PestUtilsLib
lib = PestUtilsLib() #the constructor searches for the shared lib
grid_info = lib.install_mf6_grid_from_file("grid","freyberg6.disv.grb")
easting,northing,elev = lib.get_cell_centres_mf6("grid",grid_info["ncells"])
```

The higher-level helper functions obsecure the calls the fortran library and string together multiple low-level function calls:

```
import pypestutils.helpers as helpers
grid_info = helpers.get_grid_info_from_file("freyberg6.disv.grb")
```

## Documentation

The documentation for pypestutils can be found [here](https://pypest.github.io/pypestutils/)

The documentation for the shared FORTRAN library can be found [here](docs/pestutilslib/fortran_library_documentation.md)

## Installation

### Dependencies

pypestutils requires `numpy` and `pandas`

### Upstream PyPI release

Use `pip` to install a built distribution for Windows, Linux or macOS:

    pip install pypestutils

to also include optional requirements use:

    pip install pypestutils[optional]

### From source

Installation from source requires a Fortran compiler and build tools. See [BUILD.md](BUILD.md) for details.

## Disclaimer


This software is preliminary or provisional and is subject to revision. It is
being provided to meet the need for timely best science. This software is
provided "as is" and "as-available", and makes no representations or warranties
of any kind concerning the software, whether express, implied, statutory, or
other. This includes, without limitation, warranties of title,
merchantability, fitness for a particular purpose, non-infringement, absence
of latent or other defects, accuracy, or the presence or absence of errors,
whether or not known or discoverable.
